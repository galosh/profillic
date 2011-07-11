/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profillic
##  File:
##      Train.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The profillic program.  It takes a some unaligned sequences and
##      estimates parameters of a profile HMM, and then outputs the profile HMM
##      in the galosh profile format.
##
#******************************************************************************
#*
#*    This file is part of profillic, a suite of programs for estimating parameters of
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2009, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *        http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
#*****************************************************************************/

#include "Train.hpp"

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace galosh;

int
main ( int const argc, char const ** argv )
{
  //typedef bfloat ProbabilityType;
  //typedef logspace ProbabilityType;
  //typedef floatrealspace ProbabilityType;
  typedef doublerealspace ProbabilityType;
  
  typedef bfloat ScoreType; // Preferred
  //typedef logspace ScoreType; // SLOWer than bfloat
  //typedef realspace ScoreType; // Only for very few & small sequences
  
  // if using anything other than LogProbability for the MatrixValueType,
  // params.useRabinerScaling should be set to true.
  typedef bfloat MatrixValueType;
  //typedef logspace MatrixValueType;
  //typedef doublerealspace MatrixValueType;
  //typedef floatrealspace MatrixValueType;


#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  string fasta_filename;
  string profile_filename;
  string profile_output_filename;
  string * profile_output_filename_ptr = NULL; // NULL means don't output the profile to its own file.
  int use_unconditional_bw = 0; // Use Conditional BW (or QA) unless 3rd arg is non-zero.
  unsigned int initial_profile_length = 0; // 0 means use the median length
  double lengthadjust_deletion_threshold = .5;//.15;
  double lengthadjust_insertion_threshold = lengthadjust_deletion_threshold;// / seqan::ValueSize<ResidueType>::VALUE;
  double lengthadjust_threshold_increment = .0005;
  bool use_lengthadjust = true;
  double lengthadjust_occupancy_threshold = .5; // not used unless count_seqs_exceeding_alignment_profile_thresholds is true in ProfileTrainer.hpp
  int lengthadjust_use_sensitive_thresholding = 0;
  unsigned int lengthadjust_increase_thresholds_for_length_changes_start_iteration = 500;
  double lengthadjust_increase_thresholds_for_length_changes_min_increment = 1E-4;
  unsigned long random_seed = 0; // 0 means use the time. (NOTE: Presently ignored since profile is started at even() instead of uniform().
  double even_starting_profile_multiple = 0; // If < 0, even() the positions of the starting profile; if > 0, starting profile's positions are averaged with this multiple of the even profile. (NOTE: Presently irrelevant; see above about starting at even() by default. TODO: Change so that if this argument is given, uniform() *is* used..)

  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <output profile file> <fasta sequences file>  [<input profile filename>|<initial profile length or 0 to use median or 1 to use max>|0 [<use Unconditional training>|0 [<lengthadjust insertion threshold or 0 to disable changes to the profile length>|.5 [<lengthadjust deletion threshold>|.5 [<lengthadjust threshold increment>|.0005 [<increase thresholds for length changes after iteration>|500 [<random seed>|0 [<even starting profile multiple (or startEven, if negative)>|0 ]...]" << endl;
    cout << "Example usage (typical; with DMS starting from median sequence length): " << argv[ 0 ] << " <output profile file> <fasta sequences file>" << endl;
    cout << "Example usage with Unconditional training and DMS: " << argv[ 0 ] << " <output profile file> <fasta sequences file> 0 1" << endl;
    cout << "Example usage with Unconditional training and DMS, starting from maximum sequence length: " << argv[ 0 ] << " <output profile file> <fasta sequences file> 1 1" << endl;
    cout << "Example usage with fixed profile length (no DMS): " << argv[ 0 ] << " <output profile file> <fasta sequences file> <initial profile length> 0 0" << endl;
    cout << "Example usage with a starting profile and a fixed profile length: " << argv[ 0 ] << " <output profile file> <fasta sequences file> <input profile filename> 0 0" << endl;
    exit( 1 );
  }
  // At least two arguments: profile_output_filename and fasta_filename
  profile_output_filename = argv[ 1 ];
  if( ( profile_output_filename == "" ) || ( profile_output_filename == "-" ) ) {
    profile_output_filename_ptr = NULL; // Don't output a profile.
  } else {
    profile_output_filename_ptr = &profile_output_filename;
    // TODO: REMOVE
    cout << "Output profile filename: " << profile_output_filename << endl;
  }
  fasta_filename = argv[ 2 ];

  if( argc > 3 ) {
    // At least 3 arguments: 3rd is either profile_filename or initial_profile_length
    try {
      // TODO: REMOVE
      //cout << "argv[ 3 ] is " << argv[ 3 ] << endl;
      initial_profile_length = boost::lexical_cast<int>( argv[ 3 ] );
      //cout << "Initial profile length: " << initial_profile_length << endl;
    } catch( boost::bad_lexical_cast & ) {
      // TODO: REMOVE
      //cout << "it is NOT an integer." << endl;
      // Ok maybe it is a filename.
      std::ifstream fs( argv[ 3 ] );
      if( fs.is_open() ) {
        // Okay, we'll assume it's the profile then.
        profile_filename = argv[ 3 ];
        initial_profile_length = 0;
        // TODO: REMOVE
        cout << "Starting profile filename: " << profile_filename << endl;
        //cout << "profile_filename.length() returns " << profile_filename.length() << endl;
      } else {
        std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as an integer for use as the " << ( use_lengthadjust ? "initial " : "" ) << "profile length.  It also does not seem to be a filename (if it did, we would think it's the initial profile file)." << std::endl;
        exit( 1 );
      }
    } // End try .. catch block for lexical_cast
  }
  if( argc > 4 ) {
    // At least 4 arguments: 4th indicates useUnconditionalBaumWelch if non-zero
    try {
      use_unconditional_bw =
        boost::lexical_cast<int>( argv[ 4 ] );
      if( use_unconditional_bw ) {
        std::cout << "Using Unconditional training." << std::endl;
      } else {
        std::cout << "Using Conditional training." << std::endl;
      }
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 4 ] << "' as an int value for use as a flag to indicate use of Unconditional training." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // At least 4 arguments: 4th is use_unconditional_bw
  if( argc > 5 ) {
    // At least 5 arguments: 5th is lengthadjust insertion threshold (or 0 for no lenthadjust)
    try {
      lengthadjust_insertion_threshold = boost::lexical_cast<double>( argv[ 5 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 5 ] << "' as a real value for use as the lengthadjust threshold (or as 0 to indicate no lengthadjust)." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( lengthadjust_insertion_threshold < 0 ) {
      std::cerr << "The given lengthadjust insertion threshold value, " << lengthadjust_insertion_threshold << ", is negative.  You must supply a value between 0 and 1 (or 0 to indicate no lengthadjust)." << std::endl;
      exit( 1 );
    }
    if( lengthadjust_insertion_threshold > 1 ) {
      std::cerr << "The given lengthadjust insertion threshold value, " << lengthadjust_insertion_threshold << ", is greater than 1.  You must supply a value between 0 and 1 (or 0 to indicate no lengthadjust)." << std::endl;
      exit( 1 );
    }
    if( lengthadjust_insertion_threshold == 0 ) {
      std::cout << "Lengthadjust disabled." << endl;
      use_lengthadjust = false;
    } else {
      use_lengthadjust = true;
      //std::cout << "Initial lengthadjust insertion threshold: " << lengthadjust_insertion_threshold << std::endl;
    }
  }
  if( argc > 6 ) {
    // At least 6 arguments: 6th is lengthadjust deletion threshold (or 0 to use lengthadjust_insertion_threshold)
    try {
      lengthadjust_deletion_threshold = boost::lexical_cast<double>( argv[ 6 ] );
      //if( use_lengthadjust ) {
      //  if( lengthadjust_deletion_threshold ) {
      //    std::cout << "Initial lengthadjust deletion threshold: " << lengthadjust_deletion_threshold << std::endl;
      //  } else {
      //    std::cout << "Initial lengthadjust deletion threshold: " << lengthadjust_insertion_threshold << std::endl;
      //  }
      //}
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 6 ] << "' as a real value for use as the deletion lengthadjust threshold (or as 0 to use insertion lengthadjust threshold)." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( lengthadjust_deletion_threshold < 0 ) {
      std::cerr << "The given lengthadjust deletion threshold value, " << lengthadjust_deletion_threshold << ", is negative.  You must supply a value between 0 and 1 (or 0 to use insertion lengthadjust threshold)." << std::endl;
      exit( 1 );
    }
    if( lengthadjust_deletion_threshold > 1 ) {
      std::cerr << "The given lengthadjust deletion threshold value, " << lengthadjust_deletion_threshold << ", is greater than 1.  You must supply a value between 0 and 1 (or 0 to indicate use insertion lengthadjust threshold)." << std::endl;
      exit( 1 );
    }
    if( lengthadjust_deletion_threshold == 0 ) {
      lengthadjust_deletion_threshold = lengthadjust_insertion_threshold;
    }
  }
  if( argc > 7 ) {
    // At least 7 arguments: 7th is lengthadjust threshold increment
    try {
      lengthadjust_threshold_increment = boost::lexical_cast<double>( argv[ 7 ] );
      //if( use_lengthadjust ) {
      //  std::cout << "Lengthadjust threshold increment: " << lengthadjust_threshold_increment << std::endl;
      //}
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 7 ] << "' as a real value for use as the lengthadjust threshold increment." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( lengthadjust_threshold_increment <= 0 ) {
      std::cerr << "The given lengthadjust threshold increment value, " << lengthadjust_threshold_increment << ", is not positive.  You must supply a positive value between 0 and 1." << std::endl;
      exit( 1 );
    }
    if( lengthadjust_threshold_increment > 1 ) {
      std::cerr << "The given lengthadjust threshold increment value, " << lengthadjust_threshold_increment << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
      exit( 1 );
    }
    if( lengthadjust_threshold_increment == 0 ) {
      std::cerr << "WARNING: Lengthadjust threshold incrementing disabled." << endl;
    }
  }
  if( argc > 8 ) {
    // At least 8 arguments: 8th is increase thresholds for length changes start iteration
    try {
      lengthadjust_increase_thresholds_for_length_changes_start_iteration = boost::lexical_cast<unsigned long>( argv[ 8 ] );
      std::cout << "Increase thresholds for length changes start iteration: " << lengthadjust_increase_thresholds_for_length_changes_start_iteration << std::endl;
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 8 ] << "' as an unsigend long value for use as the iteration at which we start increasing thresholds for all length changes." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  }
  if( argc > 9 ) {
    // At least 9 arguments: 9th is the random seed
    try {
      random_seed = boost::lexical_cast<unsigned long>( argv[ 9 ] );
      std::cout << "Random seed: " << random_seed << std::endl;
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 9 ] << "' as an unsigend long value for use as the random seed." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  }
  if( argc > 10 ) {
    // At least 10 arguments: 10th is even_starting_profile_multiple (If < 0, even() the positions of the starting profile; if > 0, starting profile's positions are averaged with this multiple of the even profile.).  NOTE that this is irrelevant if not starting from a profile in a file, since it'll be even always otherwise.
    try {
      even_starting_profile_multiple =
        boost::lexical_cast<double>( argv[ 10 ] );
      if( even_starting_profile_multiple < 0 ) {
        std::cout << "Starting from evenly-distributed profile positions." << std::endl;
      } else if( even_starting_profile_multiple > 0 ) {
        std::cout << "Mixing starting positions with " << even_starting_profile_multiple << " times the evenly-distributed profile." << std::endl;
      }
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 10 ] << "' as a double value for use as the even starting profile multiple." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  }
  
  Train<ProbabilityType, ScoreType, MatrixValueType, ResidueType, SequenceResidueType> train;

  if( profile_filename.length() == 0 ) {
    // TODO: REMOVE
    //cout << "Training using a randomly generated profile." << endl;
    cout << "Training using the default starting profile." << endl; // Not random; even positions..
    train.train(
      fasta_filename,
      initial_profile_length,
      profile_output_filename_ptr,
      lengthadjust_insertion_threshold,
      lengthadjust_deletion_threshold,
      lengthadjust_threshold_increment, //( lengthadjust_threshold_increment / ( lengthadjust_deletion_threshold / lengthadjust_insertion_threshold ) ),
      lengthadjust_threshold_increment,
      random_seed,
      lengthadjust_occupancy_threshold,
      ( bool )lengthadjust_use_sensitive_thresholding,
      lengthadjust_increase_thresholds_for_length_changes_start_iteration,
      lengthadjust_increase_thresholds_for_length_changes_min_increment,
      ( bool )use_unconditional_bw,
      even_starting_profile_multiple,
      false // don't train globals first
    );
  } else { // if we were given a profile_filename .. else ..
    // TODO: REMOVE
    cout << "Training using the profile in file \"" << profile_filename << "\"." << endl;
    train.train(
      fasta_filename,
      profile_filename,
      profile_output_filename_ptr,
      lengthadjust_insertion_threshold,
      lengthadjust_deletion_threshold,
      lengthadjust_threshold_increment, //( lengthadjust_threshold_increment / ( lengthadjust_deletion_threshold / lengthadjust_insertion_threshold ) ),
      lengthadjust_threshold_increment,
      lengthadjust_occupancy_threshold,
      ( bool )lengthadjust_use_sensitive_thresholding,
      lengthadjust_increase_thresholds_for_length_changes_start_iteration,
      lengthadjust_increase_thresholds_for_length_changes_min_increment,
      ( bool )use_unconditional_bw,
      even_starting_profile_multiple,
      false // no, don't /////true // do train globals first
    );
  } // End if we were given a profile_filename .. else ..

  return 0; // success
} // main (..)

