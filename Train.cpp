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
namespace po = boost::program_options;

int
main ( int const argc, char ** argv )
{
  //typedef bfloat ProbabilityType;
  //typedef logspace ProbabilityType;
  //typedef floatrealspace ProbabilityType;
  typedef doublerealspace ProbabilityType;
  
  typedef bfloat ScoreType; // Preferred
  //typedef logspace ScoreType; // SLOWer than bfloat
  //typedef realspace ScoreType; // Only for very few & small sequences
  
  // You must use either logspace or bfloat for the MatrixValueType, unless you're debugging with very few and small training sequences.
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

  try {
    string config_file;
    
    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic( "Generic options" );
    generic.add_options()
      ( "version,v", "print version string" )
      ( "help", "produce help message" )
      ( "config,c", po::value<string>( &config_file )->default_value( "profillic.cfg" ),
        "name of a file of a configuration." )
      ;

    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config( "Configuration" );
    config.add_options()
      ( "output_profile,p", 
        po::value<string>(),
        "filename: where to put the trained output profile" )
      ( "fasta,f", 
        po::value<string>(),
        "input sequences, in (unaligned) Fasta format" )
      ( "starting_profile,s", 
        po::value<string>(),
        "optional starting-profile input filename; note that you can read in a profile or start with a random profile, then reset its emission probabilities to even or mix them with some number of multiples of an even profile: see even_starting_profile_multiple." )
      ( "initial_profile_length,l",
        po::value<uint32_t>(),
        "length of the starting profile (only used if the starting_profile is not specified); special value 0 means use the median length of the input sequences and 1 means use the maximum length (default is 1); the starting profile will begin with uniformly-determined emission probabilities, but by default it'll be reset to even before training -- see even_starting_profile_multiple" )
      ( "random_seed,r",
        po::value<uint32_t>(),
        "Random seed to use for the initial profile; defaults to 0, which means use the current time; note that this is generally irrelevant unless you explicitly set even_starting_profile_multiple" )
      ( "even_starting_profile_multiple,e",
        po::value<double>(),
        "If < 0, reset the emission probabilities of the starting profile to even; if > 0, mix the starting profile with this number of multiples of an even profile (defaults to 0 when specifying a starting profile filename, -1 otherwise)." )
      ( "nseq,n",
        po::value<uint32_t>(),
        "number of sequences to use (default is ALL)" )
      ( "proposeProfileLengthChanges,use_lengthadjust,dms,d", "use the Dynamic Model Surgery (DMS) algorithm to prune unused model states and introduce states for oft-used insertions" )
      ( "useUnconditionalBaumWelch,unconditional,u", "use the unconditional variant of the algorithm to simultaneously update all positions in each iteration" )
      ( "trainGlobalsFirst,globals_first,g", "train the globals before the first positions-training step" )
      ( "maxIterations,max_iterations,i",
        po::value<uint32_t>(),
        "maximum number of iterations to use in training the profile (default 1000)" )
      ( "priorStrength,emisions_prior_strength",
        po::value<double>(),
        "strength of priors for emissions parameters (not used for Quadratic Ascent training)" )
      ( "priorStrength_internal_transitions,transitions_prior_strength",
        po::value<double>(),
        "strength of priors for transitions parameters" )
      ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description lengthadjust_opts( "Lengthadjust options" );
    lengthadjust_opts.add_options()
      ( "proposeInsertingOccupancyThreshold,dms.occupancy_threshold",
        po::value<double>(),
        "DMS threshold for fraction of sequences inserting/deleting to trigger a model edit of a given position" )
      ( "proposeInsertingThreshold,dms.insertion_threshold",
        po::value<double>(),
        "DMS per-sequence threshold on insertion-open to count that sequence as inserting at a given position" )
      ( "proposeDeletingThreshold,dms.deletion_threshold",
        po::value<double>(),
        "DMS per-sequence threshold on the probability of deletion at a position to count that sequence as deleting the position" )
      ( "proposeInsertingThreshold_increment,dms.insertion_threshold_increment",
        po::value<double>(),
        "DMS insertion threshold increment" )
      ( "proposeDeletingThreshold_increment,dms.deletion_threshold_increment",
        po::value<double>(),
        "DMS deletion threshold increment" )
      ( "increaseThresholdsForLengthChanges_startIteration,dms.increase_thresholds_for_length_changes_start_iteration",
        po::value<uint32_t>(),
        "Eventually we get impatient and want the DMS process to hurry along; this determines the iteration at which we start increasing thresholds using the threshold increments and the min_increment" )
      ( "increaseThresholdsForLengthChanges_minIncrement,dms.increase_thresholds_for_length_changes_min_increment",
        po::value<double>(),
        "When increasing DMS thresholds for length changes, apply at least this minimum increment per iteration" )
      ;

    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
    ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters params;
    
    po::options_description cmdline_options;
    cmdline_options.add( generic ).add( params.m_galosh_options_description ).add( config ).add( lengthadjust_opts );

    po::options_description config_file_options;
    config_file_options.add( params.m_galosh_options_description ).add( config ).add( lengthadjust_opts );

    po::options_description visible( "Basic options" );
    visible.add( generic ).add( config );

    po::positional_options_description p;
    p.add( "output_profile", 1 );
    p.add( "fasta", 1 );
        
    store( po::command_line_parser( argc, argv ).options( cmdline_options ).positional( p ).run(), params.m_galosh_options_map );
    notify( params.m_galosh_options_map );

    // TODO: REMOVE
    //cout << params << endl;
    //cout << endl;

#define USAGE() " " << argv[ 0 ] << " [options] <output profile file> <fasta sequences file>"

        
//    cout << "Usage: " << argv[ 0 ] << " <output profile file> <fasta sequences file>  [<input profile filename>|<initial profile length or 0 to use median or 1 to use max>|0 [<use Unconditional training>|0 [<lengthadjust insertion threshold or 0 to disable changes to the profile length>|.5 [<lengthadjust deletion threshold>|.5 [<lengthadjust threshold increment>|.0005 [<increase thresholds for length changes after iteration>|500 [<random seed>|0 [<even starting profile multiple (or startEven, if negative)>|0 ]...]" << endl;
//    cout << "Example usage (typical; with DMS starting from median sequence length): " << argv[ 0 ] << " <output profile file> <fasta sequences file>" << endl;
//    cout << "Example usage with Unconditional training and DMS: " << argv[ 0 ] << " <output profile file> <fasta sequences file> 0 1" << endl;
//    cout << "Example usage with Unconditional training and DMS, starting from maximum sequence length: " << argv[ 0 ] << " <output profile file> <fasta sequences file> 1 1" << endl;
//    cout << "Example usage with fixed profile length (no DMS): " << argv[ 0 ] << " <output profile file> <fasta sequences file> <initial profile length> 0 0" << endl;
//    cout << "Example usage with a starting profile and a fixed profile length: " << argv[ 0 ] << " <output profile file> <fasta sequences file> <input profile filename> 0 0" << endl;


    // Read in the config file.
    if( config_file.length() > 0 ) {
      ifstream ifs( config_file.c_str() );
      if( !ifs ) {
        if(!params.m_galosh_options_map["config"].defaulted()) {         //TAH 3/13 don't choke if config file was defaulted and is missing
           cout << "Can't open the config file named \"" << config_file << "\"\n";
           return 0;
        } 
      } else {
        store( parse_config_file( ifs, config_file_options ), params.m_galosh_options_map );
        notify( params.m_galosh_options_map );
      }
    }

    if( params.m_galosh_options_map.count( "help" ) > 0 ) {
      cout << "Usage: " << USAGE() << endl;
      cout << visible << "\n";
      return 0;
    }

    if( params.m_galosh_options_map.count( "version" ) ) {
      cout << "Profillic, version 1.0\n";
      return 0;
    }

    //if( params.m_galosh_options_map.count( "debug" ) ) {
    //  cout << "[DEBUGGING]\n";
    //  return 0;
    //}

    // Required options
    if( ( params.m_galosh_options_map.count( "output_profile" ) == 0 ) || ( params.m_galosh_options_map.count( "fasta" ) == 0 ) ) {
      cout << "Usage: " << USAGE() << endl;
      return 1;
    }
    
    //if( params.m_galosh_options_map.count( "include-path" ) ) {
    //  cout << "Include paths are: " 
    //       << params.m_galosh_options_map["include-path"].as< vector<string> >() << "\n";
    //}
    //
    //    if (params.m_galosh_options_map.count("input-file"))
    //    {
    //        cout << "Input files are: " 
    //             << params.m_galosh_options_map["input-file"].as< vector<string> >() << "\n";
    //    }

    // TODO: ERE I AM
    Train<ProbabilityType, ScoreType, MatrixValueType, ResidueType, SequenceResidueType> train;
    train.train( params );

    return 0;
  } catch( std::exception& e ) { /// exceptions thrown by boost stuff (etc)
    cerr << "error: " << e.what() << endl;
    return 1;
  } catch( string &err ) {      /// exceptions thrown as strings
    cerr << "error: " << err << endl;
    return 1;
  } catch( ... ) {               /// anything else
    cerr << "Strange unknown exception" << endl;
    return 1;
  }

    /// OLD:
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
  //train.train( params.m_galosh_options_map );
//
//  if( profile_filename.length() == 0 ) {
//    // TODO: REMOVE
//    //cout << "Training using a randomly generated profile." << endl;
//    cout << "Training using the default starting profile." << endl; // Not random; even positions..
//    train.train(
//      fasta_filename,
//      initial_profile_length,
//      profile_output_filename_ptr,
//      lengthadjust_insertion_threshold,
//      lengthadjust_deletion_threshold,
//      lengthadjust_threshold_increment, //( lengthadjust_threshold_increment / ( lengthadjust_deletion_threshold / lengthadjust_insertion_threshold ) ),
//      lengthadjust_threshold_increment,
//      random_seed,
//      lengthadjust_occupancy_threshold,
//      ( bool )lengthadjust_use_sensitive_thresholding,
//      lengthadjust_increase_thresholds_for_length_changes_start_iteration,
//      lengthadjust_increase_thresholds_for_length_changes_min_increment,
//      ( bool )use_unconditional_bw,
//      even_starting_profile_multiple,
//      false // don't train globals first
//    );
//  } else { // if we were given a profile_filename .. else ..
//    // TODO: REMOVE
//    cout << "Training using the profile in file \"" << profile_filename << "\"." << endl;
//    train.train(
//      fasta_filename,
//      profile_filename,
//      profile_output_filename_ptr,
//      lengthadjust_insertion_threshold,
//      lengthadjust_deletion_threshold,
//      lengthadjust_threshold_increment, //( lengthadjust_threshold_increment / ( lengthadjust_deletion_threshold / lengthadjust_insertion_threshold ) ),
//      lengthadjust_threshold_increment,
//      lengthadjust_occupancy_threshold,
//      ( bool )lengthadjust_use_sensitive_thresholding,
//      lengthadjust_increase_thresholds_for_length_changes_start_iteration,
//      lengthadjust_increase_thresholds_for_length_changes_min_increment,
//      ( bool )use_unconditional_bw,
//      even_starting_profile_multiple,
//      false // no, don't /////true // do train globals first
//    );
//  } // End if we were given a profile_filename .. else ..

  return 0; // success
} // main (..)

