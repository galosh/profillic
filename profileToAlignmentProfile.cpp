/**
 * @file profileToAlignmentProfile.cpp
 *
 * @author  Ted Holzman <tholzman@scharp.org>
 *  modified from the Align.cpp file by Paul Edlefsen <pedlefsen@gmail.com>
 * @version 0.1
 * @date 02/01/2012
 * @section LICENSE
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * \brief Create one (or more) "Alignment Profiles" from a fasta file and a
 * a profile
 *
 * @section DESCRIPTION
 * Part of the \a profuse library.
 *
 * Converts a "standard" galosh profile to an alignment profile containing
 * positional information about insertions and deletions.  It gleans this
 * information from a fasta file containing several sequences.  profileToAlignmentProfile
 * aligns these sequences, determines the indels, and combines the information
 * with the given profile to produce an Alignment Profile.
 *
 * The following options are understood:
 *<p><pre>
 * Usage:  profileToAlignmentProfile [options] <profile-file-name> <fasta-file-name>
 *
 * profileToAlignmentProfile options:
 *
 * -h [ --help ]                 output this help message
 * --verbosity arg (=0)          set verbosity level (0 - 3)
 * -i [ --individual ]           output individual alignment profiles instead of
 *                               average
 * -f [ --format ] arg (=pileup) output format of alignments (verbosity>=3)
 * -n [ --nseq ] arg             number of sequences to use (default is ALL)
 * -v [ --viterbi ]              use viterbi algorithm
 * </pre>
 *
 */

#include <boost/program_options.hpp>

#include "GenAlignmentProfiles.hpp"

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace galosh;

/**
 * \fn int main(int const argc, char const ** argv)
 * \brief main driver.  Parses command line and calls gen_alignment_profiles
 * prints alignment profile output.
 */
int
main ( int const argc, char const ** argv )
{
  typedef doublerealspace ProbabilityType;
  typedef bfloat ScoreType;
  typedef bfloat MatrixValueType;

#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..
  /**
   * This section is for parsing commandline options.
   * \var po::options_description desc
   * \brief is for \a named options.
   * \var po::positional_options_description pdesc
   * \brief is for \a positional options
   * \var po::variables_map vm
   * \brief is a map that holds all the parsed commandline info
   * \todo Add ability to read options from a file
   */
#define USAGE() " " << argv[0] << " [options] <profile-file-name> <fasta-file-name>"
  namespace po = boost::program_options;
  po::options_description desc("profileToAlignmentProfile options");
  desc.add_options()
      ("help,h", "output this help message")
      ("verbosity", po::value<int>()->default_value(0), "set verbosity level (0 - 3)")
      ("individual,i","output individual alignment profiles instead of average")
      ("format,f",po::value<std::string>()->default_value("pileup"),"output format of alignments (verbosity>=3)")
      ("nseq,n",po::value<int>(),"number of sequences to use (default is ALL)")
      ("viterbi,v","use viterbi algorithm")
      ("profile",po::value<string>(),"Name of file containing galosh profile")
      ("fasta",po::value<string>(),"Name of fasta file containing sequences")
  ;
  po::positional_options_description pdesc;
  pdesc.add("profile",1);
  pdesc.add("fasta",1);
  po::variables_map vm;
  try
  {
    store(po::command_line_parser(argc,argv).options(desc).positional(pdesc).run(),vm);
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
     
    /// if --help or -h was specified, give help and normal exit
    if( vm.count("help") )
    {
      cout << "Usage: " << USAGE() << endl;
      cout << desc << endl;
      return 0;
    }

    /// If either the profile file or the fasta file wasn't specified
    /// give short usage message and error exit
    if( !vm.count("profile") || !vm.count("fasta") )
    {
      cerr << "Usage: " << USAGE() << endl;
      return( 1 );
    }
    
    /// Do the work
    GenAlignmentProfiles<ProbabilityType, ScoreType, MatrixValueType, ResidueType, SequenceResidueType> genAlignProf;
    std::vector<DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> alignment_profiles;
    alignment_profiles = genAlignProf.gen_alignment_profiles( vm );
    
    /// Output the results
    /// \todo Make the profile scanner skip "#" lines
    /// \todo Label lines of individual profiles with something
    for( int i = 0; i < alignment_profiles.size(); i++ )
    {
      std::cout << alignment_profiles[ i ];
      if( alignment_profiles.size() > 1 ) {
      std:cout << "#" << std::endl;
      }
    }
    
    return 0; // success
  } catch( std::exception& e ) { /// exceptions thrown by boost stuff
    cerr << "error: " << e.what() << endl;
    return 1;
  } catch( string &err ) {      /// exceptions thrown by GenAlignmentProfiles, etc.
    cerr << "error: " << err << endl;
    return 1;
  } catch( ... ) {               /// anything else
    cerr << "Strange unknown exception" << endl;
    return 1;
  }
} // main (..)

