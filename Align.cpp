/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      Align.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The align program.  It takes a profile and some unaligned sequences and
##      computes the viterbi alignment.
##
#******************************************************************************
#*
#*    This file is part of profuse, a suite of programs for working with
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
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

#include "ScoreAndMaybeAlign.hpp"

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

  string profile_filename;
  string fasta_filename;
  if( argc < 3 ) {
    cout << "Usage: " << argv[ 0 ] << " <profile file> <fasta sequences file> [<number of sequences to use>]" << endl;
    exit( 1 );
  }
  // else { // ( argc >= 3 )
  // At least two arguments: profile filename and fasta filename
  profile_filename = argv[ 1 ];
  fasta_filename = argv[ 2 ];
  uint32_t sequence_count = 0; // 0 means use all of the seqs in the fasta file.
  if( argc > 3 ) {
    try {
      sequence_count = boost::lexical_cast<uint32_t>( argv[ 3 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as a count value for use as the number of sequences to use (or as 0 to that all sequences should be used)." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // End if argc > 3

  ScoreAndMaybeAlign<ProbabilityType, ScoreType, MatrixValueType, ResidueType, SequenceResidueType> score_and_maybe_align;

  score_and_maybe_align.score_and_maybe_align(
    profile_filename,
    fasta_filename,
    sequence_count,
    true // use viterbi & make alignments
  );

  return 0; // success
} // main (..)

