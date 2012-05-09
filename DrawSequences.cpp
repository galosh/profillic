/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      DrawSequences.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The drawSequences program.  It draws sequences from a distribution
##      represented in a Profile HMM (found in a profile file).
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
#*    profuse is free software: you can redistribute it and/or modify it under
#*    the terms of the GNU Lesser Public License as published by the Free
#*    Software Foundation, either version 3 of the License, or (at your option)
#*    any later version.
#*
#*    profuse is distributed in the hope that it will be useful, but WITHOUT
#*    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#*    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for
#*    more details.
#*
#*    You should have received a copy of the GNU Lesser Public License along
#*    with profuse.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************/

#include "Algebra.hpp"
#include "Profile.hpp"
#include "DynamicProgramming.hpp"
#include "Fasta.hpp"
#include "Random.hpp"

#include <iostream>

#include <seqan/basic.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

int
main ( int const argc, char const ** argv )
{
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
  typedef seqan::AminoAcid SequenceResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
  typedef seqan::Iupac SequenceResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  if( argc < 3 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (galosh Profile) filename> <num seqs> [<output (Fasta) filename> [<random seed>]]" << endl;
    exit( 1 );
  }

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

  const bool be_verbose = ( argc >= 6 );
  if( be_verbose ) {
    cout << "Reading profile from file '" << argv[ 1 ] << "'" << endl;
  }
  typedef galosh::ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
  ProfileType profile;
  profile.fromFile( argv[ 1 ] );
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << profile;
    cout << endl;
  }

  uint32_t num_draws = 0;
  try {
    num_draws = boost::lexical_cast<uint32_t>( argv[ 2 ] );
  } catch( boost::bad_lexical_cast & ) {
    std::cerr << "Unable to interpret the argument '" << argv[ 2 ] << "' as an unsigned long value for use as the number of sequences to generate." << std::endl;
    exit( 1 );
  } // End try .. catch block for lexical_cast
  
  uint32_t random_seed = static_cast<uint32_t>( std::time( NULL ) );
  if( argc >= 5 ) {
    try {
      random_seed = boost::lexical_cast<uint32_t>( argv[ 4 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 4 ] << "' as an unsigned long value for use as the random seed." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // End if argc >= 4
  galosh::Random random( random_seed );

  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters parameters;

  galosh::Fasta<ResidueType> random_seqs_fasta( num_draws );
  galosh::DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::MultipleAlignment<ProfileType, ResidueType> random_seqs_true_multiple_alignment;

  std::string name_prefix( "Randomly generated sequence from profile \"" );
  name_prefix += argv[ 1 ];
  name_prefix += "\" #";
  dp.drawSequences(
                   parameters,
                   profile,
                   num_draws,
                   name_prefix,
                   random,
                   random_seqs_fasta,
                   random_seqs_true_multiple_alignment
                   );
  if( be_verbose ) { // TODO: ? 
    cout << "Alignment paths of the training sequences are:" << endl;
    random_seqs_true_multiple_alignment.toPairwiseStream( cout );
  } // End if be_verbose

  if( argc >= 4 ) {
    if( be_verbose ) {
      cout << "Writing Fasta to file '" << argv[ 3 ] << "'" << endl;
    }
    std::ofstream fasta_stream( argv[ 3 ] );
    assert( fasta_stream.good() );
    fasta_stream << random_seqs_fasta;
    fasta_stream.close();
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }
  } else {
    if( be_verbose ) {
      cout << "Fasta is:" << endl;
    }
    cout << random_seqs_fasta;
    cout << endl;
  }

  exit( 0 );
} // main (..)
