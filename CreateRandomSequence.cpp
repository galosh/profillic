/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      CreateRandomSequence.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The createRandomSequence program.  It creates a sequence of a specified
##      length by drawing from a discrete distribution over the DNA alphabet.
##
##      NOTE: At present this supports only DNA, not AminoAcids.
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

#include "MultinomialDistribution.hpp"
#include "Fasta.hpp"
#include "Sequence.hpp"
#include "Random.hpp"
#include "Algebra.hpp"

#include <iostream>

#include <boost/lexical_cast.hpp>

#include <seqan/basic.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

namespace galosh {

// For now we're assuming Dna -- TODO: generalize
template <class ProbabilityType>
void
createRandomSequence (
  const uint32_t & length,
  const ProbabilityType & a_prob,
  const ProbabilityType & c_prob,
  const ProbabilityType & g_prob,
  const ProbabilityType & t_prob,
  Random & random,
  Sequence<Dna> & sequence
)
{
  MultinomialDistribution<Dna,ProbabilityType> residue_dist;

  Dna residue;
  residue_dist[ ( residue = 'a' ) ] = a_prob;
  residue_dist[ ( residue = 'c' ) ] = c_prob;
  residue_dist[ ( residue = 'g' ) ] = g_prob;
  residue_dist[ ( residue = 't' ) ] = t_prob;
  residue_dist.normalize(); // in case they don't add to 1..

  sequence.reinitialize( length );
  for( uint32_t i = 0; i < length; i++ ) {
    sequence[ i ] = residue_dist.draw( random );
  }

  return;
} // createRandomSequence(..)

// For now we're assuming Dna -- TODO: generalize
template <class ProbabilityType>
void
createRandomSequence (
  const uint32_t & length,
  const ProbabilityType & a_prob,
  const ProbabilityType & c_prob,
  const ProbabilityType & g_prob,
  const ProbabilityType & t_prob,
  uint32_t const & random_seed,
  Sequence<Dna> & sequence
)
{
  Random random( random_seed );
  createRandomSequence(
    length,
    a_prob,
    c_prob,
    g_prob,
    t_prob,
    random,
    sequence
  );
} // createRandomSequence(..)

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
  // For now we assume a Dna distribution.  TODO: Generalize.
  if( argc < 6 ) {
    cout << "Usage: " << argv[ 0 ] << " <length> <a_prob> <c_prob> <g_prob> <t_prob> [<output (Fasta) filename> [<sequence name> [<random seed>]]]" << endl;
    exit( 1 );
  }
  const bool be_verbose = ( argc >= 7 );

  uint32_t length;
  float a_prob, c_prob, g_prob, t_prob;
  try {
    length = boost::lexical_cast<uint32_t>( argv[ 1 ] );
  } catch( boost::bad_lexical_cast & ) {
    std::cerr << "Unable to interpret the argument '" << argv[ 1 ] << "' as an unsigned long value for use as the desired sequence length." << std::endl;
    exit( 1 );
  } // End try .. catch block for lexical_cast
  try {
    a_prob = boost::lexical_cast<float>( argv[ 2 ] );
  } catch( boost::bad_lexical_cast & ) {
    std::cerr << "Unable to interpret the argument '" << argv[ 2 ] << "' as a real value for use as the probability of the 'A' residue." << std::endl;
    exit( 1 );
  } // End try .. catch block for lexical_cast
  if( a_prob < 0 ) {
    std::cerr << "The given value for the probability of the 'A' residue, " << a_prob << ", is negative.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  if( a_prob > 1 ) {
    std::cerr << "The given value for the probability of the 'A' residue, " << a_prob << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  try {
    c_prob = boost::lexical_cast<float>( argv[ 3 ] );
  } catch( boost::bad_lexical_cast & ) {
    std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as a real value for use as the probability of the 'C' residue." << std::endl;
    exit( 1 );
  } // End try .. catch block for lexical_cast
  if( c_prob < 0 ) {
    std::cerr << "The given value for the probability of the 'C' residue, " << c_prob << ", is negative.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  if( c_prob > 1 ) {
    std::cerr << "The given value for the probability of the 'C' residue, " << c_prob << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  try {
    g_prob = boost::lexical_cast<float>( argv[ 4 ] );
  } catch( boost::bad_lexical_cast & ) {
    std::cerr << "Unable to interpret the argument '" << argv[ 4 ] << "' as a real value for use as the probability of the 'G' residue." << std::endl;
    exit( 1 );
  } // End try .. catch block for lexical_cast
  if( g_prob < 0 ) {
    std::cerr << "The given value for the probability of the 'G' residue, " << g_prob << ", is negative.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  if( g_prob > 1 ) {
    std::cerr << "The given value for the probability of the 'G' residue, " << g_prob << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  try {
    t_prob = boost::lexical_cast<float>( argv[ 5 ] );
  } catch( boost::bad_lexical_cast & ) {
    std::cerr << "Unable to interpret the argument '" << argv[ 5 ] << "' as a real value for use as the probability of the 'T' residue." << std::endl;
    exit( 1 );
  } // End try .. catch block for lexical_cast
  if( t_prob < 0 ) {
    std::cerr << "The given value for the probability of the 'T' residue, " << t_prob << ", is negative.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  if( t_prob > 1 ) {
    std::cerr << "The given value for the probability of the 'T' residue, " << t_prob << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
    exit( 1 );
  }
  string sequence_name = "Randomly generated sequence";
  if( argc >= 8 ) {
    sequence_name = argv[ 7 ];
  }
  uint32_t random_seed = static_cast<uint32_t>( std::time( NULL ) );
  if( argc >= 9 ) {
    try {
      random_seed = boost::lexical_cast<uint32_t>( argv[ 8 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 8 ] << "' as an unsigned long value for use as the random seed." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
  } // End if argc >= 9

  galosh::Fasta<Dna> fasta( 1 );
  fasta.m_descriptions[ 0 ] = sequence_name;
  galosh::createRandomSequence(
    length,
    ( floatrealspace )a_prob,
    ( floatrealspace )c_prob,
    ( floatrealspace )g_prob,
    ( floatrealspace )t_prob,
    random_seed,
    fasta[ 0 ]
  );

  if( argc >= 7 ) {
    if( be_verbose ) {
      cout << "Writing sequence to file '" << argv[ 6 ] << "'." << endl;
    }
    std::ofstream fasta_stream( argv[ 6 ] );
    assert( fasta_stream.good() );
    fasta_stream << fasta;
    fasta_stream.close();
  } else {
    if( be_verbose ) {
      cout << "Writing sequence to stdout." << endl;
    }
    cout << endl;
    cout << fasta;
  }

  exit( 0 );
} // main (..)
