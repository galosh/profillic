/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      ProfileToConsensus.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The profileToConsensus program.  It creates a consensus sequence representing
##      a Profile HMM (found in a profile file).
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
#include "Fasta.hpp"

#include <iostream>

#include <seqan/basic.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace seqan;

namespace galosh {

/////////////
/**
 * Given a pointer to a HMMer profile (created using hmmer::AllocPlan7Shell()),
 * allocate its body and populate it with data from the given galosh Profile.
 * Later, you should free the hmm with FreePlan7( &hmm ).
 */
template <class ProfileType, class ResidueType>
void
profileToConsensus (
  ProfileType const & profile,
  Sequence<ResidueType> & sequence
)
{
  uint32_t profile_length = profile.length();
  sequence.reinitialize( profile_length );
  for( uint32_t pos_i = 0; pos_i < profile_length; pos_i++ ) {
    sequence[ pos_i ] =
      profile[ pos_i ][ Emission::Match ].maximumValueType();
  }
  // That's it.
  return;
} // profileToConsensus( ProfileType const &, Sequence & )

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..

  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (galosh Profile) filename> [<output (Fasta) filename>]" << endl;
    exit( 1 );
  }
  const bool be_verbose = ( argc >= 3 );
  if( be_verbose ) {
    cout << "Reading profile from file '" << argv[ 1 ] << "'" << endl;
  }
  galosh::ProfileTreeRoot<ResidueType, floatrealspace> profile;
  profile.fromFile( argv[ 1 ] );
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << profile;
    cout << endl;
  }

  galosh::Fasta<ResidueType> fasta( 1 );
  fasta.m_descriptions[ 0 ] = "Consensus"; // TODO: MAGIC #
  profileToConsensus( profile, fasta[ 0 ] );

  if( argc >= 3 ) {
    if( be_verbose ) {
      cout << "Writing Fasta to file '" << argv[ 2 ] << "'" << endl;
    }
    std::ofstream fasta_stream( argv[ 2 ] );
    assert( fasta_stream.good() );
    fasta_stream << fasta;
    fasta_stream.close();
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }
  } else {
    if( be_verbose ) {
      cout << "Fasta is:" << endl;
    }
    cout << fasta;
    cout << endl;
  }

  exit( 0 );
} // main (..)
