/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      SequenceToProfile.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The sequenceToProfile program.  It takes a sequence (in unaligned Fasta
##      format) and creates a DNA Profile HMM from it with the given
##      conservation rate (default: .75).  The profile will have the same
##      length as the input sequence, with each position having some mass (the
##      conservation rate) on the residue corresponding to the input sequence
##      at the same position, and all remaining mass evenly divided among the
##      remaining residues.  The transition parameters are set to the same
##      defaults as are used in other profuse programs (see
##      ProlificParameters.hpp).
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
#include "ProlificParameters.hpp" // for the parameters

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
 */
template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType,
          class ProfileType>
void
setTransitionsFromParameters (
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters const & parameters,
  ProfileType & profile
)
{
  uint32_t profile_length = profile.length();

  // Use the first element of the expectedDeletionsCounts parameter vector, if
  // there is one, or 1.0.
  double expected_deletions_count =
    (
      parameters.expectedDeletionsCounts ?
      ( *parameters.expectedDeletionsCounts )[ 0 ] :
      1.0
    );
  // If useDeletionsForInsertionsParameters, use the expected_deletions_count;
  // otherwise, use the first element of the expectedInsertionsCounts parameter
  // vector, if there is one, or 1.0.
  double expected_insertions_count =
    (
      parameters.useDeletionsForInsertionsParameters ?
      expected_deletions_count :
      (
        ( parameters.expectedInsertionsCounts ?
        ( *parameters.expectedInsertionsCounts )[ 0 ] :
        1.0
        )
      )
    );
  double expected_deletion_length_as_profile_length_fraction =
    (
      parameters.expectedDeletionLengthAsProfileLengthFractions ?
      ( *parameters.expectedDeletionLengthAsProfileLengthFractions )[ 0 ] :
      1.0
    );
  double expected_insertion_length_as_profile_length_fraction =
    (
      parameters.useDeletionsForInsertionsParameters ?
      expected_deletion_length_as_profile_length_fraction :
      (
        parameters.expectedInsertionLengthAsProfileLengthFractions ?
        ( *parameters.expectedInsertionLengthAsProfileLengthFractions )[ 0 ] :
        1.0
      )
    );

  ProbabilityType deletion_open =
    ( expected_deletions_count / profile_length );
  ProbabilityType insertion_open =
    ( parameters.useDeletionsForInsertionsParameters ?
    deletion_open :
    ( expected_insertions_count / profile_length ) );
        
  // [ the EV of a geometric is 1/p, where p is prob of stopping, so if q is the prob of continuing, we want ( 1 - q ) = 1/EV. ]
  ProbabilityType deletion_extension =
    ( 1.0 - min( ( 1.0 / ( expected_deletion_length_as_profile_length_fraction * profile_length ) ), ( 1.0 / parameters.minExpectedDeletionLength ) ) );
  ProbabilityType insertion_extension =
    ( parameters.useDeletionsForInsertionsParameters ? deletion_extension : ( 1.0 - min( ( 1.0 / ( expected_insertion_length_as_profile_length_fraction * profile_length ) ), ( 1.0 / parameters.minExpectedInsertionLength ) ) ) );
                    
  // Now set up the profile(s)
  profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] =
    ( parameters.preAlignInsertion );
  profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] =
    ( 1 ) -
    profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
  profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] =
    deletion_open;
  profile[ Transition::fromBegin ][ TransitionFromBegin::toMatch ] =
    ( 1 ) -
    profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
                        
  profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
    insertion_open;
  profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
    deletion_open;
  profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
    ( 1.0 ) -
    (
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] +
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
    );
  profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
    insertion_extension;
  profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
    ( 1.0 ) -
    profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
  profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
    deletion_extension;
  profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
    ( 1.0 ) -
    profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
                    
  // For now we don't use the End distribution ..
  //profile[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ] = ( 1 );
  //profile[ Transition::fromEnd ][ TransitionFromEnd::toLoop ] = ( 0 );
                      
  profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] =
    ( parameters.postAlignInsertion );
  profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] =
    ( 1.0 ) -
    profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];

} // setTransitionsFromParameters( ProlificParameters::Parameters const &, ProfileType & profile )


/////////////
/**
 */
template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType,
          class SequenceResidueType,
          class ProfileType>
void
consensusToProfile (
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters const & parameters,
  Sequence<SequenceResidueType> const & sequence,
  ProfileType & profile,
  double const & conservation_rate
)
{
  uint32_t profile_length = sequence.length();

  // Resize it, and reinitialize while we're at it.  Note that this will also
  // even() it.
  profile.reinitialize( profile_length );

  // First calculate the appropriate indel values.
  setTransitionsFromParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType, ProfileType>(
    parameters,
    profile
  );
  // Note Insertion distribution(s) are even now.  TODO: allow an option to set
  // the Insertion distribution to something else.

  // This is a trick to get the values set correctly:
  // make an even profile, then set one of them higher than you need it, but
  // normalize to get the right thing.

  // The particular position we use here is arbitrary, since the profile starts
  // out even().
  ProbabilityType pattern_trick_value =
    ( ( conservation_rate == 1.0 ) ? ( 1.0 ) :
    ( ( ( 1.0 ) - profile[ 0 ][ Emission::Match ][ sequence[ 0 ] ] ) *
    ( ( conservation_rate / ( 1.0 - conservation_rate ) ) ) ) );
                      
  // r is current remaining value (1 - P(base)), p is target value.
  //x/(r + x) = p
  //p( r + x) = x
  // rp + px = x
  // x - px = rp
  // x( 1 - p ) = rp
  // x = r( p / ( 1 - p ) )
  //cout << "Pattern trick value " << pattern_trick_value << endl; 
                      
  for( uint32_t pos_i = 0; pos_i < profile_length; pos_i++ ) {
    if( conservation_rate == 1.0 ) {
      profile[ pos_i ][ Emission::Match ].zero();
    }
    profile[ pos_i ][ Emission::Match ][ sequence[ pos_i ] ] = pattern_trick_value;
    if( conservation_rate != 1.0 ) {
      profile[ pos_i ][ Emission::Match ].normalize( 0 );
    }
  } // End foreach position, set it up according to the pattern and conservation_rate.

  // That's it.
  return;
} // consensusToProfile( ProlificParameters::Parameters const &, Sequence const &, ProfileType &, double const & )

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
  // For now we assume a Dna distribution.  TODO: Generalize.
  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (Fasta) filename> [<output (galosh Profile) filename> [<conservation rate>]]" << endl;
    exit( 1 );
  }
  const bool be_verbose = true; //( argc >= 3 );
  if( be_verbose ) {
    cout << "Reading Fasta from file '" << argv[ 1 ] << "'" << endl;
  }
  galosh::Fasta<seqan::Dna> fasta;
  fasta.fromFile( argv[ 1 ] );
  if( fasta.size() == 0 ) {
    cout << "No sequences were found in the Fasta file '" << argv[ 1 ] << "'" << endl;
    return 1;
  } else if( fasta.size() > 1 ) {
    if( be_verbose ) {
      cout << "WARNING: Using only the first sequence in the given Fasta file." << endl;
    }
  }
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << fasta[ 0 ];
    cout << endl;
  }

  // The parameters
  // TODO: Let user choose some from command line.
  galosh::ProlificParameters<seqan::Dna, floatrealspace, floatrealspace, floatrealspace>::Parameters parameters;

  // arg 3 is the conservation rate (default .75)
  double conservation_rate = .75; // TODO: DEHACKIFY MAGIC # DEFAULT conservation_rate !!
  if( argc >= 3 ) {
    try {
      conservation_rate = boost::lexical_cast<double>( argv[ 3 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as a real value for use as the conservation rate." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( conservation_rate <= 0 ) {
      std::cerr << "The given conservation rate value, " << conservation_rate << ", is zero or negative.  You must supply a value between 0 and 1, and not 0." << std::endl;
      exit( 1 );
    }
    if( conservation_rate > 1 ) {
      std::cerr << "The given conservation rate value, " << conservation_rate << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
      exit( 1 );
    }
  } // End if argc >= 3

  galosh::ProfileTreeRoot<seqan::Dna, floatrealspace> profile;
  galosh::consensusToProfile<seqan::Dna, floatrealspace, floatrealspace, floatrealspace, seqan::Dna, galosh::ProfileTreeRoot<seqan::Dna, floatrealspace> >( parameters, fasta[ 0 ], profile, conservation_rate );

  if( argc >= 3 ) {
    if( be_verbose ) {
      cout << "Writing Profile to file '" << argv[ 2 ] << "'" << endl;
    }
    std::ofstream profile_stream( argv[ 2 ] );
    assert( profile_stream.good() );
    profile_stream << profile;
    profile_stream.close();
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }
  } else {
    if( be_verbose ) {
      cout << "Profile is:" << endl;
    }
    cout << profile;
    cout << endl;
  }

  exit( 0 );
} // main (..)
