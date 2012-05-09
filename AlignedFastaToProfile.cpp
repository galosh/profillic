/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      AlignedFastaToProfile.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The alignedFastaToProfile program.  It takes a multiple alignment and
##      creates a profile from it.
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
#include <seqan/index.h>

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



template <class SourceElementType,
          class TargetElementType>
void
ungapSequence (
               seqan::String<SourceElementType> const & source,
               seqan::String<TargetElementType> & target
)
{
  typedef typename seqan::Iterator< seqan::String<SourceElementType > >::Type TSourceIter;

  for( TSourceIter it = seqan::begin( source ); it != seqan::end( source ); ++it ) {
    if( seqan::value( it ) != '-' ) {
      target += seqan::value( it );
    }
  }

  // That's it.
  return;
} // ungapSequence( source const &, target & )

template <class SourceElementType,
          class TargetElementType>
void
ungapFasta (
            galosh::Fasta<SourceElementType> const & source,
            galosh::Fasta<TargetElementType> & target
)
{
  for( int seq_i = 0; seq_i < source.size(); ++seq_i ) {
    target.m_descriptions[ seq_i ] = source.m_descriptions[ seq_i ];
    ungapSequence( source[ seq_i ], target[ seq_i ] );
  }

  // That's it.
  return;
} // ungapFasta( source const &, target & )

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
gappedFastaAndConsensusToProfile (
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters const & parameters,
  galosh::Fasta<char> const & gapped_fasta,
  seqan::String<char> const & gapped_consensus,
 ProfileType & profile
)
{
  galosh::Sequence<SequenceResidueType> ungapped_consensus;
  galosh::ungapSequence( gapped_consensus, ungapped_consensus );
  // TODO: REMOVE
  //std::cout << "gapped consensus: " << gapped_consensus << std::endl;
  // TODO: REMOVE
  //std::cout << "ungapped consensus: " << ungapped_consensus << std::endl;

  uint32_t profile_length = ungapped_consensus.length();

  // Resize it, and reinitialize while we're at it.  Note that this will also
  // even() it.
  profile.reinitialize( profile_length );

  // First calculate the default indel values.
  // For now we just use these -- TODO: calculate these values from the input multiple alignmnt, too.
  setTransitionsFromParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType, ProfileType>(
    parameters,
    profile
  );
  // Note Insertion distribution(s) are even now.  TODO: allow an option to set
  // the Insertion distribution to something else.
  typedef typename seqan::Iterator< seqan::String<char > >::Type TGappedConsensusIter;
  const uint32_t num_seqs = gapped_fasta.size();

  uint32_t prof_pos_i = 0; // Pos in profile.
  uint32_t col_i = 0; // Pos in gapped_fasta seqs
  uint32_t gaps_in_col_i;
  for( TGappedConsensusIter it = seqan::begin( gapped_consensus ); it != seqan::end( gapped_consensus ); ++it, ++col_i ) {
    if( seqan::value( it ) == '-' ) {
      continue;
    }
    gaps_in_col_i = 0;
    profile[ prof_pos_i ][ Emission::Match ].zero();
    for( uint32_t seq_i = 0; seq_i < num_seqs; seq_i++ ) {
      if( gapped_fasta[ seq_i ][ col_i ] == '-' ) {
        gaps_in_col_i += 1;
      } else {
        profile[ prof_pos_i ][ Emission::Match ][ ResidueType( gapped_fasta[ seq_i ][ col_i ] ) ] += 1;
      }
    }
    // TODO: REMOVE
    //std::cout << "Prof pos " << prof_pos_i << ": totals are " << profile[ prof_pos_i ][ Emission::Match ] << "; there were " << gaps_in_col_i << " gaps." << std::endl;
    if( gaps_in_col_i == num_seqs ) {
      profile[ prof_pos_i ][ Emission::Match ].even();
    } else {
      profile[ prof_pos_i ][ Emission::Match ] /=
        ( double( num_seqs ) - gaps_in_col_i );
    }
    ++prof_pos_i;
  } // End foreach gapped_consensus char..

  // That's it.
  return;
} // gappedFastaAndConsensusToProfile( ProlificParameters::Parameters const &, Fasta<char> const &, String<char> const &, ProfileType & )

} // End namespace galosh

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
    cout << "Usage: " << argv[ 0 ] << " <input (Aligned Fasta) filename> <consensus (Aligned Fasta) filename> [<output (galosh Profile) filename>]" << endl;
    cout << "Note that only the first sequence in the consensus file will be used, and only its gaps matter: it must be the same length as the aligned strings in the first argument, and its gaps will determine which columns of the alignment are to be used to construct the profile (ie all non-gap positions of the 'consensus')." << endl;
    exit( 1 );
  }
  const bool be_verbose = true; //( argc > 4 );
  if( be_verbose ) {
    cout << "Reading multiple alignment from Aligned Fasta file '" << argv[ 1 ] << "'" << endl;
  }
  galosh::Fasta<char> gapped_fasta;
  gapped_fasta.fromFile( argv[ 1 ] );
  if( gapped_fasta.size() == 0 ) {
    cout << "No sequences were found in the Fasta file '" << argv[ 1 ] << "'" << endl;
    return 1;
  }
  if( be_verbose ) {
    cout << "Reading consensus from Aligned Fasta file '" << argv[ 2 ] << "'" << endl;
  }
  galosh::Fasta<char> gapped_consensus_fasta;
  gapped_consensus_fasta.fromFile( argv[ 2 ] );
  if( gapped_consensus_fasta.size() == 0 ) {
    cout << "No sequences were found in the consensus Fasta file '" << argv[ 2 ] << "'" << endl;
    return 1;
  } else if( gapped_consensus_fasta.size() > 1 ) {
    if( be_verbose ) {
      cout << "WARNING: Using only the first sequence in the given consensus Fasta file." << endl;
    }
  }
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << gapped_consensus_fasta[ 0 ];
    cout << endl;
  }

  // The parameters
  // TODO: Let user choose some from command line.
  galosh::ProlificParameters<ResidueType, floatrealspace, floatrealspace, floatrealspace>::Parameters parameters;

  // Create ungapped consensus
  seqan::String<char> gapped_consensus( *( dynamic_cast<const seqan::String<char> * const>( & gapped_consensus_fasta[ 0 ] ) ) );
  galosh::ProfileTreeRoot<ResidueType, floatrealspace> profile;
  galosh::gappedFastaAndConsensusToProfile<ResidueType, floatrealspace, floatrealspace, floatrealspace, SequenceResidueType, galosh::ProfileTreeRoot<ResidueType, floatrealspace> >( parameters, gapped_fasta, gapped_consensus, profile );

  // TODO: walk through consensus, and at non-gap positions, calculate frequencies.  Also can calculate gap open and extend probs..

 if( argc > 3 ) {
   if( be_verbose ) {
     cout << "Writing Profile to file '" << argv[ 3 ] << "'" << endl;
   }
   std::ofstream profile_stream( argv[ 3 ] );
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
