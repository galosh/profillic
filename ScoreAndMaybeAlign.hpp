/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      ScoreAndMaybeAlign.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the ScoreAndMaybeAlign class, which contains the
##      single method "score_and_maybe_align", implementing both the intended
##      behavior of "Score.cpp" and of "Align.cpp".  See those files for more
##      on the intended behavior.
##
#******************************************************************************
#*
#*    This file is part of profuse, a suite of programs for working with
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2008, 2009, 2011 by Paul T. Edlefsen, Fred Hutchinson
#*    Cancer Research Center.
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

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_SCOREANDMAYBEALIGN_HPP__
#define __GALOSH_SCOREANDMAYBEALIGN_HPP__

#include <Algebra.hpp>

#include "Ambiguous.hpp"
#include "Sequence.hpp"
#include "MultinomialDistribution.hpp"
#include "ProfileHMM.hpp"
#include "Profile.hpp"
#include "Fasta.hpp"
#include "Random.hpp"
#include "DynamicProgramming.hpp"

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find_motif.h>

#ifdef __HAVE_MUSCLE
#include "muscle/distfunc.h"
#include "muscle/clustsetdf.h"
#include "muscle/clust.h"
#include "muscle/tree.h"
#include "muscle/textfile.h"
#endif // __HAVE_MUSCLE

#include <boost/lexical_cast.hpp>

namespace galosh {
template <typename ProbabilityType,
          typename ScoreType,
          typename MatrixValueType,
          typename ResidueType,
          typename SequenceResidueType>
class ScoreAndMaybeAlign {
public:
  // read in a profile and some sequences, calculate forward score and return
  // it, or calculate viterbi score and viterbi alignment.  Returns the score.
  ScoreType
  score_and_maybe_align (
    string const & profile_filename,
    string const & fasta_filename,
    uint32_t sequence_count, // 0 to use the number of sequences in the fasta file
    bool const & use_viterbi // if false, don't align; just calc the forward score.
  ) const
  {
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;

    const bool be_verbose = false;
    const bool be_verbose_show_profiles = false;
    const bool be_verbose_show_sequences = false;

    ProfileType profile;
    if( be_verbose ) {
      cout << "Reading profile from file '" << profile_filename << "'" << endl;
    }
    profile.fromFile( profile_filename );
    if( be_verbose ) {
      if( be_verbose_show_profiles ) {
        cout << "\tgot:" << endl;
        cout << profile;
        cout << endl;
      } else {
        cout << "\tdone." << endl;
      }
    } // End if be_verbose

    Fasta<SequenceResidueType> fasta;
    if( be_verbose ) {
      cout << "Reading sequences from Fasta file '" << fasta_filename << "'" << endl;
    }
    fasta.fromFile( fasta_filename );
    if( be_verbose ) {
      if( be_verbose_show_sequences ) {
        cout << "\tgot:" << endl;
        cout << fasta;
        cout << endl;
      } else {
        cout << "\tdone." << endl;
      }
    } // End if be_verbose

    sequence_count = ( ( sequence_count == 0 ) ? fasta.size() : min( static_cast<size_t>( sequence_count ), fasta.size() ) );

    if( be_verbose ) {
      cout << "Allocating the dp matrices." << endl;
    }
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer dp_matrices(
      profile,
      fasta,
      sequence_count
    );
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }

    ScoreType score;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters parameters;

    if( use_viterbi ) {
      if( be_verbose ) {
        cout << "Calculating the viterbi score, and computing the dp matrices for the multiple alignment." << endl;
      }
      score =
        dp.forward_score_viterbi(
          parameters,
          profile,
          fasta,
          sequence_count,
          dp_matrices
        );
      if( be_verbose ) {
        cout << "\tThe total viterbi score for these sequences is: " << score << endl;
      }
    } else { // if use_viterbi .. else ..
      if( be_verbose ) {
        cout << "Calculating the forward score." << endl;
      }
      score =
        dp.forward_score(
          parameters,
          profile,
          fasta,
          sequence_count,
          dp_matrices
        );
      if( be_verbose ) {
        cout << "\tThe total probability of these sequences, given this profile model, is: " << score << endl;
      }
      return score; // Can't align unless we make viterbi matrices.
    } // End if use_viterbi .. else ..
    // End calculating viterbi score and filling the dp matrices
  
    if( be_verbose ) {
      cout << "Backtracing to compute the alignments." << endl;
    }
    // Show multiple alignment
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> ma(
      &profile,
      &fasta,
      sequence_count
    );
    dp.forward_viterbiAlign(
      parameters,
      dp_matrices,
      ma
    );
    if( be_verbose ) {
      cout << "\tThe multiple alignment is:" << endl;
    }
    // TODO: Make this a parameter to choose between output formats.
    //ma.toPileupStream( cout, &fasta.m_descriptions );
    ma.toPairwiseStream( cout, &fasta.m_descriptions );
    //ma.toAlignedFastaStream( cout, &fasta.m_descriptions );

    return score;
  } // score_and_maybe_align ( string const &, string const &, bool const & use_viterbi )

}; // End class ScoreAndMaybeAlign

} // End namespace galosh

#endif // __GALOSH_SCOREANDMAYBEALIGN_HPP__
