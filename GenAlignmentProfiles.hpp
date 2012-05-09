/**-------------------------------------------------------------------------##
 *  Library:
 *     galosh::profuse
 *  @file
 *      GenAlignmentProfiles.hpp
 *  @author
 *      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 *      Slightly modified from ScoreAndMaybeAlign class, Ted Holzman 2/1/2012
 *  @description
 *      Class definition for the GenAlignmentProfile class, which contains the
 *      single method "gen_alignment_profiles", implementing both the intended
 *      behavior of "Align.cpp" and calculateAlignmentProfiles in
 *      DynamicProgramming.hpp. See those files for more
 *      on the intended behavior.
 *
 *******************************************************************************
 *    @license
 *    This file is part of profuse, a suite of programs for working with
 *    Profile HMMs.  Please see the document CITING, which should have been
 *    included with this file.  You may use at will, subject to the license
 *    (Apache v2.0), but *please cite the relevant papers* in your documentation
 *    and publications associated with uses of this library.  Thank you!
 *
 *    Copyright (C) 2008, 2009, 2011 by Paul T. Edlefsen, Fred Hutchinson
 *    Cancer Research Center.
 *
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
 **//*******************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_GENALIGNMENTPROFILES_HPP__
#define __GALOSH_GENALIGNMENTPROFILES_HPP__

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
#include "stddef.h"

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
#include <boost/program_options.hpp>
namespace galosh {
/**
 * \class GenAlignmentProfile
 * \brief contains routines to create and manipulate Alignment Profiles
 *
 */
template <typename ProbabilityType,
          typename ScoreType,
          typename MatrixValueType,
          typename ResidueType,
          typename SequenceResidueType>
class GenAlignmentProfiles {
public:
  /**
   * \fn gen_alignment_profiles
   * \brief read in a profile and some sequences, generate one or more
   * Alignment Profiles.
   **/
  std::vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile>
  gen_alignment_profiles ( boost::program_options::variables_map vm ) const
  {
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
    /**
     * obtain program parameters from variables_map
     */
    const std::string profile_filename = vm["profile"].as<string>();
    const std::string fasta_filename = vm["fasta"].as<string>();
    int sequence_count = 0;
    if( vm.count( "nseq" ) ) {
      sequence_count = vm["nseq"].as<int>();
    }
    int verbosity = 0;
    if( vm.count( "verbosity" ) ) {
      verbosity = vm["verbosity"].as<int>();
    }
    const bool be_verbose = verbosity > 0;
    const bool be_verbose_show_profiles = verbosity > 1;
    const bool be_verbose_show_sequences = verbosity > 2;
    const bool use_viterbi = vm.count( "viterbi" ) > 0;
    const bool indiv_profiles = vm.count( "individual" ) > 0;
    std::string alignment_format = "fasta";
    if( vm.count( "format" ) ) {
      alignment_format = vm["format"].as<std::string>();
    }

    ProfileType profile;
    if( be_verbose ) {
      cerr << "Reading profile from file '" << profile_filename << "'" << endl;
    }
    if( !profile.fromFile( profile_filename ) )
    {
      throw( "Can't open profile file " + profile_filename );
    }
    if( be_verbose ) {
      if( be_verbose_show_profiles ) {
        cerr << "\tgot:" << endl;
        cerr << profile;
        cerr << endl;
      } else {
        cerr << "\tdone." << endl;
      }
    } // End if be_verbose

    Fasta<SequenceResidueType> fasta;
    if( be_verbose ) {
      cerr << "Reading sequences from Fasta file '" << fasta_filename << "'" << endl;
    }
    /// \todo Find out why Fasta.fromFile(string) returns void instead of boolean
    if( !fasta.fromFile( fasta_filename.c_str() ) ) {
      throw "Can't open fasta file " + fasta_filename;
    }
    if( be_verbose ) {
      if( be_verbose_show_sequences ) {
        cerr << "\tgot:" << endl;
        cerr << fasta;
        cerr << endl;
      } else {
        cerr << "\tdone." << endl;
      }
    } // End if be_verbose

    sequence_count = ( ( sequence_count == 0 ) ? fasta.size() : min( static_cast<size_t>( sequence_count ), fasta.size() ) );

    if( be_verbose ) {
      cerr << "Allocating the dp matrices for " << sequence_count << " sequences." << endl;
    }
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer dp_matrices(
      profile,
      fasta,
      sequence_count
    );
    if( be_verbose ) {
      cerr << "\tdone." << endl;
    }

    ScoreType score;
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters parameters;

    if( be_verbose ) {
      cerr << "Computing the dp matrices for the multiple alignment." << endl;
    }
    if( use_viterbi ) {
      score =
        dp.forward_score_viterbi(
          parameters,
          profile,
          fasta,
          sequence_count,
          dp_matrices
        );
      if( be_verbose ) {
        cerr << "\tThe total viterbi score for these sequences is: " << score << endl;
      }

      if( be_verbose_show_sequences ) {
        if( be_verbose ) {
          cerr << "Backtracing to compute the alignments." << endl;
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
        if(be_verbose) cerr << "\tdone." << endl;
        cerr << "\tThe multiple alignment is:" << endl;
        if( alignment_format == "pairwise" ) {
          ma.toPairwiseStream( cerr, &fasta.m_descriptions );
        } else if( alignment_format == "fasta" ) {
          ma.toAlignedFastaStream( cerr, &fasta.m_descriptions );
        } else { //if( alignment_format == "pileup" ) {
          ma.toPileupStream( cerr, &fasta.m_descriptions );
        }
      }

    } else { // if use_viterbi .. else ..
      score =
        dp.forward_score(
          parameters,
          profile,
          fasta,
          sequence_count,
          dp_matrices
        );
      if( be_verbose ) {
        cerr << "\tThe total probability of these sequences, given this profile model, is: " << score << endl;
      }
      //return score; // Can't align unless we make viterbi matrices. ///TAH ? is this true or is it a leftover comment? // DOPTE: You can't make traditional/standard alignments without knowing one "best" alignment per sequence, which is something viterbi() does but not forward().
    } // End if use_viterbi .. else ..
    // End calculating score and filling the dp matrices
  
    std::vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile> alignment_profiles(sequence_count);
    for( int i = 0; i < alignment_profiles.size(); i++ )
    {
      alignment_profiles[ i ].reinitialize( profile.length() + 1 );
    }
    if( be_verbose ) {
      cerr << "Calculating the alignment profiles." << endl;
    }
    dp.calculateAlignmentProfiles(
       parameters,
       profile,
       fasta,
       sequence_count,
       dp_matrices,
       alignment_profiles
    );
    if( be_verbose ) {
      cerr << "\tdone." << endl;
    }

    /// DOPTE: Normalizing each individual profile isn't quite right, because the expected counts for the multi-use parameters (the internal transition parameters and the insertion emission parameters) can be greater than one and their sums vary by sequence.
    //if( be_verbose ) {
    //  cerr << "Normalizing " << alignment_profiles.size() << " alignment profiles" << endl;
    //}
    //for( int i = 0; i < alignment_profiles.size(); i++ )
    //{
    //  alignment_profiles[ i ].normalize( 0.0 );
    //}
    //if( be_verbose ) {
    //  cerr << "\tdone." << endl;
    //}

    if( indiv_profiles ) {
      return alignment_profiles;
    }

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile averageProfile;
    averageProfile.reinitialize( profile.length() + 1 );
    averageProfile.zero();
    if( be_verbose ) {
      cerr << "Averaging " << alignment_profiles.size() << " profiles" << endl;
    }
    for( int i = 0; i < alignment_profiles.size(); i++ )
    {
      averageProfile += alignment_profiles[ i ];
    }
    if( be_verbose ) {
      cerr << "\tdone." << endl;
    }
    averageProfile.normalize( 0.0 );

    alignment_profiles.clear();
    alignment_profiles.push_back( averageProfile );
    return alignment_profiles;

  } // gen_alignment_profiles( variables_map vm )

}; // End class GenAlignmentProfiles

} // End namespace galosh

#endif // __GALOSH_GENALIGNMENTPROFILES_HPP__
