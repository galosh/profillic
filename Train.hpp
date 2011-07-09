/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profillic
##  File:
##      Train.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the Train class, which implements the profillic
##      program.  It takes a some unaligned sequences and estimates parameters
##      of a profile HMM, and then outputs the profile HMM in the galosh
##      profile format.
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

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_TRAIN_HPP__
#define __GALOSH_TRAIN_HPP__

#include "Ambiguous.hpp"
#include "Algebra.hpp"
#include "Sequence.hpp"
#include "MultinomialDistribution.hpp"
#include "ProfileHMM.hpp"
#include "Profile.hpp"
#include "Fasta.hpp"
#include "Random.hpp"
#include "DynamicProgramming.hpp"
#include "ProfileTrainer.hpp"
//#include "ProfileTreeTrainer.hpp"
//#include "ProfileGibbs.hpp"
//#include "ProlificParameters.hpp"

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
class Train {
public:
  ScoreType
  train (
    string const & fasta_filename,
    unsigned int const & initial_profile_length_arg, // 0 means use median length, 1 means use max length
    double const & lengthadjust_insertion_threshold, // 0 means don't use lengthadjust.
    double const & lengthadjust_deletion_threshold,
    double const & lengthadjust_insertion_threshold_increment,
    double const & lengthadjust_deletion_threshold_increment,
    long const & random_seed_arg, // 0 means use the time.
    double const & lengthadjust_occupancy_threshold,
    bool const & lengthadjust_use_sensitive_thresholding,
    uint32_t const & lengthadjust_increase_thresholds_for_length_changes_start_iteration,
    double const & lengthadjust_increase_thresholds_for_length_changes_min_increment,
    bool const & use_unconditional_bw,
    double const & even_starting_profile_multiple, // If < 0, even() the positions of the starting profile; if > 0, starting profile's positions are averaged with this multiple of the even profile.
    bool const & train_globals_first
  )
  {
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
  
    uint32_t initial_profile_length = initial_profile_length_arg;

    // NOTE: See below: we start with positions even, since the results are much better.
    //Random random = Random();
    //unsigned long seed = random_seed_arg;
    //if( random_seed_arg == 0 ) {
    //  seed = std::time( NULL );
    //}
    //cout << "NOTE: RANDOM SEED IS " << seed << endl;
    //random.setSeed( seed );
  
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
  
    Fasta<SequenceResidueType> fasta;
    fasta.fromFile( fasta_filename );
    uint32_t num_sequences_to_use = fasta.size();
  
    //cout << "FASTA from file:" << endl;
    //cout << fasta << endl;

    if( initial_profile_length == 1 ) {
      // Guess the length.  Use the maximum length of the sequences.
      initial_profile_length = 0;
      uint32_t seq_length;
      for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
        seq_length = fasta[ seq_i ].length();
        if( seq_length > initial_profile_length ) {
          initial_profile_length = seq_length;
        }
      } // End foreach seq, find the one with the max length
    } else if( initial_profile_length == 0 ) { // if initial_profile_length == 1 .. else 0 ..
      // Guess the length.  Use the median length of the sequences.
      vector<uint32_t> lengths( num_sequences_to_use );
      for( uint32_t seq_i = 0; seq_i < num_sequences_to_use; seq_i++ ) {
        lengths[ seq_i ] = fasta[ seq_i ].length();
      }
      // Sort them.
      std::sort( lengths.begin(), lengths.end() );
      if( ( num_sequences_to_use % 2 ) == 1 ) {
        // Odd number, use the middle value.
        initial_profile_length = lengths[ ( num_sequences_to_use - 1 ) / 2 ];
      } else {
        // Even number, use the average of the two middle values.
        initial_profile_length = lengths[ ( num_sequences_to_use / 2 ) - 1 ];
        initial_profile_length += lengths[ ( num_sequences_to_use / 2 ) ];
        // Round it, rather than truncate it:
        initial_profile_length =
          uint32_t( ( ( float )initial_profile_length / 2.0f ) + .5f );
      }
    } // End if initial_profile_length == 1 .. else 0 ..
    if( initial_profile_length == 1 ) { // We can't handle length 1 profiles.
      initial_profile_length = 2;
    }
    assert( initial_profile_length > 1 );

    // Set up the profile with some reasonable initial global values
    ProfileType profile( initial_profile_length );

#ifndef DISALLOW_FLANKING_TRANSITIONS
    profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] =
      .01;
    profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] =
      1.0 -
      profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
    profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
      .01;
    profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ] =
      .9;
    profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
      1.0 -
      (
        profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] +
        profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletionIn ]
      );
    profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ] =
      .9999;//.5;
    profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toMatch ] =
      1.0 -
      profile[ galosh::Transition::fromDeletionIn ][ galosh::TransitionFromDeletionIn::toDeletionIn ];
#else
    profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
      .01;
    profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
      1.0 -
      profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..
    
#ifdef USE_DEL_IN_DEL_OUT
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] =
      .0025;
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] =
      .0025;
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ] =
      .9; // .5;
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
      1.0 -
      (
        profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
        profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +
        profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletionOut ]
      );
    profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ] =
      .9999;//.5;
    profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toEnd ] =
      1.0 -
      profile[ galosh::Transition::fromDeletionOut ][ galosh::TransitionFromDeletionOut::toDeletionOut ];
#else
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] =
      .01;
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] =
      .01;
    profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] =
      1.0 -
      (
        profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +
        profile[ 0 ][ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ]
      );
#endif // USE_DEL_IN_DEL_OUT .. else ..
    profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] =
      .5;
    profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] =
      1.0 -
      profile[ 0 ][ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ];
    profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] =
      .5;
    profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] =
      1.0 -
      profile[ 0 ][ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ];
#ifndef DISALLOW_FLANKING_TRANSITIONS
    profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] =
      .01;
    profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] =
      1.0 -
      profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
  
    // Start from uniform profile positions
    // UPDATE: Don't!  starting from even helps a lot.
    //profile.uniformPositions( random );

    // Ensure values aren't too tiny.
    profile.normalize( 1E-5 );

    return
      train(
        fasta,
        profile,
        lengthadjust_insertion_threshold, // 0 means don't use lengthadjust.
        lengthadjust_deletion_threshold,
        lengthadjust_insertion_threshold_increment,
        lengthadjust_deletion_threshold_increment,
        lengthadjust_occupancy_threshold,
        lengthadjust_use_sensitive_thresholding,
        lengthadjust_increase_thresholds_for_length_changes_start_iteration,
        lengthadjust_increase_thresholds_for_length_changes_min_increment,
        use_unconditional_bw,
        even_starting_profile_multiple,
        train_globals_first
      );
  } // train( string const &, int const &, double const &, double const &, double const &, usigned long const &, double const &, bool const &, uint32_t const &, double const &, bool const &, double const &, bool const & )

  ScoreType
  train (
    string const & fasta_filename,
    string const & profile_filename,
    double const & lengthadjust_insertion_threshold, // 0 means don't use lengthadjust.
    double const & lengthadjust_deletion_threshold,
    double const & lengthadjust_insertion_threshold_increment,
    double const & lengthadjust_deletion_threshold_increment,
    double const & lengthadjust_occupancy_threshold,
    bool const & lengthadjust_use_sensitive_thresholding,
    uint32_t const & lengthadjust_increase_thresholds_for_length_changes_start_iteration,
    double const & lengthadjust_increase_thresholds_for_length_changes_min_increment,
    bool const & use_unconditional_bw,
    double const & even_starting_profile_multiple, // If < 0, even() the positions of the starting profile; if > 0, starting profile's positions are averaged with this multiple of the even profile.
    bool const & train_globals_first
  )
  {
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;

    ProfileType profile;
    profile.fromFile( profile_filename );
    assert( profile.length() > 1 );

    // Ensure values aren't too tiny.
    profile.normalize( 1E-5 );

    Fasta<SequenceResidueType> fasta;
    fasta.fromFile( fasta_filename );
    uint32_t num_sequences_to_use = fasta.size();
  
    //cout << "FASTA from file:" << endl;
    //cout << fasta << endl;

    return
      train(
        fasta,
        profile,
        lengthadjust_insertion_threshold, // 0 means don't use lengthadjust.
        lengthadjust_deletion_threshold,
        lengthadjust_insertion_threshold_increment,
        lengthadjust_deletion_threshold_increment,
        lengthadjust_occupancy_threshold,
        lengthadjust_use_sensitive_thresholding,
        lengthadjust_increase_thresholds_for_length_changes_start_iteration,
        lengthadjust_increase_thresholds_for_length_changes_min_increment,
        use_unconditional_bw,
        even_starting_profile_multiple,
        train_globals_first
      );
  } // train( string const &, string const &, double const &, double const &, double const &, bool const &, uint32_t const &, double const &, double const &, bool const &, double const &, bool const & )

  ScoreType
  train (
    Fasta<SequenceResidueType> const & fasta,
    ProfileTreeRoot<ResidueType, ProbabilityType> & profile,
    double const & lengthadjust_insertion_threshold, // 0 means don't use lengthadjust.
    double const & lengthadjust_deletion_threshold,
    double const & lengthadjust_insertion_threshold_increment,
    double const & lengthadjust_deletion_threshold_increment,
    double const & lengthadjust_occupancy_threshold,
    bool const & lengthadjust_use_sensitive_thresholding,
    uint32_t const & lengthadjust_increase_thresholds_for_length_changes_start_iteration,
    double const & lengthadjust_increase_thresholds_for_length_changes_min_increment,
    bool const & use_unconditional_bw,
    double const & even_starting_profile_multiple, // If < 0, even() the positions of the starting profile; if > 0, starting profile's positions are averaged with this multiple of the even profile.
    bool const & train_globals_first
  )
  {
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> ProfileType;
    typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
  
    const bool use_lengthadjust = ( lengthadjust_insertion_threshold > 0 );

    const uint32_t initial_profile_length = profile.length();

    const uint32_t train_iterations = 1;//5;
    // Calculate viterbi score after training (root profile only)?
    const bool calculate_viterbi_score = false;//true;
    // Calculate multiple alignment after training (root profile only)?
    const bool calculate_multiple_alignment = false;//true;

    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> dp;
  
    ScoreType score;
  
    if( use_lengthadjust ) {
      cout << "Initial profile length: " << profile.length() << endl;
      cout << "Lengthadjust insertion threshold: " << lengthadjust_insertion_threshold << endl;
      cout << "Lengthadjust deletion threshold: " << lengthadjust_deletion_threshold << endl;
      cout << "Lengthadjust insertion threshold increment: " << lengthadjust_insertion_threshold_increment << endl;
      cout << "Lengthadjust deletion threshold increment: " << lengthadjust_deletion_threshold_increment << endl;
    } else {
      cout << "Profile length: " << profile.length() << endl;
    }
  
    if( even_starting_profile_multiple < 0 ) {
      // Start from even profile positions
      profile.evenPositions();
    } else if( even_starting_profile_multiple > 0 ) {
      ProfileType even_profile( profile.length() );
      even_profile.evenPositions();
      even_profile.multiplyByPositions( even_starting_profile_multiple );
      profile.incrementByPositions( even_profile );
      profile.normalize( 1E-5 );
    } // End if even_starting_profile_multiple

    uint32_t num_sequences_to_use = fasta.size();
    ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> trainer( &profile, fasta, num_sequences_to_use );
  
    // To match ProfuseTest.cpp params, we set them as below:
    // We need this to get the right (default) params for the priorMtoM, etc.
    typename ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters training_parameters_template;
//    {
//      ProfuseTest<ResidueType, ProbabilityType, ScoreType, MatrixValueType, SequenceResidueType> profuse_test;
//      training_parameters_template =
//        profuse_test.m_parameters;
//    }
        
    // TODO: REMOVE!  Testing..  (This is a good setting for the small numbers of observed sequences in eg the STEP trial when training from the starting profile (pretrained from eg all lanl seqs).  Use the known globals -- this avoids converging to silly all-indel profiles, and also acknowledges that the small data set shouldn't overwhelm what we already know about gaps (but TODO: allow using starting profile as a prior, with adjustable prior strength).  Disallowing priors here is key to see the effects of the data, but of course weakening the prior would also be a possibility, and disabling it can encourage local optima.  Again, starting with the 'known' profile seems key.
    //training_parameters_template.trainProfileGlobals = false; // TODO: REMOVE.
    //training_parameters_template.usePriors = false; //true; // TODO: PUT BACK TRUE
    training_parameters_template.usePriors = false;//true;

    // TODO: Try with the lengthadjust: maxBaumWelchInverseScalar > 0
    training_parameters_template.minIterations = 1; // But see below...
    training_parameters_template.maxIterations = 1000;
    // Having maxPositionCycles_globals > 1 seems ok; takes about the same
    // number of iterations, converges to roughly the same place; takes
    // longer by virtue of having more pos cycles per iteration of course.
    training_parameters_template.maxPositionCycles = 1;
    // Having maxPositionCycles_globals > 1 seems to make convergence way
    // slower when lengthadjust is on.  Length keeps adjusting..
    training_parameters_template.maxPositionCycles_globals = 1;
    training_parameters_template.minBaumWelchInverseScalar = 0; // Straight-up bw.
    training_parameters_template.maxBaumWelchInverseScalar = 0; // Straight-up bw.
    training_parameters_template.minBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    training_parameters_template.maxBaumWelchInverseScalar_globals = 0; // Straight-up bw.
    training_parameters_template.scorePercentChangeMinimum_position_cycle = 1;
    training_parameters_template.scorePercentChangeMinimum_iteration = .01;
    training_parameters_template.alwaysAccept = false;//true;//( use_lengthadjust ? true : false );
    training_parameters_template.alwaysAccept_disallowThreshold_profileDistance_iteration = 1E-5;

    // The default is 1E-5...
    //training_parameters_template.euclideanDistanceMinimum_iteration = 5E-7; // 5E-7 is great, though slow...

    // TODO: REMOVE!  TESTING.
    //training_parameters_template.profileValueMinimum = 0; // REMOVE!
  
    training_parameters_template.proposeProfileLengthChanges = use_lengthadjust;
    training_parameters_template.useAlignmentProfiles = true;
    training_parameters_template.numIterationsBetweenLengthChanges = 0;
    training_parameters_template.proposeDeletingThreshold =
      lengthadjust_deletion_threshold; //.01;//.025;//.1;
    training_parameters_template.proposeDeletingThreshold_increment =
      lengthadjust_deletion_threshold_increment; //.0005;//.00005;//.0005;//5E-5;//.0003125;//.00625;//.025;
    training_parameters_template.proposeInsertingThreshold =
      lengthadjust_insertion_threshold;
    training_parameters_template.proposeInsertingPreAlignThreshold = //.35; //.5;
      training_parameters_template.proposeInsertingThreshold;
    training_parameters_template.proposeInsertingPostAlignThreshold = //.35;//.5;
      training_parameters_template.proposeInsertingThreshold;
    training_parameters_template.proposeInsertingThreshold_increment =
      lengthadjust_insertion_threshold_increment;

    training_parameters_template.proposeInsertingOccupancyThreshold =
      lengthadjust_occupancy_threshold;
    training_parameters_template.useSensitiveThresholding =
      lengthadjust_use_sensitive_thresholding;
    training_parameters_template.increaseThresholdsForLengthChanges_startIteration =
      lengthadjust_increase_thresholds_for_length_changes_start_iteration;
    training_parameters_template.increaseThresholdsForLengthChanges_minIncrement =
      lengthadjust_increase_thresholds_for_length_changes_min_increment;
  
    //if( have_trained_profile && start_with_trained_profile ) {
    //  // When we start with the trained profile, we need to get past the
    //  // length modification wait time (which is just one iteration).
    //  training_parameters_template.minIterations =
    //    max( training_parameters_template.minIterations, ( uint32_t )2 );
    //}
  
    if( train_globals_first ) {
      // When we start with the globals, we might not get an opportunity to
      // really train unless we force at least 2 iterations.
      training_parameters_template.minIterations =
        max( training_parameters_template.minIterations, ( uint32_t )2 );
    }

    // Use rabiner scaling? (default true)
    // NOTE: You must change the MatrixValueType to logspace or bfloat if this is false!
    // NOTE: For now if USE_DEL_IN_DEL_OUT, keep this as false.
    training_parameters_template.useRabinerScaling = false;

    // Train globals first?
    training_parameters_template.trainGlobalsFirst = train_globals_first;

    // Use Ubw?
    training_parameters_template.useUnconditionalBaumWelch =
      use_unconditional_bw;
  
#ifdef ALLOW_BOLTZMANN_GIBBS
    // NOTE about priors:  since globals are presently not updated using Baldi, you can still usePriors and they will affect the globals *but not the positions*.
    training_parameters_template.baldiLearningRate = 1; // 0 means noBaldi!
    training_parameters_template.baldiTemperature = 1;
    training_parameters_template.baldiHybrid = false;
    training_parameters_template.siegelMaxFindingThePeakAttempts_positions = 1000; // 0 means Baldi not Baldi / Siegel !!!
    training_parameters_template.siegelEpsilonScaleFactor = 1.5;
    training_parameters_template.siegelMaxRefiningThePeakSteps_positions = 1;//1000;
    training_parameters_template.siegelRefiningThePeakStepsConvergenceThreshold = 1E-5;
    training_parameters_template.minBaumWelchInverseScalar = 0;
    training_parameters_template.maxBaumWelchInverseScalar = 0;
    //training_parameters_template.maxPositionCycles = 10;
#endif //ALLOW_BOLTZMANN_GIBBS
  
    // This is copied from ProfuseTest.hpp, to test with the same priors used there.
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> matchEmissionPrior;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> globalPrior;
    if( training_parameters_template.usePriors ) {
#ifdef USE_DEL_IN_DEL_OUT
      training_parameters_template.priorMtoM = .095;
      training_parameters_template.priorMtoI = .0025;
      training_parameters_template.priorMtoD = .0025;
      double priorMtoW                       = .9;
      double priorWtoW                       = .9999;
      double priorWtoE                       = .0001;
      double priorZtoZ                       = .9999;
      double priorZtoM                       = .0001;
      double priorBtoZ                       = .9;
      double priorBtoD                       = .01;
      double priorBtoM                       = .09;
#else
      double priorBtoD                       = .01;
      double priorBtoM                       = .9;
#endif // USE_DEL_IN_DEL_OUT .. else ..
#ifndef DISALLOW_FLANKING_TRANSITIONS
      double priorCtoC                       = .01;
      double priorCtoT                       = .99;
      double priorNtoN                       = .01;
      double priorNtoB                       = .99;
#endif // !DISALLOW_FLANKING_TRANSITIONS
      double priorStrength = 1;
      double priorStrength_internal_transitions = 1;
      double priorStrength_flanking_self_transitions = 1; // 1000
      double priorStrength_flanking_other_transitions = 1; // 1000
      // Here's how it is set up in ProfuseTest.hpp:
      matchEmissionPrior.reinitializeToEven( priorStrength );
      globalPrior.reinitializeToEven( priorStrength );
      training_parameters_template.matchEmissionPrior = &matchEmissionPrior;
      training_parameters_template.globalPrior = &globalPrior;
      // Put a strong prior on a small C->C and N->N
#ifndef DISALLOW_FLANKING_TRANSITIONS
      globalPrior[ 0 ][ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] = ( priorStrength_flanking_self_transitions * priorNtoN );
      globalPrior[ 0 ][ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] = ( priorStrength_flanking_other_transitions * priorNtoB );
      globalPrior[ 0 ][ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] = ( priorStrength_flanking_self_transitions * priorCtoC );
      globalPrior[ 0 ][ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] = ( priorStrength_flanking_other_transitions * priorCtoT );
#endif // !DISALLOW_FLANKING_TRANSITIONS

      // Do additional setup of the transition priors for those transtions
      // observed many times (since the number of transitions depends on
      // the profile length).
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toMatch ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorMtoM );
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorMtoD );
#ifdef USE_DEL_IN_DEL_OUT
      // TODO: Create training_parameters_template.priorMtoW
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] = ( initial_profile_length * priorStrength_internal_transitions * priorMtoW );
#endif // USE_DEL_IN_DEL_OUT
      globalPrior[ 0 ][ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorMtoI );
#ifdef USE_DEL_IN_DEL_OUT
      // TODO: Create training_parameters_template.priorWtoW
      globalPrior[ 0 ][ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ] = ( initial_profile_length * priorStrength_flanking_other_transitions * priorWtoW );
      // TODO: Create training_parameters_template.priorWtoE
      globalPrior[ 0 ][ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] = ( initial_profile_length * priorStrength_flanking_other_transitions * priorWtoE );
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toMatch ] = ( priorStrength_flanking_other_transitions * priorBtoM );
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] = ( priorStrength_flanking_other_transitions * priorBtoD );
      // TODO: Create training_parameters_template.priorBtoZ
      globalPrior[ 0 ][ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] = ( priorStrength_flanking_other_transitions * priorBtoZ );
      // TODO: Create training_parameters_template.priorZtoZ
      globalPrior[ 0 ][ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toDeletionIn ] = ( initial_profile_length * priorStrength_flanking_other_transitions * priorZtoZ );
      // TODO: Create training_parameters_template.priorZtoM
      globalPrior[ 0 ][ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] = ( initial_profile_length * priorStrength_flanking_other_transitions * priorZtoM );
#endif // USE_DEL_IN_DEL_OUT
      globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorItoM );
      globalPrior[ 0 ][ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorItoI );
      globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorDtoM );
      globalPrior[ 0 ][ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] = ( initial_profile_length * priorStrength_internal_transitions * training_parameters_template.priorDtoD );
    } // End if usePriors
  
    //training_parameters_template.verbosity = VERBOSITY_All;
    training_parameters_template.verbosity = VERBOSITY_Low;
    //training_parameters_template.verbosity = VERBOSITY_None;
    //training_parameters_template.debug = DEBUG_All;
      
    //training_parameters_template.minBaumWelchInverseScalar = (1.0/4.0);
      
    // TODO: REMOVE.  Testing matrixRowScaleFactor
    //training_parameters_template.matrixRowScaleFactor =
    //  pow( numeric_limits<double>::max(), .25 ) - 1.0;
  
    trainer.m_parameters = training_parameters_template;
    // TODO: REMOVE
    //cout << "just before training, trainer parameters are " << trainer.m_parameters << endl;
  
    if( true ) {
      cout << "The profile (before) is:" << endl;
      cout << *( trainer.m_profile ) << endl;

    } // End if print the before-profile
    // Do the training
    cout << "Training." << endl;
    score = trainer.train();

    cout << "Now (after training), the score is " << score << ", and the profile (length " << ( trainer.m_profile )->length() << ") is:" << endl;
    cout << *trainer.m_profile << endl;
    cout << "Its length is " << trainer.m_profile->length() << endl;
    cout << "It took " << trainer.m_iteration << " iterations." << endl;

    cout << "That score, again, was " << score << endl;
  
    if( calculate_viterbi_score || calculate_multiple_alignment ) {
      // Free memory first
      trainer.reinitialize();

      // Note that (for now) this will use the root profile only.
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer viterbi_matrices(
        *( trainer.m_profile ),
        fasta,
        num_sequences_to_use
      );
      score =
        dp.forward_score_viterbi(
          trainer.m_parameters,
          *( trainer.m_profile ),
          fasta,
          num_sequences_to_use,
          viterbi_matrices // TODO: ARE THESE MATRICES SIZED CORRECTLY?
        );
      cout << "The total score for all sequences, using viterbi, is: " << score << endl;
      // End calculating viterbi score
  
      if( calculate_multiple_alignment ) {
        // Show multiple alignment
        typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> ma(
          ( trainer.m_profile ),
          &fasta,
          num_sequences_to_use
        );
        dp.forward_viterbiAlign(
          trainer.m_parameters,
          viterbi_matrices,
          ma
        );
        cout << "Multiple Alignment is:" << endl;
        //cout << ma << endl;
        //ma.toPairwiseStream( cout, &fasta.m_descriptions );
        ma.toPileupStream( cout, &fasta.m_descriptions );
      } // End if calculate_multiple_alignment
    } // End if calculate_viterbi_score || calculate_multiple_alignment
    
    return score;
  } // train( string const &, Profile &, double const &, double const &, double const &, bool const &, uint32_t const &, double const &, double const &, bool const & )

}; // End class Train

} // End namespace galosh

#endif // __GALOSH_TRAIN_HPP__
