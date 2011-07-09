/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profillic
##  File:
##      ProfileTrainer.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the Galosh Profile HMM parameter estimator class.
##
#******************************************************************************
#*
#*    This file is part of profillic, a suite of programs for estimating parameters of
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

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROFILETRAINER_HPP__
#define __GALOSH_PROFILETRAINER_HPP__

#include "Profillic.hpp"

#include "ProlificParameters.hpp"

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "Profile.hpp"
using galosh::ProfileTreeRoot;

#include "Sequence.hpp"
using galosh::Sequence;

#include "Random.hpp"
using galosh::Random;

#include "DynamicProgramming.hpp"
using galosh::DynamicProgramming;

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
//#include <cmath> // for isnan( double )
//using std::isnan; // doesn't work !  why?

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

namespace galosh {

template <class ProfileType,
          class ScoreType,
          class MatrixValueType,
          class SequenceResidueType = typename profile_traits<ProfileType>::ResidueType,
          //class InternalNodeType = ProfileTreeInternalNode<typename profile_traits<ProfileType>::ProbabilityType> >
          class InternalNodeType = ProfileTreeRoot<typename profile_traits<ProfileType>::ResidueType, typename profile_traits<ProfileType>::ProbabilityType> >
  class ProfileTrainer {
  public:
    typedef typename profile_traits<ProfileType>::ResidueType ResidueType;
    typedef typename profile_traits<ProfileType>::ProbabilityType ProbabilityType;
    typedef Sequence<SequenceResidueType> SequenceType;

    class Parameters :
    public ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters
    {
      // Boost serialization
    private:
      typedef typename ProlificParameters<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::Parameters dynamic_programming_parameters_t;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( dynamic_programming_parameters_t );

        ar & BOOST_SERIALIZATION_NVP( trainProfileGlobals );
        ar & BOOST_SERIALIZATION_NVP( trainProfilePositions );
        ar & BOOST_SERIALIZATION_NVP( trainGlobalsFirst );
        ar & BOOST_SERIALIZATION_NVP( minIterations );
        ar & BOOST_SERIALIZATION_NVP( maxIterations );
        ar & BOOST_SERIALIZATION_NVP( maxPositionCycles );
        ar & BOOST_SERIALIZATION_NVP( maxPositionCycles_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( maxPositionCycles_globals );
        ar & BOOST_SERIALIZATION_NVP( maxBaumWelchInverseScalar );
        ar & BOOST_SERIALIZATION_NVP( minBaumWelchInverseScalar );
        ar & BOOST_SERIALIZATION_NVP( baumWelchInverseScalarIncrement );
        ar & BOOST_SERIALIZATION_NVP( maxBaumWelchInverseScalar_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( minBaumWelchInverseScalar_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( baumWelchInverseScalarIncrement_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( maxBaumWelchInverseScalar_globals );
        ar & BOOST_SERIALIZATION_NVP( minBaumWelchInverseScalar_globals );
        ar & BOOST_SERIALIZATION_NVP( baumWelchInverseScalarIncrement_globals );
        ar & BOOST_SERIALIZATION_NVP( scorePercentChangeMinimum_iteration );
        ar & BOOST_SERIALIZATION_NVP( scorePercentChangeMinimum_position_cycle );
        ar & BOOST_SERIALIZATION_NVP( euclideanDistanceMinimum_iteration );
        ar & BOOST_SERIALIZATION_NVP( euclideanDistanceMinimum_position_cycle );
        ar & BOOST_SERIALIZATION_NVP( alwaysAccept );
        ar & BOOST_SERIALIZATION_NVP( useAlignmentProfiles );
        ar & BOOST_SERIALIZATION_NVP( proposeProfileLengthChanges );
        ar & BOOST_SERIALIZATION_NVP( proposeDeletingThreshold );
        ar & BOOST_SERIALIZATION_NVP( proposeInsertingThreshold );
        ar & BOOST_SERIALIZATION_NVP( proposeDeletingThreshold_increment );
        ar & BOOST_SERIALIZATION_NVP( proposeInsertingThreshold_increment );
        ar & BOOST_SERIALIZATION_NVP( proposeInsertingPreAlignThreshold );
        ar & BOOST_SERIALIZATION_NVP( proposeInsertingPostAlignThreshold );
        ar & BOOST_SERIALIZATION_NVP( proposeInsertingOccupancyThreshold );
        ar & BOOST_SERIALIZATION_NVP( useSensitiveThresholding );
        ar & BOOST_SERIALIZATION_NVP( increaseThresholdsForLengthChanges_startIteration );
        ar & BOOST_SERIALIZATION_NVP( increaseThresholdsForLengthChanges_minIncrement );
        ar & BOOST_SERIALIZATION_NVP( alwaysAccept_disallowThreshold_profileDistance_iteration );
        ar & BOOST_SERIALIZATION_NVP( numIterationsBetweenLengthChanges );
        ar & BOOST_SERIALIZATION_NVP( useUnconditionalBaumWelch );
        ar & BOOST_SERIALIZATION_NVP( baldiLearningRate );
        ar & BOOST_SERIALIZATION_NVP( baldiTemperature );
        ar & BOOST_SERIALIZATION_NVP( profileValueMinimum );
        ar & BOOST_SERIALIZATION_NVP( usePriors );
        ar & BOOST_SERIALIZATION_NVP( positionShouldBeTrained );
        ar & BOOST_SERIALIZATION_NVP( siegelMaxFindingThePeakAttempts_positions );
        ar & BOOST_SERIALIZATION_NVP( baldiHybrid );
        ar & BOOST_SERIALIZATION_NVP( siegelEpsilonScaleFactor );
        ar & BOOST_SERIALIZATION_NVP( unconditionalIsolatesGlobals );
        ar & BOOST_SERIALIZATION_NVP( siegelMaxRefiningThePeakSteps_positions );
        ar & BOOST_SERIALIZATION_NVP( siegelRefiningThePeakStepsConvergenceThreshold );
        ar & BOOST_SERIALIZATION_NVP( siegelMaxFindingTheGradientAttempts_positions );
        ar & BOOST_SERIALIZATION_NVP( siegelMinEpsilon );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// PARAMETERS
      /**
       * Should we train the global parameters at all, or just do positions?
       * If true, train the globals.  If false, don't train the globals.
       */
      bool trainProfileGlobals;
  #define DEFAULT_trainProfileGlobals true
  
      /**
       * Should we train the position parameters at all?  If true, train the
       * position-specific parameters.  If false, don't train them.
       */
      bool trainProfilePositions;
  #define DEFAULT_trainProfilePositions true
  
      /**
       * Should we start with the global parameters, then do positions?  If
       * true, start with globals.  If false, start with position params,
       * then do globals.
       */
      bool trainGlobalsFirst;
  #define DEFAULT_trainGlobalsFirst false
  
      /**
       * How many times at minimum should we iterate between position cycles and
       * global training steps before we end the training process ( one
       * iteration includes both positions and globals )?
       * Note that this must be at least 1.
       */
      uint32_t minIterations; // at least 1
  #define DEFAULT_minIterations 1
  
      /**
       * How many times at maximum should we iterate between position cycles and
       * global training steps before we end the training process ( one
       * iteration includes both positions and globals )?
       * Note that this must be at least 1.
       */
      uint32_t maxIterations; // at least 1
  #define DEFAULT_maxIterations 30
  
      /**
       * How many times should we cycle through the training steps before,
       * despite continuing improvement, we move on to refining indels?
       */
      uint32_t maxPositionCycles;
  #define DEFAULT_maxPositionCycles 4

      /**
       * How many times should we cycle through the training steps before,
       * despite continuing improvement, we move on to refining indels?
       *
       * This is for the sequence identifier params training step.
       */
      uint32_t maxPositionCycles_sequence_identifiers;
  #define DEFAULT_maxPositionCycles_sequence_identifiers 4

      /**
       * How many times should we cycle through the training steps before,
       * despite continuing improvement, we move on to refining indels?
       *
       * This is for the global params training step.
       */
      uint32_t maxPositionCycles_globals;
  #define DEFAULT_maxPositionCycles_globals 4

      /**
       * Sometimes we have to reject the conditional BW updates because the
       * score would go down.  In such cases we try again at various
       * intermediate steps.  We do this by adding some fraction of the updated
       * position values to the original values.  When this fraction is 1, for
       * instance, the old and new values will be averaged.  When the fraction
       * is less than 1, it will be a weighted average.  The value of this
       * parameter is the largest denomenator of the weight on the new values.
       * We try denomenators from
       * minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by
       * baumWelchInverseScalarIncrement.  If minBaumWelchInverseScalar is
       * greater than 0, we don't even try the usual direct update.
       */
      float maxBaumWelchInverseScalar;
  #define DEFAULT_maxBaumWelchInverseScalar 8.0F

      /**
       * See maxBaumWelchInverseScalar.
       */
      float minBaumWelchInverseScalar;
  #define DEFAULT_minBaumWelchInverseScalar 0.0F

      /**
       * See maxBaumWelchInverseScalar.
       */
      float baumWelchInverseScalarIncrement;
  #define DEFAULT_baumWelchInverseScalarIncrement 4.0F

      /**
       * Sometimes we have to reject the conditional BW updates because the
       * score would go down.  In such cases we try again at various
       * intermediate steps.  We do this by adding some fraction of the updated
       * position values to the original values.  When this fraction is 1, for
       * instance, the old and new values will be averaged.  When the fraction
       * is less than 1, it will be a weighted average.  The value of this
       * parameter is the largest denomenator of the weight on the new values.
       * We try denomenators from
       * minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by
       * baumWelchInverseScalarIncrement.  If minBaumWelchInverseScalar is
       * greater than 0, we don't even try the usual direct update.
       *
       * This is for the sequence identifier params training step.
       */
      float maxBaumWelchInverseScalar_sequence_identifiers;
  #define DEFAULT_maxBaumWelchInverseScalar_sequence_identifiers 40.0F

      /**
       * See maxBaumWelchInverseScalar_sequence_identifiers.
       *
       * This is for the sequence identifier params training step.
       */
      float minBaumWelchInverseScalar_sequence_identifiers;
  #define DEFAULT_minBaumWelchInverseScalar_sequence_identifiers 0.0F

      /**
       * See maxBaumWelchInverseScalar_sequence_identifiers.
       *
       * This is for the sequence identifier params training step.
       */
      float baumWelchInverseScalarIncrement_sequence_identifiers;
  #define DEFAULT_baumWelchInverseScalarIncrement_sequence_identifiers 10.0F

      /**
       * Sometimes we have to reject the conditional BW updates because the
       * score would go down.  In such cases we try again at various
       * intermediate steps.  We do this by adding some fraction of the updated
       * position values to the original values.  When this fraction is 1, for
       * instance, the old and new values will be averaged.  When the fraction
       * is less than 1, it will be a weighted average.  The value of this
       * parameter is the largest denomenator of the weight on the new values.
       * We try denomenators from
       * minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by
       * baumWelchInverseScalarIncrement.  If minBaumWelchInverseScalar is
       * greater than 0, we don't even try the usual direct update.
       *
       * This is for the global params training step.
       */
      float maxBaumWelchInverseScalar_globals;
  #define DEFAULT_maxBaumWelchInverseScalar_globals 100.0F

      /**
       * See maxBaumWelchInverseScalar_globals.
       *
       * This is for the global params training step.
       */
      float minBaumWelchInverseScalar_globals;
  #define DEFAULT_minBaumWelchInverseScalar_globals 0.0F

      /**
       * See maxBaumWelchInverseScalar_globals.
       *
       * This is for the global params training step.
       */
      float baumWelchInverseScalarIncrement_globals;
  #define DEFAULT_baumWelchInverseScalarIncrement_globals 20.0F

      /**
       * How much (percentage) change in forward score after an iteration is
       * close enough to 0 that we should stop training?
       */
      float scorePercentChangeMinimum_iteration;
  #define DEFAULT_scorePercentChangeMinimum_iteration 1
    
      /**
       * How much (percentage) change in forward score in a position cycle is
       * close enough to 0 that we should move on to refining indels?
       */
      float scorePercentChangeMinimum_position_cycle;
  #define DEFAULT_scorePercentChangeMinimum_position_cycle 1

      /**
       * What average squared euclidean distance (averaged over free
       * parameters) is close enough to 0 that we should stop training?
       */
      float euclideanDistanceMinimum_iteration;
  #define DEFAULT_euclideanDistanceMinimum_iteration 1E-5

      /**
       * What average squared euclidean distance (averaged over free
       * parameters) is close enough to 0 that we should move on to refining
       * indels?
       */
      float euclideanDistanceMinimum_position_cycle;
  #define DEFAULT_euclideanDistanceMinimum_position_cycle 1E-5

      /**
       * Always accept the changes to the profile, even if the score goes down?
       * Theoretically, the BW update should always lead to an improvement, but
       * because of numerical instability (or strong priors), sometimes it
       * won't.
       *
       * See alwaysAccept_disallowThreshold_profileDistance_iteration.
       */
      bool alwaysAccept;
  #define DEFAULT_alwaysAccept true

      /**
       * Calculate and use Alignment Profiles for the positions BW update?  If
       * true, will do so; if false, will use
       * calculatePositionSpecificSequenceScoreCoefficients(..) and
       * updatePositionEntenteForSequence(..) instead.
       *
       * Note that alignment profiles will be calculated if
       * proposeProfileLengthChanges is true, regardless of the value of this
       * parameter.  They will only be used for the BW update, though, if this
       * is true.
       */
      bool useAlignmentProfiles;
  #define DEFAULT_useAlignmentProfiles true

      /**
       * Resize the profile while training it?  If this is true, alignment
       * profiles will be calculated and used to delete or insert positions.
       * 'propose' is a bit of a misnomer: if the rate of deletion or insertion
       * exceeds the corresponding threshold, that position will be deleted or
       * inserted (unless a cycle has been detected: there's some fancy
       * footwork involved in preventing infinite loops here: see
       * proposeInsertingThreshold_increment).
       *
       * See also proposeDeletingThreshold, proposeInsertingThreshold, and
       * numIterationsBetweenLengthChanges.
       */
      bool proposeProfileLengthChanges;
  #define DEFAULT_proposeProfileLengthChanges false

      /**
       * If the alignment profile (for a particular sequence) at a position
       * position has higher that this threshold of deletions, then the
       * sequence will be counted towards the total number of sequences with
       * "deletes" at that position -- and if the fraction of sequences exceeds
       * one minus the proposeInsertingOccupancyThreshold, the position will be
       * deleted (if proposeProfileLengthChanges is true).
       *
       * See proposeProfileLengthChanges, proposeInsertingOccupancyThreshold.
       */
      double proposeDeletingThreshold;
#define DEFAULT_proposeDeletingThreshold .1

      /**
       * If the alignment profile (for a particular sequence) at an internal
       * position has higher that this threshold of insertions, then
       * the sequence will be counted towards the total number of sequences
       * with "inserts" at that position -- and if the fraction of sequences
       * exceeds the proposeInsertingOccupancyThreshold, a new position will be
       * inserted (if proposeProfileLengthChanges is true).
       *
       * See proposeProfileLengthChanges, proposeInsertingPreAlignThreshold,
       * proposeInsertingPostAlignThreshold, proposeInsertingOccupancyThreshold.
       */
      double proposeInsertingThreshold;
#define DEFAULT_proposeInsertingThreshold .1

      /**
       * If we detect a cycle, we increment the proposeDeletingThreshold by
       * this amount.  We also increment it by this amount if it gets to the
       * minimum or maximum profile length (the minimum is 2, the maximum is
       * the longest sequence being used to train the profile).  See also
       * proposeInsertingThreshold_increment.  If useSensitiveThresholding is
       * true, thresholds will instead be incremented to the least amount that
       * would have made a difference in which positions were inserted or
       * deleted -- and this value will be used added as a little extra boost.
       *
       * See proposeProfileLengthChanges, useSensitiveThresholding.
       */
      double proposeDeletingThreshold_increment;
#define DEFAULT_proposeDeletingThreshold_increment .025

      /**
       * If we detect a cycle, we increment the proposeInsertingThresholds by
       * this amount.  We also increment them by this amount if the length gets
       * to the minimum or maximum profile length (the minimum is 2, the
       * maximum is the longest sequence being used to train the profile).  See
       * also proposeDeletingThreshold_increment.  If useSensitiveThresholding
       * is true, thresholds will instead be incremented to the least amount
       * that would have made a difference in which positions were inserted or
       * deleted -- and this value will be used added as a little extra boost.
       *
       * See proposeProfileLengthChanges, useSensitiveThresholding.
       */
      double proposeInsertingThreshold_increment;
#define DEFAULT_proposeInsertingThreshold_increment .025

      /**
       * If the alignment profile (for a particular sequence) for the
       * first/pre-align position has higher that this threshold of pre-align
       * insertions, then the sequence will be counted towards the total number
       * of sequences with "inserts" at that position -- and if the fraction of
       * sequences exceeds the proposeInsertingOccupancyThreshold, a new
       * position will be inserted (if proposeProfileLengthChanges is true).
       *
       * See proposeProfileLengthChanges, proposeInsertingThreshold,
       * proposeInsertingPostAlignThreshold, proposeInsertingOccupancyThreshold.
       */
      double proposeInsertingPreAlignThreshold;
#define DEFAULT_proposeInsertingPreAlignThreshold .1

      /**
       * If the alignment profile (for a particular sequence) at the last
       * position has higher that this threshold of post-align insertions, then
       * the sequence will be counted towards the total number of sequences
       * with "inserts" at that position -- and if the fraction of sequences
       * exceeds the proposeInsertingOccupancyThreshold, a new position will be
       * inserted (if proposeProfileLengthChanges is true).
       *
       * See proposeProfileLengthChanges, proposeInsertingThreshold, and
       * proposeInsertingPreAlignThreshold, proposeInsertingOccupancyThreshold.
       */
      double proposeInsertingPostAlignThreshold;
#define DEFAULT_proposeInsertingPostAlignThreshold .1

      /**
       * If the fraction of sequences for which the insertion fraction exceeds
       * the relevant threshold is greater than
       * proposeInsertingOccupancyThreshold (for a particular position), then a
       * new position will be inserted.  If the fraction of sequences for which
       * the deletion fraction exceeds the proposeDeletingThreshold is greater
       * than ( 1.0 - proposeInsertingOccupancyThreshold ), then the position
       * will be deleted (if proposeProfileLengthChanges is true).
       *
       * See proposeProfileLengthChanges, proposeInsertingThreshold, and
       * proposeInsertingPreAlignThreshold, proposeInsertingPostAlignThreshold,
       * and proposeDeletingThreshold.
       */
      double proposeInsertingOccupancyThreshold;
#define DEFAULT_proposeInsertingOccupancyThreshold .5

      /**
       * If true, thresholds will not be incremented until the end of a cycle,
       * and then will be incremented to the least value that would have made a
       * difference in the positions being inserted or deleted.  Note that the
       * threshold increment will also be added to this calculated min value,
       * so when useSensitiveThresholding is true those increments should be
       * small (they provide a bit of a bump to offset the fact that after the
       * global parameters update, the indel fractions tend to increase a
       * touch).
       */
      bool useSensitiveThresholding;
#define DEFAULT_useSensitiveThresholding true

      /**
       * If things are taking too long, it might be because there's an
       * undetected cycle.  At
       * increaseThresholdsForLengthChanges_startIteration, we start treating
       * all length changes like they're part of a cycle, to force convergence.
       *
       * See increaseThresholdsForLengthChanges_minIncrement.
       */
      uint32_t increaseThresholdsForLengthChanges_startIteration;
#define DEFAULT_increaseThresholdsForLengthChanges_startIteration 500

      /**
       * If things are taking too long, it might be because there's an
       * undetected cycle.  At
       * increaseThresholdsForLengthChanges_startIteration, we start treating
       * all length changes like they're part of a cycle, to force convergence.
       * If the threshold increments are very small, this might not help much,
       * so we also enforce a minimum increment.
       *
       * See increaseThresholdsForLengthChanges_startIteration.
       */
      float increaseThresholdsForLengthChanges_minIncrement;
#define DEFAULT_increaseThresholdsForLengthChanges_minIncrement 1E-4

      /**
       * When the profile isn't changing much, numerical errors dominate.  The
       * BW and CBW algs are theoretically guaranteed to only improve the
       * score, but these numerical errors can cause the score to diminish a
       * bit, and convergence to take a very long time.  When the profile
       * distance (in an iteration) is less than
       * alwaysAccept_disallowThreshold_profileDistance_iteration, we turn
       * alwaysAccept off, so the score never goes down and the algorithm
       * converges.
       *
       * See alwaysAccept.
       */
      float alwaysAccept_disallowThreshold_profileDistance_iteration;
#define DEFAULT_alwaysAccept_disallowThreshold_profileDistance_iteration 1E-5

      /**
       * if proposeProfileLengthChanges is true, this number of iterations must
       * pass between proposals.  It may be 0.  Regardless of its value, at
       * least 1 iteration must pass at the very beginning of training before
       * length changes are allowed.
       */
      uint32_t numIterationsBetweenLengthChanges;
#define DEFAULT_numIterationsBetweenLengthChanges 0

      /**
       * Use the usual kind of Baum-Welch, in which parameter updates are
       * simultaneous?
       * If true, updates all parameters at once via Baum-Welch.
       * Otherwise, performs a conditional optimization.
       */
      bool useUnconditionalBaumWelch;
  #define DEFAULT_useUnconditionalBaumWelch false
    
      /**
       * Use Baldi's style of gradient ascent (GEM vs EM)?  If so, with what
       * learning rate?
       *
       * If non-zero, Baldi's "Smooth on-line learning algorithm" will be used,
       * with the given learning rate.
       *
       * @see baldiTemperature
       */
      double baldiLearningRate;
  #define DEFAULT_baldiLearningRate 0
    
      /**
       * If we are using Baldi's style of gradient ascent (GEM vs EM), use this
       * temperature.
       *
       * @see baldiLearningRate
       */
      double baldiTemperature;
  #define DEFAULT_baldiTemperature 1.0
    
      /**
       * Minimum profile distribution value.  Applies to everything except the
       * E (End) state.  Given as a double, despite the Profile ProbabilityType.
       */
      double profileValueMinimum;
  #define DEFAULT_profileValueMinimum 1E-5

      /**
       * Use priors?  If true, m_matchEmissionPrior and m_globalsPrior will be
       * defined and used.  See also matchEmissionPrior and globalsPrior.
       */
      bool usePriors;
  #define DEFAULT_usePriors true

    /**
     * When usePriors is true, and when this is non-null, we use this for the
     * Match-state emission prior.  If this is null, a simple Laplace prior
     * will be used.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> * matchEmissionPrior;

    /**
     * When usePriors is true, and when this is non-null, we use this for the
     * priors for everything except the match-state emissions.  If this is
     * null, a simple Laplace prior will be used.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> * globalPrior;

      /**
       * It is possible to turn off training on a position-by-position basis.
       * If this vector is non-null, it must be of the same length as the
       * profile, and the truth value of the i^th position determines whether
       * the i^th profile position should be trained (true means yes, train
       * that position).
       */
      vector<bool> * positionShouldBeTrained;
  #define DEFAULT_positionShouldBeTrained NULL

      /**
       * Set this to 0 for straight-up Baldi.  If non-zero, use Baldi / Siegel and quit trying to find a peak after this many attempts.  Note this only applies when baldiLearningRate > 0 and when ALLOW_BOLTZMANN_GIBBS is #defined.
       */
      uint32_t siegelMaxFindingThePeakAttempts_positions;
#define DEFAULT_siegelMaxFindingThePeakAttempts_positions 1000

      /**
       * Set to true to try both BW / CBW *and* Baldi/Siegel, then take whichever move is best.
       */
      bool baldiHybrid;
#define DEFAULT_baldiHybrid false

      /**
       * Length of steps for Baldi / Siegel (this is multiplied into the learning rate from eqn 3 from the Baldi paper to explore trapping a peak).  Note steps down will use the inverse of this.  This must be strictly > 1.
       */
      double siegelEpsilonScaleFactor;
#define DEFAULT_siegelEpsilonScaleFactor 2.0

      /**
       * If doing Unconditional BW or Unconditional Baldi, should we still isolate updates to the globals from updates to the position params?  This is a hybrid approach, conditional but not per-position.
       */
      bool unconditionalIsolatesGlobals;
#define DEFAULT_unconditionalIsolatesGlobals false

      /**
       * This should generally be set to a fairly high value; it determines the number of attempts to maximize the quadratic after trapping a peak in the Baldi / Siegel Newton-Raphson part.  The more relevant parameter to adjust is siegelRefiningThePeakStepsConvergenceThreshold.
       */
      uint32_t siegelMaxRefiningThePeakSteps_positions;
#define DEFAULT_siegelMaxRefiningThePeakSteps_positions 1000

      /**
       * When (the larger minus the smaller, divided by the smaller) of the last-best and the current-best score is less than this threshold, declare the Baldi / Siegel hillclimb done.
       */
      ScoreType siegelRefiningThePeakStepsConvergenceThreshold;
#define DEFAULT_siegelRefiningThePeakStepsConvergenceThreshold 1E-5

      /**
       * If the score doesn't change after this number of attempts to multiply the gradient by an epsilon, give up.  Note that this must be > 0. Note this only applies when baldiLearningRate > 0, when siegelMaxFindingThePeakAttempts_positions > 0, and when ALLOW_BOLTZMANN_GIBBS is #defined.
       */
      uint32_t siegelMaxFindingTheGradientAttempts_positions;
#define DEFAULT_siegelMaxFindingTheGradientAttempts_positions 10

      /**
       * Give up trapping the peak if the epsilon value we're multiplying by the gradient gets below this value (we multiply it by siegelEpsilonScaleFactor or 1/SiegelEpsilonScaleFactor as we go..).  Note that this must be > 0. Note this only applies when baldiLearningRate > 0, when siegelMaxFindingThePeakAttempts_positions > 0, and when ALLOW_BOLTZMANN_GIBBS is #defined.
       */
      double siegelMinEpsilon;
#define DEFAULT_siegelMinEpsilon 1E-5

      Parameters ();
      virtual ~Parameters () {};
    
      // Copy constructor
      template <class AnyParameters>
      Parameters ( const AnyParameters & copy_from );
    
      // Copy constructor/operator
      template <class AnyParameters>
      Parameters & operator= (
        AnyParameters const& copy_from
      );

      template <class AnyParameters>
      void
      copyFromNonVirtual (
        AnyParameters const & copy_from
      );

      template <class AnyParameters>
      void
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      );

      virtual void
      copyFrom ( const Parameters & copy_from );
    
      virtual void
      resetToDefaults ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        Parameters const& parameters
      )
      {
        parameters.writeParameters( os );

        return os;
      } // friend operator<< ( basic_ostream &, Parameters const& )

      template<class CharT, class Traits>
      void
      writeParameters ( 
        std::basic_ostream<CharT,Traits>& os
      ) const;

    }; // End inner class Parameters

    template <class ParametersType>
    class ParametersModifierTemplate :
      public ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template ParametersModifierTemplate<ParametersType>
    {
      typedef typename ProlificParameters<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t;

      // Boost serialization
    private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information.  This will serialize the
        // parameters too.
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( base_parameters_modifier_t );

        // Serialize the new isModified_ stuff
        ar & BOOST_SERIALIZATION_NVP( isModified_trainProfileGlobals );
        ar & BOOST_SERIALIZATION_NVP( isModified_trainProfilePositions );
        ar & BOOST_SERIALIZATION_NVP( isModified_trainGlobalsFirst );
        ar & BOOST_SERIALIZATION_NVP( isModified_minIterations );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxIterations );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxPositionCycles );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxPositionCycles_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxPositionCycles_globals );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxBaumWelchInverseScalar );
        ar & BOOST_SERIALIZATION_NVP( isModified_minBaumWelchInverseScalar );
        ar & BOOST_SERIALIZATION_NVP( isModified_baumWelchInverseScalarIncrement );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxBaumWelchInverseScalar_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( isModified_minBaumWelchInverseScalar_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( isModified_baumWelchInverseScalarIncrement_sequence_identifiers );
        ar & BOOST_SERIALIZATION_NVP( isModified_maxBaumWelchInverseScalar_globals );
        ar & BOOST_SERIALIZATION_NVP( isModified_minBaumWelchInverseScalar_globals );
        ar & BOOST_SERIALIZATION_NVP( isModified_baumWelchInverseScalarIncrement_globals );
        ar & BOOST_SERIALIZATION_NVP( isModified_scorePercentChangeMinimum_iteration );
        ar & BOOST_SERIALIZATION_NVP( isModified_scorePercentChangeMinimum_position_cycle );
        ar & BOOST_SERIALIZATION_NVP( isModified_euclideanDistanceMinimum_iteration );
        ar & BOOST_SERIALIZATION_NVP( isModified_euclideanDistanceMinimum_position_cycle );
        ar & BOOST_SERIALIZATION_NVP( isModified_alwaysAccept );
        ar & BOOST_SERIALIZATION_NVP( isModified_useAlignmentProfiles );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeProfileLengthChanges );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeDeletingThreshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeInsertingThreshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeDeletingThreshold_increment );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeInsertingThreshold_increment );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeInsertingPreAlignThreshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeInsertingPostAlignThreshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_proposeInsertingOccupancyThreshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_useSensitiveThresholding );
        ar & BOOST_SERIALIZATION_NVP( isModified_increaseThresholdsForLengthChanges_startIteration );
        ar & BOOST_SERIALIZATION_NVP( isModified_increaseThresholdsForLengthChanges_minIncrement );
        ar & BOOST_SERIALIZATION_NVP( isModified_alwaysAccept_disallowThreshold_profileDistance_iteration );
        ar & BOOST_SERIALIZATION_NVP( isModified_numIterationsBetweenLengthChanges );
        ar & BOOST_SERIALIZATION_NVP( isModified_useUnconditionalBaumWelch );
        ar & BOOST_SERIALIZATION_NVP( isModified_baldiLearningRate );
        ar & BOOST_SERIALIZATION_NVP( isModified_baldiTemperature );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileValueMinimum );
        ar & BOOST_SERIALIZATION_NVP( isModified_usePriors );
        ar & BOOST_SERIALIZATION_NVP( isModified_positionShouldBeTrained );
        ar & BOOST_SERIALIZATION_NVP( isModified_siegelMaxFindingThePeakAttempts_positions );
        ar & BOOST_SERIALIZATION_NVP( isModified_baldiHybrid );
        ar & BOOST_SERIALIZATION_NVP( isModified_siegelEpsilonScaleFactor );
        ar & BOOST_SERIALIZATION_NVP( isModified_unconditionalIsolatesGlobals );
        ar & BOOST_SERIALIZATION_NVP( isModified_siegelMaxRefiningThePeakSteps_positions );
        ar & BOOST_SERIALIZATION_NVP( isModified_siegelRefiningThePeakStepsConvergenceThreshold );
        ar & BOOST_SERIALIZATION_NVP( isModified_siegelMaxFindingTheGradientAttempts_positions );
        ar & BOOST_SERIALIZATION_NVP( isModified_siegelMinEpsilon );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      bool isModified_trainProfileGlobals;
  
      bool isModified_trainProfilePositions;
  
      bool isModified_trainGlobalsFirst;
  
      bool isModified_minIterations; // at least 1
  
      bool isModified_maxIterations; // at least 1
  
      bool isModified_maxPositionCycles;

      bool isModified_maxPositionCycles_sequence_identifiers;

      bool isModified_maxPositionCycles_globals;

      bool isModified_maxBaumWelchInverseScalar;

      bool isModified_minBaumWelchInverseScalar;

      bool isModified_baumWelchInverseScalarIncrement;

      bool isModified_maxBaumWelchInverseScalar_sequence_identifiers;

      bool isModified_minBaumWelchInverseScalar_sequence_identifiers;

      bool isModified_baumWelchInverseScalarIncrement_sequence_identifiers;

      bool isModified_maxBaumWelchInverseScalar_globals;

      bool isModified_minBaumWelchInverseScalar_globals;

      bool isModified_baumWelchInverseScalarIncrement_globals;

      bool isModified_scorePercentChangeMinimum_iteration;
    
      bool isModified_scorePercentChangeMinimum_position_cycle;

      bool isModified_euclideanDistanceMinimum_iteration;

      bool isModified_euclideanDistanceMinimum_position_cycle;

      bool isModified_alwaysAccept;

      bool isModified_useAlignmentProfiles;

      bool isModified_proposeProfileLengthChanges;

      bool isModified_proposeDeletingThreshold;

      bool isModified_proposeInsertingThreshold;

      bool isModified_proposeDeletingThreshold_increment;
      bool isModified_proposeInsertingThreshold_increment;

      bool isModified_proposeInsertingPreAlignThreshold;

      bool isModified_proposeInsertingPostAlignThreshold;

      bool isModified_proposeInsertingOccupancyThreshold;

      bool isModified_useSensitiveThresholding;

      bool isModified_increaseThresholdsForLengthChanges_startIteration;

      bool isModified_increaseThresholdsForLengthChanges_minIncrement;

      bool isModified_alwaysAccept_disallowThreshold_profileDistance_iteration;

      bool isModified_numIterationsBetweenLengthChanges;

      bool isModified_useUnconditionalBaumWelch;
    
      bool isModified_baldiLearningRate;
    
      bool isModified_baldiTemperature;
    
      bool isModified_profileValueMinimum;

      bool isModified_usePriors;

      bool isModified_positionShouldBeTrained;

      bool isModified_siegelMaxFindingThePeakAttempts_positions;

      bool isModified_baldiHybrid;

      bool isModified_siegelEpsilonScaleFactor;

      bool isModified_unconditionalIsolatesGlobals;

      bool isModified_siegelMaxRefiningThePeakSteps_positions;

      bool isModified_siegelRefiningThePeakStepsConvergenceThreshold;

      bool isModified_siegelMaxFindingTheGradientAttempts_positions;

      bool isModified_siegelMinEpsilon;

      ParametersModifierTemplate ();
    
      // Copy constructor
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from );
    
      // Copy constructor/operator
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate & operator= (
        AnyParametersModifierTemplate const& copy_from
      );

      template <class AnyParametersModifierTemplate>
      void
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );
        
      template <class AnyParametersModifierTemplate>
      void
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      void reset ();

      void
      isModified_reset ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        ParametersModifierTemplate const& parameters_modifier
      )
      {
        parameters_modifier.writeParametersModifier( os );

        return os;
      } // friend operator<< ( basic_ostream &, ParametersModifierTemplate const& )

      template<class CharT, class Traits>
      void
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      );

      template<class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters );

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename ProfileTrainer::Parameters> ParametersModifier;

    class Error {
    public:
      ProfileTrainer * m_trainer; // The enclosing state.
    
      /**
       * The different types of error that the ProfileTrainer::ProfileTrainer could be in.
       * These should be consecutive from 0, and at the corresponding index in
       * c_defaultDescription there should be some string entry (maybe "").
       */
      enum type_t {
        NoError = 0,
        IncosistentStateError = 1
      } m_type;
    
      // TODO: Figure out how to get a static, constant array of strings
      #define DEFAULT_DESCRIPTION_0 "There is no known problem with the ProfileTrainer instance."
      #define DEFAULT_DESCRIPTION_1 "The ProfileTrainer instance appears to be in an incosistent state."
      /**
       * The (string) description of this error.
       */
      string m_description;
    
      Error (
        ProfileTrainer * trainer,
        type_t type = NoError,
        string description = DEFAULT_DESCRIPTION_0 // TODO: Make this depend
                                                        // on type
      ) :
        m_trainer( trainer ),
        m_type( type ),
        m_description( description )
      {
        // Do nothing else
      } // <init>( trainer, [ type_t, [ string ] ] )
    
      ~Error () {
        // Do nothing
      } // <destroy>()
    
      // TODO: Cast as int to the enum, and as str to the desc.
    }; // End inner class Error

    static int const trainingPhaseCount = 2;
    typedef uint8_t TrainingPhase;
    static int const TRAINING_PHASE_Positions = 0;
    static int const TRAINING_PHASE_Globals = 1;

    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType> m_dynamic_programming;
    Parameters m_parameters;
    ProfileType * m_profile;
    vector<SequenceType> const & m_sequences;

    typedef Error error_t;
    
    error_t m_error;
    Parameters m_trainingParameters;
    
    /**
     * This profile is used to store globals to revert to during
     * lengthadjust. It gets altered every time there is a reversion: the
     * globals at that time are added into this, weighed by the iteration
     * number.
     */
    ProfileType m_globalsToRevertTo;

    /**
     * This is the starting score.  It is not changed once it is set.
     */
    ScoreType m_scoreBeforeTraining;
    
    /** 
     * The score after the performing the current or most recent training
     * step.
     */
    ScoreType m_endingScore;
    
    /** 
     * The score at the beginning of the current iteration.
     */
    ScoreType m_startingScore_iteration;

    /** 
     * The score at the beginning of the current position cycle.
     */
    ScoreType m_startingScore_position_cycle;

    /** 
     * The score at the beginning of the current position training step.
     */
    ScoreType m_startingScore_position_step;

    /** 
     * The score at the beginning of the current position training substep,
     * when that substep is training the position-specific mixing parameters.
     */
    ScoreType m_startingScore_position_step_mixingParameters;

    /** 
     * The current position cycle (from 0)
     */
    uint32_t m_position_cycle;

    /**
     * Sometimes we have to reject the conditional BW updates because the score
     * would go down.  In such cases we try again at various intermediate
     * steps.  We do this by adding some fraction of the updated position
     * values to the original values.  When this fraction is 1, for instance,
     * the old and new values will be averaged.  When the fraction is less than
     * 1, it will be a weighted average.  We try denomenators from
     * minBaumWelchInverseScalar..maxBaumWelchInverseScalar by
     * baumWelchInverseScalarIncrement.  If minBaumWelchInverseScalar > 0, we
     * don't even try the usual unscaled update.  This is the current
     * denomenator.
     */
    float m_bw_inverse_scalar;
    
    /** 
     * The current position of the profile (from the end, backwards)
     */
    uint32_t m_row_i;

    /**
     * The index of the current sequence
     */
    uint32_t m_seq_i;

    /**
     * The number of sequences to use in training (the index one greater than
     * that of the last one to be used; must be <= m_sequences.size().
     */
    // TODO: REMOVE CONST
    const uint32_t m_sequence_count;

    /**
     * The percent change between m_startingScore_* and m_endingScore (as a
     * percent of m_startingScore_*).  Note that this could be negative (if
     * the score decreases), and it could be a very large or a very small
     * number.  See also percentChange( a, b ).
     */
    double m_scorePercentChange;
    
    /**
     * The current iteration, from 0..m_trainingParameters.maxIterations.
     */
    uint32_t m_iteration;

    /**
     * Which set of parameters are we currently training?
     */
    TrainingPhase m_trainingPhase;

    /**
     * The profile at the beginning of the current iteration.
     */
    ProfileType m_startingProfile_iteration;

    /**
     * The profile at the beginning of the current position cycle.
     */
    ProfileType m_startingProfile_position_cycle;

    /**
     * The forward rows, indexed first by row, then by sequence.
     */
    // NOTE: No longer used.
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer m_forward_matrices;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector m_forward_rows_1;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector m_forward_rows_2;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector m_unconditional_last_forward_rows;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector m_unconditional_second_to_last_forward_rows;

    /**
     * A pointer to the forward_rows we're using ( for profile position
     * m_row_i - 1 )
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector * m_forward_rows_ptr;
    /**
     * A pointer to the forward_rows for the previous row...
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector * m_prev_forward_rows_ptr;

    /**
     * Temp forward_rows ptr
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector * m_temp_forward_rows_ptr;

    /**
     * The anchor columns, for forward_reverseCalculateRow(..).
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer m_anchor_columns;

    /**
     * The anchor rows, indexed first by row, then by sequence.  This is the
     * subset of the forward matrix that we store to account for numerical
     * drift.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer m_anchor_rows;

    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector m_backward_rows_1;
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector m_backward_rows_2;

    /**
     * A pointer to the backward_rows we're using ( for profile position
     * m_row_i - 1 )
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector * m_backward_rows_ptr;
    /**
     * A pointer to the backward_rows for the following row...
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector * m_next_backward_rows_ptr;
    /**
     * Temp backward_rows ptr
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::RowVector * m_temp_backward_rows_ptr;

    /**
     * The global_entente is the data structure that holds the intermediate
     * values we need for Baum-Welch update of the global params.
     */    
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::GlobalEntente m_global_entente;

    /**
     * This is the new way of storing entente information.  It only stores the
     * entente for one position (the current position), is only used for
     * position-specific parameters (not globals), and is derived from a
     * PositionSpecificParameters object.  This is used in conjunction with the
     * m_coefficients_vector object (it is the weighted sum of those
     * coefficients).
     *
     * @see m_coefficients_vector
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente m_position_entente;

    /**
     * This object stores the position-specific and sequence-specific sequence
     * score coefficients of the position-specific parameters.  It is a vector
     * of length (the number of sequences).
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionSpecificSequenceScoreCoefficientsVector m_coefficients_vector;

    /**
     * When using unconditional Baum-Welch, we need to store up all of the
     * position ententes and wait till the end to update the actual profile.
     * Here is where we store the position ententes while we wait.
     *
     * Only used if parameters.useUnconditionalBaumWelch is true.
     */
    vector<typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente> m_unconditional_position_ententes_vector;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  /**
   * Many of the dp methods take the scaled match distributions as optional
   * arguments.  For efficiency we keep them up-to-date here and pass them into
   * those methods.  They are necessary because (only when USE_DEL_IN_DEL_OUT
   * is defined) the Match distribution contains the M->W (DeletionOut open)
   * probability, but actually the W->W extensions need to be incorporated, and they
   * differ depending on the profile position (in particular, on its distance
   * to the end of the profile).
   */
  //vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ProbabilityType> > m_scaled_match_distributions;
  // NOTE: I've changed this to MatrixValueType because with all of those DelOut Extensions, the values can get very very teensy tiny.
  vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> > m_scaled_match_distributions;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef ALLOW_BOLTZMANN_GIBBS
  // TODO: REMOVE?  FOR DEBUGGING BALDI STUFF.
  static const bool m_baldi_be_verbose = false;//true;
  //bool m_baldi_be_verbose;

    /**
     * When using unconditional Baum-Welch, we need to store up all of the
     * position Baldi-style updates and wait till the end to update the actual
     * profile.  Here is where we store the position updates while we wait.
     *
     * Only used if parameters.useUnconditionalBaumWelch is true.
     * Only used if parameters.baldiLearningRate > 0
     */
  vector<PositionBoltzmannGibbs<ResidueType, ScoreType> > m_unconditional_baldi_position_boltzmann_gibbs_changes_vector;

    /**
     * When using Baldi-style gradient ascent, we do our changes in a
     * Boltzmann-Gibbs transformed space.  This is the current profile
     * position, transformed.
     *
     * @see m_baldi_position_boltzmann_gibbs_change
     * Only used if parameters.baldiLearningRate > 0
     */
    PositionBoltzmannGibbs<ResidueType, ScoreType> m_baldi_position_boltzmann_gibbs;

    /**
     * When using Baldi-style gradient ascent, we do our changes in a
     * Boltzmann-Gibbs transformed space.  This is the change to the
     * transformed current profile position that is induced by the Baldi update
     * step.
     *
     * @see m_baldi_position_boltzmann_gibbs
     * Only used if parameters.baldiLearningRate > 0
     */
    //typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionBoltzmannGibbsChange m_baldi_position_boltzmann_gibbs_change;
  PositionBoltzmannGibbs<ResidueType, ScoreType> m_baldi_position_boltzmann_gibbs_change;
#endif // ALLOW_BOLTZMANN_GIBBS

    /**
     * When usePriors is true, this is the prior we use for the Match-state
     * emissions.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureMatchEmissionPrior<float> m_matchEmissionPrior;

    /**
     * When usePriors is true, this is the prior we use for everything except
     * the Match-state emissions.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template DirichletMixtureGlobalPrior<float> m_globalPrior;

    /**
     * This is used to save the profile position in case we want to revert
     * after modifying it.
     */
    ProfilePosition<ResidueType, ProbabilityType> m_backupProfilePosition;

    /**
     * This is used to save the entente position because when we do the wacky
     * m_bw_inverse_scalar loop, we need to revert back to the original entente
     * after scaling it.
     */
    typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::PositionEntente m_backupEntentePosition;

    /**
     * This is used to save the global profile values in case we want to revert
     * after modifying them.
     */
    ProfileTreeRoot<ResidueType, ProbabilityType> m_backupProfileGlobals;

    /**
     * This is used to save the profile values in case we want to revert after
     * modifying them, when we are using unconditional BW.
     *
     * Only used if parameters.useUnconditionalBaumWelch is true.
     */
    ProfileTreeRoot<ResidueType, ProbabilityType> m_unconditional_backupProfile;

    /**
     * Default constructor.  See reinitialize(..).
     */  
    ProfileTrainer ();

    /**
     * Construct a profile trainer with the given profile and sequences.
     */  
    ProfileTrainer (
      ProfileType * profile,
      vector<SequenceType> const & sequences
    );

    /**
     * Construct a profile trainer with the same profile, sequences, and number
     * of sequences to use as the given ProfileTrainer reference.  Note that
     * other members are not copied.
     */  
    ProfileTrainer (
      ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> const & copy_from
    );

    /**
     * Construct a profile trainer with the given profile and sequences, and
     * the number of sequences to use (use the first num_sequences_to_use
     * sequences only).
     */  
    ProfileTrainer (
      ProfileType * profile,
      vector<SequenceType> const & sequences,
      uint32_t const & num_sequences_to_use // use only the first X sequences...
    );

    /**
     * Reset everything to initial values, except for the profile, the
     * sequences, and the number of sequences to use.
     */  
    void
    reinitialize ();

    /**
     * (Re)set the profile and reinitialize.
     */  
    void
    reinitialize (
      ProfileType * profile
    );

    /**
     * Store the given parameters, and starting profile, and set the trainer to
     * its initial values.
     *
     * Note that (for now) the profile must be the same length as the original
     * profile.
     */
    void
    restart (
      const Parameters & parameters,
      ProfileType * profile
    );
  
    ScoreType
    train ();

    // Updated to return the percent change of the log score (when we are using
    // LogProbability-type scores).
    double
    percentChange ( ScoreType const & new_val, ScoreType const & old_val );

    // TODO: move this down below with the other instantiations..
    // Return the x value at which y is maximized, if ( x1, y1 ), ( x2, y2 ), and (
    // x3, y3 ) are points on a parabola.
    //
    double
    maximizeThreePointQuadratic (
      double x1,
      ScoreType y1,
      double x2,
      ScoreType y2,
      double x3,
      ScoreType y3
    )
    {
      // They can't all be 0.
      assert( !( ( x1 == 0 ) && ( x2 == 0 ) && ( x3 == 0 ) ) );
      // Everything must be non-negative.  The ys surely are (since they are ScoreTypes).
      assert( ( x1 >= 0 ) && ( x2 >= 0 ) && ( x3 >= 0 ) );
      // The middle value must be the highest.
      assert( y2 >= y1 );
      assert( y2 >= y3 );
      // The x values must all be different and must be in order.
      assert( x1 != x2 );
      if( x1 < x2 ) {
        assert( x2 < x3 );
      } else { // x1 > x2
        assert( x2 > x3 );
      }
    
      if( ( y1 == y2 ) && ( y2 == y3 ) ) {
        // They're all the same.  This is a plateau.  Tha max is y1 == y2
        // == y3, but since we're conservative we'll return the least x.
        return // min( x1, x2, x3 ):
          ( ( x1 < x2 ) ? ( ( x1 < x3 ) ? x1 : x3 ) : ( ( x2 < x3 ) ? x2 : x3 ) );
      }

      // from http://www.chemeng.ed.ac.uk/people/jack/MSO/section5/maths/part3/handout1/index.html
      ScoreType s21 = ( y2 - y1 );
      bool s21_is_negative = false;
      if( x2 > x1 ) {
        s21 /= ( x2 - x1 );
      } else { // x1 > x2
        s21 /= ( x1 - x2 );
        s21_is_negative = true;
      }
      ScoreType s32 = ( y2 - y3 );
      bool s32_is_negative = true;
      if( x3 > x2 ) {
        s32 /= ( x3 - x2 );
      } else { // x2 > x3
        s32 /= ( x2 - x3 );
        s32_is_negative = false;
      }
      assert( s21_is_negative != s32_is_negative );
      ScoreType s32_minus_s21 = ( s32 + s21 ); // adding because of the difference in signs
      int8_t s32_minus_s21_sign = ( s32_is_negative ? -1 : 1 );
      double s21_fraction_as_double = toDouble( s21 / s32_minus_s21 );
      int8_t sign_fix = ( s21_is_negative ? -1 : 1 ) * s32_minus_s21_sign;
      double return_value = ( ( ( x1 + x2 ) / 2 ) - ( ( ( x3 - x1 ) / 2 ) * sign_fix * s21_fraction_as_double ) );
      if( x1 < x2 ) {
        assert( return_value >= ( ( x1 + x2 ) / 2 ) );
        assert( return_value <= ( ( x2 + x3 ) / 2 ) );
      } else {
        assert( return_value <= ( ( x1 + x2 ) / 2 ) );
        assert( return_value >= ( ( x2 + x3 ) / 2 ) );
      }
      return( return_value );
    } // maximizeThreePointQuadratic(..)


  }; // End class ProfileTrainer

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::ProfileTrainer::Parameters ////
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
  Parameters ()
  {
    if( DEFAULT_debug >= DEBUG_All ) {
      cout << "[debug] ProfileTrainer::Parameters::<init>()" << endl;
    } // End if DEBUG_All
    this->resetToDefaults();
  } // galosh::ProfileTrainer::Parameters::<init>()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileTrainer::Parameters::<init>( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  typename ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters &
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
      // Copy constructor/operator
      operator= (
        AnyParameters const& copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTrainer::Parameters::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileTrainer::Parameters::copyFromNonVirtual( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtualDontDelegate( copy_from );
      } // copyFromNonVirtual( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
  copyFromNonVirtualDontDelegate (
    AnyParameters const & copy_from
  )
  {
    trainProfileGlobals = copy_from.trainProfileGlobals;
    trainProfilePositions = copy_from.trainProfilePositions;
    trainGlobalsFirst = copy_from.trainGlobalsFirst;
    minIterations = copy_from.minIterations;
    maxIterations = copy_from.maxIterations;
    maxPositionCycles = copy_from.maxPositionCycles;
    maxPositionCycles_sequence_identifiers = copy_from.maxPositionCycles_sequence_identifiers;
    maxPositionCycles_globals = copy_from.maxPositionCycles_globals;
    maxBaumWelchInverseScalar = copy_from.maxBaumWelchInverseScalar;
    minBaumWelchInverseScalar = copy_from.minBaumWelchInverseScalar;
    baumWelchInverseScalarIncrement = copy_from.baumWelchInverseScalarIncrement;
    maxBaumWelchInverseScalar_sequence_identifiers = copy_from.maxBaumWelchInverseScalar_sequence_identifiers;
    minBaumWelchInverseScalar_sequence_identifiers = copy_from.minBaumWelchInverseScalar_sequence_identifiers;
    baumWelchInverseScalarIncrement_sequence_identifiers = copy_from.baumWelchInverseScalarIncrement_sequence_identifiers;
    maxBaumWelchInverseScalar_globals = copy_from.maxBaumWelchInverseScalar_globals;
    minBaumWelchInverseScalar_globals = copy_from.minBaumWelchInverseScalar_globals;
    baumWelchInverseScalarIncrement_globals = copy_from.baumWelchInverseScalarIncrement_globals;
    scorePercentChangeMinimum_iteration = copy_from.scorePercentChangeMinimum_iteration;
    scorePercentChangeMinimum_position_cycle = copy_from.scorePercentChangeMinimum_position_cycle;
    euclideanDistanceMinimum_iteration = copy_from.euclideanDistanceMinimum_iteration;
    euclideanDistanceMinimum_position_cycle = copy_from.euclideanDistanceMinimum_position_cycle;
    alwaysAccept = copy_from.alwaysAccept;
    useAlignmentProfiles = copy_from.useAlignmentProfiles;
    proposeProfileLengthChanges = copy_from.proposeProfileLengthChanges;
    proposeDeletingThreshold = copy_from.proposeDeletingThreshold;
    proposeInsertingThreshold = copy_from.proposeInsertingThreshold;
    proposeDeletingThreshold_increment = copy_from.proposeDeletingThreshold_increment;
    proposeInsertingThreshold_increment = copy_from.proposeInsertingThreshold_increment;
    proposeInsertingPreAlignThreshold = copy_from.proposeInsertingPreAlignThreshold;
    proposeInsertingPostAlignThreshold = copy_from.proposeInsertingPostAlignThreshold;
    proposeInsertingOccupancyThreshold = copy_from.proposeInsertingOccupancyThreshold;
    useSensitiveThresholding = copy_from.useSensitiveThresholding;
    increaseThresholdsForLengthChanges_startIteration = copy_from.increaseThresholdsForLengthChanges_startIteration;
    increaseThresholdsForLengthChanges_minIncrement = copy_from.increaseThresholdsForLengthChanges_minIncrement;
    alwaysAccept_disallowThreshold_profileDistance_iteration = copy_from.alwaysAccept_disallowThreshold_profileDistance_iteration;
    numIterationsBetweenLengthChanges = copy_from.numIterationsBetweenLengthChanges;
    useUnconditionalBaumWelch = copy_from.useUnconditionalBaumWelch;
    baldiLearningRate = copy_from.baldiLearningRate;
    baldiTemperature = copy_from.baldiTemperature;
    profileValueMinimum = copy_from.profileValueMinimum;
    usePriors = copy_from.usePriors;
    matchEmissionPrior = copy_from.matchEmissionPrior;
    globalPrior = copy_from.globalPrior;
    positionShouldBeTrained = copy_from.positionShouldBeTrained;
    siegelMaxFindingThePeakAttempts_positions = copy_from.siegelMaxFindingThePeakAttempts_positions;
    baldiHybrid = copy_from.baldiHybrid;
    siegelEpsilonScaleFactor = copy_from.siegelEpsilonScaleFactor;
    unconditionalIsolatesGlobals = copy_from.unconditionalIsolatesGlobals;
    siegelMaxRefiningThePeakSteps_positions = copy_from.siegelMaxRefiningThePeakSteps_positions;
    siegelRefiningThePeakStepsConvergenceThreshold = copy_from.siegelRefiningThePeakStepsConvergenceThreshold;
    siegelMaxFindingTheGradientAttempts_positions = copy_from.siegelMaxFindingTheGradientAttempts_positions;
    siegelMinEpsilon = copy_from.siegelMinEpsilon;
  } // copyFromNonVirtualDontDelegate( AnyParameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_COPY
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
      copyFrom ( const Parameters & copy_from )
      {
        copyFromNonVirtual( copy_from );
      } // copyFrom( Parameters const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
      resetToDefaults ()
      {
        DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::resetToDefaults();
        // TODO: Why isn't the compiler finding "debug" in galosh::Parameters?
        //if( debug >= DEBUG_All ) {
        //  cout << "[debug] ProfileTrainer::Parameters::resetToDefaults()" << endl;
        //} // End if DEBUG_All
        trainProfileGlobals = DEFAULT_trainProfileGlobals;
        trainProfilePositions = DEFAULT_trainProfilePositions;
        trainGlobalsFirst = DEFAULT_trainGlobalsFirst;
        minIterations = DEFAULT_minIterations;
        maxIterations = DEFAULT_maxIterations;
        maxPositionCycles = DEFAULT_maxPositionCycles;
        maxPositionCycles_sequence_identifiers = DEFAULT_maxPositionCycles_sequence_identifiers;
        maxPositionCycles_globals = DEFAULT_maxPositionCycles_globals;
        maxBaumWelchInverseScalar = DEFAULT_maxBaumWelchInverseScalar;
        minBaumWelchInverseScalar = DEFAULT_minBaumWelchInverseScalar;
        baumWelchInverseScalarIncrement = DEFAULT_baumWelchInverseScalarIncrement;
        maxBaumWelchInverseScalar_sequence_identifiers = DEFAULT_maxBaumWelchInverseScalar_sequence_identifiers;
        minBaumWelchInverseScalar_sequence_identifiers = DEFAULT_minBaumWelchInverseScalar_sequence_identifiers;
        baumWelchInverseScalarIncrement_sequence_identifiers = DEFAULT_baumWelchInverseScalarIncrement_sequence_identifiers;
        maxBaumWelchInverseScalar_globals = DEFAULT_maxBaumWelchInverseScalar_globals;
        minBaumWelchInverseScalar_globals = DEFAULT_minBaumWelchInverseScalar_globals;
        baumWelchInverseScalarIncrement_globals = DEFAULT_baumWelchInverseScalarIncrement_globals;
        scorePercentChangeMinimum_iteration = DEFAULT_scorePercentChangeMinimum_iteration;
        scorePercentChangeMinimum_position_cycle = DEFAULT_scorePercentChangeMinimum_position_cycle;
        euclideanDistanceMinimum_iteration = DEFAULT_euclideanDistanceMinimum_iteration;
        euclideanDistanceMinimum_position_cycle = DEFAULT_euclideanDistanceMinimum_position_cycle;
        alwaysAccept = DEFAULT_alwaysAccept;
        useAlignmentProfiles = DEFAULT_useAlignmentProfiles;
        proposeProfileLengthChanges = DEFAULT_proposeProfileLengthChanges;
        proposeDeletingThreshold = DEFAULT_proposeDeletingThreshold;
        proposeInsertingThreshold = DEFAULT_proposeInsertingThreshold;
        proposeDeletingThreshold_increment = DEFAULT_proposeDeletingThreshold_increment;
        proposeInsertingThreshold_increment = DEFAULT_proposeInsertingThreshold_increment;
        proposeInsertingPreAlignThreshold = DEFAULT_proposeInsertingPreAlignThreshold;
        proposeInsertingPostAlignThreshold = DEFAULT_proposeInsertingPostAlignThreshold;
        proposeInsertingOccupancyThreshold = DEFAULT_proposeInsertingOccupancyThreshold;
        useSensitiveThresholding = DEFAULT_useSensitiveThresholding;
        increaseThresholdsForLengthChanges_startIteration = DEFAULT_increaseThresholdsForLengthChanges_startIteration;
        increaseThresholdsForLengthChanges_minIncrement = DEFAULT_increaseThresholdsForLengthChanges_minIncrement;
        alwaysAccept_disallowThreshold_profileDistance_iteration = DEFAULT_alwaysAccept_disallowThreshold_profileDistance_iteration;
        numIterationsBetweenLengthChanges = DEFAULT_numIterationsBetweenLengthChanges;
        useUnconditionalBaumWelch = DEFAULT_useUnconditionalBaumWelch;
        baldiLearningRate = DEFAULT_baldiLearningRate;
        baldiTemperature = DEFAULT_baldiTemperature;
        profileValueMinimum = DEFAULT_profileValueMinimum;
        usePriors = DEFAULT_usePriors;
        positionShouldBeTrained = DEFAULT_positionShouldBeTrained;
        siegelMaxFindingThePeakAttempts_positions = DEFAULT_siegelMaxFindingThePeakAttempts_positions;
        baldiHybrid = DEFAULT_baldiHybrid;
        siegelEpsilonScaleFactor = DEFAULT_siegelEpsilonScaleFactor;
        unconditionalIsolatesGlobals = DEFAULT_unconditionalIsolatesGlobals;
        siegelMaxRefiningThePeakSteps_positions = DEFAULT_siegelMaxRefiningThePeakSteps_positions;
        siegelRefiningThePeakStepsConvergenceThreshold = DEFAULT_siegelRefiningThePeakStepsConvergenceThreshold;
        siegelMaxFindingTheGradientAttempts_positions = DEFAULT_siegelMaxFindingTheGradientAttempts_positions;
        siegelMinEpsilon = DEFAULT_siegelMinEpsilon;
      } // resetToDefaults()
    
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
      template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
      void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::Parameters::
      writeParameters ( 
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        DynamicProgramming<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::Parameters::writeParameters( os );
        os << endl;

        os << "[ProfileTrainer]" << endl;
        os << "trainProfileGlobals = " << trainProfileGlobals << endl;
        os << "trainProfilePositions = " << trainProfilePositions << endl;
        os << "trainGlobalsFirst = " << trainGlobalsFirst << endl;
        os << "minIterations = " << minIterations << endl;
        os << "maxIterations = " << maxIterations << endl;
        os << "maxPositionCycles = " << maxPositionCycles << endl;
        os << "maxPositionCycles_sequence_identifiers = " << maxPositionCycles_sequence_identifiers << endl;
        os << "maxPositionCycles_globals = " << maxPositionCycles_globals << endl;
        os << "maxBaumWelchInverseScalar = " << maxBaumWelchInverseScalar << endl;
        os << "minBaumWelchInverseScalar = " << minBaumWelchInverseScalar << endl;
        os << "baumWelchInverseScalarIncrement = " << baumWelchInverseScalarIncrement << endl;
        os << "maxBaumWelchInverseScalar_sequence_identifiers = " << maxBaumWelchInverseScalar_sequence_identifiers << endl;
        os << "minBaumWelchInverseScalar_sequence_identifiers = " << minBaumWelchInverseScalar_sequence_identifiers << endl;
        os << "baumWelchInverseScalarIncrement_sequence_identifiers = " << baumWelchInverseScalarIncrement_sequence_identifiers << endl;
        os << "maxBaumWelchInverseScalar_globals = " << maxBaumWelchInverseScalar_globals << endl;
        os << "minBaumWelchInverseScalar_globals = " << minBaumWelchInverseScalar_globals << endl;
        os << "baumWelchInverseScalarIncrement_globals = " << baumWelchInverseScalarIncrement_globals << endl;
        os << "scorePercentChangeMinimum_iteration = " << scorePercentChangeMinimum_iteration << endl;
        os << "scorePercentChangeMinimum_position_cycle = " << scorePercentChangeMinimum_position_cycle << endl;
        os << "euclideanDistanceMinimum_iteration = " << euclideanDistanceMinimum_iteration << endl;
        os << "euclideanDistanceMinimum_position_cycle = " << euclideanDistanceMinimum_position_cycle << endl;
        os << "alwaysAccept = " << alwaysAccept << endl;
        os << "useAlignmentProfiles = " << useAlignmentProfiles << endl;
        os << "proposeProfileLengthChanges = " << proposeProfileLengthChanges << endl;
        os << "proposeDeletingThreshold = " << proposeDeletingThreshold << endl;
        os << "proposeInsertingThreshold = " << proposeInsertingThreshold << endl;
        os << "proposeDeletingThreshold_increment = " << proposeDeletingThreshold_increment << endl;
        os << "proposeInsertingThreshold_increment = " << proposeInsertingThreshold_increment << endl;
        os << "proposeInsertingPreAlignThreshold = " << proposeInsertingPreAlignThreshold << endl;
        os << "proposeInsertingPostAlignThreshold = " << proposeInsertingPostAlignThreshold << endl;
        os << "proposeInsertingOccupancyThreshold = " << proposeInsertingOccupancyThreshold << endl;
        os << "useSensitiveThresholding = " << useSensitiveThresholding << endl;
        os << "increaseThresholdsForLengthChanges_startIteration = " << increaseThresholdsForLengthChanges_startIteration << endl;
        os << "increaseThresholdsForLengthChanges_minIncrement = " << increaseThresholdsForLengthChanges_minIncrement << endl;
        os << "alwaysAccept_disallowThreshold_profileDistance_iteration = " << alwaysAccept_disallowThreshold_profileDistance_iteration << endl;
        os << "numIterationsBetweenLengthChanges = " << numIterationsBetweenLengthChanges << endl;
        os << "useUnconditionalBaumWelch = " << useUnconditionalBaumWelch << endl << endl;
        os << "baldiLearningRate = " << baldiLearningRate << endl;
        os << "baldiTemperature = " << baldiTemperature << endl;
        os << "profileValueMinimum = " << profileValueMinimum << endl;
        os << "usePriors = " << usePriors << endl;
        os << "positionShouldBeTrained = " << positionShouldBeTrained << endl;
        os << "siegelMaxFindingThePeakAttempts_positions = " << siegelMaxFindingThePeakAttempts_positions << endl;
        os << "baldiHybrid = " << baldiHybrid << endl;
        os << "siegelEpsilonScaleFactor = " << siegelEpsilonScaleFactor << endl;
        os << "unconditionalIsolatesGlobals = " << unconditionalIsolatesGlobals << endl;
        os << "siegelMaxRefiningThePeakSteps_positions = " << siegelMaxRefiningThePeakSteps_positions << endl;
        os << "siegelRefiningThePeakStepsConvergenceThreshold = " << siegelRefiningThePeakStepsConvergenceThreshold << endl;
        os << "siegelMaxFindingTheGradientAttempts_positions = " << siegelMaxFindingTheGradientAttempts_positions << endl;
        os << "siegelMinEpsilon = " << siegelMinEpsilon << endl;
      } // writeParameters( basic_ostream & ) const

  ////// Class galosh::ProfileTrainer::ParametersModifierTemplate ////
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( base_parameters_modifier_t::parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTrainer::ParametersModifierTemplate::<init>()" << endl;
        } // End if DEBUG_All
        isModified_reset();
      } // <init>()

      // Copy constructor
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTrainer::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )
    
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
      // Copy constructor/operator
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType> &
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      operator= (
        AnyParametersModifierTemplate const& copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProfileTrainer::ParametersModifierTemplate::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParametersModifierTemplate const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
      void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );
        isModified_trainProfileGlobals = copy_from.isModified_trainProfileGlobals;
        isModified_trainProfilePositions = copy_from.isModified_trainProfilePositions;
        isModified_trainGlobalsFirst = copy_from.isModified_trainGlobalsFirst;
        isModified_minIterations = copy_from.isModified_minIterations;
        isModified_maxIterations = copy_from.isModified_maxIterations;
        isModified_maxPositionCycles = copy_from.isModified_maxPositionCycles;
        isModified_maxPositionCycles_sequence_identifiers = copy_from.isModified_maxPositionCycles_sequence_identifiers;
        isModified_maxPositionCycles_globals = copy_from.isModified_maxPositionCycles_globals;
        isModified_maxBaumWelchInverseScalar = copy_from.isModified_maxBaumWelchInverseScalar;
        isModified_minBaumWelchInverseScalar = copy_from.isModified_minBaumWelchInverseScalar;
        isModified_baumWelchInverseScalarIncrement = copy_from.isModified_baumWelchInverseScalarIncrement;
        isModified_maxBaumWelchInverseScalar_sequence_identifiers = copy_from.isModified_maxBaumWelchInverseScalar_sequence_identifiers;
        isModified_minBaumWelchInverseScalar_sequence_identifiers = copy_from.isModified_minBaumWelchInverseScalar_sequence_identifiers;
        isModified_baumWelchInverseScalarIncrement_sequence_identifiers = copy_from.isModified_baumWelchInverseScalarIncrement_sequence_identifiers;
        isModified_maxBaumWelchInverseScalar_globals = copy_from.isModified_maxBaumWelchInverseScalar_globals;
        isModified_minBaumWelchInverseScalar_globals = copy_from.isModified_minBaumWelchInverseScalar_globals;
        isModified_baumWelchInverseScalarIncrement_globals = copy_from.isModified_baumWelchInverseScalarIncrement_globals;
        isModified_scorePercentChangeMinimum_iteration = copy_from.isModified_scorePercentChangeMinimum_iteration;
        isModified_scorePercentChangeMinimum_position_cycle = copy_from.isModified_scorePercentChangeMinimum_position_cycle;
        isModified_euclideanDistanceMinimum_iteration = copy_from.isModified_euclideanDistanceMinimum_iteration;
        isModified_euclideanDistanceMinimum_position_cycle = copy_from.isModified_euclideanDistanceMinimum_position_cycle;
        isModified_alwaysAccept = copy_from.isModified_alwaysAccept;
        isModified_useAlignmentProfiles = copy_from.isModified_useAlignmentProfiles;
        isModified_proposeProfileLengthChanges = copy_from.isModified_proposeProfileLengthChanges;
        isModified_proposeDeletingThreshold = copy_from.isModified_proposeDeletingThreshold;
        isModified_proposeInsertingThreshold = copy_from.isModified_proposeInsertingThreshold;
        isModified_proposeDeletingThreshold_increment = copy_from.isModified_proposeDeletingThreshold_increment;
        isModified_proposeInsertingThreshold_increment = copy_from.isModified_proposeInsertingThreshold_increment;
        isModified_proposeInsertingPreAlignThreshold = copy_from.isModified_proposeInsertingPreAlignThreshold;
        isModified_proposeInsertingPostAlignThreshold = copy_from.isModified_proposeInsertingPostAlignThreshold;
        isModified_proposeInsertingOccupancyThreshold = copy_from.isModified_proposeInsertingOccupancyThreshold;
        isModified_useSensitiveThresholding = copy_from.isModified_useSensitiveThresholding;
        isModified_increaseThresholdsForLengthChanges_startIteration = copy_from.isModified_increaseThresholdsForLengthChanges_startIteration;
        isModified_increaseThresholdsForLengthChanges_minIncrement = copy_from.isModified_increaseThresholdsForLengthChanges_minIncrement;
        isModified_alwaysAccept_disallowThreshold_profileDistance_iteration = copy_from.isModified_alwaysAccept_disallowThreshold_profileDistance_iteration;
        isModified_numIterationsBetweenLengthChanges = copy_from.isModified_numIterationsBetweenLengthChanges;
        isModified_useUnconditionalBaumWelch = copy_from.isModified_useUnconditionalBaumWelch;
        isModified_baldiLearningRate = copy_from.isModified_baldiLearningRate;
        isModified_baldiTemperature = copy_from.isModified_baldiTemperature;
        isModified_profileValueMinimum = copy_from.isModified_profileValueMinimum;
        isModified_usePriors = copy_from.isModified_usePriors;
        isModified_positionShouldBeTrained = copy_from.isModified_positionShouldBeTrained;
        isModified_siegelMaxFindingThePeakAttempts_positions = copy_from.isModified_siegelMaxFindingThePeakAttempts_positions;
        isModified_baldiHybrid = copy_from.isModified_baldiHybrid;
        isModified_siegelEpsilonScaleFactor = copy_from.isModified_siegelEpsilonScaleFactor;
        isModified_unconditionalIsolatesGlobals = copy_from.isModified_unconditionalIsolatesGlobals;
        isModified_siegelMaxRefiningThePeakSteps_positions = copy_from.isModified_siegelMaxRefiningThePeakSteps_positions;
        isModified_siegelRefiningThePeakStepsConvergenceThreshold = copy_from.isModified_siegelRefiningThePeakStepsConvergenceThreshold;
        isModified_siegelMaxFindingTheGradientAttempts_positions = copy_from.isModified_siegelMaxFindingTheGradientAttempts_positions;
        isModified_siegelMinEpsilon = copy_from.isModified_siegelMinEpsilon;
      } // isModified_copyFromNonVirtual ( AnyParametersModifierTemplate const & )


  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
  reset ()
      {
        isModified_reset();
        base_parameters_modifier_t::parameters.resetToDefaults();
      } // reset()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
      void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        base_parameters_modifier_t::isModified_reset();

        isModified_trainProfileGlobals = false;
        isModified_trainProfilePositions = false;
        isModified_trainGlobalsFirst = false;
        isModified_minIterations = false;
        isModified_maxIterations = false;
        isModified_maxPositionCycles = false;
        isModified_maxPositionCycles_sequence_identifiers = false;
        isModified_maxPositionCycles_globals = false;
        isModified_maxBaumWelchInverseScalar = false;
        isModified_minBaumWelchInverseScalar = false;
        isModified_baumWelchInverseScalarIncrement = false;
        isModified_maxBaumWelchInverseScalar_sequence_identifiers = false;
        isModified_minBaumWelchInverseScalar_sequence_identifiers = false;
        isModified_baumWelchInverseScalarIncrement_sequence_identifiers = false;
        isModified_maxBaumWelchInverseScalar_globals = false;
        isModified_minBaumWelchInverseScalar_globals = false;
        isModified_baumWelchInverseScalarIncrement_globals = false;
        isModified_scorePercentChangeMinimum_iteration = false;
        isModified_scorePercentChangeMinimum_position_cycle = false;
        isModified_euclideanDistanceMinimum_iteration = false;
        isModified_euclideanDistanceMinimum_position_cycle = false;
        isModified_alwaysAccept = false;
        isModified_useAlignmentProfiles = false;
        isModified_proposeProfileLengthChanges = false;
        isModified_proposeDeletingThreshold = false;
        isModified_proposeInsertingThreshold = false;
        isModified_proposeDeletingThreshold_increment = false;
        isModified_proposeInsertingThreshold_increment = false;
        isModified_proposeInsertingPreAlignThreshold = false;
        isModified_proposeInsertingPostAlignThreshold = false;
        isModified_proposeInsertingOccupancyThreshold = false;
        isModified_useSensitiveThresholding = false;
        isModified_increaseThresholdsForLengthChanges_startIteration = false;
        isModified_increaseThresholdsForLengthChanges_minIncrement = false;
        isModified_alwaysAccept_disallowThreshold_profileDistance_iteration = false;
        isModified_numIterationsBetweenLengthChanges = false;
        isModified_useUnconditionalBaumWelch = false;
        isModified_baldiLearningRate = false;
        isModified_baldiTemperature = false;
        isModified_profileValueMinimum = false;
        isModified_usePriors = false;
        isModified_positionShouldBeTrained = false;
        isModified_siegelMaxFindingThePeakAttempts_positions = false;
        isModified_baldiHybrid = false;
        isModified_siegelEpsilonScaleFactor = false;
        isModified_unconditionalIsolatesGlobals = false;
        isModified_siegelMaxRefiningThePeakSteps_positions = false;
        isModified_siegelRefiningThePeakStepsConvergenceThreshold = false;
        isModified_siegelMaxFindingTheGradientAttempts_positions = false;
        isModified_siegelMinEpsilon = false;
      } // isModified_reset()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
      void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      )
      {
        DynamicProgramming<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType>::writeParametersModifier( os );
        os << endl;

        os << "[ProfileTrainer]" << endl;
        if( isModified_trainProfileGlobals ) {
          os << "trainProfileGlobals = " << base_parameters_modifier_t::parameters.trainProfileGlobals << endl;
        }
        if( isModified_trainProfilePositions ) {
          os << "trainProfilePositions = " << base_parameters_modifier_t::parameters.trainProfilePositions << endl;
        }
        if( isModified_trainGlobalsFirst ) {
          os << "trainGlobalsFirst = " << base_parameters_modifier_t::parameters.trainGlobalsFirst << endl;
        }
        if( isModified_minIterations ) {
          os << "minIterations = " << base_parameters_modifier_t::parameters.minIterations << endl;
        }
        if( isModified_maxIterations ) {
          os << "maxIterations = " << base_parameters_modifier_t::parameters.maxIterations << endl;
        }
        if( isModified_maxPositionCycles ) {
          os << "maxPositionCycles = " << base_parameters_modifier_t::parameters.maxPositionCycles << endl;
        }
        if( isModified_maxPositionCycles_sequence_identifiers ) {
          os << "maxPositionCycles_sequence_identifiers = " << base_parameters_modifier_t::parameters.maxPositionCycles_sequence_identifiers << endl;
        }
        if( isModified_maxPositionCycles_globals ) {
          os << "maxPositionCycles_globals = " << base_parameters_modifier_t::parameters.maxPositionCycles_globals << endl;
        }
        if( isModified_maxBaumWelchInverseScalar ) {
          os << "maxBaumWelchInverseScalar = " << base_parameters_modifier_t::parameters.maxBaumWelchInverseScalar << endl;
        }
        if( isModified_minBaumWelchInverseScalar ) {
          os << "minBaumWelchInverseScalar = " << base_parameters_modifier_t::parameters.minBaumWelchInverseScalar << endl;
        }
        if( isModified_baumWelchInverseScalarIncrement ) {
          os << "baumWelchInverseScalarIncrement = " << base_parameters_modifier_t::parameters.baumWelchInverseScalarIncrement << endl;
        }
        if( isModified_maxBaumWelchInverseScalar_sequence_identifiers ) {
          os << "maxBaumWelchInverseScalar_sequence_identifiers = " << base_parameters_modifier_t::parameters.maxBaumWelchInverseScalar_sequence_identifiers << endl;
        }
        if( isModified_minBaumWelchInverseScalar_sequence_identifiers ) {
          os << "minBaumWelchInverseScalar_sequence_identifiers = " << base_parameters_modifier_t::parameters.minBaumWelchInverseScalar_sequence_identifiers << endl;
        }
        if( isModified_baumWelchInverseScalarIncrement_sequence_identifiers ) {
          os << "baumWelchInverseScalarIncrement_sequence_identifiers = " << base_parameters_modifier_t::parameters.baumWelchInverseScalarIncrement_sequence_identifiers << endl;
        }
        if( isModified_maxBaumWelchInverseScalar_globals ) {
          os << "maxBaumWelchInverseScalar_globals = " << base_parameters_modifier_t::parameters.maxBaumWelchInverseScalar_globals << endl;
        }
        if( isModified_minBaumWelchInverseScalar_globals ) {
          os << "minBaumWelchInverseScalar_globals = " << base_parameters_modifier_t::parameters.minBaumWelchInverseScalar_globals << endl;
        }
        if( isModified_baumWelchInverseScalarIncrement_globals ) {
          os << "baumWelchInverseScalarIncrement_globals = " << base_parameters_modifier_t::parameters.baumWelchInverseScalarIncrement_globals << endl;
        }
        if( isModified_scorePercentChangeMinimum_iteration ) {
          os << "scorePercentChangeMinimum_iteration = " << base_parameters_modifier_t::parameters.scorePercentChangeMinimum_iteration << endl;
        }
        if( isModified_scorePercentChangeMinimum_position_cycle ) {
          os << "scorePercentChangeMinimum_position_cycle = " << base_parameters_modifier_t::parameters.scorePercentChangeMinimum_position_cycle << endl;
        }
        if( isModified_euclideanDistanceMinimum_iteration ) {
          os << "euclideanDistanceMinimum_iteration = " << base_parameters_modifier_t::parameters.euclideanDistanceMinimum_iteration << endl;
        }
        if( isModified_euclideanDistanceMinimum_position_cycle ) {
          os << "euclideanDistanceMinimum_position_cycle = " << base_parameters_modifier_t::parameters.euclideanDistanceMinimum_position_cycle << endl;
        }
        if( isModified_alwaysAccept ) {
          os << "alwaysAccept = " << base_parameters_modifier_t::parameters.alwaysAccept << endl;
        }
        if( isModified_useAlignmentProfiles ) {
          os << "useAlignmentProfiles = " << base_parameters_modifier_t::parameters.useAlignmentProfiles << endl;
        }
        if( isModified_proposeProfileLengthChanges ) {
          os << "proposeProfileLengthChanges = " << base_parameters_modifier_t::parameters.proposeProfileLengthChanges << endl;
        }
        if( isModified_proposeDeletingThreshold ) {
          os << "proposeDeletingThreshold = " << base_parameters_modifier_t::parameters.proposeDeletingThreshold << endl;
        }
        if( isModified_proposeInsertingThreshold ) {
          os << "proposeInsertingThreshold = " << base_parameters_modifier_t::parameters.proposeInsertingThreshold << endl;
        }
        if( isModified_proposeDeletingThreshold_increment ) {
          os << "proposeDeletingThreshold_increment = " << base_parameters_modifier_t::parameters.proposeDeletingThreshold_increment << endl;
        }
        if( isModified_proposeInsertingThreshold_increment ) {
          os << "proposeInsertingThreshold_increment = " << base_parameters_modifier_t::parameters.proposeInsertingThreshold_increment << endl;
        }
        if( isModified_proposeInsertingPreAlignThreshold ) {
          os << "proposeInsertingPreAlignThreshold = " << base_parameters_modifier_t::parameters.proposeInsertingPreAlignThreshold << endl;
        }
        if( isModified_proposeInsertingPostAlignThreshold ) {
          os << "proposeInsertingPostAlignThreshold = " << base_parameters_modifier_t::parameters.proposeInsertingPostAlignThreshold << endl;
        }
        if( isModified_proposeInsertingOccupancyThreshold ) {
          os << "proposeInsertingOccupancyThreshold = " << base_parameters_modifier_t::parameters.proposeInsertingOccupancyThreshold << endl;
        }
        if( isModified_useSensitiveThresholding ) {
          os << "useSensitiveThresholding = " << base_parameters_modifier_t::parameters.useSensitiveThresholding << endl;
        }
        if( isModified_increaseThresholdsForLengthChanges_startIteration ) {
          os << "increaseThresholdsForLengthChanges_startIteration = " << base_parameters_modifier_t::parameters.increaseThresholdsForLengthChanges_startIteration << endl;
        }
        if( isModified_increaseThresholdsForLengthChanges_minIncrement ) {
          os << "increaseThresholdsForLengthChanges_minIncrement = " << base_parameters_modifier_t::parameters.increaseThresholdsForLengthChanges_minIncrement << endl;
        }
        if( isModified_alwaysAccept_disallowThreshold_profileDistance_iteration ) {
          os << "alwaysAccept_disallowThreshold_profileDistance_iteration = " << base_parameters_modifier_t::parameters.alwaysAccept_disallowThreshold_profileDistance_iteration << endl;
        }
        if( isModified_numIterationsBetweenLengthChanges ) {
          os << "numIterationsBetweenLengthChanges = " << base_parameters_modifier_t::parameters.numIterationsBetweenLengthChanges << endl;
        }
        if( isModified_useUnconditionalBaumWelch ) {
          os << "useUnconditionalBaumWelch = " << base_parameters_modifier_t::parameters.useUnconditionalBaumWelch << endl;
        }
        if( isModified_baldiLearningRate ) {
          os << "baldiLearningRate = " << base_parameters_modifier_t::parameters.baldiLearningRate << endl;
        }
        if( isModified_baldiTemperature ) {
          os << "baldiTemperature = " << base_parameters_modifier_t::parameters.baldiTemperature << endl;
        }
        if( isModified_profileValueMinimum ) {
          os << "profileValueMinimum = " << base_parameters_modifier_t::parameters.profileValueMinimum << endl;
        }
        if( isModified_usePriors ) {
          os << "usePriors = " << base_parameters_modifier_t::parameters.usePriors << endl;
        }
        if( isModified_positionShouldBeTrained ) {
          os << "positionShouldBeTrained = " << base_parameters_modifier_t::parameters.positionShouldBeTrained << endl;
        }
        if( isModified_siegelMaxFindingThePeakAttempts_positions ) {
          os << "siegelMaxFindingThePeakAttempts_positions = " << base_parameters_modifier_t::parameters.siegelMaxFindingThePeakAttempts_positions << endl;
        }
        if( isModified_baldiHybrid ) {
          os << "baldiHybrid = " << base_parameters_modifier_t::parameters.baldiHybrid << endl;
        }
        if( isModified_siegelEpsilonScaleFactor ) {
          os << "siegelEpsilonScaleFactor = " << base_parameters_modifier_t::parameters.siegelEpsilonScaleFactor << endl;
        }
        if( isModified_unconditionalIsolatesGlobals ) {
          os << "unconditionalIsolatesGlobals = " << base_parameters_modifier_t::parameters.unconditionalIsolatesGlobals << endl;
        }
        if( isModified_siegelMaxRefiningThePeakSteps_positions ) {
          os << "siegelMaxRefiningThePeakSteps_positions = " << base_parameters_modifier_t::parameters.siegelMaxRefiningThePeakSteps_positions << endl;
        }
        if( isModified_siegelRefiningThePeakStepsConvergenceThreshold ) {
          os << "siegelRefiningThePeakStepsConvergenceThreshold = " << base_parameters_modifier_t::parameters.siegelRefiningThePeakStepsConvergenceThreshold << endl;
        }
        if( isModified_siegelMaxFindingTheGradientAttempts_positions ) {
          os << "siegelMaxFindingTheGradientAttempts_positions = " << base_parameters_modifier_t::parameters.siegelMaxFindingTheGradientAttempts_positions << endl;
        }
        if( isModified_siegelMinEpsilon ) {
          os << "siegelMinEpsilon = " << base_parameters_modifier_t::parameters.siegelMinEpsilon << endl;
        }
      } // writeParametersModifier ( basic_ostream & ) const

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  template <class ParametersType>
  template<class AnyParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
      void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParameters & target_parameters )
      {
        base_parameters_modifier_t::applyModifications( target_parameters );

        if( isModified_trainProfileGlobals ) {
          target_parameters.trainProfileGlobals =
            base_parameters_modifier_t::parameters.trainProfileGlobals;
        }
        if( isModified_trainProfilePositions ) {
          target_parameters.trainProfilePositions =
            base_parameters_modifier_t::parameters.trainProfilePositions;
        }
        if( isModified_trainGlobalsFirst ) {
          target_parameters.trainGlobalsFirst =
            base_parameters_modifier_t::parameters.trainGlobalsFirst;
        }
        if( isModified_minIterations ) {
          target_parameters.minIterations =
            base_parameters_modifier_t::parameters.minIterations;
        }
        if( isModified_maxIterations ) {
          target_parameters.maxIterations =
            base_parameters_modifier_t::parameters.maxIterations;
        }
        if( isModified_maxPositionCycles ) {
          target_parameters.maxPositionCycles =
            base_parameters_modifier_t::parameters.maxPositionCycles;
        }
        if( isModified_maxPositionCycles_sequence_identifiers ) {
          target_parameters.maxPositionCycles_sequence_identifiers =
            base_parameters_modifier_t::parameters.maxPositionCycles_sequence_identifiers;
        }
        if( isModified_maxPositionCycles_globals ) {
          target_parameters.maxPositionCycles_globals =
            base_parameters_modifier_t::parameters.maxPositionCycles_globals;
        }
        if( isModified_maxBaumWelchInverseScalar ) {
          target_parameters.maxBaumWelchInverseScalar =
            base_parameters_modifier_t::parameters.maxBaumWelchInverseScalar;
        }
        if( isModified_minBaumWelchInverseScalar ) {
          target_parameters.minBaumWelchInverseScalar =
            base_parameters_modifier_t::parameters.minBaumWelchInverseScalar;
        }
        if( isModified_baumWelchInverseScalarIncrement ) {
          target_parameters.baumWelchInverseScalarIncrement =
            base_parameters_modifier_t::parameters.baumWelchInverseScalarIncrement;
        }
        if( isModified_maxBaumWelchInverseScalar_sequence_identifiers ) {
          target_parameters.maxBaumWelchInverseScalar_sequence_identifiers =
            base_parameters_modifier_t::parameters.maxBaumWelchInverseScalar_sequence_identifiers;
        }
        if( isModified_minBaumWelchInverseScalar_sequence_identifiers ) {
          target_parameters.minBaumWelchInverseScalar_sequence_identifiers =
            base_parameters_modifier_t::parameters.minBaumWelchInverseScalar_sequence_identifiers;
        }
        if( isModified_baumWelchInverseScalarIncrement_sequence_identifiers ) {
          target_parameters.baumWelchInverseScalarIncrement_sequence_identifiers =
            base_parameters_modifier_t::parameters.baumWelchInverseScalarIncrement_sequence_identifiers;
        }
        if( isModified_maxBaumWelchInverseScalar_globals ) {
          target_parameters.maxBaumWelchInverseScalar_globals =
            base_parameters_modifier_t::parameters.maxBaumWelchInverseScalar_globals;
        }
        if( isModified_minBaumWelchInverseScalar_globals ) {
          target_parameters.minBaumWelchInverseScalar_globals =
            base_parameters_modifier_t::parameters.minBaumWelchInverseScalar_globals;
        }
        if( isModified_baumWelchInverseScalarIncrement_globals ) {
          target_parameters.baumWelchInverseScalarIncrement_globals =
            base_parameters_modifier_t::parameters.baumWelchInverseScalarIncrement_globals;
        }
        if( isModified_scorePercentChangeMinimum_iteration ) {
          target_parameters.scorePercentChangeMinimum_iteration =
            base_parameters_modifier_t::parameters.scorePercentChangeMinimum_iteration;
        }
        if( isModified_scorePercentChangeMinimum_position_cycle ) {
          target_parameters.scorePercentChangeMinimum_position_cycle =
            base_parameters_modifier_t::parameters.scorePercentChangeMinimum_position_cycle;
        }
        if( isModified_euclideanDistanceMinimum_iteration ) {
          target_parameters.euclideanDistanceMinimum_iteration =
            base_parameters_modifier_t::parameters.euclideanDistanceMinimum_iteration;
        }
        if( isModified_euclideanDistanceMinimum_position_cycle ) {
          target_parameters.euclideanDistanceMinimum_position_cycle =
            base_parameters_modifier_t::parameters.euclideanDistanceMinimum_position_cycle;
        }
        if( isModified_alwaysAccept ) {
          target_parameters.alwaysAccept =
            base_parameters_modifier_t::parameters.alwaysAccept;
        }
        if( isModified_useAlignmentProfiles ) {
          target_parameters.useAlignmentProfiles =
            base_parameters_modifier_t::parameters.useAlignmentProfiles;
        }
        if( isModified_proposeProfileLengthChanges ) {
          target_parameters.proposeProfileLengthChanges =
            base_parameters_modifier_t::parameters.proposeProfileLengthChanges;
        }
        if( isModified_proposeDeletingThreshold ) {
          target_parameters.proposeDeletingThreshold =
            base_parameters_modifier_t::parameters.proposeDeletingThreshold;
        }
        if( isModified_proposeInsertingThreshold ) {
          target_parameters.proposeInsertingThreshold =
            base_parameters_modifier_t::parameters.proposeInsertingThreshold;
        }
        if( isModified_proposeDeletingThreshold_increment ) {
          target_parameters.proposeDeletingThreshold_increment =
            base_parameters_modifier_t::parameters.proposeDeletingThreshold_increment;
        }
        if( isModified_proposeInsertingThreshold_increment ) {
          target_parameters.proposeInsertingThreshold_increment =
            base_parameters_modifier_t::parameters.proposeInsertingThreshold_increment;
        }
        if( isModified_proposeInsertingPreAlignThreshold ) {
          target_parameters.proposeInsertingPreAlignThreshold =
            base_parameters_modifier_t::parameters.proposeInsertingPreAlignThreshold;
        }
        if( isModified_proposeInsertingPostAlignThreshold ) {
          target_parameters.proposeInsertingPostAlignThreshold =
            base_parameters_modifier_t::parameters.proposeInsertingPostAlignThreshold;
        }
        if( isModified_proposeInsertingOccupancyThreshold ) {
          target_parameters.proposeInsertingOccupancyThreshold =
            base_parameters_modifier_t::parameters.proposeInsertingOccupancyThreshold;
        }
        if( isModified_useSensitiveThresholding ) {
          target_parameters.useSensitiveThresholding =
            base_parameters_modifier_t::parameters.useSensitiveThresholding;
        }
        if( isModified_increaseThresholdsForLengthChanges_startIteration ) {
          target_parameters.increaseThresholdsForLengthChanges_startIteration =
            base_parameters_modifier_t::parameters.increaseThresholdsForLengthChanges_startIteration;
        }
        if( isModified_increaseThresholdsForLengthChanges_minIncrement ) {
          target_parameters.increaseThresholdsForLengthChanges_minIncrement =
            base_parameters_modifier_t::parameters.increaseThresholdsForLengthChanges_minIncrement;
        }
        if( isModified_alwaysAccept_disallowThreshold_profileDistance_iteration ) {
          target_parameters.alwaysAccept_disallowThreshold_profileDistance_iteration =
            base_parameters_modifier_t::parameters.alwaysAccept_disallowThreshold_profileDistance_iteration;
        }
        if( isModified_numIterationsBetweenLengthChanges ) {
          target_parameters.numIterationsBetweenLengthChanges =
            base_parameters_modifier_t::parameters.numIterationsBetweenLengthChanges;
        }
        if( isModified_useUnconditionalBaumWelch ) {
          target_parameters.useUnconditionalBaumWelch =
            base_parameters_modifier_t::parameters.useUnconditionalBaumWelch;
        }
        if( isModified_baldiLearningRate ) {
          target_parameters.baldiLearningRate =
            base_parameters_modifier_t::parameters.baldiLearningRate;
        }
        if( isModified_baldiTemperature ) {
          target_parameters.baldiTemperature =
            base_parameters_modifier_t::parameters.baldiTemperature;
        }
        if( isModified_profileValueMinimum ) {
          target_parameters.profileValueMinimum =
            base_parameters_modifier_t::parameters.profileValueMinimum;
        }
        if( isModified_usePriors ) {
          target_parameters.usePriors =
            base_parameters_modifier_t::parameters.usePriors;
        }
        if( isModified_positionShouldBeTrained ) {
          target_parameters.positionShouldBeTrained =
            base_parameters_modifier_t::parameters.positionShouldBeTrained;
        }
        if( isModified_siegelMaxFindingThePeakAttempts_positions ) {
          target_parameters.siegelMaxFindingThePeakAttempts_positions =
            base_parameters_modifier_t::parameters.siegelMaxFindingThePeakAttempts_positions;
        }
        if( isModified_baldiHybrid ) {
          target_parameters.baldiHybrid =
            base_parameters_modifier_t::parameters.baldiHybrid;
        }
        if( isModified_siegelEpsilonScaleFactor ) {
          target_parameters.siegelEpsilonScaleFactor =
            base_parameters_modifier_t::parameters.siegelEpsilonScaleFactor;
        }
        if( isModified_unconditionalIsolatesGlobals ) {
          target_parameters.unconditionalIsolatesGlobals =
            base_parameters_modifier_t::parameters.unconditionalIsolatesGlobals;
        }
        if( isModified_siegelMaxRefiningThePeakSteps_positions ) {
          target_parameters.siegelMaxRefiningThePeakSteps_positions =
            base_parameters_modifier_t::parameters.siegelMaxRefiningThePeakSteps_positions;
        }
        if( isModified_siegelRefiningThePeakStepsConvergenceThreshold ) {
          target_parameters.siegelRefiningThePeakStepsConvergenceThreshold =
            base_parameters_modifier_t::parameters.siegelRefiningThePeakStepsConvergenceThreshold;
        }
        if( isModified_siegelMaxFindingTheGradientAttempts_positions ) {
          target_parameters.siegelMaxFindingTheGradientAttempts_positions =
            base_parameters_modifier_t::parameters.siegelMaxFindingTheGradientAttempts_positions;
        }
        if( isModified_siegelMinEpsilon ) {
          target_parameters.siegelMinEpsilon =
            base_parameters_modifier_t::parameters.siegelMinEpsilon;
        }
      } // applyModifications( AnyParameters & )

  ////// Class galosh::ProfileTrainer ////
  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    /**
     * Default constructor.  See reinitialize(..).
     */  
    ProfileTrainer (
    ) :
      m_dynamic_programming(),
      m_profile( 0 ),
      m_sequences( 0 ),
      m_error( this ),
      m_sequence_count( 0 )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::<init>()" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    /**
     * Construct a profile trainer with the given profile and sequences.
     */  
    ProfileTrainer (
      ProfileType * profile,
      vector<SequenceType> const & sequences
    ) :
      m_dynamic_programming(),
      m_profile( profile ),
      m_sequences( sequences ),
      m_error( this ),
      m_sequence_count( sequences.size() )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::<init>( ProfileType, vector<SequenceType> )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( ProfileType *, vector<SequenceType> const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    /**
     * Construct a profile trainer with the given profile and sequences, and
     * the number of sequences to use (use the first num_sequences_to_use
     * sequences only).
     */  
    ProfileTrainer (
      ProfileType * profile,
      vector<SequenceType> const & sequences,
      uint32_t const & num_sequences_to_use // use only the first X sequences...
    ) :
      m_dynamic_programming(),
      m_profile( profile ),
      m_sequences( sequences ),
      m_error( this ),
      m_sequence_count( min( ( size_t )( ( num_sequences_to_use > 0 ) ? num_sequences_to_use : sequences.size() ), sequences.size() ) )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::<init>( ProfileType, vector<SequenceType>, " << num_sequences_to_use << " )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( ProfileType *, vector<SequenceType> const &, uint32_t const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_INIT
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    /**
     * Construct a profile trainer with the same profile, sequences, and number
     * of sequences to use as the given ProfileTrainer reference.  Note that
     * other members are not copied.
     */  
    ProfileTrainer (
      ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType> const & copy_from
    ) :
      m_dynamic_programming(),
      m_profile( copy_from.m_profile ),
      m_sequences( copy_from.m_sequences ),
      m_error( this ),
      m_sequence_count( copy_from.m_sequence_count )
  {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::<init>( ProfileTrainer const & )" << endl;
      } // End if DEBUG_All
      // Do nothing else
    } // <init>( ProfileTrainer const & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    /**
     * (Re)set the profile and reinitialize.
     */  
    reinitialize (
      ProfileType * profile
    )
    {
      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::reinitialize( ProfileType )" << endl;
      } // End if DEBUG_All

      m_profile = profile;
      reinitialize();
    } // reinitialize( ProfileType * )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    reinitialize (
    )
    {
      if( m_trainingParameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::reinitialize()" << endl;
      } // End if DEBUG_All

      m_globalsToRevertTo.reinitialize();
      m_startingProfile_iteration.reinitialize();
      m_startingProfile_position_cycle.reinitialize();

#ifndef NDEBUG
      // TODO :REMOVE 
      m_forward_matrices.reinitialize();
#endif // !NDEBUG

      m_forward_rows_1.reinitialize();
      m_forward_rows_2.reinitialize();

      m_unconditional_last_forward_rows.reinitialize();
      m_unconditional_last_forward_rows.reinitialize();
      m_unconditional_second_to_last_forward_rows.reinitialize();

      m_anchor_columns.reinitialize();
      m_anchor_rows.reinitialize();

      m_backward_rows_1.reinitialize();
      m_backward_rows_2.reinitialize();

      m_backupProfilePosition.reinitialize();

      m_position_entente.reinitialize();
      m_global_entente.reinitialize();

      m_coefficients_vector.reinitialize();

      m_unconditional_position_ententes_vector.resize( 0 );

#ifdef ALLOW_BOLTZMANN_GIBBS
      m_unconditional_baldi_position_boltzmann_gibbs_changes_vector.resize( 0 );
      m_baldi_position_boltzmann_gibbs.reinitialize();
      m_baldi_position_boltzmann_gibbs_change.reinitialize();
#endif // ALLOW_BOLTZMANN_GIBBS

      m_matchEmissionPrior.reinitialize();
      m_globalPrior.reinitialize();

      return;
    } // reinitialize()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    /**
     * Store the given parameters, and starting profile, and set the trainer to
     * its initial values.
     *
     * Note that (for now) the profile must be the same length as the original
     * profile.
     */
    restart (
      const Parameters & parameters,
      ProfileType * profile
    )
    {
      if( parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::restart( parameters, profile )" << endl;
      } // End if DEBUG_All
      m_trainingParameters = parameters;

      // TODO: REMOVE
      //cout << "IN RESTART (1): PROFILE IS " << *m_profile << endl;
      m_profile = profile;
      // Make sure each distribution's probabilities sum to 1!
      m_profile->normalize( 0 ); // TODO: Should we use m_trainingParameters.profileValueMinimum ?
      // TODO: Why is this necessary?!
      m_profile->ensurePositionsKnowTheirRoot();

      m_globalsToRevertTo.copyExceptPositions( *m_profile );
      // TODO: REMOVE! MAGIC #
      m_globalsToRevertTo.multiplyByExceptPositions( 10.0 );
#ifndef DISALLOW_FLANKING_TRANSITIONS
      // TODO: REMOVE! TESTING.  FLOOB!
      // Make an even stronger pull toward the initial flanking transitions.
      //m_globalsToRevertTo[ Transition::fromPreAlign ] *= 10.0;
      //m_globalsToRevertTo[ Transition::fromPostAlign ] *= 10.0;
#endif // !DISALLOW_FLANKING_TRANSITIONS
      // TODO: REMOVE
      //cout << "IN RESTART (2): PROFILE IS " << *m_profile << endl;

      m_startingProfile_iteration.reinitialize( *m_profile );
      // TODO: REMOVE
      //cout << "IN RESTART (2.1): PROFILE IS " << *m_profile << endl;
      m_startingProfile_position_cycle.reinitialize( *m_profile );
      // TODO: REMOVE
      //cout << "IN RESTART (2.2): PROFILE IS " << *m_profile << endl;

      // TODO: Add a parameter to determine whether we use anchor rows/cols or
      // use the full forward matrices.

#ifndef NDEBUG
      // TODO :REMOVE 
      m_forward_matrices.reinitialize(
        m_profile->length(),
        m_sequences,
        m_sequence_count
      );
#endif // !NDEBUG

      m_forward_rows_1.reinitialize(
        m_sequences,
        m_sequence_count
      );
      
      m_forward_rows_2.reinitialize(
        m_sequences,
        m_sequence_count
      );

      if( m_trainingParameters.useUnconditionalBaumWelch ) {
        m_unconditional_last_forward_rows.reinitialize(
          m_sequences,
          m_sequence_count
        );
        m_unconditional_last_forward_rows.reinitialize(
          m_sequences,
          m_sequence_count
        );
        m_unconditional_second_to_last_forward_rows.reinitialize(
          m_sequences,
          m_sequence_count
        );
      } // End if useUnconditionalBaumWelch
      
      m_forward_rows_ptr = &m_forward_rows_1;
      m_prev_forward_rows_ptr = &m_forward_rows_2;

      // TODO: MAGIC # (store_every_Nth_row)
      // TODO: Make this a parameter
      uint32_t store_every_nth_row = 1;//5;
      // TODO: REMOVE!
      //cout << "calling m_anchor_rows.reinitialize(..)" << endl;
      m_anchor_rows.reinitialize(
        m_profile->length(),
        m_sequences,
        m_sequence_count,
        numeric_limits<uint32_t>::max(),
        store_every_nth_row
      );
      // TODO: REMOVE!
      //cout << "RETURNED from m_anchor_rows.reinitialize(..)" << endl;

      // TODO: MAGIC # (store_every_Nth_column)
      // TODO: Make this a parameter
      uint32_t store_every_nth_column =
        ( ( m_anchor_rows.m_storeEveryNthRow == 1 ) ? 0 : 5 );
      if(
        ( store_every_nth_column > 0 )
      ) {
        m_anchor_columns.reinitialize(
          m_profile->length(),
          m_sequences,
          m_sequence_count,
          store_every_nth_column,
          ( ( store_every_nth_column == 0 ) ? 0 : numeric_limits<uint32_t>::max() )
        );
      }

      m_backward_rows_1.reinitialize(
        m_sequences,
        m_sequence_count
      );
      
      m_backward_rows_2.reinitialize(
        m_sequences,
        m_sequence_count
      );
      
      m_backward_rows_ptr = &m_backward_rows_1;
      m_next_backward_rows_ptr = &m_backward_rows_2;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      m_scaled_match_distributions.resize( m_profile->length() - 1 );
      for( uint32_t tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
        m_profile->createScaledMatchDistributionForPosition(
          tmp_pos_i,
          m_scaled_match_distributions[ tmp_pos_i ]
        );
        // TODO: REMOVE
        //cout << "scaled match distribution for pos " << tmp_pos_i << ": " << m_scaled_match_distributions[ tmp_pos_i ] << endl;
        // TODO: REMOVE
        assert( m_scaled_match_distributions[ tmp_pos_i ][ TransitionFromMatch::toMatch ] > 0.0 );
      } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      // Set the profile ptr of the backup profile position.
      m_backupProfilePosition.setProfileTreeRoot( m_profile->getProfileTreeRoot() );

      // Start out training position-specific params.  Globals after.  Or the
      // other way 'round, depending on the parameters.
      m_trainingPhase =
        ( m_trainingParameters.trainGlobalsFirst ?
          TRAINING_PHASE_Globals :
          TRAINING_PHASE_Positions );
      //m_trainingGlobals = m_trainingParameters.trainGlobalsFirst;

      m_position_entente.reinitialize();
      m_global_entente.reinitialize();

      m_coefficients_vector.reinitialize(
        m_sequence_count
      );

      if( m_trainingParameters.useUnconditionalBaumWelch ) {
        m_unconditional_position_ententes_vector.resize( m_profile->length() ); // One for each row except the first
      } // End if useUnconditionalBaumWelch

#ifdef ALLOW_BOLTZMANN_GIBBS      
      if( m_trainingParameters.baldiLearningRate > 0 ) {
        if( m_trainingParameters.useUnconditionalBaumWelch ) {
          m_unconditional_baldi_position_boltzmann_gibbs_changes_vector.resize( m_profile->length() ); // One for each row except the first
        } // End if useUnconditionalBaumWelch

        m_baldi_position_boltzmann_gibbs.reinitialize(
          m_trainingParameters.baldiTemperature
        );
        m_baldi_position_boltzmann_gibbs_change.reinitialize(
          m_trainingParameters.baldiTemperature
        );
      } // End if baldiLearningRate > 0
#endif //ALLOW_BOLTZMANN_GIBBS

      if( m_trainingParameters.usePriors ) {
        if( m_trainingParameters.matchEmissionPrior != NULL ) {
          m_matchEmissionPrior =
            *( m_trainingParameters.matchEmissionPrior );
        } else {
          m_matchEmissionPrior.reinitializeToLaplace();
        }

        // TODO: REMOVE
        //cout << "MatchEmission prior is " << m_matchEmissionPrior << endl;

        if( m_trainingParameters.globalPrior != NULL ) {
          m_globalPrior =
            *( m_trainingParameters.globalPrior );
        } else {
          m_globalPrior.reinitializeToLaplace();
        }

        // TODO: REMOVE
        //cout << "Global prior is " << m_globalPrior << endl;

      } // End if usePriors

      return;
    } // restart( const Parameters &, const ProfileType & )

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_TRAIN
    ScoreType
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    train ()
    {

      if( m_parameters.debug >= DEBUG_All ) {
        cout << "[debug] ProfileTrainer::train()" << endl;
      } // End if DEBUG_All
      // TODO: REMOVE
      //cout << "BEFORE RESTARTING, PROFILE IS " << *m_profile << endl;
      // TODO: REMOVE
      //cout << "calling restart( .. )" << endl;
      restart( m_parameters, m_profile );
      // TODO: REMOVE
      //cout << "RETURNED from restart( .. )" << endl;

      uint32_t last_row = m_profile->length();
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator anchor_columns_iter;
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator anchor_rows_iter;
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator tmp_matrices_iter;
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::SequentialAccessContainer::iterator tmp_matrices_iter2;

      // TODO: Move these to be member values (m_sequence_scores, etc)
      vector<ScoreType> sequence_scores( m_sequence_count );
      // TODO: Move these to be member values (m_sequence_scores, etc)
      vector<ScoreType> backup_sequence_scores( m_sequence_count );
      // For scaling the ententes to avoid overflow
      ScoreType largest_sequence_score;

      bool forward_rows_are_already_rotated_and_calculated = false;
      bool backward_rows_are_already_rotated_and_calculated = false;

#ifdef ALLOW_BOLTZMANN_GIBBS
      register uint8_t temp_index;
      register ScoreType score_after_modification;
      register double epsilon;
      register double epsilon_scale_value;
      register double epsilon_scale_value_scale_factor; // TODO: Incorporate into Unconditional BS.
      register double baldi_change_maxed_out_epsilon; // TODO: REMOVE.  TESTING.
      register double baldi_change_had_an_effect; // TODO: REMOVE.  TESTING.
      ProfilePosition<ResidueType, ProbabilityType> baldi_last_profile_position; // TODO: REMOVE.  TESTING.
      register uint32_t finding_the_peak_attempts;
      register uint32_t refining_the_peak_steps;
      vector<uint8_t> sample_indices( 3 );
      vector<ScoreType> sample_scores( 3 );
      vector<double> sample_epsilons( 3 );
#endif // ALLOW_BOLTZMANN_GIBBS


      // TODO: REMOVE.  TESTING
      // If doing UBW, only do positions on even iters and only globals on odd iters?
      static const bool ubw_cbw_hybrid = m_trainingParameters.unconditionalIsolatesGlobals && m_trainingParameters.useUnconditionalBaumWelch;

      // TODO: REMOVE
      //cout << "PROFILE IS " << *m_profile << endl;

      bool swap =
        m_dynamic_programming.forward_score(
          m_trainingParameters,
          false, // don't use viterbi
          *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          m_sequences,
          m_sequence_count,
          &m_anchor_columns,
          &m_anchor_rows,
          *m_prev_forward_rows_ptr,
          *m_forward_rows_ptr,
          &sequence_scores,
          &largest_sequence_score,
          m_scoreBeforeTraining
        );
      m_endingScore = m_scoreBeforeTraining;
      if( swap ) {
        // Rotate the forward_rows
        m_temp_forward_rows_ptr = m_forward_rows_ptr;
        m_forward_rows_ptr = m_prev_forward_rows_ptr;
        m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
      }
      forward_rows_are_already_rotated_and_calculated = true;

      // TODO: REMOVE
      //cout << "Ending score, calculated the new way, is " << m_endingScore << endl;

      // TODO: REMOVE?
      assert( !isnan( m_endingScore ) );

      if( m_endingScore == 0 ) {
        std::cerr << "WARNING: Underflow detected (forward_score() returns 0).  You should try setting the MatrixValueType to bfloat, which is the type with the largest range available." << endl;
        assert( m_endingScore != 0 );
      }

      //cout << "final row for seq 0, calculated the new way, is " << ( *m_forward_rows_ptr )[ 0 ] << endl;
      //cout << "second-to-last row for seq 0, calculated the new way, is " << ( *m_prev_forward_rows_ptr )[ 0 ] << endl;

      //m_scoreBeforeTraining = m_endingScore =
      //  m_dynamic_programming.forward_score( 
      //    m_trainingParameters,
      //    *m_profile,
      //    m_sequences,
      //    m_sequence_count,
      //    m_forward_matrices,
      //    &sequence_scores,
      //    &largest_sequence_score
      //  );
      //cout << "Ending score, calculated the OLD way, is " << m_endingScore << endl;
      //cout << "final row for seq 0, calculated the OLD way, is " << m_forward_matrices[ last_row ][ 0 ] << endl;
      //cout << "second-to-last row for seq 0, calculated the OLD way, is " << m_forward_matrices[ last_row - 1 ][ 0 ] << endl;

      if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
        cout << "\nThe total forward score before training is "
             << m_scoreBeforeTraining << endl;
        cout << "Training using ";
        if( m_trainingParameters.useUnconditionalBaumWelch ) {
          cout << "Unconditional ";
          if( ubw_cbw_hybrid ) {
            cout << "(Hybrid) ";
          }
        } else {
          cout << "Conditional ";
        }
#ifdef ALLOW_BOLTZMANN_GIBBS
        if( m_trainingParameters.baldiLearningRate > 0 ) {
          cout << "Baldi";
          if( m_trainingParameters.siegelMaxFindingThePeakAttempts_positions > 0 ) {
            cout << " / Siegel, using siegelEpsilonScaleFactor " << m_trainingParameters.siegelEpsilonScaleFactor;
          } else {
            cout << ", using baldiLearningRate " << m_trainingParameters.baldiLearningRate;
          }
          if( m_trainingParameters.baldiHybrid ) {
            cout << ", hybrid with Baum-Welch";
          }
        } else {
          cout << "Baum-Welch";
        }
#else // !ALLOW_BOLTZMANN_GIBBS
        cout << "Baum-Welch";
#endif ALLOW_BOLTZMANN_GIBBS
        if( m_trainingParameters.proposeProfileLengthChanges ) {
          cout << ", using Dynamic Model Surgery with initial ";
          if(
             ( m_trainingParameters.proposeInsertingThreshold == m_trainingParameters.proposeInsertingPreAlignThreshold ) &&
             ( m_trainingParameters.proposeInsertingThreshold == m_trainingParameters.proposeInsertingPostAlignThreshold )
          ) {
            if( m_trainingParameters.proposeInsertingThreshold == m_trainingParameters.proposeDeletingThreshold ) {
              cout << "threshold " << m_trainingParameters.proposeInsertingThreshold << " ";
            } else {
              cout << "thresholds (ins=" << m_trainingParameters.proposeInsertingThreshold << ",del=" << m_trainingParameters.proposeDeletingThreshold << ") ";
            }
          } else {
            cout << "thresholds (ins=" << m_trainingParameters.proposeInsertingThreshold << ",pre=" << m_trainingParameters.proposeInsertingPreAlignThreshold << ",post=" << m_trainingParameters.proposeInsertingPostAlignThreshold << ",del=" << m_trainingParameters.proposeDeletingThreshold << ") ";
          }
          if( m_trainingParameters.proposeInsertingThreshold_increment == m_trainingParameters.proposeDeletingThreshold_increment ) {
            cout << "and threshold increment " << m_trainingParameters.proposeInsertingThreshold_increment;
          } else {
            cout << "and threshold increments (ins=" << m_trainingParameters.proposeInsertingThreshold_increment << ",del=" << m_trainingParameters.proposeDeletingThreshold_increment;
          }
        } // End if proposeProfileLengthChanges
        cout << "." << endl;
       } // End if verbose

      // Local vars
      ScoreType sequence_score = 0;

      ScoreType total_sequence_score = 0;
      ScoreType scaled_total_sequence_score = 0;
      
      ProbabilityType tmp_prob = 0;

      double profile_distance_iteration = 0.0;
      double profile_distance_position_cycle = 0.0;

      // Alignment profiles (and profile length changes) stuff

      // TODO: Why are pre-align insertions happening so much,
      // especially when the profile is too *long*? (TODO: Still?)

      // 
      //bool increment_pre_align_insertion_threshold_with_nonconsecutive_uses = true;
      //bool increment_post_align_insertion_threshold_with_nonconsecutive_uses = true;

      static const bool count_seqs_exceeding_alignment_profile_thresholds = false;

      const double fraction_seqs_exceeding_insertion_threshold_insert_threshold = m_trainingParameters.proposeInsertingOccupancyThreshold;
      const double fraction_seqs_exceeding_deletion_threshold_delete_threshold = ( 1.0 - fraction_seqs_exceeding_insertion_threshold_insert_threshold );
      //double fraction_seqs_exceeding_deletion_threshold_maybe_delete_threshold = .5;//.25;
      //double fraction_seqs_exceeding_insertion_threshold_maybe_insert_threshold = .5;//.25;

      double propose_inserting_threshold = m_trainingParameters.proposeInsertingThreshold;
      double propose_inserting_prealign_threshold = m_trainingParameters.proposeInsertingPreAlignThreshold;
      double propose_inserting_postalign_threshold = m_trainingParameters.proposeInsertingPostAlignThreshold;
      double propose_deleting_threshold = m_trainingParameters.proposeDeletingThreshold;

      double propose_inserting_threshold_increment =
        m_trainingParameters.proposeInsertingThreshold_increment;
      double propose_deleting_threshold_increment =
        m_trainingParameters.proposeDeletingThreshold_increment;

      // What, if anything, to print to cout
      const bool cout_profile_length_changes = ( m_trainingParameters.verbosity >= VERBOSITY_High ); //false;//true;
      const bool cout_indel_fractions = ( m_trainingParameters.verbosity >= VERBOSITY_All ); //false;//true;
      static const bool cout_profile_length_changes_seqs = false;
      static const bool cout_indel_fractions_seqs = false;

      // TODO: REMOVE?  For debugging...
      static const bool cout_anchor_rows_iter_changes = false;
      uint32_t anchor_rows_iter_i = 0; // Used only if cout_anchor_rows_iter_changes is true.
      static const bool cout_anchor_columns_iter_changes = false;
      uint32_t anchor_columns_iter_i = 0; // Used only if cout_anchor_columns_iter_changes is true.

      static const uint32_t min_profile_length = 2; // magic #, sort of (our algorithms can't handle profiles smaller than size 2).
      uint32_t max_profile_length = 0; // see below
      if( m_trainingParameters.proposeProfileLengthChanges ) {
        // for now, make max_profile_length twice the largest sequence size.
        for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
          if( max_profile_length < m_sequences[ m_seq_i ].length() ) {
            max_profile_length = m_sequences[ m_seq_i ].length();
          }
        } // End foreach seq.
        // TODO: MAGIC # (max_profile_length multiple of longest seq len)
        max_profile_length *= 2;
        if( cout_profile_length_changes ) {
          // TODO: REMOVE?
          cout << "Allowing profile length changes in the range " << min_profile_length << " to " << max_profile_length << "." << endl;
        } // End if cout_profile_length_changes
      } // End if m_trainingParameters.proposeProfileLengthChanges


      /// TODO: Also replace the global entente with a single AlignmentProfilePosition. Don't forget to account for the tricksy storage of pre-align and post-align insertion opens and extensions.
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfilePosition alignment_profile_position;
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfilePosition alignment_profile_position_allseqs;
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfilePosition lastpos_alignment_profile_position_allseqs;

      // Profile length changes stuff
      uint32_t num_iterations_until_next_length_change_proposal =
        ( m_trainingParameters.numIterationsBetweenLengthChanges + 1 );
      bool just_grew_profile = false;
      bool just_shrunk_profile = false;
      bool just_caught_preAlign_cycle = false;
      bool last_length_change_was_insertion = false;
      bool last_length_change_was_deletion = false;
      bool prev_to_last_length_change_was_insertion = false;
      bool profile_length_changed_this_iteration = false;
      bool profile_length_changed_last_iteration = false;
      bool prev_profile_length_changed_this_iteration = false;
      bool cycle_detected_this_iteration = false;
      bool cycle_detected_last_iteration = false;
      bool enable_reverting = true;
      bool revert_globals_for_changed_profile_length = false;
      bool just_reverted_globals_for_changed_profile_length = false;
      bool globals_are_at_starting_values = true;
      bool cycle_detected_deletions = false;
      bool cycle_detected_insertions = false;
      bool fancy_cycle_detected_deletions = false;
      bool fancy_cycle_detected_insertions = false;
      double deletion_fraction = 0.0;
      double insertion_fraction = 0.0;
      double last_pos_deletion_fraction = 0.0;
      double last_pos_insertion_fraction = 0.0;
      double largest_insertion_fraction = 0.0;
      double largest_deletion_fraction = 0.0;
      uint32_t largest_insertion_fraction_row = 0;
      uint32_t largest_deletion_fraction_row = 0;
      uint32_t pos_i = 0;
      uint32_t rpos_i = 0;  // measured in reverse: 0 is ( last_row - 1 )...
      uint32_t num_seqs_exceeding_last_pos_insertion_threshold = 0;
      uint32_t num_seqs_exceeding_last_pos_deletion_threshold = 0;
      uint32_t num_seqs_exceeding_deletion_threshold = 0;
      uint32_t num_seqs_exceeding_insertion_threshold = 0;
      bool delete_last_pos = false;
      uint32_t propose_inserting_pos_i = 0;
      uint32_t propose_deleting_pos_i = 0;
      uint32_t propose_inserting_rpos_i = 0;
      uint32_t propose_deleting_rpos_i = 0;
      // These are set to true the first time we delete/insert, and stay true.
      // It's so we know whether we can test last_deleted_pos_i, etc.  I could
      // have used a special value of those to indicate that it's invalid, but
      // instead I did it this way.
      bool have_deleted_a_pos = false;
      bool have_inserted_a_pos = false;
      uint32_t last_deleted_pos_i = 0;
      uint32_t last_inserted_pos_i = 0;
      uint32_t last_deleted_rpos_i = 0;
      uint32_t last_inserted_rpos_i = 0;
      uint32_t length_before_last_insertion = 0;
      uint32_t length_before_last_deletion = 0;
      uint32_t prev_to_last_deleted_pos_i = 0;
      uint32_t prev_to_last_inserted_pos_i = 0;
      uint32_t prev_to_last_deleted_rpos_i = 0;
      uint32_t prev_to_last_inserted_rpos_i = 0;
      uint32_t prev_length_before_last_insertion = 0;
      uint32_t prev_length_before_last_deletion = 0;
      bool insertion_fraction_exceeds_threshold = false;
      bool last_pos_insertion_fraction_exceeds_threshold = false;
      bool deletion_fraction_exceeds_threshold = false;
      bool last_pos_deletion_fraction_exceeds_threshold = false;
      uint32_t length_at_last_revert = 0;
      uint32_t length_at_prev_to_last_revert = 0;
      uint32_t length_at_last_iteration_in_which_length_changed = 0;
      // To keep track of the iteration in which each position was created:
      vector<uint32_t> lengthadjust_position_creation_iters( m_profile->length(), numeric_limits<uint32_t>::max() );
      /////// Parameters //////
      // If this is false, increase thresholds only when cycles are detected.
      // If it is true, also increase thresholds when anything is inserted or
      // deleted.
      // NOTE: in either case, there is some danger of bias being introduced
      // because the first rows processed (which are actually the
      // higher-numbered rows) have a lower threshold than the later-processed
      // rows.
      // TODO: REMOVE?
      bool increase_thresholds_for_length_changes = false;//true;
      // Eventually we should increase thresholds for length changes, in case
      // we get stuck in a long cycle.
      // TODO: MAGIC # (200)
      const uint32_t start_increasing_thresholds_for_length_changes_after_iteration = m_trainingParameters.increaseThresholdsForLengthChanges_startIteration;//500;
      const double min_threshold_increment_for_increasing_threshold_for_length_changes = m_trainingParameters.increaseThresholdsForLengthChanges_minIncrement;//1E-4;
      // TODO: PUT BACK.  TESTING.
      //uint32_t length_change_min_iteration = 1;
      // TODO: REMOVE.  TESTING.
      const uint32_t length_change_min_iteration =
 	( m_trainingParameters.trainGlobalsFirst ? 0 : 1 );
      // TODO: Make this a parameter
      static const bool ensure_even_indels_after_reverting = true;//false;

      // TODO: REMOVE? TESTING.
      bool use_sensitive_threshold_increments =
        m_trainingParameters.useSensitiveThresholding;//true;//false;

      double min_insertion_fraction_passing_threshold;
      double min_last_pos_insertion_fraction_passing_threshold;
      double min_insertion_fraction_passing_threshold_causing_a_cycle;

      double min_deletion_fraction_passing_threshold;
      double min_last_pos_deletion_fraction_passing_threshold;
      double min_deletion_fraction_passing_threshold_causing_a_cycle;

      // TODO: REMOVE.  TESTING.  UNFINISHED!  HARD because when changing the profile to the best so far, everything changes length, which means we need to update a lot of things.  Too complicated for me right now.
      static const bool revert_to_best_profile_when_disabling_reversion = false;
      ProfileType best_profile( *m_profile );
      ScoreType best_profile_score;

      // TODO: REMOVE.  TESTING.

      // If the profile_distance_position_cycle is less than this value, turn
      // disallow_alwaysAccept on.  If it is greater than or equal to this
      // value, turn it off.
      const double disallow_alwaysAccept_profile_distance_iteration_threshold = m_trainingParameters.alwaysAccept_disallowThreshold_profileDistance_iteration;//1E-5;
      bool disallow_alwaysAccept = false;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      // TODO: REMOVE? TESTING.  The right thing to do is to
      // recalculate all forward rows preceding an inserted or deleted
      // position during lengthadjust when USE_DEL_IN_DEL_OUT is true,
      // since the del-out probs depend on the distance a position is
      // to the end of the profile.  But it slows it way down and it
      // doesn't really seem necessary.  This should remain false
      // unless you are experimenting...
      static const bool do_not_recalculate_all_forward_rows_during_lengthadjust = false;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      uint32_t tmp_pos_i;
      // TODO: REMOVE!
      ScoreType tmp_score;

      // When using priors, we expect the score to get worse sometimes, so we
      // should only look at the convergence of the profile.
      // Update: actually the euclidean convergence metric works just fine, so
      // we'll save a tiny bit of work by not computing the score's percent
      // change.
      static const bool use_score_percent_change = false; //!m_trainingParameters.usePriors;
      // TODO: if not using percent_change, can we save time by not recomputing
      // the score?

      // Here's where the refining begins (above here is just setup
      // stuff)..
      if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
          ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
          m_trainingParameters.verbosity >= VERBOSITY_All ) {
        cout << "Training..";
        cout.flush();
      }

      bool just_changed_iterations = true;
      for( m_iteration = 0;
           (
            ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) < m_trainingParameters.minIterations ) ||
             (
               ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) < m_trainingParameters.maxIterations ) &&
               (
                 revert_globals_for_changed_profile_length ||
                 just_reverted_globals_for_changed_profile_length ||
                 profile_length_changed_this_iteration ||
                 cycle_detected_this_iteration ||
                 // TODO: REMOVE. TESTING. Note this is after updating m_iteration..
                 ( ubw_cbw_hybrid && ( ( m_iteration % 2 ) == 1 ) ) ? true :
                 (
                   ( 
                     !use_score_percent_change ||
                     (
                       // Note that this will only be evaluated if m_iteration
                       // >= minIterations, so it's ok that
                       // m_startingScore_iteration is undefined until after
                       // the first iteration.
                       (
                         m_scorePercentChange =
                         percentChange(
                           m_endingScore,
                           m_startingScore_iteration
                         )
                       ) >
                       m_trainingParameters.scorePercentChangeMinimum_iteration
                     )
                   ) &&
                   (
                     profile_distance_iteration >
                     m_trainingParameters.euclideanDistanceMinimum_iteration
                   )
                 )
               )
             )
           );
           m_iteration++
      ) {

        // TODO: REMOVE
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
//          // TODO: REMOVE
//          cout << "[ " << tmp_pos_i << " ] Scaled match distribution: " << m_scaled_match_distributions[ tmp_pos_i ] << endl;
//        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        // TODO: REMOVE
        assert( m_endingScore < 1 );

        // TODO: REMOVE
        //if( m_iteration == 30 ) {
        //  cout_profile_length_changes = true;
        //  cout_indel_fractions = true;
        //  cout_profile_length_changes_seqs = true;
        //  cout_indel_fractions_seqs = true;
        //}

        if( true ) {
        
          // TODO: REMOVE TEST. TESTING.
          if( !ubw_cbw_hybrid || ( ( m_iteration % 2 ) == 0 ) ) {

            // Start over for this iteration..
            m_startingScore_iteration = m_endingScore;
            // Remember the profile as it was.
            m_startingProfile_iteration.copyFrom( *m_profile );

          // TODO: REMOVE TEST. TESTING.
          }

          if(
            m_trainingParameters.proposeProfileLengthChanges
          ) {

                // TODO: REMOVE.  TESTING
                if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                  if( m_iteration > 1 ) {
                    cout << endl;
                    cout << "largest_insertion_fraction: " << largest_insertion_fraction << " at row " << largest_insertion_fraction_row << endl;
                    cout << "largest_deletion_fraction: " << largest_deletion_fraction << " at row " << largest_deletion_fraction_row << endl;
                  }
                }
                // TODO: REMOVE

                if( false && ( m_iteration > 1 ) && !just_grew_profile && !just_shrunk_profile && !just_reverted_globals_for_changed_profile_length ) {
                  // TODO: REMOVE.  REALLY ..
                  //if( largest_insertion_fraction > 0 ) {
                  //  ( ( m_row_i == 0 ) ? propose_inserting_prealign_threshold : ( ( m_row_i == last_row ) ? propose_inserting_postalign_threshold : propose_inserting_threshold ) ) = largest_insertion_fraction;
                  //}
                  //if( largest_deletion_fraction > 0 ) {
                  //  propose_deleting_threshold = largest_deletion_fraction;
                  //}
                  //use_sensitive_threshold_increments = false;
                  //enable_reverting = false;

                  // TODO: REMOVE!!!
                  if( ( largest_insertion_fraction > 0 ) && ( largest_deletion_fraction > 0 ) &&
                      // if we're to delete the position on either side of the inserted one, don't bother.
                      ( largest_insertion_fraction_row != largest_deletion_fraction_row ) &&
                      ( largest_deletion_fraction_row != ( largest_insertion_fraction_row + 1 ) )
                  ) {
                    if(
                       //( largest_insertion_fraction < propose_inserting_threshold ) &&
                       ( largest_insertion_fraction > ( propose_inserting_threshold / 2 ) ) &&
                       //( largest_deletion_fraction < propose_deleting_threshold ) &&
                       ( largest_deletion_fraction > ( propose_deleting_threshold / 2 ) )
                    ) {
                      // TODO: REMOVE ALL OF THIS
                      if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                        cout << "Swapping..." << endl;
                      }
                    //if( largest_insertion_fraction_row != 0 ) {
                      // shift them around.  eg delete the largest_deletion_fraction_row and insert after the largest_insertion_fraction_row.
                      if( largest_insertion_fraction_row < largest_deletion_fraction_row ) {
                        ( *m_profile )[ largest_insertion_fraction_row ].even(); // Start the new pos even().
                        for( tmp_pos_i = largest_insertion_fraction_row + 1; tmp_pos_i < largest_deletion_fraction_row; tmp_pos_i++ ) {
                          ( *m_profile )[ tmp_pos_i ] = m_startingProfile_position_cycle[ tmp_pos_i - 1 ];
                        }
                      } else { // largest_deletion_fraction_row < largest_insertion_fraction_row
                        for( tmp_pos_i = largest_deletion_fraction_row - 1; tmp_pos_i < largest_insertion_fraction_row - 1; tmp_pos_i++ ) {
                          ( *m_profile )[ tmp_pos_i ] = m_startingProfile_position_cycle[ tmp_pos_i + 1 ];
                        }
                        ( *m_profile )[ largest_insertion_fraction_row - 1 ].even(); // Start the new pos even().
                      } // End if largest_insertion_fraction_row < largest_deletion_fraction_row .. else ..
                      // Now recalc the score & matrices..
                      bool swap =
                        m_dynamic_programming.forward_score(
                          m_trainingParameters,
                          false, // don't use viterbi
                          *m_profile,
                #if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          m_scaled_match_distributions,
                #endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          m_sequences,
                          m_sequence_count,
                          &m_anchor_columns,
                          &m_anchor_rows,
                          *m_prev_forward_rows_ptr,
                          *m_forward_rows_ptr,
                          &sequence_scores,
                          &largest_sequence_score,
                          m_endingScore
                        );
                      if( swap ) {
                        // Rotate the forward_rows
                        m_temp_forward_rows_ptr = m_forward_rows_ptr;
                        m_forward_rows_ptr = m_prev_forward_rows_ptr;
                        m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                      }
                      forward_rows_are_already_rotated_and_calculated = true;
                    } // End if the fractions are in-range..
                  } // End if we can swap the rows..
                } // End if m_iteration > 1
                // TODO: REMOVE.  TESTING.
                largest_insertion_fraction = 0;
                largest_deletion_fraction = 0;
                largest_insertion_fraction_row = numeric_limits<uint32_t>::max();
                largest_deletion_fraction_row = numeric_limits<uint32_t>::max();



            // TODO: REMOVE?
            if(
              !increase_thresholds_for_length_changes &&
              ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) > start_increasing_thresholds_for_length_changes_after_iteration )
            ) {
              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                cout << "Henceforth the thresholds will be increased for every length change." << endl;
              }
              increase_thresholds_for_length_changes = true;
            }

            // TODO: REMOVE.  TESTING.
            if(
              revert_to_best_profile_when_disabling_reversion &&
              enable_reverting
            ) {
              if( m_endingScore < best_profile_score ) {
                // Less is better..
                best_profile.copyFrom( *m_profile );
                best_profile_score = m_endingScore;
              }
            }
          
            // TODO: REMOVE.  TESTING.
            if(
              profile_length_changed_this_iteration
            ) {
              // If the profile length "changed" to the same length as it was
              // before, then increase the thresholds.  Note that we do this
              // even if use_sensitive_threshold_increments is enabled, but only
              // when we have no cycles (and therefore no threshold increase
              // due to sensitive thresholding).
              if(
                (
                  use_sensitive_threshold_increments ?
                  !cycle_detected_this_iteration :
                  true
                ) &&
                (
                  length_at_last_iteration_in_which_length_changed ==
                  m_profile->length()
                )
              ) {
                // Length didn't change; increase thresholds.
                propose_inserting_threshold +=
                  propose_inserting_threshold_increment;
                propose_inserting_prealign_threshold +=
                  propose_inserting_threshold_increment;
                propose_inserting_postalign_threshold +=
                  propose_inserting_threshold_increment;
                propose_deleting_threshold +=
                  propose_deleting_threshold_increment;
                          
                if( cout_profile_length_changes ) {
                  cout << "Thresholds are being increased because the length changes didn't actually change the length.  propose_deleting_threshold increased to " << propose_deleting_threshold << " and propose_inserting_threshold increased to " << propose_inserting_threshold << "." << endl;
                } // End if cout_profile_length_changes
                if(
                  (
                    ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                    ( m_trainingParameters.verbosity < VERBOSITY_High )
                  ) ||
                  ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                ) {
                  if(
                    ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
                    ( propose_inserting_threshold == propose_inserting_postalign_threshold )
                  ) {
                    if( propose_inserting_threshold == propose_deleting_threshold ) {
                      cout << "threshold:" << propose_inserting_threshold;
                    } else {
                      cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
                    }
                  } else {
                    cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
                  }
                  if(
                    m_trainingParameters.verbosity >= VERBOSITY_Medium
                  ) {
                    cout << endl;
                  } else {
                    cout.flush();
                  }
                } // End if verbose

              } // End if the profile length "changed" to the same length as it was before, then increase the thresholds.
              length_at_last_iteration_in_which_length_changed =
                m_profile->length();
            } // End if profile length changed this iteration.
            if(
              cycle_detected_this_iteration &&
              (
                use_sensitive_threshold_increments ||
                !profile_length_changed_this_iteration
              )
            ) {
              // Ok there was a cycle. Increase the thresholds if we're doing
              // use_sensitive_threshold_increments or if there was no length
              // change...
              if( use_sensitive_threshold_increments && ( min_insertion_fraction_passing_threshold_causing_a_cycle < 1 ) ) {
                if(
                   propose_inserting_threshold <
                   min_insertion_fraction_passing_threshold_causing_a_cycle
                ) {
                  propose_inserting_threshold =
                    min_insertion_fraction_passing_threshold_causing_a_cycle;
                  // TODO: REMOVE
                  //cout << endl << "propose_inserting_threshold:" << propose_inserting_threshold << endl;
                } else if( !increase_thresholds_for_length_changes ) {
                  // TODO: REMOVE / FIX.  This assumes that the thresholds are in lock-step.
                  //cout << "HUH? Somehow the threshold was already increased beyond what the use_sensitive_threshold_increments value would set it to?!" << endl;
                  //assert( false );
                  //exit( 1 );
                }
                if(
                   propose_inserting_prealign_threshold <
                  min_insertion_fraction_passing_threshold_causing_a_cycle
                ) {
                  propose_inserting_prealign_threshold =
                    min_insertion_fraction_passing_threshold_causing_a_cycle;
                } else if( !increase_thresholds_for_length_changes ) {
                  // TODO: REMOVE / FIX.  This assumes that the thresholds are in lock-step.
                  //cout << "HUH? Somehow the threshold was already increased beyond what the use_sensitive_threshold_increments value would set it to?!" << endl;
                  //assert( false );
                  //exit( 1 );
                }
                if(
                  propose_inserting_postalign_threshold <
                  min_insertion_fraction_passing_threshold_causing_a_cycle
                ) {
                  propose_inserting_postalign_threshold =
                    min_insertion_fraction_passing_threshold_causing_a_cycle;
                } else if( !increase_thresholds_for_length_changes ) {
                  // TODO: REMOVE / FIX.  This assumes that the thresholds are in lock-step.
                  //cout << "HUH? Somehow the threshold was already increased beyond what the use_sensitive_threshold_increments value would set it to?!" << endl;
                  //assert( false );
                  //exit( 1 );
                }
              } // End if use_sensitive_threshold_increments && min_insertion_fraction_passing_threshold_causing_a_cycle < 1
              if( use_sensitive_threshold_increments && ( min_deletion_fraction_passing_threshold_causing_a_cycle < 1 ) ) {
                if(
                  propose_deleting_threshold <
                  min_deletion_fraction_passing_threshold_causing_a_cycle
                ) {
                  propose_deleting_threshold =
                    min_deletion_fraction_passing_threshold_causing_a_cycle;
                  // TODO: REMOVE
                  //cout << endl << "propose_deleting_threshold:" << propose_deleting_threshold << endl;
                } else if( !increase_thresholds_for_length_changes ) {
                  // TODO: REMOVE / FIX.  This assumes that the thresholds are in lock-step.
                  //cout << "HUH? Somehow the threshold was already increased beyond what the use_sensitive_threshold_increments value would set it to?!" << endl;
                  //assert( false );
                  //exit( 1 );
                }
              } // End if use_sensitive_threshold_increments && min_deletion_fraction_passing_threshold_causing_a_cycle < 1

              propose_inserting_threshold +=
                propose_inserting_threshold_increment;
              propose_inserting_prealign_threshold +=
                propose_inserting_threshold_increment;
              propose_inserting_postalign_threshold +=
                propose_inserting_threshold_increment;
              propose_deleting_threshold +=
                propose_deleting_threshold_increment;

              if( cout_profile_length_changes ) {
                if( use_sensitive_threshold_increments ) {
                  cout << "Thresholds are being increased to the lowest fraction affecting the count of sequences passing the threshold, and then increased by the threshold_increment...";
                } else {
                  cout << "Thresholds are being increased because cycles were detected last iteration, but the length didn't change.";
                }
                cout << "propose_deleting_threshold increased to " << propose_deleting_threshold << " and propose_inserting_threshold increased to " << propose_inserting_threshold << "." << endl;
              } // End if cout_profile_length_changes
              if(
                (
                  ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                  ( m_trainingParameters.verbosity < VERBOSITY_High )
                ) ||
                ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
              ) {
                if(
                  ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
                  ( propose_inserting_threshold == propose_inserting_postalign_threshold )
                ) {
                  if( propose_inserting_threshold == propose_deleting_threshold ) {
                    cout << "threshold:" << propose_inserting_threshold;
                  } else {
                    cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
                  }
                } else {
                  cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
                }
                if(
                  m_trainingParameters.verbosity >= VERBOSITY_Medium
                ) {
                  cout << endl;
                } else {
                  cout.flush();
                }
              } // End if verbose

              // If we've had this situation for two consecutive iterations,
              // reset cycle detection, to make sure we're allowing any pending
              // changes.
              if(
                cycle_detected_last_iteration &&
                !profile_length_changed_last_iteration
              ) {
                if(
                  ( have_deleted_a_pos && last_length_change_was_deletion ) ||
                  ( have_inserted_a_pos && last_length_change_was_insertion )
                ) {
                  // TODO: REMOVE
                  if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                    cout << endl << "Resetting cycle detection." << endl;
                  }
                  // reset
                  have_deleted_a_pos = false;
                  have_inserted_a_pos = false;
                  last_length_change_was_deletion = false;
                  last_length_change_was_insertion = false;
                } else
                if( !use_sensitive_threshold_increments ) {
                  // If we've already reset, then switch to sensitive
                  // thresholding.
                  use_sensitive_threshold_increments = true;
                  // TODO: REMOVE
                  if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                    cout << endl << "Sensitive thresholding enabled." << endl;
                  }
                } // End if we've already reset
              } // End if cycle_detected_last_iteration && !prev_profile_length_changed_this_iteration
            } // End if !profile_length_changed_this_iteration && cycle_detected_this_iteration
            
            // TODO: REMOVE.  TESTING.
            // If we didn't change lengths last iteration then reset cycle
            // detection stuffs.
            //if(
            //  ( num_iterations_until_next_length_change_proposal == 0 ) &&
            //  !profile_length_changed_this_iteration
            //) {
            //  // TODO: REMOVE
            //  cout << endl << "RESETTING" << endl;
            //  // reset
            //  have_deleted_a_pos = false;
            //  have_inserted_a_pos = false;
            //  last_length_change_was_deletion = false;
            //  last_length_change_was_insertion = false;
            //}

            if( num_iterations_until_next_length_change_proposal > 0 ) {
              num_iterations_until_next_length_change_proposal -= 1;
            }
            just_reverted_globals_for_changed_profile_length = false;

            //if( num_iterations_until_next_length_change_proposal == 0 ) {
              profile_length_changed_last_iteration =
                profile_length_changed_this_iteration;
              profile_length_changed_this_iteration = false;
            //}

            cycle_detected_last_iteration = cycle_detected_this_iteration;
            cycle_detected_this_iteration = false;
            // TOOD: REMOVE.  TESTING.
            // Use the trained-once globals for the initial global values to
            // revert to..
            if( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) == 1 ) {
              m_globalsToRevertTo.copyExceptPositions( *m_profile );
              // TODO: REMOVE! MAGIC #
              m_globalsToRevertTo.multiplyByExceptPositions( 10.0 );
            }

            // TODO: REMOVE.  TESTING.
            //if( m_iteration > 10 ) {
            //  globals_are_at_starting_values = true;
            //  m_trainingParameters.numIterationsBetweenLengthChanges = 1;
            //}

            // TODO: REMOVE.  TESTING.
            //if( m_iteration > 12 ) {
            //  m_trainingParameters.proposeProfileLengthChanges = false;
            //}

            // TODO: REMOVE.  TESTING.  MAGIC #s
            // This is very useful, it seems, for the Tigger2a data, but messes up the ProfuseTest results.
            // This ensures an earlier convergence:
            //const uint32_t scale_threshold_increments_every_X_iterations = 50;//10;
            //const double threshold_increment_scalar = 10.0;
            //const double threshold_increment_maximum = .05;
            //if(
            //  false &&
            //  (
            //    (
            //      ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) > 0 ) &&
            //      ( ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) % scale_threshold_increments_every_X_iterations ) == 0 )
            //    )
            //    // TODO: REMOVE.  TESTING.
            //    //||
            //    //( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) == 10 )
            //  )
            //) {
            //  if(
            //    false &&
            //    ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) == 10 ) 
            //  ) {
            //    // jump
            //    //  propose_inserting_threshold = .5;
            //    //  propose_inserting_prealign_threshold = .5;
            //    //  propose_inserting_postalign_threshold = .5;
            //    //  propose_deleting_threshold = .5;
            //
            //    // TODO: REMOVE!  TESTING trying to get a longer profile for the Tigger2a data.
            //    // favor insertions: lower the insertion threshold (will get incremented below)
            //    propose_inserting_threshold = 0;
            //    propose_inserting_prealign_threshold = 0;
            //    propose_inserting_postalign_threshold = 0;
            //    // Also raise the deletion rate, just to make sure.
            //    propose_deleting_threshold = 1.0;
            //  }
            //
            //  // extra boost
            //  propose_inserting_threshold_increment *=
            //    threshold_increment_scalar;
            //  if(
            //    propose_inserting_threshold_increment >
            //    threshold_increment_maximum
            //  ) {
            //    propose_inserting_threshold_increment =
            //      threshold_increment_maximum;
            //  }
            //  propose_deleting_threshold_increment *=
            //    threshold_increment_scalar;
            //  if(
            //    propose_deleting_threshold_increment >
            //    threshold_increment_maximum
            //  ) {
            //    propose_deleting_threshold_increment =
            //      threshold_increment_maximum;
            //  }
            //
            //  // Apply it once
            //  propose_inserting_threshold +=
            //    propose_inserting_threshold_increment;
            //  propose_inserting_prealign_threshold +=
            //    propose_inserting_threshold_increment;
            //  propose_inserting_postalign_threshold +=
            //    propose_inserting_threshold_increment;
            //  propose_deleting_threshold +=
            //    propose_deleting_threshold_increment;
            //            
            //  if( cout_profile_length_changes ) {
            //    cout << "Thresholds are being increased because the iteration is a multiple of " << scale_threshold_increments_every_X_iterations << ".  propose_deleting_threshold increased to " << propose_deleting_threshold << " and propose_inserting_threshold increased to " << propose_inserting_threshold << "." << endl;
            //  } // End if cout_profile_length_changes
            //  if(
            //    (
            //      ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
            //      ( m_trainingParameters.verbosity < VERBOSITY_High )
            //    ) ||
            //    ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
            //  ) {
            //    if(
            //      ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
            //      ( propose_inserting_threshold == propose_inserting_postalign_threshold )
            //    ) {
            //      if( propose_inserting_threshold == propose_deleting_threshold ) {
            //        cout << "threshold:" << propose_inserting_threshold;
            //      } else {
            //        cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
            //      }
            //    } else {
            //      cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
            //    }
            //    if(
            //      m_trainingParameters.verbosity >= VERBOSITY_Medium
            //    ) {
            //      cout << endl;
            //    } else {
            //      cout.flush();
            //    }
            //  } // End if verbose
            //
            //  // reset
            //  have_deleted_a_pos = false;
            //  have_inserted_a_pos = false;
            //  last_length_change_was_deletion = false;
            //  last_length_change_was_insertion = false;
            //} // End if this iteration is 10, 20, etc.

            if( use_sensitive_threshold_increments ) {
              // These are per-iteration.
              min_insertion_fraction_passing_threshold_causing_a_cycle = 1;
              min_deletion_fraction_passing_threshold_causing_a_cycle = 1;
            }
          } // End if m_trainingParameters.proposeProfileLengthChanges
          
          // TODO: REMOVE.  TESTING.
          //cout << "percentChange: " << m_scorePercentChange;
          
          if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
              ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
              m_trainingParameters.verbosity >= VERBOSITY_All ) {
            cout << "{Iteration#" << ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) + 1 );
            cout.flush();
          }
          if( m_trainingParameters.verbosity >=
              VERBOSITY_High ) {
            cout << endl;
          }
          just_changed_iterations = false;
        } // End if ( true )

        for( m_trainingPhase =
               ( m_parameters.trainGlobalsFirst ?
                 TRAINING_PHASE_Globals :
                 TRAINING_PHASE_Positions );
             ( m_trainingPhase < trainingPhaseCount );
             ( m_parameters.trainGlobalsFirst ?
               ( ( m_trainingPhase == TRAINING_PHASE_Globals ) ?
                 ( m_trainingPhase = TRAINING_PHASE_Positions ) :
                 ( ( m_trainingPhase == TRAINING_PHASE_Positions ) ?
                   ( m_trainingPhase = ( TRAINING_PHASE_Globals + 1 ) ) :
                   ++m_trainingPhase ) ) :
               ++m_trainingPhase )
        ) {

          // Don't train positions if the user has requested that we don't.
          if( ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
              !( m_trainingParameters.trainProfilePositions && ( m_trainingParameters.maxPositionCycles > 0 ) ) ) {
            continue;
          }

          // Don't train globals if the user has requested that we don't.  Note
          // that if useUnconditionalBaumWelch is true, we use the Globals
          // phase even if trainProfileGlobals is false (since we update the
          // position-specific params during the globals phase in this case).
          if( ( m_trainingPhase == TRAINING_PHASE_Globals ) &&
              !( ( m_trainingParameters.trainProfileGlobals && ( m_trainingParameters.maxPositionCycles_globals > 0 ) ) || m_trainingParameters.useUnconditionalBaumWelch ) ) {
            continue;
          }

          // TODO: REMOVE
          //cout << "iteration is " << ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) << "; training phase is " << static_cast<int>( m_trainingPhase ) << "." << endl;
         
          for( m_position_cycle = 0;
               (
                ( m_trainingParameters.useUnconditionalBaumWelch ? ( m_position_cycle == 0 ) : true ) && // This sez only do one position cycle when using unconditional BW.
                 // Only use the percentChange after the first cycle...
                 (
                   ( m_position_cycle == 0 ) ||
                   (
                    ( m_position_cycle < ( ( ( m_trainingPhase != TRAINING_PHASE_Positions ) ? m_trainingParameters.maxPositionCycles_globals : m_trainingParameters.maxPositionCycles ) ) ) &&
                    // Note that we won't get to this point if
                    // revert_globals_for_changed_profile_length or
                    // just_reverted_globals_for_changed_profile_length is
                    // true, since we set the position cycle to be its max in
                    // those cases.
                    ( !use_score_percent_change || (
                     ( m_scorePercentChange =
                       percentChange(
                         m_endingScore,
                         m_startingScore_position_cycle
                       ) )
                     >
                     m_trainingParameters.scorePercentChangeMinimum_position_cycle
                    ) ) &&
                    (
                     profile_distance_position_cycle >
                     m_trainingParameters.euclideanDistanceMinimum_position_cycle
                    )
                   )
                 )
               );
               m_position_cycle++
          ) {

            // Start over for this position training cycle..
            m_startingScore_position_cycle = m_endingScore;
          
            // Remember the profile as it was.
            m_startingProfile_position_cycle.copyFrom( *m_profile );

            if( m_trainingPhase == TRAINING_PHASE_Globals ) {
              m_global_entente.zero(); // Don't forget to reset the entente each time around.
#ifdef DISALLOW_FLANKING_TRANSITIONS
              // Note that zero() won't actually zero these when DISALLOW_FLANKING_TRANSITIONS is set, so we do it explicitly:
              m_global_entente[
                Transition::fromPreAlign
              ].zero();
              m_global_entente[
                Transition::fromPostAlign
              ].zero();
#endif // DISALLOW_FLANKING_TRANSITIONS
            }
          
            if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
                m_trainingParameters.verbosity >= VERBOSITY_All ) {
              if( ( m_trainingPhase != TRAINING_PHASE_Positions ) ) {
                if( m_trainingParameters.useUnconditionalBaumWelch ) {
                  cout << "[Unconditional ";
                  if( m_trainingParameters.baldiLearningRate > 0 ) {
                    if( m_trainingParameters.siegelMaxFindingThePeakAttempts_positions > 0 ) {
                      cout << "BS";
                    } else {
                      cout << "Baldi";
                    }
                  } else {
                    cout << "BW";
                  }
                  cout << " Update:" << ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) + 1 ) << ( ubw_cbw_hybrid ? ( ( ( m_iteration % 2 ) == 0 ) ? "-Positions" : "-Globals" ) : "" );
                } else {
                  cout << "[Globals:" << ( m_iteration + 1 ) << "#" << ( m_position_cycle + 1 );
                }
              } else { // phase == TRAINING_PHASE_Positions
                if( m_trainingParameters.useUnconditionalBaumWelch ) {
                  if( ubw_cbw_hybrid ) {
                    if( ( m_iteration % 2 ) == 0 ) {
                      cout << "[Positions:" << ( static_cast<int>( m_iteration / 2 ) + 1 ) << "#" << ( m_position_cycle + 1 );
                    } else {
                      cout << "[Globals:" << ( static_cast<int>( m_iteration / 2 ) + 1 ) << "#" << ( m_position_cycle + 1 );
                    }
                  } else {
                    cout << "[Positions&Globals:" << ( m_iteration + 1 ) << "#" << ( m_position_cycle + 1 );
                  }
                } else {
                  cout << "[Positions:" << ( m_iteration + 1 ) << "#" << ( m_position_cycle + 1 );
                }
              }
              if( m_trainingParameters.verbosity >=
                         VERBOSITY_Medium ) {
                cout << endl;
              } else {
                cout.flush();
              }
            }
          
            // Foreach position of the profile, do the training steps.  The
            // forward matrix rows are offset by one from the profile positions.
            m_row_i = last_row;
            anchor_rows_iter =
              m_anchor_rows.end();
            anchor_rows_iter--;
            if( cout_anchor_rows_iter_changes ) {
              anchor_rows_iter_i = m_anchor_rows.size() - 1;
              cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
            }
            if(
              ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
              ( m_anchor_rows.m_storeEveryNthRow > 1 )
            ) {
              anchor_columns_iter =
                m_anchor_columns.end();
              anchor_columns_iter--;
              if( cout_anchor_columns_iter_changes ) {
                anchor_columns_iter_i = m_anchor_columns.size() - 1;
                cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
              }
            } // End if m_anchor_columns.m_storeEveryNthColumn > 0
          
            // Actually if we are using unconditional bw and not training
            // globals, then we don't need to go through all the rows.
            if( ( m_trainingPhase == TRAINING_PHASE_Globals ) &&
                  !m_trainingParameters.trainProfileGlobals ) {
              m_row_i = 0;
            }

            if( m_trainingParameters.proposeProfileLengthChanges ) {
              if( m_trainingPhase == TRAINING_PHASE_Positions ) {
                just_grew_profile = false;
                just_shrunk_profile = false;
                just_caught_preAlign_cycle = false;
              }
            } // End if m_trainingParameters.proposeProfileLengthChanges


            ///////////////////////////////////////////
            ///// Begin foreach m_row_i ////////////////
            ///////////////////////////////////////////
            do { // iterating on m_row_i, which is unsigned.
              // Break after after 1 for positions, unless looking to change
              // the length of the profile..

              if(
                ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
                !( m_trainingParameters.proposeProfileLengthChanges || ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) == 0 ) ) &&
                ( m_row_i == 0 )
              ) {
                break;
              } // End if this is row 0 and we're currently training positions (etc.)


              if(
                ( m_row_i != last_row ) &&
                // The anchor row is stored at the beginning of the block of
                // rows, so if the next row is the first of a block, then we
                // are in the block of a new anchor row, as we progress
                // backwards.
                ( ( ( m_row_i + 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 )
              ) {
                // TODO: REMOVE
                assert( anchor_rows_iter != m_anchor_rows.begin() );
                anchor_rows_iter--;
                if( cout_anchor_rows_iter_changes ) {
                  anchor_rows_iter_i--;
                  cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                }
              }
              if(
                ( m_row_i != last_row ) &&
                ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                ( m_anchor_rows.m_storeEveryNthRow > 1 )
              ) {
                // TODO: REMOVE
                assert( anchor_columns_iter != m_anchor_columns.begin() );
                anchor_columns_iter--;
                if( cout_anchor_columns_iter_changes ) {
                  anchor_columns_iter_i--;
                  cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                }
              } // End if ( not the last row ) and m_anchor_columns.m_storeEveryNthColumn > 0 

              // Remember the score before we do any modifications to the profile.
              m_startingScore_position_step = m_endingScore;

          
              if( m_trainingParameters.verbosity >=
                  VERBOSITY_High ) {
                if( ( m_trainingPhase != TRAINING_PHASE_Positions ) ) {
                  if( m_trainingParameters.useUnconditionalBaumWelch ) {
                    // Globals as part of unconditional bw
                    if( m_row_i == last_row ) {
                      cout <<
                        "===================================================" << endl;
                      cout << "Calculating entente values for profile globals for training iteration " <<
                        ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) + 1 ) << "." << endl;
                    } else if( m_row_i == 0 ) {
                      cout << ".done.  Using the entente values to update the profile using unconditional Baum-Welch." << endl;
                    } else {
                      cout << ".";
                      //cout.flush();
                    } // End if m_row_i == last_row .. else if m_row_i == 0 .. else ..
                  } else {
                    if( m_row_i == last_row ) {
                      cout <<
                        "===================================================" << endl;
                      cout << "Training profile globals for training iteration " <<
                        ( m_iteration + 1 ) <<
                        ", globals training cycle " << ( m_position_cycle + 1 ) << "." <<
                        endl;
                    } else if( m_row_i == 0 ) {
                      cout << ".done." << endl;
                    } else {
                      cout << ".";
                      //cout.flush();
                    } // End if m_row_i == last_row .. else if m_row_i == 0 .. else ..
                  } // End if useUnconditionalBaumWelch .. else ..
                } else if( m_row_i > 0 ) { // if ( m_trainingPhase != TRAINING_PHASE_Positions ) .. else .. 
                  cout <<
                    "===================================================" << endl;
                  if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                      !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                    cout << "Skipping profile position " << ( m_row_i - 1 ) << endl;
                  } else if( m_trainingParameters.useUnconditionalBaumWelch ) {
                    cout << "Calculating entente values for profile position " << ( m_row_i - 1 ) <<
                      ", via forward matrix row " << m_row_i <<
                      " ( from " << last_row << " through 1 )" <<
                      " for training iteration " <<
                      ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) + 1 ) << "." << endl;
                  } else {
                    cout << "Training profile position " << ( m_row_i - 1 ) <<
                      ", via forward matrix row " << m_row_i <<
                      " ( from " << last_row << " through 1 )" <<
                      " for training iteration " <<
                      ( m_iteration + 1 ) <<
                      ", position cycle " << ( m_position_cycle + 1 ) << "." <<
                      endl;
                  }
                  cout << endl;
                } // End if ( m_trainingPhase != TRAINING_PHASE_Positions ) .. else .. 
              } // End if verbose

              if(
                 ( m_trainingPhase == TRAINING_PHASE_Globals ) &&
                 revert_globals_for_changed_profile_length
              ) {
                m_row_i = 0;
                if( m_trainingParameters.verbosity >= VERBOSITY_High ) {
                  cout << "_done." << endl;
                }
              } // End if the profile length just changed and this is the
                // Globals phase, skip all the intermediate rows so we can get straight to reverting..



                if( m_trainingPhase == TRAINING_PHASE_Positions ) {
                  if(
                    ( m_row_i > 0 ) &&
                    ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                    !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ]
                  ) {
                    // We're skipping the (parent profile) training for this
                    // position ..
                  } else {
                    // Reset the position entente for this new position...
                    m_position_entente.zero();
                    // Setting m_scalar = 0 signals to the
                    // ScalableParameterCollection::operator+=(
                    // ScalableParameterCollection ) code that it should just
                    // copy the entente update into the entente rather than
                    // scale the update to match the entente's scalar.  This is
                    // safer, since the right scale might be tiny and we don't
                    // know that until doing the fist update.
                    m_position_entente.m_scalar = 0;
                  
                    if(
                      m_trainingParameters.useAlignmentProfiles ||
                      m_trainingParameters.proposeProfileLengthChanges
                    ) {
                      alignment_profile_position_allseqs.zero();
                      // Setting m_scalar = 0 signals to the
                      // ScalableParameterCollection::operator+=(
                      // ScalableParameterCollection ) code that it should just
                      // copy the entente update into the entente rather than
                      // scale the update to match the entente's scalar.  This is
                      // safer, since the right scale might be tiny and we don't
                      // know that until doing the first update.
                      alignment_profile_position_allseqs.m_scalar = 0;
                      if(
                        m_trainingParameters.proposeProfileLengthChanges &&
                        count_seqs_exceeding_alignment_profile_thresholds &&
                        ( ubw_cbw_hybrid ? ( ( m_iteration % 2 ) == 0 ) : true )
                      ) {
                        num_seqs_exceeding_insertion_threshold = 0;
                        num_seqs_exceeding_deletion_threshold = 0;
                        num_seqs_exceeding_last_pos_insertion_threshold = 0;
                        num_seqs_exceeding_last_pos_deletion_threshold = 0;
                      } // End if count_seqs_exceeding_alignment_profile_thresholds
                      if( use_sensitive_threshold_increments ) {
                        min_insertion_fraction_passing_threshold = 1;
                        min_last_pos_insertion_fraction_passing_threshold = 1;
                        min_deletion_fraction_passing_threshold = 1;
                        min_last_pos_deletion_fraction_passing_threshold = 1;
                      } // End if use_sensitive_threshold_increments
                    } // end if m_trainingParameters.useAlignmentProfiles || m_trainingParameters.proposeProfileLengthChanges 
                  
                  } // End if we are skipping (parent profile) training for this position .. else ..
                } // End if trainingPhase == TRAINING_PHASE_Positions



                // Update the modifier ( The - 1 is because profile indices are
                // offset by 1 from the forward matrix indices. )
                //m_modifier.position( ( m_row_i - 1 ) );
              
                // Update the positions and backward_rows for this new position
                if( m_trainingParameters.debug > DEBUG_None ) {
                  cout <<
                    "[debug] Updating the positions and backward_rows for row " << m_row_i << "..";
                  cout.flush();
                } // End if debug

                // TODO: Add parameter to determine whether we use forward matrices or anchor rows/cols..
                //m_forward_rows_ptr = &m_forward_matrices[ m_row_i ];
                //if( m_row_i > 0 ) {
                //  m_prev_forward_rows_ptr = &m_forward_matrices[ m_row_i - 1 ];
                //}
                if(
                  m_trainingParameters.useUnconditionalBaumWelch &&
                  ( m_row_i == last_row )
                ) {
                  if( m_trainingPhase == TRAINING_PHASE_Positions ) {
                    m_unconditional_last_forward_rows =
                      ( *m_forward_rows_ptr );
                    m_unconditional_second_to_last_forward_rows =
                      ( *m_prev_forward_rows_ptr );
                    if( m_trainingParameters.useRabinerScaling ) {
                      for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                        // Copy the rabiner scalars.
                        m_unconditional_second_to_last_forward_rows[ m_seq_i ].m_rabinerInverseScalar =
                          ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                        m_unconditional_second_to_last_forward_rows[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                          ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                        m_unconditional_last_forward_rows[ m_seq_i ].m_rabinerInverseScalar =
                          ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                        m_unconditional_last_forward_rows[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                          ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                      } // End foreach sequence
                    } // End if useRabinerScaling
                  } else if( m_trainingPhase == TRAINING_PHASE_Globals ) {
                    m_forward_rows_1 =
                      m_unconditional_last_forward_rows;
                    m_forward_rows_2 =
                      m_unconditional_second_to_last_forward_rows;
                    if( m_trainingParameters.useRabinerScaling ) {
                      for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                        // Copy the rabiner scalars.
                        m_forward_rows_1[ m_seq_i ].m_rabinerInverseScalar =
                          m_unconditional_last_forward_rows[ m_seq_i ].m_rabinerInverseScalar;
                        m_forward_rows_1[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                          m_unconditional_last_forward_rows[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                        m_forward_rows_2[ m_seq_i ].m_rabinerInverseScalar =
                          m_unconditional_second_to_last_forward_rows[ m_seq_i ].m_rabinerInverseScalar;
                        m_forward_rows_2[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                          m_unconditional_second_to_last_forward_rows[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                      } // End foreach sequence
                    } // End if useRabinerScaling
                    m_forward_rows_ptr = &m_forward_rows_1;
                    m_prev_forward_rows_ptr = &m_forward_rows_2;
                    forward_rows_are_already_rotated_and_calculated = true;
                  } // End if this is the beginning of the first pass of an
                    // unconditional update, else if it is the beginning of the
                    // second pass .. (save/restore final two rows to/from
                    // backups).

                } // End if this is the beginning of an unconditional update,
                  // either save or restore the final two rows to/from backups.

                // TODO: REMOVE
                if( false && ( m_row_i == last_row ) ) {
                  cout << "last row is " << ( *m_forward_rows_ptr )[ 0 ] << endl;
                  cout << "last row m_rabinerInverseScalar is " << ( *m_forward_rows_ptr )[ 0 ].m_rabinerInverseScalar << ", and m_rabinerCumulativeInverseScalar is " << ( *m_forward_rows_ptr )[ 0 ].m_rabinerCumulativeInverseScalar << endl;
                } // End if testing


                if( !forward_rows_are_already_rotated_and_calculated ) {
                  // Rotate the forward_rows
                  m_temp_forward_rows_ptr = m_forward_rows_ptr;
                  m_forward_rows_ptr = m_prev_forward_rows_ptr;
                  m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;

                  // TODO: REMOVE
                  //cout << "After rotating the forward rows, forward row for row " << m_row_i << ", for seq 0, is " << ( *m_forward_rows_ptr )[ 0 ] << endl;
                  //cout << "OLD: " << m_forward_matrices[ m_row_i ][ 0 ] << endl;

                  // Now recalculate the previous forward row.
                  if( m_row_i != 0 ) {
                    if( ( ( m_row_i - 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                      // Restore the anchor row.
                      if( ( m_row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                        // Then the anchor row for the previous row is not the
                        // same as for this one.

                        // TODO: REMOVE
                        assert( anchor_rows_iter != m_anchor_rows.begin() );

                        tmp_matrices_iter = anchor_rows_iter;
                        tmp_matrices_iter--;
                        // TODO: Don't copy the data!  Copy the pointer..
                        ( *m_prev_forward_rows_ptr ) =
                          *tmp_matrices_iter;

                        // TODO: REMOVE
                        //assert( ( *m_prev_forward_rows_ptr )[ 0 ] == ( *tmp_matrices_iter )[ 0 ] );

                        if( m_trainingParameters.useRabinerScaling ) {
                          for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                            // Copy the rabiner scalars.
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                              ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar;
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                              ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                          } // End foreach sequence
                        } // End if useRabinerScaling
                        // TODO: REMOVE
                        //cout << "Restored anchor row for prev row ( " << ( m_row_i - 1 ) << " ): " << ( *m_prev_forward_rows_ptr )[ 0 ] << endl;
                        //cout << "Restored anchor row for prev row ( " << ( m_row_i - 1 ) << " ): m_rabinerInverseScalar is " << ( *m_prev_forward_rows_ptr )[ 0 ].m_rabinerInverseScalar << ", and m_rabinerCumulativeInverseScalar is " << ( *m_prev_forward_rows_ptr )[ 0 ].m_rabinerCumulativeInverseScalar << endl;
                        //cout << "OLD: " << m_forward_matrices[ m_row_i - 1 ][ 0 ].m_rabinerInverseScalar << " and " << m_forward_matrices[ m_row_i - 1 ][ 0 ].m_rabinerCumulativeInverseScalar << endl;
                      } else { // if the current row is also an anchor row  .. else ..
                        // Then the anchor row for the previous row is the
                        // same as for this one.

                        // TODO: Don't copy the data!  Copy the pointer..
                        ( *m_prev_forward_rows_ptr ) =
                          *anchor_rows_iter;
                        if( m_trainingParameters.useRabinerScaling ) {
                          for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                            // Copy the rabiner scalars.
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                              ( *anchor_rows_iter )[ m_seq_i ].m_rabinerInverseScalar;
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                              ( *anchor_rows_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                          } // End foreach sequence
                        } // End if useRabinerScaling
                        // TODO: REMOVE
                        //cout << "Restored anchor row for prev row ( " << ( m_row_i - 1 ) << " ): " << ( *m_prev_forward_rows_ptr )[ 0 ] << endl;
                        //cout << "Restored anchor row for prev row ( " << ( m_row_i - 1 ) << " ): m_rabinerInverseScalar is " << ( *m_prev_forward_rows_ptr )[ 0 ].m_rabinerInverseScalar << ", and m_rabinerCumulativeInverseScalar is " << ( *m_prev_forward_rows_ptr )[ 0 ].m_rabinerCumulativeInverseScalar << endl;
                        //cout << "OLD: " << m_forward_matrices[ m_row_i - 1 ][ 0 ].m_rabinerInverseScalar << " and " << m_forward_matrices[ m_row_i - 1 ][ 0 ].m_rabinerCumulativeInverseScalar << endl;
                      } // End if the current row is also an anchor row  .. else ..
                    } else { // if we can just restore the anchor row, then
                             // do it. .. else ..
                      for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                        if(
                          ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                          ( m_anchor_rows.m_storeEveryNthRow > 1 ) // redundant here
                        ) {
                          // TODO: REMOVE
                          assert( anchor_columns_iter != m_anchor_columns.begin() );
                          tmp_matrices_iter = anchor_columns_iter;
                          tmp_matrices_iter--;
                          ( *tmp_matrices_iter )[ m_seq_i ].restoreColumns(
                            m_anchor_columns.m_storeEveryNthColumn,
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ]
                          );
                          
                          if( m_trainingParameters.useRabinerScaling ) {
                            // Copy the rabiner scalars.
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                              ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar;
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                              ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                          } // End if useRabinerScaling
                        } // End if m_anchor_columns.m_storeEveryNthColumn > 0
                    
                        m_dynamic_programming.forward_reverseCalculateRow(
                          m_trainingParameters,
                          *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          (
                            ( m_row_i == 0 ) ?
                            m_scaled_match_distributions[ 0 ] : // ignored
                            m_scaled_match_distributions[ m_row_i - 1 ]
                          ),
                          (
                            ( m_row_i <= 1 ) ?
                            m_scaled_match_distributions[ 0 ] : // ignored
                            m_scaled_match_distributions[ m_row_i - 2 ]
                          ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          m_sequences[ m_seq_i ],
                          ( m_row_i - 1 ),
                          m_anchor_columns.m_storeEveryNthColumn,
                          ( *m_forward_rows_ptr )[ m_seq_i ],
                          ( *m_prev_forward_rows_ptr )[ m_seq_i ]
                        );
                    
                        if( false || ( m_parameters.debug >= DEBUG_All ) ) {
                          cout << endl;
                          // TODO: REMOVE
                          //cout << "Sequence " << m_seq_i << " length is " << m_sequences[ m_seq_i ].length() << endl;
                          cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] ";
                          cout << "NEW: " << ( *m_prev_forward_rows_ptr )[ m_seq_i ] << endl;
                          cout << endl;
                        } // End if DEBUG_All
                        // Now use the new one.
                      } // End foreach m_seq_i
                    } // End if ( m_row_i - 1 ) is an anchor row .. else .. 
                  } // End if ( row_i != 0 )
                } // End if !forward_rows_are_already_rotated_and_calculated

                if( !backward_rows_are_already_rotated_and_calculated ) {
                  // Rotate the backward_rows
                  m_temp_backward_rows_ptr = m_backward_rows_ptr;
                  m_backward_rows_ptr = m_next_backward_rows_ptr;
                  m_next_backward_rows_ptr = m_temp_backward_rows_ptr;

                  for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                    if( m_trainingParameters.useRabinerScaling ) {
                      ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                        ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                    } // End if useRabinerScaling
                    m_dynamic_programming.backward_calculateRow(
                       m_trainingParameters,
                       false, // don't use viterbi
                       *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                       (
                         ( m_row_i == 0 ) ?
                         m_scaled_match_distributions[ m_row_i ] : // ignored
                         m_scaled_match_distributions[ m_row_i - 1 ]
                       ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                       m_sequences[ m_seq_i ],
                       m_row_i,
                       ( *m_next_backward_rows_ptr )[ m_seq_i ],
                       ( *m_backward_rows_ptr )[ m_seq_i ]
                     );
                    if( m_trainingParameters.useRabinerScaling ) {
                      if( m_row_i == last_row ) {
                        ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                          ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                      } else {
                        ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                          (
                           ( *m_next_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar *
                           ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar
                          );
                      } // End if m_row_i == last_row .. else ..
                    } // End if useRabinerScaling
                  } // End foreach m_seq_i
                } // End if !backward_rows_are_already_rotated_and_calculated

                if( forward_rows_are_already_rotated_and_calculated ) {
                  forward_rows_are_already_rotated_and_calculated = false;
                }
                if( backward_rows_are_already_rotated_and_calculated ) {
                  backward_rows_are_already_rotated_and_calculated = false;
                }


                ////////////////////////////////////////////////////////
                // Update GEFS and alignmentProfilePosition for each seq
                // Also, Lengthadjust / Dynamic Model Surgery Stuff for each seq
                //
#ifdef DISALLOW_FLANKING_TRANSITIONS
      //  TODO: REMOVE
                assert( ( *m_profile )[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] == 0 );
                assert( ( *m_profile )[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] == 0 );
#endif // DISALLOW_FLANKING_TRANSITIONS
                for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                  if( false && ( m_parameters.debug >= DEBUG_All ) ) {
                    cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] " <<
                      "forward:" << ( *m_forward_rows_ptr )[ m_seq_i ] << endl;
                    cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] " <<
                      "backward:" << m_backward_rows_ptr->operator[]( m_seq_i ) << endl;
                  } // End if DEBUG_All

                  if( m_trainingPhase == TRAINING_PHASE_Globals ) {

                    // We could have ( m_trainingPhase ==
                    // TRAINING_PHASE_Globals ) true for some reason other than
                    // that we are actually going to change the global params:
                    // maybe we are using unconditional Baum-Welch and have
                    // turned global parameter training off. Thus the
                    // additional test:
                    if(
                      m_trainingParameters.trainProfileGlobals &&
                      !(
                        m_trainingParameters.proposeProfileLengthChanges &&
                        revert_globals_for_changed_profile_length
                      )
                    ) {
                      /// TODO: Also replace the global entente with a single AlignmentProfilePosition. Don't forget to account for the tricksy storage of pre-align and post-align insertion opens and extensions. SEE BELOW FOR CODE TESTING EQUIVALENCE
                      m_dynamic_programming.updateGlobalEntenteForSequence(
                        m_trainingParameters,
                        *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        (
                          ( m_row_i <= 1 ) ?
                          m_scaled_match_distributions[ m_row_i ] : // ignored
                          m_scaled_match_distributions[ m_row_i - 2 ]
                        ),
                        (
                          ( m_row_i == 0 ) ?
                          m_scaled_match_distributions[ m_row_i ] : // ignored
                          m_scaled_match_distributions[ m_row_i - 1 ]
                        ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        m_sequences[ m_seq_i ],
                        m_row_i,
                        ( ( m_row_i == 0 ) ?
                          ( *m_forward_rows_ptr )[ m_seq_i ] /*ignored*/ :
                        ( *m_prev_forward_rows_ptr )[ m_seq_i ] ),
                        ( *m_forward_rows_ptr )[ m_seq_i ],
                        m_backward_rows_ptr->operator[](  m_seq_i ),
                        &sequence_scores[ m_seq_i ],
                        m_global_entente
                      );
                      //cout << "[ " << m_row_i << ", " << m_seq_i << " ] Now m_global_entente is " << m_global_entente << endl;
                      //cout << "\tunscaled, m_global_entente is " << m_global_entente.createUnscaledCopy() << endl;

                      // TEST EQUIVALENCE of updateGEFS and calculateAPP
                      // TODO: REMOVE.  Testing.
//                      if( false ) {
//                        // Test that updateGEFS and
//                        // calculateAlignmentProfilePosition(..) do the same
//                        // thing -- but note that updateGEFS counts transitions
//                        // INTO a row, whereas
//                        // calculateAlignmentProfilePosition(..) counts
//                        // transitions OUT OF a row, so they won't be the same
//                        // thing until all added together.  We are using the
//                        // alignment_profile_position_allseqs as if it were a
//                        // global entente..
//
//                        m_dynamic_programming.calculateAlignmentProfilePosition(
//                          m_trainingParameters,
//                          *m_profile,
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//                          (
//                            ( m_row_i == 1 ) ?
//                            m_scaled_match_distributions[ m_row_i ] : // ignored
//                            m_scaled_match_distributions[ m_row_i - 2 ]
//                          ),
//                          m_scaled_match_distributions[ m_row_i - 1 ],
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
//                          m_sequences[ m_seq_i ],
//                          m_row_i,
//                          (
//                            ( m_row_i == 0 ) ?
//                            ( *m_forward_rows_ptr )[ m_seq_i ] : // ignored
//                            ( *m_prev_forward_rows_ptr )[ m_seq_i ]
//                          ),
//                          ( *m_forward_rows_ptr )[ m_seq_i ],
//                          m_backward_rows_ptr->operator[](  m_seq_i ),
//                          m_next_backward_rows_ptr->operator[](  m_seq_i ),
//                          // TODO: pass NULL unless !alwaysAccept ?
//                          &sequence_scores[ m_seq_i ],
//                          alignment_profile_position,
//                          &m_coefficients_vector[ m_seq_i ]
//                        );
//                        // TODO: REMOVE
//                        cout << "(1) alignment profile position is " << alignment_profile_position << endl;
//                        //if( m_row_i > 0 ) {
//                        //  cout << "sequence " << m_seq_i << " score, calculated using these coefficients, is " << m_coefficients_vector[ m_seq_i ].calculateScore( ( *m_profile )[ m_row_i - 1 ] ) << endl;
//                        //}
//
//#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
//                        if( m_row_i == 0 ) {
//                          // We've been tricksy and have stored the pre-align
//                          // opens separately from the extensions.  Now we want
//                          // to put them back to how they 'should' be.
//                          alignment_profile_position[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] +=
//                            alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
//                          alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = 0.0;
//                          alignment_profile_position[ Emission::PreAlignInsertion ] +=
//                            alignment_profile_position[ Emission::Insertion ];
//                          alignment_profile_position[ Emission::Insertion ].zero();
//                        } // End if m_row_i == 0
//                        if( m_row_i == last_row ) {
//                          // We've been tricksy and have stored the post-align
//                          // opens separately from the extensions.  Now we want
//                          // to put them back to how they 'should' be.
//                          alignment_profile_position[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] +=
//                            alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
//                          alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = 0.0;
//                          alignment_profile_position[ Emission::PostAlignInsertion ] +=
//                            alignment_profile_position[ Emission::Insertion ];
//                          alignment_profile_position[ Emission::Insertion ].zero();
//                        } // End if m_row_i == last_row
//#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
//                        if( ( m_seq_i == 0 ) && ( m_row_i == last_row ) ) {
//                          alignment_profile_position_allseqs.zero();
//                          alignment_profile_position_allseqs.m_scalar = 0;
//                        } // End if starting over, reset alignment_profile_position_allseqs...
//                        // 
//                        //cout << "About to call (alignment_profile_position_allseqs:" << alignment_profile_position_allseqs << ") += " << alignment_profile_position << endl;
//                        alignment_profile_position_allseqs +=
//                          alignment_profile_position;
//                        //cout << "[ " << m_row_i << ", " << m_seq_i << " ] Now alignment_profile_position_allseqs is " << alignment_profile_position_allseqs << endl;
//                        //cout << "\tunscaled, alignment_profile_position_allseqs is " << alignment_profile_position_allseqs.createUnscaledCopy() << endl;
//                      } // End if testing

                    } // End if trainProfileGlobals
                  } else { // End if TRAINING_PHASE_Globals .. else ..
                    // ( trainingPhase == TRAINING_PHASE_Positions )
                    if(
                      m_trainingParameters.useAlignmentProfiles ||
                      m_trainingParameters.proposeProfileLengthChanges
                    ) {
                      // TODO: REMOVE tmp_score
                      tmp_score =
                      m_dynamic_programming.calculateAlignmentProfilePosition(
                        m_trainingParameters,
                        *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        (
                          ( m_row_i <= 1 ) ?
                          m_scaled_match_distributions[ m_row_i ] : // ignored
                          m_scaled_match_distributions[ m_row_i - 2 ]
                        ),
                        (
                          ( m_row_i == 0 ) ?
                          m_scaled_match_distributions[ m_row_i ] : // ignored
                          m_scaled_match_distributions[ m_row_i - 1 ]
                        ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        m_sequences[ m_seq_i ],
                        m_row_i,
                        (
                          ( m_row_i == 0 ) ?
                          ( *m_forward_rows_ptr )[ m_seq_i ] : // ignored
                          ( *m_prev_forward_rows_ptr )[ m_seq_i ]
                        ),
                        ( *m_forward_rows_ptr )[ m_seq_i ],
                        m_backward_rows_ptr->operator[](  m_seq_i ),
                        m_next_backward_rows_ptr->operator[](  m_seq_i ),
                        // TODO: I think we should be able to put NULL here (unless useMaximumValue (in DynamicProgramming params) is true, maybe?
                        &sequence_scores[ m_seq_i ],
                        alignment_profile_position,
                        &m_coefficients_vector[ m_seq_i ]
                      );
                      // TODO: REMOVE
                      //cout << "(2) alignment profile position is " << alignment_profile_position << endl;
                      //cout << "coefficients, calculated using calculateAlignmentProfilePosition(..), are " << m_coefficients_vector[ m_seq_i ] << endl;
                      // TODO: REMOVE
                      // This check to see if the scores are right only makes sense if we recalculate the scores, which we only do when !alwaysAccept..
#ifndef NDEBUG // if assert(..) is defined:
                      if( true && !( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) && ( m_row_i > 0 ) ) { // debugging...
                        sequence_score =
                          m_coefficients_vector[ m_seq_i ].calculateScore( ( *m_profile )[ m_row_i - 1 ] );
                        //  cout << "sequence " << m_seq_i << " score, calculated using these coefficients, is " << sequence_score << endl;
                        // TODO: REMOVE?
                        if( sequence_score > sequence_scores[ m_seq_i ] ) {
                          if( toDouble( sequence_score - sequence_scores[ m_seq_i ] ) >= 1E-5 ) {
                            cout << "UH-OH: sequence " << m_seq_i << " score, calculated using the coefficients for row " << m_row_i << ", is " << sequence_score << ", but it should be " << sequence_scores[ m_seq_i ] << endl;
                            // TODO: REMOVE
                            cout << "score, calculated using calculateAlignmentProfilePosition(..), is " << tmp_score << endl;
                            cout << "( *m_prev_forward_rows_ptr )[ m_seq_i ] is " << ( *m_prev_forward_rows_ptr )[ m_seq_i ] << endl;
                            cout << "m_backward_rows_ptr->operator[](  m_seq_i ) is " << m_backward_rows_ptr->operator[](  m_seq_i ) << endl;
                            cout << "( *m_forward_rows_ptr )[ m_seq_i ] is " << ( *m_forward_rows_ptr )[ m_seq_i ] << endl;

                            cout << "alignment profile position is " << alignment_profile_position << endl;
                            alignment_profile_position.unscale();
                            cout << "unscaled, alignment profile position is " << alignment_profile_position << endl;
                            cout << "coefficients, calculated using calculateAlignmentProfilePosition(..), are " << m_coefficients_vector[ m_seq_i ] << endl;
                            cout << "( *m_profile )[ m_row_i - 1 ] is " <<  ( *m_profile )[ m_row_i - 1 ] << endl;
                            cout << "profile globals are ";
                            m_profile->writeExceptPositions( cout );
                            cout << endl;
                            m_dynamic_programming.calculatePositionSpecificSequenceScoreCoefficients(
                              m_trainingParameters,
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              (
                                ( m_row_i <= 1 ) ?
                                m_scaled_match_distributions[ m_row_i ] : // ignored
                                m_scaled_match_distributions[ m_row_i - 2 ]
                              ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences[ m_seq_i ],
                              m_row_i,
                              ( *m_prev_forward_rows_ptr )[ m_seq_i ],
                              m_backward_rows_ptr->operator[](  m_seq_i ),
                              m_coefficients_vector[ m_seq_i ]
                            );
                            cout << "coefficients, calculated using calculatePositionSpecificSequenceScoreCoefficients(..), are " << m_coefficients_vector[ m_seq_i ] << endl;
                            cout << "score, recalculated using forward_score( Parameters, bool, Row, Row ), is " <<
                              m_dynamic_programming.forward_score(
                                m_trainingParameters,
                                *m_profile,
                                m_row_i,
                                ( *m_forward_rows_ptr )[ m_seq_i ],
                                m_backward_rows_ptr->operator[](  m_seq_i )
                              ) << endl;
                            cout << "score, recaculated using the most complex forward_score(..), is ";
                            m_dynamic_programming.forward_score(
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences,
                              m_sequence_count,
                              &m_anchor_columns,
                              &m_anchor_rows,
                              *m_prev_forward_rows_ptr,
                              *m_forward_rows_ptr,
                              &sequence_scores,
                              &largest_sequence_score,
                              m_endingScore
                            );
                            cout << sequence_scores[ m_seq_i ] << endl;
                            assert( false );
                            assert( toDouble( sequence_score - sequence_scores[ m_seq_i ] ) < 1E-5 );
                          }
                        } else if( sequence_score < sequence_scores[ m_seq_i ] ) {
                          if( toDouble( sequence_scores[ m_seq_i ] - sequence_score ) >= 1E-5 ) {
                            cout << "UH-OH: sequence " << m_seq_i << " score, calculated using the coefficients for row " << m_row_i << ", is " << sequence_score << ", but it should be " << sequence_scores[ m_seq_i ] << endl;
                            // TODO: REMOVE
                            cout << "score, calculated using calculateAlignmentProfilePosition(..), is " << tmp_score << endl;
                            cout << "( *m_prev_forward_rows_ptr )[ m_seq_i ] is " << ( *m_prev_forward_rows_ptr )[ m_seq_i ] << endl;
                            cout << "m_backward_rows_ptr->operator[](  m_seq_i ) is " << m_backward_rows_ptr->operator[](  m_seq_i ) << endl;
                            cout << "( *m_forward_rows_ptr )[ m_seq_i ] is " << ( *m_forward_rows_ptr )[ m_seq_i ] << endl;
                            cout << "m_next_backward_rows_ptr->operator[](  m_seq_i ) is " << m_next_backward_rows_ptr->operator[](  m_seq_i ) << endl;

                            cout << "alignment profile position is " << alignment_profile_position << endl;
                            alignment_profile_position.unscale();
                            cout << "unscaled, alignment profile position is " << alignment_profile_position << endl;
                            cout << "coefficients, calculated using calculateAlignmentProfilePosition(..), are " << m_coefficients_vector[ m_seq_i ] << endl;
                            cout << "( *m_profile )[ m_row_i - 1 ] is " <<  ( *m_profile )[ m_row_i - 1 ] << endl;
                            cout << "profile globals are ";
                            m_profile->writeExceptPositions( cout );
                            cout << endl;
                            m_dynamic_programming.calculatePositionSpecificSequenceScoreCoefficients(
                              m_trainingParameters,
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              (
                                ( m_row_i <= 1 ) ?
                                m_scaled_match_distributions[ m_row_i ] : // ignored
                                m_scaled_match_distributions[ m_row_i - 2 ]
                              ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences[ m_seq_i ],
                              m_row_i,
                              ( *m_prev_forward_rows_ptr )[ m_seq_i ],
                              m_backward_rows_ptr->operator[](  m_seq_i ),
                              m_coefficients_vector[ m_seq_i ]
                            );
                            cout << "coefficients, calculated using calculatePositionSpecificSequenceScoreCoefficients(..), are " << m_coefficients_vector[ m_seq_i ] << endl;
                            cout << "score, recalculated using forward_score( Parameters, bool, Row, Row ), is " <<
                              m_dynamic_programming.forward_score(
                                m_trainingParameters,
                                *m_profile,
                                m_row_i,
                                ( *m_forward_rows_ptr )[ m_seq_i ],
                                m_backward_rows_ptr->operator[](  m_seq_i )
                              ) << endl;
                            cout << "score, recaculated using the most complex forward_score(..), is ";
                            m_dynamic_programming.forward_score(
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences,
                              m_sequence_count,
                              &m_anchor_columns,
                              &m_anchor_rows,
                              *m_prev_forward_rows_ptr,
                              *m_forward_rows_ptr,
                              &sequence_scores,
                              &largest_sequence_score,
                              m_endingScore
                            );
                            cout << sequence_scores[ m_seq_i ] << endl;
                            assert( false );
                            assert( toDouble( sequence_scores[ m_seq_i ] - sequence_score ) < 1E-5 );
                          }
                        }
                      } // End DEBUGGING: if !alwaysAccept and m_row_i > 0, do some debug checking..
#endif // !NDEBUG

                      // if we aren't just using the
                      // alignment_profile_position_allseqs to determine adding
                      // or removing a position, then we also need to calculate
                      // the insertion_fraction and deletion_fraction for
                      // individual sequences to count how many exceed the
                      // threshold.
                      if(
                        m_trainingParameters.proposeProfileLengthChanges &&
                        count_seqs_exceeding_alignment_profile_thresholds &&
                        ( ubw_cbw_hybrid ? ( ( m_iteration % 2 ) == 0 ) : true )
                      ) {
                        // In the last row we store deletions in the Match state,
                        // so we have to wait until the row is ( last_row - 1 )
                        // to find out when we deleted there.
                        if(
                          ( m_row_i < last_row ) &&
                          ( m_row_i > 0 )
                        ) {
                          delete_last_pos = false;
                          if( m_row_i == ( last_row - 1 ) ) {
                            // First check up on deletions of the final position,
                            // since we couldn't do it when m_row_i was last_row
                            // (because we don't store deletions in the deletion
                            // state there).
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          last_pos_deletion_fraction =
                            toDouble( (  alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] + alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] + alignment_profile_position[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] + alignment_profile_position[ Transition::fromDeletionOut ].total() ) / ( ( alignment_profile_position[ Transition::fromMatch ].total() + alignment_profile_position[ Transition::fromDeletion ].total() + alignment_profile_position[ Transition::fromDeletionOut ].total() ) ) );
#else // !USE_DEL_IN_DEL_OUT
                          last_pos_deletion_fraction =
                            toDouble( (  alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] + alignment_profile_position[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] ) / ( ( alignment_profile_position[ Transition::fromMatch ].total() + alignment_profile_position[ Transition::fromDeletion ].total() ) ) );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                            if( cout_indel_fractions_seqs ) {
                              cout << "[ " << m_row_i << ", " << m_seq_i << " ] (last pos) Deletion fraction (sequence " << m_seq_i << ") " << last_pos_deletion_fraction << endl;
                            }
#ifndef DISALLOW_FLANKING_TRANSITIONS
                            // TODO: REMOVE
                            if( last_pos_deletion_fraction == 0 ) {
                              cout << "[ " << m_row_i << ", " << m_seq_i << " ] Somehow the last_pos_deletion_fraction is 0?!  The alignment_profile_position is " << alignment_profile_position << ", " << " the profile position is " << ( *m_profile )[ m_row_i + 1 - 1 ] << ", the coefficients vector is " << m_coefficients_vector[ m_seq_i ] << endl;
                            }
#endif // !DISALLOW_FLANKING_TRANSITIONS
                            // TODO: REMOVE
                            //assert( last_pos_deletion_fraction > 0 );

                            // TODO: REMOVE
                            assert( last_pos_deletion_fraction >= 0.0 );
                            assert( last_pos_deletion_fraction <= 1.0 );

                            if(
                              ( last_pos_deletion_fraction > propose_deleting_threshold )
                            ) {
                              // TODO: REMOVE
                              //cout << "DELETING LAST POS" << endl;
                              pos_i =
                                ( last_row - 1 );
                              num_seqs_exceeding_last_pos_deletion_threshold += 1;
                              if(
                                use_sensitive_threshold_increments &&
                                (
                                  last_pos_deletion_fraction <
                                  min_last_pos_deletion_fraction_passing_threshold
                                )
                              ) {
                                min_last_pos_deletion_fraction_passing_threshold =
                                  last_pos_deletion_fraction;
                              }
                              if( cout_profile_length_changes_seqs ) {
                                cout << "[ " << m_row_i << ", " << m_seq_i << " ]\t\tProfile position " << pos_i << " exceeds deletion threshold." << endl;
                                if( !cout_indel_fractions_seqs ) { // We've already printed it.
                                  cout << "\t Last pos Deletion fraction (sequence " << m_seq_i << ") " << last_pos_deletion_fraction << endl;
                                  //cout << "[ " << m_row_i << ", " << m_seq_i << " ] alignment profile position is " << alignment_profile_position << endl;
                                } // End if !cout_indel_fractions_seqs
                              } // End if cout_profile_length_changes_seqs
                            } // End if last_pos_deletion_fraction > propose_deleting_threshold
                          } // End if ( m_row_i == ( last_row - 1 ) )
                          // TODO: In the final row, will this work?
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          deletion_fraction =
                            toDouble( ( alignment_profile_position[ Transition::fromDeletion ].total() + alignment_profile_position[ Transition::fromDeletionIn ].total() + alignment_profile_position[ Transition::fromDeletionOut ].total() ) / ( ( alignment_profile_position[ Transition::fromMatch ].total() + alignment_profile_position[ Transition::fromDeletion ].total() + alignment_profile_position[ Transition::fromDeletionIn ].total() + alignment_profile_position[ Transition::fromDeletionOut ].total() ) ) );
#else // !USE_DEL_IN_DEL_OUT
                          deletion_fraction =
                            toDouble( ( alignment_profile_position[ Transition::fromDeletion ].total() ) / ( ( alignment_profile_position[ Transition::fromMatch ].total() + alignment_profile_position[ Transition::fromDeletion ].total() ) ) );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                          if( cout_indel_fractions_seqs ) {
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ] Deletion fraction (sequence " << m_seq_i << ") " << deletion_fraction << endl;
                          }
#ifndef DISALLOW_FLANKING_TRANSITIONS
                          // TODO: REMOVE
                          //if( deletion_fraction == 0 ) {
                          //  cout << "[ " << m_row_i << ", " << m_seq_i << " ] Somehow the deletion_fraction is 0?!  The alignment_profile_position is " << alignment_profile_position << endl;
                          //}
#endif // !DISALLOW_FLANKING_TRANSITIONS
                          // TODO: REMOVE
                          //assert( deletion_fraction > 0 );

                          // TODO: REMOVE
                          assert( deletion_fraction >= 0.0 );
                          assert( deletion_fraction <= 1.0 );

                          if(
                            deletion_fraction > propose_deleting_threshold
                          ) {
                            pos_i =
                              ( m_row_i - 1 );
                            num_seqs_exceeding_deletion_threshold += 1;
                            if(
                              use_sensitive_threshold_increments &&
                              (
                                deletion_fraction <
                                min_deletion_fraction_passing_threshold
                              )
                            ) {
                              min_deletion_fraction_passing_threshold =
                                deletion_fraction;
                            }

                            if( cout_profile_length_changes_seqs ) {
                              cout << "[ " << m_row_i << ", " << m_seq_i << " ]\t\tProfile position " << pos_i << " exceeds deletion threshold." << endl;
                              if( !cout_indel_fractions_seqs ) { // We've already printed it.
                                cout << "\t Deletion fraction (sequence " << m_seq_i << ") " << deletion_fraction << endl;
                                //cout << "[ " << m_row_i << ", " << m_seq_i << " ] alignment profile position is " << alignment_profile_position << endl;
                              } // End if !cout_indel_fractions_seqs
                            } // End if cout_profile_length_changes_seqs
                          } // End if deletion_fraction > propose_deleting_threshold
                        } // End if it's not the first nor the last row

                        if( m_row_i == 0 ) {
                          // Note that we're only using the insertion _open_, which we've tricksily stored separately from the insertion extension.
                          // Sanity check
                          // TODO: Why is there apparent numeric error in the scalar?
                          //assert( ( alignment_profile_position[ Transition::fromBegin ].total() / alignment_profile_position.m_scalar ) == 1 );
                          assert( abs( toDouble( alignment_profile_position[ Transition::fromBegin ].total() / alignment_profile_position.m_scalar ) - 1 ) < 1E-5 ); // MAGIC #!

                          // Pre-align..
                          insertion_fraction =
                            toDouble( alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] / alignment_profile_position.m_scalar );
                          //toDouble( alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] / alignment_profile_position[ Transition::fromBegin ].total() );
                          if( cout_indel_fractions_seqs ) {
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ]  Pre-align insertion fraction (sequence " << m_seq_i << "): " << insertion_fraction << endl;
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ] alignment profile position is " << alignment_profile_position << endl;
                          }
                          // TODO: REMOVE
                          assert( insertion_fraction >= -1E-5 );
                          assert( ( 1.0 - insertion_fraction ) > -1E-5 );
                        } else if( m_row_i == last_row ) {
                          // Note that we're only using the insertion _open_, which we've tricksily stored separately from the insertion extension.
                          assert( ( alignment_profile_position[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] / alignment_profile_position.m_scalar ) == 1 );
                          last_pos_insertion_fraction =
                            //toDouble( ( alignment_profile_position[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] ) / ( alignment_profile_position[ Transition::fromPostAlign ].total() ) );
                            toDouble( ( alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / alignment_profile_position.m_scalar );
                            //toDouble( ( alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] ) );
                          // TODO: Put back.
                          //if( cout_indel_fractions_seqs ) {
                          if( ( 1.0 - last_pos_insertion_fraction ) <= -1E-5 ) {
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ]  Post-align insertion fraction (sequence " << m_seq_i << "): " << last_pos_insertion_fraction << endl;
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ] alignment profile position is " << alignment_profile_position << endl;
                            cout << "( 1.0 - last_pos_insertion_fraction ) is " << ( 1.0 - last_pos_insertion_fraction ) << endl;
                            cout << "-1E-5 is " << -1E-5 << endl;
                          }
                          // TODO: REMOVE
                          assert( last_pos_insertion_fraction >= -1E-5 );
                          assert( ( 1.0 - last_pos_insertion_fraction ) > -1E-5 );
                        } else {
                          insertion_fraction =
                            // Denom should include paths that skip this position.
                            toDouble( ( alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position[ Transition::fromMatch ].total() + alignment_profile_position[ Transition::fromDeletion ].total() ) );
                          // TODO: PUT BACK!
                            //toDouble( ( alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position[ Transition::fromMatch ].total() ) );
                          if( cout_indel_fractions_seqs ) {
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ]  Insertion fraction (sequence " << m_seq_i << "): " << insertion_fraction << endl;
                            //cout << "[ " << m_row_i << ", " << m_seq_i << " ] alignment profile position is " << alignment_profile_position << endl;
                          }
                          // TODO: REMOVE
                          assert( insertion_fraction >= -1E-5 );
                          assert( ( 1.0 - insertion_fraction ) > -1E-5 );
                        } // End switch (first row, last row, else)

                        //cout << "[ " << m_row_i << ", " << m_seq_i << " ] Insertion fraction (sequence " << m_seq_i << ") " << insertion_fraction << endl;
                        if(
                           ( ( m_row_i == last_row ) ? ( last_pos_insertion_fraction > propose_inserting_postalign_threshold ) : ( insertion_fraction > ( ( m_row_i == 0 ) ? propose_inserting_prealign_threshold : propose_inserting_threshold ) ) )
                        ) {
                          if( m_row_i == last_row ) {
                            num_seqs_exceeding_last_pos_insertion_threshold += 1;
                          } else {
                            num_seqs_exceeding_insertion_threshold += 1;
                          }
                          if(
                             use_sensitive_threshold_increments &&
                             ( m_row_i == last_row ) &&
                             (
                              last_pos_insertion_fraction <
                              min_last_pos_insertion_fraction_passing_threshold
                              ) 
                          ) {
                            min_last_pos_insertion_fraction_passing_threshold =
                              last_pos_insertion_fraction;
                          }

                          if( use_sensitive_threshold_increments &&
                              ( m_row_i < last_row ) &&
                              (
                               insertion_fraction <
                                min_insertion_fraction_passing_threshold
                              )
                          ) {
                            min_insertion_fraction_passing_threshold =
                              insertion_fraction;
                          }

                          if( cout_profile_length_changes_seqs ) {
                            cout << "[ " << m_row_i << ", " << m_seq_i << " ]\t\tProfile position " << m_row_i << " exceeds insertion threshold." << endl;
                            if( !cout_indel_fractions_seqs ) { // we've already printed it.
                              cout << "\t Insertion fraction (sequence " << m_seq_i << ") " << ( ( m_row_i == last_row ) ? last_pos_insertion_fraction : insertion_fraction ) << endl;
                              //cout << "[ " << m_row_i << ", " << m_seq_i << " ] alignment profile position is " << alignment_profile_position << endl;
                            }
                          } // End if cout_profile_length_changes_seqs
                        } // End if insertion_fraction > propose_inserting_threshold
                      } // End if count_seqs_exceeding_alignment_profile_thresholds

                      // 
                      //cout << "About to call (alignment_profile_position_allseqs:" << alignment_profile_position_allseqs << ") += " << alignment_profile_position << endl;
                      alignment_profile_position_allseqs +=
                        alignment_profile_position;
                      //cout << "Now alignment_profile_position_allseqs is " << alignment_profile_position_allseqs << endl;
                      //cout << "\tunscaled, alignment_profile_position_allseqs is " << alignment_profile_position_allseqs.createUnscaledCopy() << endl;
                      if( m_parameters.debug >= DEBUG_All ) {
                        cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] " <<
                          "Alignment profile position (just this sequence):" << alignment_profile_position << endl;
                        cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] " <<
                          "Alignment profile position (all sequences so far):" << alignment_profile_position_allseqs << endl;
                      } // End if DEBUG_ALL
                    } // end if m_trainingParameters.useAlignmentProfiles || m_trainingParameters.proposeProfileLengthChanges

                    if( m_row_i > 0 ) {
                      if( m_trainingParameters.useAlignmentProfiles ) {
                        if( m_row_i == 0 ) { // NOTE! NEVER TRUE!!!!
                          /// TODO: Also replace the global entente with a single AlignmentProfilePosition. Don't forget to account for the tricksy storage of pre-align and post-align insertion opens and extensions.

                          // We've been tricksy and have stored the pre-align
                          // opens separately from the extensions.  Now we want
                          // to put them back to how they 'should' be.

                          alignment_profile_position[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] +=
                            alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
                          alignment_profile_position[ Transition::fromMatch ] = 0.0;

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
                          alignment_profile_position[ Emission::PreAlignInsertion ] +=
                            alignment_profile_position[ Emission::Insertion ];
                          alignment_profile_position[ Emission::Insertion ].zero();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
                        } // End if m_row_i == 0
                        if( m_row_i == last_row ) {
                          // We've been tricksy and have stored the post-align
                          // opens separately from the extensions.  Now we want
                          // to put them back to how they 'should' be.
                          alignment_profile_position[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] +=
                            alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
                          alignment_profile_position[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = 0.0;
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
                          alignment_profile_position[ Emission::PostAlignInsertion ] +=
                            alignment_profile_position[ Emission::Insertion ];
                          alignment_profile_position[ Emission::Insertion ].zero();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
                        } // End if m_row_i == last_row
                        // This will just add the part of it that is relevant to
                        // the current profile position's parameters.
                        //cout << "About to call (m_position_entente:" << m_position_entente << ") += " << alignment_profile_position << endl;
                        m_position_entente += alignment_profile_position;
                        //cout << "Now m_position_entente is " << m_position_entente << endl;
                        //cout << "\tunscaled, m_position_entente is " << m_position_entente.createUnscaledCopy() << endl;
                        //m_dynamic_programming.updatePositionEntente(
                        //  m_trainingParameters,
                        //  alignment_profile_position,
                        //  m_position_entente
                        //);
                      
                        // TODO: REMOVE
                        //cout << "[Row " << m_row_i << ", Sequence " << m_seq_i << "] ";
                        //cout << "m_position_entente (using alignment profile position) is now " << m_position_entente << endl;
                      } else { // if m_trainingParameters.useAlignmentProfiles .. else ..
                      
                        // Note that we can't skip this, even if we're supposed
                        // to skip training the parent profile for this
                        // position, since we need the coefficients and the
                        // score..
                        m_dynamic_programming.calculatePositionSpecificSequenceScoreCoefficients(
                          m_trainingParameters,
                          *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          (
                            ( m_row_i <= 1 ) ?
                            m_scaled_match_distributions[ m_row_i ] : // ignored
                            m_scaled_match_distributions[ m_row_i - 2 ]
                          ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          m_sequences[ m_seq_i ],
                          m_row_i,
                          ( *m_prev_forward_rows_ptr )[ m_seq_i ],
                          m_backward_rows_ptr->operator[](  m_seq_i ),
                          m_coefficients_vector[ m_seq_i ]
                        );
                        // TODO: REMOVE
                        //cout << "[Row " << m_row_i << ", Sequence " << m_seq_i << "] ";
                        //cout << "Coefficients (just this sequence and position):" << m_coefficients_vector[ m_seq_i ] << endl;
                        //cout << "sequence " << m_seq_i << " score, calculated using these coefficients, is " << m_coefficients_vector[ m_seq_i ].calculateScore( ( *m_profile )[ m_row_i - 1 ] ) << endl;

                        // TODO: This will be for bw and conditional bw but not
                        // for the gauss-newton maximization...  NOTE: we
                        // really don't need to store the vector of the
                        // coefficients if we are using gauss-newton, since
                        // here is the only place we access it.
                        m_dynamic_programming.updatePositionEntenteForSequence(
                          m_trainingParameters,
                          ( *m_profile )[ m_row_i - 1 ],
                          m_coefficients_vector[ m_seq_i ],
                          NULL,
                          m_position_entente
                        );
                        // TODO: REMOVE
                        //cout << "[Row " << m_row_i << ", Sequence " << m_seq_i << "] ";
                        //cout << "m_position_entente is now:" << m_position_entente << endl;

                      } // End if m_trainingParameters.useAlignmentProfiles .. else ..
                    } // End if m_row_i > 0
                    
                    if( m_parameters.debug >= DEBUG_All ) {
                      cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] " <<
                        "Coeffs (just this sequence):" << m_coefficients_vector[ m_seq_i ] << endl;
                      cout << "[debug] [Row " << m_row_i << ", Sequence " << m_seq_i << "] " <<
                        "Entente Pos (so far):" << m_position_entente << endl;
                    } // End if DEBUG_ALL
                    // TODO: REMOVE!
                    //exit( 0 );
                  } // End if trainingGlobals .. else if trainingPositions ..
                } // End foreach m_seq_i


                  
                if( m_trainingPhase == TRAINING_PHASE_Globals ) {
                  // TODO: REMOVE
                  //cout << "[ " << m_row_i << " ] Now m_global_entente is " << m_global_entente << endl;
                  //cout << "\tunscaled, m_global_entente is " << m_global_entente.createUnscaledCopy() << endl;
                  // TODO: REMOVE?
                  assert( !isnan( m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] ) );
#ifndef NDEBUG
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                  // Check that the global entente has sane values for del-in and del-out (after all rows and sequences have been accounted for).
                  if( m_row_i == 0 ) {
                    // TODO: DEHACKIFY MAGIC #.  Note allowance for some drift due to anchor rows & cols, so we have a fairly loose determination of equality as just being on the same scale.
                    if( m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] !=
                        m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] ) {
                      bool MW_gt_WE = ( m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ]  > m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] );
                      MatrixValueType diff = ( MW_gt_WE ? ( m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] - m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] ) : ( m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] - m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] ) );
                      // TODO: REMOVE
                      cout << "[ " << m_row_i << " ] Now m_global_entente is " << m_global_entente << endl;
                      cout << "\tunscaled, m_global_entente is " << m_global_entente.createUnscaledCopy() << endl;
                      // TODO: REMOVE
                      cout << "abs( m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ]  - m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] ) / max(..): " << ( diff / ( MW_gt_WE ? m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] : m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] ) ) << endl;
                      assert(  toDouble( ( diff / ( MW_gt_WE ? m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] : m_global_entente[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ] ) ) < 1E-5 ) ); // TODO: DEHACIFY MAGIC # 1E-5 !!!
                    }
                    if( m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] !=
                        m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] ) {
                      bool BZ_gt_ZM = ( m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ]  > m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] );
                      MatrixValueType diff = ( BZ_gt_ZM ? ( m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] - m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] ) : ( m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] - m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] ) );
                      // TODO: REMOVE
                      cout << "[ " << m_row_i << " ] Now m_global_entente is " << m_global_entente << endl;
                      cout << "\tunscaled, m_global_entente is " << m_global_entente.createUnscaledCopy() << endl;
                      // TODO: REMOVE
                      cout << "abs( m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ]  - m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] ) / max(..): " << ( diff / ( BZ_gt_ZM ? m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] : m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] ) ) << endl;
                      assert(  toDouble( ( diff / ( BZ_gt_ZM ? m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] : m_global_entente[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ] ) ) < 1E-5 ) ); // TODO: DEHACIFY MAGIC # 1E-5 !!!
                    }
                  } // End if m_row_i == 0
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#endif // !NDEBUG
                } // End if TRAINING_PHASE_Globals
                
                //////////////////////////////////////////////////////
                // Lengthadjust / Dynamic Model Surgery Stuff for allseqs
                //
                  
                if(
                  m_trainingParameters.proposeProfileLengthChanges &&
                  ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
                  (
                    ( m_row_i == 0 ) ||
                    !(
                      ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                      !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ]
                    )
                  ) &&
                  ( ubw_cbw_hybrid ? ( ( m_iteration % 2 ) == 0 ) : true ) &&
                  ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) >= length_change_min_iteration ) &&
                  ( num_iterations_until_next_length_change_proposal == 0 )
                ) {

                  //cout << "[ " << m_row_i << " ] alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs << endl;
                  //cout << "[ " << m_row_i << " ] unscaled alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs.createUnscaledCopy() << endl;
                  //if( m_row_i > 0 ) {
                  //  cout << "[ " << m_row_i << " ] m_position_entente is " << m_position_entente << endl;
                  //  cout << "[ " << m_row_i << " ] unscaled m_position_entente is " << m_position_entente.createUnscaledCopy() << endl;
                  //}
                
                  insertion_fraction_exceeds_threshold = false;
                  if( m_row_i != ( last_row - 1 ) ) {
                    // We need to remember what we learned from last_row when
                    // m_row_i is ( last_row - 1 ).
                    last_pos_insertion_fraction_exceeds_threshold = false;
                  }
                  if(
                     ( m_profile->length() < max_profile_length ) // if we should consider adding a new position after this one..
                  ) {
                    if( count_seqs_exceeding_alignment_profile_thresholds ) {
                      if( m_row_i == last_row ) {
                        if( cout_indel_fractions ) {
                          cout << "\t Num seqs exceeding last pos insertion threshold: " << num_seqs_exceeding_last_pos_insertion_threshold << endl;
                          cout << "\t\t fraction is " << ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / ( double )m_sequence_count ) << endl;
                        }
                        // TODO: REMOVE.
                        if( ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / ( double )m_sequence_count ) > largest_insertion_fraction ) {
                          largest_insertion_fraction = ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / ( double )m_sequence_count );
                          largest_insertion_fraction_row = m_row_i;
                        }
                        // Should we try adding a profile position?
                        if(
                            ( ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / ( double )m_sequence_count ) > ( double )fraction_seqs_exceeding_insertion_threshold_insert_threshold )
                           ) {
                          last_pos_insertion_fraction_exceeds_threshold = true;
                          if(
                             use_sensitive_threshold_increments &&
                             (
                              last_pos_insertion_fraction <
                              min_last_pos_insertion_fraction_passing_threshold
                              )
                          ) {
                            min_last_pos_insertion_fraction_passing_threshold =
                              last_pos_insertion_fraction;
                          }
                        } // End checking if we can insert here.
                      } else { // m_row_i < last_row
                        if( cout_indel_fractions ) {
                          cout << "\t Num seqs exceeding insertion threshold: " << num_seqs_exceeding_insertion_threshold << endl;
                          cout << "\t\t fraction is " << ( ( double )num_seqs_exceeding_insertion_threshold / ( double )m_sequence_count ) << endl;
                        }
                        // TODO: REMOVE.
                        if( ( ( double )num_seqs_exceeding_insertion_threshold / ( double )m_sequence_count ) > largest_insertion_fraction ) {
                          largest_insertion_fraction = ( ( double )num_seqs_exceeding_insertion_threshold / ( double )m_sequence_count );
                          largest_insertion_fraction_row = m_row_i;
                        }
                        // Should we try adding a profile position?
                        if(
                            ( ( ( double )num_seqs_exceeding_insertion_threshold / ( double )m_sequence_count ) > ( double )fraction_seqs_exceeding_insertion_threshold_insert_threshold )
                           ) {
                          insertion_fraction_exceeds_threshold = true;
                          if(
                             use_sensitive_threshold_increments &&
                             (
                              insertion_fraction <
                              min_insertion_fraction_passing_threshold
                              )
                          ) {
                            min_insertion_fraction_passing_threshold =
                              insertion_fraction;
                          }
                        } // End checking if we can insert here.
                      } // End if m_row_i == last_row .. else ..
                    } else { // if count_seqs_exceeding_alignment_profile_thresholds .. else ..
                      // Then instead use the alignment_profile_position_allseqs 
                
                      if( m_row_i == 0 ) {
                        // Note that we're only using the insertion _open_, which we've tricksily stored separately from the insertion extension.
                        // Sanity check.  Since there are m_sequence_count sequences, each of which has a 1.0 in its B->M cell, the total of B->M should be that count.
                        //assert( ( alignment_profile_position_allseqs[ Transition::fromBegin ].total() / alignment_profile_position_allseqs.m_scalar ) == m_sequence_count );
                        //assert( abs( toDouble( alignment_profile_position_allseqs[ Transition::fromBegin ].total() / alignment_profile_position_allseqs.m_scalar ) - m_sequence_count ) < 1E-5 ); // MAGIC #!
                        insertion_fraction =
                          toDouble( ( alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position_allseqs.m_scalar * m_sequence_count ) );
                        if( cout_indel_fractions ) {
                          cout << "[ " << m_row_i << " ]  Pre-align insertion fraction: " << insertion_fraction << endl;
                          cout << "[ " << m_row_i << " ] alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs << endl;
                        }
                        // TODO: REMOVE
                        assert( insertion_fraction >= -1E-5 );
                        assert( ( 1.0 - insertion_fraction ) > -1E-5 );
                      } else if( m_row_i == last_row ) {
                        // Note that we're only using the insertion _open_, which we've tricksily stored separately from the insertion extension.
                        // Sanity check.  Since there are m_sequence_count sequences, each of which has a 1.0 in its C->T cell, the total of C->T should be that count.
                        // TODO: REMOVE
                        //double tmp_double = toDouble( alignment_profile_position_allseqs[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] / alignment_profile_position_allseqs.m_scalar );
                        //cout << "tmp_double: " << tmp_double << endl;
                        //assert( abs( toDouble( alignment_profile_position_allseqs[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] / alignment_profile_position_allseqs.m_scalar ) - m_sequence_count ) < 1E-5 ); // MAGIC #!
                        //assert( ( alignment_profile_position_allseqs[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] / alignment_profile_position_allseqs.m_scalar ) == m_sequence_count );
                        last_pos_insertion_fraction =
                          toDouble( ( alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position_allseqs.m_scalar * m_sequence_count ) );
                          //toDouble( ( alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position_allseqs[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] ) );
                        if( cout_indel_fractions ) {
                          cout << "[ " << m_row_i << " ]  Post-align insertion fraction: " << last_pos_insertion_fraction << endl;
                          cout << "[ " << m_row_i << " ] alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs << endl;
                        }
                        // TODO: REMOVE
                        assert( last_pos_insertion_fraction >= -1E-5 );
                        assert( ( 1.0 - last_pos_insertion_fraction ) > -1E-5 );
                      } else {
                        insertion_fraction =
                          // Denom should include paths that skip this position.
                          toDouble( ( alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position_allseqs[ Transition::fromMatch ].total() + alignment_profile_position_allseqs[ Transition::fromDeletion ].total() ) );
                          //toDouble( ( alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] ) / ( alignment_profile_position_allseqs[ Transition::fromMatch ].total() ) );
                        if( cout_indel_fractions ) {
                          cout << "[ " << m_row_i << " ]  Insertion fraction: " << insertion_fraction << endl;
                          //cout << "[ " << m_row_i << " ] alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs << endl;
                        }
                        // TODO: REMOVE
                        assert( insertion_fraction >= -1E-5 );
                        assert( ( 1.0 - insertion_fraction ) > -1E-5 );
                      }
                
                      if( m_row_i == last_row ) {
                        // TODO: REMOVE.
                        if( last_pos_insertion_fraction > largest_insertion_fraction ) {
                          largest_insertion_fraction = last_pos_insertion_fraction;
                          largest_insertion_fraction_row = m_row_i;
                        }
                        //cout << "[ " << m_row_i << " ]  Last pos insertion fraction: " << last_pos_insertion_fraction << endl;
                        // Should we try adding a profile position?
                        if(
                           last_pos_insertion_fraction >
                           propose_inserting_postalign_threshold
                        ) {
                          last_pos_insertion_fraction_exceeds_threshold = true;
                          if(
                             use_sensitive_threshold_increments &&
                             (
                              last_pos_insertion_fraction <
                              min_last_pos_insertion_fraction_passing_threshold
                              )
                          ) {
                            min_last_pos_insertion_fraction_passing_threshold =
                              last_pos_insertion_fraction;
                          }
                        } // End checking if we can insert here.
                      } else { // if m_row_i == last_row .. else ..
                        // TODO: REMOVE.
                        if( insertion_fraction > largest_insertion_fraction ) {
                          largest_insertion_fraction = insertion_fraction;
                          largest_insertion_fraction_row = m_row_i;
                        }
                        //cout << "[ " << m_row_i << " ]  Insertion fraction: " << insertion_fraction << endl;
                        // Should we try adding a profile position?
                        if(
                           insertion_fraction >
                           ( ( m_row_i == 0 ) ? propose_inserting_prealign_threshold : propose_inserting_threshold )
                        ) {
                          insertion_fraction_exceeds_threshold = true;
                          if(
                             use_sensitive_threshold_increments &&
                             (
                              insertion_fraction <
                              min_insertion_fraction_passing_threshold
                              )
                          ) {
                            min_insertion_fraction_passing_threshold =
                              insertion_fraction;
                          }
                        } // End checking if we can insert here.
                      } // End if m_row_i == last_row .. else ..

                    } // End if count_seqs_exceeding_alignment_profile_thresholds .. else ..
  
                    // Wait to do last row indels until last_row - 1.
                    // NOTE: Last row changes take precedence over
                    // second-to-last-row changes.
                    if( m_row_i == last_row ) {
                      if( last_pos_insertion_fraction_exceeds_threshold ) {
                        // Create backup of the alignment profile position so we
                        // can use it when inserting the new last position...
                        lastpos_alignment_profile_position_allseqs =
                          alignment_profile_position_allseqs;
                      }
                    }
  
                  } // End if ( m_profile->length() < max_profile_length ) (if we should consider adding a new position here ..)
  
                  deletion_fraction_exceeds_threshold = false;
                  if( m_row_i != ( last_row - 1 ) ) {
                    // We need to remember what we learned from last_row when
                    // m_row_i is ( last_row - 1 ).
                    last_pos_deletion_fraction_exceeds_threshold = false;
                  }
                  if( m_profile->length() > min_profile_length ) { // Should we bother considering deleting this position?
                    // In the last row we store deletions in the Match state,
                    // so we have to wait until the row is ( last_row - 1 )
                    // to find out when we deleted there.
                    if(
                       ( m_row_i < last_row ) &&
                       ( m_row_i > 0 ) // also we can't delete "row 0"
                     ) {
                      if( count_seqs_exceeding_alignment_profile_thresholds ) {
                        if( m_row_i == ( last_row - 1 ) ) {
                          if( cout_indel_fractions ) {
                            cout << "\t Num seqs exceeding (last pos) deletion threshold: " << num_seqs_exceeding_last_pos_deletion_threshold << endl;
                            cout << "\t\t fraction is " << ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / m_sequence_count ) << endl;
                          }
                          // TODO: REMOVE.
                          if( ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / ( double )m_sequence_count ) > largest_deletion_fraction ) {
                            largest_deletion_fraction = ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / ( double )m_sequence_count );
                            largest_deletion_fraction_row = last_row;
                          }

                          if(
                              ( ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / ( double )m_sequence_count ) > ( double )fraction_seqs_exceeding_deletion_threshold_delete_threshold )
                          ) {
                            // TODO: REMOVE
                            //cout << "DELETING LAST POS" << endl;
                            last_pos_deletion_fraction_exceeds_threshold = true;
                          } // End if last_pos_deletion_fraction > propose_deleting_threshold
                        } // End if ( m_row_i == ( last_row - 1 ) )
                      
                        if( cout_indel_fractions ) {
                          cout << "\t Num seqs exceeding deletion threshold: " << num_seqs_exceeding_deletion_threshold << endl;
                          cout << "\t\t fraction is " << ( ( double )num_seqs_exceeding_deletion_threshold / m_sequence_count ) << endl;
                        }
                        // TODO: REMOVE.
                        if( ( ( double )num_seqs_exceeding_deletion_threshold / ( double )m_sequence_count ) > largest_deletion_fraction ) {
                          largest_deletion_fraction = ( ( double )num_seqs_exceeding_deletion_threshold / ( double )m_sequence_count );
                          largest_deletion_fraction_row = m_row_i;
                        }
                        if(
                            ( ( ( double )num_seqs_exceeding_deletion_threshold / ( double )m_sequence_count ) > ( double )fraction_seqs_exceeding_deletion_threshold_delete_threshold )
                        ) {
                          deletion_fraction_exceeds_threshold = true;
                        } // End if deletion_fraction > propose_deleting_threshold
                      } else { // if count_seqs_exceeding_alignment_profile_thresholds .. else ..
                        if( m_row_i == ( last_row - 1 ) ) {
                          // First check up on deletions of the final position,
                          // since we couldn't do it when m_row_i was last_row
                          // (because we don't store deletions in the deletion
                          // state there).
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          last_pos_deletion_fraction =
                            toDouble( (  alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] + alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] + alignment_profile_position_allseqs[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] + alignment_profile_position_allseqs[ Transition::fromDeletionOut ].total() ) / ( ( alignment_profile_position_allseqs[ Transition::fromMatch ].total() + alignment_profile_position_allseqs[ Transition::fromDeletion ].total() + alignment_profile_position_allseqs[ Transition::fromDeletionOut ].total() ) ) );
#else // !USE_DEL_IN_DEL_OUT
                            last_pos_deletion_fraction =
                              toDouble( (  alignment_profile_position_allseqs[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] + alignment_profile_position_allseqs[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] ) / ( ( alignment_profile_position_allseqs[ Transition::fromMatch ].total() + alignment_profile_position_allseqs[ Transition::fromDeletion ].total() ) ) );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                            // TODO: REMOVE.
                            if( last_pos_deletion_fraction > largest_deletion_fraction ) {
                              largest_deletion_fraction = last_pos_deletion_fraction;
                              largest_deletion_fraction_row = last_row;
                            }
                            if( cout_indel_fractions ) {
                              cout << "[ " << m_row_i << " ] (last pos) Deletion fraction: " << last_pos_deletion_fraction << endl;
                            }
                            if(
                               ( last_pos_deletion_fraction > propose_deleting_threshold )
                               ) {
                              // TODO: REMOVE
                              //cout << "DELETING LAST POS" << endl;
                              last_pos_deletion_fraction_exceeds_threshold = true;
                            } // End if last_pos_deletion_fraction > propose_deleting_threshold
                          if(
                             use_sensitive_threshold_increments &&
                             (
                              last_pos_deletion_fraction <
                              min_last_pos_deletion_fraction_passing_threshold
                              )
                          ) {
                            min_last_pos_deletion_fraction_passing_threshold =
                              last_pos_deletion_fraction;
                          }

                        } // End if ( m_row_i == ( last_row - 1 ) )
                    
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        deletion_fraction =
                          toDouble( ( alignment_profile_position_allseqs[ Transition::fromDeletion ].total() + alignment_profile_position_allseqs[ Transition::fromDeletionIn ].total() + alignment_profile_position_allseqs[ Transition::fromDeletionOut ].total() ) / ( ( alignment_profile_position_allseqs[ Transition::fromMatch ].total() + alignment_profile_position_allseqs[ Transition::fromDeletion ].total() + alignment_profile_position_allseqs[ Transition::fromDeletionIn ].total() + alignment_profile_position_allseqs[ Transition::fromDeletionOut ].total() ) ) );
#else // !USE_DEL_IN_DEL_OUT
                        deletion_fraction =
                          toDouble( ( alignment_profile_position_allseqs[ Transition::fromDeletion ].total() ) / ( ( alignment_profile_position_allseqs[ Transition::fromMatch ].total() + alignment_profile_position_allseqs[ Transition::fromDeletion ].total() ) ) );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..

                        // TODO: REMOVE.
                        if( deletion_fraction > largest_deletion_fraction ) {
                          largest_deletion_fraction = deletion_fraction;
                          largest_deletion_fraction_row = m_row_i;
                        }
  
                        if( cout_indel_fractions ) {
                          cout << "[ " << m_row_i << " ] Deletion fraction: " << deletion_fraction << endl;
                        }
                        if(
                           ( deletion_fraction > propose_deleting_threshold )
                           ) {
                          deletion_fraction_exceeds_threshold = true;

                          if(
                             use_sensitive_threshold_increments &&
                             (
                              deletion_fraction <
                              min_deletion_fraction_passing_threshold
                              )
                          ) {
                            min_deletion_fraction_passing_threshold =
                              deletion_fraction;
                          }

                        } // End if deletion_fraction > propose_deleting_threshold
                      } // End if count_seqs_exceeding_alignment_profile_thresholds .. else ..
                    } // ( m_row_i < last_row ) && ( m_row_i > 0 )
                  } // End if m_profile->length() > min_profile_length
                  
                  if( m_row_i != last_row ) { // postpone length changes until row last_pos - 1


//                    /// TODO: REMOVE
//                    if( m_row_i == ( last_row - 1 ) ) {
//                      last_pos_deletion_fraction_exceeds_threshold = true;
//                      last_pos_insertion_fraction_exceeds_threshold = false;//true;
//                      enable_reverting = false;
//                    }

                    assert( last_pos_insertion_fraction_exceeds_threshold ? ( m_row_i == ( last_row - 1 ) ) : true );
                    assert( last_pos_deletion_fraction_exceeds_threshold ? ( m_row_i == ( last_row - 1 ) ) : true );
                    // Are we dealing with the last pos?  It always takes precedent.
                    if(
                      ( m_row_i == ( last_row - 1 ) ) &&
                      (
                        last_pos_insertion_fraction_exceeds_threshold ||
                        last_pos_deletion_fraction_exceeds_threshold
                      )
                    ) {
                      insertion_fraction_exceeds_threshold = false;
                      deletion_fraction_exceeds_threshold = false;
                    }
                    
                    // Set up the *_pos_i values (the position-to-be-removed or
                    // the position-to-be-inserted)
                    if(
                      ( m_row_i == ( last_row - 1 ) ) &&
                      (
                        last_pos_insertion_fraction_exceeds_threshold
                      )
                    ) {
                      propose_inserting_pos_i = last_row;
                      propose_inserting_rpos_i = numeric_limits<uint32_t>::max();
                    } else if( insertion_fraction_exceeds_threshold ) {
                      propose_inserting_pos_i = m_row_i;
                      propose_inserting_rpos_i = 
                        ( m_profile->length() - propose_inserting_pos_i );
                    } else {
                      // Unused if nothing exceeds threshold.
                      propose_inserting_pos_i = numeric_limits<uint32_t>::max();
                      propose_inserting_rpos_i = numeric_limits<uint32_t>::max();
                    }
                    if(
                      ( m_row_i == ( last_row - 1 ) ) &&
                      (
                        last_pos_deletion_fraction_exceeds_threshold
                      )
                    ) {
                      propose_deleting_pos_i = ( last_row - 1 );
                      propose_deleting_rpos_i = 0;
                    } else if( deletion_fraction_exceeds_threshold ) {
                      propose_deleting_pos_i = ( m_row_i - 1 );
                      propose_deleting_rpos_i = 
                        ( m_profile->length() - 1 - propose_deleting_pos_i );
                    } else {
                      // Unused if nothing exceeds threshold.
                      propose_deleting_pos_i = numeric_limits<uint32_t>::max();
                      propose_deleting_rpos_i = numeric_limits<uint32_t>::max();
                    }
                    
                    // Check for cycles
                    cycle_detected_insertions = false;
                    cycle_detected_deletions = false;
                    fancy_cycle_detected_insertions = false;
                    fancy_cycle_detected_deletions = false;
                    if(
                       ( ( m_row_i == ( last_row - 1 ) ) && last_pos_deletion_fraction_exceeds_threshold ) ||
                      deletion_fraction_exceeds_threshold
                    ) {
                      if(
                         ( ( m_row_i == ( last_row - 1 ) ) && last_pos_insertion_fraction_exceeds_threshold ) ||
                        insertion_fraction_exceeds_threshold
                      ) {
                        // If we're asked to *both* insert and delete here,
                        // call it a cycle.
                        // TODO: Increase the threshold by more in this instance?
                        cycle_detected_deletions = true;
                        cycle_detected_insertions = true;
                      } // End if both ins and del are called for.
                      if(
                        (
                          // We deleted this same position the last time we
                          // deleted something -- and the length is the same as
                          // the last time we deleted it.
                          (
                            have_deleted_a_pos &&
                            (
                              ( last_deleted_pos_i == propose_deleting_pos_i ) ||
                              ( last_deleted_rpos_i == propose_deleting_rpos_i )
                            )
                          ) &&
                          ( m_profile->length() == length_before_last_deletion )
                        ) ||
                        // Disallow deletion of just-inserted positions (give
                        // 'em a chance!).  By not allowing it, we force a wait
                        // until the threshold increases to the point where
                        // inserting it doesn't result in an immediate
                        // deletion.  If it is the 0th position, and if we just
                        // inserted it, then it was a pre-align insertion.
                        // Since it is special (there is no Match position
                        // there so there's opportunity for it to have been
                        // deleted earlier), we let it be immediately deleted,
                        // but we take note to avoid forever inserting and
                        // removing it.
                        (
                          //just_grew_profile &&
                          last_length_change_was_insertion &&
                          (
                            ( last_inserted_pos_i == propose_deleting_pos_i ) ||
                            ( last_inserted_rpos_i == propose_deleting_rpos_i )
                          ) &&
                          (
                            ( propose_deleting_pos_i > 0 ) ||
                            ( lengthadjust_position_creation_iters[ propose_deleting_pos_i ] != ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) ) // if we're about to delete the pos 0 that we inserted in a previous iteration, call it a cycle too.  it's only not a cycle (yet) if we inserted pos 0 on this iteration, and are now supposed to delete it.
                          )
                        ) ||
                        (
                          // Don't delete something that was created this
                          // iteration, (unless it's the pre-align pos).
                          ( lengthadjust_position_creation_iters[ propose_deleting_pos_i ] == ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) ) &&
                          ( propose_deleting_pos_i > 0 )
                        )
                      ) {
                        cycle_detected_deletions = true;
                        fancy_cycle_detected_deletions = true;
                    
                        // We should never be calling it a cycle if we're
                        // already at max threshold, right?
                        // TODO: REMOVE.
                        assert( propose_deleting_threshold <= 1.0 );
                      }
                    } // End if last_pos_deletion_fraction_exceeds_threshold || deletion_fraction_exceeds_threshold
                    if(
                       //!cycle_detected_insertions &&
                       (
                        ( m_row_i == ( last_row - 1 ) && last_pos_insertion_fraction_exceeds_threshold ) ||
                        insertion_fraction_exceeds_threshold
                      )
                    ) {
                      if(
                        just_caught_preAlign_cycle ||
                        (
                          // We inserted this same position the last time we
                          // inserted something -- and the length was the same
                          // then as it is now.
                          (
                            have_inserted_a_pos &&
                            (
                              ( last_inserted_pos_i == propose_inserting_pos_i ) ||
                              ( last_inserted_rpos_i == propose_inserting_rpos_i )
                            )
                          ) &&
                          ( m_profile->length() == length_before_last_insertion )
                        ) ||
                        // Disallow insertion of just-deleted positions.  By
                        // not allowing it, we force a wait until the threshold
                        // increases to the point where deleting it doesn't
                        // result in an immediate insertion.
                        (
                          //just_shrunk_profile &&
                          last_length_change_was_deletion &&
                          (
                            ( last_deleted_pos_i == propose_inserting_pos_i ) ||
                            ( last_deleted_rpos_i == propose_inserting_rpos_i )
                          )
                        )
                      ) {
                        cycle_detected_insertions = true;
                        fancy_cycle_detected_insertions = true;
                    
                        // We should never be calling it a cycle if we're
                        // already at max threshold, right?
                        // TODO: REMOVE.
                        assert( ( ( propose_inserting_pos_i == 0 ) ? propose_inserting_prealign_threshold : ( ( propose_inserting_pos_i == last_row ) ? propose_inserting_postalign_threshold : propose_inserting_threshold ) ) <= 1.0 );
                      } // End if( last_inserted_pos_i == propose_inserting_pos_i )
                    } // End if last_pos_insertion_fraction_exceeds_threshold || insertion_fraction_exceeds_threshold
                    
                    if( cycle_detected_insertions || cycle_detected_deletions ) {
                      if( !just_caught_preAlign_cycle ) {
                        if( cout_profile_length_changes ) {
                            cout << "Cycle detected.";
                          if( cycle_detected_insertions ) {
                            cout << "  Not inserting";
                          }
                          if( cycle_detected_deletions ) {
                            if( cycle_detected_insertions ) {
                              cout << " and not deleting";
                            } else {
                              cout << "  Not deleting";
                            }
                          }
                          cout << "." << endl;
                        } // End if cout_profile_length_changes
                        if(
                          (
                            ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                            ( m_trainingParameters.verbosity < VERBOSITY_High )
                          ) ||
                          ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                        ) {
                          cout << "cycle@";
                          if( cycle_detected_insertions ) {
                            cout << "i" << propose_inserting_pos_i;
                          }
                          if( cycle_detected_deletions ) {
                            if( cycle_detected_insertions ) {
                              cout << "&";
                            }
                            cout << "d" << propose_deleting_pos_i;
                          }
                          cout << "!";

                          if( count_seqs_exceeding_alignment_profile_thresholds ) {
                            if( last_pos_deletion_fraction_exceeds_threshold ) {
                              cout << "D:" << ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / m_sequence_count );
                            }
                            if( deletion_fraction_exceeds_threshold ) {
                              cout << "d:" << ( ( double )num_seqs_exceeding_deletion_threshold / m_sequence_count );
                            }
                            if( last_pos_insertion_fraction_exceeds_threshold ) {
                              cout << "I:" << ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / m_sequence_count );
                            }
                            if( insertion_fraction_exceeds_threshold ) {
                              cout << "i:" << ( ( double )num_seqs_exceeding_insertion_threshold / m_sequence_count );
                            }
                          } else { // if count_seqs_exceeding_alignment_profile_thresholds .. else ..
                            if( last_pos_deletion_fraction_exceeds_threshold ) {
                              cout << "D:" << last_pos_deletion_fraction;
                            }
                            if( deletion_fraction_exceeds_threshold ) {
                              cout << "d:" << deletion_fraction;
                            }
                            if( last_pos_insertion_fraction_exceeds_threshold ) {
                              cout << "I:" << last_pos_insertion_fraction;
                            }
                            if( insertion_fraction_exceeds_threshold ) {
                              cout << "i:" << insertion_fraction;
                            }
                          } // End if count_seqs_exceeding_alignment_profile_thresholds .. else ..
                        } // End if verbose  
                        if( use_sensitive_threshold_increments ) {
                          if(
                            cycle_detected_insertions &&
                            last_pos_insertion_fraction_exceeds_threshold &&
                            (
                              min_last_pos_insertion_fraction_passing_threshold <
                              min_insertion_fraction_passing_threshold_causing_a_cycle
                            )
                          ) {
                            min_insertion_fraction_passing_threshold_causing_a_cycle =
                              min_last_pos_insertion_fraction_passing_threshold;
                            if(
                              (
                                ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                                ( m_trainingParameters.verbosity < VERBOSITY_High )
                              ) ||
                              ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                            ) {
                              cout << "minCycleThreshold@I:" << min_insertion_fraction_passing_threshold_causing_a_cycle;
                            } // End if verbose
                          } else if(
                            cycle_detected_insertions &&
                            insertion_fraction_exceeds_threshold &&
                            (
                              min_insertion_fraction_passing_threshold <
                              min_insertion_fraction_passing_threshold_causing_a_cycle
                            )
                          ) {
                            min_insertion_fraction_passing_threshold_causing_a_cycle =
                              min_insertion_fraction_passing_threshold;
                            if(
                              (
                                ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                                ( m_trainingParameters.verbosity < VERBOSITY_High )
                              ) ||
                              ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                            ) {
                              cout << "minCycleThreshold@i:" << min_insertion_fraction_passing_threshold_causing_a_cycle;
                            } // End if verbose
                          }
                          if( cycle_detected_deletions ) {
                            if(
                              last_pos_deletion_fraction_exceeds_threshold &&
                              (
                                min_last_pos_deletion_fraction_passing_threshold <
                                min_deletion_fraction_passing_threshold_causing_a_cycle
                              )
                            ) {
                              min_deletion_fraction_passing_threshold_causing_a_cycle =
                                min_last_pos_deletion_fraction_passing_threshold;
                              if(
                                (
                                  ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                                  ( m_trainingParameters.verbosity < VERBOSITY_High )
                                ) ||
                                ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                              ) {
                                cout << "minCycleThreshold@D:" << min_deletion_fraction_passing_threshold_causing_a_cycle;
                              } // End if verbose
                            } else if(
                              deletion_fraction_exceeds_threshold &&
                              (
                                min_deletion_fraction_passing_threshold <
                                min_deletion_fraction_passing_threshold_causing_a_cycle
                              )
                            ) {
                              min_deletion_fraction_passing_threshold_causing_a_cycle =
                                min_deletion_fraction_passing_threshold;
                              if(
                                (
                                  ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                                  ( m_trainingParameters.verbosity < VERBOSITY_High )
                                ) ||
                                ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                              ) {
                                cout << "minCycleThreshold@d:" << min_deletion_fraction_passing_threshold_causing_a_cycle;
                              } // End if verbose
                            }
                          } // End if cycle_detected_deletions
                        } else { // if use_sensitive_threshold_increments .. else ..
                          propose_inserting_threshold +=
                            propose_inserting_threshold_increment;
                          propose_inserting_prealign_threshold +=
                            propose_inserting_threshold_increment;
                          propose_inserting_postalign_threshold +=
                            propose_inserting_threshold_increment;
                          propose_deleting_threshold +=
                            propose_deleting_threshold_increment;
                          
                          if( cout_profile_length_changes ) {
                            cout << "  propose_deleting_threshold increased to " << propose_deleting_threshold << " and propose_inserting_threshold increased to " << propose_inserting_threshold << ".";
                          } // End if cout_profile_length_changes
                          if(
                            (
                              ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                              ( m_trainingParameters.verbosity < VERBOSITY_High )
                            ) ||
                            ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                          ) {
                            if(
                              ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
                              ( propose_inserting_threshold == propose_inserting_postalign_threshold )
                            ) {
                              if( propose_inserting_threshold == propose_deleting_threshold ) {
                                cout << "threshold:" << propose_inserting_threshold;
                              } else {
                                cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
                              }
                            } else {
                              cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
                            }
                          } // End if verbose
                        } // End if use_sensitive_threshold_increments .. else ..
                        if(
                          m_trainingParameters.verbosity >= VERBOSITY_Medium
                        ) {
                          cout << endl;
                        } else {
                          cout.flush();
                        }
                        cycle_detected_this_iteration = true;
                      } // End if !just_caught_preAlign_cycle
                    } //else // End if cycle_detected_insertions || cycle_detected_deletions
                    if( 
                    // TODO:REMOVE
                       !fancy_cycle_detected_deletions &&
                       ( just_caught_preAlign_cycle || !cycle_detected_deletions || ( cycle_detected_deletions && cycle_detected_insertions && ( count_seqs_exceeding_alignment_profile_thresholds ? ( ( last_pos_deletion_fraction_exceeds_threshold ? ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / m_sequence_count ) : ( ( double )num_seqs_exceeding_deletion_threshold / m_sequence_count ) ) > ( ( last_pos_insertion_fraction_exceeds_threshold ? ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / m_sequence_count ) : ( ( double )num_seqs_exceeding_insertion_threshold / m_sequence_count ) ) ) ) : ( ( last_pos_deletion_fraction_exceeds_threshold ? last_pos_deletion_fraction : deletion_fraction ) > ( ( last_pos_insertion_fraction_exceeds_threshold ? last_pos_insertion_fraction : insertion_fraction ) ) ) ) ) ) &&

                      (
                        last_pos_deletion_fraction_exceeds_threshold ||
                        deletion_fraction_exceeds_threshold
                      )
                       // TODO: PUT BACK.
                      //&&
                      //(
                      //  !(
                      //    profile_length_changed_this_iteration &&
                      //    last_length_change_was_insertion
                      //  ) ||
                      //  (
                      //    ( propose_deleting_pos_i == 0 ) && 
                      //    // redundant: ( last_inserted_pos_i == 0 ) &&
                      //    ( lengthadjust_position_creation_iters[ 0 ] == ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) )
                      //  )
                      //)
                    ) {
                      have_deleted_a_pos = true;
                      if(
                        ( propose_deleting_pos_i == 0 ) && 
                        // redundant: ( last_inserted_pos_i == 0 ) &&
                        ( lengthadjust_position_creation_iters[ 0 ] == ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) )
                      ) {
                        // Special case of deleting the just-inserted preAlign
                        // position: This is a cycle, though we couldn't catch
                        // it through the usual channels, and we need to go
                        // ahead and do the deletion since we did the
                        // insertion.
                        just_caught_preAlign_cycle = true;
                        if( cout_profile_length_changes ) {
                          cout << "Cycle detected.";
                          cout << "  Deleting just-inserted preAlign position.";
                        } // End if cout_profile_length_changes
                        if(
                          (
                            ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                            ( m_trainingParameters.verbosity < VERBOSITY_High )
                          ) ||
                          ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                        ) {
                          cout << "&d0";
                          // TODO: REMOVE
                          cout << "@";
                          if( deletion_fraction_exceeds_threshold ) {
                            if( count_seqs_exceeding_alignment_profile_thresholds ) {
                              cout << "d:" << ( ( double )num_seqs_exceeding_deletion_threshold / m_sequence_count );
                            } else {
                              cout << "d:" << deletion_fraction;
                            }
                          }
                          cout << "cycle!";
                        } // End if verbose
                        if( use_sensitive_threshold_increments ) {
                          if(
                            insertion_fraction_exceeds_threshold &&
                            (
                              min_insertion_fraction_passing_threshold <
                              min_insertion_fraction_passing_threshold_causing_a_cycle
                            )
                          ) {
                            min_insertion_fraction_passing_threshold_causing_a_cycle =
                              min_insertion_fraction_passing_threshold;
                            if(
                              (
                                ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                                ( m_trainingParameters.verbosity < VERBOSITY_High )
                              ) ||
                              ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                            ) {
                              cout << "minCycleThreshold@i:" << min_insertion_fraction_passing_threshold_causing_a_cycle;
                            } // End if verbose
                          }
                          if(
                            deletion_fraction_exceeds_threshold &&
                            (
                              min_deletion_fraction_passing_threshold <
                              min_deletion_fraction_passing_threshold_causing_a_cycle
                            )
                          ) {
                            min_deletion_fraction_passing_threshold_causing_a_cycle =
                              min_deletion_fraction_passing_threshold;
                            if(
                              (
                                ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                                ( m_trainingParameters.verbosity < VERBOSITY_High )
                              ) ||
                              ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                            ) {
                              cout << "minCycleThreshold@d:" << min_deletion_fraction_passing_threshold_causing_a_cycle;
                            } // End if verbose
                          }
                        } else { // if use_sensitive_threshold_increments .. else ..
                          propose_inserting_threshold +=
                            propose_inserting_threshold_increment;
                          propose_inserting_prealign_threshold +=
                            propose_inserting_threshold_increment;
                          propose_inserting_postalign_threshold +=
                            propose_inserting_threshold_increment;
                          propose_deleting_threshold +=
                            propose_deleting_threshold_increment;
                    
                          if( cout_profile_length_changes ) {
                            cout << "  propose_deleting_threshold increased to " << propose_deleting_threshold << " and propose_inserting_threshold increased to " << propose_inserting_threshold << ".";
                          } // End if cout_profile_length_changes
                          if(
                            (
                              ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                              ( m_trainingParameters.verbosity < VERBOSITY_High )
                            ) ||
                            ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                          ) {
                            if(
                              ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
                              ( propose_inserting_threshold == propose_inserting_postalign_threshold )
                            ) {
                              if( propose_inserting_threshold == propose_deleting_threshold ) {
                                cout << "threshold:" << propose_inserting_threshold;
                              } else {
                                cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
                              }
                            } else {
                              cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
                            }
                          } // End if verbose
                        } // End if use_sensitive_threshold_increments .. else ..
                        if(
                          m_trainingParameters.verbosity >= VERBOSITY_Medium
                        ) {
                          cout << endl;
                        } else {
                          cout.flush();
                        }
                        cycle_detected_this_iteration = true;
                    
                        // Revert to previous last_inserted_pos_i,
                        // length_before_last_insertion,
                        // last_length_change_was_*, and
                        // profile_length_changed_this_iteration, since this is
                        // effectively a cycle we detected late, not a real
                        // change.
                        last_inserted_pos_i = prev_to_last_inserted_pos_i;
                        last_inserted_rpos_i = prev_to_last_inserted_rpos_i;
                        length_before_last_insertion = prev_length_before_last_insertion;
                        last_length_change_was_insertion =
                          prev_to_last_length_change_was_insertion;
                        last_length_change_was_deletion =
                          !prev_to_last_length_change_was_insertion;
                        profile_length_changed_this_iteration =
                          prev_profile_length_changed_this_iteration;
                    
                        // And keep the old last_deleted_pos_i stuff.
                        // And keep the old last_deleted_pos_i and
                        // length_before_last_deletion stuff.
                      } else { // if we're catching a pre-align cycle .. else ..
                        prev_to_last_deleted_pos_i = last_deleted_pos_i;
                        prev_to_last_deleted_rpos_i = last_deleted_rpos_i;
                        prev_length_before_last_deletion = length_before_last_deletion;
                        last_deleted_pos_i = propose_deleting_pos_i;
                        last_deleted_rpos_i = propose_deleting_rpos_i;
                        length_before_last_deletion = m_profile->length();
                      } // End if we're catching a pre-align cycle .. else ..
                      if( cout_profile_length_changes ) {
                        cout << "\t\tDeleting profile position " << propose_deleting_pos_i << endl;
                        if( propose_deleting_pos_i == ( last_row - 1 ) ) {
                          if( !cout_indel_fractions ) { // We've already printed it.
                            if( count_seqs_exceeding_alignment_profile_thresholds ) {
                              cout << "\t Num seqs exceeding (last pos) deletion threshold: " << num_seqs_exceeding_last_pos_deletion_threshold << endl;
                              cout << "\t\t fraction is " << ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / ( double )m_sequence_count ) << endl;
                            } else {
                              cout << "\t Last pos Deletion fraction: " << last_pos_deletion_fraction << endl;
                            }
                          }
                        } else {
                          if( !cout_indel_fractions ) { // We've already printed it.
                            if( count_seqs_exceeding_alignment_profile_thresholds ) {
                              cout << "\t Num seqs exceeding deletion threshold: " << num_seqs_exceeding_deletion_threshold << endl;
                              cout << "\t\t fraction is " << ( ( double )num_seqs_exceeding_deletion_threshold / ( double )m_sequence_count ) << endl;
                            } else {
                              cout << "\t Deletion fraction: " << deletion_fraction << endl;
                            }
                          }
                        }
                        //cout << "[ " << m_row_i << " ] alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs << endl;
                      } // End if cout_profile_length_changes
                      if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                        if( !just_caught_preAlign_cycle ) {
                          cout << "d" << propose_deleting_pos_i;
                          // TODO: REMOVE
                          cout << "@";
                          if( last_pos_deletion_fraction_exceeds_threshold ) {
                            if( count_seqs_exceeding_alignment_profile_thresholds ) {
                              cout << "D:" << ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / m_sequence_count );
                            } else {
                              cout << "D:" << last_pos_deletion_fraction;
                            }
                          }
                          if( deletion_fraction_exceeds_threshold ) {
                            if( count_seqs_exceeding_alignment_profile_thresholds ) {
                              cout << "d:" << ( ( double )num_seqs_exceeding_deletion_threshold / m_sequence_count );
                            } else {
                              cout << "d:" << deletion_fraction;
                            }
                          }
                        } // else we already printed that it's a cycle.
                        if(
                          m_trainingParameters.verbosity >= VERBOSITY_Medium
                        ) {
                          cout << endl;
                        } else {
                          cout.flush();
                        }
                      } // End if verbose
                    
                      // First shorten the profile.
                      //cout << "\t\tOld profile: " << *m_profile << endl;
                      m_profile->erase(
                        m_profile->begin() + propose_deleting_pos_i
                      );
                      //cout << "\t\tNew profile: " << *m_profile << endl;
                      if( cout_profile_length_changes ) {
                        cout << "\t\tNew profile length: " << m_profile->length() << endl;
                      }
                    
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                      // Also shorten the scaled match distributions.
                      // Technically the match distributions after this
                      // position don't need to be recalculated, so maybe TODO:
                      // optimize more.
                      m_scaled_match_distributions.resize( m_profile->length() - 1 );
                      for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                        m_profile->createScaledMatchDistributionForPosition(
                          tmp_pos_i,
                          m_scaled_match_distributions[ tmp_pos_i ]
                        );
                      } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          
                      if( m_trainingParameters.useUnconditionalBaumWelch ) {
                        m_unconditional_position_ententes_vector.erase(
                          m_unconditional_position_ententes_vector.begin() + propose_deleting_pos_i
                        );
#ifdef ALLOW_BOLTZMANN_GIBBS
                        if( m_trainingParameters.baldiLearningRate > 0 ) {
                          m_unconditional_baldi_position_boltzmann_gibbs_changes_vector.erase(
                            m_unconditional_baldi_position_boltzmann_gibbs_changes_vector.begin() + propose_deleting_pos_i
                          );
                        } // End if baldiLearningRate > 0
#endif // ALLOW_BOLTZMANN_GIBBS
                     } // End if useUnconditionalBaumWelch
                         
                     // Also resize the starting profiles? No.  They'll be
                     // reset when they need to be.
                         
                     // Resize the forward matrices?  Nah, it doesn't hurt
                     // to have them be too long..
                   
                     // UPDATE: If we ever do put the forward matrices back,
                     // then we should shrink them here.
                   
                     // Shrink the anchor columns.
                     if(
                       ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                       ( m_anchor_rows.m_storeEveryNthRow > 1 )
                     ) {
                       if( propose_deleting_pos_i == ( last_row - 1 ) ) {
                         // Remove the last element.
                         anchor_columns_iter = m_anchor_columns.end();
                         anchor_columns_iter--;
                   
                         m_anchor_columns.erase(
                           anchor_columns_iter
                         );
                         // Keep the iterator up-to-date
                         anchor_columns_iter =
                           m_anchor_columns.end();
                         if( cout_anchor_columns_iter_changes ) {
                           anchor_columns_iter_i = m_anchor_columns.size();
                           cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                         }
                       } else { // if propose_deleting_pos_i == ( last_row - 1 ) .. else ..
                         // Remove the last row of anchor columns,
                         // since it doesn't matter which one we
                         // delete.  (All rows of anchor columns after
                         // the one associated with the deleted
                         // position are bogus anyway.)
                   
                         tmp_matrices_iter = m_anchor_columns.end();
                         tmp_matrices_iter--;
                         assert( anchor_columns_iter != tmp_matrices_iter );
                   
                         m_anchor_columns.erase(
                           anchor_columns_iter
                         );
                         // The iterator was just invalidated..
                   
                         anchor_columns_iter = m_anchor_columns.begin();
                         advance( anchor_columns_iter, ( propose_deleting_pos_i + 1 ) );
                         if( cout_anchor_columns_iter_changes ) {
                           anchor_columns_iter_i = propose_deleting_pos_i + 1;
                           cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                         }
                       } // End if propose_deleting_pos_i == ( last_row - 1 ) .. else ..
                       // After deleting, the anchor_columns_iter points to
                       // the RowVector just after the one we've
                       // deleted.
                   
                       // That one -- and all subsequent RowVectors -- are
                       // bogus now, but it doesn't matter since we're moving
                       // backwards through the rows.
                     }  // End if m_anchor_columns.m_storeEveryNthColumn > 0
                   
                     // Shrink the anchor rows if there's one fewer now.
                     if(
                       // We remove an anchor row if the new length changes
                       // the number of anchor rows: that is, if the old last
                       // index, last_row, divides evenly by
                       // m_anchor_rows.m_storeEveryNthRow.
                       ( ( last_row % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                     ) {
                       if( propose_deleting_pos_i == ( last_row - 1 ) ) {
                         // Remove the last anchor row.
                         anchor_rows_iter =
                           m_anchor_rows.end();
                         anchor_rows_iter--;
                   
                         m_anchor_rows.erase(
                           anchor_rows_iter
                         );
                         // Point after the deleted position
                         anchor_rows_iter =
                           m_anchor_rows.end();
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i = m_anchor_rows.size();
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                       } else { // if propose_deleting_pos_i == ( m_row_i - 1 ) .. else ..
                         // Delete the last anchor row, since it doesn't
                         // matter which one we delete. (All anchor rows after
                         // the one associated with the deleted position are
                         // bogus anyway.)
                   
                         // TODO: REMOVE Unless we're removing the last
                         // position, the number of anchor rows shouldn't
                         // change when the current anchor row is the last
                         // anchor row (that is, if we get to this point, we
                         // know that the number of anchor rows is changing,
                         // so we know that it's not the last anchor row):
                         tmp_matrices_iter = m_anchor_rows.end();
                         tmp_matrices_iter--;
                         assert( anchor_rows_iter != tmp_matrices_iter );
                   
                         // Note that the anchor_rows_iter won't change.  We
                         // want to be pointing to the anchor row associated
                         // with the position just after the deleted position,
                         // which is always the current anchor_rows_iter: if
                         // the deleted pos was associated with an anchor row,
                         // then (despite that we actually are removing the
                         // _last_ anchor row) the new role for
                         // *anchor_rows_iter is to be the anchor row
                         // associated with the position after the removed
                         // one.
                         m_anchor_rows.erase(
                           tmp_matrices_iter
                         );
                         // The iterator was just invalidated..
                         anchor_rows_iter = m_anchor_rows.begin();
                         advance( anchor_rows_iter, ( propose_deleting_pos_i + 1 ) / m_anchor_rows.m_storeEveryNthRow );
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i = ( ( propose_deleting_pos_i + 1 ) / m_anchor_rows.m_storeEveryNthRow );
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                         // After deleting, the anchor_rows_iter points to the
                         // anchor row associated with the position just after
                         // the position-to-be-deleted .
  
                         // But see below where we increment the anchor_rows_iter...
                       } // End if propose_deleting_pos_i == ( m_row_i - 1 ) .. else ..
                     } else { // if we should remove an anchor row .. else ..
                       // Make sure that the anchor_rows_iter points to the
                       // anchor row associated with the position just after
                       // the position-to-be-deleted.
                       if( propose_deleting_pos_i == ( last_row - 1 ) ) {
                         anchor_rows_iter = m_anchor_rows.end();
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i = m_anchor_rows.size();
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                         /// NOTE that I've commented this out because it
                         /// doesn't really matter, since we would have to
                         /// undo it anyway, and we never use the
                         /// ancho_rows_iter in this circumstance, since we
                         /// only need to recalc the forward rows if we're
                         /// deleting the last row (see below) -- also when I
                         /// uncomment this and the corresponding undo code
                         /// (see below right before "just_shrunk_profile =
                         /// true"), it breaks everything.  So screw it.
                       //} else if( ( ( propose_deleting_pos_i + 2 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                       //  // The position after the one that we are deleting is
                       //  // itself associated with an anchor row, so it has a
                       //  // different anchor row than the one we're deleting.
                       //  // It must be the subsequent one.
                       //  anchor_rows_iter++;
                       //  if( cout_anchor_rows_iter_changes ) {
                       //    anchor_rows_iter_i++;
                       //    cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                       //  }
                       } // else it's the same as the current one, so do nothing.
                   
                     } // End if we should remove an anchor row .. else ..

  #if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        if( !do_not_recalculate_all_forward_rows_during_lengthadjust ) {
                          // Since we've changed the distance of the rows to the
                          // end, we need to recalculate all of the forward rows up
                          // to and including the one for the row just before
                          // propose_deleting_pos_i.
                          tmp_matrices_iter = m_anchor_rows.begin();
                          if(
                            ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                            ( m_anchor_rows.m_storeEveryNthRow > 1 )
                          ) {
                            tmp_matrices_iter2 = m_anchor_columns.begin();
                          }
                      
                          swap = false;
                          for( uint32_t row_i = 0; row_i <= propose_deleting_pos_i; row_i++ ) {
                            if(
                              ( row_i != 0 ) &&
                              ( ( row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                            ) {
                              tmp_matrices_iter++;
                              assert( tmp_matrices_iter != m_anchor_rows.end() );
                            } // End if storing anchor rows too and this is an anchor row
                            if(
                              ( row_i != 0 ) &&
                              ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                              ( m_anchor_rows.m_storeEveryNthRow > 1 )
                            ) {
                              tmp_matrices_iter2++;
                            }
                      
                            // Every other time, we must swap which forward row is current.
                            swap = !swap;
                      
                            for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                              m_dynamic_programming.forward_calculateRow(
                                m_trainingParameters,
                                false, // don't use viterbi
                                *m_profile,
    //#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                                (
                                  ( ( row_i <= 1 ) || ( ( row_i - 2 ) >= ( last_row - 1 ) ) ) ?
                                  m_scaled_match_distributions[ row_i ] : // ignored
                                  m_scaled_match_distributions[ row_i - 2 ]
                                ),
                                (
                                  ( ( row_i == 0 ) || ( ( row_i - 1 ) >= ( last_row - 1 ) ) ) ?
                                  m_scaled_match_distributions[ row_i ] : // ignored
                                  m_scaled_match_distributions[ row_i - 1 ]
                                ),
    //#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                                m_sequences[ m_seq_i ],
                                row_i,
                                ( swap ? ( *m_forward_rows_ptr )[ m_seq_i ] : ( *m_prev_forward_rows_ptr )[ m_seq_i ] ),
                                ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] )
                              );
                              if(
                                ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                                ( m_anchor_rows.m_storeEveryNthRow > 1 )
                              ) {
                                // This won't copy the rabiner scalars -- see below for that.
                                ( *tmp_matrices_iter2 )[ m_seq_i ].storeColumns(
                                  ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ),
                                  m_anchor_columns.m_storeEveryNthColumn
                                );
                              } // End if storing anchor columns too.
                              if(
                                ( m_anchor_rows.m_storeEveryNthRow != 0 ) &&
                                ( ( row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                              ) {
                                // This won't copy rabiner scalars.  See below for that.
                                ( *tmp_matrices_iter )[ m_seq_i ] =
                                  ( swap ? ( *m_prev_forward_rows_ptr ) : ( *m_forward_rows_ptr ) )[ m_seq_i ];
                              } // End if storing anchor rows too and this is an anchor row
                      
                              if( m_trainingParameters.useRabinerScaling ) {
                                if( row_i == 0 ) {
                                  ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerInverseScalar;
                                } else {
                                  ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar =
                                    (
                                      ( swap ? ( *m_forward_rows_ptr )[ m_seq_i ] : ( *m_prev_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar *
                                      ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerInverseScalar
                                    );
                                }
                                if(
                                  ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                                  ( m_anchor_rows.m_storeEveryNthRow > 1 )
                                ) {
                                  ( *tmp_matrices_iter2 )[ m_seq_i ].m_rabinerInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerInverseScalar;
                                  ( *tmp_matrices_iter2 )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar;
                                } // End if anchor_columns != 0;
                                if(
                                  ( m_anchor_rows.m_storeEveryNthRow != 0 ) &&
                                  ( ( row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                                ) {
                                  ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr ) : ( *m_forward_rows_ptr ) )[ m_seq_i ].m_rabinerInverseScalar;
                                  ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr ) : ( *m_forward_rows_ptr ) )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                                } // End if storing anchor rows too and this is an anchor row
                              } // End if useRabinerScaling
                      
                            } // End foreach sequence.
                          } // End foreach row up to the one just before the
                            // deleted one,, recalculate the forward rows.
                          if( swap ) {
                            // Rotate the forward_rows
                            m_temp_forward_rows_ptr = m_forward_rows_ptr;
                            m_forward_rows_ptr = m_prev_forward_rows_ptr;
                            m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                          }
                          forward_rows_are_already_rotated_and_calculated = true;
                        } else { // !do_not_recalculate_all_forward_rows_during_lengthadjust .. else ..
#else // !USE_DEL_IN_DEL_OUT
                     if( true ) {
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..                    
                        // But we may need to recalculate forward rows whose
                        // calculation depended on their row position.  This
                        // applies to the first, second, and last rows only;
                        // but of course we didn't delete at the first row (row
                        // 0).  So we really only have to do it if we just
                        // deleted the last position or the first position.
                        // Technically we don't even need to do it if it's the
                        // first pos, since it would get done anyway on the
                        // next iteration's forward row rotation.
                        if(
                          ( propose_deleting_pos_i == ( last_row - 1 ) )
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                            || true
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        ) {
                          // When deleting the last pos, the current row is
                          // actually one less than the last row.  Thus we
                          // already have the correct previous row to
                          // recalculate the current row, which is becoming the
                          // last row (since we are deleting the *old* last
                          // row).  Here we just recalculate the now-final row.
                          for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                            m_dynamic_programming.forward_calculateRow(
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              (
                                ( propose_deleting_pos_i == 1 ) ?
                                m_scaled_match_distributions[ 0 ] : // ignored
                                m_scaled_match_distributions[ propose_deleting_pos_i - 2 ]
                              ),
                              m_scaled_match_distributions[ propose_deleting_pos_i - 1 ],
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                               m_sequences[ m_seq_i ],
                               propose_deleting_pos_i,
                               ( *m_prev_forward_rows_ptr )[ m_seq_i ],
                               ( *m_forward_rows_ptr )[ m_seq_i ]
                             );
                             if( m_trainingParameters.useRabinerScaling ) {
                               // Note that we know that propose_deleting_pos_i is
                               // last_row - 1.
                               if( propose_deleting_pos_i == 0 ) {
                                 ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                               } else {
                                 ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                                 ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar *=
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                               } // End if propose_deleting_pos_i == 0 .. else ..
                               // NOTE: The backward rows' rabiner scalars are
                               // updated below.
                             } // End if useRabinerScaling
                             // Keep the anchor columns up-to-date
                             if(
                               ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                               ( m_anchor_rows.m_storeEveryNthRow > 1 )
                             ) {
                               if( propose_deleting_pos_i == ( last_row - 1 ) ) { 
                                 // Note that we know that propose_deleting_pos_i is
                                 // last_row - 1, which is what we set it to above
                                 // -- so the anchor_columns_iter points to the
                                 // right place.
                                 // TODO: REMOVE
                                 assert( anchor_columns_iter == m_anchor_columns.end() );
                               } else {
                                 // It points to just after the deleted row of anchor columns.
                               }
                               tmp_matrices_iter = anchor_columns_iter;
                               tmp_matrices_iter--;
                               ( *tmp_matrices_iter )[ m_seq_i ].storeColumns(
                                 ( *m_forward_rows_ptr )[ m_seq_i ],
                                 m_anchor_columns.m_storeEveryNthColumn
                               );
                               if( m_trainingParameters.useRabinerScaling ) {
                                 // Copy the rabiner scalars.
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                               } // End if useRabinerScaling
                             } // End if m_anchor_columns.m_storeEveryNthColumn > 0
                             // Keep the anchor rows up-to-date
                             if( ( propose_deleting_pos_i % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                               // Note that we know & utilize the fact that the
                               // anchor row associated with the now-final
                               // position is the final anchor row, always.
                     
                               // temporarily decrement it...
                               if( propose_deleting_pos_i == ( last_row - 1 ) ) {
                                 // TODO: REMOVE
                                 assert( anchor_rows_iter == m_anchor_rows.end() );
                               } else {
                                 // It points to just after the deleted row..
                               }
                               tmp_matrices_iter = anchor_rows_iter;
                               tmp_matrices_iter--;
                     
                               ( *tmp_matrices_iter )[ m_seq_i ] =
                                 ( *m_forward_rows_ptr )[ m_seq_i ];
                               if( m_trainingParameters.useRabinerScaling ) {
                                 // Copy the rabiner scalars.
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                               } // End if useRabinerScaling
                             } // End if this is an anchor row, store it.
                           } // End foreach m_seq_i
                           // At this moment, the m_forward_rows_ptr and
                           // m_prev_forward_rows_ptr are correctly set up for
                           // the next iteration.
                           forward_rows_are_already_rotated_and_calculated = true;
                         } // End if propose_deleting_pos_i == ( last_row - 1 ) (or maybe true)
                       } // End if true (or maybe if do_not_recalculate_all_forward_rows_during_lengthadjust)
                           
                       if( propose_deleting_pos_i == ( last_row - 1 ) ) {
                         // These won't be decremented, since m_row_i will be
                         // last_row the next time around.  Pointing to the
                         // last element is what we want.
                         assert( anchor_rows_iter == m_anchor_rows.end() );
                         anchor_rows_iter--;
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i = m_anchor_rows.size() - 1;
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                         if(
                           ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                           ( m_anchor_rows.m_storeEveryNthRow > 1 )
                         ) {
                           assert( anchor_columns_iter == m_anchor_columns.end() );
                           anchor_columns_iter--;
                           if( cout_anchor_columns_iter_changes ) {
                             anchor_columns_iter_i--;
                             cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                           }
                         }  // End if m_anchor_columns.m_storeEveryNthColumn > 0
                     
                         // Both of these will be decremented, but it'll end
                         // up doing the right thing: starting the whole
                         // position cycle over...
                         m_row_i = last_row;
                       } // End if propose_deleting_pos_i == ( last_row - 1 )
                     
                       // Lastly, we need to revert to the previous
                       // backward_rows..
                     
                       // We set this because we'll reuse the current backward_row and prev_backward_row (but see below)
                       backward_rows_are_already_rotated_and_calculated = true;

                       // Note that we don't have to recalculate any
                       // backward rows -- we can use the backward row that
                       // was calculated for ( m_row_i + 1 ), that we've
                       // just rotated into the m_backward_rows_ptr
                       // position.  The calculation of the current backward
                       // row doesn't depend on any profile values in the
                       // just-deleted position, or on the row location
                       // except as last_row or ( last_row - 1 ) --
                       // distances from the end -- which is unaffected by
                       // removing this position.  We're now treating it as
                       // the backward row for m_row_i.
                     
                       // The only time we need to *actually* recompute the backward
                       // row is when we just deleted the last position --
                       // because then there is nothing in the backward row
                       // (it's just garbage), but we need it (its last
                       // position -- and only if not using rabiner
                       // scaling!) to calculate the score using the
                       // forward_score(..) method, below.
                       if( ( propose_deleting_pos_i == 0 ) || ( propose_deleting_pos_i == ( last_row - 1 ) ) ) {
                         // We don't have to do it again next time 'round
                     
                         for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                           m_dynamic_programming.backward_calculateRow(
                             m_trainingParameters,
                             false, // don't use viterbi
                             *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                             (
                              ( propose_deleting_pos_i == 0 ) ?
                              m_scaled_match_distributions[ 0 ] : // ignored
                              m_scaled_match_distributions[ propose_deleting_pos_i - 1 ]
                             ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                             m_sequences[ m_seq_i ],
                             ( ( propose_deleting_pos_i + 1 ) - 1 ),
                             ( *m_next_backward_rows_ptr )[ m_seq_i ], // ignored when propose_deleting_pos_i == ( last_row - 1 )
                             ( *m_backward_rows_ptr )[ m_seq_i ]
                            );
                         } // End foreach seq_i
                       } // End if ( propose_deleting_pos_i == 0 || ( last_row - 1 ) )
                      
                      // If we don't alwaysAccept, then we need to
                      // recalculate the score, so we know if we have an
                      // improvement next time (Note that we can't reject the
                      // length change; this is for the next iteration...).
                    
                      // TODO: PUT BACK.  TESTING.
                      //if( !( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) ) {
                      if( true ) {
                        // TODO: REMOVE
                        //cout << "(1) Ending score, before recalculating, is " << m_endingScore << endl;
                        // The calculation of the current backward row
                        // depends on the *next* profile position, so we
                        // don't need to recalculate it after removing *this*
                        // profile position, before calculating the score.
                        m_endingScore =
                          m_dynamic_programming.forward_score( 
                            m_trainingParameters,
                            *m_profile,
                            m_sequence_count,
                            propose_deleting_pos_i, // row_i - 1 or last_row - 1
                            //m_forward_matrices[ ( propose_deleting_pos_i + 1 ) - 1 ],
                            (
                              forward_rows_are_already_rotated_and_calculated ?
                              *m_forward_rows_ptr :
                              *m_prev_forward_rows_ptr
                            ),
                            *m_backward_rows_ptr,
                            &sequence_scores,
                            &largest_sequence_score
                          );
                        // TODO: REMOVE
                        //cout << "(1) Ending score, calculated the new way, is " << m_endingScore << endl;
                    
                        // TODO: REMOVE?
                        assert( !isnan( m_endingScore ) );
                      } // End if !alwaysAccept, recalculate the score(s)
                    
                      if( !backward_rows_are_already_rotated_and_calculated ) {
                        // Rotate the backward rows to account for the
                        // deletion.  Note that they will be rotated again on
                        // the next iteration, so this is just enforcing that
                        // the next_backward_rows that were used to calculate
                        // the backward_rows will be used again on the next
                        // iteration.  From this point on the current
                        // backward_rows (hereafter called
                        // next_backward_rows) are garbage.
                        m_temp_backward_rows_ptr = m_backward_rows_ptr;
                        m_backward_rows_ptr = m_next_backward_rows_ptr;
                        m_next_backward_rows_ptr = m_temp_backward_rows_ptr;
                    
                        // TODO: Do I need to do this when propose_deleting_pos_i == 0 too?
                        if( m_trainingParameters.useRabinerScaling ) {
                          // The backward rows were scaled by the
                          // rabinerInverseScalar of the corresponding forward
                          // matrix rows, so now the scale factor needs to be
                          // updated.
                          ScoreType scale_ratio;
                          for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                            // This is the old scale:
                            scale_ratio =
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            // TODO: REMOVE. Note to self.
                            //pos_i ==
                            //( forward_rows_are_already_rotated_and_calculated ? ( last_row - 1 ) : ( m_row_i - 1 ) );
                            scale_ratio /=
                              // The new scale:
                              //m_forward_matrices[ ( propose_deleting_pos_i + 1 ) - 1 ][ m_seq_i ].m_rabinerInverseScalar;
                              //m_forward_matrices[ ( forward_rows_are_already_rotated_and_calculated ? m_row_i : ( m_row_i - 1 ) ) ][ m_seq_i ].m_rabinerInverseScalar;
                              ( forward_rows_are_already_rotated_and_calculated ? ( *m_forward_rows_ptr ) : ( *m_prev_forward_rows_ptr ) )[ m_seq_i ].m_rabinerInverseScalar;
                    
                            // TODO: Why not just change the scalar (and not
                            // actually scale the data): isn't this more likely
                            // to result in underflow or overflow?
                    
                            // Rescale it to the new scale.
                            ( *m_backward_rows_ptr )[ m_seq_i ] *= scale_ratio;
                            
                            // Note the new scale.
                            ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                            //////m_forward_matrices[ ( propose_deleting_pos_i + 1 ) - 1 ][ m_seq_i ].m_rabinerInverseScalar;
                              ( forward_rows_are_already_rotated_and_calculated ? ( *m_forward_rows_ptr ) : ( *m_prev_forward_rows_ptr ) )[ m_seq_i ].m_rabinerInverseScalar;
                            if( forward_rows_are_already_rotated_and_calculated ) {
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            } else {
                                ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar /=
                                scale_ratio;
                            } // End if forward_rows_are_already_rotated_and_calculated .. else ..
                          } // End foreach sequence, fix the backward row's rabinerInverseScalars.. 
                        } // End if useRabinerScaling
                    
                      } // End if !backward_rows_are_already_rotated_and_calculated
                    
                      // Keep track of the creation time. (In this case,
                      // un-keep track).
                      lengthadjust_position_creation_iters.erase(
                        lengthadjust_position_creation_iters.begin() + propose_deleting_pos_i
                      );

                      // Make sure that the anchor_rows_iter (&anchor_columns_iter) will point to the right thing after decrementing it next 'round.
                      //if(
                      //   ( m_row_i != last_row ) && // it won't be decremented if m_row_i == last_row.
                      //  // The anchor row is stored at the beginning of the block of
                      //  // rows, so if the next row is the first of a block, then we
                      //  // are in the block of a new anchor row, as we progress
                      //  // backwards.
                      //  ( ( ( m_row_i + 1 - 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                      //) {
                      //  anchor_rows_iter++;
                      //  if( cout_anchor_rows_iter_changes ) {
                      //    anchor_rows_iter_i++;
                      //    cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                      //  }
                      //}
                      //if(
                      //   ( m_row_i != last_row ) &&
                      //   ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                      //   ( m_anchor_rows.m_storeEveryNthRow > 1 )
                      //) {
                      //  anchor_columns_iter++;
                      //  if( cout_anchor_columns_iter_changes ) {
                      //    anchor_columns_iter_i++;
                      //    cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                      //  }
                      //} // End if ( not the last row ) and m_anchor_columns.m_storeEveryNthColumn > 0 

                      /// NOTE that I've commented this out.  See above (look for "screw it").
                      // Okay so at this point we've got the anchor_rows_iter
                      // pointing to the anchor row associated with the
                      // position after the one we just deleted.  We need it to
                      // point instead to the one associated with the position
                      // we've deleted, which will only be different if the pos
                      // after this one is the first in its block.
                      //if(
                      //  ( propose_deleting_pos_i != ( last_row - 1 ) ) &&
                      //  ( ( ( propose_deleting_pos_i + 2 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                      //  ) {
                      //  // The position after the one that we are deleting is
                      //  // itself associated with an anchor row, so it has a
                      //  // different anchor row than the one we're deleting.
                      //  // It must be the subsequent one.
                      //  anchor_rows_iter--;
                      //  if( cout_anchor_rows_iter_changes ) {
                      //    anchor_rows_iter_i--;
                      //    cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                      //  }
                      //}
                          
                      // Oh yes, and there is one fewer row...
                      last_row--;
                    
                      // Note what we're doing now.
                      just_shrunk_profile = true;
                      just_grew_profile = false;
                    
                      if( increase_thresholds_for_length_changes ) {
                        if(
                          propose_inserting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {
                          propose_inserting_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_inserting_threshold +=
                            propose_inserting_threshold_increment;
                        }
                        if(
                          propose_inserting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {
                          propose_inserting_prealign_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_inserting_prealign_threshold +=
                            propose_inserting_threshold_increment;
                        }
                        if(
                          propose_inserting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {
                          propose_inserting_postalign_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_inserting_postalign_threshold +=
                            propose_inserting_threshold_increment;
                        }
                        if(
                          propose_deleting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {            
                          propose_deleting_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_deleting_threshold +=
                            propose_deleting_threshold_increment;
                        }
                        if(
                          (
                            ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                            ( m_trainingParameters.verbosity < VERBOSITY_High )
                          ) ||
                          ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                        ) {
                          if(
                            ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
                            ( propose_inserting_threshold == propose_inserting_postalign_threshold )
                          ) {
                            if( propose_inserting_threshold == propose_deleting_threshold ) {
                              cout << "threshold:" << propose_inserting_threshold;
                            } else {
                              cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
                            }
                          } else {
                            cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
                          }
                          if(
                            m_trainingParameters.verbosity >= VERBOSITY_Medium
                          ) {
                            cout << endl;
                          } else {
                            cout.flush();
                          }
                        } // End if verbose
                      } // End if increase_thresholds_for_length_changes
                    
                      // Note that the profile length did change this
                      // iteration.
                      if( !just_caught_preAlign_cycle ) { // Unless it's a preAlign cycle.
                        // Don't reinsert what we just added
                        prev_to_last_length_change_was_insertion =
                          last_length_change_was_insertion;
                        last_length_change_was_deletion = true;
                        last_length_change_was_insertion = false;
                    
                        prev_profile_length_changed_this_iteration =
                          profile_length_changed_this_iteration;
                        profile_length_changed_this_iteration = true;
                      }

                      if(
                         globals_are_at_starting_values ||
                         !enable_reverting
                       ) {
                        if( m_trainingParameters.useUnconditionalBaumWelch && enable_reverting ) {
                            // Don't try to do any unconditional bw updating
                            // this time, since the change in length messes
                            // with the forward calculation, and therefore
                            // makes all of our updates conditional (and on the
                            // wrong values).
                            revert_globals_for_changed_profile_length = true;
                        }
                        // NOTE: If we just_caught_preAlign_cycle, we won't
                        // break now, we will break after retraining the 0th
                        // position, after the next iteration.  Why is this
                        // necessary?  Not quite sure, but if we just break
                        // here, the subsequent Globals phase lowers the score
                        // significantly.
                    
                        // Move along
                        continue; // decrement m_row_i .. 
                      } else if(
                        // TODO: REMOVE!  Testing.
                                //false &&
                        ( length_at_last_revert == m_profile->length() ) &&
                        ( length_at_last_revert == length_at_prev_to_last_revert )
                      ) {
                        // Okay the length doesn't seem to be changing any
                        // more.  Stop reverting.
                        enable_reverting = false;
                        revert_globals_for_changed_profile_length = false;
                          
                        // TODO: REMOVE
                        if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                          cout << endl << "REVERTING DISABLED." << endl;
                        }
                    
                        continue;
                      } else {
                        // There is no point in continuing with the rows.
                        // We're going to have to change the globals to
                        // account for the change in the profile length.
                        // Skip past all the remaining position cycles.
                        revert_globals_for_changed_profile_length = true;

                        m_position_cycle = m_trainingParameters.maxPositionCycles;
                        num_iterations_until_next_length_change_proposal =
                          ( m_trainingParameters.numIterationsBetweenLengthChanges + 1 );
                        break;
                      } // End if ( globals_are_at_starting_values || !enable_reverting ) .. else ..
                    } // End if we should do the deletion
                    // Should we do the insertion?
                    else if(
                    // TODO: REMOVE
                       !fancy_cycle_detected_insertions &&
                       ( !just_caught_preAlign_cycle &&
                         ( !cycle_detected_insertions || ( cycle_detected_deletions && cycle_detected_insertions && ( count_seqs_exceeding_alignment_profile_thresholds ? ( ( last_pos_deletion_fraction_exceeds_threshold ? ( ( double )num_seqs_exceeding_last_pos_deletion_threshold / m_sequence_count ) : ( ( double )num_seqs_exceeding_deletion_threshold / m_sequence_count ) ) < ( ( last_pos_insertion_fraction_exceeds_threshold ? ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / m_sequence_count ) : ( ( double )num_seqs_exceeding_insertion_threshold / m_sequence_count ) ) ) ) : ( ( last_pos_deletion_fraction_exceeds_threshold ? last_pos_deletion_fraction : deletion_fraction ) < ( ( last_pos_insertion_fraction_exceeds_threshold ? last_pos_insertion_fraction : insertion_fraction ) ) ) ) ) )
                       ) &&
                      (
                        last_pos_insertion_fraction_exceeds_threshold ||
                        insertion_fraction_exceeds_threshold
                      )
                       // TODO: PUT BACK.
                      // &&
                      //!(
                      //  profile_length_changed_this_iteration &&
                      //  last_length_change_was_deletion
                      //)
                    ) {
                      // TODO: Update for delaying the last pos insertion..
                      have_inserted_a_pos = true;
                      prev_to_last_inserted_pos_i = last_inserted_pos_i;
                      prev_to_last_inserted_rpos_i = last_inserted_rpos_i;
                      prev_length_before_last_insertion = length_before_last_insertion;
                      last_inserted_pos_i = propose_inserting_pos_i;
                      last_inserted_rpos_i = propose_inserting_rpos_i;
                      length_before_last_insertion = m_profile->length();
                      if( cout_profile_length_changes ) {
                        cout << "\t\tInserting profile position " << propose_inserting_pos_i << endl;
                        if( !cout_indel_fractions ) { // We've already printed it.
                          if( count_seqs_exceeding_alignment_profile_thresholds ) {
                            cout << "\t Num seqs exceeding insertion threshold: " << ( ( propose_inserting_pos_i == last_row ) ? num_seqs_exceeding_last_pos_insertion_threshold : num_seqs_exceeding_insertion_threshold ) << endl;
                            cout << "\t\t fraction is " << ( ( double )( ( propose_inserting_pos_i == last_row ) ? num_seqs_exceeding_last_pos_insertion_threshold : num_seqs_exceeding_insertion_threshold ) / ( double )m_sequence_count ) << endl;
                          } else {
                            cout << "\t Insertion fraction: " << ( ( propose_inserting_pos_i == last_row ) ? last_pos_insertion_fraction : insertion_fraction ) << endl;
                          }
                        }
                        //cout << "[ " << m_row_i << " ] alignment profile position (for all sequences together) is " << alignment_profile_position_allseqs << endl;
                      } // End if cout_profile_length_changes
                      if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                        cout << "i" << propose_inserting_pos_i;
                        // TODO: REMOVE
                        cout << "@";
                        if( last_pos_insertion_fraction_exceeds_threshold ) {
                          if( count_seqs_exceeding_alignment_profile_thresholds ) {
                            cout << "I:" << ( ( double )num_seqs_exceeding_last_pos_insertion_threshold / m_sequence_count );
                          } else {
                            cout << "I:" << last_pos_insertion_fraction;
                          }
                        }
                        if( insertion_fraction_exceeds_threshold ) {
                          if( count_seqs_exceeding_alignment_profile_thresholds ) {
                            cout << "i:" << ( ( double )num_seqs_exceeding_insertion_threshold / m_sequence_count );
                          } else {
                            cout << "i:" << insertion_fraction;
                          }
                        }
                        if(
                          m_trainingParameters.verbosity >= VERBOSITY_Medium
                        ) {
                          cout << endl;
                        } else {
                          cout.flush();
                        }
                      } // End if verbose
                    
                      // First lengthen the profile.
                      //cout << "\t\tOld profile: " << *m_profile << endl;
                    
                      // TODO: Note that if profile were an internalNode,
                      // we'd have to do a lot of upkeep, not just insert a
                      // new position.
                      m_profile->insert(
                        m_profile->begin() + propose_inserting_pos_i,
                        ( *m_profile )[ 0 ]
                      );
                      ( *m_profile )[ propose_inserting_pos_i ].reinitialize(
                        m_profile
                      );
                      // TODO: REMOVE
                      //cout << "Just after insertion, the new position (at " << propose_inserting_pos_i << " is " << ( *m_profile )[ propose_inserting_pos_i ] << endl;
                      // TODO: REMOVE?
                      //alignment_profile_position_allseqs.normalize();
                    
                      if( propose_inserting_pos_i == 0 ) {
                        ( *m_profile )[ propose_inserting_pos_i ][ Emission::Match ] =
                          alignment_profile_position_allseqs[ Emission::Insertion ]; // This is just the initial insertion column.  TODO: Use the extensions too?
                      } else if( propose_inserting_pos_i == last_row ) {
                        ( *m_profile )[ propose_inserting_pos_i ][ Emission::Match ] =
                          // TODO: Something different if there's flanking insertions?
                          // Note that we use the backed-up allseqs...
                          lastpos_alignment_profile_position_allseqs[ Emission::Insertion ]; // This is just the initial insertion column.  TODO: Use the extensions too?
                      } else {
                        // An internal row
                        ( *m_profile )[ propose_inserting_pos_i ][ Emission::Match ] =
                          alignment_profile_position_allseqs[ Emission::Insertion ];
                      }
                      // TODO: REMOVE
                      if( false && cout_profile_length_changes ) {
                        cout << "\t\tNew profile position, before priors and normalization: " << ( *m_profile )[ propose_inserting_pos_i ] << endl;
                        cout << "\t\t lastpos_alignment_profile_position_allseqs is " << lastpos_alignment_profile_position_allseqs << endl;
                        cout << "\t\t alignment_profile_position_allseqs is " << alignment_profile_position_allseqs << endl;
                        cout << "\t\t m_row_i is " << m_row_i << endl;
                      }
                      if( !m_trainingParameters.usePriors ) {
                        // Magic #!  This just weakens the signal a bit.
                        ( *m_profile )[ propose_inserting_pos_i ][ Emission::Match ] += 1.0;
                      } // End if !usePriors
                      if( m_trainingParameters.usePriors ) {
                        m_matchEmissionPrior.incorporatePrior( ( *m_profile )[ propose_inserting_pos_i ] );
                      }
                      ( *m_profile )[ propose_inserting_pos_i ].normalize( m_trainingParameters.profileValueMinimum );
                    
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                      // Also lengthen the scaled match distributions.
                      // Technically the match distributions after this
                      // position don't need to be recalculated, so maybe TODO:
                      // optimize more.
                      m_scaled_match_distributions.resize( m_profile->length() - 1 );
                      for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                        m_profile->createScaledMatchDistributionForPosition(
                          tmp_pos_i,
                          m_scaled_match_distributions[ tmp_pos_i ]
                        );
                      } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          
                      // TODO: REMOVE
                      //cout << "\t\tNew profile: " << *m_profile << endl;
                      if( cout_profile_length_changes ) {
                        cout << "\t\tNew profile length: " << m_profile->length() << endl;
                        // TODO: REMOVE
                        //cout << "\t\tNew profile position: " << ( *m_profile )[ propose_inserting_pos_i ] << endl;
                        //cout << "\t\tlast_row is " << last_row << endl;
                      }
                    
                      if( m_trainingParameters.useUnconditionalBaumWelch ) {
                        m_unconditional_position_ententes_vector.insert(
                          m_unconditional_position_ententes_vector.begin() + propose_inserting_pos_i,
                          m_unconditional_position_ententes_vector[ 0 ]
                        );
                        // We shouldn't have to reinitialize the entente
                        // there, because it'll be copied in when we come
                        // back to this row..
#ifdef ALLOW_BOLTZMANN_GIBBS
                        if( m_trainingParameters.baldiLearningRate > 0 ) {
                          m_unconditional_baldi_position_boltzmann_gibbs_changes_vector.insert(
                            m_unconditional_baldi_position_boltzmann_gibbs_changes_vector.begin() + propose_inserting_pos_i,
                            m_unconditional_baldi_position_boltzmann_gibbs_changes_vector[ 0 ]
                          );
                          // We shouldn't have to reinitialize the object
                          // there, because it'll be copied in when we come
                          // back to this row..
                        } // End if baldiLearningRate > 0
#endif // ALLOW_BOLTZMANN_GIBBS
                     } // End if useUnconditionalBaumWelch
                   
                     // Grow the anchor columns and rows if we're exceeding
                     // the largest we've had in the past.
                     if(
                       ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                       ( m_anchor_rows.m_storeEveryNthRow > 1 )
                     ) {
                       if( propose_inserting_pos_i == last_row ) {
                         // Insert a copy of the last anchor columns
                         // RowVector before it.  Point to the second one.
                         tmp_matrices_iter = m_anchor_columns.end();
                         tmp_matrices_iter--;
                   
                         anchor_columns_iter = // points to the inserted row
                           m_anchor_columns.insert(
                             tmp_matrices_iter, // insert just before end
                             *tmp_matrices_iter // a copy of the last RowVector
                           );
                         // we copied behind the element, so actually we're
                         // now pointing to the second-to-last RowVector, not
                         // to the last one.
                         anchor_columns_iter++;

                         // TODO: REMOVE
                         tmp_matrices_iter = m_anchor_columns.end();
                         tmp_matrices_iter--;
                         assert( anchor_columns_iter == tmp_matrices_iter );

                         // Note that anchor_columns_iter_i does change here
                         // because we don't process the post-aligns until the
                         // second-to-last row (so we're changing it to be the
                         // last row).
                         if( cout_anchor_columns_iter_changes ) {
                           anchor_columns_iter_i = m_anchor_columns.size() - 1;
                           cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                         }
                       } else { // if propose_inserting_pos_i == last_row .. else ..
                         // Insert a copy of the current anchor columns
                         // RowVector before it.  Point to the second one.
                   
                         // TODO: REMOVE
                         tmp_matrices_iter = anchor_columns_iter;
                   
                         anchor_columns_iter = // points to the inserted row
                           m_anchor_columns.insert(
                             tmp_matrices_iter, // insert before here
                             *tmp_matrices_iter // a copy of the current RowVector
                           );
                         // we copied behind the element, so actually we're
                         // now pointing to the copy that is associated with
                         // the current one.
                         anchor_columns_iter++;

                         if( cout_anchor_columns_iter_changes ) {
                           anchor_columns_iter_i++;
                           cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                         }
                       } // End if propose_inserting_pos_i == last_row .. else ..

                       // At this point, the anchor_rows_iter points to the
                       // RowVector we've just inserted to m_anchor_columns.
                       // TODO: REMOVE
#ifndef NDEBUG
                       tmp_matrices_iter = m_anchor_columns.begin();
                       advance( tmp_matrices_iter, ( propose_inserting_pos_i + 1 ) );
                       assert( anchor_columns_iter == tmp_matrices_iter );
#endif // !NDEBUG
                     } // End if m_anchor_columns.m_storeEveryNthColumn > 0
                   
                     // Grow the anchor rows if necessary; either way, make
                     // sure we've got the anchor_rows_iter pointing to the
                     // anchor row associated with the position just before
                     // the inserted position.
                     if(
                       // We need to add a new anchor row if the new length
                       // changes the number of anchor rows: that is, if the
                       // last index, m_profile->length(), divides evenly by
                       // m_anchor_rows.m_storeEveryNthRow.
                       ( ( m_profile->length() % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                     ) {
                       if( propose_inserting_pos_i == last_row ) {
                         // Insert a duplicate copy of the last anchor row.
                         tmp_matrices_iter = m_anchor_rows.end();
                         tmp_matrices_iter--;
                   
                         anchor_rows_iter = // iter will point to new element.
                           m_anchor_rows.insert(
                             tmp_matrices_iter, // inserts just before the end
                             *tmp_matrices_iter // a copy of the last element
                           );
                         // we copied behind the element, so actually we're
                         // now pointing to the second-to-last RowVector, not
                         // to the last one.
                         anchor_rows_iter++;

                         // TODO: REMOVE
                         tmp_matrices_iter = m_anchor_rows.end();
                         tmp_matrices_iter--;
                         assert( anchor_rows_iter == tmp_matrices_iter );
                   
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i = m_anchor_rows.size() - 1;
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                       } else { // if propose_inserting_pos_i == last_row .. else ..
                         // Insert a duplicate copy of the current anchor row.
                   
                         // TODO: REMOVE
                         tmp_matrices_iter = anchor_rows_iter;
                   
                         anchor_rows_iter = // iter will point to new element.
                           m_anchor_rows.insert(
                             anchor_rows_iter, // inserts before this point.
                             *anchor_rows_iter
                           );
                         // we copied behind the element, so actually we're
                         // now pointing to the copy that is associated with
                         // the current one.
                         anchor_rows_iter++;
                   
                         // TODO: REMOVE
                         // REMOVED because technically, tmp_matrices_iter is presently invalid.
                         //assert( anchor_rows_iter == tmp_matrices_iter );
                   
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i++;
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                       } // End if propose_inserting_pos_i == last_row .. else ..

                       // At this point, the anchor_rows_iter points to the
                       // RowVector associated with the new position.
                       // This might not be the right anchor row for the new
                       // position.
                       if( ( ( propose_inserting_pos_i + 1 ) % m_anchor_rows.m_storeEveryNthRow ) != 0 ) {
                         // Then the new row is not an anchor row, meaning we
                         // should point back one row earlier.
                         assert( anchor_rows_iter != m_anchor_rows.begin() );
                         anchor_rows_iter--;
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i--;
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                       }
                   
                       // All anchor rows after anchor_rows_iter are now
                       // invalid, but the only ones that matter are those
                       // for propose_inserting_pos_i and for (
                       // propose_inserting_pos_i + 1 ) .. see below.
                     } else { // if we should add another anchor row .. else ..
                   
                       // Ok we don't need to add a new row, but do we need
                       // to move the iter?  We want it pointing to the anchor row
                       // corresponding to the position we're adding.
                       if( propose_inserting_pos_i == last_row ) {
                         // Since we didn't add a new anchor row, the right one must
                         // must be the existing last anchor row.
                         anchor_rows_iter = m_anchor_rows.end();
                         anchor_rows_iter--;
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i = m_anchor_rows.size() - 1;
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                       } else if( ( ( propose_inserting_pos_i + 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                         // At this point, the anchor_rows_iter points
                         // to the RowVector for the current position,
                         // not for the position we're adding.  Since
                         // the new position does correspond to an
                         // anchor row, we must be currently pointing
                         // to the one before it.  Thus we advance, to
                         // point to the anchor row for the new
                         // position.
                         assert( anchor_rows_iter != m_anchor_rows.end() );
                         anchor_rows_iter++;
                         if( cout_anchor_rows_iter_changes ) {
                           anchor_rows_iter_i++;
                           cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                         }
                       }
                     } // End if we should add another anchor row .. else ..

                     // Ok, so now the anchor_rows_iter points to the
                     // anchor row associated with the newly created
                     // position.
#ifndef NDEBUG
                     tmp_matrices_iter = m_anchor_rows.begin();
                     advance( tmp_matrices_iter, ( ( propose_inserting_pos_i + 1 ) / m_anchor_rows.m_storeEveryNthRow ) );
                     assert( anchor_rows_iter == tmp_matrices_iter );
#endif //!NDEBUG

  #if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        if( !do_not_recalculate_all_forward_rows_during_lengthadjust ) {
                          // Since we've changed the distance of the rows to the
                          // end, we need to recalculate all of the forward rows up
                          // to and including the one corresponding to the newly
                          // added position.
                          tmp_matrices_iter = m_anchor_rows.begin();
                          if(
                            ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                            ( m_anchor_rows.m_storeEveryNthRow > 1 )
                          ) {
                            tmp_matrices_iter2 = m_anchor_columns.begin();
                          }
                      
                          swap = false;
                          for( uint32_t row_i = 0; row_i <= ( propose_inserting_pos_i + 1 ); row_i++ ) {
                            if(
                              ( row_i != 0 ) &&
                              ( ( row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                            ) {
                              tmp_matrices_iter++;
                              assert( tmp_matrices_iter != m_anchor_rows.end() );
                            } // End if storing anchor rows too and this is an anchor row
                            if(
                              ( row_i != 0 ) &&
                              ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                              ( m_anchor_rows.m_storeEveryNthRow > 1 )
                            ) {
                              tmp_matrices_iter2++;
                            }
                      
                            // Every other time, we must swap which forward row is current.
                            swap = !swap;
                      
                            for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                              m_dynamic_programming.forward_calculateRow(
                                m_trainingParameters,
                                false, // don't use viterbi
                                *m_profile,
    //#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                                (
                                  ( ( row_i <= 1 ) || ( ( row_i - 2 ) >= ( last_row - 1 ) ) ) ?
                                  m_scaled_match_distributions[ row_i ] : // ignored
                                  m_scaled_match_distributions[ row_i - 2 ]
                                ),
                                (
                                  ( ( row_i == 0 ) || ( ( row_i - 1 ) >= ( last_row - 1 ) ) ) ?
                                  m_scaled_match_distributions[ row_i ] : // ignored
                                  m_scaled_match_distributions[ row_i - 1 ]
                                ),
    //#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                                m_sequences[ m_seq_i ],
                                row_i,
                                ( swap ? ( *m_forward_rows_ptr )[ m_seq_i ] : ( *m_prev_forward_rows_ptr )[ m_seq_i ] ),
                                ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] )
                              );
                              if(
                                ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                                ( m_anchor_rows.m_storeEveryNthRow > 1 )
                              ) {
                                // This won't copy the rabiner scalars -- see below for that.
                                ( *tmp_matrices_iter2 )[ m_seq_i ].storeColumns(
                                  ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ),
                                  m_anchor_columns.m_storeEveryNthColumn
                                );
                              } // End if storing anchor columns too.
                              if(
                                ( m_anchor_rows.m_storeEveryNthRow != 0 ) &&
                                ( ( row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                              ) {
                                // This won't copy rabiner scalars.  See below for that.
                                ( *tmp_matrices_iter )[ m_seq_i ] =
                                  ( swap ? ( *m_prev_forward_rows_ptr ) : ( *m_forward_rows_ptr ) )[ m_seq_i ];
                              } // End if storing anchor rows too and this is an anchor row
                      
                              if( m_trainingParameters.useRabinerScaling ) {
                                if( row_i == 0 ) {
                                  ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerInverseScalar;
                                } else {
                                  ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar =
                                    (
                                      ( swap ? ( *m_forward_rows_ptr )[ m_seq_i ] : ( *m_prev_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar *
                                      ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerInverseScalar
                                    );
                                }
                                if(
                                  ( m_anchor_columns.m_storeEveryNthColumn != 0 ) &&
                                  ( m_anchor_rows.m_storeEveryNthRow > 1 )
                                ) {
                                  ( *tmp_matrices_iter2 )[ m_seq_i ].m_rabinerInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerInverseScalar;
                                  ( *tmp_matrices_iter2 )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr )[ m_seq_i ] : ( *m_forward_rows_ptr )[ m_seq_i ] ).m_rabinerCumulativeInverseScalar;
                                } // End if anchor_columns != 0;
                                if(
                                  ( m_anchor_rows.m_storeEveryNthRow != 0 ) &&
                                  ( ( row_i % m_anchor_rows.m_storeEveryNthRow ) == 0 )
                                ) {
                                  ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr ) : ( *m_forward_rows_ptr ) )[ m_seq_i ].m_rabinerInverseScalar;
                                  ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                    ( swap ? ( *m_prev_forward_rows_ptr ) : ( *m_forward_rows_ptr ) )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                                } // End if storing anchor rows too and this is an anchor row
                              } // End if useRabinerScaling
                            } // End foreach sequence.
                          } // End foreach row up to and including the newly
                            // inserted one, recalculate the forward rows.
                          if( swap ) {
                            // Rotate the forward_rows
                            m_temp_forward_rows_ptr = m_forward_rows_ptr;
                            m_forward_rows_ptr = m_prev_forward_rows_ptr;
                            m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                          }
                          forward_rows_are_already_rotated_and_calculated = true;
                        } else { // if !do_not_recalculate_all_forward_rows_during_lengthadjust .. else ..
#else // !USE_DEL_IN_DEL_OUT
                     if( true ) {
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                        // Now we need to recalculate the new forward
                        // row, since we're going to be using it, and
                        // the profile has changed at the relevant
                        // position.  We might also need to
                        // recalculate the current forward row, which
                        // we will check first (or to be clearer: the
                        // "current" one in scare quotes, meaning the
                        // one before the newly inserted one -- a
                        // necessary distinction because we delay
                        // post-align insertions).
                      
                        // Only some forward row calculations depend on the index
                        // of the row: when row_i is 0, 1, or last_row.  (unless
                        // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT, in which case they all depend on
                        // the index of the row).
                       if( // here "propose_inserting_pos_i" is standing in for the row associated with the position before propose_inserting_pos_i ..
                          ( propose_inserting_pos_i == 0 ) ||
                          ( propose_inserting_pos_i == 1 ) ||
                          ( propose_inserting_pos_i == last_row )
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          || true
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        ) {
                          if( propose_inserting_pos_i == last_row ) {
                            // First, rotate the forward_rows to put the
                            // forward_row in the prev_forward_row place where it
                            // belongs (remember, m_row_i is actually ( last_row
                            // - 1 ) because we delayed the post-align insertion).
                            m_temp_forward_rows_ptr = m_forward_rows_ptr;
                            m_forward_rows_ptr = m_prev_forward_rows_ptr;
                            m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                          }
                          for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                            // Then we also need to recalculate the current
                            // forward row..
                            m_dynamic_programming.forward_calculateRow(
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              m_scaled_match_distributions[ propose_inserting_pos_i - 2 ],
                              m_scaled_match_distributions[ propose_inserting_pos_i - 1 ],
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                               m_sequences[ m_seq_i ],
                               propose_inserting_pos_i,
                              ( *m_prev_forward_rows_ptr )[ m_seq_i ], // ignored if propose_inserting_pos_i == 0
                              ( *m_forward_rows_ptr )[ m_seq_i ]
                             );
                             if( m_trainingParameters.useRabinerScaling ) {
                               if( propose_inserting_pos_i == 0 ) {
                                 ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                               } else {
                                 ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                                 ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar *=
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                               } // End if propose_inserting_pos_i == 0 .. else ..
                             } // End if useRabinerScaling
                             // Keep the anchor columns up-to-date
                             if(
                               ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                               ( m_anchor_rows.m_storeEveryNthRow > 1 )
                             ) {
                               // Remember, the anchor_columns_iter is
                               // currently pointing to the
                               // newly-inserted one.
                               tmp_matrices_iter = anchor_columns_iter;
                               tmp_matrices_iter--;

                               ( *tmp_matrices_iter )[ m_seq_i ].storeColumns(
                                 ( *m_forward_rows_ptr )[ m_seq_i ],
                                 m_anchor_columns.m_storeEveryNthColumn
                               );
                               if( m_trainingParameters.useRabinerScaling ) {
                                 // Copy the rabiner scalars.
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                               } // End if useRabinerScaling
                             }  // End if m_anchor_columns.m_storeEveryNthColumn > 0
                     
                             // Keep the anchor rows up-to-date
                             if( ( propose_inserting_pos_i % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                               // Remember that at this point,
                               // anchor_rows_iter points to the
                               // anchor row for the newly created
                               // position...  What we want to replace
                               // here is the anchor row associated
                               // with the position just before the
                               // newly created one...
                               if( ( ( propose_inserting_pos_i + 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                                 // TODO: REMOVE
                                 assert( anchor_rows_iter != m_anchor_rows.begin() );
                                 tmp_matrices_iter = anchor_rows_iter;
                                 tmp_matrices_iter--;
                               }
                   
                               // TODO: Don't copy the data!  Copy the pointer!
                               ( *tmp_matrices_iter )[ m_seq_i ] =
                                 ( *m_forward_rows_ptr )[ m_seq_i ];
                               if( m_trainingParameters.useRabinerScaling ) {
                                 // Copy the rabiner scalars.
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                                 ( *tmp_matrices_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                   ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                               } // End if useRabinerScaling
                             } // End if this is an anchor row, store it.
                           } // End foreach sequence
                         } // End if we're inserting after the first, second,
                         // or last row, recalc the current forward row.
                       
                         // Now calculate the new forward row.
                       
                         // First, rotate the forward_rows
                         m_temp_forward_rows_ptr = m_forward_rows_ptr;
                         m_forward_rows_ptr = m_prev_forward_rows_ptr;
                         m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                     
                         for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                           m_dynamic_programming.forward_calculateRow(
                             m_trainingParameters,
                             false, // don't use viterbi
                             *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                            (
                              ( propose_inserting_pos_i == 0 ) ?
                              m_scaled_match_distributions[ 0 ] : // ignored
                              m_scaled_match_distributions[ propose_inserting_pos_i - 1 ]
                            ),
                            m_scaled_match_distributions[ propose_inserting_pos_i ],
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                            m_sequences[ m_seq_i ],
                            propose_inserting_pos_i + 1,
                            ( *m_prev_forward_rows_ptr )[ m_seq_i ],
                            ( *m_forward_rows_ptr )[ m_seq_i ]
                          );
                          if( m_trainingParameters.useRabinerScaling ) {
                            ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                              ( *m_prev_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                            ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar *=
                              ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                          } // End if useRabinerScaling
                          // Keep the anchor columns up-to-date
                          if(
                            ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                            ( m_anchor_rows.m_storeEveryNthRow > 1 )
                          ) {
                            // We've already set the iterator to the right place (see above).
                            ( *anchor_columns_iter )[ m_seq_i ].storeColumns(
                              ( *m_forward_rows_ptr )[ m_seq_i ],
                              m_anchor_columns.m_storeEveryNthColumn
                            );
                            if( m_trainingParameters.useRabinerScaling ) {
                              // Copy the rabiner scalars.
                              ( *anchor_columns_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                              ( *anchor_columns_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                            } // End if useRabinerScaling
                          } // End if m_anchor_columns.m_storeEveryNthColumn > 0
                    
                          // Keep the anchor rows up-to-date
                          if( ( ( propose_inserting_pos_i + 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 ) {
                            // We've already set the iterator to the right place (see above).
                    
                            // TODO: Don't copy the data!  Copy the pointer!
                            ( *anchor_rows_iter )[ m_seq_i ] =
                              ( *m_forward_rows_ptr )[ m_seq_i ];
                            if( m_trainingParameters.useRabinerScaling ) {
                              // Copy the rabiner scalars.
                              ( *anchor_rows_iter )[ m_seq_i ].m_rabinerInverseScalar =
                                ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                              ( *anchor_rows_iter )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar;
                            } // End if useRabinerScaling
                          } // End if this is an anchor row, store it.
                        } // End foreach m_seq_i
                      
                        // And remember that we already rotated the rows.
                        forward_rows_are_already_rotated_and_calculated = true;
                      } // End if true (or maybe if do_not_recalculate_all_forward_rows_during_lengthadjust)
                    
                      // At this point the anchor rows and cols iters point to
                      // the row for the newly-inserted position.
                    
                      if( propose_inserting_pos_i != last_row ) {
                        // Rotate the backward_rows, since we're effectively
                        // undoing our last move and doing it again with the
                        // new row.  If we just inserted the new final row,
                        // we're already off (but it doesn't matter since in
                        // calculating backward rows for the final row, the
                        // next_backward_rows_ptr is ignored.
                        m_temp_backward_rows_ptr = m_backward_rows_ptr;
                        m_backward_rows_ptr = m_next_backward_rows_ptr;
                        m_next_backward_rows_ptr = m_temp_backward_rows_ptr;
                        // On the next time 'round (when we decrement m_row_i),
                        // it will be rotated again, making the
                        // m_backward_rows_ptr become the
                        // m_next_backward_rows_ptr for the next time thru.
                      
                        if( m_trainingParameters.useRabinerScaling ) {
                          // The backward rows were scaled by the
                          // rabinerInverseScalar of the corresponding forward
                          // matrix rows, so now the scale factor needs to be
                          // updated.
                          ScoreType scale_ratio;
                          for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                            // This is the old scale:
                            scale_ratio =
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            scale_ratio /=
                              // The new scale:
                              //m_forward_matrices[ propose_inserting_pos_i + 1 ][ m_seq_i ].m_rabinerInverseScalar;
                              ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            // Rescale it to the new scale.
                            ( *m_backward_rows_ptr )[ m_seq_i ] *= scale_ratio;
                      
                            // Note the new scale.
                            ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                              //////m_forward_matrices[ propose_inserting_pos_i + 1 ][ m_seq_i ].m_rabinerInverseScalar;
                            ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            if( ( propose_inserting_pos_i + 1 ) == ( last_row + 1 ) ) {
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            } else {
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar /=
                                scale_ratio;
                            } // End if ( propose_inserting_pos_i + 1 ) == ( last_row + 1 ) .. else ..
                          } // End foreach sequence, fix the backward row's rabinerInverseScalars.. 
                        } // End if useRabinerScaling
                      } // End if propose_inserting_pos_i != last_row, rotate
                        // backward rows.
                    
                      // If we don't alwaysAccept, then we need to
                      // recalculate the score, so we know if we have an
                      // improvement next time (Note that we can't reject the
                      // length change; this is for the next time around,
                      // when we come back to the newly created row...)
                      // TODO: PUT BACK.. TESTING.
                      //if( !( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) ) {
                      if( true ) {
                    
                        // To efficiently recalculate the score, we can
                        // calculate the backward_rows for the newly created
                        // row now...
                        // Rotate the backward_rows
                        m_temp_backward_rows_ptr = m_backward_rows_ptr;
                        m_backward_rows_ptr = m_next_backward_rows_ptr;
                        m_next_backward_rows_ptr = m_temp_backward_rows_ptr;
                    
                        backward_rows_are_already_rotated_and_calculated = true;
                    
                        for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                    
                          if( m_trainingParameters.useRabinerScaling ) {
                            ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar =
                              /////m_forward_matrices[ ( propose_inserting_pos_i + 1 ) ][ m_seq_i ].m_rabinerInverseScalar;
                              ( *m_forward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                          } // End if useRabinerScaling
                          m_dynamic_programming.backward_calculateRow(
                            m_trainingParameters,
                            false, // don't use viterbi
                            *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                            m_scaled_match_distributions[ propose_inserting_pos_i ],
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                            m_sequences[ m_seq_i ],
                            ( propose_inserting_pos_i + 1 ),
                            ( *m_next_backward_rows_ptr )[ m_seq_i ],
                            ( *m_backward_rows_ptr )[ m_seq_i ]
                          );
                          if( m_trainingParameters.useRabinerScaling ) {
                            if( ( propose_inserting_pos_i + 1 ) == ( last_row + 1 ) ) {
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar;
                            } else {
                              ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar =
                                (
                                  ( *m_next_backward_rows_ptr )[ m_seq_i ].m_rabinerCumulativeInverseScalar *
                                  ( *m_backward_rows_ptr )[ m_seq_i ].m_rabinerInverseScalar
                                );
                            } // End if ( propose_inserting_pos_i + 1 ) == ( last_row + 1 ) .. else ..
                          } // End if useRabinerScaling
                        } // End foreach sequence, calculate the backward
                        // rows for the new row.
                    
                        // TODO: REMOVE
                        //cout << "(2) Ending score, before recalculating, is " << m_endingScore << endl;
                        m_endingScore =
                          m_dynamic_programming.forward_score( 
                            m_trainingParameters,
                            *m_profile,
                            m_sequence_count,
                            ( propose_inserting_pos_i + 1 ),
                            //m_forward_matrices[ propose_inserting_pos_i + 1 ],
                            ( *m_forward_rows_ptr ),
                            ( *m_backward_rows_ptr ),
                            &sequence_scores,
                            &largest_sequence_score
                          );
                        // TODO: REMOVE!
                        //cout << "sequence 1's forward row for row " << ( propose_inserting_pos_i + 1 ) << ", used to recalculate the score, is " << ( *m_forward_rows_ptr )[ 1 ] << endl;
                        //cout << "sequence 1's backward row for row " << ( propose_inserting_pos_i + 1 ) << ", used to recalculate the score, is " << ( *m_backward_rows_ptr )[ 1 ] << endl;
                        // TODO: REMOVE!
                        //cout << "sequence 1's score is " << sequence_scores[ 1 ] << endl;
                        // TODO: REMOVE
                        //cout << "(2) Ending score, calculated the new way, is " << m_endingScore << endl;
                        // TODO: REMOVE
                    
                        // TODO: REMOVE?
                        assert( !isnan( m_endingScore ) );
                    
                      } // End if !alwaysAccept, recalculate the score(s)
                    
                      // Keep track of the creation time.
                      lengthadjust_position_creation_iters.insert(
                        lengthadjust_position_creation_iters.begin() + propose_inserting_pos_i,
                        ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration )
                      );
                    
                      // Now increment the m_row_i and last_row
                      last_row++;
                    
                      m_row_i = ( propose_inserting_pos_i + 2 ); // it will be
                                                                 // decremented
                                                                 // on
                                                                 // continue,
                                                                 // so
                      // it will be propose_inserting_pos_i + 1 when we come
                      // around again.
                    
                      // This too will be decremented when we come around again.
                      if(
                        ( ( ( ( m_row_i - 1 ) + 1 ) % m_anchor_rows.m_storeEveryNthRow ) == 0 ) &&
                        ( ( m_row_i - 1 ) != last_row ) // Note this is the *new* m_row_i & last_row!
                      ) {
                        // TODO: REMOVE
                        assert( anchor_rows_iter != m_anchor_rows.end() );
                        anchor_rows_iter++;
                        if( cout_anchor_rows_iter_changes ) {
                          anchor_rows_iter_i++;
                          cout << "[anchor_rows_iter_i=" << anchor_rows_iter_i << "]";
                        }
                        // TODO: REMOVE
                        //assert( anchor_rows_iter == ( m_anchor_rows.begin() + ( propose_inserting_pos_i + 2 ) ) );
                      } // End if we need to increment the anchor rows iterator
                    
                      // This too will be decremented when we come around again.
                      if(
                        ( ( m_row_i - 1 ) != last_row ) &&
                        ( m_anchor_columns.m_storeEveryNthColumn > 0 ) &&
                        ( m_anchor_rows.m_storeEveryNthRow > 1 )
                      ) {
                        // TODO: REMOVE
                        assert( anchor_columns_iter != m_anchor_columns.end() );
                        anchor_columns_iter++;
                        if( cout_anchor_columns_iter_changes ) {
                          anchor_columns_iter_i++;
                          cout << "[anchor_columns_iter_i=" << anchor_columns_iter_i << "]";
                        }
                      } // End if ( m_row_i - 1 ) != last_row and m_anchor_columns.m_storeEveryNthColumn > 0
                    
                      // Note what we're doing now.
                      just_grew_profile = true;
                      just_shrunk_profile = false;
                    
                      if( increase_thresholds_for_length_changes ) {
                        if(
                          propose_inserting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {
                          propose_inserting_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_inserting_threshold +=
                            propose_inserting_threshold_increment;
                        }
                        if(
                          propose_inserting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {
                          propose_inserting_prealign_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_inserting_prealign_threshold +=
                            propose_inserting_threshold_increment;
                        }
                        if(
                          propose_inserting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {
                          propose_inserting_postalign_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_inserting_postalign_threshold +=
                            propose_inserting_threshold_increment;
                        }
                        if(
                          propose_deleting_threshold_increment <
                          min_threshold_increment_for_increasing_threshold_for_length_changes
                        ) {            
                          propose_deleting_threshold +=
                            min_threshold_increment_for_increasing_threshold_for_length_changes;
                        } else {
                          propose_deleting_threshold +=
                            propose_deleting_threshold_increment;
                        }
                        if(
                          (
                            ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                            ( m_trainingParameters.verbosity < VERBOSITY_High )
                          ) ||
                          ( m_trainingParameters.verbosity >= VERBOSITY_All ) 
                        ) {
                          if(
                            ( propose_inserting_threshold == propose_inserting_prealign_threshold ) &&
                            ( propose_inserting_threshold == propose_inserting_postalign_threshold )
                          ) {
                            if( propose_inserting_threshold == propose_deleting_threshold ) {
                              cout << "threshold:" << propose_inserting_threshold;
                            } else {
                              cout << "thresholds:(ins=" << propose_inserting_threshold << ",del=" << propose_deleting_threshold << ")";
                            }
                          } else {
                            cout << "thresholds:(ins=" << propose_inserting_threshold << ",pre=" << propose_inserting_prealign_threshold << ",post=" << propose_inserting_postalign_threshold << ",del=" << propose_deleting_threshold << ")";
                          }
                          if(
                            m_trainingParameters.verbosity >= VERBOSITY_Medium
                          ) {
                            cout << endl;
                          } else {
                            cout.flush();
                          }
                        } // End if verbose
                      } // End if increase_thresholds_for_length_changes
  
                      // Don't redelete what we just added
                      prev_to_last_length_change_was_insertion =
                        last_length_change_was_insertion;
                      last_length_change_was_insertion = true;
                      last_length_change_was_deletion = false;
                      // Note that the profile length did change this
                      // iteration.
                      prev_profile_length_changed_this_iteration =
                        profile_length_changed_this_iteration;
                      profile_length_changed_this_iteration = true;
                    
                      if( globals_are_at_starting_values || !enable_reverting ) {
                        if(
                          m_trainingParameters.useUnconditionalBaumWelch
                        ) {
                          if( enable_reverting ) {
                            // Don't try to do any unconditional bw updating
                            // this time, since the change in length messes
                            // with the forward calculation, and therefore
                            // makes all of our updates conditional (and on the
                            // wrong values).
                            revert_globals_for_changed_profile_length = true;
                          }
                        }
                        // move along
                        continue;
                      } else if(
                        // TODO: REMOVE!  Testing.
                                //false &&
                        ( length_at_last_revert == m_profile->length() ) &&
                        ( length_at_last_revert == length_at_prev_to_last_revert )
                      ) {
                        // Okay the length doesn't seem to be changing any
                        // more.  Stop reverting.
                        enable_reverting = false;
                        revert_globals_for_changed_profile_length = false;
                          
                        // TODO: REMOVE
                        if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                          cout << endl << "REVERTING DISABLED." << endl;
                        }
  
                        continue;
                      } else {
                        // There is no point in continuing with the rows.
                        // We're going to have to change the globals to
                        // account for the change in the profile length.
                        // Skip past all the remaining position cycles.
                        revert_globals_for_changed_profile_length = true;
                        m_position_cycle =
                          m_trainingParameters.maxPositionCycles;
                        num_iterations_until_next_length_change_proposal =
                          ( m_trainingParameters.numIterationsBetweenLengthChanges + 1 );
                        break;
                      } // End if ( globals_are_at_starting_values || !enable_reverting ) .. else ..
                    } // End if there's a cycle .. else if we should delete .. else if we should insert ..
                  } // End if m_row_i < last_row
                } // End if we are proposing profile length changes, do it.

                if(
                  ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
                  ( m_row_i == 0 )
                ) {
                  if( num_iterations_until_next_length_change_proposal == 0 ) {
                    num_iterations_until_next_length_change_proposal =
                      ( m_trainingParameters.numIterationsBetweenLengthChanges + 1 );
                  }

                  // And now break, since we can't train in row 0 ( no more
                  // profile positions! )
                  break;
                }

                // Incorporate priors OR do Baldi / Siegel calcs
                if( m_trainingPhase == TRAINING_PHASE_Positions ) {
#ifdef ALLOW_BOLTZMANN_GIBBS
                  // Don't use priors if using Baldi..
                  if( m_trainingParameters.baldiLearningRate > 0 ) {
                    m_baldi_position_boltzmann_gibbs_change.zero();
                    // NOTE: This would be messed up by usePriors, since it is
                    // based on computing the contribution of the
                    // (un-altered-by-priors) sequence-specific
                    // position_entente_update to change in the total
                    // position_entente.
                    for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                      m_dynamic_programming.updatePositionBoltzmannGibbsChangeForSequence(
                        m_trainingParameters,
                        ( *m_profile )[ m_row_i - 1 ],
                        m_coefficients_vector[ m_seq_i ],
                        NULL,
                        m_baldi_position_boltzmann_gibbs_change
                      );
                      if( false && m_baldi_be_verbose ) {
                        cout << "[seq " << m_seq_i << "] m_baldi_position_boltzmann_gibbs_change is now " << m_baldi_position_boltzmann_gibbs_change << endl;
                        cout << "\tprof pos is " << ( *m_profile )[ m_row_i - 1 ] << endl;
                        cout << "\tcoeffs vector is " << m_coefficients_vector[ m_seq_i ] << endl;
                        cout << "\tpos entente is " << m_position_entente << endl;
                      }
                    } // End foreach sequence
                  } else
#endif // ALLOW_BOLTZMANN_GIBBS
                  if( m_trainingParameters.usePriors ) {
                    // TODO: REMOVE
                    assert( m_position_entente.m_scalar != 0 );
                    // TODO: REMOVE
                    //cout << "BEFORE usePriors, m_position_entente is " << m_position_entente << endl;
                    //m_position_entente.unscale(); // TODO: REMOVE.  TESTING.
                    m_matchEmissionPrior.incorporatePrior( m_position_entente );
                    // TODO: REMOVE
                    //cout << "AFTER usePriors, m_position_entente is " << m_position_entente << endl;
                  } // End if ( baldiLearningRate > 0 ) .. else if usePriors ..
                } // End if TRAINING_PHASE_Positions
          
                if( m_trainingParameters.debug > DEBUG_None ) {
                  cout << ".done." << endl;
                } // End if debug
          

              if( ( m_trainingPhase != TRAINING_PHASE_Positions ) &&
                  ( m_row_i == 0 ) ) {

                  //TODO: REMOVE
                  //if( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) == 5 ) {
                  //  cout << "[ " << m_row_i << " ]proposed new profile globals (before normalizing): " << m_global_entente << endl;
                  //  cout << "[ " << m_row_i << " ]\t unscaled, proposed new profile globals (before normalizing): " << m_global_entente.createUnscaledCopy() << endl;
                  //}
                  ////TODO: REMOVE
                  //if( false ) {
                  //  cout << "alignment_profile_position_allseqs is " << alignment_profile_position_allseqs << endl;
                  //  cout << "\t unscaled, alignment_profile_position_allseqs is " << alignment_profile_position_allseqs.createUnscaledCopy() << endl;
                  //} // End if testing
          
                  if( m_trainingParameters.usePriors ) {
                    // TODO: REMOVE
                    //m_global_entente.unscale();
                    //cout << "BEFORE usePriors, m_global_entente is " << m_global_entente << endl;
                    //cout << "m_globalPrior is " << m_globalPrior << endl;

                    m_globalPrior.incorporatePrior( m_global_entente );

                    // TODO: REMOVE
                    //cout << "AFTER usePriors, m_global_entente is " << m_global_entente << endl;
                  } // End if usePriors

                  // We are currently training the profile global params
                  m_global_entente.normalize(
                    m_trainingParameters.profileValueMinimum
                  );

                  // TODO: REMOVE! TESTING theory about C->C being the problem.
                  //m_global_entente[
                  //  ProfileType::Transition::fromPostAlign
                  //][
                  //  ProfileType::TransitionFromPostAlign::toPostAlign
                  //] = .5;
                  //m_global_entente[
                  //  ProfileType::Transition::fromPostAlign
                  //][
                  //  ProfileType::TransitionFromPostAlign::toTerminal
                  //] = .5;

                  //TODO: REMOVE
                  //if( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) == 5 ) {
                  //  cout << "[ " << m_row_i << " ] proposed new profile globals (after normalizing): " << m_global_entente << endl;
                  //}

                  if( m_trainingParameters.useUnconditionalBaumWelch ) {

                    if(
                      !(
                        m_trainingParameters.proposeProfileLengthChanges &&
                        revert_globals_for_changed_profile_length
                      )
                    ) {
                      m_unconditional_backupProfile.copyFrom( *m_profile );
                    } // Back up if there's any chance we'll wish we did.

#ifdef ALLOW_BOLTZMANN_GIBBS
                    if( m_trainingParameters.baldiLearningRate > 0 ) { // Do Baldi or Baldi / Siegel ..

                      // TODO: REMOVE
                      //  if(
                      //     m_trainingParameters.proposeProfileLengthChanges &&
                      //     revert_globals_for_changed_profile_length
                      //     ) {
                      //    cout << "SKIPPING BALDI BECAUSE WE'RE REVERTING GLOBALS." << endl;
                      //}
                      
                      if( 
                         !(
                           m_trainingParameters.proposeProfileLengthChanges &&
                           revert_globals_for_changed_profile_length
                           )
                          // TODO: REMOVE.  TESTING.
                          && ( ubw_cbw_hybrid ? ( ( m_iteration % 2 ) == 0 ) : true ) // only train positions every *other* iter
                      ) {
                        m_startingScore_position_step = m_endingScore;
                        
                        if( m_trainingParameters.baldiHybrid ) {

                          // Also make sure we couldn't do better just using the UBW update.
                          for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                            if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                                !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                              // Then don't train it.
                              continue;
                            }
                            m_position_entente =
                              m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                            m_position_entente.normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                            ( *m_profile )[ m_row_i - 1 ] = m_position_entente; // already normalized.
                          } // End foreach row_i
                              
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          // Make sure the scaled distributions are up-to-date.
                          for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                            m_profile->createScaledMatchDistributionForPosition(
                              tmp_pos_i,
                              m_scaled_match_distributions[ tmp_pos_i ]
                            );
                          } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          bool swap =
                            m_dynamic_programming.forward_score( 
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences,
                              m_sequence_count,
                              &m_anchor_columns,
                              &m_anchor_rows,
                              *m_prev_forward_rows_ptr,
                              *m_forward_rows_ptr,
                              &sequence_scores,
                              &largest_sequence_score,
                              m_endingScore
                            );
                          //m_endingScore =
                          //  m_dynamic_programming.forward_score( 
                          //    m_trainingParameters,
                          //    *m_profile,
                          //    m_sequences,
                          //    m_sequence_count,
                          //    m_forward_matrices,
                          //    &sequence_scores,
                          //    &largest_sequence_score
                          //  );
                          // TODO: REMOVE?
                          assert( !isnan( m_endingScore ) );

                          // Return it to what it was before.
                          m_profile->copyPositions( m_unconditional_backupProfile );
  
                          if( m_baldi_be_verbose ) {
                            cout << "The score if we just use standard UBW is " << m_endingScore << endl;
                          }

                        } // End if m_trainingParameters.baldiHybrid

                        if( m_trainingParameters.siegelMaxFindingThePeakAttempts_positions == 0 ) {
                          // Just use straight-up Baldi.
                          for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                            if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                                !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                              // Then don't train it.
                              continue;
                            }

                            m_position_entente =
                              m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                            m_position_entente.normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                        
                            if( false && m_baldi_be_verbose ) {
                              cout << "The profile position was " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                              
                            m_baldi_position_boltzmann_gibbs_change =
                              m_unconditional_baldi_position_boltzmann_gibbs_changes_vector[ m_row_i - 1 ];
                              
                            // Convert the profile position to Boltzmann-Gibbs.
                            m_baldi_position_boltzmann_gibbs.fromProfilePosition(
                              ( *m_profile )[ m_row_i - 1 ]
                            );
                            m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                              m_trainingParameters.baldiLearningRate;

                            // In the Baldi paper, this corresponds to equation
                            // 3-prime (the learning rate should incorporate the
                            // current score):
                            //m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                            //  m_startingScore_position_step; // the current score...
                            // This is wrong, but seems to do better: ! Even still, not so good.
                            //m_baldi_position_boltzmann_gibbs_change.m_scalar *=
                            //  m_startingScore_position_step; // the current score...

                            if( false && m_baldi_be_verbose ) {
                              cout << "weighted by the learning-rate and by current score m_startingScore_position_step, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                            }
                          
                            m_baldi_position_boltzmann_gibbs +=
                              m_baldi_position_boltzmann_gibbs_change;
                          
                            m_baldi_position_boltzmann_gibbs.toProfilePosition(
                              ( *m_profile )[ m_row_i - 1 ]
                            );
                            if( false && m_baldi_be_verbose ) {
                              cout << "... so the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                            ( *m_profile )[ m_row_i - 1 ].normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                            if( false && m_baldi_be_verbose ) {
                              cout << "Normalized, the profile position is now " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                          } // End foreach row_i

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          // Make sure the scaled distributions are up-to-date.
                          for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                            m_profile->createScaledMatchDistributionForPosition(
                              tmp_pos_i,
                              m_scaled_match_distributions[ tmp_pos_i ]
                            );
                          } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          bool swap =
                            m_dynamic_programming.forward_score( 
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences,
                              m_sequence_count,
                              &m_anchor_columns,
                              &m_anchor_rows,
                              *m_prev_forward_rows_ptr,
                              *m_forward_rows_ptr,
                              &sequence_scores,
                              &largest_sequence_score,
                              score_after_modification
                            );
                          //score_after_modification =
                          //  m_dynamic_programming.forward_score( 
                          //    m_trainingParameters,
                          //    *m_profile,
                          //    m_sequences,
                          //    m_sequence_count,
                          //    m_forward_matrices,
                          //    &sequence_scores,
                          //    &largest_sequence_score
                          //  );
                          // TODO: REMOVE?
                          assert( !isnan( score_after_modification ) );
  
                          // Return it to what it was before.
                          m_profile->copyPositions( m_unconditional_backupProfile );

                        } else { // End if use Baldi .. else use Baldi / Siegel ..

                          // Get three sample points for which the
                          // score of the third is less than that of
                          // the second.  (Capture the peak.)
                          sample_indices[ 0 ] = 0;
                          sample_indices[ 1 ] = 1;
                          sample_indices[ 2 ] = 2;
                          sample_scores[ 0 ] = 0;
                          sample_scores[ 1 ] = 0;
                          sample_scores[ 2 ] = 0;
                          sample_epsilons[ 0 ] = 0;
                          sample_epsilons[ 1 ] = 0;
                          sample_epsilons[ 2 ] = 0;
                          epsilon_scale_value_scale_factor = m_trainingParameters.siegelEpsilonScaleFactor;
                          for( finding_the_peak_attempts = 0,
                                 epsilon = ( 1.0 / ( epsilon_scale_value_scale_factor * epsilon_scale_value_scale_factor * epsilon_scale_value_scale_factor ) ),
                                 epsilon_scale_value = epsilon_scale_value_scale_factor;
                                 (
                                  ( epsilon >= m_trainingParameters.siegelMinEpsilon ) &&
                                  ( finding_the_peak_attempts <=
                                    m_trainingParameters.siegelMaxFindingThePeakAttempts_positions ) &&
                                  !(
                                    ( sample_epsilons[ 2 ] == 0 ) &&
                                    ( finding_the_peak_attempts >
                                      m_trainingParameters.siegelMaxFindingTheGradientAttempts_positions )
                                  ) &&
                                  (
                                   ( sample_epsilons[ 2 ] == 0 ) ||
                                    ( sample_scores[ sample_indices[ 2 ] ] >=
                                      sample_scores[ sample_indices[ 1 ] ] )
                                  )
                                 );
                               ++finding_the_peak_attempts,
                                 ( epsilon *= epsilon_scale_value )
                          ) { // for( finding_the_peak_attempts ... )
                            
                            if( m_baldi_be_verbose ) {
                              cout << " Unconditional Baldi / Siegel with finding-the-peak with epsilon = " << epsilon << endl;
                            }

                            for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                              if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                                  !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                                // Then don't train it.
                                continue;
                              }

                              m_position_entente =
                                m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                              m_position_entente.normalize(
                                m_trainingParameters.profileValueMinimum
                              );
                        
                              if( false && m_baldi_be_verbose ) {
                                cout << "The profile position was " << ( *m_profile )[ m_row_i - 1 ] << endl;
                              }
                              
                              m_baldi_position_boltzmann_gibbs_change =
                                m_unconditional_baldi_position_boltzmann_gibbs_changes_vector[ m_row_i - 1 ];
                              
                              // Convert the profile position to Boltzmann-Gibbs.
                              m_baldi_position_boltzmann_gibbs.fromProfilePosition(
                                ( *m_profile )[ m_row_i - 1 ]
                              );
                              if( false && m_baldi_be_verbose ) {
                                cout << ".. so the corresponding m_baldi_position_boltzmann_gibbs is " << m_baldi_position_boltzmann_gibbs << endl;
                                cout << "Scaled, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                              }

                              m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                                epsilon;
                            
                              m_baldi_position_boltzmann_gibbs +=
                                m_baldi_position_boltzmann_gibbs_change;
                            
                              m_baldi_position_boltzmann_gibbs.toProfilePosition(
                                ( *m_profile )[ m_row_i - 1 ]
                              );
                              if( false && m_baldi_be_verbose ) {
                                cout << "... so the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                              }
                              ( *m_profile )[ m_row_i - 1 ].normalize(
                                m_trainingParameters.profileValueMinimum
                              );
                            } // End foreach row_i (in trying to find the peak using epsilon)

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                            // Make sure the scaled distributions are up-to-date.
                            for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                              m_profile->createScaledMatchDistributionForPosition(
                                tmp_pos_i,
                                m_scaled_match_distributions[ tmp_pos_i ]
                              );
                            } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                            bool swap =
                              m_dynamic_programming.forward_score( 
                                m_trainingParameters,
                                false, // don't use viterbi
                                *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                                m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                                m_sequences,
                                m_sequence_count,
                                &m_anchor_columns,
                                &m_anchor_rows,
                                *m_prev_forward_rows_ptr,
                                *m_forward_rows_ptr,
                                &sequence_scores,
                                &largest_sequence_score,
                                score_after_modification
                              );
                            //score_after_modification =
                            //  m_dynamic_programming.forward_score( 
                            //    m_trainingParameters,
                            //    *m_profile,
                            //    m_sequences,
                            //    m_sequence_count,
                            //    m_forward_matrices,
                            //    &sequence_scores,
                            //    &largest_sequence_score
                            //  );
                            // TODO: REMOVE?
                            assert( !isnan( score_after_modification ) );
                            
                            // Return it to what it was before.
                            m_profile->copyPositions( m_unconditional_backupProfile );
                            
                            if( m_baldi_be_verbose ) {
                              cout << "Got score diff " << ( ( m_startingScore_position_step > score_after_modification ) ? "(negative) " : "" ) << ( ( m_startingScore_position_step > score_after_modification ) ? ( m_startingScore_position_step / score_after_modification ) : ( score_after_modification / m_startingScore_position_step ) ) << ": was " << m_startingScore_position_step << "; is now " << score_after_modification << endl;
                            }
                            
                            if( epsilon_scale_value < 1 ) {
                              assert( epsilon_scale_value == ( 1.0 / epsilon_scale_value_scale_factor ) );
                              if( finding_the_peak_attempts == 0 ) {
                                assert( false ); // The code expects that epsilon_scale_value > 1 to start with...
                                cout << "ERROR: finding_the_peak_attempts == 0 but epsilon_scale_value < 1!" << endl;
                                exit( 1 );
                                //assert( sample_indices[ 0 ] == 0 );
                                //assert( sample_indices[ 1 ] == 1 );
                                //assert( sample_indices[ 2 ] == 2 );
                                //
                                //sample_scores[ 0 ] =
                                //  score_after_modification;
                                //sample_epsilons[ 0 ] =
                                //  epsilon;
                              } else if( sample_epsilons[ 1 ] == 0 ) {
                                // Set the scale value.
                                assert( sample_indices[ 0 ] == 0 );
                                assert( sample_indices[ 1 ] == 1 );
                                assert( sample_indices[ 2 ] == 2 );
                        
                                if(
                                   score_after_modification < sample_scores[ 0 ]
                                  ) {
                                  // Score went down.  Flip.
                                  epsilon_scale_value = epsilon_scale_value_scale_factor;
                        
                                  // Put the new score in the front..
                                  sample_scores[ 1 ] = sample_scores[ 0 ];
                                  sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                  sample_scores[ 0 ] = score_after_modification;
                                  sample_epsilons[ 0 ] = epsilon;
                        
                                  // We've already calculated epsilon*2 (it's in sample_scores[ 1 ])
                                 epsilon *= epsilon_scale_value;
                                } else if( score_after_modification == sample_scores[ 0 ] ) {
                                  // Hrm.  No change.
                                  if( m_baldi_be_verbose ) {
                                    cout << "Couldn't get a score change.  ";
                                  }
  
                                  // Ok.  Try again with a bigger epsilon scale factor.
                                  epsilon_scale_value_scale_factor = 1 + ( ( epsilon_scale_value_scale_factor - 1 ) * 2 ); // TODO: DEHACKIFY MAGIC # 2.
                                  epsilon /= epsilon_scale_value;
                                  epsilon_scale_value = ( 1 / epsilon_scale_value_scale_factor );
                                  if( m_baldi_be_verbose ) {
                                    cout << "Trying again with epsilon_scale_value_scale_factor = " << epsilon_scale_value_scale_factor << ", and an epsilon_scale_value the inverse of that." << endl;
                                  }
                                } else { // score_after_modification > sample_scores[ 0 ]
                                  // Score's going up.  Good!
                                  sample_scores[ 1 ] = score_after_modification;
                                  sample_epsilons[ 1 ] = epsilon;
                                }
                        
                              } else if( sample_epsilons[ 2 ] == 0 ) {
                                assert( sample_indices[ 0 ] == 0 );
                                assert( sample_indices[ 1 ] == 1 );
                                assert( sample_indices[ 2 ] == 2 );
  
                                sample_scores[ 2 ] =
                                  score_after_modification;
                                sample_epsilons[ 2 ] =
                                  epsilon;
                              } else { // if we're still working on getting the first samples .. else ..
                        
                                assert( finding_the_peak_attempts > 2 );
  
                                // Put the new score at the end..
                                temp_index          = sample_indices[ 0 ];
                                sample_indices[ 0 ] = sample_indices[ 1 ];
                                sample_indices[ 1 ] = sample_indices[ 2 ];
                                sample_indices[ 2 ] = temp_index;
                        
                                sample_scores[ temp_index ] =
                                  score_after_modification;
                                sample_epsilons[ temp_index ] =
                                  epsilon;
                        
                              }  // End if we're still working on getting the first samples .. else ..
                        
                            } else { // if the epsilon scale value is < 1 .. else ..
                              assert( epsilon_scale_value == epsilon_scale_value_scale_factor );
  
                              if( finding_the_peak_attempts == 0 ) { // The first sample.
                                assert( sample_indices[ 0 ] == 0 );
                                assert( sample_indices[ 1 ] == 1 );
                                assert( sample_indices[ 2 ] == 2 );
  
                                sample_scores[ 0 ] =
                                  score_after_modification;
                                sample_epsilons[ 0 ] =
                                  epsilon;
                              } else if( sample_epsilons[ 1 ] == 0 ) {
                                // Set the scale value.
                                assert( sample_indices[ 0 ] == 0 );
                                assert( sample_indices[ 1 ] == 1 );
                                assert( sample_indices[ 2 ] == 2 );
                        
                                if(
                                   score_after_modification < sample_scores[ 0 ]
                                  ) {
                                  // Score went down..
                                  epsilon_scale_value =
                                    ( 1.0 / epsilon_scale_value_scale_factor );
                        
                                  // Put the new score in the front..
                                  sample_scores[ 1 ] = sample_scores[ 0 ];
                                  sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                  sample_scores[ 0 ] = score_after_modification;
                                  sample_epsilons[ 0 ] = epsilon;
                        
                                  // We've already calculated epsilon/2 (it's in sample_scores[ 1 ])
                                  epsilon *= epsilon_scale_value;
                                } else if( score_after_modification == sample_scores[ 0 ] ) {
                                  // Hrm.  No change.
                                  if( m_baldi_be_verbose ) {
                                    cout << "Couldn't get a score change.  ";
                                  }
                                  assert( epsilon_scale_value_scale_factor > 1 );
                                  if( epsilon == baldi_change_maxed_out_epsilon ) {
                                    // Erp! There's no change
                                    // because we've hit a prob=1,
                                    // maxing out what can be
                                    // accomplished by going bigger.
                                    // Try going < 1...
                                    if( m_baldi_be_verbose ) {
                                      cout << "Trying again with a flipped epsilon_scale_value." << endl;
                                    }
                                    epsilon_scale_value =
                                      ( 1.0 / epsilon_scale_value_scale_factor );
                                    // We've already calculated epsilon/2 (it's in sample_scores[ 0 ])
                                    epsilon *= epsilon_scale_value;
                                  } else {
                                    // Ok.  Try again with a bigger epsilon scale factor.
                                    epsilon_scale_value_scale_factor = 1 + ( ( epsilon_scale_value_scale_factor - 1 ) * 2 ); // TODO: DEHACKIFY MAGIC # 2.
                                    epsilon /= epsilon_scale_value;
                                    epsilon_scale_value = epsilon_scale_value_scale_factor;
                                    if( m_baldi_be_verbose ) {
                                      cout << "Trying again with epsilon_scale_value_scale_factor = " << epsilon_scale_value_scale_factor << endl;
                                    }
                                  }
                                } else { // score_after_modification > sample_scores[ 0 ]
                                  // Score's going up.  Good!
                                  sample_scores[ 1 ] = score_after_modification;
                                  sample_epsilons[ 1 ] = epsilon;
                                }
                        
                              } else if( sample_epsilons[ 2 ] == 0 ) {
                                assert( sample_indices[ 0 ] == 0 );
                                assert( sample_indices[ 1 ] == 1 );
                                assert( sample_indices[ 2 ] == 2 );
                        
                                sample_scores[ 2 ] =
                                  score_after_modification;
                                sample_epsilons[ 2 ] =
                                  epsilon;
                        
                              } else { // If still working on filling the 3 samples for the first time .. else ..
                        
                                // Put the new score at the end..
                                temp_index          = sample_indices[ 0 ];
                                sample_indices[ 0 ] = sample_indices[ 1 ];
                                sample_indices[ 1 ] = sample_indices[ 2 ];
                                sample_indices[ 2 ] = temp_index;
                        
                                sample_scores[ temp_index ] =
                                  score_after_modification;
                                sample_epsilons[ temp_index ] =
                                  epsilon;
                        
                              } // End if still working on filling the 3 samples for the first time .. else ..
                        
                            } // End if the epsilon scale value is < 1 .. else ..
                        
                            if( m_baldi_be_verbose ) {
                              cout <<
                                "CONTINUE: The three sample points are ( [ " <<
                                sample_epsilons[ sample_indices[ 0 ] ] << ", " <<
                                sample_scores[ sample_indices[ 0 ] ] << " ], [ " <<
                                sample_epsilons[ sample_indices[ 1 ] ] << ", " <<
                                sample_scores[ sample_indices[ 1 ] ] << " ], [ " <<
                                sample_epsilons[ sample_indices[ 2 ] ] << ", " <<
                                sample_scores[ sample_indices[ 2 ] ] << " ] )" << endl;
                              cout << "Scale value is " << epsilon_scale_value << 
                                ".  Local epsilon was " << epsilon << endl;
                            }
                            if(
                               ( finding_the_peak_attempts ==
                                 m_trainingParameters.siegelMaxFindingThePeakAttempts_positions )
                            ) {
                              if( m_baldi_be_verbose ) {
                                cout << "Giving up.  We've tried too many times." << endl;
                              }
                              if( sample_epsilons[ 2 ] == 0 ) {
                                if( sample_epsilons[ 1 ] == 0 ) {
                                  sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                  sample_epsilons[ 2 ] = sample_epsilons[ 0 ];
                                  sample_scores[ 1 ] = sample_scores[ 0 ];
                                  sample_scores[ 2 ] = sample_scores[ 0 ];
                                } else {
                                  sample_epsilons[ 2 ] = sample_epsilons[ 1 ];
                                  sample_scores[ 2 ] = sample_scores[ 1 ];
                                }
                              }
                              // TODO: REMOVE
                              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                cout << "?";
                              } // End if verbose
                              break;
                            } else if( ( epsilon * epsilon_scale_value ) < m_trainingParameters.siegelMinEpsilon ) {
                              if( m_baldi_be_verbose ) {
                                cout << "Giving up.  Epsilon is too small." << endl;
                              }
                              if( sample_epsilons[ 2 ] == 0 ) {
                                if( sample_epsilons[ 1 ] == 0 ) {
                                  sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                  sample_epsilons[ 2 ] = sample_epsilons[ 0 ];
                                  sample_scores[ 1 ] = sample_scores[ 0 ];
                                  sample_scores[ 2 ] = sample_scores[ 0 ];
                                } else {
                                  sample_epsilons[ 2 ] = sample_epsilons[ 1 ];
                                  sample_scores[ 2 ] = sample_scores[ 1 ];
                                }
                              }
                              // TODO: REMOVE
                              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                cout << "'";
                              } // End if verbose
                              break;
                            } else if( 
                              ( sample_epsilons[ 2 ] == 0 ) &&
                              ( finding_the_peak_attempts >=
                                m_trainingParameters.siegelMaxFindingTheGradientAttempts_positions )
                            ) {
                              if( m_baldi_be_verbose ) {
                                cout << "Giving up.  Couldn't get a score change." << endl;
                              }
                              if( sample_epsilons[ 1 ] == 0 ) {
                                sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                sample_epsilons[ 2 ] = sample_epsilons[ 0 ];
                                sample_scores[ 1 ] = sample_scores[ 0 ];
                                sample_scores[ 2 ] = sample_scores[ 0 ];
                              } else {
                                sample_epsilons[ 2 ] = sample_epsilons[ 1 ];
                                sample_scores[ 2 ] = sample_scores[ 1 ];
                              }
                              // TODO: REMOVE
                              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                cout << ";";
                              } // End if verbose
                              break;
                            } else if( ( sample_epsilons[ 2 ] != 0 ) && ( sample_scores[ 0 ] == sample_scores[ 1 ] ) && ( sample_scores[ 1 ] == sample_scores[ 2 ] ) ) {
                              // A plateau
                              if( m_baldi_be_verbose ) {
                                cout << "Plateau detected." << endl;
                              }
                              // TODO: REMOVE
                              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                cout << "`";
                              } // End if verbose
                              break;
                            } else if(
                               ( sample_epsilons[ 2 ] == 0 ) ||
                               ( sample_scores[ sample_indices[ 2 ] ] >=
                                 sample_scores[ sample_indices[ 1 ] ] )
                              ) {
                              if( m_baldi_be_verbose ) {
                                cout << "continuing to try to find the peak..." << endl;
                              }
                            } else {
                              if( m_baldi_be_verbose ) {
                                cout << "Aha!  Trapped the peak in the three sample points." << endl;
                              }
                              break;
                            }
                        
                          } // End foreach finding_the_peak_attempts : getting three sample points (trapping the peak)
                          // TODO: REMOVE
                          if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                            cout.flush();
                          } // End if verbose

                          // If we got caught on a plateau, just take whatever the last epsilon was...
                          if(
                             //( sample_scores[ 0 ] == sample_scores[ 1 ] ) &&
                             ( sample_scores[ 1 ] == sample_scores[ 2 ] )
                          ) {
                            // Do nothing.  Already, epsilon and score_after_modification are what we want.
                            if( m_baldi_be_verbose ) {
                              cout << "estimated plateau epsilon is " << epsilon << endl;
                            }
                            if( m_baldi_be_verbose ) {
                              cout << "The score using the estimated plateau epsilon is " << score_after_modification << endl;
                            }
  
                          } else if(
                             sample_scores[ sample_indices[ 2 ] ] <=
                             sample_scores[ sample_indices[ 1 ] ]
                          ) {
                            // If we trapped the peak, find the max of the quadratic..
  
                            // But first make sure they're ordered
                            // consistently low to high, for ease of
                            // coding below..
                            if( sample_epsilons[ sample_indices[ 0 ] ] > sample_epsilons[ sample_indices[ 1 ] ] ) {
                              assert( sample_epsilons[ sample_indices[ 1 ] ] >= sample_epsilons[ sample_indices[ 2 ] ] );
                              // This is the reverse order. Swap the ends..
                              temp_index          = sample_indices[ 0 ];
                              sample_indices[ 0 ] = sample_indices[ 2 ];
                              sample_indices[ 2 ] = temp_index;
                            }
                            assert( sample_epsilons[ sample_indices[ 2 ] ] >= sample_epsilons[ sample_indices[ 1 ] ] );
                            assert( sample_epsilons[ sample_indices[ 1 ] ] >= sample_epsilons[ sample_indices[ 0 ] ] );
                          
                            if( m_baldi_be_verbose ) {
                              cout <<
                                "Trapped the peak: The three sample points are ( [ " <<
                                sample_epsilons[ sample_indices[ 0 ] ] << ", " <<
                                sample_scores[ sample_indices[ 0 ] ] << " ], [ " <<
                                sample_epsilons[ sample_indices[ 1 ] ] << ", " <<
                                sample_scores[ sample_indices[ 1 ] ] << " ], [ " <<
                                sample_epsilons[ sample_indices[ 2 ] ] << ", " <<
                                sample_scores[ sample_indices[ 2 ] ] << " ] )" << endl;
                            }
                      
                            for( refining_the_peak_steps = 0;
                                  ( refining_the_peak_steps <
                                    m_trainingParameters.siegelMaxRefiningThePeakSteps_positions );
                                 ++refining_the_peak_steps
                            ) {
                          
                              // Find the epsilon at which the score is maximized, under the
                              // assumption that the landscape is locally quadratic.
                              epsilon =
                                maximizeThreePointQuadratic(
                                  sample_epsilons[ sample_indices[ 0 ] ],
                                  sample_scores[ sample_indices[ 0 ] ],
                                  sample_epsilons[ sample_indices[ 1 ] ],
                                  sample_scores[ sample_indices[ 1 ] ],
                                  sample_epsilons[ sample_indices[ 2 ] ],
                                  sample_scores[ sample_indices[ 2 ] ]
                                );
                                
                              if( m_baldi_be_verbose ) {
                                cout << "estimated peak epsilon is " << epsilon << endl;
                              }
                                  
                              for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                                if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                                    !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                                  // Then don't train it.
                                  continue;
                                }
                              
                                m_position_entente =
                                  m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                                m_position_entente.normalize(
                                  m_trainingParameters.profileValueMinimum
                                );
                              
                                if( false && m_baldi_be_verbose ) {
                                  cout << "The profile position was " << ( *m_profile )[ m_row_i - 1 ] << endl;
                                }
                                
                                m_baldi_position_boltzmann_gibbs_change =
                                  m_unconditional_baldi_position_boltzmann_gibbs_changes_vector[ m_row_i - 1 ];
                                
                                // Convert the profile position to Boltzmann-Gibbs.
                                m_baldi_position_boltzmann_gibbs.fromProfilePosition(
                                  ( *m_profile )[ m_row_i - 1 ]
                                );
                                if( false && m_baldi_be_verbose ) {
                                  cout << ".. so the corresponding m_baldi_position_boltzmann_gibbs is " << m_baldi_position_boltzmann_gibbs << endl;
                                  cout << "Scaled, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                                }
                                
                                m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                                  epsilon;
                                
                                // In the Baldi paper, this corresponds to equation
                                // 3-prime (the learning rate should incorporate the
                                // current score):
                                //m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                                //  m_startingScore_position_step; // the current score...
                                // This is wrong, but seems to do better: ! Even still, not so good.
                                //m_baldi_position_boltzmann_gibbs_change.m_scalar *=
                                //  m_startingScore_position_step; // the current score...
                                
                                if( false && m_baldi_be_verbose ) {
                                  cout << "weighted by the learning-rate, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                                }
                                
                                m_baldi_position_boltzmann_gibbs +=
                                  m_baldi_position_boltzmann_gibbs_change;
                                
                                m_baldi_position_boltzmann_gibbs.toProfilePosition(
                                  ( *m_profile )[ m_row_i - 1 ]
                                );
                                if( false && m_baldi_be_verbose ) {
                                  cout << "... so, using the estimated peak epsilon, the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                                }
                                ( *m_profile )[ m_row_i - 1 ].normalize(
                                  m_trainingParameters.profileValueMinimum
                                );
                              } // End foreach row_i
                              
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              // Make sure the scaled distributions are up-to-date.
                              for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                                m_profile->createScaledMatchDistributionForPosition(
                                  tmp_pos_i,
                                  m_scaled_match_distributions[ tmp_pos_i ]
                                );
                              } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              bool swap =
                                m_dynamic_programming.forward_score( 
                                  m_trainingParameters,
                                  false, // don't use viterbi
                                  *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                                  m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                                  m_sequences,
                                  m_sequence_count,
                                  &m_anchor_columns,
                                  &m_anchor_rows,
                                  *m_prev_forward_rows_ptr,
                                  *m_forward_rows_ptr,
                                  &sequence_scores,
                                  &largest_sequence_score,
                                  score_after_modification
                                );
                              //score_after_modification =
                              //  m_dynamic_programming.forward_score( 
                              //    m_trainingParameters,
                              //    *m_profile,
                              //    m_sequences,
                              //    m_sequence_count,
                              //    m_forward_matrices,
                              //    &sequence_scores,
                              //    &largest_sequence_score
                              //  );
                              // TODO: REMOVE?
                              assert( !isnan( score_after_modification ) );
    
                              // Return it to what it was before.
                              m_profile->copyPositions( m_unconditional_backupProfile );

                              if( m_baldi_be_verbose ) {
                                cout << "The score using the estimated peak epsilon is " << score_after_modification << endl;
                              }
                                
                              if( score_after_modification < sample_scores[ sample_indices[ 1 ] ] ) {
                                // See if we're converged..
                                if( ( ( sample_scores[ sample_indices[ 1 ] ] - score_after_modification ) / sample_scores[ sample_indices[ 1 ] ] ) < m_trainingParameters.siegelRefiningThePeakStepsConvergenceThreshold ) {
                                  break;
                                }
                                if( epsilon < sample_epsilons[ sample_indices[ 1 ] ] ) {
                                  // Put the new score at the beginning, replacing the current beginning.
                                  temp_index          = sample_indices[ 0 ];
                                } else { // the new epsilon is below the old peak .. else above..
                                  // Put the new score at the end, replacing the current end.
                                  temp_index          = sample_indices[ 2 ];
                                }
                              } else if( score_after_modification == sample_scores[ sample_indices[ 1 ] ] ) {
                                  // Then we're certainly within the convergence threhsold.  Done.
                                  break;
                              } else { // The new score is better, so should go in the middle.
                                // See if we're converged..
                                if( ( ( score_after_modification - sample_scores[ sample_indices[ 1 ] ] ) / score_after_modification ) < m_trainingParameters.siegelRefiningThePeakStepsConvergenceThreshold ) {
                                  break;
                                }
                                if( epsilon < sample_epsilons[ sample_indices[ 1 ] ] ) {
                                  // Put the new score in the middle, move the current peak up, and leave the lower value where it is.
                                  temp_index          = sample_indices[ 2 ];
                                  sample_indices[ 2 ] = sample_indices[ 1 ];
                                  sample_indices[ 1 ] = temp_index;
                                } else { // the new epsilon is below the old peak .. else above..
                                  // Put the new score in the middle, move the current peak down..
                                  temp_index          = sample_indices[ 0 ];
                                  sample_indices[ 0 ] = sample_indices[ 1 ];
                                  sample_indices[ 1 ] = temp_index;
                                }
                              } // End if the new score is worse .. else the same .. else better ..
                        
                              sample_scores[ temp_index ] =
                                score_after_modification;
                              sample_epsilons[ temp_index ] =
                                epsilon;
                            } // End foreach refining_the_peak_steps
  
                            if( refining_the_peak_steps ==
                                m_trainingParameters.siegelMaxRefiningThePeakSteps_positions ) {
                              // Then we ended because we ran out of steps, not because we actually converged.
                              // TODO: REMOVE!
                              //cout << "DIDN'T CONVERGE TO A PEAK." << endl;
                            }
                            epsilon = sample_epsilons[ sample_indices[ 1 ] ];
                            score_after_modification = sample_scores[ sample_indices[ 1 ] ];
                          } // End if we got caught on a plateau .. else if we trapped a peak ..
                        
                          // Even if we did manage to trap the peak, our new score might
                          // not be better than the old score (stranger things have
                          // happened)...
                          if( 
                             ( sample_scores[ sample_indices[ 2 ] ] >
                               sample_scores[ sample_indices[ 1 ] ] )
                              ||
                             ( score_after_modification <
                               m_startingScore_position_step )
                            ) {
                        
                            if( m_baldi_be_verbose ) {
                              cout << "Darn.  Since we couldn't find an adequate peak, we'll use the best epsilon we've seen so far." << endl;
                            }
                        
                            // In case we did find the peak and modified for it but it did no good,
                            if( 
                               sample_scores[ sample_indices[ 2 ] ] >
                               sample_scores[ sample_indices[ 1 ] ]
                              ) {
                              epsilon = sample_epsilons[ sample_indices[ 2 ] ];
                              score_after_modification = sample_scores[ sample_indices[ 2 ] ];
                            } else {
                              epsilon = sample_epsilons[ sample_indices[ 1 ] ];
                              score_after_modification = sample_scores[ sample_indices[ 1 ] ];
                            }
                        
                            if( m_baldi_be_verbose ) {
                              cout << "epsilon is now " << epsilon << endl;
                              cout << "score_after_modification is now " << score_after_modification << endl;
                            }
                          } // End if we have to fall back on the best we've found so far
                          
                        } // End if use Baldi .. else use Baldi / Siegel ..

                        if( m_baldi_be_verbose ) {
                          cout << "Score diff using espilon = " << epsilon << " is " << ( ( m_startingScore_position_step > score_after_modification ) ? "(negative) " : "" ) << ( ( m_startingScore_position_step > score_after_modification ) ? ( m_startingScore_position_step / score_after_modification ) : ( score_after_modification / m_startingScore_position_step ) ) << ": was " << m_startingScore_position_step << "; is now " << score_after_modification << endl;
                        }

                        if( !m_trainingParameters.baldiHybrid || ( m_endingScore < score_after_modification ) ) {
                          if( m_trainingParameters.baldiHybrid ) {
                            // Using Baldi or Baldi / Siegel gave us a better score than using UBW!
                            if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                              cout << "~";
                            } // End if verbose
                          }
                          
                          // Ok reapply the epsilon that worked best.
                          for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                            if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                                !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                              // Then don't train it.
                              continue;
                            }
                          
                            m_position_entente =
                              m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                            m_position_entente.normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                          
                            if( false && m_baldi_be_verbose ) {
                              cout << "The profile position was " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                            
                            m_baldi_position_boltzmann_gibbs_change =
                              m_unconditional_baldi_position_boltzmann_gibbs_changes_vector[ m_row_i - 1 ];
                            
                            // Convert the profile position to Boltzmann-Gibbs.
                            m_baldi_position_boltzmann_gibbs.fromProfilePosition(
                              ( *m_profile )[ m_row_i - 1 ]
                            );
                            if( false && m_baldi_be_verbose ) {
                              cout << ".. so the corresponding m_baldi_position_boltzmann_gibbs is " << m_baldi_position_boltzmann_gibbs << endl;
                              cout << "Scaled, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                            }
                            
                            m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                              epsilon;
                            
                            // In the Baldi paper, this corresponds to equation
                            // 3-prime (the learning rate should incorporate the
                            // current score):
                            //m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                            //  m_endingScore; // the current score...
                            // This is wrong, but seems to do better: ! Even still, not so good.
                            //m_baldi_position_boltzmann_gibbs_change.m_scalar *=
                            //  m_endingScore; // the current score...
                            
                            if( false && m_baldi_be_verbose ) {
                              cout << "weighted by epsilon, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                            }
                            
                            m_baldi_position_boltzmann_gibbs +=
                              m_baldi_position_boltzmann_gibbs_change;
                            
                            m_baldi_position_boltzmann_gibbs.toProfilePosition(
                              ( *m_profile )[ m_row_i - 1 ]
                            );
                            if( false && m_baldi_be_verbose ) {
                              cout << "... so the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                            ( *m_profile )[ m_row_i - 1 ].normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                            if( false && m_baldi_be_verbose ) {
                              cout << "Normalized, the profile position is now " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                          } // End foreach row_i
  
                          // No need to recalculate the score.  We have it.
                          m_endingScore = score_after_modification;
                        } else { // if !baldiHybird or if Baldi / Siegel gave us a better result than UBW .. else ..
                          for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                            if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                                !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                              // Then don't train it.
                              continue;
                            }
                            m_position_entente =
                              m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                            m_position_entente.normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                            ( *m_profile )[ m_row_i - 1 ] = m_position_entente; // already normalized.
                          } // End foreach row_i

                          // No need to recalculate the score; we already have it.
                          // We need score_after_modification to be updated for later use..
                          score_after_modification = m_endingScore;
                         
                          if( m_trainingParameters.baldiHybrid ) {
                            if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                              cout << " ";
                            } // End if verbose
                          }
                        } // End if !baldiHybrid or if Baldi / Siegel gave us a better result than UBW .. else ..
  
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // Make sure the scaled distribution is up-to-date.
                        if( ( m_row_i - 1 ) < ( m_profile->length() - 1 ) ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            m_row_i - 1,
                            m_scaled_match_distributions[ m_row_i - 1 ]
                          );
                        }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                      } // End if !revert_globals_for_changed_profile_length
                    } // End if baldiLearningRate > 0 ..
#endif // ALLOW_BOLTZMANN_GIBBS
                    // Note that we do the following block even if baldiLearningRate > 0, because it does the globals updating..
#ifdef ALLOW_BOLTZMANN_GIBBS
                    //if( m_trainingParameters.baldiLearningRate > 0 ) {
                    //  m_startingScore_position_step = score_after_modification;
                    //}
#endif // ALLOW_BOLTZMANN_GIBBS
                      for( m_bw_inverse_scalar =
                             m_trainingParameters.minBaumWelchInverseScalar;
                           (
                             ( m_bw_inverse_scalar <=
                               m_trainingParameters.maxBaumWelchInverseScalar ) &&
                             ( (
                                m_trainingParameters.proposeProfileLengthChanges &&
                                revert_globals_for_changed_profile_length
                                ) ? true : 
#ifdef ALLOW_BOLTZMANN_GIBBS
                               ( ( m_trainingParameters.baldiLearningRate > 0 ) ? ( m_endingScore <= score_after_modification ) : ( m_endingScore <= m_startingScore_position_step ) )
#else // !ALLOW_BOLTZMANN_GIBBS
                             ( m_endingScore <= m_startingScore_position_step )
#endif // ALLOW_BOLTZMANN_GIBBS
                               )
                           );
                           m_bw_inverse_scalar +=
                             m_trainingParameters.baumWelchInverseScalarIncrement
                      ) {
#ifdef ALLOW_BOLTZMANN_GIBBS
                        if( m_trainingParameters.baldiLearningRate > 0 ) {
                          m_endingScore = score_after_modification;
                        } else {
#endif // ALLOW_BOLTZMANN_GIBBS
                          m_endingScore = m_startingScore_position_step;
#ifdef ALLOW_BOLTZMANN_GIBBS
                        }
#endif // ALLOW_BOLTZMANN_GIBBS
                      
                        if( m_bw_inverse_scalar !=
                            m_trainingParameters.minBaumWelchInverseScalar ) {
                          m_profile->copyFrom( m_unconditional_backupProfile );
                        } // End if this is not the first attempt, restore from
                          // backups.
                      
                        if(
                          !(
                            m_trainingParameters.proposeProfileLengthChanges &&
                            revert_globals_for_changed_profile_length
                            )
#ifdef ALLOW_BOLTZMANN_GIBBS
                          && !( m_trainingParameters.baldiLearningRate > 0 )  // if we did Baldi, we've already done positions ...
#endif // ALLOW_BOLTZMANN_GIBBS
                          // TODO: REMOVE.  TESTING.
                          && ( ubw_cbw_hybrid ? ( ( m_iteration % 2 ) == 0 ) : true ) // only train positions every *other* iter
                        ) {
                          for( m_row_i = last_row; m_row_i > 0; m_row_i-- ) {
                            if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                            !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                              // Then don't train it.
                              continue;
                            }
                      
                            //if(
                            //  revert_globals_for_changed_profile_length &&
                            //  (
                            //    just_grew_profile ?
                            //    ( m_row_i == ( last_inserted_pos_i + 1 ) ) :
                            //    ( m_row_i < ( last_deleted_pos_i + 1 ) )
                            //  )
                            //) {
                            //  // Stop updating at the same place we'd stop if we
                            //  // were doing conditional bw.
                            //  break;
                            //}
                            m_position_entente =
                              m_unconditional_position_ententes_vector[ m_row_i - 1 ];
                            m_position_entente.normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                      
                            if( m_bw_inverse_scalar > 0 ) {
                              m_position_entente.m_scalar *=
                                m_bw_inverse_scalar;
                              m_position_entente.unscale();
                      
                              ( *m_profile )[ m_row_i - 1 ] += m_position_entente;
                              ( *m_profile )[ m_row_i - 1 ].normalize(
                                m_trainingParameters.profileValueMinimum
                              );
                              // We wait to update the scaled_match_distributions
                              // until we've updated the globals.
                            } else { // if m_bw_inverse_scalar > 0 .. else 
                              // TODO: REMOVE
                              //cout << "not using baldi." << endl;
                              ( *m_profile )[ m_row_i - 1 ] = m_position_entente;
                              // We wait to update the
                              // scaled_match_distributions until we've updated
                              // the globals.
                            } // End if m_bw_inverse_scalar > 0 .. else ..
                          } // End foreach row_i
                        } // End if trainProfilePositions
                      
                        // Okay, now update the globals
                        if( m_trainingParameters.trainProfileGlobals
                            // TODO: REMOVE.  TESTING.
                            && ( ubw_cbw_hybrid ? ( ( m_iteration % 2 ) == 1 ) : true ) // only train globals every *other* iter
                            ) {
                          if(
                            m_trainingParameters.proposeProfileLengthChanges &&
                            revert_globals_for_changed_profile_length
                          ) {
                            // HERE IS WHERE WE REVERT GLOBALS
                            length_at_prev_to_last_revert = length_at_last_revert;
                            length_at_last_revert = m_profile->length();
                            // For Uncondional BW, when
                            // proposeProfileLengthChanges is true, sometimes
                            // we've already reverted but we set
                            // revert_globals_for_changed_profile_length to true
                            // anyway, so that there is no BW update on the
                            // iteration in which a length change occurred.
                            if( !globals_are_at_starting_values ) {
                              // TODO: REMOVE!  Testing.
                              // TODO: REMOVE
                              //cout << "Before incorporating current globals, which are ";
                              //m_profile->writeExceptPositions( cout );
                              //cout << ", m_globalsToRevertTo's globals are ";
                              //m_globalsToRevertTo.writeExceptPositions( cout );
                              //cout << endl;
                              // Weigh it more as we progress...
                              m_profile->multiplyByExceptPositions( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) );
                              m_globalsToRevertTo.incrementByExceptPositions( *m_profile );
                              // TODO: REMOVE
                              //cout << "After, they are ";
                              //m_globalsToRevertTo.writeExceptPositions( cout );
                              //cout << endl;
                      
                              m_profile->copyExceptPositions( m_globalsToRevertTo );
                              m_profile->normalizeExceptPositions( m_trainingParameters.profileValueMinimum );
                      
                              // TODO: REMOVE
                              //cout << "globals are now " << endl;
                              //m_profile->writeExceptPositions( cout );
                              //cout << endl;
                      
                              // TODO: REMOVE?
                              if( ensure_even_indels_after_reverting ) {
                                ProbabilityType average_indel_open = 1.0;
                                average_indel_open -= ( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
                                average_indel_open /= 2;
                                ( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
                                  average_indel_open;
                                ( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
                                  average_indel_open;
#ifndef DISALLOW_FLANKING_TRANSITIONS
                                // TODO: REMOVE?!
                                ( *m_profile )[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] = average_indel_open;
                                ( *m_profile )[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] = ( 1.0 - average_indel_open );
                                ( *m_profile )[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] = average_indel_open;
                                ( *m_profile )[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] = ( 1.0 - average_indel_open );
#endif // !DISALLOW_FLANKING_TRANSITIONS
                                // TODO: REMOVE
                                //cout << "After reverting and setting indel opens to " << average_indel_open << ", globals are " << endl;
                                //m_profile->writeExceptPositions( cout );
                                //cout << endl;
                              } // End if ensure_even_indels_after_reverting
                            } else {  // If !globals_are_at_starting_values .. else ..
                              // Force an iteration without length changes.
                              if( num_iterations_until_next_length_change_proposal < 2 ) {
                                // This will be immediately decremented to 1.
                                num_iterations_until_next_length_change_proposal = 2;
                              }
                            } // End if !globals_are_at_starting_values .. else ..
                            if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                              cout << "*"; // indicating that we've reverted.
                            } // End if verbose
                            revert_globals_for_changed_profile_length = false;
                            // Don't stop training yet:
                            just_reverted_globals_for_changed_profile_length = true;
                            globals_are_at_starting_values = true;
                      
                            // only update the one time.
                            m_bw_inverse_scalar =
                              m_trainingParameters.maxBaumWelchInverseScalar;
                            if( m_trainingParameters.debug >= DEBUG_All ) {
                            //if( m_trainingParameters.verbosity >= VERBOSITY_All ) {
                              cout <<
                                "[debug] Restarting profile globals from their original values, since the profile length changed." << endl;
                              m_profile->writeExceptPositions( cout );
                              cout << endl;
                            } // End if DEBUG_All
                          } else { // if revert_globals_for_changed_profile_length .. else ..
                            if( m_bw_inverse_scalar == 0 ) {
                              m_profile->copyExceptPositions( m_global_entente );
                            } else {
                              m_global_entente.m_scalar *= m_bw_inverse_scalar;
                              m_global_entente.unscale();
                            
                              // TODO: REMOVE
                              //cout << "m_bw_inverse_scalar is " << m_bw_inverse_scalar << endl;
                              //cout << "proposed new profile globals will add ";
                              //cout << m_global_entente << " to ";
                              //m_profile->writeExceptPositions( cout );
                              //cout << ".";
                              //cout << endl;
                              m_profile->incrementByExceptPositions( m_global_entente );
                              m_profile->normalizeExceptPositions(
                                m_trainingParameters.profileValueMinimum
                              );
                              // TODO: REMOVE
                              //cout << "proposed new profile globals: ";
                              //m_profile->writeExceptPositions( cout );
                              //cout << endl;
                            } // End if m_bw_inverse_scalar == 0 .. else ..
                          } // End if revert_globals_for_changed_profile_length .. else ..
                        } // End if trainProfileGlobals
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        if( m_trainingParameters.debug >= DEBUG_All ) {
                          cout <<
                            "[debug] Before profile update, the total score is " << m_endingScore << endl;
                        } // End if DEBUG_All
                      
                          bool swap =
                            m_dynamic_programming.forward_score( 
                              m_trainingParameters,
                              false, // don't use viterbi
                              *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                              m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                              m_sequences,
                              m_sequence_count,
                              &m_anchor_columns,
                              &m_anchor_rows,
                              *m_prev_forward_rows_ptr,
                              *m_forward_rows_ptr,
                              &sequence_scores,
                              &largest_sequence_score,
                              m_endingScore
                            );
                          if( swap ) {
                            // Rotate the forward_rows
                            m_temp_forward_rows_ptr = m_forward_rows_ptr;
                            m_forward_rows_ptr = m_prev_forward_rows_ptr;
                            m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                          }
                          forward_rows_are_already_rotated_and_calculated = true;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          // TODO: REMOVE.  TESTING.
                          if( false ) {
                            cout << "using scaled_match_distributions, m_endingScore is " << m_endingScore << endl;
                            bool swap =
                              m_dynamic_programming.forward_score( 
                                m_trainingParameters,
                                false, // don't use viterbi
                                *m_profile,
                                m_sequences,
                                m_sequence_count,
                                &m_anchor_columns,
                                &m_anchor_rows,
                                *m_prev_forward_rows_ptr,
                                *m_forward_rows_ptr,
                                &sequence_scores,
                                &largest_sequence_score,
                                m_endingScore
                              );
                            if( swap ) {
                              // Rotate the forward_rows
                              m_temp_forward_rows_ptr = m_forward_rows_ptr;
                              m_forward_rows_ptr = m_prev_forward_rows_ptr;
                              m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                            }
                            forward_rows_are_already_rotated_and_calculated = true;
                            cout << "WITHOUT scaled_match_distributions, m_endingScore is " << m_endingScore << endl;
                          } // End false
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          //m_endingScore =
                          //  m_dynamic_programming.forward_score( 
                          //    m_trainingParameters,
                          //    *m_profile,
                          //    m_sequences,
                          //    m_sequence_count,
                          //    m_forward_matrices,
                          //    &sequence_scores,
                          //    &largest_sequence_score
                          //  );
                      
                          // TODO: REMOVE?
                          assert( !isnan( m_endingScore ) );
                      
                        if( ( m_trainingParameters.debug >= DEBUG_All ) ) {
                          cout <<
                            "[debug] [Profile] " <<
                            "After Unconditional Baum-Welch update, profile" <<
                            " is " << *m_profile << "." << endl;
                        } // End if DEBUG_All
                        if( ( m_trainingParameters.debug >= DEBUG_All ) ) {
                          cout <<
                            "[debug] After profile update, the total score is " << m_endingScore <<
                            endl;
                        } // End if DEBUG_All
                      } // End foreach m_bw_inverse_scalar..
                  } else { // if useUnconditionalBaumWelch .. else .. 
                    if(
                      !revert_globals_for_changed_profile_length
                    ) {
                      m_backupProfileGlobals.copyExceptPositions( *m_profile );
                    }

                    for( m_bw_inverse_scalar =
                           m_trainingParameters.minBaumWelchInverseScalar_globals;
                         (
                           ( m_bw_inverse_scalar <=
                           m_trainingParameters.maxBaumWelchInverseScalar_globals ) &&
                           ( m_endingScore <= m_startingScore_position_step )
                         );
                         m_bw_inverse_scalar +=
                           m_trainingParameters.baumWelchInverseScalarIncrement_globals
                    ) {
                      m_endingScore = m_startingScore_position_step;

                      if( m_bw_inverse_scalar !=
                          m_trainingParameters.minBaumWelchInverseScalar_globals ) {
                        // Restore old values
                        m_profile->copyExceptPositions( m_backupProfileGlobals );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // TODO: Is it necessary to recompute these here?
                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                      }

                      if(
                        revert_globals_for_changed_profile_length
                      ) {

                        length_at_prev_to_last_revert = length_at_last_revert;
                        length_at_last_revert = m_profile->length();
                        // TODO: REMOVE
                        //cout << "Before incorporating current globals, which are ";
                        //m_profile->writeExceptPositions( cout );
                        //cout << ", m_globalsToRevertTo's globals are ";
                        //m_globalsToRevertTo.writeExceptPositions( cout );
                        //cout << endl;
                        // Weigh it more as we progress...
                        m_profile->multiplyByExceptPositions( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) );
                        m_globalsToRevertTo.incrementByExceptPositions( *m_profile );
                        // TODO: REMOVE
                        //cout << "After, they are ";
                        //m_globalsToRevertTo.writeExceptPositions( cout );
                        //cout << endl;
                        m_profile->copyExceptPositions( m_globalsToRevertTo );
                        m_profile->normalizeExceptPositions( m_trainingParameters.profileValueMinimum );
                        // TODO: REMOVE
                        //cout << "globals are now " << endl;
                        //m_profile->writeExceptPositions( cout );
                        //cout << endl;

                        // TODO: REMOVE?
                        if( ensure_even_indels_after_reverting ) {
                          ProbabilityType average_indel_open = 1.0;
                          average_indel_open -= ( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
                          average_indel_open /= 2.0;
                          // TODO: REMOVE?
                          assert( average_indel_open > 0.0 );
                          ( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
                            average_indel_open;
                          ( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
                            average_indel_open;
                          // TODO: REMOVE?!
#ifndef DISALLOW_FLANKING_TRANSITIONS
                          ( *m_profile )[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] = average_indel_open;
                          ( *m_profile )[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ] = ( 1.0 - average_indel_open );
                          ( *m_profile )[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] = average_indel_open;
                          ( *m_profile )[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ] = ( 1.0 - average_indel_open );
#endif // !DISALLOW_FLANKING_TRANSITIONS
                          // TODO: REMOVE
                          //cout << "After reverting and setting indel opens to " << average_indel_open << ", globals are " << endl;
                          //m_profile->writeExceptPositions( cout );
                          //cout << endl;
                        } // End if ensure_even_indels_after_reverting
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                          cout << "*"; // indicating that we've reverted.
                        } // End if verbose

                        revert_globals_for_changed_profile_length = false;
                        // Don't stop training yet:
                        just_reverted_globals_for_changed_profile_length = true;
                        globals_are_at_starting_values = true;
                        // And don't keep going with these globals iters,
                        // just start over.
                        m_bw_inverse_scalar =
                          m_trainingParameters.maxBaumWelchInverseScalar_globals;                          
                        m_position_cycle =
                          m_trainingParameters.maxPositionCycles_globals;
                        if( m_trainingParameters.debug >= DEBUG_All ) {
                          //if( m_trainingParameters.verbosity >= VERBOSITY_All ) {
                          cout <<
                            "[debug] Restarting profile globals from their original values, since the profile length changed." << endl;
                          m_profile->writeExceptPositions( cout );
                          cout << endl;
                        } // End if DEBUG_All
                      } else if( m_bw_inverse_scalar == 0 ) {
                        if(// true ||
                        ( m_trainingParameters.debug >= DEBUG_All ) ) {
                          cout << "[debug] m_bw_inverse_scalar is 0" << endl;
                          cout << "[debug] proposed new profile globals will be ";
                          cout << m_global_entente << " (versus ";
                          m_profile->writeExceptPositions( cout );
                          cout << " ).";
                          cout << endl;
                        } // End if DEBUG_All

//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//                  // TODO: REMOVE!!!! TESTING!!!!
//                  m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] = 0;
//                  m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] = 0;
//                  //m_global_entente[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] = 0;
//                  m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ] = 0;
//                  //m_global_entente[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ] = 0;
//                  m_global_entente.normalize( 0 );
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        m_profile->copyExceptPositions( m_global_entente );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // TODO: REMOVE!!!! TESTING!!!!
                        //( *m_profile )[ Transition::fromMatch ] = m_backupProfileGlobals[ Transition::fromMatch ];

                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                      } else { // if revert_globals_for_changed_profile_length .. else if m_bw_inverse_scalar == 0 .. else ..
                        m_global_entente.m_scalar *=
                          ( m_bw_inverse_scalar );
                        m_global_entente.unscale();

                        if( //true ||
                        ( m_trainingParameters.debug >= DEBUG_All ) ) {
                          cout << "[debug] m_bw_inverse_scalar is " << m_bw_inverse_scalar << endl;
                          cout << "[debug] proposed new profile globals will add ";
                          cout << m_global_entente << " to ";
                          m_profile->writeExceptPositions( cout );
                          cout << ".";
                          cout << endl;
                        } // End if DEBUG_All

                        m_profile->incrementByExceptPositions( m_global_entente );
                        m_profile->normalizeExceptPositions(
                          m_trainingParameters.profileValueMinimum
                        );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        if( //true ||
                          ( m_trainingParameters.debug >= DEBUG_All ) ) {
                          cout << "[debug ] proposed new profile globals: ";
                          m_profile->writeExceptPositions( cout );
                          cout << endl;
                        } // End if DEBUG_All
                      } // End if revert_globals_for_changed_profile_length .. else if m_bw_inverse_scalar == 0 .. else ..

                      if( m_trainingParameters.debug >= DEBUG_All ) {
                        cout <<
                          "[debug] Before profile globals update, the total score is " << m_endingScore << endl;
                      } // End if DEBUG_All
                      
                        // TODO: REMOVE
                        //cout << "calling m_dynamic_programming.forward_score(..)" << endl;
                        //cout << "m_trainingParameters are " << m_trainingParameters << endl;
                        //cout << "*m_profile is " << *m_profile << endl;
                        //cout << "m_sequences are:" << endl;
                        //for( uint32_t i = 0; i < m_sequences.size(); i++ ) {
                        //  cout << m_sequences[ i ] << endl;
                        //}
                        //cout << "m_sequence_count is " << m_sequence_count << endl;
                        //cout << "[skipping forward matrices]" << endl;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                      // TODO: REMOVE! TESTING!
                      // Revert Del-in values to the backup values
                      //( *m_profile )[ Transition::fromBegin ] =
                      //  m_backupProfileGlobals[ Transition::fromBegin ];
                      //( *m_profile )[ Transition::fromDeletionIn ] =
                      //  m_backupProfileGlobals[ Transition::fromDeletionIn ];
                      // TODO: REMOVE! TESTING!
                      // Revert Del-out values to the backup values
                      //( *m_profile )[ Transition::fromMatch ] =
                      //  m_backupProfileGlobals[ Transition::fromMatch ];
                      // TODO: REMOVE! Testing
                      // Revert just the fromMatch::toDeletionOut prob.
                      //( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] = 0;
                      //( *m_profile )[ Transition::fromMatch ].normalize( 0 );
                      //( *m_profile )[ Transition::fromMatch ] *=
                      //  ( 1.0 - m_backupProfileGlobals[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] );
                      //( *m_profile )[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ] =
                      //  m_backupProfileGlobals[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ];

                      //( *m_profile )[ Transition::fromDeletionOut ] =
                      //  m_backupProfileGlobals[ Transition::fromDeletionOut ];

                      // Make sure the scaled distributions are up-to-date.
                      //for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                      //  m_profile->createScaledMatchDistributionForPosition(
                      //    tmp_pos_i,
                      //    m_scaled_match_distributions[ tmp_pos_i ]
                      //  );
                      //} // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

                        bool swap =
                          m_dynamic_programming.forward_score( 
                            m_trainingParameters,
                            false, // don't use viterbi
                            *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                            m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                            m_sequences,
                            m_sequence_count,
                            &m_anchor_columns,
                            &m_anchor_rows,
                            *m_prev_forward_rows_ptr,
                            *m_forward_rows_ptr,
                            &sequence_scores,
                            &largest_sequence_score,
                            m_endingScore
                          );
                        if( swap ) {
                          // Rotate the forward_rows
                          m_temp_forward_rows_ptr = m_forward_rows_ptr;
                          m_forward_rows_ptr = m_prev_forward_rows_ptr;
                          m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                        }
                        forward_rows_are_already_rotated_and_calculated = true;
                        //m_endingScore =
                        //  m_dynamic_programming.forward_score( 
                        //    m_trainingParameters,
                        //    *m_profile,
                        //    m_sequences,
                        //    m_sequence_count,
                        //    m_forward_matrices,
                        //    &sequence_scores,
                        //    &largest_sequence_score
                        //  );

                        // TODO: REMOVE?
                        assert( !isnan( m_endingScore ) );
                        if( //true ||
                          ( m_trainingParameters.debug >= DEBUG_All ) ) {
                        cout <<
                          "[debug] [Profile Globals] " <<
                          "After Baum-Welch update, profile" <<
                          " is " << *m_profile << "." << endl;
                      } // End if DEBUG_All
                      if( m_trainingParameters.debug >= DEBUG_All ) {
                        cout <<
                          "[debug] After profile globals update, the total score is " << m_endingScore <<
                          endl;
                      } // End if DEBUG_All

                    } // End foreach m_bw_inverse_scalar..
                  } // End if useUnconditionalBaumWelch .. else ..
              } else if(
                ( m_row_i > 0 ) &&
                ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
                m_trainingParameters.trainProfilePositions
              ) {
                if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                    !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                  // We're skipping (parent profile) the training for this position ..
                  if( m_trainingParameters.useUnconditionalBaumWelch ) {
                    // TODO: REMOVE?  Unnecessary.
                    m_unconditional_position_ententes_vector[ m_row_i - 1 ].zero();
                  } // End if useUnconditionalBaumWelch
                } else if( m_trainingParameters.useUnconditionalBaumWelch ) {
                  // TODO: REMOVE
                  //cout << "SAVING for row " << m_row_i << ": pos entente is " << m_position_entente << endl;
                  m_unconditional_position_ententes_vector[ m_row_i - 1 ] =
                    m_position_entente;
#ifdef ALLOW_BOLTZMANN_GIBBS
                  if( m_trainingParameters.baldiLearningRate > 0 ) {
                    m_unconditional_baldi_position_boltzmann_gibbs_changes_vector[ m_row_i - 1 ] =
                      m_baldi_position_boltzmann_gibbs_change;
                  } // End if baldiLearningRate > 0
#endif // ALLOW_BOLTZMANN_GIBBS
                } else {
                  m_backupProfilePosition = ( *m_profile )[ m_row_i - 1 ];
                  m_backupEntentePosition = m_position_entente;
                  // Also save the sequence scores.
                  for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                    backup_sequence_scores[ m_seq_i ] =
                      sequence_scores[ m_seq_i ];
                  } // End foreach m_seq_i, copy its score to the backup.
                } // End if we are skipping (parent profile) training for this position .. else if useUnconditionalBaumWelch .. else ..
                if( m_trainingParameters.useUnconditionalBaumWelch ) {
                  // Then don't actually modify the profile
                } else { // if useUnconditionalBaumWelch .. else ..
                  for( m_bw_inverse_scalar =
                         m_trainingParameters.minBaumWelchInverseScalar;
                       ( ( m_bw_inverse_scalar <=
                           m_trainingParameters.maxBaumWelchInverseScalar ) &&
                         ( m_endingScore <= m_startingScore_position_step ) );
                       m_bw_inverse_scalar +=
                         m_trainingParameters.baumWelchInverseScalarIncrement
                  ) {
                    m_endingScore = m_startingScore_position_step;
                    if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                        !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                      // We're skipping (parent profile) the training for this position ..
                      // TODO: REMOVE
                      //cout << "SKIPPING PROFLE at position " << ( m_row_i - 1 ) << endl;
                    } else {
                      ( *m_profile )[ m_row_i - 1 ] = m_backupProfilePosition;
                      m_position_entente = m_backupEntentePosition;
                      // TODO: REMOVE
                      //m_position_entente /= m_position_entente.m_scalar;
                      //cout << "Unscaled, m_position_entente is " << m_position_entente << endl;
                      m_position_entente.normalize(
                        m_trainingParameters.profileValueMinimum
                      );
                      if( m_bw_inverse_scalar > 0 ) {
                          m_position_entente.m_scalar *=
                            m_bw_inverse_scalar;
                          m_position_entente.unscale();
                          ( *m_profile )[ m_row_i - 1 ] += m_position_entente;
                          ( *m_profile )[ m_row_i - 1 ].normalize(
                            m_trainingParameters.profileValueMinimum
                          );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          // Make sure the scaled distribution is up-to-date.
                          if( ( m_row_i - 1 ) < ( m_profile->length() - 1 ) ) {
                            m_profile->createScaledMatchDistributionForPosition(
                              m_row_i - 1,
                              m_scaled_match_distributions[ m_row_i - 1 ]
                            );
                          }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

                      } else { // if m_bw_inverse_scalar > 0 .. else (it's 0)
#ifdef ALLOW_BOLTZMANN_GIBBS
                        if( m_trainingParameters.baldiLearningRate > 0 ) { // Do Baldi or Baldi / Siegel ..

                          m_startingScore_position_step = m_endingScore;

                          if( m_baldi_be_verbose ) {
                            cout << "The profile position was " << ( *m_profile )[ m_row_i - 1 ] << endl;
                          }
                          // Convert the profile position to Boltzmann-Gibbs.
                          m_baldi_position_boltzmann_gibbs.fromProfilePosition(
                            ( *m_profile )[ m_row_i - 1 ]
                          );
                          if( m_baldi_be_verbose ) {
                            cout << ".. so the corresponding m_baldi_position_boltzmann_gibbs is " << m_baldi_position_boltzmann_gibbs << endl;
                            cout << "Scaled, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                          }

                          if( m_trainingParameters.baldiHybrid ) {
                            // First determine the result if we were just using straight-up CBW.
                            ( *m_profile )[ m_row_i - 1 ] = m_position_entente; // already normalized.
                            
                            m_endingScore = 1.0;
                            for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                              m_endingScore *=
                                m_coefficients_vector[ m_seq_i ].calculateScore( 
                                  ( *m_profile )[ m_row_i - 1 ]
                                );
                            } // End foreach seq_i, recompute the score.

                            // Now we have the score to beat, plus a good starting place to determine the optimal epsilon.
                            // TODO: Could I get a starting value for epsilon from this?
                          } // End if m_trainingParameters.baldiHybrid

                          if( m_trainingParameters.siegelMaxFindingThePeakAttempts_positions == 0 ) {
                            // Just use straight-up Baldi.
                            m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                              m_trainingParameters.baldiLearningRate;
                            if( m_baldi_be_verbose ) {
                              cout << "weighted by the learning-rate, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                            }

                            // In the Baldi paper, this corresponds to equation
                            // 3-prime (the learning rate should incorporate the
                            // current score):
                            //m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                            //  m_startingScore_position_step; // the current score...
                            // This is wrong, but seems to do better: ! Even still, not so good.
                            //m_baldi_position_boltzmann_gibbs_change.m_scalar *=
                            //  m_startingScore_position_step; // the current score...

                            if( m_baldi_be_verbose ) {
                              cout << "weighted by the learning-rate and by current score m_startingScore_position_step, m_baldi_position_boltzmann_gibbs_change is " << m_baldi_position_boltzmann_gibbs_change << endl;
                            }
                          
                            m_baldi_position_boltzmann_gibbs +=
                              m_baldi_position_boltzmann_gibbs_change;
                          
                            m_baldi_position_boltzmann_gibbs.toProfilePosition(
                              ( *m_profile )[ m_row_i - 1 ]
                            );
                            if( m_baldi_be_verbose ) {
                              cout << "... so the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                            ( *m_profile )[ m_row_i - 1 ].normalize(
                              m_trainingParameters.profileValueMinimum
                            );

                            score_after_modification = 1.0;
                            for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                              score_after_modification *=
                                m_coefficients_vector[ m_seq_i ].calculateScore( 
                                  ( *m_profile )[ m_row_i - 1 ]
                                );
                            } // End foreach seq_i, recompute the score.
                          } else { // End if use Baldi .. else use Baldi / Siegel ..

                            // Get three sample points for which the
                            // score of the third is less than that of
                            // the second.  (Capture the peak.)
                            sample_indices[ 0 ] = 0;
                            sample_indices[ 1 ] = 1;
                            sample_indices[ 2 ] = 2;
                            sample_scores[ 0 ] = 0;
                            sample_scores[ 1 ] = 0;
                            sample_scores[ 2 ] = 0;
                            sample_epsilons[ 0 ] = 0;
                            sample_epsilons[ 1 ] = 0;
                            sample_epsilons[ 2 ] = 0;
                            epsilon_scale_value_scale_factor = m_trainingParameters.siegelEpsilonScaleFactor;
                            baldi_change_maxed_out_epsilon = 0; // TODO: REMOVE. TESTING.
                            for( finding_the_peak_attempts = 0,
                                   epsilon = ( 1.0 / ( epsilon_scale_value_scale_factor * epsilon_scale_value_scale_factor * epsilon_scale_value_scale_factor ) ),
                                   epsilon_scale_value = epsilon_scale_value_scale_factor;
                                   (
                                    ( epsilon >= m_trainingParameters.siegelMinEpsilon ) &&
                                    ( finding_the_peak_attempts <=
                                      m_trainingParameters.siegelMaxFindingThePeakAttempts_positions ) &&
                                    !(
                                      ( sample_epsilons[ 2 ] == 0 ) &&
                                      ( finding_the_peak_attempts >
                                        m_trainingParameters.siegelMaxFindingTheGradientAttempts_positions )
                                    ) &&
                                    (
                                     ( sample_epsilons[ 2 ] == 0 ) ||
                                      ( sample_scores[ sample_indices[ 2 ] ] >=
                                        sample_scores[ sample_indices[ 1 ] ] )
                                    )
                                   );
                                 ++finding_the_peak_attempts,
                                   ( epsilon *= epsilon_scale_value )
                            ) { // for( finding_the_peak_attempts ... )
                          
                              if( m_baldi_be_verbose ) {
                                cout << " Conditional Baldi / Siegel with finding-the-peak with epsilon = " << epsilon << endl;
                              }
                              m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                                epsilon;
                            
                              m_baldi_position_boltzmann_gibbs +=
                                m_baldi_position_boltzmann_gibbs_change;

                              // Save the position first so we can see if we got a change.
                              baldi_last_profile_position = ( *m_profile )[ m_row_i - 1 ];
                            
                              m_baldi_position_boltzmann_gibbs.toProfilePosition(
                                ( *m_profile )[ m_row_i - 1 ]
                              );
                              if( m_baldi_be_verbose ) {
                                cout << "... so the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                              }

                              // Check to see if we've maxed out our ability to change the profile position.
                              if( ( ( *m_profile )[ m_row_i - 1 ][ Emission::Match ].maximumValue() == 1 ) &&
                                  // Also check that the change to that residue is to increase it...
                                  ( m_baldi_position_boltzmann_gibbs_change[ Emission::Match ][ ( *m_profile )[ m_row_i - 1 ][ Emission::Match ].maximumValueType() ] > 0 )
                              ) {
                                baldi_change_maxed_out_epsilon = epsilon;
                                if( m_baldi_be_verbose ) {
                                  cout << "Baldi/Siegel change maxed out the profile position with prob 1 on " << ( *m_profile )[ m_row_i - 1 ][ Emission::Match ].maximumValueType() << "!" << endl;
                                }
                              }
                              ( *m_profile )[ m_row_i - 1 ].normalize(
                                m_trainingParameters.profileValueMinimum
                              );

                              // Check to see if we actually changed the profile
                              if( ( *m_profile )[ m_row_i - 1 ][ Emission::Match ] == baldi_last_profile_position[ Emission::Match ] ) {
                                // We don't need to recalculate the score.
                                baldi_change_had_an_effect = false;
                                if( m_baldi_be_verbose ) {
                                  cout << "Baldi/Siegel change had no effect on the profile position!" << endl;
                                }
                              } else {
                                baldi_change_had_an_effect = true;
                                score_after_modification = 1.0;
                                for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                                  score_after_modification *=
                                    m_coefficients_vector[ m_seq_i ].calculateScore( 
                                      ( *m_profile )[ m_row_i - 1 ]
                                    );
                                } // End foreach seq_i, recompute the score.
                              } // End if the profile position didn't change at all .. else ..

                              // Return it to what it was before.
                              m_baldi_position_boltzmann_gibbs -=
                                m_baldi_position_boltzmann_gibbs_change;
                              m_baldi_position_boltzmann_gibbs_change.m_scalar *=
                                epsilon;
                          
                              if( m_baldi_be_verbose ) {
                                cout << "Got score diff " << ( ( m_startingScore_position_step > score_after_modification ) ? "(negative) " : "" ) << ( ( m_startingScore_position_step > score_after_modification ) ? ( m_startingScore_position_step / score_after_modification ) : ( score_after_modification / m_startingScore_position_step ) ) << " after changing the position by " <<
                                  m_baldi_position_boltzmann_gibbs_change.createUnscaledCopy() << ": was " << m_startingScore_position_step << "; is now " << score_after_modification << endl;
                              }

                              if( epsilon_scale_value < 1 ) {
                                assert( epsilon_scale_value == ( 1.0 / epsilon_scale_value_scale_factor ) );
                                if( finding_the_peak_attempts == 0 ) {
                                  assert( false ); // The code expects that epsilon_scale_value > 1 to start with...
                                  cout << "ERROR: finding_the_peak_attempts == 0 but epsilon_scale_value < 1!" << endl;
                                  exit( 1 );
                                  //assert( sample_indices[ 0 ] == 0 );
                                  //assert( sample_indices[ 1 ] == 1 );
                                  //assert( sample_indices[ 2 ] == 2 );
                                  //
                                  //sample_scores[ 0 ] =
                                  //  score_after_modification;
                                  //sample_epsilons[ 0 ] =
                                  //  epsilon;
                                } else if( sample_epsilons[ 1 ] == 0 ) {
                                  // Set the scale value.
                                  assert( sample_indices[ 0 ] == 0 );
                                  assert( sample_indices[ 1 ] == 1 );
                                  assert( sample_indices[ 2 ] == 2 );
                          
                                  if(
                                     score_after_modification < sample_scores[ 0 ]
                                    ) {
                                    // Score went down.  Flip.
                                    epsilon_scale_value = epsilon_scale_value_scale_factor;
                          
                                    // Put the new score in the front..
                                    sample_scores[ 1 ] = sample_scores[ 0 ];
                                    sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                    sample_scores[ 0 ] = score_after_modification;
                                    sample_epsilons[ 0 ] = epsilon;
                          
                                    // We've already calculated epsilon*2 (it's in sample_scores[ 1 ])
                                   epsilon *= epsilon_scale_value;
                                  } else if( score_after_modification == sample_scores[ 0 ] ) {
                                    // Hrm.  No change.
                                    if( m_baldi_be_verbose ) {
                                      cout << "Couldn't get a score change.  ";
                                    }

                                    // Ok.  Try again with a bigger epsilon scale factor.
                                    epsilon_scale_value_scale_factor = 1 + ( ( epsilon_scale_value_scale_factor - 1 ) * 2 ); // TODO: DEHACKIFY MAGIC # 2.
                                    epsilon /= epsilon_scale_value;
                                    epsilon_scale_value = ( 1 / epsilon_scale_value_scale_factor );
                                    if( m_baldi_be_verbose ) {
                                      cout << "Trying again with epsilon_scale_value_scale_factor = " << epsilon_scale_value_scale_factor << ", and an epsilon_scale_value the inverse of that." << endl;
                                    }
                                  } else { // score_after_modification > sample_scores[ 0 ]
                                    // Score's going up.  Good!
                                    sample_scores[ 1 ] = score_after_modification;
                                    sample_epsilons[ 1 ] = epsilon;
                                  }
                          
                                } else if( sample_epsilons[ 2 ] == 0 ) {
                                  assert( sample_indices[ 0 ] == 0 );
                                  assert( sample_indices[ 1 ] == 1 );
                                  assert( sample_indices[ 2 ] == 2 );

                                  sample_scores[ 2 ] =
                                    score_after_modification;
                                  sample_epsilons[ 2 ] =
                                    epsilon;
                                } else { // if we're still working on getting the first samples .. else ..
                          
                                  assert( finding_the_peak_attempts > 2 );

                                  // Put the new score at the end..
                                  temp_index          = sample_indices[ 0 ];
                                  sample_indices[ 0 ] = sample_indices[ 1 ];
                                  sample_indices[ 1 ] = sample_indices[ 2 ];
                                  sample_indices[ 2 ] = temp_index;
                          
                                  sample_scores[ temp_index ] =
                                    score_after_modification;
                                  sample_epsilons[ temp_index ] =
                                    epsilon;
                          
                                }  // End if we're still working on getting the first samples .. else ..
                          
                              } else { // if the epsilon scale value is < 1 .. else ..
                                assert( epsilon_scale_value == epsilon_scale_value_scale_factor );

                                if( finding_the_peak_attempts == 0 ) { // The first sample.
                                  assert( sample_indices[ 0 ] == 0 );
                                  assert( sample_indices[ 1 ] == 1 );
                                  assert( sample_indices[ 2 ] == 2 );

                                  sample_scores[ 0 ] =
                                    score_after_modification;
                                  sample_epsilons[ 0 ] =
                                    epsilon;
                                } else if( sample_epsilons[ 1 ] == 0 ) {
                                  // Set the scale value.
                                  assert( sample_indices[ 0 ] == 0 );
                                  assert( sample_indices[ 1 ] == 1 );
                                  assert( sample_indices[ 2 ] == 2 );
                          
                                  if(
                                     score_after_modification < sample_scores[ 0 ]
                                    ) {
                                    // Score went down..
                                    epsilon_scale_value =
                                      ( 1.0 / epsilon_scale_value_scale_factor );
                          
                                    // Put the new score in the front..
                                    sample_scores[ 1 ] = sample_scores[ 0 ];
                                    sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                    sample_scores[ 0 ] = score_after_modification;
                                    sample_epsilons[ 0 ] = epsilon;
                          
                                    // We've already calculated epsilon/2 (it's in sample_scores[ 1 ])
                                    epsilon *= epsilon_scale_value;
                                  } else if( score_after_modification == sample_scores[ 0 ] ) {
                                    // Hrm.  No change.
                                    if( m_baldi_be_verbose ) {
                                      cout << "Couldn't get a score change.  ";
                                    }
                                    assert( epsilon_scale_value_scale_factor > 1 );
                                    if( epsilon == baldi_change_maxed_out_epsilon ) {
                                      // Erp! There's no change
                                      // because we've hit a prob=1,
                                      // maxing out what can be
                                      // accomplished by going bigger.
                                      // Try going < 1...
                                      if( m_baldi_be_verbose ) {
                                        cout << "Trying again with a flipped epsilon_scale_value." << endl;
                                      }
                                      epsilon_scale_value =
                                        ( 1.0 / epsilon_scale_value_scale_factor );
                                      // We've already calculated epsilon/2 (it's in sample_scores[ 0 ])
                                      epsilon *= epsilon_scale_value;
                                    } else {
                                      // Ok.  Try again with a bigger epsilon scale factor.
                                      epsilon_scale_value_scale_factor = 1 + ( ( epsilon_scale_value_scale_factor - 1 ) * 2 ); // TODO: DEHACKIFY MAGIC # 2.
                                      epsilon /= epsilon_scale_value;
                                      epsilon_scale_value = epsilon_scale_value_scale_factor;
                                      if( m_baldi_be_verbose ) {
                                        cout << "Trying again with epsilon_scale_value_scale_factor = " << epsilon_scale_value_scale_factor << endl;
                                      }
                                    }
                                  } else { // score_after_modification > sample_scores[ 0 ]
                                    // Score's going up.  Good!
                                    sample_scores[ 1 ] = score_after_modification;
                                    sample_epsilons[ 1 ] = epsilon;
                                  }
                          
                                } else if( sample_epsilons[ 2 ] == 0 ) {
                                  assert( sample_indices[ 0 ] == 0 );
                                  assert( sample_indices[ 1 ] == 1 );
                                  assert( sample_indices[ 2 ] == 2 );
                          
                                  sample_scores[ 2 ] =
                                    score_after_modification;
                                  sample_epsilons[ 2 ] =
                                    epsilon;
                          
                                } else { // If still working on filling the 3 samples for the first time .. else ..
                          
                                  // Put the new score at the end..
                                  temp_index          = sample_indices[ 0 ];
                                  sample_indices[ 0 ] = sample_indices[ 1 ];
                                  sample_indices[ 1 ] = sample_indices[ 2 ];
                                  sample_indices[ 2 ] = temp_index;
                          
                                  sample_scores[ temp_index ] =
                                    score_after_modification;
                                  sample_epsilons[ temp_index ] =
                                    epsilon;
                          
                                } // End if still working on filling the 3 samples for the first time .. else ..
                          
                              } // End if the epsilon scale value is < 1 .. else ..
                          
                              if( m_baldi_be_verbose ) {
                                cout <<
                                  "CONTINUE: The three sample points are ( [ " <<
                                  sample_epsilons[ sample_indices[ 0 ] ] << ", " <<
                                  sample_scores[ sample_indices[ 0 ] ] << " ], [ " <<
                                  sample_epsilons[ sample_indices[ 1 ] ] << ", " <<
                                  sample_scores[ sample_indices[ 1 ] ] << " ], [ " <<
                                  sample_epsilons[ sample_indices[ 2 ] ] << ", " <<
                                  sample_scores[ sample_indices[ 2 ] ] << " ] )" << endl;
                                cout << "Scale value is " << epsilon_scale_value << 
                                  ".  Local epsilon was " << epsilon << endl;
                              }
                              if(
                                 ( finding_the_peak_attempts ==
                                   m_trainingParameters.siegelMaxFindingThePeakAttempts_positions )
                              ) {
                                if( m_baldi_be_verbose ) {
                                  cout << "Giving up.  We've tried too many times." << endl;
                                }
                                if( sample_epsilons[ 2 ] == 0 ) {
                                  if( sample_epsilons[ 1 ] == 0 ) {
                                    sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                    sample_epsilons[ 2 ] = sample_epsilons[ 0 ];
                                    sample_scores[ 1 ] = sample_scores[ 0 ];
                                    sample_scores[ 2 ] = sample_scores[ 0 ];
                                  } else {
                                    sample_epsilons[ 2 ] = sample_epsilons[ 1 ];
                                    sample_scores[ 2 ] = sample_scores[ 1 ];
                                  }
                                }
                                // TODO: REMOVE
                                if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                  cout << "?";
                                } // End if verbose
                                break;
                              } else if( ( epsilon * epsilon_scale_value ) < m_trainingParameters.siegelMinEpsilon ) {
                                if( m_baldi_be_verbose ) {
                                  cout << "Giving up.  Epsilon is too small." << endl;
                                }
                                if( sample_epsilons[ 2 ] == 0 ) {
                                  if( sample_epsilons[ 1 ] == 0 ) {
                                    sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                    sample_epsilons[ 2 ] = sample_epsilons[ 0 ];
                                    sample_scores[ 1 ] = sample_scores[ 0 ];
                                    sample_scores[ 2 ] = sample_scores[ 0 ];
                                  } else {
                                    sample_epsilons[ 2 ] = sample_epsilons[ 1 ];
                                    sample_scores[ 2 ] = sample_scores[ 1 ];
                                  }
                                }
                                // TODO: REMOVE
                                if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                  cout << "'";
                                } // End if verbose
                                break;
                              } else if( 
                                ( sample_epsilons[ 2 ] == 0 ) &&
                                ( finding_the_peak_attempts >=
                                  m_trainingParameters.siegelMaxFindingTheGradientAttempts_positions )
                              ) {
                                if( m_baldi_be_verbose ) {
                                  cout << "Giving up.  Couldn't get a score change." << endl;
                                }
                                if( sample_epsilons[ 1 ] == 0 ) {
                                  sample_epsilons[ 1 ] = sample_epsilons[ 0 ];
                                  sample_epsilons[ 2 ] = sample_epsilons[ 0 ];
                                  sample_scores[ 1 ] = sample_scores[ 0 ];
                                  sample_scores[ 2 ] = sample_scores[ 0 ];
                                } else {
                                  sample_epsilons[ 2 ] = sample_epsilons[ 1 ];
                                  sample_scores[ 2 ] = sample_scores[ 1 ];
                                }
                                // TODO: REMOVE
                                if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                  cout << ";";
                                } // End if verbose
                                break;
                              } else if( ( sample_epsilons[ 2 ] != 0 ) && ( sample_scores[ 0 ] == sample_scores[ 1 ] ) && ( sample_scores[ 1 ] == sample_scores[ 2 ] ) ) {
                                // A plateau
                                if( m_baldi_be_verbose ) {
                                  cout << "Plateau detected." << endl;
                                }
                                // TODO: REMOVE
                                if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                  cout << "`";
                                } // End if verbose
                                break;
                              } else if(
                                 ( sample_epsilons[ 2 ] == 0 ) ||
                                 ( sample_scores[ sample_indices[ 2 ] ] >=
                                   sample_scores[ sample_indices[ 1 ] ] )
                                ) {
                                if( m_baldi_be_verbose ) {
                                  cout << "continuing to try to find the peak..." << endl;
                                }
                              } else {
                                if( m_baldi_be_verbose ) {
                                  cout << "Aha!  Trapped the peak in the three sample points." << endl;
                                }
                                break;
                              }
                          
                            } // End foreach finding_the_peak_attempts : getting three sample points (trapping the peak)
                            // TODO: REMOVE
                            if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                              cout.flush();
                            } // End if verbose

                            // If we got caught on a plateau, just take whatever the last epsilon was...
                            if(
                               //( sample_scores[ 0 ] == sample_scores[ 1 ] ) &&
                               ( sample_scores[ 1 ] == sample_scores[ 2 ] )
                            ) {
                              // Do nothing.  Already, epsilon and score_after_modification are what we want.
                              if( m_baldi_be_verbose ) {
                                cout << "estimated plateau epsilon is " << epsilon << endl;
                              }
                              if( m_baldi_be_verbose ) {
                                cout << "The score using the estimated plateau epsilon is " << score_after_modification << endl;
                              }
  
                            } else if(
                               sample_scores[ sample_indices[ 2 ] ] <=
                               sample_scores[ sample_indices[ 1 ] ]
                            ) {
                              // If we trapped the peak, find the max of the quadratic..

                              // But first make sure they're ordered
                              // consistently low to high, for ease of
                              // coding below..
                              if( sample_epsilons[ sample_indices[ 0 ] ] > sample_epsilons[ sample_indices[ 1 ] ] ) {
                                assert( sample_epsilons[ sample_indices[ 1 ] ] >= sample_epsilons[ sample_indices[ 2 ] ] );
                                // This is the reverse order. Swap the ends..
                                temp_index          = sample_indices[ 0 ];
                                sample_indices[ 0 ] = sample_indices[ 2 ];
                                sample_indices[ 2 ] = temp_index;
                              }
                              assert( sample_epsilons[ sample_indices[ 2 ] ] >= sample_epsilons[ sample_indices[ 1 ] ] );
                              assert( sample_epsilons[ sample_indices[ 1 ] ] >= sample_epsilons[ sample_indices[ 0 ] ] );
                            
                              if( m_baldi_be_verbose ) {
                                cout <<
                                  "Trapped the peak: The three sample points are ( [ " <<
                                  sample_epsilons[ sample_indices[ 0 ] ] << ", " <<
                                  sample_scores[ sample_indices[ 0 ] ] << " ], [ " <<
                                  sample_epsilons[ sample_indices[ 1 ] ] << ", " <<
                                  sample_scores[ sample_indices[ 1 ] ] << " ], [ " <<
                                  sample_epsilons[ sample_indices[ 2 ] ] << ", " <<
                                  sample_scores[ sample_indices[ 2 ] ] << " ] )" << endl;
                              }

                              for( refining_the_peak_steps = 0;
                                    ( refining_the_peak_steps <
                                      m_trainingParameters.siegelMaxRefiningThePeakSteps_positions );
                                   ++refining_the_peak_steps
                              ) {
                            
                                // Find the epsilon at which the score is maximized, under the
                                // assumption that the landscape is locally quadratic.
                                epsilon =
                                  maximizeThreePointQuadratic(
                                    sample_epsilons[ sample_indices[ 0 ] ],
                                    sample_scores[ sample_indices[ 0 ] ],
                                    sample_epsilons[ sample_indices[ 1 ] ],
                                    sample_scores[ sample_indices[ 1 ] ],
                                    sample_epsilons[ sample_indices[ 2 ] ],
                                    sample_scores[ sample_indices[ 2 ] ]
                                  );
                                  
                                if( m_baldi_be_verbose ) {
                                  cout << "estimated peak epsilon is " << epsilon << endl;
                                }
                                  
                                //if( isnan( epsilon ) ) { // Doesn't work!  Why?  Workaround:
                                assert( !( epsilon != epsilon ) );
                                
                                m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                                  epsilon;
                                  
                                m_baldi_position_boltzmann_gibbs +=
                                  m_baldi_position_boltzmann_gibbs_change;
                                  
                                m_baldi_position_boltzmann_gibbs.toProfilePosition(
                                  ( *m_profile )[ m_row_i - 1 ]
                                );
                                if( m_baldi_be_verbose ) {
                                  cout << "... so, with the estimated peak epsilon, the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                                }
                                ( *m_profile )[ m_row_i - 1 ].normalize(
                                  m_trainingParameters.profileValueMinimum
                                );
                                  
                                score_after_modification = 1.0;
                                for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                                  score_after_modification *=
                                    m_coefficients_vector[ m_seq_i ].calculateScore( 
                                      ( *m_profile )[ m_row_i - 1 ]
                                    );
                                } // End foreach seq_i, recompute the score.
                                  
                                // Return it to what it was before.
                                m_baldi_position_boltzmann_gibbs -=
                                  m_baldi_position_boltzmann_gibbs_change;
                                m_baldi_position_boltzmann_gibbs_change.m_scalar *=
                                  epsilon;
                                  
                                if( m_baldi_be_verbose ) {
                                  cout << "The score using the estimated peak epsilon is " << score_after_modification << endl;
                                }
                                  
                                if( score_after_modification < sample_scores[ sample_indices[ 1 ] ] ) {
                                  // See if we're converged..
                                  if( ( ( sample_scores[ sample_indices[ 1 ] ] - score_after_modification ) / sample_scores[ sample_indices[ 1 ] ] ) < m_trainingParameters.siegelRefiningThePeakStepsConvergenceThreshold ) {
                                    break;
                                  }
                                  if( epsilon < sample_epsilons[ sample_indices[ 1 ] ] ) {
                                    // Put the new score at the beginning, replacing the current beginning.
                                    temp_index          = sample_indices[ 0 ];
                                  } else { // the new epsilon is below the old peak .. else above..
                                    // Put the new score at the end, replacing the current end.
                                    temp_index          = sample_indices[ 2 ];
                                  }
                                } else if( score_after_modification == sample_scores[ sample_indices[ 1 ] ] ) {
                                    // Then we're certainly within the convergence threhsold.  Done.
                                    break;
                                } else { // The new score is better, so should go in the middle.
                                  // See if we're converged..
                                  if( ( ( score_after_modification - sample_scores[ sample_indices[ 1 ] ] ) / score_after_modification ) < m_trainingParameters.siegelRefiningThePeakStepsConvergenceThreshold ) {
                                    break;
                                  }
                                  if( epsilon < sample_epsilons[ sample_indices[ 1 ] ] ) {
                                    // Put the new score in the middle, move the current peak up, and leave the lower value where it is.
                                    temp_index          = sample_indices[ 2 ];
                                    sample_indices[ 2 ] = sample_indices[ 1 ];
                                    sample_indices[ 1 ] = temp_index;
                                  } else { // the new epsilon is below the old peak .. else above..
                                    // Put the new score in the middle, move the current peak down..
                                    temp_index          = sample_indices[ 0 ];
                                    sample_indices[ 0 ] = sample_indices[ 1 ];
                                    sample_indices[ 1 ] = temp_index;
                                  }
                                } // End if the new score is worse .. else the same .. else better ..
                          
                                sample_scores[ temp_index ] =
                                  score_after_modification;
                                sample_epsilons[ temp_index ] =
                                  epsilon;
                              } // End foreach refining_the_peak_steps

                              if( refining_the_peak_steps ==
                                  m_trainingParameters.siegelMaxRefiningThePeakSteps_positions ) {
                                // Then we ended because we ran out of steps, not because we actually converged.
                                // TODO: REMOVE!
                                //cout << "DIDN'T CONVERGE TO A PEAK." << endl;
                              }
                              epsilon = sample_epsilons[ sample_indices[ 1 ] ];
                              score_after_modification = sample_scores[ sample_indices[ 1 ] ];
                            } // End if we got caught on a plateau .. else if we trapped a peak ..
                          
                            // Even if we did manage to trap the peak, our new score might
                            // not be better than the old score (stranger things have
                            // happened)...
                            if( 
                               ( sample_scores[ sample_indices[ 2 ] ] >
                                 sample_scores[ sample_indices[ 1 ] ] )
                                ||
                               ( score_after_modification <
                                 m_startingScore_position_step )
                              ) {
                          
                              if( m_baldi_be_verbose ) {
                                cout << "Darn.  Since we couldn't find an adequate peak, we'll use the best epsilon we've seen so far." << endl;
                              }
                          
                              // In case we did find the peak and modified for it but it did no good,
                              if( 
                                 sample_scores[ sample_indices[ 2 ] ] >
                                 sample_scores[ sample_indices[ 1 ] ]
                                ) {
                                epsilon = sample_epsilons[ sample_indices[ 2 ] ];
                                score_after_modification = sample_scores[ sample_indices[ 2 ] ];
                              } else {
                                epsilon = sample_epsilons[ sample_indices[ 1 ] ];
                                score_after_modification = sample_scores[ sample_indices[ 1 ] ];
                              }
                          
                              if( m_baldi_be_verbose ) {
                                cout << "epsilon is now " << epsilon << endl;
                                cout << "score_after_modification is now " << score_after_modification << endl;
                              }
                            } // End if we have to fall back on the best we've found so far
                          
                            m_baldi_position_boltzmann_gibbs_change.m_scalar /=
                              epsilon;
                          
                            m_baldi_position_boltzmann_gibbs +=
                              m_baldi_position_boltzmann_gibbs_change;
                              
                            m_baldi_position_boltzmann_gibbs.toProfilePosition(
                              ( *m_profile )[ m_row_i - 1 ]
                            );
                            if( false && m_baldi_be_verbose ) {
                              cout << "... so the profile position becomes " << ( *m_profile )[ m_row_i - 1 ] << endl;
                            }
                            ( *m_profile )[ m_row_i - 1 ].normalize(
                              m_trainingParameters.profileValueMinimum
                            );
                          } // End if use Baldi .. else use Baldi / Siegel ..

                          if( m_baldi_be_verbose ) {
                            cout << "Score diff using espilon = " << epsilon << ": after changing the position by " <<
                              m_baldi_position_boltzmann_gibbs_change.createUnscaledCopy() << " is " << ( ( m_startingScore_position_step > score_after_modification ) ? "(negative) " : "" ) << ( ( m_startingScore_position_step > score_after_modification ) ? ( m_startingScore_position_step / score_after_modification ) : ( score_after_modification / m_startingScore_position_step ) ) << ": was " << m_startingScore_position_step << "; is now " << score_after_modification << endl;
                          }

                          if( m_trainingParameters.baldiHybrid ) {
                            // Also make sure we couldn't do better just using the CBW update.
                            if( m_endingScore < score_after_modification ) {
                              // Using Baldi / Siegel gave us a better score than using CBW!
                              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                cout << "~";
                              } // End if verbose
                            
                              m_endingScore = score_after_modification;
                            } else { // if Baldi / Siegel gave us a better result than CBW .. else ..
                              // We are better off just using the usual CBW update.
                              ( *m_profile )[ m_row_i - 1 ] = m_position_entente;
                              
                              if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                                cout << " ";
                              } // End if verbose
                            } // End if Baldi / Siegel gave us a better result than CBW .. else ..
                          } else { // if m_trainingParameters.baldiHybrid .. else ..
                            m_endingScore = score_after_modification;
                          } // End if m_trainingParameters.baldiHybrid .. else ..

                          // TODO: REMOVE
//#ifndef NDEBUG
//                          // TODO: The problem with trying to call this is that the forward_score methods need to be updated to take in the scaled distros..
//                          // Make sure the score is what we say it is.
//                          score_after_modification =
//                         m_dynamic_programming.forward_score( 
//                           m_trainingParameters,
//                           *m_profile,
//                           m_sequences,
//                           m_sequence_count,
//                           m_forward_matrices,
//                           NULL,
//                           NULL
//                         );
//                          
//#endif // !NDEBUG

                          if( m_baldi_be_verbose ) {
                            cout << "Normalized, the profile position is now " << ( *m_profile )[ m_row_i - 1 ] << endl;
                          }
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          // Make sure the scaled distribution is up-to-date.
                          if( ( m_row_i - 1 ) < ( m_profile->length() - 1 ) ) {
                            m_profile->createScaledMatchDistributionForPosition(
                              m_row_i - 1,
                              m_scaled_match_distributions[ m_row_i - 1 ]
                            );
                          }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        } else // if baldiLearningRate > 0 .. else ..
#endif // ALLOW_BOLTZMANN_GIBBS
                        {
                          // TODO: REMOVE
                          //cout << "not using baldi." << endl;
                          //cout << "Before, profile position " << ( m_row_i - 1 ) << " is " << ( *m_profile )[ m_row_i - 1 ] << endl;
                          //if( ( m_row_i - 1 ) >= 1 ) {
                          //  cout << "Prev profile position " << ( m_row_i - 2 ) << " is " << ( *m_profile )[ m_row_i - 2 ] << endl;
                          //}
                          //if( m_row_i < m_profile->length() ) {
                          //  cout << "Next profile position " << ( m_row_i ) << " is " << ( *m_profile )[ m_row_i ] << endl;
                          //}
                          //cout << "Setting profile position " << ( m_row_i - 1 ) << " to " << m_position_entente << endl;
                          ( *m_profile )[ m_row_i - 1 ] = m_position_entente;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          // Make sure the scaled distribution is up-to-date.
                          if( ( m_row_i - 1 ) < ( m_profile->length() - 1 ) ) {
                            m_profile->createScaledMatchDistributionForPosition(
                              m_row_i - 1,
                              m_scaled_match_distributions[ m_row_i - 1 ]
                            );
                          }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                        } // End if baldiLearningRate > 0 .. else ..
                      } // End if m_bw_inverse_scalar > 0 .. else 
                    } // End if we are skipping (parent profile) training for this position .. else ..
                  
                    if( m_trainingParameters.debug >= DEBUG_All ) {
                      if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                          !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                        // We're skipping the (parent profile) training for this
                        // position ..
                      } else {
                        cout <<
                          "[debug] [Row " << m_row_i << " (using scalar = (1/" << m_bw_inverse_scalar << ")] " <<
                          "After Baum-Welch update, profile position " << (m_row_i - 1) <<
                          " is " << ( *m_profile )[ m_row_i - 1 ] << "." << endl;
                      } // End if we are skipping (parent profile) training for this position .. else ..
                    } // End if DEBUG_All

                    // Recalculate the score using the coefficients vectors.
                    // Note that (even if we accept the change) we don't need
                    // to change the forward rows (since we're moving backwards
                    // through the rows), nor do we need to change the backward
                    // rows (since we do that on the next round, after
                    // decrementing m_row_i).
                    if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                        !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                      // We're skipping the (parent profile) training for this
                      // position ..
                      //m_position_entente.zero();
                    } else {
                      largest_sequence_score = numeric_limits<double>::min();
                      m_endingScore = 1.0;
                      for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                        sequence_score =
                          m_coefficients_vector[ m_seq_i ].calculateScore( 
                            ( *m_profile )[ m_row_i - 1 ]
                          );
  
                        sequence_scores[ m_seq_i ] = sequence_score;
                        if( sequence_score > largest_sequence_score ) {
                          largest_sequence_score = sequence_score;
                        }
                          
                        if( m_trainingParameters.debug >= DEBUG_All ) {
                          cout <<
                            "[debug] [Row " << m_row_i << " (using scalar = (1/" << m_bw_inverse_scalar << ")] " <<
                            "After Baum-Welch update, sequence " << m_seq_i <<
                            " score is " << sequence_score << "." << endl;
                        } // End if DEBUG_All
  
                        m_endingScore *= sequence_score;
                      } // End foreach sequence, update forward, backward, entente,
                        // and recalculate ending score.
  
                      // TODO: REMOVE?
                      assert( !isnan( m_endingScore ) );
                    } // End if we are skipping (parent profile) training for this position .. else ..
                  } // End foreach m_bw_inverse_scalar
                } // End if useUnconditionalBaumWelch .. else ..
              } // End if ( m_trainingPhase != TRAINING_PHASE_Positions ) .. else ..
              
              if(
                (
                  ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
                  m_trainingParameters.trainProfilePositions &&
                  !m_trainingParameters.useUnconditionalBaumWelch
                ) ||
                ( m_row_i == 0 )
              ) {
                if( ( m_trainingParameters.debug > DEBUG_None ) ) {

                  cout << "[debug] Before this training step, the score was " << m_startingScore_position_step << endl;
                  cout << "[debug] After this training step, the score is " << m_endingScore << endl;
                } // End if debug
                if(
                  !(
                    ( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) ||
                    // Always accept after reverting:
                    just_reverted_globals_for_changed_profile_length
                  )
                  &&
                  ( m_endingScore < m_startingScore_position_step )
                ) {
                  // Reject it.

                  // TODO: REMOVE?  Debugging problem where the score goes down by a clearly non-neglible amount, violating the proof of the EM.
                  //assert( percentChange( m_endingScore, m_startingScore_position_step ) > -100 );
                  // TODO: REMOVE
                  //if( percentChange( m_endingScore, m_startingScore_position_step ) <= -100 ) {
                  // TODO: REMOVE
                  //cout << "REJECTING: (" << m_endingScore << "<" << m_startingScore_position_step << ")" << endl;
                  //cout << "           Percent Change is " << percentChange( m_endingScore, m_startingScore_position_step ) << endl;
                  // TODO: REMOVE
                  //if(
                  //  ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
                  //  m_trainingParameters.trainProfilePositions &&
                  //  !m_trainingParameters.useUnconditionalBaumWelch
                  //) {
                  //  cout << "m_position_entente was " << m_position_entente << endl;
                  //}

                  //  if( m_trainingParameters.useUnconditionalBaumWelch ) {
                  //    cout << "We tried this profile:\n" << ( *m_profile ) << endl;
                  //    cout << "We will revert to\n" << m_unconditional_backupProfile << endl;
                  //    //m_backupProfileGlobals.copyExceptPositions( m_unconditional_backupProfile );
                  //    //cout << "We will revert to " << m_backupProfileGlobals << endl;
                  //  } else if( m_trainingPhase == TRAINING_PHASE_Globals ) {
                  //    cout << "We tried these globals: ";
                  //    m_profile->writeExceptPositions( cout );
                  //    cout << endl;
                  //    cout << "We will revert to " << m_backupProfileGlobals << endl;
                  //  }
                  // TODO: REMOVE?
                  //assert( percentChange( m_endingScore, m_startingScore_position_step ) > -100 );
                  //}
                    //}
          
                  // TODO: REMOVE (?)
                  if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                    // TODO: REMOVE
                    if( ( m_scorePercentChange = percentChange( m_endingScore, m_startingScore_position_step ) ) > -100.0 ) {
                      cout << "=";
                    } else {
                      cout << "-";
                    }
                    // TODO: REMOVE
                    //cout << "(diff:" << ( m_startingScore_position_step - m_endingScore ) << "[" << m_scorePercentChange << "%])";
                    // TODO: REMOVE
                    if( m_scorePercentChange <= -100.0 ) {
                      assert( m_scorePercentChange >= -100.0 );
                      cout << "(diff(" << m_endingScore << "," << m_startingScore_position_step << "):-" << ( m_startingScore_position_step - m_endingScore ) << "/" << m_startingScore_position_step << "[" << m_scorePercentChange << "%])";
                      // TODO: PUT BACK
                      //exit( 1 );
                      // TODO: REMOVE
                      if( !m_trainingParameters.usePriors && ( m_trainingPhase == TRAINING_PHASE_Globals ) && !m_trainingParameters.useUnconditionalBaumWelch ) {
                        cout << "Attempted new profile globals are ";
                        m_profile->writeExceptPositions( cout );
                        cout << endl;
                        cout << "m_global_entente is " << m_global_entente << endl;
                        cout << "\tunscaled, m_global_entente is " << m_global_entente.createUnscaledCopy() << endl;
                        cout << "Old profile globals were ";
                        cout << m_backupProfileGlobals;
                        cout << endl;
                        // TODO:  REMOVE.  Making sure that the new score was right.
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//                        // Make sure the scaled distributions are up-to-date.
//                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
//                          m_profile->createScaledMatchDistributionForPosition(
//                            tmp_pos_i,
//                            m_scaled_match_distributions[ tmp_pos_i ]
//                          );
//                        } // End foreach tmp_pos_i
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
//                        bool swap =
//                          m_dynamic_programming.forward_score( 
//                            m_trainingParameters,
//                            false, // don't use viterbi
//                            *m_profile,
//  #if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//                            m_scaled_match_distributions,
//  #endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
//                            m_sequences,
//                            m_sequence_count,
//                            NULL,
//                            NULL,
//                            *m_prev_forward_rows_ptr,
//                            *m_forward_rows_ptr,
//                            NULL,
//                            NULL,
//                            score_after_modification
//                          );
//                        cout << "Recalculated new score is " << score_after_modification << endl;
//                        // TODO: REMOVE.  Making sure that the prev score was right.
//                        m_profile->copyExceptPositions( m_backupProfileGlobals );
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//                        // Make sure the scaled distributions are up-to-date.
//                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
//                          m_profile->createScaledMatchDistributionForPosition(
//                            tmp_pos_i,
//                            m_scaled_match_distributions[ tmp_pos_i ]
//                          );
//                        } // End foreach tmp_pos_i
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
//                        swap =
//                          m_dynamic_programming.forward_score( 
//                            m_trainingParameters,
//                            false, // don't use viterbi
//                            *m_profile,
//  #if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//                            m_scaled_match_distributions,
//  #endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
//                            m_sequences,
//                            m_sequence_count,
//                            NULL,
//                            NULL,
//                            *m_prev_forward_rows_ptr,
//                            *m_forward_rows_ptr,
//                            NULL,
//                            NULL,
//                            score_after_modification
//                          );
//                        cout << "Recalculated old score is " << score_after_modification << endl;
                      } // End if !usePriors & TRAINING_PHASE_Globals & !useUnconditionalBaumWelch
                      assert( false );
                    }
                  } // End if verbose
                  
                  if( m_trainingParameters.debug > DEBUG_None ) {
                    cout << "[debug] Rejecting and reverting, since the score went down." << endl;
                  } // End if debug
                  if( m_trainingPhase != TRAINING_PHASE_Positions ) {
                      if( m_trainingParameters.useUnconditionalBaumWelch ) {
                        if( m_trainingParameters.debug >= DEBUG_All ) {
                          cout <<
                            "[debug] Before undoing profile update, the total score is " << m_endingScore <<
                            endl;
                        } // End if DEBUG_All
              
                        m_profile->copyFrom( m_unconditional_backupProfile );
                      } else { // if useUnconditionalBaumWelch .. else ..
                        if( m_trainingParameters.debug >= DEBUG_All ) {
                          cout <<
                            "[debug] Before undoing profile globals update, the total score is " << m_endingScore <<
                            endl;
                        } // End if DEBUG_All
              
                        m_profile->copyExceptPositions( m_backupProfileGlobals );
                      } // End if useUnconditionalBaumWelch .. else ..
              
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                      bool swap =
                        m_dynamic_programming.forward_score( 
                          m_trainingParameters,
                          false, // don't use viterbi
                          *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                          m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                          m_sequences,
                          m_sequence_count,
                          &m_anchor_columns,
                          &m_anchor_rows,
                          *m_prev_forward_rows_ptr,
                          *m_forward_rows_ptr,
                          &sequence_scores,
                          &largest_sequence_score,
                          m_endingScore
                        );
                      if( swap ) {
                        // Rotate the forward_rows
                        m_temp_forward_rows_ptr = m_forward_rows_ptr;
                        m_forward_rows_ptr = m_prev_forward_rows_ptr;
                        m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                      }
                      forward_rows_are_already_rotated_and_calculated = true;
                      //m_endingScore =
                      //  m_dynamic_programming.forward_score( 
                      //    m_trainingParameters,
                      //    *m_profile,
                      //    m_sequences,
                      //    m_sequence_count,
                      //    m_forward_matrices,
                      //    &sequence_scores,
                      //    &largest_sequence_score
                      //  );

                      // TODO: REMOVE?
                      assert( !isnan( m_endingScore ) );

                      if( m_trainingParameters.debug >= DEBUG_All ) {
                        cout <<
                          "[debug] After undoing profile globals update, the total score is " << m_endingScore <<
                          endl;
                      } // End if DEBUG_All
              
                  } else { // ( m_trainingPhase == TRAINING_PHASE_Positions )
                    if( m_trainingParameters.debug >= DEBUG_All ) {
                      cout << "[debug] m_bw_inverse_scalar was " << ( m_bw_inverse_scalar - m_trainingParameters.baumWelchInverseScalarIncrement ) << endl;
                    } // End if DEBUG_All
                    if( ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
                        !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
                      // We're skipping the (parent profile) training for this
                      // position ..
                    } else {
                      ( *m_profile )[ m_row_i - 1 ] = m_backupProfilePosition;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                      // Make sure the scaled distribution is up-to-date.
                      if( ( m_row_i - 1 ) < ( m_profile->length() - 1 ) ) {
                        m_profile->createScaledMatchDistributionForPosition(
                          m_row_i - 1,
                          m_scaled_match_distributions[ m_row_i - 1 ]
                        );
                      }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

                      // Also restore the sequence scores.
                      for( m_seq_i = 0; m_seq_i < m_sequence_count; m_seq_i++ ) {
                        sequence_scores[ m_seq_i ] =
                          backup_sequence_scores[ m_seq_i ];
                      } // End foreach m_seq_i, restore its score from the backup.
                    } // End if we are skipping training for this position .. else ..
                    m_endingScore = m_startingScore_position_step;
                  } // End if ( m_trainingPhase != TRAINING_PHASE_Positions ) .. else .. 
                } else { // ( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) || just_reverted_globals_for_changed_profile_length || ( m_endingScore >= m_startingScore_position_step )
                  // Accept changes
                  if(
                    m_trainingParameters.trainProfileGlobals &&
                    ( m_trainingPhase == TRAINING_PHASE_Globals ) &&
                    !just_reverted_globals_for_changed_profile_length
                  ) {
                    globals_are_at_starting_values = false;
                  }
          
                  // TODO: REMOVE
                  //if( m_trainingPhase == TRAINING_PHASE_Globals ) {
                  //  cout << "ACCEPTING: (" << m_endingScore << ">" << m_startingScore_position_step << ")" << endl;
                  //}
          
                  // TODO: REMOVE
                  //if( isnan( ( *m_profile )[ ( m_row_i - 1 ) ][ Emission::Match ].total() ) ) {
                  //  cout << "====== total is nan!" << endl;
                  //  cout << "====== Profile at row " << ( m_row_i - 1 ) << " is " << ( *m_profile )[ ( m_row_i - 1 ) ] << endl;
                  //  exit( 0 );
                  //}

                  // TODO: REMOVE (?)
                  if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
                    if( !( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) || ( m_endingScore >= m_startingScore_position_step ) ) {
                      cout << "+";
                    } // End if !( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) || ( m_endingScore >= m_startingScore_position_step )
                  } // End if verbose
                  
                  if( m_trainingParameters.debug > DEBUG_None ) {
                    if( just_reverted_globals_for_changed_profile_length ) {
                      cout << "[debug] Accepting, since we just reverted the globals." << endl;
                    } else if( !( !disallow_alwaysAccept && m_trainingParameters.alwaysAccept ) || ( m_endingScore >= m_startingScore_position_step ) ) {
                      cout << "[debug] Accepting, since the score went up." << endl;
                    } else {
                      cout << "[debug] Accepting, since we always accept." << endl;
                    }
                  } // End if debug
                  if( m_trainingParameters.debug >= DEBUG_All ) {
                    if( m_trainingPhase == TRAINING_PHASE_Positions ) {
                      cout << "[debug] m_bw_inverse_scalar was " << ( m_bw_inverse_scalar - m_trainingParameters.baumWelchInverseScalarIncrement ) << endl;
                    } else if( m_trainingPhase == TRAINING_PHASE_Globals ) {
                      cout << "[debug] m_bw_inverse_scalar was " << ( m_bw_inverse_scalar - m_trainingParameters.baumWelchInverseScalarIncrement_globals ) << endl;
                    } else {
                      // Erp!?  This assertion will fail
                      assert( m_trainingPhase < trainingPhaseCount );
                    }
                  } // End if DEBUG_All

                } // End if rejecting .. else accepting ...
                if( m_trainingParameters.debug > DEBUG_None ) {
                  if( ( m_trainingPhase != TRAINING_PHASE_Positions ) ) {
                      cout <<
                        "[debug] [Profile Globals] ";
                  } else {
                    cout <<
                      "[debug] [Row " << m_row_i << "] ";
                  }
                  cout <<
                    "After Baum-Welch update, total" <<
                    " score is " << m_endingScore << "." << endl;
                } // End if debug
              } // End if this isn't a global training step with an internal
                // row, update the profile.
          
              if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                  ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
                  m_trainingParameters.verbosity >= VERBOSITY_All ) {
                cout << ".";
                cout.flush();
              }
          
              // TODO: REMOVE!
              //if( ( m_row_i > 0 ) &&
              //    ( m_trainingParameters.positionShouldBeTrained != NULL ) &&
              //    !( *m_trainingParameters.positionShouldBeTrained )[ m_row_i - 1 ] ) {
              //  cout << "THE SKIPPED PROFILE POSITION ( " << ( m_row_i - 1 ) << " ) is:" << endl;
              //  cout << ( *m_profile )[ m_row_i - 1 ] << endl;
              //}
          
          
            } while( m_row_i-- > 0 ); // End for each forward row, downto 0 (incl. 0)

                  //// TODO: REMOVE!!  TESTING!!
                  //cout << "The profile is now: " << endl << *m_profile << endl;
          
            if(
              ( m_trainingPhase == TRAINING_PHASE_Positions ) &&
              !m_trainingParameters.useUnconditionalBaumWelch &&
              !revert_globals_for_changed_profile_length
            ) {
              // Recalculate the forward matrices.
              if( m_trainingParameters.debug >= DEBUG_All ) {
                cout <<
                  "[debug] Before recalculating, the total score is " << m_endingScore <<
                  endl;
              } // End if DEBUG_All
              // TODO: Is recalculation necessary?

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                        // Make sure the scaled distributions are up-to-date.
                        for( tmp_pos_i = 0; tmp_pos_i < ( m_profile->length() - 1 ); tmp_pos_i++ ) {
                          m_profile->createScaledMatchDistributionForPosition(
                            tmp_pos_i,
                            m_scaled_match_distributions[ tmp_pos_i ]
                          );
                        } // End foreach tmp_pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                bool swap =
                  m_dynamic_programming.forward_score( 
                    m_trainingParameters,
                    false, // don't use viterbi
                    *m_profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                    m_scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
                    m_sequences,
                    m_sequence_count,
                    &m_anchor_columns,
                    &m_anchor_rows,
                    *m_prev_forward_rows_ptr,
                    *m_forward_rows_ptr,
                    &sequence_scores,
                    &largest_sequence_score,
                    m_endingScore
                  );
                if( swap ) {
                  // Rotate the forward_rows
                  m_temp_forward_rows_ptr = m_forward_rows_ptr;
                  m_forward_rows_ptr = m_prev_forward_rows_ptr;
                  m_prev_forward_rows_ptr = m_temp_forward_rows_ptr;
                }
                forward_rows_are_already_rotated_and_calculated = true;
                //m_endingScore =
                //  m_dynamic_programming.forward_score( 
                //    m_trainingParameters,
                //    *m_profile,
                //    m_sequences,
                //    m_sequence_count,
                //    m_forward_matrices,
                //    &sequence_scores,
                //    &largest_sequence_score
                //  );

                // TODO: REMOVE?
                assert( !isnan( m_endingScore ) );
              if( m_trainingParameters.debug >= DEBUG_All ) {
                cout <<
                  "[debug] After recalculating, the total score is " << m_endingScore <<
                  endl;
              } // End if DEBUG_All
            } // End if we need to recalculate as this is the end of the
              // position-specifics in this position cycle
          
            if( m_trainingParameters.useUnconditionalBaumWelch &&
                ( m_trainingPhase == TRAINING_PHASE_Positions ) ) {
              if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                  ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
                  m_trainingParameters.verbosity >= VERBOSITY_All ) {
                cout << "]";
                if( m_trainingParameters.verbosity >=
                           VERBOSITY_Medium ) {
                  cout << endl;
                } else {
                  cout.flush();
                }
              }
            } else { // if useUnconditionalBaumWelch && TRAINING_PHASE_Positions .. else .. 

              // profile distance is the squared euclidean distance divided by
              // the number of free parameters. Note that the "euclidean
              // distance" is calculated in each multinomial distribution using
              // all but the last probability; this is an arbitrary way to
              // compensate for the fact that those probabilities live in a
              // simplex (so the dimensions are not orthogonal), but note that
              // the result would be slightly different if we used a different
              // subset of the parameters in each MultinomialDistribution.
              
              if( m_trainingPhase == TRAINING_PHASE_Globals ) {
                profile_distance_position_cycle =
                  ( m_startingProfile_position_cycle.euclideanDistanceSquaredExceptPositions( *m_profile ) / m_profile->freeParameterCountExceptPositions() );
              } else {
                  if(
                    m_trainingParameters.proposeProfileLengthChanges &&
                    ( m_startingProfile_position_cycle.length() != m_profile->length() )
                  ) {
                    profile_distance_position_cycle =
                      numeric_limits<double>::max();
                  } else {
                    // profile_distance is the squared euclidean distance /
                    // m_profile->length.
                    if( m_trainingParameters.trainProfileGlobals ) {
                      profile_distance_position_cycle =
                        ( m_startingProfile_position_cycle.euclideanDistanceSquared( *m_profile ) / m_profile->freeParameterCount() );
                    } else { // if trainProfileGlobals .. else ..
                      profile_distance_position_cycle =
                        ( m_startingProfile_position_cycle.euclideanDistanceSquaredPositions( *m_profile ) / m_profile->freeParameterCountPositions() );
                    } // End if trainProfileGlobals .. else ..
                  } // End if the profile has changed length .. else .. 

              } // End if( m_trainingPhase == TRAINING_PHASE_Globals ) .. else 
              // TODO: REMOVE
              //cout << "AT THE END OF THE POSITION CYCLE, THE PROFILE DISTANCE IS " << profile_distance_position_cycle << endl;

              if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
                  ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
                  m_trainingParameters.verbosity >= VERBOSITY_All ) {
                if( m_trainingParameters.useUnconditionalBaumWelch ) {
                  if(
                    profile_length_changed_this_iteration
                  ) {
                    cout << "score:" << m_endingScore << ";";
                    cout << "length:" << m_profile->length() << "]";
                  } else {
                    cout << "score:" << m_endingScore << "]";
                  }
                } else if(
                  profile_length_changed_this_iteration &&
                  ( m_trainingPhase == TRAINING_PHASE_Positions )
                ) {
                  if( !revert_globals_for_changed_profile_length ) {
                    cout << "score:" << m_endingScore << ";";
                  }
                  cout << "length:" << m_profile->length() << "]";
                } else {
                  cout << "score:" << m_endingScore << ";distance:" << profile_distance_position_cycle << "]";
                }
                if( m_trainingParameters.verbosity >=
                           VERBOSITY_Medium ) {
                  cout << endl;
                } else {
                  cout.flush();
                }
              }
            } // End if we are doing unconditional bw and we have just
              // finished with the positions (not having come to the globals
              // yet, thus not having changed the score at all) .. else ..
          } // End for each position_cycle

        } // End for each training phase (positions or globals or seq ids)

        // profile distance is the squared euclidean distance divided by
        // the number of free parameters. Note that the "euclidean
        // distance" is calculated in each multinomial distribution using
        // all but the last probability; this is an arbitrary way to
        // compensate for the fact that those probabilities live in a
        // simplex (so the dimensions are not orthogonal), but note that
        // the result would be slightly different if we used a different
        // subset of the parameters in each MultinomialDistribution.
        
          // profile_distance is the squared euclidean distance /
          // m_profile->length.
          if(
            m_trainingParameters.proposeProfileLengthChanges &&
            ( m_startingProfile_iteration.length() != m_profile->length() )
          ) {
            profile_distance_iteration =
              numeric_limits<double>::max();
          } else {
            if( m_trainingParameters.trainProfileGlobals ) {
              profile_distance_iteration =
                ( m_startingProfile_iteration.euclideanDistanceSquared( *m_profile ) / m_profile->freeParameterCount() );
            } else {
              profile_distance_iteration =
                ( m_startingProfile_iteration.euclideanDistanceSquaredPositions( *m_profile ) / m_profile->freeParameterCountPositions() );
            } // End if trainProfileGlobals .. else ..

            // When the distance gets very small, numerical drift takes over.
            // To prevent it, don't alwaysAccept when the distance is small.
            // TODO: REMOVE? Testing.
            if(
              profile_distance_iteration <
              disallow_alwaysAccept_profile_distance_iteration_threshold
            ) {
              disallow_alwaysAccept = true;
            } else {
              disallow_alwaysAccept = false;
            }
          } // End if the profile length has changed .. else ..

        // TODO: REMOVE
        //cout << "AT THE END OF ITERATION " << ( ( ( ubw_cbw_hybrid ) ? static_cast<int>( m_iteration / 2 ) : m_iteration ) + 1 ) << ", THE PROFILE DISTANCE IS " << profile_distance_iteration << endl;

        if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
            ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
            m_trainingParameters.verbosity >= VERBOSITY_All ) {
          cout << "distance:" << profile_distance_iteration << "}";
          cout.flush();
        }
      } // End for each iteration
  
      // Yay!  Done!
      if( ( m_trainingParameters.verbosity >= VERBOSITY_Low ) &&
          ( m_trainingParameters.verbosity <  VERBOSITY_High ) ||
          m_trainingParameters.verbosity >= VERBOSITY_All ) {
        cout << "done." << endl;
      }

      // TODO: REMOVE.  TESTING.
      //cout << "percentChange: " << m_scorePercentChange << endl;

      //m_scoreImprovement =
      //  m_endingScore - m_scoreBeforeTraining;
    
      if( m_trainingParameters.verbosity >= VERBOSITY_Low ) {
        cout << "Done training.  Final score is " << m_endingScore << ".  Starting score was " << m_scoreBeforeTraining << "." << endl;
        if( m_trainingParameters.proposeProfileLengthChanges ) {
          cout << "\tFinal length is " << m_profile->length() << "." << endl;
        }
      } // End if verbose

      return m_endingScore;
    } // train()

  template <class ProfileType,
            class ScoreType,
            class MatrixValueType,
            class SequenceResidueType,
            class InternalNodeType>
  GALOSH_INLINE_TRAIN_PERCENT_CHANGE
  double
  ProfileTrainer<ProfileType, ScoreType, MatrixValueType, SequenceResidueType, InternalNodeType>::
    // Updated to return the percent change of the log score (when we are using
    // LogProbability-type scores).
  percentChange ( ScoreType const & new_val, ScoreType const & old_val )
  {
      if( m_trainingParameters.debug >= DEBUG_All ) {
       
          const double pct_change =
            ( 100.0 * ( toDouble( new_val / old_val ) - 1.0 ) );
        cout << "PERCENT_CHANGE(" << new_val << "," << old_val << "): " << pct_change << endl;
        return pct_change;
      } else {
        return ( 100.0 * ( toDouble( new_val / old_val ) - 1.0 ) );
      } // End if DEBUG_All .. else ..
    } // percentChange( ScoreType const &, ScoreType const & )

} // End namespace galosh

#endif // __GALOSH_PROFILETRAINER_HPP__
