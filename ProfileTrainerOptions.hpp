/**
 * \file ProfileTrainerOptions.hpp
 * \author Ted Holzman 
 * \par Library:
 * Galosh profillic
 * \brief Options for ProfileTrainer
 * \copyright &copy; 2008, 2011, 2012, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/

///Useful for modifying GALOSH_DEF_OPT behavior.  See vector-valued options below
#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

/**
 * \note Note the lack of "protective" \#define symbols in this file.  We may
 * want to include it several times, redefining GALOSH_DEF_OPT as needed.
 */

      /**
       * Should we train the global parameters at all, or just do positions?
       * If true, train the globals.  If false, don't train the globals.
       */
GALOSH_DEF_OPT(trainProfileGlobals,bool,true,"Should we train the global parameters at all, or just do positions?" );

      /**
       * Should we train the position parameters at all?  If true, train the
       * position-specific parameters.  If false, don't train them.
       */
GALOSH_DEF_OPT(trainProfilePositions,bool,true,"Should we train the position-specific parameters at all?" );

      /**
       * Should we start with the global parameters, then do positions?  If
       * true, start with globals.  If false, start with position params,
       * then do globals.
       */
GALOSH_DEF_OPT(trainGlobalsFirst,bool,false,"Should we start with the global parameters, then do positions?  If true, start with globals.  If false, start with position params, then do globals." );

      /**
       * How many times at minimum should we iterate between position cycles and
       * global training steps before we end the training process ( one
       * iteration includes both positions and globals )?
       * Note that this must be at least 1.
       */
GALOSH_DEF_OPT(minIterations,uint32_t,1,"How many times at minimum should we iterate between position cycles and global training steps before we end the training process ( one iteration includes both positions and globals )? Note that this must be at least 1." );

      /**
       * How many times at maximum should we iterate between position cycles and
       * global training steps before we end the training process ( one
       * iteration includes both positions and globals )?
       * Note that this must be at least 1.
       */
GALOSH_DEF_OPT(maxIterations,uint32_t,1000,"How many times at maximum should we iterate between position cycles and global training steps before we end the training process ( one iteration includes both positions and globals )? Note that this must be at least 1." );

      /**
       * How many times should we cycle through the training steps before,
       * despite continuing improvement, we move on to refining indels?
       */
GALOSH_DEF_OPT(maxPositionCycles,uint32_t,1,"How many times should we cycle through the training steps before, despite continuing improvement, we move on to the next iteration?" );

      /**
       * How many times should we cycle through the training steps before,
       * despite continuing improvement, we move on to refining indels?
       *
       * This is for the sequence identifier params training step.
       */
GALOSH_DEF_OPT(maxPositionCycles_sequence_identifiers,uint32_t,4,"How many times should we cycle through the training steps before, despite continuing improvement, we move on to refining indels?");

      /**
       * How many times should we cycle through the training steps before,
       * despite continuing improvement, we move on to refining indels?
       *
       * This is for the global params training step.
       */
GALOSH_DEF_OPT(maxPositionCycles_globals,uint32_t,4,"How many times should we cycle through the training steps before, despite continuing improvement, we move on the next iteration? This is for the global params training step." );

      /**
       * Set this to 0 for straight-up BW/CBW.  Non-zero values allow
       * intermediate step sizes in the direction of the BW/CBW update.
       * Sometimes we have to reject the conditional BW updates because the
       * score would go down.  In such cases we can try again at various
       * intermediate steps.  We do this by adding some fraction of the updated
       * position values to the original values.  When this fraction is 1, for
       * instance, the old and new values will be averaged.  When the fraction
       * is not 1, it will be a weighted average.  The value of this
       * parameter is the largest denomenator of the weight on the new values.
       * We try denomenators from
       * minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by
       * baumWelchInverseScalarIncrement.
       */
GALOSH_DEF_OPT(maxBaumWelchInverseScalar,float,0.0f,"Set this to 0 for straight-up BW/CBW.  Non-zero values allow intermediate step sizes in the direction of the BW/CBW update.\nDetails: Sometimes we have to reject the conditional BW updates because the score would go down.  In such cases we can try again at various intermediate steps.  We do this by adding some fraction of the updated position values to the original values.  When this fraction is 1, for instance, the old and new values will be averaged.  When the fraction is not 1, it will be a weighted average.  The value of this parameter is the largest denomenator of the weight on the new values. We try denomenators from minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by baumWelchInverseScalarIncrement." );

      /**
       * See maxBaumWelchInverseScalar.
       */
GALOSH_DEF_OPT(minBaumWelchInverseScalar,float,0.0F,"See maxBaumWelchInverseScalar.  If minBaumWelchInverseScalar is greater than 0, we don't even try the usual direct update." );

      /**
       * See maxBaumWelchInverseScalar.
       */
GALOSH_DEF_OPT(baumWelchInverseScalarIncrement,float,4.0f,"See maxBaumWelchInverseScalar.  This should be between minBaumWelchInverseScalar and maxBaumWelchInverseScalar." );

      /**
       * Sometimes we have to reject the conditional BW updates because the
       * score would go down.  In such cases we try again at various
       * intermediate steps.  We do this by adding some fraction of the updated
       * position values to the original values.  When this fraction is 1, for
       * instance, the old and new values will be averaged.  When the fraction
       * is less than 1, it will be a weighted average.  The value of this
       * parameter is the largest denominator of the weight on the new values.
       * We try denominators from
       * minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by
       * baumWelchInverseScalarIncrement.  If minBaumWelchInverseScalar is
       * greater than 0, we don't even try the usual direct update.
       *
       * This is for the sequence identifier params training step.
       */
GALOSH_DEF_OPT(maxBaumWelchInverseScalar_sequence_identifiers,float,40.0F,"This is the max sequence identifier for the params training step");

      /**
       * See maxBaumWelchInverseScalar_sequence_identifiers.
       *
       * This is for the sequence identifier params training step.
       */
GALOSH_DEF_OPT(minBaumWelchInverseScalar_sequence_identifiers,float,0.0F,"This is the min sequence identifier for the params training step");

      /**
       * See maxBaumWelchInverseScalar_sequence_identifiers.
       *
       * This is for the sequence identifier params training step.
       */
GALOSH_DEF_OPT(baumWelchInverseScalarIncrement_sequence_identifiers,float,10.0F,"This is increment for the for the params training step");

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
GALOSH_DEF_OPT(maxBaumWelchInverseScalar_globals,float,100.0,"See maxBaumWelchInverseScalar.  This is for the global params training step." );

      /**
       * See maxBaumWelchInverseScalar_globals.
       *
       * This is for the global params training step.
       */
GALOSH_DEF_OPT(minBaumWelchInverseScalar_globals,float,0.0F,"See maxBaumWelchInverseScalar_globals.  If minBaumWelchInverseScalar_globals is greater than 0, we don't even try the usual direct update." );

      /**
       * See maxBaumWelchInverseScalar_globals.
       *
       * This is for the global params training step.
       */
GALOSH_DEF_OPT(baumWelchInverseScalarIncrement_globals,float,20.0F,"See maxBaumWelchInverseScalar_globals.  This should be between minBaumWelchInverseScalar_globals and maxBaumWelchInverseScalar_globals." );

      /**
       * How much (percentage) change in forward score after an iteration is
       * close enough to 0 that we should stop training?
       */
GALOSH_DEF_OPT(scorePercentChangeMinimum_iteration,float,1.0F,"How much (percentage) change in forward score after an iteration is close enough to 0 that we should stop training?" );

      /**
       * How much (percentage) change in forward score in a position cycle is
       * close enough to 0 that we should move on to refining indels?
       */
GALOSH_DEF_OPT(scorePercentChangeMinimum_position_cycle,float,1.0F,"How much (percentage) change in forward score in a position cycle is close enough to 0 that we should move on to refining indels?" );

      /**
       * What average squared euclidean distance (averaged over free
       * parameters) is close enough to 0 that we should stop training?
       */
GALOSH_DEF_OPT(euclideanDistanceMinimum_iteration,float,1E-5F,"What average squared euclidean distance (averaged over free parameters) is close enough to 0 that we should stop training?" );

      /**
       * What average squared euclidean distance (averaged over free
       * parameters) is close enough to 0 that we should move on to refining
       * indels?
       */
GALOSH_DEF_OPT(euclideanDistanceMinimum_position_cycle,float,1E-5F,"What average squared euclidean distance (averaged over free parameters) is close enough to 0 that we should move on to refining indels?" );

      /**
       * Always accept the changes to the profile, even if the score goes down?
       * Theoretically, the BW update should always lead to an improvement, but
       * because of numerical instability (or strong priors), sometimes it
       * won't.
       *
       * See alwaysAccept_disallowThreshold_profileDistance_iteration.
       */
GALOSH_DEF_OPT(alwaysAccept,bool,true,"Always accept the changes to the profile, even if the score goes down? Theoretically, the BW update should always lead to an improvement, but because of numerical instability (or strong priors), sometimes it won't.  Note that we recommend that you set alwaysAccept to false if usePriors=true and you are using lengthadjust." );

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
GALOSH_DEF_OPT(useAlignmentProfiles,bool,true,"Calculate and use Alignment Profiles for the positions BW update");

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
GALOSH_DEF_OPT(proposeProfileLengthChanges,bool,false,"Resize the profile while training it");

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
GALOSH_DEF_OPT(proposeDeletingThreshold,double,0.1,"If the alignment profile (for a particular sequence) at a position position has higher that this threshold of deletions, then the sequence will be counted towards the total number of sequences with 'deletes' at that position -- and if the fraction of sequences exceeds one minus the proposeInsertingOccupancyThreshold, the position will be deleted (if proposeProfileLengthChanges is true)");

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
GALOSH_DEF_OPT(proposeInsertingThreshold,double,0.1,"If the alignment profile (for a particular sequence) at an internal position has higher that this threshold of insertions, then the sequence will be counted towards the total number of sequences with 'inserts' at that position -- and if the fraction of sequences exceeds the proposeInsertingOccupancyThreshold, a new position will be inserted (if proposeProfileLengthChanges is true)");

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
GALOSH_DEF_OPT(proposeDeletingThreshold_increment,double,0.025,"If we detect a cycle, we increment the proposeDeletingThreshold by this amount.  We also increment it by this amount if it gets to the minimum or maximum profile length (the minimum is 2, the maximum is the longest sequence being used to train the profile).  See also proposeInsertingThreshold_increment.  If useSensitiveThresholding is true, thresholds will instead be incremented to the least amount that would have made a difference in which positions were inserted or deleted -- and this value will be used added as a little extra boost");

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
GALOSH_DEF_OPT(proposeInsertingThreshold_increment,double,0.025,"If we detect a cycle, we increment the proposeInsertingThresholds by this amount.  We also increment them by this amount if the length gets to the minimum or maximum profile length (the minimum is 2, the maximum is the longest sequence being used to train the profile).  See also proposeDeletingThreshold_increment.  If useSensitiveThresholding is true, thresholds will instead be incremented to the least amount that would have made a difference in which positions were inserted or deleted -- and this value will be used added as a little extra boost");

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
GALOSH_DEF_OPT(proposeInsertingPreAlignThreshold,double,0.1,"If the alignment profile (for a particular sequence) for the first/pre-align position has higher than this threshold of pre-align insertions, then the sequence will be counted towards the total number of sequences with 'inserts' at that position -- and if the fraction of sequences exceeds the proposeInsertingOccupancyThreshold, a new position will be inserted (if proposeProfileLengthChanges is true)");

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
GALOSH_DEF_OPT(proposeInsertingPostAlignThreshold,double,0.1,"If the alignment profile (for a particular sequence) at the last position has higher than this threshold of post-align insertions, then the sequence will be counted towards the total number of sequences with 'inserts' at that position -- and if the fraction of sequences exceeds the proposeInsertingOccupancyThreshold, a new position will be inserted (if proposeProfileLengthChanges is true)");

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
GALOSH_DEF_OPT(proposeInsertingOccupancyThreshold,double,0.5,"If the fraction of sequences for which the insertion fraction exceeds the relevant threshold is greater than proposeInsertingOccupancyThreshold (for a particular position), then a new position will be inserted.  If the fraction of sequences for which the deletion fraction exceeds the proposeDeletingThreshold is greater than ( 1.0 - proposeInsertingOccupancyThreshold ), then the position will be deleted (if proposeProfileLengthChanges is true)");

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
GALOSH_DEF_OPT(useSensitiveThresholding,bool,true,"If true, thresholds will not be incremented until the end of a cycle, and then will be incremented to the least value that would have made a difference in the positions being inserted or deleted");

      /**
       * If things are taking too long, it might be because there's an
       * undetected cycle.  At
       * increaseThresholdsForLengthChanges_startIteration, we start treating
       * all length changes like they're part of a cycle, to force convergence.
       *
       * See increaseThresholdsForLengthChanges_minIncrement.
       */
GALOSH_DEF_OPT(increaseThresholdsForLengthChanges_startIteration,uint32_t,500,"If things are taking too long, it might be because there's an undetected cycle.  At increaseThresholdsForLengthChanges_startIteration, we start treating all length changes like they're part of a cycle, to force convergence");

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
GALOSH_DEF_OPT(increaseThresholdsForLengthChanges_minIncrement,float,1E-4F,"If things are taking too long, it might be because there's an undetected cycle.  At increaseThresholdsForLengthChanges_startIteration, we start treating all length changes like they're part of a cycle, to force convergence. If the threshold increments are very small, this might not help much, so we also enforce a minimum increment");

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
GALOSH_DEF_OPT(alwaysAccept_disallowThreshold_profileDistance_iteration,float,1E-5F,"When the profile distance (in an iteration) is less than alwaysAccept_disallowThreshold_profileDistance_iteration, we turn alwaysAccept off, so the score never goes down and the algorithm converges");

      /**
       * if proposeProfileLengthChanges is true, this number of iterations must
       * pass between proposals.  It may be 0.  Regardless of its value, at
       * least 1 iteration must pass at the very beginning of training before
       * length changes are allowed.
       */
GALOSH_DEF_OPT(numIterationsBetweenLengthChanges,uint32_t,0,"If proposeProfileLengthChanges is true, this number of iterations must pass between proposals");

      /**
       * Use the usual kind of Baum-Welch, in which parameter updates are
       * simultaneous?
       * If true, updates all parameters at once via Baum-Welch.
       * Otherwise, performs a conditional optimization.
       */
GALOSH_DEF_OPT(useUnconditionalBaumWelch,bool,false,"Use the usual kind of Baum-Welch");

      /**
       * Use Baldi's style of gradient ascent (GEM vs EM)?  If so, with what
       * learning rate?
       *
       * If non-zero, Baldi's "Smooth on-line learning algorithm" will be used,
       * with the given learning rate.
       *
       * @see baldiTemperature
       */
GALOSH_DEF_OPT(baldiLearningRate,double,0.0,"Use Baldi's style of gradient ascent (GEM vs EM)?  If so, with what learning rate?");

      /**
       * If we are using Baldi's style of gradient ascent (GEM vs EM), use this
       * temperature.
       *
       * @see baldiLearningRate
       */
GALOSH_DEF_OPT(baldiTemperature,double,1.0,"If we are using Baldi's style of gradient ascent (GEM vs EM), use this temperature");

      /**
       * Minimum profile distribution value.  Applies to everything except the
       * E (End) state.  Given as a double, despite the Profile ProbabilityType.
       */
GALOSH_DEF_OPT(profileValueMinimum,double,1E-5,"Minimum profile distribution value");

      /**
       * Use priors?  If true, m_matchEmissionPrior and m_globalsPrior will be
       * defined and used.  See also matchEmissionPrior and globalsPrior.
       */
GALOSH_DEF_OPT(usePriors,bool,true,"Use priors?  If true, m_matchEmissionPrior and m_globalsPrior will be defined and used.  Note that we recommend that you set alwaysAccept to false if using priors with lengthadjust." );
      /**
       * Set this to 0 for straight-up Baldi.  If non-zero, use Baldi / Siegel and 
       * quit trying to find a peak after this many attempts.  Note this only applies 
       * when baldiLearningRate > 0 and when ALLOW_BOLTZMANN_GIBBS is #defined.
       */
GALOSH_DEF_OPT(siegelMaxFindingThePeakAttempts_positions,uint32_t,1000,"Set this to 0 for straight-up Baldi.  If non-zero, use Baldi / Siegel and quit trying to find a peak after this many attempts");

      /**
       * Set to true to try both BW / CBW *and* Baldi/Siegel, then take whichever move is best.
       */
GALOSH_DEF_OPT(baldiHybrid,bool,false,"Set to true to try both BW / CBW *and* Baldi/Siegel, then take whichever move is best");

      /**
       * Length of steps for Baldi / Siegel (this is multiplied into the learning rate from eqn 3
       * from the Baldi paper to explore trapping a peak).  Note steps down will use the inverse of this.  
       * This must be strictly > 1.
       */
GALOSH_DEF_OPT(siegelEpsilonScaleFactor,double,2.0,"Length of steps for Baldi / Siegel (this is multiplied into the learning rate from eqn 3 from the Baldi paper to explore trapping a peak).  Note steps down will use the inverse of this.  This must be strictly > 1.");

      /**
       * If doing Unconditional BW or Unconditional Baldi, should we still isolate updates to the
       * globals from updates to the position params?  This is a hybrid approach, conditional but 
       * not per-position.
       */
GALOSH_DEF_OPT(unconditionalIsolatesGlobals,bool,false,"If doing Unconditional BW or Unconditional Baldi, should we still isolate updates to the globals from updates to the position params?  This is a hybrid approach, conditional but not per-position.");

      /**
       * This should generally be set to a fairly high value; it determines the number of attempts 
       * to maximize the quadratic after trapping a peak in the Baldi / Siegel Newton-Raphson part.  
       * The more relevant parameter to adjust is siegelRefiningThePeakStepsConvergenceThreshold.
       */
GALOSH_DEF_OPT(siegelMaxRefiningThePeakSteps_positions,uint32_t,1000,"This should generally be set to a fairly high value; it determines the number of attempts to maximize the quadratic after trapping a peak in the Baldi / Siegel Newton-Raphson part.  The more relevant parameter to adjust is siegelRefiningThePeakStepsConvergenceThreshold.");

GALOSH_DEF_OPT(siegelRefiningThePeakStepsConvergenceThreshold,double,1E-5,"When (the larger minus the smaller, divided by the smaller) of the last-best and the current-best score is less than this threshold, declare the Baldi / Siegel hillclimb done.");

      /**
       * If the score doesn't change after this number of attempts to multiply the gradient by an epsilon,
       * give up.  Note that this must be > 0. Note this only applies when baldiLearningRate > 0, when 
       * siegelMaxFindingThePeakAttempts_positions > 0, and when ALLOW_BOLTZMANN_GIBBS is #defined.
       */
GALOSH_DEF_OPT(siegelMaxFindingTheGradientAttempts_positions,uint32_t,10,"If the score doesn't change after this number of attempts to multiply the gradient by an epsilon, give up.");

      /**
       * Give up trapping the peak if the epsilon value we're multiplying by the gradient gets below
       * this value (we multiply it by siegelEpsilonScaleFactor or 1/SiegelEpsilonScaleFactor as we go..).  
       * Note that this must be > 0. Note this only applies when baldiLearningRate > 0, when 
       * siegelMaxFindingThePeakAttempts_positions > 0, and when ALLOW_BOLTZMANN_GIBBS is #defined.
       */
GALOSH_DEF_OPT(siegelMinEpsilon,double,1E-5,"Give up trapping the peak if the epsilon value we're multiplying by the gradient gets below this value (we multiply it by siegelEpsilonScaleFactor or 1/SiegelEpsilonScaleFactor as we go..).  Note that this must be > 0.");

      /**
       * When (the larger minus the smaller, divided by the smaller) of the last-best and the current-best 
       * score is less than this threshold, declare the Baldi / Siegel hillclimb done.
       */

/// This is how we're handling vectors.  It is a work-around because vectors are handled specially
/// by boost::program_options.  It allows the command line to look something like
///
/// --profileLengths 10 20 30
///
/// The TMP_EXTRA_STUFF must be set to include (at least) the ->multitoken() thing.
/// It should also be unset at the bottom of the vector initializations.


#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF ->multitoken()
///Vector definitions go here.

/** do this after the vector definition section */
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
