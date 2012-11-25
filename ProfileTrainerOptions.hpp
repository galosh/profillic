/**
 * \file ProfileTrainerOptions.hpp
 * \author Paul T Edlefsen (pedlefse@scharp.org)
 * \par Library:
 * Galosh profillic
 * \brief Options for ProfileTrainer
 * \copyright &copy; 2008, 2011, 2012 by Paul T. Edlefsen, Fred Hutchinson Cancer
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
 * The configuration file contains the non-default, non-command-line options
 *
 */
//GALOSH_DEF_OPT(configFile,string,"Train.cfg","File path for the configuration file");

/**
 * The RNG seed..
 */
//GALOSH_DEF_OPT(seed,uint32_t,0,"The seed to use, or 0 to use the current time.");


GALOSH_DEF_OPT(unconditionalIsolatesGlobals,bool,false,"If doing Unconditional BW or Unconditional Baldi, should we still isolate updates to the globals from updates to the position params?  This is a hybrid approach, conditional but not per-position.");

GALOSH_DEF_OPT(trainProfileGlobals,bool,true,"Should we train the global parameters at all, or just do positions?" );
GALOSH_DEF_OPT(trainProfilePositions,bool,true,"Should we train the position-specific parameters at all?" );
GALOSH_DEF_OPT(trainGlobalsFirst,bool,false,"Should we start with the global parameters, then do positions?  If true, start with globals.  If false, start with position params, then do globals." );

GALOSH_DEF_OPT(minIterations,uint32_t,1,"How many times at minimum should we iterate between position cycles and global training steps before we end the training process ( one iteration includes both positions and globals )? Note that this must be at least 1." );
GALOSH_DEF_OPT(maxIterations,uint32_t,1000,"How many times at maximum should we iterate between position cycles and global training steps before we end the training process ( one iteration includes both positions and globals )? Note that this must be at least 1." );
GALOSH_DEF_OPT(maxPositionCycles,uint32_t,1,"How many times should we cycle through the training steps before, despite continuing improvement, we move on to the next iteration?" );
GALOSH_DEF_OPT(maxPositionCycles_globals,uint32_t,1,"How many times should we cycle through the training steps before, despite continuing improvement, we move on the next iteration? This is for the global params training step." );

GALOSH_DEF_OPT(maxBaumWelchInverseScalar,float,8.0,"Set this to 0 for straight-up BW/CBW.  Non-zero values allow intermediate step sizes in the direction of the BW/CBW update.\nDetails: Sometimes we have to reject the conditional BW updates because the score would go down.  In such cases we can try again at various intermediate steps.  We do this by adding some fraction of the updated position values to the original values.  When this fraction is 1, for instance, the old and new values will be averaged.  When the fraction is not 1, it will be a weighted average.  The value of this parameter is the largest denomenator of the weight on the new values. We try denomenators from minBaumWelchInverseScalar..maxBaumWelchInverseScalar, by baumWelchInverseScalarIncrement." );
GALOSH_DEF_OPT(minBaumWelchInverseScalar,float,0,"See maxBaumWelchInverseScalar.  If minBaumWelchInverseScalar is greater than 0, we don't even try the usual direct update." );
GALOSH_DEF_OPT(baumWelchInverseScalarIncrement,float,4.0,"See maxBaumWelchInverseScalar.  This should be between minBaumWelchInverseScalar and maxBaumWelchInverseScalar." );

GALOSH_DEF_OPT(maxBaumWelchInverseScalar_globals,float,100.0,"See maxBaumWelchInverseScalar.  This is for the global params training step." );
GALOSH_DEF_OPT(minBaumWelchInverseScalar_globals,float,0,"See maxBaumWelchInverseScalar_globals.  If minBaumWelchInverseScalar_globals is greater than 0, we don't even try the usual direct update." );
GALOSH_DEF_OPT(baumWelchInverseScalarIncrement_globals,float,20.0,"See maxBaumWelchInverseScalar_globals.  This should be between minBaumWelchInverseScalar_globals and maxBaumWelchInverseScalar_globals." );

GALOSH_DEF_OPT(scorePercentChangeMinimum_iteration,float,1.0,"How much (percentage) change in forward score after an iteration is close enough to 0 that we should stop training?" );

GALOSH_DEF_OPT(scorePercentChangeMinimum_position_cycle,float,1.0,"How much (percentage) change in forward score in a position cycle is close enough to 0 that we should move on to refining indels?" );

GALOSH_DEF_OPT(usePriors,bool,true,"Use priors?  If true, m_matchEmissionPrior and m_globalsPrior will be defined and used.  Note that we recommend that you set alwaysAccept to false if using priors with lengthadjust." );
GALOSH_DEF_OPT(alwaysAccept,bool,true,"Always accept the changes to the profile, even if the score goes down? Theoretically, the BW update should always lead to an improvement, but because of numerical instability (or strong priors), sometimes it won't.  Note that we recommend that you set alwaysAccept to false if usePriors=true and you are using lengthadjust." );


GALOSH_DEF_OPT(euclideanDistanceMinimum_iteration,float,1E-5,"What average squared euclidean distance (averaged over free parameters) is close enough to 0 that we should stop training?" );

GALOSH_DEF_OPT(euclideanDistanceMinimum_position_cycle,float,1E-5,"What average squared euclidean distance (averaged over free parameters) is close enough to 0 that we should move on to refining indels?" );


GALOSH_DEF_OPT(baldiHybrid,bool,false,"Set to true to try both BW / CBW *and* Baldi/Siegel, then take whichever move is best.");

GALOSH_DEF_OPT(siegelEpsilonScaleFactor,double,1.5,"Length of steps for Baldi / Siegel (this is multiplied into the learning rate from eqn 3 from the Baldi paper to explore trapping a peak).  Note steps down will use the inverse of this.  This must be strictly > 1.");

GALOSH_DEF_OPT(siegelMaxFindingThePeakAttempts_positions,uint32_t,1000,"Set this to 0 for straight-up Baldi.  If non-zero, use Baldi / Siegel and quit trying to find a peak after this many attempts.");
GALOSH_DEF_OPT(siegelMaxRefiningThePeakSteps_positions,uint32_t,1000,"This should generally be set to a fairly high value; it determines the number of attempts to maximize the quadratic after trapping a peak in the Baldi / Siegel Newton-Raphson part.  The more relevant parameter to adjust is siegelRefiningThePeakStepsConvergenceThreshold.");
GALOSH_DEF_OPT(siegelRefiningThePeakStepsConvergenceThreshold,double,1E-5,"When (the larger minus the smaller, divided by the smaller) of the last-best and the current-best score is less than this threshold, declare the Baldi / Siegel hillclimb done.");

GALOSH_DEF_OPT(siegelMaxFindingTheGradientAttempts_positions,uint32_t,10,"If the score doesn't change after this number of attempts to multiply the gradient by an epsilon, give up.");

GALOSH_DEF_OPT(siegelMinEpsilon,double,1E-5,"Give up trapping the peak if the epsilon value we're multiplying by the gradient gets below this value (we multiply it by siegelEpsilonScaleFactor or 1/SiegelEpsilonScaleFactor as we go..).  Note that this must be > 0.");




/// This is how we're handling vectors.  It is a work-around because vectors are handled specially
/// by boost::program_options.  It allows the command line to look something like
///
/// --profileLengths 10 20 30
///
/// The TMP_EXTRA_STUFF must be set to include (at least) the ->multitoken() thing.
/// It should also be unset at the bottom of the vector initializations.

#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF ->multitoken()

//GALOSH_DEF_OPT(profileLengths,myVector<int>,myVector<int>(1,100) BOOST_PP_COMMA() string("100"),"Lengths of the profiles/fasta seqs");

/** do this after the vector definition section */
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
