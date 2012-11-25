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

GALOSH_DEF_OPT(baldiHybrid,bool,false,"Set to true to try both BW / CBW *and* Baldi/Siegel, then take whichever move is best.");

GALOSH_DEF_OPT(siegelEpsilonScaleFactor,double,1.5,"Length of steps for Baldi / Siegel (this is multiplied into the learning rate from eqn 3 from the Baldi paper to explore trapping a peak).  Note steps down will use the inverse of this.  This must be strictly > 1.");

GALOSH_DEF_OPT(siegelMaxFindingThePeakAttempts_positions,uint32_t,1000,"Set this to 0 for straight-up Baldi.  If non-zero, use Baldi / Siegel and quit trying to find a peak after this many attempts.");
GALOSH_DEF_OPT(siegelMaxRefiningThePeakSteps_positions,uint32_t,1000,"This should generally be set to a fairly high value; it determines the number of attempts to maximize the quadratic after trapping a peak in the Baldi / Siegel Newton-Raphson part.  The more relevant parameter to adjust is siegelRefiningThePeakStepsConvergenceThreshold.");
GALOSH_DEF_OPT(siegelRefiningThePeakStepsConvergenceThreshold,double,1E-5,"When (the larger minus the smaller, divided by the smaller) of the last-best and the current-best score is less than this threshold, declare the Baldi / Siegel hillclimb done.");

GALOSH_DEF_OPT(siegelMaxFindingTheGradientAttempts_positions,uint32_t,10,"If the score doesn't change after this number of attempts to multiply the gradient by an epsilon, give up.");

GALOSH_DEF_OPT(siegelMinEpsilon,double,1E-5,"Give up trapping the peak if the epsilon value we're multiplying by the gradient gets below this value (we multiply it by siegelEpsilonScaleFactor or 1/SiegelEpsilonScaleFactor as we go..).  Note that this must be > 0.");



GALOSH_DEF_OPT(unconditionalIsolatesGlobals,bool,false,"If doing Unconditional BW or Unconditional Baldi, should we still isolate updates to the globals from updates to the position params?  This is a hybrid approach, conditional but not per-position.");



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
