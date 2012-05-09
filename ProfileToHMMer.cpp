/// #### TODO: ERE I AM.  Need to adapt to work with aminoacids.  For now this supports only DNA.

#include "Algebra.hpp"
#include "Profile.hpp"

#include <iostream>

#include <seqan/basic.h>

// For HMMer support
namespace hmmer {
extern "C" {
#include "hmmer/src/config.h"		/* compile-time configuration constants */
#include "hmmer/squid/squidconf.h"
#include "hmmer/src/structs.h"		/* data structures, macros, #define's   */
#include "hmmer/src/funcs.h"		/* function declarations                */
#include "hmmer/src/globals.h"		/* alphabet global variables            */
#include "hmmer/squid/squid.h"		/* general sequence analysis library    */
#include "hmmer/squid/msa.h"            /* squid's multiple alignment i/o       */
} // End extern "C"
} // End namespace hmmer

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
template <class ProfileType>
void
profileToHMMerProfile (
  //ProfileType const & profile,
  ProfileType & profile, // for balanced_indels and min_indels, which change the profile.
  struct hmmer::plan7_s * hmm,
  const double swentry = .5, // ignored if USE_DEL_IN_DEL_OUT
  const double swexit = .5, // ignored if USE_DEL_IN_DEL_OUT
  const double expected_distance_between_hits = 0, // 0 means use insertion length
  const double null_model_expected_length = 0, // 0 means use insertion length
  const bool be_verbose = false
)
{
  // If false, use the profile's insertion distribution for all insertions and for the null.  If true, use an even distribution for all of those.
  // TODO: ?
  static const bool use_even_insertion_distribution = false;//true;

  // TODO: REMOVE.  TESTING.
  static const bool balanced_indels = false;
  // TODO: REMOVE.  TESTING.
  static const bool min_indels = false;//true;
  // TODO: REMOVE.  TESTING.
  static const bool max_indels = false;//true;
  // TODO: REMOVE.  TESTING.
  static const double hack_indel_init_prob = 0;//.1;//.05; // 0 means use the one in the profile (or the min, max, or avg, depending on the other hack options)
  // TODO: REMOVE.  TESTING.
  static const double hack_expected_indel_length = 0;//2; // 0 means use the one in the profile (or the min, max, or avg, depending on the other hack options)

  // Make sure it's non-null.
  assert( hmm != 0 );
  // Make sure it's clean & new.
  //assert( hmm->dna2 == -( hmmer::INFTY ) );

  // For now, assuming Dna.  TODO: Generalize
  hmmer::SetAlphabet( hmmNUCLEIC );

  hmmer::AllocPlan7Body( hmm, profile.length() );
    
  if( be_verbose ) {
    cout << "Multiple local (hmmfs)" << endl;
#ifndef USE_DEL_IN_DEL_OUT
    // Lower entry/exit probs favor larger/longer hits.
    cout << "S/W aggregate entry probability:   " << swentry << endl;
    cout << "S/W aggregate exit probability:    " << swexit << endl;
#endif
    if( expected_distance_between_hits == 0 ) {
      cout << "Expected distance between hits:    " << ( 1.0 / profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] ) << endl;
    } else {
      cout << "Expected distance between hits:    " << expected_distance_between_hits << endl;
    }
    if( null_model_expected_length == 0 ) {
      cout << "Null model expected length:        " << ( 1.0 / profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] ) << endl;
    } else {
      cout << "Null model expected length:        " << null_model_expected_length << endl;
    }
    cout << "Converting galosh Profile to hmmer hmm:";
    cout.flush();
  }

  if( hack_indel_init_prob > 0 ) {
    profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
      hack_indel_init_prob;
    profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
    profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
      (
        1.0 -
        ( 2.0 * profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] )
      );
    cout << "Hacked indel init probabilities: M->I and M->D are " << profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] << "." << endl;
  } // End if hack_indel_init_prob
  if( hack_expected_indel_length > 0 ) {
    profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
      1.0 - ( 1.0 / hack_expected_indel_length );
    profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
    profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
      (
        1.0 -
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ]
      );
    profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];

    cout << "Hacked indel extension probabilities: I->I and D->D are " << profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] << "." << endl;
  } // End if hack_expected_indel_length

#ifdef USE_DEL_IN_DEL_OUT
      vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, doublerealspace> > scaled_match_distributions( profile.length() - 1 );
      for( uint32_t pos_i = 0; pos_i < ( profile.length() - 1 ); pos_i++ ) {
        profile.createScaledMatchDistributionForPosition(
          pos_i,
          scaled_match_distributions[ pos_i ]
        );
      } // End foreach pos_i
#endif // USE_DEL_IN_DEL_OUT

  // TODO: REMOVe.  TESTING.
  if( balanced_indels ) {
    // Use the average of the del-open and ins-open for both.
#ifdef USE_DEL_IN_DEL_OUT
    for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
      scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] +=
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ];
      scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] /=
        2.0;
      scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ] =
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ];

      cout << "[pos " << pos_i << "]: Using average indel probabilities: M->I and M->D are " << scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] << endl;
    } // End foreach pos_i
#else
    profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] +=
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ];
    profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] /=
      2.0;
    profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
#endif // USE_DEL_IN_DEL_OUT .. else ..

    // Use the average of del-ext and ins-ext for both.
    profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] +=
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
    profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] /=
      2.0;
    profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
      (
        1.0 -
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ]
      );
    profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
    profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
#ifndef USE_DEL_IN_DEL_OUT
    cout << "Using average indel probabilities: M->I and M->D are " << profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] << "; ";
#endif // !USE_DEL_IN_DEL_OUT
    cout << "I->I and D->D are " << profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] << "." << endl;
  } // End if balanced_indels

  // TODO: REMOVe.  TESTING.
  else if( min_indels ) {
    // Use the minimum of the del-open and ins-open for both, and use the -ext
    // values from that minimum-open type.
    if(
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] <
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
    ) {
#ifdef USE_DEL_IN_DEL_OUT
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ] =
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ];
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toMatch ] =
          (
            1.0 -
            ( 2.0 * scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] )
          );
        cout << "[pos " << pos_i << "]: Using minimum indel probabilities: M->I and M->D are " << scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] << endl;
      } // End foreach pos_i
#else
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
        profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
      profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
        (
          1.0 -
          ( 2.0 * profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] )
        );
      cout << "Using minimum indel probabilities: M->I and M->D are " << profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] << "; ";
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
    } else { // del-open <= ins-open
#ifdef USE_DEL_IN_DEL_OUT
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] =
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ];
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toMatch ] =
          (
            1.0 -
            ( 2.0 * scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] )
          );
        cout << "[pos " << pos_i << "]: Using minimum indel probabilities: M->I and M->D are " << scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] << endl;
      } // End foreach pos_i
#else
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
        profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ];
      profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
        (
          1.0 -
          ( 2.0 * profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] )
        );
      cout << "Using minimum indel probabilities: M->I and M->D are " << profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] << "; ";
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
        profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ];
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
        profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
    } // End if ins-open < del-open .. else ..
    cout << "I->I and D->D are " << profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] << "." << endl;
  } // End if min_indels

  // TODO: REMOVe.  TESTING.
  else if( max_indels ) {
    // Use the maximum of the del-open and ins-open for both, and use the -ext
    // values from that maximum-open type.
    if(
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] >
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
    ) {
#ifdef USE_DEL_IN_DEL_OUT
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ] =
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ];
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toMatch ] =
          (
            1.0 -
            ( 2.0 * scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] )
          );
        cout << "[pos " << pos_i << "]: Using maximum indel probabilities: M->I and M->D are " << scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] << endl;
      } // End foreach pos_i
#else
      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ] =
        profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
      profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
        (
          1.0 -
          ( 2.0 * profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] )
        );
      cout << "Using maximum indel probabilities: M->I and M->D are " << profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] << "; ";
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ] =
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
      profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ] =
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
    } else { // del-open <= ins-open
#ifdef USE_DEL_IN_DEL_OUT
      for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] =
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ];
        scaled_match_distributions[ pos_i ][ TransitionFromMatch::toMatch ] =
          (
            1.0 -
            ( 2.0 * scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] )
          );
        cout << "[pos " << pos_i << "]: Using maximum indel probabilities: M->I and M->D are " << scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ] << endl;
      } // End foreach pos_i
#else
      profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] =
        profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ];
      profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ] =
        (
          1.0 -
          ( 2.0 * profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] )
        );
      cout << "Using maximum indel probabilities: M->I and M->D are " << profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ] << "; ";
#endif // USE_DEL_IN_DEL_OUT .. else ..
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ] =
        profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ];
      profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] =
        profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
    } // End if ins-open < del-open .. else ..
    cout << "I->I and D->D are " << profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ] << "." << endl;
  } // End if max_indels
  // END HACKS

  for( uint32_t pos_i = 0; pos_i < profile.length(); pos_i++ ) {
    if( be_verbose ) {
      cout << '.';
      cout.flush();
    }
    for( uint32_t res_i = 0; res_i < seqan::ValueSize<Dna>::VALUE; res_i++ ) {
      hmm->mat[ pos_i + 1 ][ hmmer::SymbolIndex( ( char )Dna( res_i ) ) ] =
        toDouble(
          profile[ pos_i ][ Emission::Match ][ res_i ]
        );
      if( use_even_insertion_distribution ) {
        hmm->ins[ pos_i ][ hmmer::SymbolIndex( ( char )Dna( res_i ) ) ] =
          .25;
      } else {
        hmm->ins[ pos_i ][ hmmer::SymbolIndex( ( char )Dna( res_i ) ) ] =
          toDouble(
            profile[ Emission::Insertion ][ res_i ]
          );
      } // End if use_even_insertion_distribution
    } // End foreach res_i
    if( pos_i != ( hmm->M - 1 ) ) { // Last pos has no transitions
#ifdef USE_DEL_IN_DEL_OUT
      hmm->t[ pos_i + 1 ][ 0 ] =
        toDouble(
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toMatch ]
        );
      hmm->t[ pos_i + 1 ][ 1 ] =
        toDouble(
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toInsertion ]
        );
      hmm->t[ pos_i + 1 ][ 2 ] =
        toDouble(
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletion ]
        );
#else
      hmm->t[ pos_i + 1 ][ 0 ] =
        toDouble(
          profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ]
        );
      hmm->t[ pos_i + 1 ][ 1 ] =
        toDouble(
          profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ]
        );
      hmm->t[ pos_i + 1 ][ 2 ] =
        toDouble(
          profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
        );
#endif // USE_DEL_IN_DEL_OUT .. else ..
      hmm->t[ pos_i + 1 ][ 3 ] =
        toDouble(
          profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ]  
        );
      hmm->t[ pos_i + 1 ][ 4 ] =
        toDouble(
          profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ]
        );
      hmm->t[ pos_i + 1 ][ 5 ] =
        toDouble(
          profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ]
        );
      hmm->t[ pos_i + 1 ][ 6 ] =
        toDouble(
          profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ]
        );
    } // End if this is not the last position  
  } // End foreach pos_i

#ifdef USE_DEL_IN_DEL_OUT
  hmm->tbd1 =
    toDouble(
      profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ]
    );
#else
  // TODO: REMOVE? TESTING.
  // Use the M->D prob for B->D.
  hmm->tbd1 =
    (
      toDouble(
        profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
      ) /
      (
        1.0 -
        toDouble(
          profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ]
        )
      )
    );
  // TODO: REMOVE.
  //cout << "M->D: " << hmm->tbd1 << endl; 
  //cout << "\t = ( " << 
  //    toDouble(
  //      profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ]
  //    ) << " / ( 1.0 - " << toDouble(
  //        profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ]
  //    ) << " )" << endl;
#endif // USE_DEL_IN_DEL_OUT .. else ..
  
  // Plan7FSConfig sets the XT values using the p1 (null model extension prob)
  // for the pre- and post-align insertion extension probs, and for the loop
  // continuation prob -- it governs the amount of non-model stuff to expect in
  // between hits.  Note we'll set the p1 to a different value below (based on
  // the null_model_expected_length), because it also governs how many hits
  // hmmsearch returns.
  if( expected_distance_between_hits == 0 ) {
    hmm->p1 =
      toDouble(
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ]
      );
  } else {
    hmm->p1 = 
      1.0 - ( 1.0 / expected_distance_between_hits );
  }

  // Configure it
  //hmmer::Plan7LSConfig( hmm ); // "For repeat detection"
  //hmmer::Plan7GlobalConfig( hmm );
  hmmer::Plan7FSConfig( hmm, swentry, swexit );

  // Plan7FSConfig sets the XT values using the p1 (null model extension prob)
  // for the pre- and post-align insertion extension probs.

  // If we want, we can use our trained values of pre-align and post-align
  // insertion extension probs instead of the p1 value..
#ifndef DISALLOW_FLANKING_TRANSITIONS
  if( false ) {
    hmm->xt[ XTN ][ LOOP ] = toDouble( 
      profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ]
    );
    hmm->xt[ XTN ][ MOVE ] =
      toDouble( 
        profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ]
      );
  
    hmm->xt[ XTC ][ LOOP ] =
      toDouble(
        profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ]
      );
  
    hmm->xt[ XTC ][ MOVE ] =
      toDouble(
        profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ]
      );
  } // End if false
#endif // !DISALLOW_FLANKING_TRANSITIONS

#ifdef USE_DEL_IN_DEL_OUT
    /* Configure entry.  Use the Del-In values.
     */
    hmm->begin[1] =
      toDouble(
        profile[ Transition::fromBegin ][ TransitionFromBegin::toMatch ]
      );
    double cum_extend =
      toDouble(
        profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ]
      );
    for( uint32_t k = 2; k <= hmm->M; k++ ) {
      hmm->begin[ k ] =
        ( 
          cum_extend *
          toDouble(
            profile[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ]
          )
        );
      cum_extend *=
        toDouble(
          profile[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toDeletionIn ]
        );
    }

    /* Configure exit.  Use the Del-Out probs -- already calculated in the
     * scaled_match_distributions.
     */
    hmm->end[hmm->M] = 1.0;
    for( uint32_t pos_i = 0; pos_i < ( profile.length() - 1 ); pos_i++ ) {
      hmm->end[ pos_i + 1 ] =
        toDouble(
          scaled_match_distributions[ pos_i ][ TransitionFromMatch::toDeletionOut ]
        );
    }
    // Note that we don't do Plan7RenormalizeExits(hmm) since they are already
    // normalized (by createScaledMatchDistributionForPosition).
#else
  // Another option is to use a geometric del-in and del-out
  if( 1 ) {
    // params:
    // swentry is (B->Z * hmm->M), given !(B->D)
    // swentry_extend is Z->Z
    float swentry_extend = ( 1.0f - ( 1.0f / ( hmm->M / 4.0f ) ) ); // for now, so expected length of del-in is 1/4 of the profile...
    // swexit is (M->U * hmm->M)
    // swexit_extend is U->U
    float swexit_extend = ( 1.0f - ( 1.0f / ( hmm->M / 4.0f ) ) ); // for now, so expected length of del-out is 1/4 of the profile...

    int   k;			/* counter over states      */

    /* Configure entry.  Given no B->D or B->M, there's a geometric del-in
     * length, max M-1.
     */
    hmm->begin[1] = (1. - swentry) * (1. - hmm->tbd1);
    float cum_extend = ( 1.0 - hmm->tbd1 ) * ( swentry / hmm->M );
    for( k = 2; k < hmm->M; k++ ) {
      hmm->begin[ k ] = cum_extend * ( 1.0 - swentry_extend );
      cum_extend *= swentry_extend;
    }
    hmm->begin[ hmm->M ] = cum_extend; // it always starts at or before the end!

    /* Configure exit.  del-out lengths are also geometric, but note there's
     * never a path that leaves del-out once it enters..
     */
    hmm->end[hmm->M] = 1.0;
    cum_extend = ( swexit / hmm->M );
    for( k = ( hmm->M - 1 ); k >= 1; k-- ) {
      hmm->end[k] = cum_extend;
      cum_extend *= swexit_extend;
    }
    Plan7RenormalizeExits(hmm);
  } // End if true
#endif // USE_DEL_IN_DEL_OUT .. else ..

  // Null model
  float            p1;		/* null sequence model p1 transition       */
  float  null_model[ MAXABET ];	/* null sequence model                     */
  hmmer::P7DefaultNullModel( null_model, &p1 );
  if( use_even_insertion_distribution ) {
    // Set null model to the insertion distribution, which is even.
    null_model[ 0 ] = .25;
    null_model[ 1 ] = .25;
    null_model[ 2 ] = .25;
    null_model[ 3 ] = .25;
  } else {
    // Set null model to the insertion distribution.
    null_model[ 0 ] =
      toDouble(
        profile[ Emission::Insertion ][ seqan::Dna( 'a' ) ]
      );
    null_model[ 1 ] =
      toDouble(
        profile[ Emission::Insertion ][ seqan::Dna( 'c' ) ]
      );
    null_model[ 2 ] =
      toDouble(
        profile[ Emission::Insertion ][ seqan::Dna( 'g' ) ]
      );
    null_model[ 3 ] =
      toDouble(
        profile[ Emission::Insertion ][ seqan::Dna( 't' ) ]
      );
  }

  if( null_model_expected_length == 0 ) {
    // By default, use insertion length.
    p1 =
      toDouble(
        profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ]
      );
  } else {
    p1 = 1.0 - ( 1.0 / null_model_expected_length );
  }

  hmmer::Plan7SetNullModel( hmm, null_model, p1 );

  // NOTE: I've commented this out.  It is for converting pseudocount vectors
  // to a profile, but we have no pseudocounts...
  //// Priors
  ////struct hmmer::p7prior_s *pri;        /* Dirichlet priors to use */
  ////pri = hmmer::P7DefaultPrior();
  ////P7PriorifyHMM( hmm, pri );

  // Name it
  hmmer::Plan7SetName( hmm, "GaloshProfile" ); // TODO: MAGIC # (name)
  // Add the date
  hmmer::Plan7SetCtime( hmm );

  // Now renormalize it 
  hmmer::Plan7Renormalize( hmm );

  if( be_verbose ) {
    cout << ".done." << endl;
  }
} // copyFromHMMerProfile( ProfileType &, hmmer::hmm const *, SupportedTypes const & )
} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
  // For now we assume a Dna distribution.  TODO: Generalize.
  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input (galosh Profile) filename> [<output (HMMer profile) filename> [<sw aggregate entry prob> [<sw aggregate exit prob> [<expected distance between hits> [<null model expected length>]]]]" << endl;
    exit( 1 );
  }
  const bool be_verbose = true; //( argc >= 5 );
  if( be_verbose ) {
    cout << "Reading profile from file '" << argv[ 1 ] << "'" << endl;
  }
  galosh::ProfileTreeRoot<seqan::Dna, floatrealspace> profile;
  profile.fromFile( argv[ 1 ] );
  if( be_verbose ) {
    cout << "\tgot:" << std::endl;
    cout << profile;
    cout << endl;
  }
  if( profile.length() <= 1 ) {
    if( be_verbose ) {
      cout << "ERROR: Profile is too short!" << endl;
    }
    exit( 1 );
  }

  double swentry = .5;
  double swexit = .5;
  if( argc >= 4 ) {
    try {
      swentry = boost::lexical_cast<double>( argv[ 3 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 3 ] << "' as a real value for use as the aggregate S/W entry probability." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( swentry <= 0 ) {
      std::cerr << "The given aggregate S/W entry probability value, " << swentry << ", is zero or negative.  You must supply a value between 0 and 1, and not 0." << std::endl;
      exit( 1 );
    }
    if( swentry > 1 ) {
      std::cerr << "The given aggregate S/W entry probability value, " << swentry << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
      exit( 1 );
    }
    if( argc >= 5 ) {
      try {
        swexit = boost::lexical_cast<double>( argv[ 4 ] );
      } catch( boost::bad_lexical_cast & ) {
        std::cerr << "Unable to interpret the argument '" << argv[ 4 ] << "' as a real value for use as the aggregate S/W exit probability." << std::endl;
        exit( 1 );
      } // End try .. catch block for lexical_cast
#ifndef USE_DEL_IN_DEL_OUT
      if( swexit <= 0 ) {
        std::cerr << "The given aggregate S/W exit probability value, " << swexit << ", is zero or negative.  You must supply a value between 0 and 1, and not 0." << std::endl;
        exit( 1 );
      }
      if( swexit > 1 ) {
        std::cerr << "The given aggregate S/W exit probability value, " << swexit << ", is greater than 1.  You must supply a value between 0 and 1." << std::endl;
        exit( 1 );
      }
#endif // USE_DEL_IN_DEL_OUT
    } else {
      swexit = swentry;
    } // End if argc >= 5
  } // End if argc >= 4
  double expected_distance_between_hits = 0;
  if( argc >= 6 ) {
    try {
      expected_distance_between_hits = boost::lexical_cast<double>( argv[ 5 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 5 ] << "' as a real value for use as the expected distance between hits." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( expected_distance_between_hits < 0 ) {
      std::cerr << "The given expected distance between hits, " << expected_distance_between_hits << ", is negative.  You must supply a value greater than 0, or 0 to indicate that the insertion length should be used." << std::endl;
      exit( 1 );
    }
  } // End if argc >= 6
  double null_model_expected_length = 0;
  if( argc >= 7 ) {
    try {
      null_model_expected_length = boost::lexical_cast<double>( argv[ 6 ] );
    } catch( boost::bad_lexical_cast & ) {
      std::cerr << "Unable to interpret the argument '" << argv[ 6 ] << "' as a real value for use as the null model expected length." << std::endl;
      exit( 1 );
    } // End try .. catch block for lexical_cast
    if( null_model_expected_length < 0 ) {
      std::cerr << "The given null model expected length, " << null_model_expected_length << ", is negative.  You must supply a value greater than 0, or 0 to indicate that the insertion length should be used." << std::endl;
      exit( 1 );
    }
  } // End if argc >= 7

  struct hmmer::plan7_s *hmm = hmmer::AllocPlan7Shell();
  profileToHMMerProfile( profile, hmm, swentry, swexit, expected_distance_between_hits, null_model_expected_length, be_verbose );

  if( ( argc >= 3 ) && be_verbose ) {
    cout << "Writing HMMer profile to file '" << argv[ 2 ] << "'" << endl;
  } else if( be_verbose ) {
    cout << "Writing HMMer profile to stdout " << endl;
  }
  FILE *hmmfp;
  if( argc >= 3 ) {
    if( ( hmmfp = fopen( argv[ 2 ], "w" ) ) == NULL ) {
      cout << "Failed to open HMMer profile file '" << argv[ 2 ] << "' for writing." << endl;
    }
  } else {
    hmmfp = stdout;
  }
  WriteAscHMM( hmmfp, hmm );
  fflush( hmmfp );

  if( argc >= 3 ) {
    // Close the output file.
    fclose( hmmfp );
  }

  FreePlan7( hmm );
  exit( 0 );
} // main (..)
