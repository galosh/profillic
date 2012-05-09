#========================================
# USER RESERVED -- The following are reserved for users to set on the
# command line.  Makefiles should not set these.  These variables are
# for C/C++ compilation, and linking.
## -pg is for use with gprof.  Use on CFLAGS and on LDFLAGS
#CFLAGS		= -Wall -pg
#CFLAGS         = -Winline
# TAH 2/12 for debugging within Eclipse
CFLAGS		= -g3 -gdwarf-2  
#JFLAGS		=
#LDFLAGS	= -pg

# OPTIMIZE with the -O option.  Override from the command line for
# building debug versions.
#
#OPTFLAGS	= -O3 -funroll-loops -DNDEBUG=1

#========================================
### For Seqan support (required):
SEQAN_CFLAGS 	= -I./seqan
SEQAN_LDFLAGS 	=
SEQAN_LIBS	=

#========================================
## For the HMMoC-BFloat-Algebra library (required):
ALGEBRA_CFLAGS 	= -I./HMMoC-BFloat-Algebra
ALGEBRA_LDFLAGS = -L./HMMoC-BFloat-Algebra
ALGEBRA_LIBS	= -lHMMoC-BFloat-Algebra

#========================================
## For the prolific library (required):
PROLIFIC_CFLAGS 	= -I./prolific
PROLIFIC_LDFLAGS	=
PROLIFIC_LIBS		=

#========================================
## For Boost support (required):
BOOST_CFLAGS 	= -I./boost-include
BOOST_LDFLAGS 	= -L./boost-lib
# TAH 2/12 added program options
BOOST_LIBS	= -lboost_serialization -lboost_graph -lboost_filesystem -lboost_system -lboost_program_options

#========================================
### Paul added these for HMMer (and squid) support (optional)
HMMER_CFLAGS 	= -I./hmmer/src -I./hmmer/squid
HMMER_LDFLAGS 	= -L./hmmer/src -L./hmmer/squid
HMMER_LIBS	= -lsquid -lhmmer
#HMMER_CFLAGS 	=
#HMMER_LDFLAGS 	=
#HMMER_LIBS	=

###==============================================
#
# This section modified by TAH 2/12.  Note that we assume that the
# prolific library is below the current (profuse) directory.
#

# In case any of these hpp files change, rebuild programs
PROLIFIC_LIB = ./prolific
INCS = Profuse.hpp \
       $(PROLIFIC_LIB)/DynamicProgramming.hpp \
       $(PROLIFIC_LIB)/MultinomialDistribution.hpp \
       $(PROLIFIC_LIB)/Profile.hpp

ALIGN_INCS = $(INCS) \
ScoreAndMaybeAlign.hpp

PROF2ALPROF_INCS = $(INCS) \
GenAlignmentProfiles.hpp

SCORE_INCS = $(INCS) \
ScoreAndMaybeAlign.hpp

CREATERANDOMSEQUENCE_INCS = $(INCS)

DRAWSEQUENCES_INCS = $(INCS)

PROFILETOSEQUENCE_INCS = $(INCS)

SEQUENCETOPROFILE_INCS = $(INCS)

ALIGNEDFASTATOPROFILE_INCS = $(INCS)

PROFILETREETOPROFILE_INCS = $(INCS)

PROFILETOHMMER_INCS = $(INCS)

PROF2ALPROF_INCS = $(INCS) \
GenAlignmentProfiles.hpp

ALIGN_OBJS = Align.o

SCORE_OBJS = Score.o

CREATERANDOMSEQUENCE_OBJS = CreateRandomSequence.o

DRAWSEQUENCES_OBJS = DrawSequences.o

PROFILETOSEQUENCE_OBJS = ProfileToConsensus.o

SEQUENCETOPROFILE_OBJS = SequenceToProfile.o

ALIGNEDFASTATOPROFILE_OBJS = AlignedFastaToProfile.o

PROFILETREETOPROFILE_OBJS = ProfileTreeToProfile.o

PROFILETOHMMER_OBJS = ProfileToHMMer.o

PROF2ALPROF_OBJS = profileToAlignmentProfile.o

ALIGN_SOURCES = Align.cpp

SCORE_SOURCES = Score.cpp

CREATERANDOMSEQUENCE_SOURCES = CreateRandomSequence.cpp

DRAWSEQUENCES_SOURCES = DrawSequences.cpp

PROFILETOSEQUENCE_SOURCES = ProfileToConsensus.cpp

SEQUENCETOPROFILE_SOURCES = SequenceToProfile.cpp

ALIGNEDFASTATOPROFILE_SOURCES = AlignedFastaToProfile.cpp

PROFILETREETOPROFILE_SOURCES = ProfileTreeToProfile.cpp

PROFILETOHMMER_SOURCES = ProfileToHMMer.cpp

PROF2ALPROF_SOURCES = profileToAlignmentProfile.cpp

default: all

align: $(ALIGN_SOURCES) $(ALIGN_INCS) $(ALIGN_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o align $(ALIGN_OBJS) $(MUSCLE_CPPOBJ)

score: $(SCORE_SOURCES) $(SCORE_INCS) $(SCORE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o score $(SCORE_OBJS) $(MUSCLE_CPPOBJ)

createRandomSequence: $(CREATERANDOMSEQUENCE_SOURCES) $(CREATERANDOMSEQUENCE_INCS) $(CREATERANDOMSEQUENCE_OBJS)
	     $(CXX_LINK) -o createRandomSequence $(CREATERANDOMSEQUENCE_OBJS)

drawSequences: $(DRAWSEQUENCES_SOURCES) $(DRAWSEQUENCES_INCS) $(DRAWSEQUENCES_OBJS)
	     $(CXX_LINK) -o drawSequences $(DRAWSEQUENCES_OBJS)

profileToSequence: $(PROFILETOSEQUENCE_SOURCES) $(PROFILETOSEQUENCE_INCS) $(PROFILETOSEQUENCE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profileToSequence $(PROFILETOSEQUENCE_OBJS) $(MUSCLE_CPPOBJ)

sequenceToProfile: $(SEQUENCETOPROFILE_SOURCES) $(SEQUENCETOPROFILE_INCS) $(SEQUENCETOPROFILE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o sequenceToProfile $(SEQUENCETOPROFILE_OBJS) $(MUSCLE_CPPOBJ)

alignedFastaToProfile: $(ALIGNEDFASTATOPROFILE_SOURCES) $(ALIGNEDFASTATOPROFILE_INCS) $(ALIGNEDFASTATOPROFILE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o alignedFastaToProfile $(ALIGNEDFASTATOPROFILE_OBJS) $(MUSCLE_CPPOBJ)

profileTreeToProfile: $(PROFILETREETOPROFILE_SOURCES) $(PROFILETREETOPROFILE_INCS) $(PROFILETREETOPROFILE_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profileTreeToProfile $(PROFILETREETOPROFILE_OBJS)

profileToHMMer: $(PROFILETOHMMER_SOURCES) $(PROFILETOHMMER_INCS) $(PROFILETOHMMER_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) $(HMMER_LDFLAGS) $(HMMER_LIBS) -o profileToHMMer $(PROFILETOHMMER_OBJS)

profileToAlignmentProfile: $(PROF2ALPROF_SOURCES) $(PROF2ALPROF_INCS) $(PROF2ALPROF_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o profileToAlignmentProfile $(PROF2ALPROF_OBJS) $(MUSCLE_CPPOBJ)

converters: sequenceToProfile alignedFastaToProfile profileTreeToProfile profileToSequence profileToAlignmentProfile

progs: align score createRandomSequence drawSequences 

all: progs converters

## Recompile if the includes are modified ...
$(CREATERANDOMSEQUENCE_OBJS): $(CREATERANDOMSEQUENCE_SOURCES) $(CREATERANDOMSEQUENCE_INCS)
$(DRAWSEQUENCES_OBJS): $(DRAWSEQUENCES_SOURCES) $(DRAWSEQUENCES_INCS)
$(PROFILETOSEQUENCE_OBJS): $(PROFILETOSEQUENCE_SOURCES) $(PROFILETOSEQUENCE_INCS)
$(SEQUENCETOPROFILE_OBJS): $(SEQUENCETOPROFILE_SOURCES) $(SEQUENCETOPROFILE_INCS)
$(ALIGNEDFASTATOPROFILE_OBJS): $(ALIGNEDFASTATOPROFILE_SOURCES) $(ALIGNEDFASTATOPROFILE_INCS)
$(PROFILETREETOPROFILE_OBJS): $(PROFILETREETOPROFILE_SOURCES) $(PROFILETREETOPROFILE_INCS)
$(PROFILETOHMMER_OBJS): $(PROFILETOHMMER_SOURCES) $(PROFILETOHMMER_INCS)
$(ALIGN_OBJS): $(ALIGN_SOURCES) $(ALIGN_INCS)
$(SCORE_OBJS): $(SCORE_SOURCES) $(SCORE_INCS)

.PHONY: clean
clean:
	rm -f align score createRandomSequence drawSequences profileToSequence sequenceToProfile alignedFastaToProfile profileTreeToProfile profileToHMMer profileToAlignmentProfile $(ALIGN_OBJS) $(SCORE_OBJS) $(CREATERANDOMSEQUENCE_OBJS) $(DRAWSEQUENCES_OBJS) $(PROFILETOSEQUENCE_OBJS) $(SEQUENCETOPROFILE_OBJS) $(PROFILETREETOPROFILE_OBJS) $(ALIGNEDFASTATOPROFILE_OBJS) $(PROFILETOHMMER_OBJS) $(PROF2ALPROF_OBJS) 

#========================================
# FILE EXTENSIONS.  Extensions and prefixes for different types of
# files change from platform to platform.  Hide these in macros so
# that we can more easily cut and paste between makefiles.
o		= .o
EXE_SFX		= 
SCRIPT_SFX 	= 
LIB_PFX		= lib
LIB_SFX		= .a
LIB_SHARED_SFX	= .so
TMPLIB		= libtemp.a

# FILE TOOLS
AR 	= ar qv
CHMOD 	= chmod
CP	= cp
GREP	= grep
MKDIR 	= mkdir
MUNCH 	= stepmunch
MV	= mv
NM 	= nm
RANLIB	= ranlib
RM 	= rm -f
RMDIR 	= rm -rf
STRIP	= strip
UNZIP 	= unzip
ZIP 	= zip


#========================================
# ANSI C Compile and Link
#
CC		= gcc
CC_COMPILE	= $(CC) -c $(OPTFLAGS) $(CFLAGS) $(CC_CFLAGS) $(CC_SYSCFLAGS)
CC_LINK		= $(CC) $(LDFLAGS) $(CC_LDFLAGS) $(CC_SYSLDFLAGS) $(CC_LIBS)
CC_CFLAGS 	= $(ALGEBRA_CFLAGS) $(PROLIFIC_CFLAGS) $(BOOST_CFLAGS) $(SEQAN_CFLAGS) $(HMMER_CFLAGS)
CC_LDFLAGS	= $(ALGEBRA_LDFLAGS) $(PROLIFIC_LDFLAGS) $(BOOST_LDFLAGS) $(SEQAN_LDFLAGS)
CC_LIBS		= $(ALGEBRA_LIBS) $(PROLIFIC_LIBS) $(BOOST_LIBS) $(SEQAN_CLIBS)

# Global system things used for compilation, static linking, etc.
CC_SYSCFLAGS 	= -I.
CC_SYSLDFLAGS 	=
CC_SYSLIBS	=

#========================================
# C++ Compile and Link
#
CXX		= g++
CXX_COMPILE	= $(CXX) -c  $(OPTFLAGS) $(CFLAGS) $(CXX_CFLAGS) $(CXX_SYSCFLAGS)
CXX_LINK	= $(CXX) $(LDFLAGS) $(CXX_LDFLAGS) $(CXX_SYSLDFLAGS) $(CXX_LIBS)
CXX_CFLAGS 	= $(ALGEBRA_CFLAGS) $(PROLIFIC_CFLAGS) $(BOOST_CFLAGS) $(SEQAN_CFLAGS) $(HMMER_CFLAGS)
CXX_LDFLAGS	= $(ALGEBRA_LDFLAGS) $(PROLIFIC_LDFLAGS) $(BOOST_LDFLAGS) $(SEQAN_LDFLAGS)
CXX_LIBS	= $(ALGEBRA_LIBS) $(PROLIFIC_LDFLAGS) $(BOOST_LIBS) $(SEQAN_LIBS) 

# The force flags are used for C/C++ compilers that select the
# language based on the file naming conventions.  Some C++ source
# may be in files with C naming conventions.
CXX_FORCE	= 

# System Flags -- Things for static linking or making sure that the
# compiler understands that a file is a C++ file or whatever.  These
# usually change from platform to platform.
CXX_SYSCFLAGS 	= -I.
CXX_SYSLDFLAGS 	= 
CXX_SYSLIBS	= 

# Compilation Rules -- Repeat the rules for all of the different
# naming conventions.
#
.cxx.o:	; $(CXX_COMPILE) $<
.cpp.o:	; $(CXX_COMPILE) $<
.cc.o:	; $(CXX_COMPILE) $<
.C.o:	; $(CXX_COMPILE) $<

.cxx:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cpp:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cc:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.C:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)

# for legacy reasons also compile .c as c++
.c.o:	; $(CXX_COMPILE) $(CXX_FORCE) $<
.c:	
	$(CXX_COMPILE) $(CXX_FORCE) $<

