# Makefile for the ROOT test programs.  # This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

#------------------------------------------------------------------------------
#brings in the few minossoft things i need
#BOOSTFLAGS = -I boost_1_48_0

# commented out for kingbee and older versions of gcc
#CPPSTD = c++11

#CC = gcc-4.8
#CXX = g++-4.8
#LD = $(CC)

ifeq ($(CPPSTD), c++11)
CPPSTD_FLAGS = -std=c++11

#ifeq ($(CXX), clang++)
#CPPSTD_FLAGS += -stdlib=libc++
#endif

else
# If not compiling with C++11 support, all occurrences of "constexpr"
# must be replaced with "const", because "constexpr" is a keyword
# which pre-C++11 compilers do not support.
# ("constexpr" is needed in the code to perform in-class initialization
# of static non-integral member objects, i.e.:
#		static const double c_light = 2.99e8;
# which works in C++03 compilers, must be modified to:
#		static constexpr double c_light = 2.99e8;
# to work in C++11, but adding "constexpr" breaks C++03 compatibility.
# The following compiler flag defines a preprocessor macro which is
# simply:
#		#define constexpr const
# which replaces all instances of the text "constexpr" and replaces it
# with "const".
# This preserves functionality while only affecting very specific semantics.
CPPSTD_FLAGS = -Dconstexpr=const
endif


GENERAL_FLAGS = -pipe
OPTIMIZE_FLAGS = -O2
DEBUG_FLAGS = -g -ggdb
PROFILING_FLAGS =
ARCHITECTURE_FLAGS = -m64 -pthread
WARN_FLAGS = -W -Wall -Wextra -Woverloaded-virtual# -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable

CXXFLAGS += $(GENERAL_FLAGS) $(CPPSTD_FLAGS) $(ARCHITECTURE_FLAGS) $(OPTIMIZE_FLAGS) $(WARN_FLAGS) $(ROOTCFLAGS) $(INC_ANITA_UTIL)

#CXXFLAGS += $(CPPSTD_FLAGS) -g -O2 $(INC_ANITA_UTIL) $(BOOSTFLAGS) $(WARN_FLAGS)

OPTFLAGS  = -O2
DBGFLAGS  = -pipe -Wall -W -Woverloaded-virtual -g -ggdb -O0 -fno-inline

DBGCXXFLAGS = $(DBGFLAGS) $(ROOTCFLAGS) $(BOOSTFLAGS)

LDFLAGS  += $(CPPSTD_FLAGS) -g $(LD_ANITA_UTIL) -I$(BOOST_ROOT) -L.
#LDFLAGS  += $(CPPSTD_FLAGS) -g $(LD_ANITA_UTIL) -I$(BOOST_ROOT) $(ROOTLDFLAGS) -L. 


LIBS += 

HEADERS	  = rx.hpp Taumodel.hh
##ANITA_DATA_HEADERS = include/RawAnitaEvent.h include/UsefulAnitaEvent.h include/RawAnitaHeader.h include/AnitaConventions.h include/AnitaGeomTool.h include/AnitaPacketUtil.h include/simpleStructs.h
ICEMCO    = icemc.o vector.o position.o earthmodel.o balloon.o icemodel.o trigger.o signal.o ray.o Spectra.o anita.o roughness.o secondaries.o Primaries.o Tools.o counting.o Settings.o classdict.o Taumodel.o screen.o
ICEMCS    = icemc.cc vector.cc position.cc earthmodel.cc balloon.cc icemodel.cc trigger.cc signal.cc ray.cc Spectra.cc anita.cc roughness.cc secondaries.cc Primaries.cc Tools.cc counting.cc Settings.cc classdict.C Taumodel.cc screen.cc

ICEMC     = icemc$(ExeSuf)

OBJS          = $(CONDTRKO) $(ICEMCO) 

PROGRAMS      = $(ICEMC)


#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

##ANITADATALIB = libAnitaEvent.$(DllSuf)

##$(ANITADATALIB):
	@cd anita_data_format; make all; make install

all:            $(PROGRAMS)

$(ICEMC):       $(ICEMCO)
		$(LD) $(LDFLAGS) $(ICEMCO) $(LIBS) $(OutPutOpt) $(ICEMC)
		@echo "$@ done"


.PHONY: clean
clean:
		@rm -f $(OBJS) core classdict.* icemc
##@cd anita_data_format; make clean

distclean:      clean
		@rm -f $(PROGRAMS) $(ICEMCSO) $(ICEMCLIB) *dict.* *.def *.exp \
		   *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository core* classdict.* icemc

.PHONY: debug
debug: CXXFLAGS = $(DBGCXXFLAGS)
debug: LDFLAGS = -O0 -g -ggdb
debug: 		$(ICEMC)
		@echo "Compile in $@ mode done"


.PHONY: run
run:
	  ./$(ICEMC)



.SUFFIXES: .$(SrcSuf)

###

icemc.$(ObjSuf): 

classdict.C:	$(HEADERS)
	@echo "Generating dictionary…"
	@rm -f classdict*
	rootcint classdict.C -c $(ANITA_DATA_HEADERS) $(HEADERS) LinkDef.h

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
