# Makefile for the ROOT test programs.  # This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

#------------------------------------------------------------------------------

################################################################################
# Site specific flags
################################################################################
# Toggle these as needed to get things to install

#BOOSTFLAGS = -I boost_1_48_0
# commented out for kingbee and older versions of gcc
ANITA3_EVENTREADER=1
ANITA3_EVENTCORRELATOR=1

# Uncomment to enable healpix 
#USE_HEALPIX=1

# Uncomment to disable explicit vectorization (but will do nothing if ANITA_UTIL is not available) 
#VECTORIZE=1


# The ROOT flags are added to the CXXFLAGS in the .arch file
# so this should be simpler...
ifeq (,$(findstring -std=, $(CXXFLAGS)))
ifeq ($(shell test $(GCC_MAJOR) -lt 5; echo $$?),0)
ifeq ($(shell test $(GCC_MINOR) -lt 5; echo $$?),0)
CXXFLAGS += -std=c++0x
else
CXXFLAGS += -std=c++11
endif
endif
endif

################################################################################

ifeq (,$(findstring -std=c++1, $(CXXFLAGS)))
# If not compiling with C++11 (or later) support, all occurrences of "constexpr"
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



# Uses the standard ANITA environment variable to figure
# out if ANITA libs are installed
ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_EXISTS=1

ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR)
LIBS_ANITA_UTIL=-lRootFftwWrapper 
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
ANITA_UTIL_ETC_DIR=$(ANITA_UTIL_INSTALL_DIR)/etc
endif

ifdef ANITA_UTIL_EXISTS
CXXFLAGS += -DANITA_UTIL_EXISTS
endif

ifdef VECTORIZE
CXXFLAGS += -DVECTORIZE -march=native -fabi-version=0
endif

ifdef ANITA3_EVENTREADER
CXXFLAGS += -DANITA3_EVENTREADER
LIBS_ANITA_UTIL:=-lAnitaEvent $(LIBS_ANITA_UTIL)
endif

ifdef ANITA3_EVENTCORRELATOR
CXXFLAGS += -DANITA3_EVENTCORRELATOR
LIBS_ANITA_UTIL:=-lAnitaCorrelator $(LIBS_ANITA_UTIL)
endif

ifdef USE_HEALPIX
	CXXFLAGS += -DUSE_HEALPIX `pkg-config --cflags healpix_cxx`
	LIBS  += `pkg-config --libs healpix_cxx` 
endif



GENERAL_FLAGS = -g -O2 -pipe -m64 -pthread
WARN_FLAGS = -W -Wall -Wextra -Woverloaded-virtual
# -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable

CXXFLAGS += $(GENERAL_FLAGS) $(CPPSTD_FLAGS) $(WARN_FLAGS) $(ROOTCFLAGS) $(INC_ANITA_UTIL)

DBGFLAGS  = -pipe -Wall -W -Woverloaded-virtual -g -ggdb -O0 -fno-inline

DBGCXXFLAGS = $(DBGFLAGS) $(ROOTCFLAGS) $(BOOSTFLAGS)
LDFLAGS  += $(CPPSTD_FLAGS) $(LD_ANITA_UTIL) -I$(BOOST_ROOT) -L.
LIBS += $(LIBS_ANITA_UTIL)

# Mathmore not included in the standard ROOT libs
LIBS += -lMathMore

CLASS_HEADERS = rx.hpp Taumodel.hh Settings.h blazars/fava.h icemc_random.h 
DICT = classdict

OBJS = vector.o position.o earthmodel.o balloon.o icemodel.o signal.o ray.o Spectra.o anita.o roughness.o secondaries.o Primaries.o Tools.o counting.o $(DICT).o Settings.o Taumodel.o screen.o GlobalTrigger.o ChanTrigger.o SimulatedSignal.o EnvironmentVariable.o source.o  random.o


BINARIES = icemc$(ExeSuf) testTrigger$(ExeSuf) testSettings$(ExeSuf) 

ifdef ANITA_UTIL_EXISTS
BINARIES += testEAS$(ExeSuf) testWAIS$(ExeSuf) testInputAfterAntenna$(ExeSuf) testThermalNoise$(ExeSuf)
endif


#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(BINARIES)

$(BINARIES): %: %.$(SrcSuf) $(OBJS)
		$(LD) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $< $(LIBS) $(OutPutOpt) $@
		@echo "$@ done"

libicemc.so: $(OBJS) 
		$(LD) $(CXXFLAGS) -shared $(LDFLAGS) $(OBJS) $(LIBS) $(OutPutOpt) $@
		@echo "$@ done"


.PHONY: clean
clean:
		@rm -f $(OBJS) classdict.* $(BINARIES) *.so 

distclean:      clean
		@rm -f $(OBJS) $(BINAIRES) $(DICT)* *.def *.exp \
		   *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository core* classdict.* icemc

$(DICT).C : $(CLASS_HEADERS) LinkDef.h
		@echo "<**And here's the dictionary...**>" $<
		@rm -f *Dict*
		rootcint -f $@  -c -p -I./ $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h


%.$(ObjSuf) : %.$(SrcSuf) %.h
	@echo "<**Compiling**> "$<
	$(LD) $(CXXFLAGS) -c $< -o $@

%.$(ObjSuf) : %.C
	@echo "<**Compiling**> "$<
	$(LD) $(CXXFLAGS) $ -c $< -o  $@
