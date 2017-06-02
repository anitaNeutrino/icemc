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
















################################################################################
ifeq (,$(findstring -std=c++1, $(ROOTCFLAGS)))
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
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR) -lAnitaEvent -lRootFftwWrapper
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
ANITA_UTIL_ETC_DIR=$(ANITA_UTIL_INSTALL_DIR)/etc
endif

ifdef ANITA_UTIL_EXISTS
CXXFLAGS += -DANITA_UTIL_EXISTS
endif

ifdef ANITA3_EVENTREADER
CXXFLAGS += -DANITA3_EVENTREADER
endif

GENERAL_FLAGS = -g -O2 -pipe -m64 -pthread
WARN_FLAGS = -W -Wall -Wextra -Woverloaded-virtual
# -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable

CXXFLAGS += $(GENERAL_FLAGS) $(CPPSTD_FLAGS) $(WARN_FLAGS) $(ROOTCFLAGS) $(INC_ANITA_UTIL)

DBGFLAGS  = -pipe -Wall -W -Woverloaded-virtual -g -ggdb -O0 -fno-inline

DBGCXXFLAGS = $(DBGFLAGS) $(ROOTCFLAGS) $(BOOSTFLAGS)
LDFLAGS  += $(CPPSTD_FLAGS) $(LD_ANITA_UTIL) -I$(BOOST_ROOT) -L.

# Mathmore not included in the standard ROOT libs
LIBS += -lMathMore

CLASS_HEADERS = rx.hpp Taumodel.hh
DICT = classdict

OBJS = vector.o position.o earthmodel.o balloon.o icemodel.o signal.o ray.o Spectra.o anita.o roughness.o secondaries.o Primaries.o Tools.o counting.o Settings.o $(DICT).o Taumodel.o screen.o GlobalTrigger.o ChanTrigger.o SimulatedSignal.o


BINARIES = icemc$(ExeSuf) testTrigger$(ExeSuf) testSettings$(ExeSuf) testEAS$(ExeSuf)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(BINARIES)

$(BINARIES): %: %.$(SrcSuf) $(OBJS)
		$(LD) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $< $(OutPutOpt) $@
		@echo "$@ done"

.PHONY: clean
clean:
		@rm -f $(OBJS) classdict.* $(BINARIES)

distclean:      clean
		@rm -f $(OBJS) $(BINAIRES) $(DICT)* *.def *.exp \
		   *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository core* classdict.* icemc

$(DICT).C : $(HEADERS)
		@echo "<**And here's the dictionary...**>" $<
		@rm -f *Dict*
		rootcint $@ -c -p -I./ $(INC_ANITA_UTIL) $(CLASS_HEADERS) LinkDef.h

%.$(ObjSuf) : %.$(SrcSuf) %.h
	@echo "<**Compiling**> "$<
	$(LD) $(CXXFLAGS) -c $< -o $@

%.$(ObjSuf) : %.C
	@echo "<**Compiling**> "$<
	$(LD) $(CXXFLAGS) $ -c $< -o  $@
