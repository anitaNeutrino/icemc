# Makefile for the ROOT test programs.  # This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include so-dep_.inc # Updated by emacs: SO_DEP := canvas.cpp
SO_DEP_STEM = $(basename $(SO_DEP))
SO_TARGET := $(SO_DEP_STEM).so
HASH = $(shell md5sum $(SO_DEP) | gawk '{print $$1}')
SO_HOT_TARGET := SO/$(SO_DEP_STEM)_$(HASH)_.so

hot: $(SO_HOT_TARGET)

include Makefile.arch

#------------------------------------------------------------------------------

################################################################################
# Site specific flags
################################################################################
# Toggle these as needed to get things to install

#BOOSTFLAGS = -I boost_1_48_0
# commented out for kingbee and older versions of gcc
# dummy change to test changed git origin.
ANITA3_EVENTREADER=1





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
LD_ANITA_UTIL=-L$(ANITA_UTIL_LIB_DIR) -lAnitaEvent -lfftw3 -lRootFftwWrapper
INC_ANITA_UTIL=-I$(ANITA_UTIL_INC_DIR)
ANITA_UTIL_ETC_DIR=$(ANITA_UTIL_INSTALL_DIR)/etc
endif

ifdef ANITA_UTIL_EXISTS
CXXFLAGS += -DANITA_UTIL_EXISTS
endif

ifdef ANITA3_EVENTREADER
CXXFLAGS += -DANITA3_EVENTREADER
endif

GENERAL_FLAGS = -g3 -O0 -pipe -m64 -pthread
WARN_FLAGS = -W -Wall -Wextra -Woverloaded-virtual
# -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable

CXXFLAGS += $(GENERAL_FLAGS) $(CPPSTD_FLAGS) $(WARN_FLAGS) $(ROOTCFLAGS) $(INC_ANITA_UTIL)

DBGFLAGS  = -pipe -Wall -W -Woverloaded-virtual -g3 -O0 -fno-inline

DBGCXXFLAGS = $(DBGFLAGS) $(ROOTCFLAGS) $(BOOSTFLAGS)
LDFLAGS  += $(CPPSTD_FLAGS) $(LD_ANITA_UTIL) -I$(BOOST_ROOT) -L.

# Mathmore not included in the standard ROOT libs
LIBS += -lMathMore  -lX11

CLASS_HEADERS = rx.hpp Taumodel.hh Settings.h
DICT = classdict

OBJS = vector.o position.o earthmodel.o balloon.o icemodel.o signal.o ray.o Spectra.o anita.o roughness.o secondaries.o Primaries.o Tools.o counting.o $(DICT).o Settings.o Taumodel.o screen.o GlobalTrigger.o ChanTrigger.o SimulatedSignal.o EnvironmentVariable.o hot-loop.o


# BINARIES = icemc$(ExeSuf) testTrigger$(ExeSuf) testSettings$(ExeSuf) testEAS$(ExeSuf) testInputAfterAntenna$(ExeSuf) testThermalNoise$(ExeSuf)
BINARIES = testEAS$(ExeSuf) main
# BINARIES = main skel

BVVCPPFLAGS =         -pedantic -Wall -O0 -g3 -fPIC $(shell root-config --cflags)
BVVLDLIBS  = -ldl $(shell root-config --glibs)


#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(BINARIES) $(SO_TARGET)

# testEAS : testEAS.cc hot-api.h hot-loop.o
# 		g++ -pedantic -Wall -O0 -g3 -fPIC `root-config --cflags` -rdynamic -o $@ $< -ldl `root-config --glibs` hot-loop.o $(LIBS)


$(SO_HOT_TARGET): $(SO_DEP)
	g++ $(CXXFLAGS) -shared $(LDFLAGS) -ldl -fvisibility=hidden -o $@ $<
	nm $@ > $(basename $(SO_HOT_TARGET)).txt
	./notify.sh testEAS SO_LOCATION $@ || rm $(SO_HOT_TARGET) $(basename $(SO_HOT_TARGET)).txt # if default, "hot" target was called by mistake when program is not running.

$(SO_TARGET) : $(SO_DEP) hot-api.h
# g++ $(CXXFLAGS) -shared $(LDFLAGS) -fvisibility=hidden -o $@ $<
	g++ $(CXXFLAGS) -g3 -shared $(LDFLAGS) -fvisibility=hidden -o $@ $<


hot-loop.o : hot-loop.cpp hot-api.h
		g++ -c $(BVVCPPFLAGS) $(BVVLDFLAGS) -rdynamic -o $@ $< $(BVVLDLIBS)

# main : main.cpp hot-api.h hot-loop.o
# 		c++ $(BVVCPPFLAGS) $(BVVLDFLAGS) -rdynamic -o $@ $< $(BVVLDLIBS) hot-loop.o
# 
$(BINARIES): %: %.$(SrcSuf) $(OBJS)
		$(LD) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) $< $(OutPutOpt) $@
		@echo "$@ done"

.PHONY: clean
clean:
		@rm -f $(OBJS) classdict.* $(BINARIES) *.so

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
