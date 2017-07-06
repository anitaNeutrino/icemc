# Makefile for the ROOT test programs.  # This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

#------------------------------------------------------------------------------
#brings in the few minossoft things i need
CXXFLAGS += -g -O2 #-I. #$(TEST_INCLUDES)

#LIBS += -stdlib=libstdc++ -lstdc++

ICEMCO    = makelimitplot.o
ICEMCS    =  makelimitplot.cc
ICEMC     = makelimitplot$(ExeSuf)

OBJS          = $(CONDTRKO) $(ICEMCO) 

PROGRAMS      = $(ICEMC)


#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(ICEMC):       $(ICEMCO)
		$(LD) $(LDFLAGS) $(ICEMCO) $(LIBS) \
		   $(OutPutOpt)$(ICEMC)
		@echo "$@ done"


clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS) $(ICEMCSO) $(ICEMCLIB) *Dict.* *.def *.exp \
		   *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository core*

.SUFFIXES: .$(SrcSuf)

###

#icemc.$(ObjSuf): 
makelimitplot.$(ObjSuf): 



.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
