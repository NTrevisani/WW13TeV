########################################################
# Makefile for compiling a TObject derived ROOT class
#
# from the example Makefile in ROOT's test directory

########################################################
#Make sure that your Environment variable LD_LIBRARY_PATH
#has been correctly set.

# Package name
TARGET = MuonIsoMVA

########################################################
# The rest of the makefile does not need to be modified

#ROOTLIBS      = -L$(ROOTSYS)/lib -lNew -lBase -lCint -lClib -lCont -lFunc \
#                -lGraf -lGraf3d -lHist -lHtml -lMatrix -lMeta -lMinuit -lNet \
#                -lPostscript -lProof -lTree -lUnix -lZip
#ROOTGLIBS     = -lGpad -lGui -lGX11 -lX3d
ROOTCONFIG   := root-config
ROOTLIBS      = $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)

TARGETHEAD=$(TARGET).h
TARGETSRC=$(TARGET).C
TARGETCINT=$(TARGET)Cint.cxx
TARGETCINTO=$(TARGET)Cint.o
TARGETO = $(TARGET).o
TARGETSO = $(TARGET)_C.so

# Linux/gcc-2.96 version
CXX           = g++
CXXINCLUDES   = -I$(ROOTSYS)/include -I/nfs/fanae/user/calderon/PAF
CXXFLAGS      = -O -Wall -fPIC $(CXXINCLUDES)
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -Wl,-soname,$(TARGETSO) -shared
LIBS          = $(ROOTLIBS) -lm -ldl -rdynamic
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib \
                -lXpm -lX11 -lm -ldl -rdynamic

all: $(TARGETSO)

$(TARGETCINT): $(TARGETHEAD)
	rootcint -f $(TARGETCINT) -c $(CXXINCLUDES) $(TARGETHEAD)

$(TARGETCINTO): $(TARGETCINT)
	$(CXX) $(CXXFLAGS) -c -o $(TARGETCINTO) $(TARGETCINT)

$(TARGETO): $(TARGETHEAD) $(TARGETSRC)

$(TARGETSO): $(TARGETO) $(TARGETCINTO)
	$(LD) $(SOFLAGS) $(LDFLAGS) $(TARGETO) $(TARGETCINTO) -o $(TARGETSO) $(GLIBS)

clean:
	rm -f $(TARGETO) $(TARGETCINTO) $(TARGETCINT) $(TARGET)Cint.h $(TARGETSO)

