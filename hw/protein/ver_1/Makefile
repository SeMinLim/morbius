LIBPATH=../../
BOARD=vc707
BUILDTOOLS=$(LIBPATH)/buildtools/

BLIBPATH=$(LIBPATH)/../bluelib/src/

GIBBSPATH=gibbssampler

CUSTOMBSV= -p +:$(LIBPATH)/src/:$(LIBPATH)/dram/src:$(BLIBPATH)/:$(BLIBPATH)/msfp/:$(GIBBSPATH) -steps-max-intervals 3000000
CUSTOMCPP_BSIM= $(BLIBPATH)/bdpi.cpp

include $(BUILDTOOLS)/Makefile.base
