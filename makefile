STD = -std=c++0x

CXX = g++ $(STD) -g
F77 = gfortran

FFLAGS   += -O3 -fPIC 
CXXFLAGS += -O3 -fPIC -fpermissive 
LDFLAGS +=  -fno-automatic -fno-f2c -O3 -g # -I./src/Inc -Iobj

# root
ROOTINCS = $(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 
ROOTARCH = $(findstring -m64, $(ROOTINCS) )

# LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# applgrid
APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLCLIBS    = $(shell applgrid-config --ldcflags)
APPLFLIBS    = $(shell applgrid-config --ldflags)

# hoppet
HOPPETLIBS =  $(shell hoppet-config --libs)

# fastjet
FASTJET     = $(shell fastjet-config --prefix)
#FASTJETLIBS = -L$(FASTJET)/lib -lfastjet
FASTJETLIBS = $(shell fastjet-config --libs)
FASTJETINCS = $(shell fastjet-config --cxxflags)
CXXFLAGS   += $(FASTJETINCS)

# get the fotran runtime library for linking fortran 
FRTLLIB = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLIB  = -L$(subst /libgfortran.a, ,$(FRTLLIB) ) -lgfortran

HOPPETINCS=$(shell hoppet-config --cxxflags)
HOPPETLIBS=$(shell hoppet-config --libs)

# now set up the compile and link flags and libs
#CXXFLAGS += $(ROOTARCH) $(ROOTINCS) $(APPLCXXFLAGS) $(LHAPDFINCS) $(FASTJETINCS) $(HOPPETINCS)
CXXFLAGS += $(ROOTARCH) $(ROOTINCS) $(APPLCXXFLAGS) $(LHAPDFINCS) $(HOPPETINCS)

CLIBS += $(APPLCLIBS) $(HOPPETLIBS)  $(FASTJETLIBS) $(ROOTLIBS) $(FRTLIB) $(LHAPDFLIBS)
FLIBS += $(LHAPDFLIBS) $(HOPPETLIBS) $(APPLFLIBS) $(ROOTLIBS) $(FRTLIB)


OFILE= MyFrame.o MyFrameData.o MyData.o MyCrossSection.o AtlasStyle.o lhapdf_string_interface.o theory_error_info.o OptionHandler.o MyPDF.o
OFILETEST= MyFrame.o MyFrameData.o MyData.o MyCrossSection.o AtlasStyle.o lhapdf_string_interface.o OptionHandler.o MyPDF.o


install : all
all : plot_old


plot_gpp: plot_gpp.o $(OFILE) 
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

plot_old: plot_old.o $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

plot_new: plot_new.o $(OFILETEST)
	$(CXX) $(LDFLAGS) -o $@ $(OFILETEST) $<  $(CLIBS)

plot_test: plot_old.o $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

stand: stand.o $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS) $(LHAPDFLIBS) 
	
makegridfromsherpa: makegridfromsherpa.o  $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

testgridfromsherpa: testgridfromsherpa.o  $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

testmypdf: testmypdf.o $(OFILETEST)
	$(CXX) $(LDFLAGS) -o $@ $(OFILETEST) $<  $(CLIBS) 


.SUFFIXES : .cxx .o .f .c

.f.o :
	$(F77) $(FFLAGS)   -c $<

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 


clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la  *~
