
# PLEASE HERE INSERT YOUR OWN PATH TO YOUR HEALPIX LIBRARY
HEALPIXPATH=/home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc

LIB=$(HEALPIXPATH)/lib/libhealpix_cxx.a $(HEALPIXPATH)/lib/libpsht.a  -lgsl -lgslcblas $(HEALPIXPATH)/lib/libcxxsupport.a -lm $(HEALPIXPATH)/lib/libc_utils.a  $(HEALPIXPATH)/lib/libfftpack.a  -lifcore -lstdc++ -liomp5 -lpthread -lcfitsio

# In some case you need to link to a different cfitsio distribution if you have compiled your Healpix with the non-standard cfitsio library

# PLEASE HERE INSERT YOUR OWN PATH TO YOUR GSL LIBRARY

GSLPATH=/usr/include

PATH1=-B$(HEALPIXPATH)/bin -I$(HEALPIXPATH)/include -L$(GSLPATH)
COMP=icpc

OPT=-openmp -O3 -parallel  


$(codename): $(codename).cpp header.h class_func.cpp 
	$(COMP) $(OPT) -c class_func.cpp $(PATH1) -o class_func.o
	$(COMP) $(OPT) -c $(codename).cpp $(PATH1) -o $(codename).o
	$(COMP) $(OPT) -o $(codename)  $(codename).o class_func.o $(PATH1) $(LIB) 


