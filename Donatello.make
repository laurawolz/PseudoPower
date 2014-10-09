LIB=/home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc/lib/libhealpix_cxx.a /home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc/lib/libpsht.a  -lgsl -lgslcblas /home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc/lib/libcxxsupport.a -lm /home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc/lib/libc_utils.a    /home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc/lib/libfftpack.a  -lifcore -lstdc++ -liomp5 -lpthread -lcfitsio

HEALPIXPATH=/home/lwolz/Libraries/Healpix_donatello/Healpix_2.20a/src/cxx/linux_icc
GSLPATH=/usr/include

PATH1=-B$(HEALPIXPATH)/bin -I$(HEALPIXPATH)/include -L$(GSLPATH)
COMP=icpc

OPT=-openmp -O3 -parallel  #-par_report2


$(codename): $(codename).cpp header.h class_func.cpp 
	$(COMP) $(OPT) -c class_func.cpp $(PATH1) -o class_func.o
	$(COMP) $(OPT) -c $(codename).cpp $(PATH1) -o $(codename).o
	$(COMP) $(OPT) -o $(codename)  $(codename).o class_func.o $(PATH1) $(LIB) 