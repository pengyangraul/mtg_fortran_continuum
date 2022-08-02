# This is an commentary line in a makefile
# Start of the makefile
#
# Defining variables
programpath=${HOME}
objects = mat_module.o indnth.o  io_module.o main.o
fc = ifort
compile_op = -heap-arrays -O3 -qopenmp  -parallel -I${programpath}/include -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include 
link_op = -heap-arrays -O3 -parallel \
 ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -ldl -L${programpath}/lib -lslatec

# makefile
execname: $(objects)
	$(fc) $(objects) $(link_op) -o main 
%.o: %.f90
	$(fc) $(compile_op) -c $<

clean:
	rm *.mod
	rm $(objects)
# End of the makefile
