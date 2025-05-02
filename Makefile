

compiler = ifort #gfortran

OBJ = mod_parameters.o    \
	  mod_read_file.o     \
	  mod_utility.o       \
	  mod_grid.o		  \
	  mod_field.o		  \
	  mod_solver.o		  \
   	  main.o
# Has to be compiled in the above sequence since later are dependent on the previous ones.

compile: $(OBJ)
	$(compiler) -qopenmp -o run.exe *.o
# use -qopenmp for ifort compiler

%.o: %.F90
	$(compiler) -qopenmp -c $< 
# use -qopenmp for ifort compiler here also to use openmp for parallelization.


run:
	sbatch ./submit_job.sh
#	./run.exe >output.log 2>error.log

clean:
	rm *.o *.mod