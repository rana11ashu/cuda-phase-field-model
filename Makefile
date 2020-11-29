SHELL = /bin/sh
CC = gcc -std=c99
compops =  -O3 -fopenmp 
linkops =  -L/usr/local/lib -lfftw3_omp -lfftw3 -lm
objs = main.o get_input.o init_conf.o evolve.o out_conf.o
headers = binary.h

%.o: %.c
	 $(CC) $(compops) -c $<

all:spVM3D.out
spVM3D.out: $(objs) $(headers)
	 $(CC) -o spVM3D.out $(objs) $(compops) $(linkops)

clean:
	-\rm  *.o *.out
# End of the makefile
