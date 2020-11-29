#! /bin/sh

nvcc -Xptxas -O3 -arch=sm_37 -I/usr/local/cuda-9.2/samples/common/inc main.cu -o spino.out -lcufft -lcuda -lcudart -lcurand
