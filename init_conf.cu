__global__ void kernelSetRandomNumbers(double *randnum_d)
{
    int idx= (threadIdx.x+threadIdx.y*blockDim.x)+
                     (blockIdx.x*blockDim.x*blockDim.y);
    randnum_d[idx] = (2.0*randnum_d[idx]- 1.0);
} 

//Kernel Init_Conf to FIND INITIAL COMP
__global__ void kernelSetInitialComp(double* randnum_d,
                                        cuDoubleComplex* comp_d,cuDoubleComplex* dfdc_d,
                                        double mean,double alloycomp,double noise_level)
{
    int idx=(threadIdx.x+threadIdx.y*blockDim.x)+
                    (blockIdx.x*blockDim.x*blockDim.y);
    randnum_d[idx] = randnum_d[idx] - mean;
    comp_d[idx].x  = alloycomp + noise_level*randnum_d[idx];
    comp_d[idx].y  = 0.0;
    dfdc_d[idx].x  = comp_d[idx].x;
    dfdc_d[idx].y  = comp_d[idx].y;
}

void Init_Conf(){
    
    //DECLARATION OF LOCAL VARIABLE
        double *randnum_d;
        double  sum=0.0, mean=0.0;
        cudaMalloc((void **)&randnum_d, nx*ny*nz*sizeof(double));
        double *sum_d;
        cudaMalloc((void **)&sum_d, sizeof(double));
        double *d_temp_storage = NULL;
        size_t  temp_storage_bytes = 0;
    
    //GENERATE RANDOM NUMBERS USING cuRAND
    curandGenerator_t  gen;
    curandCreateGenerator (&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed (gen,1234ULL);
    curandGenerateUniformDouble(gen, randnum_d, nx*ny*nz);
    
    //KERNEL TO SCALE RANDOM NUMBERS
    kernelSetRandomNumbers<<<dim3(blocks,1,1),
                                  dim3(32,32,1)>>>
                                    (randnum_d);

    //REDUCTION(SUM) USING CUB LIBRARY                               
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, randnum_d, sum_d, nx*ny*nz );
    cudaMalloc((void**)&d_temp_storage, temp_storage_bytes);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, randnum_d, sum_d, nx*ny*nz );
    
    //FIND MEAN
    cudaMemcpy(&sum,sum_d,sizeof(double),cudaMemcpyDeviceToHost);    
    printf("sum = %lf\n", sum);
    mean = sum /(double)(nx*ny*nz);
    printf("Mean = %lf\n", mean);

    //KERNEL TO FIND INTIAL COMP
    kernelSetInitialComp<<<dim3(blocks,1,1),
                            dim3(32,32,1)>>>
                                (randnum_d,comp_d,dfdc_d,mean,alloycomp,noise_level);
    //REDUCTION(SUM) USING CUB LIBRARY TO VERIFY MEAN    
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, randnum_d, sum_d, nx*ny*nz);
    cudaMalloc((void**)&d_temp_storage, temp_storage_bytes);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, randnum_d, sum_d, nx*ny*nz);

    //FIND MEAN
    cudaMemcpy(&sum,sum_d,sizeof(double),cudaMemcpyDeviceToHost);
    mean = sum /(double)(nx*ny*nz);
    printf("Mean = %f\n", mean);

    //DESTROY LOCAL VARIABLES
    cudaFree(randnum_d);
    cudaFree(d_temp_storage);
}

/*{
    FILE *fp;
    char fn[100];
    double *random_num;
    double sum, mean;

    random_num = (double *) malloc(sizeof(double) * nx * ny * nz);
    srand(time(NULL));
    sum = 0.0;
    for (int k = 0; k < nz; k++) {
	for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
		random_num[i + nx * (j + ny * k)] =
		    (double) rand() / RAND_MAX;
		random_num[i + nx * (j + ny * k)] =
		    2.0 * random_num[i + nx * (j + ny * k)] - 1.0;
		random_num[i + nx * (j + ny * k)] *= noise_level;
		sum += random_num[i + nx * (j + ny * k)];
	    }
	}
    }


    mean = sum * one_by_nxnynz;

    for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
		comp[i + nx * (j + ny * k)] =
		    alloycomp + random_num[i + nx * (j + ny * k)] - mean;
	    }
	  }
    }

    for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
		dfdc[i + nx * (j + ny * k)] = comp[i + nx * (j + ny * k)];
	    }
	  }
    }

    sprintf(fn, "profile.in");
    if (!(fp = fopen(fn, "w"))) {
	printf("File:%s could not be opened \n", fn);
	exit(1);
    }
    for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
		fprintf(fp, "%d\t%d\t%d\t%lf\n", i, j, k, creal(comp[i + nx * (j + ny * k)]));
	    }
	    fprintf(fp, "\n");
	  }
	fprintf(fp, "\n");
    }
    fclose(fp);
    free(random_num);
}*/

/*void Read_Restart(void)
{
    FILE *fpread;
    char fr[100];

    sprintf(fr, "conf.%06d", initcount);
    fpread = fopen(fr, "r");
    if (fread(&comp[0], sizeof(double), 2 * nx * ny * nz, fpread));
    fclose(fpread);

    for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
		dfdc[i + nx * (j + ny * k)] = comp[i + nx * (j + ny * k)];
	    }
	  }
    }
}*/
