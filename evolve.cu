__global__ void kernelSetLogTable(double *logTable_d)
{
    int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                    (blockIdx.x*blockDim.x*blockDim.y);
    double temp= idx/1000000.0;
    if(idx==0)
        logTable_d[0]=log(0.0000001);
    else 
        logTable_d[idx]=log(temp);
}
__global__ void kernelSet_kx(double *kx_d,double dkx,int nx)
{
    int idx=threadIdx.x;
    int temp=nx/2;
    if(idx<temp)
        {kx_d[idx] = idx*dkx;}
    else 
        {kx_d[idx] = (idx - nx)*dkx;}
}
__global__ void kernelSet_ky(double *ky_d,double dky,int ny)
{
    int idx=threadIdx.x;
    int temp=ny/2;
    if(idx<temp)
        {ky_d[idx] = idx*dky;}
    else 
        {ky_d[idx] = (idx - ny)*dky;}
}
__global__ void kernelSet_kz(double *kz_d,double dkz,int nz)
{
    int idx=threadIdx.x;
    int temp=nz/2;
    if(idx<temp)
        {kz_d[idx] = idx*dkz;}
    else 
        {kz_d[idx] = (idx - nz)*dkz;}
}

__global__ void kernelComplexToReal(cuDoubleComplex* dfdc_d,double* tempreal_d)
        {
            int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                            (blockIdx.x*blockDim.x*blockDim.y);
            tempreal_d[idx]=dfdc_d[idx].x;   
		}
		

/////////////		
//INSIDE TIME LOOP
__global__ void kernelFind_dfdc_d( cuDoubleComplex* dfdc_d,double* logTable_d,double A )
   {
       int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                (blockIdx.x*blockDim.x*blockDim.y);
       double compTemp = dfdc_d[idx].x; 
       int logTemp1 = compTemp*1000000.0;
       int logTemp2 = (1.0-compTemp)*1000000.0;
       //dfdc_d[idx].x = 2.0*A*compTemp*(1.0-compTemp)*(1.0-2*compTemp);
       dfdc_d[idx].x = A*(1.0-2*compTemp) + logTable_d[logTemp1] - logTable_d[logTemp2];
       dfdc_d[idx].y = 0.0;
   }   
__global__ void kernelFind_dfdcPlus2kkkc(cuDoubleComplex* dfdc_d,cuDoubleComplex* comp_d,
									double* kx_d,double* ky_d,double* kz_d,double kappa_c,int nx,int ny,int nz)
	{
        int i= threadIdx.x + blockIdx.x*blockDim.x;
        int j= threadIdx.y + blockIdx.y*blockDim.y;
        int k= threadIdx.z + blockIdx.z*blockDim.z;
        int idx = k+nz*(j+i*ny);
			
			double kpow2 = kx_d[i] * kx_d[i] + ky_d[j] * ky_d[j] + kz_d[k] * kz_d[k];
		    dfdc_d[idx].x = dfdc_d[idx].x + 2.0 * kappa_c * kpow2 * comp_d[idx].x;
		    dfdc_d[idx].y = dfdc_d[idx].y + 2.0 * kappa_c * kpow2 * comp_d[idx].y;
	}		
			
//DERIVATIVE DFDCY
__global__ void kernelFind_dfdcx_d(cuDoubleComplex* dfdc_d, cuDoubleComplex* temp1,double* kx_d,int nx,int ny,int nz)
	{
			int i= threadIdx.x + blockIdx.x*blockDim.x;
        	int j= threadIdx.y + blockIdx.y*blockDim.y;
        	int k= threadIdx.z + blockIdx.z*blockDim.z;
        	int idx = k+nz*(j+i*ny);
		
		    temp1[idx].x = -1.0 * kx_d[i] * dfdc_d[idx].y;
		    temp1[idx].y     =    kx_d[i] * dfdc_d[idx].x;
	}


__global__ void kernelNormalize_dfdcx_multiplyWithMobility(cuDoubleComplex* temp1, double* tempreal_d, double one_by_nxnynz)
	{
		int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                            (blockIdx.x*blockDim.x*blockDim.y);
        	temp1[idx].x *= one_by_nxnynz;
        	temp1[idx].y *= one_by_nxnynz;

			double c = tempreal_d[idx];
		   	temp1[idx].x *= (c* (1.0 -c) );
           	temp1[idx].y *= (c* (1.0 -c) );
           
    }
__global__ void kernelStore_dfdcx_to_temp2(cuDoubleComplex* temp1,cuDoubleComplex* temp2, double *kx_d,int nx,int ny,int nz)	
	{
			int i= threadIdx.x + blockIdx.x*blockDim.x;
        	int j= threadIdx.y + blockIdx.y*blockDim.y;
        	int k= threadIdx.z + blockIdx.z*blockDim.z;
        	int idx = k+nz*(j+i*ny);
		
	 		temp2[idx].x  = (-1.0 * kx_d[i] * temp1[idx].y) ;
            temp2[idx].y  =         kx_d[i] * temp1[idx].x  ;
	}

//DERIVATIVE DFDCY
__global__ void kernelFind_dfdcy_d(cuDoubleComplex* dfdc_d, cuDoubleComplex* temp1,double* ky_d,int nx,int ny,int nz)
	{
			int i= threadIdx.x + blockIdx.x*blockDim.x;
        	int j= threadIdx.y + blockIdx.y*blockDim.y;
        	int k= threadIdx.z + blockIdx.z*blockDim.z;
        	int idx = k+nz*(j+i*ny);
		
		    temp1[idx].x = -1.0 * ky_d[j] * dfdc_d[idx].y;
		    temp1[idx].y     =    ky_d[j] * dfdc_d[idx].x;
	}

__global__ void kernelNormalize_dfdcy_multiplyWithMobility(cuDoubleComplex* temp1, double* tempreal_d, double one_by_nxnynz)
	{
		int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                            (blockIdx.x*blockDim.x*blockDim.y);
        	temp1[idx].x *= one_by_nxnynz;
        	temp1[idx].y *= one_by_nxnynz;

			double c = tempreal_d[idx];
		   	temp1[idx].x *= (c* (1.0 -c) );
           	temp1[idx].y *= (c* (1.0 -c) );
    }

__global__ void kernelStore_dfdcy_to_temp2(cuDoubleComplex* temp1,cuDoubleComplex* temp2, double *ky_d,int nx,int ny,int nz)	
	{
			int i= threadIdx.x + blockIdx.x*blockDim.x;
        	int j= threadIdx.y + blockIdx.y*blockDim.y;
        	int k= threadIdx.z + blockIdx.z*blockDim.z;
        	int idx = k+nz*(j+i*ny);
		
	 		temp2[idx].x  = temp2[idx].x + (-1.0 * ky_d[j] * temp1[idx].y) ;
            temp2[idx].y  = temp2[idx].y +  ky_d[j] * temp1[idx].x  ;
	}

//DERIVATIVE DFDCZ
__global__ void kernelFind_dfdcz_d(cuDoubleComplex* dfdc_d, cuDoubleComplex* temp1,double* kz_d,int nx,int ny,int nz)
	{
			int i= threadIdx.x + blockIdx.x*blockDim.x;
        	int j= threadIdx.y + blockIdx.y*blockDim.y;
        	int k= threadIdx.z + blockIdx.z*blockDim.z;
        	int idx = k+nz*(j+i*ny);
		
		    temp1[idx].x = -1.0 * kz_d[k] * dfdc_d[idx].y;
		    temp1[idx].y     =    kz_d[k] * dfdc_d[idx].x;
	}

__global__ void kernelNormalize_dfdcz_multiplyWithMobility(cuDoubleComplex* temp1, double* tempreal_d, double one_by_nxnynz)
	{
		int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                            (blockIdx.x*blockDim.x*blockDim.y);
        	temp1[idx].x *= one_by_nxnynz;
        	temp1[idx].y *= one_by_nxnynz;

			double c = tempreal_d[idx];
		   	temp1[idx].x *= (c* (1.0 -c) );
           	temp1[idx].y *= (c* (1.0 -c) );
    }

__global__ void kernelStore_dfdcz_to_temp2(cuDoubleComplex* temp1,cuDoubleComplex* temp2, double *kz_d,int nx,int ny,int nz)	
	{
			int i= threadIdx.x + blockIdx.x*blockDim.x;
        	int j= threadIdx.y + blockIdx.y*blockDim.y;
        	int k= threadIdx.z + blockIdx.z*blockDim.z;
        	int idx = k+nz*(j+i*ny);
		
	 		temp2[idx].x  = temp2[idx].x + (-1.0 * kz_d[k] * temp1[idx].y) ;
            temp2[idx].y  = temp2[idx].y +  kz_d[k] * temp1[idx].x  ;
	}


__global__ void kernelStep4(cuDoubleComplex* dfdc_d,cuDoubleComplex* comp_d,cuDoubleComplex* temp2,
								double* kx_d,double* ky_d,double* kz_d,
								 double P,double kappa_c,double dt,int nx,int ny,int nz)
	{
        int i= threadIdx.x + blockIdx.x*blockDim.x;
        int j= threadIdx.y + blockIdx.y*blockDim.y;
        int k= threadIdx.z + blockIdx.z*blockDim.z;
		int idx = k+nz*(j+i*ny);
			
			double kpow2 = kx_d[i] * kx_d[i] + ky_d[j] * ky_d[j] + kz_d[k] * kz_d[k];
		    double kpow4 = kpow2 * kpow2;

           dfdc_d[idx].x = temp2[idx].x;
           dfdc_d[idx].y = temp2[idx].y;

           
			double lhs = 1.0 + P * kappa_c * kpow4 * dt;
			cuDoubleComplex rhs;
			cuDoubleComplex ExPart;	
			ExPart.x = P * kappa_c * kpow4 * comp_d[idx].x + dfdc_d[idx].x;			
			ExPart.y = P * kappa_c * kpow4 * comp_d[idx].y + dfdc_d[idx].y;
		    // (P * kappa_c * kpow4 * comp[] + dfdc[]) is non-linear term (implicit)
		    rhs.x = comp_d[idx].x + dt * ExPart.x;
		    rhs.y = comp_d[idx].y + dt * ExPart.y;

		    comp_d[idx].x = rhs.x / lhs;
		    comp_d[idx].y = rhs.y / lhs;
		    
		    dfdc_d[idx].x = comp_d[idx].x;
		    dfdc_d[idx].y = comp_d[idx].y;
	}

__global__ void kernelNormalizingGroup(cuDoubleComplex* dfdc_d,double one_by_nxnynz,int* exit_variable_d)
{
    int idx= (threadIdx.x+threadIdx.y*blockDim.x) +
                        (blockIdx.x*blockDim.x*blockDim.y);
    dfdc_d[idx].x *= one_by_nxnynz;
    dfdc_d[idx].y *= one_by_nxnynz;
    if (dfdc_d[idx].x < -0.2    ||  dfdc_d[idx].x > 1.2) {
        *exit_variable_d=1;}
}
void Evolve()
{
	void Output_Conf(int );

    int count;
    double dkx, dky, dkz;
    double err, maxerror;
    double total;
    double *kx_d,*ky_d,*kz_d;
    double *tempreal_d;
	double *logTable_d;
	double *tempreal;
	int  exit_variable;
	int *exit_variable_d;

	cuDoubleComplex *temp1,*temp2;

	//GPU TIMER VARIABLES
	float gpuElapsedTime;
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);


	//cudaMalloc ((void **)&dfdcx_d, nx*ny*nz*sizeof(cuDoubleComplex));
    //cudaMalloc ((void **)&dfdcy_d, nx*ny*nz*sizeof(cuDoubleComplex));
    //cudaMalloc ((void **)&dfdcz_d, nx*ny*nz*sizeof(cuDoubleComplex));
 
    cudaMalloc ((void **)&temp1, nx*ny*nz*sizeof(cuDoubleComplex));
    cudaMalloc ((void **)&temp2, nx*ny*nz*sizeof(cuDoubleComplex));



	cudaMalloc ((void**)&kx_d,sizeof(double)*nx);
    cudaMalloc ((void**)&ky_d,sizeof(double)*ny);
    cudaMalloc ((void**)&kz_d,sizeof(double)*nz);
	cudaMalloc ((void**)&tempreal_d,sizeof(double)*nx*ny*nz);
	cudaMalloc ((void**)&logTable_d,sizeof(double)*1024000);
 	cudaMalloc ((void**)&exit_variable_d,sizeof(int));
    
    tempreal = (double *) malloc(sizeof(double) * nx * ny * nz);

//Computation starts from here
	
	alloycomp=0.5;
	exit_variable=0;
    cudaMemcpy(exit_variable_d,&exit_variable,sizeof(int),cudaMemcpyHostToDevice);

    dkx = 2.0 * PI / ((double) nx * dx);
    dky = 2.0 * PI / ((double) ny * dy);
    dkz = 2.0 * PI / ((double) nz * dz);
//TO MAKE LOG TABLE

kernelSetLogTable<<<dim3(1000,1,1),dim3(32,32,1)>>>(logTable_d);

// TO FIND KX,KY,KZ 
//KERNELS TO FIND kx_d,ky_d,kz_d
    kernelSet_kx<<<dim3(1,1,1),dim3(nx,1,1)>>>(kx_d,dkx,nx);
    kernelSet_ky<<<dim3(1,1,1),dim3(ny,1,1)>>>(ky_d,dky,ny);
    kernelSet_kz<<<dim3(1,1,1),dim3(nz,1,1)>>>(kz_d,dkz,nz);

    
	cufftExecZ2Z(plan,comp_d,comp_d,CUFFT_FORWARD);
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  	

//TIME LOOP STARTS
for (count = 0; count <= num_steps; count++) {

	printf("Iteration= %d\n",count+1);

	if(  (count % 10000) == 0 ) 
		{Output_Conf(count);}
	
    kernelComplexToReal<<<dim3(blocks,1,1),dim3(32,32,1)>>>(dfdc_d,tempreal_d);
	
    
	
   
	kernelFind_dfdc_d<<<dim3(blocks,1,1),dim3(32,32,1)>>>(dfdc_d,logTable_d,A);
		cufftExecZ2Z(plan,dfdc_d,dfdc_d,CUFFT_FORWARD);
	kernelFind_dfdcPlus2kkkc<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(dfdc_d,comp_d,kx_d,ky_d,kz_d,kappa_c,nx,ny,nz);
	

	kernelFind_dfdcx_d<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(dfdc_d, temp1,kx_d,nx,ny,nz);
	cufftExecZ2Z(plan,temp1,temp1,CUFFT_INVERSE);
	kernelNormalize_dfdcx_multiplyWithMobility<<<dim3(blocks,1,1),dim3(32,32,1)>>>(temp1,tempreal_d,one_by_nxnynz);
	cufftExecZ2Z(plan,temp1,temp1,CUFFT_FORWARD);
	kernelStore_dfdcx_to_temp2<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(temp1,temp2, kx_d,nx,ny,nz);
    

	kernelFind_dfdcy_d<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(dfdc_d, temp1,ky_d,nx,ny,nz);
	cufftExecZ2Z(plan,temp1,temp1,CUFFT_INVERSE);
	kernelNormalize_dfdcy_multiplyWithMobility<<<dim3(blocks,1,1),dim3(32,32,1)>>>(temp1,tempreal_d,one_by_nxnynz);
	cufftExecZ2Z(plan,temp1,temp1,CUFFT_FORWARD);
	kernelStore_dfdcy_to_temp2<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(temp1,temp2, ky_d,nx,ny,nz);
    
    kernelFind_dfdcz_d<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(dfdc_d, temp1,kz_d,nx,ny,nz);
	cufftExecZ2Z(plan,temp1,temp1,CUFFT_INVERSE);
	kernelNormalize_dfdcz_multiplyWithMobility<<<dim3(blocks,1,1),dim3(32,32,1)>>>(temp1,tempreal_d,one_by_nxnynz);
	cufftExecZ2Z(plan,temp1,temp1,CUFFT_FORWARD);
	kernelStore_dfdcz_to_temp2<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(temp1,temp2, kz_d,nx,ny,nz);
    


	kernelStep4<<<dim3(nx/8,ny/8,nz/8),dim3(8,8,8)>>>(dfdc_d,comp_d,temp2,kx_d,ky_d,kz_d,P,kappa_c,dt,nx,ny,nz);

	cudaMemcpy(dfdc,dfdc_d,sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);
	
	/* Check for conservation of mass */
	total = dfdc[0].x * one_by_nxnynz;
	err = fabs(total - alloycomp);
	if (err > COMPERR) {
	    printf("ELEMENTS ARE NOT CONSERVED,SORRY!!!!\n");
	    printf("error=%lf\n", err);
	    exit(1);
	}


	cufftExecZ2Z(plan,dfdc_d,dfdc_d, CUFFT_INVERSE);
    
	kernelNormalizingGroup<<<dim3(blocks,1,1),dim3(32,32,1)>>>(dfdc_d,one_by_nxnynz,exit_variable_d);
	
	cudaMemcpy(&exit_variable,exit_variable_d,sizeof(int),cudaMemcpyDeviceToHost);
        if(exit_variable==1)
        	{printf("Compositions are out of bound");
        		Output_Conf(count);
              exit(1);}
 //Check for convergence 
		maxerror = 0.0;

		cudaMemcpy(dfdc,dfdc_d,sizeof(cuDoubleComplex)*nx*ny*nz,cudaMemcpyDeviceToHost);
		cudaMemcpy(tempreal,tempreal_d,sizeof(double)*nx*ny*nz,cudaMemcpyDeviceToHost);

      	int i,j,k,idx;
        	for (k = 0; k < nz; k++) {
          	   for ( j = 0; j < ny; j++) {
           		for ( i = 0; i < nx; i++) 
           {
			   	idx= i + nx * (j + ny*k);
		    	err = fabs( tempreal[idx] - dfdc[idx].x );
				if (err > maxerror)
				maxerror = err;
			}}}

		if (maxerror <= Tolerance) {
	    	printf("maxerror=%lf\tnumbersteps=%d\n", maxerror, count);
	   		 break;}
	cudaThreadSynchronize();
  }


  cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuElapsedTime,start,stop);
    printf("Elapsed time on Gpu %f ms\n",gpuElapsedTime);

  //time loop ends
 cudaFree(dfdcx_d);
 cudaFree(dfdcy_d);
 cudaFree(dfdcz_d);
 cudaFree(kx_d);
 cudaFree(ky_d);
 cudaFree(kz_d);
 cudaFree(tempreal_d);
 cudaFree(logTable_d);
 cudaFree(exit_variable_d);
 
}
