#include"binary.h"
#include"get_input.cu"
#include"init_conf.cu"
#include"evolve.cu"
#include"out_conf.cu"

int main(void)
{
    void Get_Input_Parameters(char *fnin, char *fnout);
    void Init_Conf(void);
    void Read_Restart(void);
    void Evolve(void);
    void Output_Conf(int );

    char finput[15] = "bin1ary";
    char fnin[30], fnout[30];

    FILE *fp;

    if (!(fp = fopen(finput, "r"))) {
	printf("File:%s could not be opened\n", finput);
	exit(EXIT_FAILURE);
    }
    if (fscanf(fp, "%s", fnin) == 1) {
	printf("Input Parameters Filename:%s\n", fnin);
    }
    if (fscanf(fp, "%s", fnout) == 1) {
	printf("Output Parameters Filename:%s\n", fnout);
    }
    if (!(fpout = fopen(fnout, "w"))) {
	printf("File:%s could not be opened\n", fnout);
	exit(EXIT_FAILURE);
    }
    fclose(fp);

    Get_Input_Parameters(fnin, fnout);

    //comp   = (cuDoubleComplex*)malloc( (sizeof (cuDoubleComplex)*nx*ny*nz) );
    //dfdc   = (cuDoubleComplex*)malloc( (sizeof (cuDoubleComplex)*nx*ny*nz) );
   
    

    cudaMalloc ((void **)&comp_d, nx*ny*nz*sizeof(cuDoubleComplex));
    cudaMalloc ((void **)&dfdc_d, nx*ny*nz*sizeof(cuDoubleComplex));
    
    
    
    one_by_nxnynz = 1.0 / (double) (nx * ny * nz);
    blocks=(nx*ny*nz)/1024;
    
    cufftPlan3d(&plan, nx, ny, nz, CUFFT_Z2Z);


    
	Init_Conf();
   

    //cudaMemcpy(dfdc_d,dfdc,sizeof(cuDoubleComplex)*nx*ny*nz,cudaMemcpyHostToDevice);
    cudaMemcpy(dfdc_d,dfdc,sizeof(cuDoubleComplex)*nx*ny*nz,cudaMemcpyHostToDevice);
    
    Evolve();


    fclose(fpout);

    cufftDestroy(plan);
    cudaFree(comp_d);
    cudaFree(dfdc_d);
    
    return 0;
}
