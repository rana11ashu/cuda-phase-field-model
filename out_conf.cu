void Output_Conf(int step){

    FILE  *fp;
    cuDoubleComplex *realComp;
    char  fname[100];
    
    sprintf(fname,"prof_%d",step);
    fp  = fopen(fname,"wb");

    realComp  = (cuDoubleComplex*)malloc((sizeof (cuDoubleComplex)*nx*ny*nz));
    cudaMemcpy(realComp,dfdc_d,sizeof(cuDoubleComplex)*nx*ny*nz,cudaMemcpyDeviceToHost);
    
    int i,j,k;
    for (i = 0; i < nx; i++){
        for (j = 0; j < ny; j++){
          for (k = 0; k < nz; k++){
          	
            fprintf(fp,"%d\t%d\t%d\t%f\n",i,j,k,realComp[k+nz*(j+i*ny)].x );}
        }
      }
      


/*

    for (i = 0; i < nx; i++){
        for (j = 0; j < ny; j++){
          for (k = 0; k < nz; k++){
  
            fprintf(fp,"%d\t%d\t%d\t%e\n",
                    i,j,k,realComp[k+nz*(j+i*ny)]);}
        }
      }*/
   
    //fwrite (realComp,sizeof(double)*nx*ny*nz,1,fp);
    fclose(fp);
    free(realComp);
}

/*
void Output_Conf(int steps)
{
    FILE *fpt;

    char fn[100];

    sprintf(fn, "conf.%06d", steps);

    fpt = fopen(fn, "w");
    fwrite(&dfdc[0], sizeof(double), 2 * nx * ny * nz, fpt);
    fclose(fpt);

    sprintf(fn, "prof_gp.%06d", steps);
    fpt = fopen(fn, "w");

    for (int k = 0; k < nz; k++) {
	  for (int j = 0; j < ny; j++) {
	    for (int i = 0; i < nx; i++) {
		fprintf(fpt, "%d\t%d\t%d\t%lf\n", i, j, k, dfdc[i + nx * (j + ny * k)].x );
	    }
	    fprintf(fpt, "\n");
	  }
	  fprintf(fpt, "\n");
    }
    fclose(fpt);
}*/
