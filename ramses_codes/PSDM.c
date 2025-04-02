//Compile as gcc PSDM.c -std=gnu99 -I/mnt/exports/data/install/fftw-3.3.8/include -L/mnt/exports/data/install/fftw-3.3.8/lib -lfftw3 -I${GSL_INCDIR} -L${GSL_LIBDIR} -lm -lgsl -lgslcblas -o psdm
//Run as ./psdm

#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <fftw3.h> 
#include <gsl/gsl_sf_trig.h>
#include <fcntl.h> // for open
#include <unistd.h> // for close
#define MAX_PATH_LENGTH 512

FILE *fp; 

float kn, kf, coeff; 

double fourier_space_grid_function(double k1, double k2, double k3, int nbox)
{

    double w; 
    double k1arg, k2arg, k3arg; 
    
    k1arg = k1/(2.0*(double)kn); 
    k2arg = k2/(2.0*(double)kn);
    k3arg = k3/(2.0*(double)kn);

    w = gsl_sf_sinc(k1arg)*gsl_sf_sinc(k2arg)*gsl_sf_sinc(k3arg); 
    w = w*w; 

    return w;
    
}

void myps(int nbox, float **grid, float box_size, int level, int outs)
{

    fftw_complex *in, *out; 
    fftw_plan p; 
    int N0, N1, N2, N;
    int i, j, k, sign, nhalf, iw, i2, i1, i3, nkweights, m; 
    unsigned flags;
    float g ; 
    double *w, *deltasqk, *powspec; 
    int *iweights;				
    double tpisq;
    float *local_grid; 
    long int N_gridpoints; 
    double shot_noise_correction, c1; 
    double error, tolerance, k1, k2, p1, p2, alpha0, *powspec_dummy, *c2, g_local; 
    int i_dummy, i_c2, ii, jj, kk, i_local, j_local, k_local;
    int counter; 
    double old_alpha0, contrib, sum;

    local_grid = *grid;

    N0 = nbox;
    N1 = nbox; 
    N2 = nbox; 

    nhalf = nbox / 2;
    N = pow(nbox, 3);
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N); 
    out = in; 
    // out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N); 
    w = (double *) malloc(nbox*sizeof(double));
    p = fftw_plan_dft_3d(N0, N1, N2, in, out, FFTW_FORWARD, FFTW_MEASURE); 

    // ----------------------

    kf = 2*M_PI/box_size; // h/Mpc 
    kn = M_PI*(float)nbox/box_size; // h/Mpc; Nyquist wavenumber
    printf("kf=%f\n",kf);
    printf("kn=%f\n",kn);
    coeff = pow(box_size/(2.0*M_PI),2.0); // Mpc/h 
    printf("coeff=%f\n",coeff);
    for(i=0; i<nbox; i++)
    {
	if (i > nhalf) 
	    iw = i-nbox; 
	else 
	    iw = i; 
	w[i] = kf*(double)iw;
    }
    printf("kmax=%f\n", kf*(double)nhalf);
  
    // ----------------------

    for(i=0; i<N; i++)
    {
	in[i][0] = local_grid[i];
	in[i][1] = 0.0;
    }

    fftw_execute(p); 

    // ---------------------
  
    contrib = 0.0;
    powspec = (double *) calloc(nbox, sizeof(double));
    for (i=0; i<nbox; i++) 
	for (j=0; j<nbox; j++) 
	    for (k=0; k<nbox; k++)
	    {

		g = w[i]*w[i] + w[j]*w[j] + w[k]*w[k]; 

		if (g != 0)
		{
		    i1 = (int)(0.5+sqrt(g*coeff));
		    contrib = pow(out[k+nbox*(j+nbox*i)][0],2.0) + 
			pow(out[k+nbox*(j+nbox*i)][1],2.0);
		    contrib /= pow(fourier_space_grid_function(w[i],w[j],w[k],nbox),2.0);
		    powspec[i1] += contrib; 
		}
	    }

    fftw_free(in); 
    //fftw_free(out); 

    printf("Power calculation done\n");

    // ----------------------

    iweights=(int *) malloc(nbox*sizeof(int));
    for(i=0;i<nbox;i++)
	iweights[i]=0;
    tpisq=2.0*M_PI*M_PI;

    for(i=0;i<nbox;i++)
    {
	i1=i;
	if(i1>=nhalf)
	    i1=nbox-i1;
	for(j=0;j<nbox;j++)
	{
	    i2=j;
	    if(i2>=nhalf)
		i2=nbox-i2;
	    for(k=0;k<nbox;k++)
	    {
		i3=k;
		if(i3>=nhalf)
		    i3=nbox-i3;
		m=0.5+sqrt(i1*i1+i2*i2+i3*i3);
		iweights[m]+=1;
	    }//for...k
	}//for...j
    }//for...i

    // ----------------------
   printf("\nCalculating Delta^2 ...\n");
   deltasqk = (double *) calloc(nbox, sizeof(double));

   for (i = 0; i < nhalf; i++)
   {
     powspec[i] = powspec[i]*pow(box_size,3.0)/pow((float)nbox,6.0);
     powspec[i] /= (double)iweights[i];
     deltasqk[i] = pow(w[i],3.0)*powspec[i]/(2.0*M_PI*M_PI);
   }
  
  // ----------------------
  printf("Skipping shot noise correction ...\n");
  /* Correct shot noise effect.  See Equations 19 and 21 of Jing (2005
   * ApJ 620 559). 

   N_gridpoints = nbox*nbox*nbox; 
   for (i = 0; i < nhalf; i++)
   {
       c1 = 1.0 - (2.0/3.0)*pow(sin(M_PI*w[i]/(2*kn)),2);
       shot_noise_correction = c1/((double)N_gridpoints);
       powspec[i] -= shot_noise_correction; 
       deltasqk[i] = pow(w[i],3.0)*powspec[i]/(2.0*M_PI*M_PI);
   }*/

  // ----------------------
  
   free(iweights); 

  // ----------------------
   //fp is the file in which the output will be stored.
   char save_path[MAX_PATH_LENGTH];
    snprintf(save_path, MAX_PATH_LENGTH, "/user1/poojarani/Lya_Comparison/ramses_analysis/lev%03d_len%d/PS_DM/psdm_%03d.txt", level, (int)box_size, outs);
    fp = fopen(save_path, "w");

   for (i = 0; i < nhalf; i++) 
     fprintf(fp, "%f  %f  %f\n", w[i], powspec[i], deltasqk[i]); 
   fclose(fp);

   /*sum = 0.0;
   for (i = 0; i < nhalf-1; i++)
 	  sum += (w[i+1]-w[i])*deltasqk[i]/((w[i]+w[i+1])/2.0);
  
   printf("Variance from power spectrum: %e\n",sum);*/
   
   printf("PS written to file.\n");
   fftw_destroy_plan(p); 
   free(w); 
  
   return; 
  
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <lev> <length> <output>\n", argv[0]);
        return EXIT_FAILURE;
    }
    int lev = atoi(argv[1]);
    float boxsize = atof(argv[2]); //boxsize in Mpc/h
    int output = atoi(argv[3]);

   int ncell = pow(2,lev);		//no.of cells in each direction
   int tcells = ncell*ncell*ncell;	//total no.of cells
   //double arr[tcells];		//1st read the 21-cm data in double
   // Adding lines 
    double *arr = (double *)malloc(tcells * sizeof(double));
    if (arr == NULL) {
        perror("Memory allocation failed for arr");
        return EXIT_FAILURE;
    }
   
   //float del21[tcells];		//store the same data in float
   float *del21;
   del21 = (float *) malloc (tcells * sizeof (float));
   FILE *ptr;

   //Enter here the path to the binary file which was generated by npy2bin.py code.
   char file_path[MAX_PATH_LENGTH];
    snprintf(file_path, MAX_PATH_LENGTH, "/user1/poojarani/Lya_Comparison/ramses_analysis/lev%03d_len%d/DM_grid/dm_grid_%03d.bin", lev, (int)boxsize, output);
    ptr = fopen(file_path, "rb");

   if (ptr == NULL) {
    perror("Error opening file");
    return EXIT_FAILURE;
    }
   fread(arr, sizeof(double), tcells, ptr);
   fclose(ptr);
   
   //As the function 'myps' accepts float array, here we convert elements of (double)arr to float
   for(int i=0;i<tcells;i++)
   {
	del21[i] = (float) arr[i];
   }
   
   myps(ncell,&del21,boxsize,lev,output);	//This is the main command where the PS is calculated.

   /*for(int i=tcells-10;i<tcells;i++)
   {
	printf("%lf\n",del21[i]);
   }*/

}