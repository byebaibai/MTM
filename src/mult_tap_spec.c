
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979
#define ABS(a) ((a) < (0) ? -(a) : (a))

#define NRANSI
/* #include "jl.h" */
#include "nrutil.h"
/* #include "nr.h" */


#define perr(x,y)  (fprintf(stderr, x , y))
#define prbl (fprintf(stderr,"\n"))

void four1(float data[], unsigned long nn, int isign);

void zero_pad(float  output[], int start , int olength);
int 
adwait(double *tapers,  double *dcf,
       double *el, int nwin, int nf, double *ares, double *degf, double avar);
int
hires(double *spreal,  double *el, int nwin, int nf, double *ares);

int  multitap(int n, int nwin, double *el,  float npi, double *tapers, double *tapsum);

void  get_F_values(double *sr, double *si, int nf, int nwin,float *Fvalue, double *b);

void jrealft(float data[], unsigned long n, int isign);



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** FUNC DEF */ void  mt_get_spec(float *series, int inum, int klength, float *amp)
{
/*    series = input time series
      inum   = length of time series
      klength = number of elements in power spectrum (a power of 2)
      amp = returned power spectrum
*/

	int             i, j, isign = 1;

	unsigned long   nn;
	float           tsv;


    	nn = klength;

        

	/* copy amp onto series and apply zero padding to  klength */
       
	for (i = 0; i < inum; i++) {
		
		amp[i] = series[i];
 	
	}

       zero_pad(amp, inum, klength);


/*  Fast Fourier Transform Routine:  here we are using the Numerical Recipes
     routine jrealft which returns the fft in the 1-D input array
     packed as pairs of real numbers.
     The jrealft routine requires the input array to start at index=1
      so we must decrement the index of amp
*/


 
          jrealft(amp-1, nn, isign);


}

/** FUNC DEF */  void  do_mtap_spec(float *data, int npoints, int kind,
	    int nwin, float npi, int inorm, float dt, float *ospec, float *dof, float *Fvalues, int klen)
{
/*
     data = floating point input time series
    npoints = number of points in data
     kind = flag for choosing hires or adaptive weighting coefficients
     nwin = number of taper windows to calculate
     npi = order of the slepian functions
     inorm = flag for choice of normalization
     dt = sampling interval (time)
     ospec = output spctrum
    dof = degrees of freedom at each frequency
    Fvalues = Ftest value at each frequency estimate
    klen = number of frequecies calculated (power of 2)

*/

	int             i, j, k;
	double         *lambda, *tapers;
	long            len, longlen;
	float          *xt;
	FILE           *fopen(), *inf, *tap_file;
        FILE            *dof_file;

	int             logg;
	int             nn;
	float          *b;
	int             iwin, kk;

	/*************/
	double          anrm, norm;
        double            *ReSpec, *ImSpec;
	double         *sqr_spec,  *amu;
	float          *amp, *fv;
	double          avamp, temp, sqramp;
	double          sum, *tapsum;
	/************/
         int num_freqs;
         int len_taps, num_freq_tap;
   
	double         *dcf, *degf, avar;
	int             n1, n2, kf;
	int             flag;
	int             one = 1;

         double tem1, tem2;

/* lambda = vector of eigenvalues   
   tapsum = sum of each taper, saved for use in adaptive weighting  
   tapers =  matrix of slepian tapers, packed in a 1D double array    
*/


/*fprintf(stderr, "From do_mtap_spec:  npoints=%d kind=%d nwin=%d npi=%f inorm=%d dt=%f klen=%d\n",npoints, kind, nwin, npi, inorm, dt, klen );
               for(i=0; i< npoints; i++)
               fprintf(stderr, "%d %e\n",i, data[i]);
*/

	lambda = dvector((long)0, (long)nwin);
        tapsum=dvector((long)0,(long)nwin);
	len_taps = npoints * nwin;
	tapers = dvector((long)0,(long) len_taps);

             num_freqs = 1+klen/2;
             num_freq_tap = num_freqs*nwin;



       /* get a slepian taper  */

	k = multitap(npoints, nwin, lambda,  npi, tapers, tapsum);
#if 0

        tap_file = fopen("taper_file", "w");
        /* print out tapers for curiosity  */
        for(i=0; i<npoints; i++){
         for(j=0; j<nwin; j++)fprintf(tap_file,"%15.10f ",tapers[i+j*npoints]);
          fprintf(tap_file,"\n");
          }
          fclose(tap_file);
#endif



	
     /* choose normalization based on inorm flag  */

	anrm = 1.;

	switch (inorm) {
	case 0:
		anrm = 1.;
		break;

	case 1:
		anrm = npoints;
		break;
	case 2:
		anrm = 1 / dt;
		break;
	case 3:
		anrm = sqrt((double) npoints);
		break;
	default:
		anrm = 1.;
		break;
	}

	
	/* apply the taper in the loop.  do this nwin times  */

	b = vector((long)0, (long)npoints);
	amu = dvector((long)0,(long) num_freqs);
	sqr_spec = dvector((long)0,(long) num_freq_tap);
        ReSpec = dvector((long)0,(long) num_freq_tap);
        ImSpec = dvector((long)0,(long) num_freq_tap);


	for (iwin = 0; iwin < nwin; iwin++) {
		kk = iwin * npoints;
                kf = iwin * num_freqs;

		for (j = 0; j < npoints; j++)
			b[j] = data[j] * tapers[kk + j];   /*  application of  iwin-th taper   */
		
		amp = vector((long)0,(long) klen);

		mt_get_spec(b, npoints, klen, amp);  /* calculate the eigenspectrum */

	       
          
		
		sum = 0.0;


/* get spectrum from real fourier transform    */


          norm = 1.0/(anrm*anrm);

     



            for(i=1; i<num_freqs-1; i++){
       if(2*i+1 > klen) fprintf(stderr,"error in index\n");
       if(i+kf > num_freq_tap ) fprintf(stderr,"error in index\n");

            sqramp = SQR(amp[2*i+1])+SQR(amp[2*i]);

            ReSpec[i+kf] = amp[2*i];
            ImSpec[i+kf] = amp[2*i+1];



            sqr_spec[i+kf] =    norm*(sqramp);

             sum += sqr_spec[i+kf];
            }
          sqr_spec[0+kf] = norm*SQR(fabs(amp[0]));
          sqr_spec[num_freqs-1+kf] = norm*SQR(fabs(amp[1]));

            ReSpec[0+kf] = amp[0];
            ImSpec[0+kf] = 0.0;

            ReSpec[num_freqs-1+kf] = amp[1];
            ImSpec[num_freqs-1+kf] = 0.0;

             sum += sqr_spec[0+kf] + sqr_spec[num_freqs-1+kf];

        if(num_freqs-1+kf>num_freq_tap )fprintf(stderr,"error in index\n");

		temp = sum / (double) num_freqs;
		if (temp > 0.0)
			avamp = sqrt(temp) / anrm;
		else {
			avamp = 0.0;
			 /* fprintf(stderr," avamp = 0.0! \n"); */ 
		}


		free_vector(amp,(long) 0,(long) klen);

	}

                 free_vector(b, (long)0, (long)npoints);
		fv = vector((long)0,(long) num_freqs);

        /* choice of hi-res or adaptive weighting for spectra    */


#if 0
		if ((inf = fopen("mspec.file", "w")) == NULL) {
			fprintf(stderr, "mspec.file unable to open\n");
			return;
		}
	
		for (i = 0; i < num_freqs; i++) {
	for (iwin = 0; iwin < nwin; iwin++) {
	
		kf = iwin * num_freqs;

			fprintf(inf, "%f %f ", ReSpec[i + kf], ImSpec[i + kf]);
		}
             	fprintf(inf, "\n");
	}
	
	fclose(inf);

#endif 



	switch (kind) {
	case 1:
	
		hires(sqr_spec,  lambda, nwin, num_freqs, amu);
	        get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);
 

 
           for (i = 0; i < num_freqs; i++) {
		ospec[i] =amu[i];
                 dof[i] = nwin-1;
                 Fvalues[i] = fv[i];
                  }

		break;

	case 2:


		/* get avar = variance*/

		n1 = 0;
		n2 = npoints;


                  avar = 0.0;

		for (i = n1; i < n2; i++)
			avar += (data[i]) * (data[i]);


		switch (inorm) {
		case 0:
			avar = avar / npoints;
			break;

		case 1:
			avar = avar / (npoints * npoints);
			break;

		case 2:
			avar = avar * dt * dt;
			break;

		case 3:

			avar = avar / npoints;
			break;

		default:
			break;
		}

		 


		dcf = dvector((long)0,(long) num_freq_tap);
		degf = dvector((long)0,(long) num_freqs);

	


		adwait(sqr_spec, dcf, lambda, nwin, num_freqs, amu, degf, avar);

                get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);

#if 0
           /* dump out the degrees of freedom to a file for later inspection  */
              	if ((dof_file = fopen("dof_file", "w")) == NULL) {
		fprintf(stderr, "dof unable to open\n");
		return;}
               	for (i = 0; i < num_freqs; i++) {
	                fprintf(dof_file,"%f\n",degf[i]);
                       	}

                   fclose(dof_file);
#endif

                 /* rap up   */

           for (i = 0; i < num_freqs; i++) {
		ospec[i] =amu[i];
                 dof[i] = degf[i];
                 Fvalues[i] = fv[i];
                  }


                
		free_dvector(dcf,(long)0,(long) num_freq_tap);
		free_dvector(degf,(long)0,(long) num_freqs);
		free_vector(fv,(long)0,(long) num_freqs);


		break;
	}

/*  free up memory and return  */

        free_dvector(amu,(long)0,(long) num_freqs);
        

	free_dvector(sqr_spec, (long)0,(long) num_freq_tap);
	free_dvector(ReSpec, (long)0,(long) num_freq_tap);

	free_dvector(ImSpec, (long)0,(long) num_freq_tap);

	free_dvector(lambda,(long) 0,(long) nwin);

	free_dvector(tapers,(long) 0, (long)len_taps);
	free_dvector(tapsum,(long) 0, (long)nwin);

}
