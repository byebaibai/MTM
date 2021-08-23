
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <string.h>
#include <malloc.h>



#define PI 3.141592654
#define ABS(a) ((a) < (0) ? -(a) : (a))

#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(A, B) ((A) < (B) ? (A) : (B))

/* #include "nr.h" */
#include "nrutil.h"


void jrealft(float data[], unsigned long n, int isign);



void find_max_min(void *p, int n, float *max, float *min, int  flag)
 {
  int i,*pi=p;
  float *pf=p;
  
  if(flag){
    (*max) = (*min) = *pf;
    for(i=1; i<n; i++){
      *max = MAX(*max, *(pf+i));
      *min = MIN(*min, *(pf+i));
    }
  }
  else{
    (*max) = (*min) = *pi;
    for(i=1; i<n; i++){
      *max = MAX(*max, *(pi+i));
      *min = MIN(*min, *(pi+i));
    }
  }
}
/**************************************************************/

/** FUNC DEF **/ int get_pow_2(int inum)
{
	int             j, klength;
	/* find smallest power of 2 that encompasses the data */

	for (j = 1; pow((double) 2, (double) j) < inum; j++);
	return klength = pow((double) 2, (double) j);
}
/*******************************************************************/
/** FUNC DEF */ void  Db_scale(float *spec1, float *spec2, int num_freqs)
{
	/* make a copy of spec2 onto spec1 with the Db scale  */
	int             i;
	for (i = 0; i < num_freqs; i++) {
		if (spec2[i] <= 0.0) {
			fprintf(stderr, "negative or zero spectrum: %d\n", i);
			fprintf(stderr, "%g \n", spec2[i]);
			exit(0);
		}
		spec1[i] = 10. * log10((double) spec2[i]);
	}
}

/**FUNC DEF */ void  Log_scale(float *spec1, double *spec2, int num_freqs)
{
	/* make a copy of spec2 onto spec1 with the Db scale  */
	int             i;
	for (i = 0; i < num_freqs; i++) {
		if (spec2[i] <= 0.0) {spec1[i]=0.0;
			fprintf(stderr, "negative or zero spectrum: %d %g \n", i, spec2[i]);
			}
		spec1[i] =  log10( spec2[i]);
	}
}



/**FUNC DEF */ void Scale_Trace2(float *spec1, int num1, float *spec2, int num2, float *spec3)
{
	/*
	 * this routine scales one vector to the other for plotting purposes
	 */
	float           diff1, diff2, max1, min1, max2, min2;
	int             float_or_int = 1;
	int             i;
	find_max_min(spec1, num1, &max1, &min1, float_or_int);
	find_max_min(spec2, num2, &max2, &min2, float_or_int);
	diff1 = max2 - min2;
	diff2 = max1 - min1;
	/* fprintf(stderr, "Scale2: max1=%e, min1=%e ,max2=%e ,min2=%e\n", max1, min1, max2, min2); */
	for (i = 0; i < num2; i++) {
		spec3[i] = ((spec2[i] - min2) / diff1) * diff2 + min1;
	}

}


/*********************************************************/

/*********************************************************/
/** FUNC DEF */ double  scale_trace_RMS(float x[], int lx)
   {

     int k;
      double mean;
      double std;
      mean = 0.;
     if(lx < 2 ) return mean;
 
      for( k=0; k<lx; k++)
        {
        mean += x[k];
        }

       mean = mean/ (float)lx;

      for( k=0; k<lx; k++)
        {
        x[k] = x[k] - mean;
        }

         std = 0.;
      for( k=0; k<lx; k++)
        {
        std += x[k] * x[k];
        }
        std = sqrt(std)/(lx-1);

      for( k=0; k<lx; k++)
        {
        x[k] = x[k] / std ;
        }

    /*  fprintf(stderr," %e %e \n", mean, std);*/
      return  mean;
   }
/*********************************************************/

/** FUNC DEF */ double  remove_mean(float x[], int lx)
   {

     int k;
      double mean;

      mean = 0.;
     if(lx < 2 ) return mean;
 
      for( k=0; k<lx; k++)
        {
        mean += x[k];
        }

       mean = mean/ (float)lx;

      for( k=0; k<lx; k++)
        {
        x[k] = x[k] - mean;
        }

      /*fprintf(stderr,"MEAN = %20.10f\n",mean);*/
      return  mean;
   }
/*********************************************************/

/** FUNC DEF */void complex_array( float input[], int inlength,float  output[], int olength)
 
 {
    int i,k,j;

            for(i=0; i< inlength ;  i++)
       {
             k = 2*i;
             j = k+1;

           if(i>olength) break;
          output[k] = input[i];
           
          output[j] = 0.0;
         }

 }
/*********************************************************/
/** FUNC DEF */ void zero_pad(float  output[], int start , int olength)
{
    int i;
    for( i= start ; i< olength; i++) 
         {   
            output[i] = 0.0; 
               }
}


/*********************************************************/

/** FUNC DEF */void copy_trace(float a[], float b[], int num)
{
/* copies vector a ONTO vector b  */
     int i;
    for(i=0; i< num ;  i++)
       {
           b[i] = a[i];
         }
         
}

/*********************************************************/

/** FUNC DEF */void window_trace( float input[], float output[], int start, int num)
{
   int i;
    for(i=0; i< num ;  i++)
       {
           output[i] = input[i+start];
         }
}
/*********************************************************/

/** FUNC DEF */void print_array( float array[], int inum)
 {
     int i;
           for(i=0; i<inum; i++)
     {
         printf("%d %f\n", i, array[i]);
      }
}


/************************************************************/
/** FUNC DEF */void rm_lintrend(float *x,  float *y, int length, float a, float b)
   {   /* program removes a linear trend of the data in vector y
           returns a de-trended vector y  */
      int i;
       /* fprintf(stderr, "in rm_lintrend %f %f\n",a,b); */
      for (i = 0; i < length; i++)
      
        { y[i] = y[i] - x[i]*a - b; }

    /*  fprintf(stderr, "done in rm trend....\n"); */
     }
/************************************************************/
/** FUNC DEF */void get_abfit(float *x, float *y, int length, float *slope, float *intercept)
{   

      float s=0.0, sx=0.0, sy=0.0, sxx=0.0, sxy=0.0;
      float del;
      int i;
      
      for (i = 0; i < length; i++)
      
        {
        
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        sxy += x[i]*y[i];
        } 
      
       s = length;
       
       del = s*sxx - sx*sx;
       
       if(del != 0.0)
       {
            *intercept = (sxx*sy - sx*sxy)/del;
            *slope = (s*sxy - sx*sy)/del;
         }
       
}
/************************************************************/
/** FUNC DEF */void rm_lin_sig_trend(float *y, int n, float dt, float *slope, float *cept)
   {   /* program removes a linear trend of a time series vector y
           returns a de-trended vector y  */
      int i;
      float *x;
      float a, b;

   /* create an x vector of time */

     fprintf(stderr, "starting rm_lin_sig_trend....\n");

       x = (float *)malloc( ((n)*sizeof(float)) ); 
          for(i=0; i<n; i++)
	      { 
	      x[i] = i*dt;
              }

    /*  find the line of the data  */
     get_abfit(x, y, n, &a, &b);

/*   remove the trend from the data  */
     rm_lintrend(x, y, n, a, b);
     fprintf(stderr, "fixing slope and cept....\n");

#if 0
       *slope = a;
       *cept  = b;


     fprintf(stderr, "done in rm_lin_sig_trend %f %f \n", *slope, *cept);
#endif
      free(x);

     }

/***********************/
/** FUNC DEF */void get_indies(float t1,float t2,float dt,float tref, int numtot, int *ibeg, int *inum)
{


            *inum = (int)((t2 - t1)/dt )+1 ;  
               
            *ibeg = (int)((t1 - tref)/dt);

                     if( *ibeg < 0 ) *ibeg = 0;
             if( (*ibeg + *inum)  >  numtot ) *inum = numtot- *ibeg;

}


/***********************/
/** FUNC DEF */void get_indtim(float *t1,float *t2,float dt,float tref, int numtot, int ibeg, int inum)
{

           *t1 = (float)ibeg*dt+tref;
           *t2 = (inum-1)*dt+ *t1;

}

float get_taper(int itype,int n, int k, float percent)
{
/*
c-this function generates a single sample of a data window.
c-itype=1(rectangular), 2(tapered rectangular), 3(triangular),
c-      4(hanning), 5(hamming), or 6(blackman).
c-      (note:  tapered rectangular has cosine-tapered 10% ends.)
c-n=size (total no. samples) of window.
c-k=sample number within window, from 0 through n-1.
c-  (if k is outside this range, spwndo is set to 0.)
*/
int l;
float vwin;
vwin = 0.0;
if(itype < 1 || itype > 6) return vwin;
if(k<0 || k > n) return vwin;
vwin = 1.0;
switch(itype)
  {
     case 1:
        return vwin;
        break;
     case 2:
      l=(n-2)*percent;
      if(k<=l) vwin=0.5*(1.0-cos(k*PI/(l+1)));
      if(k>=n-l-2) vwin=0.5*(1.0-cos((n-k-1)*PI/(l+1)));
        return vwin;
        break;
      case 3:
          vwin=1.0-ABS(1.0-2*k/(n-1.0));
       return vwin;
        break;
      case 4:
          vwin=0.5*(1.0-cos(2*k*PI/(n-1))); 
       return vwin;    
        break;
      case 5:
          vwin=0.54-0.46*cos(2*k*PI/(n-1));
       return vwin;
        break;
      case 6:
          vwin=0.42-0.5*cos(2*k*PI/(n-1))+0.08*cos(4*k*PI/(n-1));
       return vwin;
        break;
     }
     return vwin;
}
/*********************************************************/
/** FUNC DEF */int  apply_taper(float x[],int lx,int itype,float tsv)
   {
     int k, ierror;
       float w;
      ierror=1;
      if(itype<1  || itype > 6) return ierror;
      tsv=0.0;
      for( k=0; k<=lx; k++)
        {
        w=get_taper(itype,lx+1,k, 0.05);
        x[k]=x[k]*w;
        tsv=tsv+w*w;
        }
      ierror=0;
      return ierror;
   }

float get_cos_taper(int n, int k, float percent)
{
/* n = number of data points in window. k = index of data point  */
int l;
float vwin;
vwin = 0.0;

if(k<0 || k > n) return vwin;
vwin = 1.0;

      l=(n-2)*percent;
      if(k<=l) vwin=0.5*(1.0-cos(k*PI/(l+1)));
      if(k>=n-l-2) vwin=0.5*(1.0-cos((n-k-1)*PI/(l+1)));


        return vwin;
}


/** FUNC DEF */void  smooth_fft(float *data, int npoints, float dt, float *naive_spec, int klen, float fWidth)
{
	int             isign = 1;
	float          *dtemp, df, freqwin, tem;
	int             num_freqs;
	int             i, k, j;
	num_freqs = 1 + klen / 2;

	dtemp = vector(0, klen);

	copy_trace(data, dtemp, npoints);


	zero_pad(dtemp, npoints, klen);
	jrealft(dtemp - 1, (unsigned long) klen, isign);




	for (i = 1; i < num_freqs - 1; i++) {
		naive_spec[i] = (SQR(dtemp[2 * i + 1]) + SQR(dtemp[2 * i]));

	}
	naive_spec[0] = SQR(fabs(dtemp[0]));
	naive_spec[num_freqs - 1] = SQR(fabs(dtemp[1]));



	df = 2 * (0.5 / dt) / klen;
	freqwin = (int) (fWidth / df) / 2;

#if 1

	/* smooth the periodogram   

	fprintf(stderr, "smooth the periodogram 4, freqwin=%d\n", freqwin);
         */


	for (i = 0; i < num_freqs; i++) {
		tem = 0.0;
		k = 0;
		for (j = i - freqwin; j <= i + freqwin; j++) {
			if (j > 0 && j < num_freqs - 1) {
				tem += naive_spec[j];
				k++;
			}
		}

		if (k > 0) {
			naive_spec[i] = tem / (float) k;
		} else
			naive_spec[i] = naive_spec[i];


	}

#endif
         free_vector(dtemp, 0 , klen);
}

