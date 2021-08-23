
#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979
#define ABS(a) ((a) < (0) ? -(a) : (a))

/* #include "jl.h" */

#define   NRANSI
#include "nrutil.h"
/* #include "nr.h" */



#define DIAG1 0
#define MAX(a,b) ((a) >= (b) ? (a) : (b))




/** FUNC DEF */ int adwait(double *sqr_spec,  double *dcf,
            double *el, int nwin, int num_freq, double *ares, double *degf, double avar)
{
/*
c  this version uses thomson's algorithm for calculating 
c  the adaptive spectrum estimate
*/
      double as,das,tol,a1,scale,ax,fn,fx;
      double *spw, *bias;
      double test_tol, dif;
      int jitter, i,j,k, kpoint, jloop;
       float df;
       /*
c  set tolerance for iterative scheme exit */


#if 0
  fprintf(stderr,"test input\n adwait: %d %d %f\n",nwin, num_freq, avar);
       fprintf(stderr,"\n Data=\n");
    for( i =0; i<num_freq; i++)
       {
            fprintf(stderr,"%d %f \n",i,sqr_spec[i]);
          }
#endif


      tol=3.0e-4;
      jitter=0;
      scale=avar;
                 /***********************************
                 c  we scale the bias by the total variance of the frequency transform
                 c  from zero freq to the nyquist
                 c  in this application we scale the eigenspectra by the bias in order to avoid
                 c  possible floating point overflow
                 ************************************/
      spw=dvector(0,nwin);
      bias=dvector(0,nwin);
      for( i=0;i<nwin; i++)
          {
            
            bias[i]=(1.00-el[i]);
             }

    /*
         for( i=1;i<=nwin; i++) fprintf(stderr,"%f %f\n",el[i], bias[i]);
         fprintf(stderr,"\n"); 
     */

       /* START do 100 */
    for( jloop=0; jloop<num_freq; jloop++)
    {   
        
       for( i=0;i<nwin; i++)
         {  kpoint=jloop+i*num_freq;
            spw[i]=(sqr_spec[kpoint])/scale ;
           }
                          /********************************************
                          c  first guess is the average of the two 
                              lowest-order eigenspectral estimates
                           ********************************************/
       as=(spw[0]+spw[1])/2.00;

                              /* START do 300 */
                              /* c  find coefficients */

        for( k=0; k<20 ; k++) 
        {
          fn=0.00;
          fx=0.00;

          for( i=0;i<nwin; i++)
           {
               a1=sqrt(el[i])*as/(el[i]*as+bias[i]);
               a1=a1*a1;
               fn=fn+a1*spw[i];
               fx=fx+a1;
           }
  

         ax=fn/fx;
         dif = ax-as;
         das=ABS(dif);
      /* fprintf(stderr,"adwait: jloop = %d k=%d %g %g %g %g\n",jloop,k, fn,fx,ax,das);*/
         test_tol = das/as;
         if( test_tol < tol )
            { 
                   break;
               }

         as=ax;
        }

        /* fprintf(stderr,"adwait: k=%d test_tol=%f\n",k, test_tol);*/
                            /* end  300  */

                           /* c  flag if iteration does not converge */

      if(k>=20)  jitter++;
  
       ares[jloop]=as*scale;
                            /* c   calculate degrees of freedom */
      df=0.0;
      for( i=0;i< nwin; i++)
       {
          kpoint=jloop+i*num_freq;
          dcf[kpoint]=sqrt(el[i])*as/(el[i]*as+bias[i]);
          df=df+dcf[kpoint]*dcf[kpoint];
       }
 			/*
			 * we normalize degrees of freedom by the weight of
			 * the first eigenspectrum this way we never have
			 * fewer than two degrees of freedom
			 */

       degf[jloop]=df*2./(dcf[jloop]*dcf[jloop]);

  }                                       /* end 100 */

     /*fprintf(stderr,"%d failed iterations\n",jitter);*/
      free_dvector(spw,0,nwin);
      free_dvector(bias,0,nwin);

     return jitter;
}
