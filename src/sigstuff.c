#include "jl.h"

int
get_pow_2(int inum)
{
	int             j, klength;
	/* find smallest power of 2 that encompasses the data */

	for (j = 1; pow((double) 2, (double) j) < inum; j++);
	return klength = pow((double) 2, (double) j);
}


double  remove_mean(float x[], int lx)
   {

     int k;
      double mean;
      mean = 0.;
     if(lx < 2 ) return mean;
 
      for( k=0; k<=lx; k++)
        {
        mean = x[k]+mean;
        }

       mean = mean/ (float)lx;

      for( k=0; k<=lx; k++)
        {
        x[k] = x[k] - mean;
        }

      
      return  mean;
   }
/*********************************************************/
void zero_pad(float  output[], int start , int olength)
{
    int i;
    for( i= start ; i< olength; i++) 
         {   
            output[i] = 0.0; 
               }
}

float get_cos_taper(int n, int k)
{
int l;
float vwin;
vwin = 0.0;

if(k<0 || k > n) return vwin;
vwin = 1.0;

      l=(n-2)/10;
      if(k<=l) vwin=0.5*(1.0-cos(k*PI/(l+1)));
      if(k>=n-l-2) vwin=0.5*(1.0-cos((n-k-1)*PI/(l+1)));


        return vwin;
}
