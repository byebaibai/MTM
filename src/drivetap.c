
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"


#define PI 3.14159265358979

#define ABS(a) ((a) < (0) ? -(a) : (a))

#define NRANSI

int getline(FILE *input,char s[],int lim);

int get_pow_2(int inum);

float get_cos_taper(int n, int k, float percent);
 
void do_mtap_spec(float *data, int npoints, int kind,
	    int nwin, float npi, int inorm, float dt, float *ospec, float *dof, float *Fvalues, int klen);


double  remove_mean(float x[], int lx);

void zero_pad(float  output[], int start , int olength);

void jrealft(float data[], unsigned long n, int isign);



void main(int argc, char **argv)
{

  int             i, j, k,l;
  int             npoints, nwin,  flag2, kdof;
  float           npi,flag1;
  float          *xt, vwin;
  float rex[2000], rey[2000];
  FILE           *fopen(), *inf, *fp, *record, *srecord;
  int             logg,lspec;
  int             num_points;
  float          *data, *dtemp , dt, tem, *dtap, *tap_spec, *autoreg;
  int iwin, kk;
  int klen;
  float          *spec, *naive_spec;
  float           *dof, *Fvalues, xline[4][2], yline[4][2];
  char  in_file[100];
 

  /************/
  int n1,n2, kind, num_freqs;
  int inorm, K , freqwin;
  float norm, fWidth, anrm;
  int increase;
  float f0, df, frq1, nyquist, *freq, *ex;
  int isign=1;
  double mean;
  float t1, t2;
  char  name[100];

  char line[101];
  int num_cof;


/* f-values by row for 
     50%, 90%, 95%, 99%   confidence with 2,K dof
   from abramowitz and stegan, p. 987-988*/

  float fvals[4][7]={ { 
      1.,0.828,.780,.757,.743,.735,.729 } ,
		      { 9.,4.32,3.46,3.11,2.92,2.81,2.73},
		      {19.0,6.94,5.14,4.46,4.1,3.89,3.74},
		      {99.,18.,10.92,8.65,7.56,6.93,6.51} };

  kind=1;
  inorm=1;




  /******* Data and Parameter I/O  
	   We need to read in the time series, and the sampling interval,
	   (which may be set to 1 if unimportant)
	   The data in this case is arranged such that the first
	   line is 
	   num_points = number of points in the time series.
	   dt = sampling rate
	   the next num_points floats are the actual time series.

	   The parameters required are 
	   the number of pi-prolate functions, npi
	   and the number of summing windows, nwin




  ******/


  npi = 3.0;
  nwin = 5;
  kind = 1;
  inorm = 1;
  
  fprintf(stderr,"argc = %d\n",argc);
  
  if(argc < 6) {


    fprintf(stderr,"HELLO:\n");
    
    fprintf(stderr,"need three args: file npi nwin kind inorm \n");
    fprintf(stderr,"example: testmt file 3 5 1 1\n");
    fprintf(stderr,"kind = 1 : hires \n");
    fprintf(stderr,"kind = 2 : adwait \n");
    fprintf(stderr,"kind = 3 : naive periodogram \n");
    
    fprintf(stderr,"inorm = 1 : standard \n");
    fprintf(stderr,"inorm = 2 : other \n");
    exit(0);
    
  }


  for(i=0; i<argc; i++) {   fprintf(stderr,"%d %s\n", i, argv[i]); }
  
  strcpy(in_file ,argv[1]);
  npi = atof(argv[2]);
  nwin = atoi(argv[3]);
  kind = atoi(argv[4]);
  inorm = atoi(argv[5]);

     


  fprintf(stderr,"\n\nfilename=%s npi=%f nwin=%d, kind=%d inorm=%d num_cof=%d\n\n",
	  in_file,npi,nwin,kind,inorm, num_cof);
  
  
  if ((inf = fopen(in_file, "r")) == NULL) {
    fprintf(stderr, "file not found\n");
    exit(0);
  }
  
  
  /*k = fscanf(inf, "%s %d %f %f %f", &name, &num_points, &dt, &t1, &t2);*/
  
  k = getline(inf, line, 100);
  
  k = sscanf(line, "%d %f", &num_points, &dt);
  
  
  
  /*  
      p. 335, Percival and Walden, choose npi=2,3,4 some small integer
      W = npi/(num_points*dt);
   or num_points*W = npi/dt ;

   K < 2*num_points*W*dt

   nwin = 0...K-1

*/


         fWidth =  npi/((float)num_points*dt);
         K = (int) 2*num_points*fWidth*dt;
        printf("fWidth = %f   K = %d \n", fWidth, K);



        nyquist = 0.5/dt;
        
	increase = 0;
	klen = get_pow_2(num_points);

#if 0

    fprintf(stderr," klen = %d , want to increase it?\n",klen);
    fprintf(stderr," type in power of increase: 0=none, 1 = double, 2=2^2...etc\n");
     scanf("%d", &increase);
           klen = klen*pow( (double) 2, (double) increase);
#endif


           /*klen = 1024;*/

           fprintf(stderr," klen = %d num_points=%d \n",klen,num_points );

     num_freqs = 1+klen/2;

	       data = (float *) malloc(num_points * sizeof(float));
	
		ex = (float *) malloc(num_points * sizeof(float));

	i = 0;
	while ((k = fscanf(inf, "%f", &data[i])) > 0)
	{
         ex[i] = i*dt;
	i++;
	}
	
	if (i != num_points) {
		fprintf(stderr, "wrong number of data points i = %d num_points = %d\n",i,  num_points);
		exit(0);
	}
	npoints = num_points;
	k = 1;
fprintf(stderr, "done getting data...\n");
	

/*  exit(0); */


    mean = remove_mean(data, npoints); 

  


/*----------------------------    do simple (naive) periodogram ------------ */
    


    naive_spec=vector( (long)0, (long) num_freqs );

    dtemp=vector( (long)0, (long) klen );

          /* 10% cosine  taper */

	for (i = 0; i < num_points; i++) {
		vwin = get_cos_taper(num_points, i, .05); 
	      dtemp[i] = vwin*data[i];

            /*   if(i<10 || i > num_points-10) printf("%d %f %f %f\n", i, vwin, dtemp[i], data[i]);*/
	}

anrm = num_points;

#if 1

           	switch (inorm) {
	case 0:
		anrm = 1.;
		break;

	case 1:
		anrm = num_points;
		break;
	case 2:
		anrm = 1 / dt;
		break;
	case 3:
		anrm = sqrt((double) num_points);
		break;
	default:
		anrm = 1.;
		break;
	}
#endif

           norm = 1./(anrm*anrm);

      fprintf(stderr,"NORM = %g  inorm= %d\n",norm,inorm);
          zero_pad(dtemp, num_points, klen);
          jrealft(dtemp-1, (unsigned long) klen, isign);




            for(i=1; i<num_freqs-1; i++){
            naive_spec[i] =    norm*(SQR(dtemp[2*i+1])+SQR(dtemp[2*i]));

            }
          naive_spec[0] = norm*SQR(fabs(dtemp[0]));
          naive_spec[num_freqs-1] = norm*SQR(fabs(dtemp[1]));



         df = 2*nyquist/klen;
         freqwin = (int) ( fWidth/df) /2 ;      

#if 1
          /* smooth the periodogram   */
        fprintf(stderr, "smooth the periodogram 4, freqwin=%d\n", freqwin);

	for (i = 0; i < num_freqs ; i++) {
		tem = 0.0;
		k = 0;
		for (j = i - freqwin; j <= i + freqwin; j++) {
			if (j > 0 && j < num_freqs - 1) {
				tem += naive_spec[j];
				k++;
			}
		}
             
       if(k>0) { dtemp[i] = tem/(float)k; } else  dtemp[i] =naive_spec[i];


	}

	/*for (i = 1; i < num_freqs - 1; i++) naive_spec[i] = dtemp[i];*/
#endif


	for (i = 0; i < num_freqs ; i++) {

     if(  naive_spec[i]<0.0 || dtemp[i] < 0.0){
            fprintf(stderr,"negative or zero spectrum: %d\n",i);
            fprintf(stderr,"%g  %g\n", naive_spec[i], dtemp[i]);
          exit(0);
           }
	
		naive_spec[i] = 10.*log10(naive_spec[i]);
                dtemp[i] = 10.*log10(dtemp[i]);

		}

 /**********************************************/

      spec=vector( (long)0, (long) klen );
      dof=vector( (long)0, (long) klen );

      Fvalues=vector( (long)0, (long) klen );
  

      do_mtap_spec(data, npoints, kind,  nwin,  npi, inorm, dt, spec, dof, Fvalues, klen);
      fprintf(stderr, " done with do_mtap_spec: %d\n",  num_freqs);

      inf = fopen("spec.out", "w");


      fprintf(stderr, " writing to file:\n");


      for (i = 0; i < num_freqs; i++) {


	frq1 =  df*i;
/* 	fprintf(stderr, "%d %f\n",i, frq1); */



	
	fprintf(stderr, "%d %f %f %f %f %f %f\n", i,frq1,
		spec[i], naive_spec[i], dtemp[i],
		dof[i], Fvalues[i]);
	




 	fprintf(inf, "%d %f %f %f %f %f %f\n", i,frq1, 
 		spec[i], naive_spec[i], dtemp[i], 
 		dof[i], Fvalues[i]); 


      }
      fclose(inf);


      /*
      free_vector(ex,0,npoints);
      free_vector(freq,0,num_freqs);
      free_vector(spec,0,klen);

      free_vector( dof,0,klen);
      free_vector(Fvalues,0,klen);

      */


     
}




