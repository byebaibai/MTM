#include "jl.h"




/*  prototypes  */

int             get_pow_2(int inum); 
float           get_cos_taper(int n, int k);

int adwait(double *sqr_spec,  double *dcf,
            double *el, int nwin, int num_freq, 
            double *ares, double *degf, double avar);
 
void get_F_values(double *sr, double *si, int nf, 
                      int nwin, float *Fvalue, double *b);

void do_mtap_spec(float *data, int npoints, int kind,
	     int nwin, float npi, int inorm, float dt, 
               float *ospec, float *dof, float *Fvalues, int klen);

int hires(double *sqr_spec,  double *el, int nwin, int num_freq, double *ares);

void jfour1(float data[], unsigned long nn, int isign);

void jrealft(float data[], unsigned long n, int isign);

int jtridib_(int *n, double *eps1, double *d, double *e, double *e2, 
        double *lb, double *ub, int *m11, int *m, 
        double *w, int *ind, int *ierr, 
	double *rv4, double *rv5);

int jtinvit_(int *nm, int *n, double *d, double *e, double *e2,
             int *m, double *w, int *ind, double *z, int *ierr, 
             double *rv1, double *rv2, 
	     double *rv3, double *rv4, double *rv6);

void  mt_get_spec(float *series, int inum, int klength, float *amp);

int  multitap(int n, int nwin, double *el,  float npi, double *tapers, double *tapsum);

void zero_pad(float  output[], int start , int olength);

double  remove_mean(float x[], int lx);

float get_cos_taper(int n, int k);

/***************************************************************/

void
main(int argc, char **argv)
{

	int             i, j, k, l;
	int             npoints, nwin, flag = 0;
	float           npi;
	float          *xt, vwin;
	float           rex[2000], rey[2000];
	FILE           *fopen(), *inf, *fp;
	int             logg, lspec;
	int             num_points;
	float          *data, *dtemp, dt, tem, *dtap, *tap_spec, *autoreg;
	int             iwin, kk;
	int             klen;
	float          *spec, *naive_spec;
	float          *dof, *Fvalues, xline[4][2], yline[4][2];
	char            in_file[100];
	/************/
	int             n1, n2, kind, num_freqs;
	int             inorm, K, freqwin;
	float           norm, fWidth;
	int             increase;
	float           f0, df, nyquist, *freq, *ex;
	int             isign = 1;
	double          mean;


	int             num_cof;



	kind = 1;
	inorm = 1;




	/*
	 * Data and Parameter I/O  We need to read in the time series, and
	 * the sampling interval, (which may be set to 1 if unimportant) The
	 * data in this case is arranged such that the first line is
	 * num_points = number of points in the time series. dt = sampling
	 * rate the next num_points floats are the actual time series.
	 * 
	 * The parameters required are the number of pi-prolate functions, npi
	 * and the number of summing windows, nwin
	 * 
	 */


	npi = 3.0;
	nwin = 5;
	kind = 1;
	inorm = 1;

	fprintf(stderr, "argc = %d\n", argc);

	if (argc < 5) {


		fprintf(stderr, "HELLO:\n");

		fprintf(stderr, "need three args: file npi nwin [ kind inorm ] \n");
		fprintf(stderr, "example: testmt file 3 5 1 1\n");
		fprintf(stderr, "kind = 1 : hires \n");
		fprintf(stderr, "kind = 2 : adwait \n");
		fprintf(stderr, "kind = 3 : naive periodogram \n");

		fprintf(stderr, "inorm = 1 : standard \n");
		fprintf(stderr, "inorm = 2 : other \n");
		exit(0);

	}
	strcpy(in_file, argv[1]);
	npi = atof(argv[2]);
	nwin = atoi(argv[3]);
	kind = atoi(argv[4]);
	inorm = atoi(argv[5]);
	



#if 0
	fprintf(stderr, "need three args: file npi nwin [ kind inorm ] \n");
	fprintf(stderr, "example: testmt file 3 5 1 1\n");
	fprintf(stderr, "kind = 1 : hires \n");
	fprintf(stderr, "kind = 2 : adwait \n");
	fprintf(stderr, "kind = 3 : naive periodogram \n");

	fprintf(stderr, "inorm = 1 : standard \n");
	fprintf(stderr, "inorm = 2 : other \n");




	fprintf(stderr, "\n\nType in the input file name, npi, nwin, kind and inorm:\n\n");
	scanf("%s", in_file);
	scanf("%f", &npi);
	scanf("%d", &nwin);
	scanf("%d", &kind);
	scanf("%d", &inorm);
#endif

	fprintf(stderr, "\n\nfilename=%s npi=%f nwin=%d, kind=%d inorm=%d\n\n",
		in_file, npi, nwin, kind, inorm);


	if ((inf = fopen(in_file, "r")) == NULL) {
		fprintf(stderr, "file not found\n");
		exit(0);
	}
	k = fscanf(inf, "%d %f", &num_points, &dt);


	/*
	 * p. 335, Percival and Walden, choose npi=2,3,4 some small integer W
	 * = npi/(num_points*dt); or num_points*W = npi/dt ;
	 * 
	 * K < 2*num_points*W*dt
	 * 
	 * nwin = 0...K-1
	 * 
	 */


	fWidth = npi / ((float) num_points * dt);
	K = (int) 2    *num_points * fWidth * dt;
	printf("fWidth = %f   K = %d \n", fWidth, K);



	nyquist = 0.5 / dt;


	klen = get_pow_2(num_points);


	increase=1;
	klen = klen * pow((double) 2, (double) increase);


	/* klen = 1024; */

	fprintf(stderr, " klen = %d num_points=%d \n", klen, num_points);

	num_freqs = 1 + klen / 2;








	/*
	 * READ IN THE TIME SERIES: floating point array
	 * 
	 */


	data = (float *) malloc(num_points * sizeof(float));



	i = 0;
	while ((k = fscanf(inf, "%f", &data[i])) > 0) {
		i++;
	}


	npoints = i;
	k = 1;
	fprintf(stderr, "done getting data...\n");





	fprintf(stderr, "INPUT: %d %d %d %f %d %f %d\n", npoints, kind, nwin, npi, inorm, dt, klen);


	mean = remove_mean(data, npoints);

	fprintf(stderr, " mean = %f \n", mean);



	/*----------------------------    do simple (naive) periodogram ------------ */



	naive_spec = (float *) malloc(klen * sizeof(float));
	dtemp = (float *) malloc(klen * sizeof(float));


	/* 10% cosine  taper */

	for (i = 0; i < num_points; i++) {
		vwin = get_cos_taper(num_points, i);




		dtemp[i] = vwin * data[i];

		/*
		 * if(i<10 || i > num_points-10) printf("%d %f %f %f\n", i,
		 * vwin, dtemp[i], data[i]);
		 */




	}
	norm = 1. / (num_points * num_points);

	zero_pad(dtemp, num_points, klen);
	jrealft(dtemp - 1, (unsigned long) klen, isign);




	for (i = 1; i < num_freqs - 1; i++) {
		naive_spec[i] = norm * (SQR(dtemp[2 * i + 1]) + SQR(dtemp[2 * i]));

	}
	naive_spec[0] = norm * SQR(fabs(dtemp[0]));
	naive_spec[num_freqs - 1] = norm * SQR(fabs(dtemp[1]));


	df = 2 * nyquist / klen;
	freqwin = (int) (fWidth / df);


	/* smooth the periodogram   */
	fprintf(stderr, "smooth the periodogram 4, freqwin=%d\n", freqwin);

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
			dtemp[i] = tem / (float) k;
		} else
			dtemp[i] = naive_spec[i];


	}





	/**********************************************/




	spec = (float *) malloc(klen * sizeof(float));
	dof = (float *) malloc(klen * sizeof(float));

	Fvalues = (float *) malloc(klen * sizeof(float));





	do_mtap_spec(data, npoints, kind, nwin, npi, inorm, dt, spec, dof, Fvalues, klen);
	fprintf(stderr, " done with do_mtap_spec\n");

	freq = (float *) malloc(klen * sizeof(float));




	inf = fopen("spec.out", "w");

	for (i = 0; i < num_freqs; i++) {
                 freq[i] = df*i;   
		fprintf(inf, "%d %g %g %g %g %g %g\n", i, freq[i],
			spec[i], naive_spec[i], dtemp[i],
			dof[i], Fvalues[i]);
	}


}
