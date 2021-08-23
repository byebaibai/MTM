#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979
#define ABS(a) ((a) < (0) ? (-a) : (a))
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define DIAG1 0

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
#include "jl.h"

/* prototypes  */

int             get_pow_2(int inum);
float           get_cos_taper(int n, int k);

int 
adwait(double *sqr_spec, double *dcf,
       double *el, int nwin, int num_freq,
       double *ares, double *degf, double avar);

void 
get_F_values(double *sr, double *si, int nf,
	     int nwin, float *Fvalue, double *b);

void 
do_mtap_spec(float *data, int npoints, int kind,
	     int nwin, float npi, int inorm, float dt,
	     float *ospec, float *dof, float *Fvalues, int klen);

int             hires(double *sqr_spec, double *el, int nwin, int num_freq, double *ares);

void            jfour1(float data[], unsigned long nn, int isign);

void            jrealft(float data[], unsigned long n, int isign);

int 
jtridib_(int *n, double *eps1, double *d, double *e, double *e2,
	 double *lb, double *ub, int *m11, int *m,
	 double *w, int *ind, int *ierr,
	 double *rv4, double *rv5);

int 
jtinvit_(int *nm, int *n, double *d, double *e, double *e2,
	 int *m, double *w, int *ind, double *z, int *ierr,
	 double *rv1, double *rv2,
	 double *rv3, double *rv4, double *rv6);

void            mt_get_spec(float *series, int inum, int klength, float *amp);

int             multitap(int n, int nwin, double *el, float npi, double *tapers, double *tapsum);

void            zero_pad(float output[], int start, int olength);

double          remove_mean(float x[], int lx);

float           get_cos_taper(int n, int k);

/***************************************************************/

void
main()
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



#if 1
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
		exit;
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


	increase = 1;
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
	freqwin = (int) (fWidth / df)/2;


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
		freq[i] = df * i;
		fprintf(inf, "%d %g %g %g %g %g %g\n", i, freq[i],
			spec[i], naive_spec[i], dtemp[i],
			dof[i], Fvalues[i]);
	}

}
/*--------------------------------------------------------*/
/*----------------mt_get_spec---------------------------*/
void
mt_get_spec(float *series, int inum, int klength, float *amp)
{
	/*
	 * series = input time series inum   = length of time series klength
	 * = number of elements in power spectrum (a power of 2) amp =
	 * returned power spectrum
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

	/*
	 * Fast Fourier Transform Routine:  here we are using the Numerical
	 * Recipes routine jrealft which returns the fft in the 1-D input
	 * array packed as pairs of real numbers. The jrealft routine
	 * requires the input array to start at index=1 so we must decrement
	 * the index of amp
	 */
	jrealft(amp - 1, nn, isign);

}
/*------------------------------------------------------------*/
/*----------------do_mtap_spec---------------------------*/
void
do_mtap_spec(float *data, int npoints, int kind,
	     int nwin, float npi, int inorm, float dt, float *ospec, float *dof, float *Fvalues, int klen)
{
	/*
	 * data = floating point input time series npoints = number of points
	 * in data kind = flag for choosing hires or adaptive weighting
	 * coefficients nwin = number of taper windows to calculate npi =
	 * order of the slepian functions inorm = flag for choice of
	 * normalization dt = sampling interval (time) ospec = output spctrum
	 * dof = degrees of freedom at each frequency Fvalues = Ftest value
	 * at each frequency estimate klen = number of frequecies calculated
	 * (power of 2)
	 * 
	 */

	int             i, j, k;
	double         *lambda, *tapers;
	long            len, longlen;
	float          *xt;
	FILE           *fopen(), *inf, *tapfile;
	FILE           *dof_file;

	int             logg;
	int             nn;
	float          *b;
	int             iwin, kk;

	/*************/
	double          anrm, norm;
	double         *ReSpec, *ImSpec;
	double         *sqr_spec, *amu;
	float          *amp, *fv;
	double          avamp, temp, sqramp;
	double          sum, *tapsum;
	/************/
	int             num_freqs;
	int             len_taps, num_freq_tap;

	double         *dcf, *degf, avar;
	int             n1, n2, kf;
	int             flag;
	int             one = 1;

	double          tem1, tem2;

	/*
	 * lambda = vector of eigenvalues   tapsum = sum of each taper, saved
	 * for use in adaptive weighting  tapers =  matrix of slepian tapers,
	 * packed in a 1D double array
	 */
	lambda = (double *) malloc((size_t) nwin * sizeof(double));
	tapsum = (double *) malloc((size_t) nwin * sizeof(double));

	len_taps = npoints * nwin;

	tapers = (double *) malloc((size_t) len_taps * sizeof(double));

	num_freqs = 1 + klen / 2;
	num_freq_tap = num_freqs * nwin;

	/* get a slepian taper  */

	k = multitap(npoints, nwin, lambda, npi, tapers, tapsum);
#if 0
	/* print out tapers for curiosity  */
	for (i = 0; i < npoints; i++) {
		for (j = 0; j < nwin; j++)
			fprintf(stderr, "%d %15.10f ", i, tapers[i + j * npoints]);
		prbl;
	}
#endif

	/* choose normalization based on inorm flag  */

	anrm = 1.;

	switch (inorm) {
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
	b = (float *) malloc((size_t) npoints * sizeof(float));

	amu = (double *) malloc((size_t) num_freqs * sizeof(double));
	sqr_spec = (double *) malloc((size_t) num_freq_tap * sizeof(double));
	ReSpec = (double *) malloc((size_t) num_freq_tap * sizeof(double));
	ImSpec = (double *) malloc((size_t) num_freq_tap * sizeof(double));

	for (iwin = 0; iwin < nwin; iwin++) {
		kk = iwin * npoints;
		kf = iwin * num_freqs;

		for (j = 0; j < npoints; j++)
			b[j] = data[j] * tapers[kk + j];	/* application of
								 * iwin-th taper   */

		amp = (float *) malloc((size_t) klen * sizeof(float));

		mt_get_spec(b, npoints, klen, amp);	/* calculate the
							 * eigenspectrum */

		free(b);

		sum = 0.0;

		/* get spectrum from real fourier transform    */

		norm = 1.0 / (anrm * anrm);

		for (i = 1; i < num_freqs - 1; i++) {
			if (2 * i + 1 > klen)
				fprintf(stderr, "error in index\n");
			if (i + kf > num_freq_tap)
				fprintf(stderr, "error in index\n");

			sqramp = SQR(amp[2 * i + 1]) + SQR(amp[2 * i]);

			ReSpec[i + kf] = amp[2 * i];
			ImSpec[i + kf] = amp[2 * i + 1];
			sqr_spec[i + kf] = norm * (sqramp);

			sum += sqramp;
		}
		sqr_spec[0 + kf] = norm * SQR(fabs(amp[0]));
		sqr_spec[num_freqs - 1 + kf] = norm * SQR(fabs(amp[1]));

		ReSpec[0 + kf] = amp[0];
		ImSpec[0 + kf] = 0.0;

		ReSpec[num_freqs - 1 + kf] = amp[1];
		ImSpec[num_freqs - 1 + kf] = 0.0;

		sum += sqr_spec[0 + kf] + sqr_spec[num_freqs - 1 + kf];

		if (num_freqs - 1 + kf > num_freq_tap)
			fprintf(stderr, "error in index\n");

		temp = sum / (double) num_freqs;
		if (temp > 0.0)
			avamp = sqrt(temp) / anrm;
		else {
			avamp = 0.0;
			/* fprintf(stderr," avamp = 0.0! \n"); */
		}


		free(amp);

	}

	fv = (float *) malloc((size_t) num_freqs * sizeof(float));

	/* choice of hi-res or adaptive weighting for spectra    */

	switch (kind) {
	case 1:

		hires(sqr_spec, lambda, nwin, num_freqs, amu);
		get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);

		for (i = 0; i < num_freqs; i++) {
			ospec[i] = amu[i];
			dof[i] = nwin - 1;
			Fvalues[i] = fv[i];
		}
		break;

	case 2:

		/* get avar = variance */

		n1 = 0;
		n2 = npoints;


		avar = 0.0;

		for (i = n1; i < n2; i++)
			avar += (data[i]) * (data[i]);


		switch (inorm) {
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

		dcf = (double *) malloc((size_t) num_freq_tap * sizeof(double));
		degf = (double *) malloc((size_t) num_freqs * sizeof(double));

		adwait(sqr_spec, dcf, lambda, nwin, num_freqs, amu, degf, avar);

		get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);

#if 1
		/*
		 * dump out the degrees of freedom to a file for later
		 * inspection
		 */
		if ((dof_file = fopen("dof_file", "w")) == NULL) {
			fprintf(stderr, "dof unable to open\n");
			return;
		}
		for (i = 0; i < num_freqs; i++) {
			fprintf(dof_file, "%f\n", degf[i]);
		}

		fclose(dof_file);
#endif

		/* rap up   */

		for (i = 0; i < num_freqs; i++) {
			ospec[i] = amu[i];
			dof[i] = degf[i];
			Fvalues[i] = fv[i];
		}

		free(dcf);
		free(degf);
		free(fv);

		break;
	}

	/* free up memory and return  */

	free(amu);
	free(sqr_spec);
	free(ReSpec);
	free(ImSpec);
	free(lambda);
	free(tapers);

}
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------multitap--------------------------*/
int
multitap(int num_points, int nwin, double *lam, float npi, double *tapers, double *tapsum)
{
	/*
	 * get the multitaper slepian functions: num_points = number of
	 * points in data stream nwin = number of windows lam= vector of
	 * eigenvalues npi = order of slepian functions tapsum = sum of each
	 * taper, saved for use in adaptive weighting  tapers =  matrix of
	 * slepian tapers, packed in a 1D double array
	 */

	int             i, j, k, kk;
	double         *z, ww, cs, ai, an, eps, rlu, rlb, aa;
	double          dfac, drat, gamma, bh, tapsq, TWOPI, DPI;
	double         *diag, *offdiag, *offsq;
	char           *k1[4];

	char            name[81];
	double         *scratch1, *scratch2, *scratch3, *scratch4, *scratch6;


	/* need to initialize iwflag = 0 */
	double          anpi;
	double         *ell;
	int             key, nbin, npad;
	int            *ip;
	double         *evecs;
	double         *zee;

	long            len;
	int             ierr;
	int             m11;
	DPI = (double) PI;
	TWOPI = (double) 2 *DPI;

	anpi = npi;
	an = (double) (num_points);
	ww = (double) (anpi) / an;	/* this corresponds to P&W's W value  */
	cs = cos(TWOPI * ww);


	ell = (double *) malloc((size_t) nwin * sizeof(double));

	diag = (double *) malloc((size_t) num_points * sizeof(double));

	offdiag = (double *) malloc((size_t) num_points * sizeof(double));
	offsq = (double *) malloc((size_t) num_points * sizeof(double));

	scratch1 = (double *) malloc((size_t) num_points * sizeof(double));
	scratch2 = (double *) malloc((size_t) num_points * sizeof(double));
	scratch3 = (double *) malloc((size_t) num_points * sizeof(double));
	scratch4 = (double *) malloc((size_t) num_points * sizeof(double));
	scratch6 = (double *) malloc((size_t) num_points * sizeof(double));

	/* make the diagonal elements of the tridiag matrix  */

	for (i = 0; i < num_points; i++) {
		ai = (double) (i);
		diag[i] = -cs * (((an - 1.) / 2. - ai)) * (((an - 1.) / 2. - ai));
		offdiag[i] = -ai * (an - ai) / 2.;
		offsq[i] = offdiag[i] * offdiag[i];
	}

	eps = 1.0e-13;
	m11 = 1;

	ip = (int *) malloc((size_t) nwin * sizeof(int));

	/* call the eispac routines to invert the tridiagonal system */

	jtridib_(&num_points, &eps, diag, offdiag, offsq, &rlb, &rlu, &m11, &nwin, lam,
		 ip, &ierr, scratch1, scratch2);
#if DIAG1
	fprintf(stderr, "ierr=%d rlb=%.8f rlu=%.8f\n", ierr, rlb, rlu);

	fprintf(stderr, "eigenvalues for the eigentapers\n");

	for (k = 0; k < nwin; k++)
		fprintf(stderr, "%.20f ", lam[k]);
	fprintf(stderr, "\n");
#endif


	len = num_points * nwin;

	evecs = (double *) malloc((size_t) len * sizeof(double));



	jtinvit_(&num_points, &num_points, diag, offdiag, offsq, &nwin, lam, ip, evecs, &ierr,
		 scratch1, scratch2, scratch3, scratch4, scratch6);

	free(scratch1);
	free(scratch2);
	free(scratch3);
	free(scratch4);
	free(scratch6);

	/*
	 * we calculate the eigenvalues of the dirichlet-kernel problem i.e.
	 * the bandwidth retention factors from slepian 1978 asymptotic
	 * formula, gotten from thomson 1982 eq 2.5 supplemented by the
	 * asymptotic formula for k near 2n from slepian 1978 eq 61 more
	 * precise values of these parameters, perhaps useful in adaptive
	 * spectral estimation, can be calculated explicitly using the
	 * rayleigh-quotient formulas in thomson (1982) and park et al (1987)
	 * 
	 */
	dfac = (double) an *DPI * ww;
	drat = (double) 8. *dfac;


	dfac = (double) 4. *sqrt(DPI * dfac) * exp((double) (-2.0) * dfac);


	for (k = 0; k < nwin; k++) {
		lam[k] = (double) 1.0 - (double) dfac;
		dfac = dfac * drat / (double) (k + 1);



		/* fails as k -> 2n */
	}


	gamma = log((double) 8. * an * sin((double) 2. * DPI * ww)) + (double) 0.5772156649;



	for (k = 0; k < nwin; k++) {
		bh = -2. * DPI * (an * ww - (double) (k) /
				  (double) 2. - (double) .25) / gamma;
		ell[k] = (double) 1. / ((double) 1. + exp(DPI * (double) bh));

	}

	for (i = 0; i < nwin; i++)
		lam[i] = MAX(ell[i], lam[i]);

	/************************************************************
        c   normalize the eigentapers to preserve power for a white process
        c   i.e. they have rms value unity
        c  tapsum is the average of the eigentaper, should be near zero for
        c  antisymmetric tapers
        ************************************************************/

	for (k = 0; k < nwin; k++) {
		kk = (k) * num_points;
		tapsum[k] = 0.;
		tapsq = 0.;
		for (i = 0; i < num_points; i++) {
			aa = evecs[i + kk];
			tapers[i + kk] = aa;
			tapsum[k] = tapsum[k] + aa;
			tapsq = tapsq + aa * aa;
		}
		aa = sqrt(tapsq / (double) num_points);
		tapsum[k] = tapsum[k] / aa;

		for (i = 0; i < num_points; i++) {
			tapers[i + kk] = tapers[i + kk] / aa;

		}
	}

	/* Free Memory */

	free(ell);
	free(diag);
	free(offdiag);
	free(offsq);
	free(ip);

	free(evecs);

	return 1;
}
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*----------------adwait--------------------------------*/
int 
adwait(double *sqr_spec, double *dcf,
double *el, int nwin, int num_freq, double *ares, double *degf, double avar)
{
	/*
	 * c  this version uses thomson's algorithm for calculating c  the
	 * adaptive spectrum estimate
	 */
	double          as, das, tol, a1, scale, ax, fn, fx;
	double         *spw, *bias;
	double          test_tol, dif;
	int             jitter, i, j, k, kpoint, jloop;
	float           df;
	/*
	 * c  set tolerance for iterative scheme exit
	 */

#if 0
	fprintf(stderr, "test input\n adwait: %d %d %f\n", nwin, num_freq, avar);
	fprintf(stderr, "\n Data=\n");
	for (i = 0; i < num_freq; i++) {
		fprintf(stderr, "%d %f \n", i, sqr_spec[i]);
	}
#endif


	tol = 3.0e-4;
	jitter = 0;
	scale = avar;
	/***********************************
        c  we scale the bias by the total variance of the frequency transform
        c  from zero freq to the nyquist
        c  in this application we scale the eigenspectra by the bias in order to avoid
        c  possible floating point overflow
        ************************************/
	spw = (double *) malloc((size_t) nwin * sizeof(double));
	bias = (double *) malloc((size_t) nwin * sizeof(double));


	for (i = 0; i < nwin; i++) {

		bias[i] = (1.00 - el[i]);
	}

	/*
	 * for( i=1;i<=nwin; i++) fprintf(stderr,"%f %f\n",el[i], bias[i]);
	 * fprintf(stderr,"\n");
	 */

	/* START do 100 */
	for (jloop = 0; jloop < num_freq; jloop++) {

		for (i = 0; i < nwin; i++) {
			kpoint = jloop + i * num_freq;
			spw[i] = (sqr_spec[kpoint]) / scale;
		}
		/********************************************
                c  first guess is the average of the two
                    lowest-order eigenspectral estimates
                 ********************************************/
		as = (spw[0] + spw[1]) / 2.00;

		/* START do 300 */
		/* c  find coefficients */

		for (k = 0; k < 20; k++) {
			fn = 0.00;
			fx = 0.00;

			for (i = 0; i < nwin; i++) {
				a1 = sqrt(el[i]) * as / (el[i] * as + bias[i]);
				a1 = a1 * a1;
				fn = fn + a1 * spw[i];
				fx = fx + a1;
			}


			ax = fn / fx;
			dif = ax - as;
			das = ABS(dif);
			/*
			 * fprintf(stderr,"adwait: jloop = %d k=%d %g %g %g
			 * %g\n",jloop,k, fn,fx,ax,das);
			 */
			test_tol = das / as;
			if (test_tol < tol) {
				break;
			}
			as = ax;
		}

		/* fprintf(stderr,"adwait: k=%d test_tol=%f\n",k, test_tol); */
		/* end  300  */

		/* c  flag if iteration does not converge */

		if (k >= 20)
			jitter++;

		ares[jloop] = as * scale;
		/* c   calculate degrees of freedom */
		df = 0.0;
		for (i = 0; i < nwin; i++) {
			kpoint = jloop + i * num_freq;
			dcf[kpoint] = sqrt(el[i]) * as / (el[i] * as + bias[i]);
			df = df + dcf[kpoint] * dcf[kpoint];
		}
		/*
		 * we normalize degrees of freedom by the weight of the first
		 * eigenspectrum this way we never have fewer than two
		 * degrees of freedom
		 */

		degf[jloop] = df * 2. / (dcf[jloop] * dcf[jloop]);

	}			/* end 100 */

	fprintf(stderr, "%d failed iterations\n", jitter);
	free(spw);
	free(bias);

	return jitter;
}
/*--------------------------------------------------------*/
/*------------get_F_values----------------------------*/
/*--------------------------------------------------------*/
void
get_F_values(double *sr, double *si, int nf, int nwin, float *Fvalue, double *b)
{
	/*
	 * b is fft of slepian eigentapers at zero freq sr si are the
	 * eigenspectra amu contains line frequency estimates and f-test
	 * parameter
	 */
	double          sum, sumr, sumi, sum2;
	int             i, j, k;
	double         *amur, *amui;
	sum = 0.;

	amur = (double *) malloc((size_t) nf * sizeof(double));
	amui = (double *) malloc((size_t) nf * sizeof(double));




	for (i = 0; i < nwin; i++) {

		sum = sum + b[i] * b[i];
	}
	for (i = 0; i < nf; i++) {
		amur[i] = 0.;
		amui[i] = 0.;
		for (j = 0; j < nwin; j++) {
			k = i + j * nf;
			amur[i] = amur[i] + sr[k] * b[j];
			amui[i] = amui[i] + si[k] * b[j];
		}
		amur[i] = amur[i] / sum;
		amui[i] = amui[i] / sum;
		sum2 = 0.;
		for (j = 0; j < nwin; j++) {
			k = i + j * nf;
			sumr = sr[k] - amur[i] * b[j];
			sumi = si[k] - amui[i] * b[j];
			sum2 = sum2 + sumr * sumr + sumi * sumi;
		}
		Fvalue[i] = (float) (nwin - 1) * (SQR(amui[i]) + SQR(amur[i])) * sum / sum2;
	}
	free(amui);
	free(amur);
	return;
}

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
/*-------------  HIRES  ----------------------------------*/
/*-------------------------------------------------------*/
int
hires(double *sqr_spec, double *el, int nwin, int num_freq, double *ares)
{
	int             i, j, k, kpoint;
	float           a;

	for (j = 0; j < num_freq; j++)
		ares[j] = 0.;

	for (i = 0; i < nwin; i++) {
		k = i * num_freq;
		a = 1. / (el[i] * nwin);
		for (j = 0; j < num_freq; j++) {
			kpoint = j + k;
			ares[j] = ares[j] +
				a * (sqr_spec[kpoint]);
		}
	}

	for (j = 0; j < num_freq; j++) {
		if (ares[j] > 0.0)
			ares[j] = sqrt(ares[j]);
		else
			printf("sqrt problem in hires pos=%d %f\n", j, ares[j]);
	}

	return 1;
}
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*---------------jtinvit----------------------------------*/
#include <math.h>
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)


/* ./ add name=tinvit */

/* ------------------------------------------------------------------ */

int 
jtinvit_(int *nm, int *n, double *d, double *e, double *e2, int *m, double *w, int *ind, double *z, int *ierr, double *rv1, double *rv2,
	 double *rv3, double *rv4, double *rv6)
{
	/* Initialized data */

	static double   machep = 1.25e-15;

	/* System generated locals */
	int             z_dim1, z_offset, i1, i2, i3;
	double          d1, d2;

	/* Builtin functions */
	double          sqrt();

	/* Local variables */
	static double   norm;
	static int      i, j, p, q, r, s;
	static double   u, v, order;
	static int      group;
	static double   x0, x1;
	static int      ii, jj, ip;
	static double   uk, xu;
	static int      tag, its;
	static double   eps2, eps3, eps4;

	static double   rtem;

	/* this subroutine is a translation of the inverse iteration tech- */
	/* nique in the algol procedure tristurm by peters and wilkinson. */
	/* handbook for auto. comp., vol.ii-linear algebra, 418-439(1971). */

	/* this subroutine finds those eigenvectors of a tridiagonal */
	/* symmetric matrix corresponding to specified eigenvalues, */
	/* using inverse iteration. */

	/* on input: */

	/* nm must be set to the row dimension of two-dimensional */
	/* array parameters as declared in the calling program */
	/* dimension statement; */

	/* n is the order of the matrix; */

	/* d contains the diagonal elements of the input matrix; */

	/* e contains the subdiagonal elements of the input matrix */
	/* in its last n-1 positions.  e(1) is arbitrary; */

	/* e2 contains the squares of the corresponding elements of e, */
	/* with zeros corresponding to negligible elements of e. */
	/* e(i) is considered negligible if it is not larger than */
	/* the product of the relative machine precision and the sum */
	/* of the magnitudes of d(i) and d(i-1).  e2(1) must contain */
	/* 0.0d0 if the eigenvalues are in ascending order, or 2.0d0 */
	/* if the eigenvalues are in descending order.  if  bisect, */
	/* tridib, or  imtqlv  has been used to find the eigenvalues, */
	/* their output e2 array is exactly what is expected here; */

	/* m is the number of specified eigenvalues; */

	/*
	 * w contains the m eigenvalues in ascending or descending order;
	 */

	/* ind contains in its first m positions the submatrix indices */
	/* associated with the corresponding eigenvalues in w -- */
	/* 1 for eigenvalues belonging to the first submatrix from */
	/*
	 * the top, 2 for those belonging to the second submatrix, etc.
	 */

	/* on output: */

	/* all input arrays are unaltered; */

	/* z contains the associated set of orthonormal eigenvectors. */
	/* any vector which fails to converge is set to zero; */

	/* ierr is set to */
	/* zero       for normal return, */
	/* -r         if the eigenvector corresponding to the r-th */
	/* eigenvalue fails to converge in 5 iterations; */

	/* rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays. */

	/* questions and comments should be directed to b. s. garbow, */
	/* applied mathematics division, argonne national laboratory */

	/*
	 * ------------------------------------------------------------------
	 */

	/* :::::::::: machep is a machine dependent parameter specifying */
	/* the relative precision of floating point arithmetic. */
	/* machep = 16.0d0**(-13) for long form arithmetic */
	/* on s360 :::::::::: */
	/* for f_floating dec fortran */
	/* data machep/1.1d-16/ */
	/* for g_floating dec fortran */
	/* Parameter adjustments */
	--rv6;
	--rv4;
	--rv3;
	--rv2;
	--rv1;
	--e2;
	--e;
	--d;
	z_dim1 = *nm;
	z_offset = z_dim1 + 1;
	z -= z_offset;
	--ind;
	--w;

	/* Function Body */

	*ierr = 0;
	if (*m == 0) {
		goto L1001;
	}
	tag = 0;
	order = 1. - e2[1];
	q = 0;
	/* :::::::::: establish and process next submatrix :::::::::: */
L100:
	p = q + 1;

	i1 = *n;
	for (q = p; q <= i1; ++q) {
		if (q == *n) {
			goto L140;
		}
		if (e2[q + 1] == 0.) {
			goto L140;
		}
		/* L120: */
	}
	/* :::::::::: find vectors by inverse iteration :::::::::: */
L140:
	++tag;
	s = 0;

	i1 = *m;
	for (r = 1; r <= i1; ++r) {
		if (ind[r] != tag) {
			goto L920;
		}
		its = 1;
		x1 = w[r];
		if (s != 0) {
			goto L510;
		}
		/* :::::::::: check for isolated root :::::::::: */
		xu = 1.;
		if (p != q) {
			goto L490;
		}
		rv6[p] = 1.;
		goto L870;
L490:
		norm = (d1 = d[p], abs(d1));
		ip = p + 1;

		i2 = q;
		for (i = ip; i <= i2; ++i) {
			/* L500: */
			norm = norm + (d1 = d[i], abs(d1)) + (d2 = e[i], abs(d2));
		}
		/* :::::::::: eps2 is the criterion for grouping, */
		/* eps3 replaces zero pivots and equal */
		/* roots are modified by eps3, */
		/*
		 * eps4 is taken very small to avoid overflow ::::::::: :
		 */
		eps2 = norm * .001;
		eps3 = machep * norm;
		uk = (double) (q - p + 1);
		eps4 = uk * eps3;
		uk = eps4 / sqrt(uk);
		s = p;
L505:
		group = 0;
		goto L520;
		/* :::::::::: look for close or coincident roots :::::::::: */
L510:
		if ((d1 = x1 - x0, abs(d1)) >= eps2) {
			goto L505;
		}
		++group;
		if (order * (x1 - x0) <= 0.) {
			x1 = x0 + order * eps3;
		}
		/* :::::::::: elimination with interchanges and */
		/* initialization of vector :::::::::: */
L520:
		v = 0.;

		i2 = q;
		for (i = p; i <= i2; ++i) {
			rv6[i] = uk;
			if (i == p) {
				goto L560;
			}
			if ((d1 = e[i], abs(d1)) < abs(u)) {
				goto L540;
			}
			/*
			 * :::::::::: warning -- a divide check may occur
			 * here if
			 */
			/*
			 * e2 array has not been specified correctly ::::::
			 * ::::
			 */
			xu = u / e[i];
			rv4[i] = xu;
			rv1[i - 1] = e[i];
			rv2[i - 1] = d[i] - x1;
			rv3[i - 1] = 0.;
			if (i != q) {
				rv3[i - 1] = e[i + 1];
			}
			u = v - xu * rv2[i - 1];
			v = -xu * rv3[i - 1];
			goto L580;
	L540:
			xu = e[i] / u;
			rv4[i] = xu;
			rv1[i - 1] = u;
			rv2[i - 1] = v;
			rv3[i - 1] = 0.;
	L560:
			u = d[i] - x1 - xu * v;
			if (i != q) {
				v = e[i + 1];
			}
	L580:
			;
		}

		if (u == 0.) {
			u = eps3;
		}
		rv1[q] = u;
		rv2[q] = 0.;
		rv3[q] = 0.;
		/* :::::::::: back substitution */
		/* for i=q step -1 until p do -- :::::::::: */
L600:
		i2 = q;
		for (ii = p; ii <= i2; ++ii) {
			i = p + q - ii;
			rtem = rv6[i] - u * rv2[i] - v * rv3[i];
			rv6[i] = (rtem) / rv1[i];
			v = u;
			u = rv6[i];
			/* L620: */
		}
		/* :::::::::: orthogonalize with respect to previous */
		/* members of group :::::::::: */
		if (group == 0) {
			goto L700;
		}
		j = r;

		i2 = group;
		for (jj = 1; jj <= i2; ++jj) {
	L630:
			--j;
			if (ind[j] != tag) {
				goto L630;
			}
			xu = 0.;

			i3 = q;
			for (i = p; i <= i3; ++i) {
				/* L640: */
				xu += rv6[i] * z[i + j * z_dim1];
			}

			i3 = q;
			for (i = p; i <= i3; ++i) {
				/* L660: */
				rv6[i] -= xu * z[i + j * z_dim1];
			}

			/* L680: */
		}

L700:
		norm = 0.;

		i2 = q;
		for (i = p; i <= i2; ++i) {
			/* L720: */
			norm += (d1 = rv6[i], abs(d1));
		}

		if (norm >= 1.) {
			goto L840;
		}
		/* :::::::::: forward substitution :::::::::: */
		if (its == 5) {
			goto L830;
		}
		if (norm != 0.) {
			goto L740;
		}
		rv6[s] = eps4;
		++s;
		if (s > q) {
			s = p;
		}
		goto L780;
L740:
		xu = eps4 / norm;

		i2 = q;
		for (i = p; i <= i2; ++i) {
			/* L760: */
			rv6[i] *= xu;
		}
		/* :::::::::: elimination operations on next vector */
		/* iterate :::::::::: */
L780:
		i2 = q;
		for (i = ip; i <= i2; ++i) {
			u = rv6[i];
			/*
			 * :::::::::: if rv1(i-1) .eq. e(i), a row
			 * interchange
			 */
			/* was performed earlier in the */
			/* triangularization process :::::::::: */
			if (rv1[i - 1] != e[i]) {
				goto L800;
			}
			u = rv6[i - 1];
			rv6[i - 1] = rv6[i];
	L800:
			rv6[i] = u - rv4[i] * rv6[i - 1];
			/* L820: */
		}

		++its;
		goto L600;
		/*
		 * :::::::::: set error -- non-converged eigenvector
		 * ::::::::::
		 */
L830:
		*ierr = -r;
		xu = 0.;
		goto L870;
		/* :::::::::: normalize so that sum of squares is */
		/* 1 and expand to full order :::::::::: */
L840:
		u = 0.;

		i2 = q;
		for (i = p; i <= i2; ++i) {
			/* L860: */
			/* Computing 2nd power */
			d1 = rv6[i];
			u += d1 * d1;
		}

		xu = 1. / sqrt(u);

L870:
		i2 = *n;
		for (i = 1; i <= i2; ++i) {
			/* L880: */
			z[i + r * z_dim1] = 0.;
		}

		i2 = q;
		for (i = p; i <= i2; ++i) {
			/* L900: */
			z[i + r * z_dim1] = rv6[i] * xu;
		}

		x0 = x1;
L920:
		;
	}

	if (q < *n) {
		goto L100;
	}
L1001:
	return 0;
	/* :::::::::: last card of tinvit :::::::::: */
}				/* tinvit_ */

/*--------------------------------------------------------*/
/*----------------jtridib---------------------------------*/
/*--------------------------------------------------------*/
int 
jtridib_(int *n, double *eps1, double *d, double *e, double *e2, double *lb, double *ub, int *m11, int *m, double *w, int *ind, int *ierr,
	 double *rv4, double *rv5)
{
	/* Initialized data */

	static double   machep = 1.25e-15;

	/* System generated locals */
	int             i1, i2;
	double          d1, d2, d3;

	/* Local variables */
	static int      i, j, k, l, p, q, r, s;
	static double   u, v;
	static int      m1, m2;
	static double   t1, t2, x0, x1;
	static int      m22, ii;
	static double   xu;
	static int      isturm, tag;



	/* this subroutine is a translation of the algol procedure bisect, */
	/* num. math. 9, 386-393(1967) by barth, martin, and wilkinson. */
	/* handbook for auto. comp., vol.ii-linear algebra, 249-256(1971). */

	/* this subroutine finds those eigenvalues of a tridiagonal */
	/* symmetric matrix between specified boundary indices, */
	/* using bisection. */

	/* on input: */

	/* n is the order of the matrix; */

	/* eps1 is an absolute error tolerance for the computed */
	/* eigenvalues.  if the input eps1 is non-positive, */
	/* it is reset for each submatrix to a default value, */
	/* namely, minus the product of the relative machine */
	/* precision and the 1-norm of the submatrix; */

	/* d contains the diagonal elements of the input matrix; */

	/* e contains the subdiagonal elements of the input matrix */
	/* in its last n-1 positions.  e(1) is arbitrary; */

	/* e2 contains the squares of the corresponding elements of e. */
	/* e2(1) is arbitrary; */

	/* m11 specifies the lower boundary index for the desired */
	/* eigenvalues; */

	/* m specifies the number of eigenvalues desired.  the upper */
	/* boundary index m22 is then obtained as m22=m11+m-1. */

	/* on output: */

	/* eps1 is unaltered unless it has been reset to its */
	/* (last) default value; */

	/* d and e are unaltered; */

	/* elements of e2, corresponding to elements of e regarded */
	/* as negligible, have been replaced by zero causing the */
	/* matrix to split into a direct sum of submatrices. */
	/* e2(1) is also set to zero; */

	/* lb and ub define an interval containing exactly the desired */
	/* eigenvalues; */

	/* w contains, in its first m positions, the eigenvalues */
	/* between indices m11 and m22 in ascending order; */

	/* ind contains in its first m positions the submatrix indices */
	/* associated with the corresponding eigenvalues in w -- */
	/* 1 for eigenvalues belonging to the first submatrix from */
	/*
	 * the top, 2 for those belonging to the second submatrix, etc.;
	 */

	/* ierr is set to */
	/* zero       for normal return, */
	/* 3*n+1      if multiple eigenvalues at index m11 make */
	/* unique selection impossible, */
	/* 3*n+2      if multiple eigenvalues at index m22 make */
	/* unique selection impossible; */

	/* rv4 and rv5 are temporary storage arrays. */

	/* note that subroutine tql1, imtql1, or tqlrat is generally faster */
	/* than tridib, if more than n/4 eigenvalues are to be found. */

	/* questions and comments should be directed to b. s. garbow, */
	/* applied mathematics division, argonne national laboratory */

	/*
	 * ------------------------------------------------------------------
	 */

	/* :::::::::: machep is a machine dependent parameter specifying */
	/* the relative precision of floating point arithmetic. */
	/* machep = 16.0d0**(-13) for long form arithmetic */
	/* on s360 :::::::::: */
	/* for f_floating dec fortran */
	/* data machep/1.1d-16/ */
	/* for g_floating dec fortran */
	/* Parameter adjustments */
	--rv5;
	--rv4;
	--e2;
	--e;
	--d;
	--ind;
	--w;

	/* Function Body */

	*ierr = 0;
	tag = 0;
	xu = d[1];
	x0 = d[1];
	u = 0.;
	/* :::::::::: look for small sub-diagonal entries and determine an */
	/* interval containing all the eigenvalues :::::::::: */
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
		x1 = u;
		u = 0.;
		if (i != *n) {
			u = (d1 = e[i + 1], abs(d1));
		}
		/* Computing MIN */
		d1 = d[i] - (x1 + u);
		xu = min(d1, xu);
		/* Computing MAX */
		d1 = d[i] + (x1 + u);
		x0 = max(d1, x0);
		if (i == 1) {
			goto L20;
		}
		if ((d1 = e[i], abs(d1)) > machep * ((d2 = d[i], abs(d2)) + (
						 d3 = d[i - 1], abs(d3)))) {
			goto L40;
		}
L20:
		e2[i] = 0.;
L40:
		;
	}

	/* Computing MAX */
	d1 = abs(xu), d2 = abs(x0);
	x1 = max(d1, d2) * machep * (double) (*n);
	xu -= x1;
	t1 = xu;
	x0 += x1;
	t2 = x0;
	/* :::::::::: determine an interval containing exactly */
	/* the desired eigenvalues :::::::::: */
	p = 1;
	q = *n;
	m1 = *m11 - 1;
	if (m1 == 0) {
		goto L75;
	}
	isturm = 1;
L50:
	v = x1;
	x1 = xu + (x0 - xu) * .5;
	if (x1 == v) {
		goto L980;
	}
	goto L320;
L60:
	if ((i1 = s - m1) < 0) {
		goto L65;
	} else if (i1 == 0) {
		goto L73;
	} else {
		goto L70;
	}
L65:
	xu = x1;
	goto L50;
L70:
	x0 = x1;
	goto L50;
L73:
	xu = x1;
	t1 = x1;
L75:
	m22 = m1 + *m;
	if (m22 == *n) {
		goto L90;
	}
	x0 = t2;
	isturm = 2;
	goto L50;
L80:
	if ((i1 = s - m22) < 0) {
		goto L65;
	} else if (i1 == 0) {
		goto L85;
	} else {
		goto L70;
	}
L85:
	t2 = x1;
L90:
	q = 0;
	r = 0;
	/* :::::::::: establish and process next submatrix, refining */
	/* interval by the gerschgorin bounds :::::::::: */
L100:
	if (r == *m) {
		goto L1001;
	}
	++tag;
	p = q + 1;
	xu = d[p];
	x0 = d[p];
	u = 0.;

	i1 = *n;
	for (q = p; q <= i1; ++q) {
		x1 = u;
		u = 0.;
		v = 0.;
		if (q == *n) {
			goto L110;
		}
		u = (d1 = e[q + 1], abs(d1));
		v = e2[q + 1];
L110:
		/* Computing MIN */
		d1 = d[q] - (x1 + u);
		xu = min(d1, xu);
		/* Computing MAX */
		d1 = d[q] + (x1 + u);
		x0 = max(d1, x0);
		if (v == 0.) {
			goto L140;
		}
		/* L120: */
	}

L140:
	/* Computing MAX */
	d1 = abs(xu), d2 = abs(x0);
	x1 = max(d1, d2) * machep;
	if (*eps1 <= 0.) {
		*eps1 = -x1;
	}
	if (p != q) {
		goto L180;
	}
	/* :::::::::: check for isolated root within interval :::::::::: */
	if (t1 > d[p] || d[p] >= t2) {
		goto L940;
	}
	m1 = p;
	m2 = p;
	rv5[p] = d[p];
	goto L900;
L180:
	x1 *= (double) (q - p + 1);
	/* Computing MAX */
	d1 = t1, d2 = xu - x1;
	*lb = max(d1, d2);
	/* Computing MIN */
	d1 = t2, d2 = x0 + x1;
	*ub = min(d1, d2);
	x1 = *lb;
	isturm = 3;
	goto L320;
L200:
	m1 = s + 1;
	x1 = *ub;
	isturm = 4;
	goto L320;
L220:
	m2 = s;
	if (m1 > m2) {
		goto L940;
	}
	/* :::::::::: find roots by bisection :::::::::: */
	x0 = *ub;
	isturm = 5;

	i1 = m2;
	for (i = m1; i <= i1; ++i) {
		rv5[i] = *ub;
		rv4[i] = *lb;
		/* L240: */
	}
	/* :::::::::: loop for k-th eigenvalue */
	/* for k=m2 step -1 until m1 do -- */
	/*
	 * (-do- not used to legalize -computed go to-) ::::::::::
	 */
	k = m2;
L250:
	xu = *lb;
	/* :::::::::: for i=k step -1 until m1 do -- :::::::::: */
	i1 = k;
	for (ii = m1; ii <= i1; ++ii) {
		i = m1 + k - ii;
		if (xu >= rv4[i]) {
			goto L260;
		}
		xu = rv4[i];
		goto L280;
L260:
		;
	}

L280:
	if (x0 > rv5[k]) {
		x0 = rv5[k];
	}
	/* :::::::::: next bisection step :::::::::: */
L300:
	x1 = (xu + x0) * .5;
	if (x0 - xu <= machep * 2. * (abs(xu) + abs(x0)) + abs(*eps1)) {
		goto L420;
	}
	/* :::::::::: in-line procedure for sturm sequence :::::::::: */
L320:
	s = p - 1;
	u = 1.;

	i1 = q;
	for (i = p; i <= i1; ++i) {
		if (u != 0.) {
			goto L325;
		}
		v = (d1 = e[i], abs(d1)) / machep;
		if (e2[i] == 0.) {
			v = 0.;
		}
		goto L330;
L325:
		v = e2[i] / u;
L330:
		u = d[i] - x1 - v;
		if (u < 0.) {
			++s;
		}
		/* L340: */
	}

	switch ((int) isturm) {
	case 1:
		goto L60;
	case 2:
		goto L80;
	case 3:
		goto L200;
	case 4:
		goto L220;
	case 5:
		goto L360;
	}
	/* :::::::::: refine intervals :::::::::: */
L360:
	if (s >= k) {
		goto L400;
	}
	xu = x1;
	if (s >= m1) {
		goto L380;
	}
	rv4[m1] = x1;
	goto L300;
L380:
	rv4[s + 1] = x1;
	if (rv5[s] > x1) {
		rv5[s] = x1;
	}
	goto L300;
L400:
	x0 = x1;
	goto L300;
	/* :::::::::: k-th eigenvalue found :::::::::: */
L420:
	rv5[k] = x1;
	--k;
	if (k >= m1) {
		goto L250;
	}
	/* :::::::::: order eigenvalues tagged with their */
	/* submatrix associations :::::::::: */
L900:
	s = r;
	r = r + m2 - m1 + 1;
	j = 1;
	k = m1;

	i1 = r;
	for (l = 1; l <= i1; ++l) {
		if (j > s) {
			goto L910;
		}
		if (k > m2) {
			goto L940;
		}
		if (rv5[k] >= w[l]) {
			goto L915;
		}
		i2 = s;
		for (ii = j; ii <= i2; ++ii) {
			i = l + s - ii;
			w[i + 1] = w[i];
			ind[i + 1] = ind[i];
			/* L905: */
		}

L910:
		w[l] = rv5[k];
		ind[l] = tag;
		++k;
		goto L920;
L915:
		++j;
L920:
		;
	}

L940:
	if (q < *n) {
		goto L100;
	}
	goto L1001;
	/* :::::::::: set error -- interval cannot be found containing */
	/* exactly the desired eigenvalues :::::::::: */
L980:
	*ierr = *n * 3 + isturm;
L1001:
	*lb = t1;
	*ub = t2;
	return 0;
	/* :::::::::: last card of tridib :::::::::: */
}				/* tridib_ */

/*--------------------------------------------------------*/
/*---------------get_pow_2----------------------------*/
/*--------------------------------------------------------*/
int
get_pow_2(int inum)
{
	int             j, klength;
	/* find smallest power of 2 that encompasses the data */

	for (j = 1; pow((double) 2, (double) j) < inum; j++);
	return klength = pow((double) 2, (double) j);
}

/*--------------------------------------------------------*/
/*-----------remove_mean----------------------------*/
/*--------------------------------------------------------*/

double 
remove_mean(float x[], int lx)
{

	int             k;
	double          mean;
	mean = 0.;
	if (lx < 2)
		return mean;

	for (k = 0; k <= lx; k++) {
		mean = x[k] + mean;
	}

	mean = mean / (float) lx;

	for (k = 0; k <= lx; k++) {
		x[k] = x[k] - mean;
	}


	return mean;
}
/*********************************************************/
/*--------------------------------------------------------*/
/*-----------zero_pad----------------------------------*/
/*--------------------------------------------------------*/

void 
zero_pad(float output[], int start, int olength)
{
	int             i;
	for (i = start; i < olength; i++) {
		output[i] = 0.0;
	}
}
/*--------------------------------------------------------*/
/*-----------get_cos_taper----------------------------*/
/*--------------------------------------------------------*/

float 
get_cos_taper(int n, int k)
{
	int             l;
	float           vwin;
	vwin = 0.0;

	if (k < 0 || k > n)
		return vwin;
	vwin = 1.0;

	l = (n - 2) / 10;
	if (k <= l)
		vwin = 0.5 * (1.0 - cos(k * PI / (l + 1)));
	if (k >= n - l - 2)
		vwin = 0.5 * (1.0 - cos((n - k - 1) * PI / (l + 1)));
	return vwin;
}
