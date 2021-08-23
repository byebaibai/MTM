#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265358979
#define ABS(a) ((a) < (0) ? -(a) : (a))

#define NRANSI
/* #include "jl.h" */
#include "nrutil.h"
/* #include "nr.h" */

/** FUNC DEF */ void  get_F_values(double *sr, double *si, int nf, int nwin, float *Fvalue, double *b)
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
	amur = dvector((long) 0, (long) nf);
	amui = dvector((long) 0, (long) nf);



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
                /* percival and walden, eq 499c, p499 */
              /* sum = Hk(0) squared  */
	}
	return;
}
