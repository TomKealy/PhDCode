/*
* interp1_table_adj_mex.c
* Mex file for *adjoint* of 1D periodic interpolation using table lookup.
*
* forward direction: (for m = 1,...,M)
* f(t_m) = \sum_{k=0}^{K-1} c_k h( (t_m - k) mod K )
*
* adjoint direction: (for k=0,...,K-1) (note complex conjugate!)
* c_k = \sum_{m=1}^M f(t_m) h^*( (t_m - k) mod K )
*
* The interpolator h is nonzero (and tabulated) for -J/2 <= t <= J/2.
*
* Copyright 2004-4-1 Jeff Fessler and Yingying Zhang, The University of Michigan
*/
#include "mex.h"
#include "math.h"
#include "string.h"
#include "def,table.h"


static void interp1_table_complex_adj(
double *r_ck,		/* [K1,1] out */
double *i_ck,
const int K1,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,	/* imaginary part of complex kernel */
const int J1,
const int L1,
const double *pt,	/* [M,1] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm, jj1;
	/* Precompute some params: Note index begins from 0 in C */
	const int ncenter = floor(J1 * L1/2);	/* ? */
	const int J_shift = (J1 % 2) ? (J1+1)/2 : J1/2;	/* nufft_offset */

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double tval = *pt++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;

	/* put t in range [0,K1-1] */
	const double tm = tval - K1 * floor(tval / K1);
	const int koff = (J1%2==1)
		? ( round(tm) - J_shift )
		: ( floor(tm) - J_shift );

	for (jj1=0; jj1<J1; jj1++) {
		const int k1 = koff + jj1 + 1;
		const int n1 = ncenter + round((tm - k1) * L1);
		register double coefr = r_h1[n1];
		register double coefi = i_h1[n1];
		const int kk = (k1 + K1) % K1;

		/* instead of f = h c, we have c += h^* f */
		r_ck[kk] += coefr * fmr + coefi * fmi;
		i_ck[kk] += coefr * fmi - coefi * fmr;
	}
    }
}


static void interp1_table_real_adj(
double *r_ck,		/* [K1,1] out */
double *i_ck,
const int K1,
const double *r_h1,	/* [J1*L1+1,1] in (real) */
const int J1,
const int L1,
const double *pt,	/* [M,1] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm, jj1;
	/* Precompute some params: Note index begins from 0 in C */
	const int ncenter = floor(J1 * L1/2);	/* ? */
	const int J_shift = (J1%2) ? (J1+1)/2 : J1/2;	/* nufft_offset */

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double tval = *pt++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;

	/* put t in range [0,K1-1] */
	const double tm = tval - K1 * floor(tval / K1);
	const int koff = (J1%2==1)
		? ( round(tm) - J_shift )
		: ( floor(tm) - J_shift );

	for (jj1=0; jj1<J1; jj1++) {
		const int k1 = koff + jj1 + 1;
		const int n1 = ncenter + round((tm - k1) * L1);
		register double coefr = r_h1[n1];
		const int kk = (k1 + K1) % K1;

		/* instead of f = h c, we have c += h^* f */
		r_ck[kk] += coefr * fmr;
		i_ck[kk] += coefr * fmi;
	}
    }
}


/*
* Usage: ck = function(fm, h_table, J, L, tm, K)
*/
static int interp1_table_adj_mex(
mxArray *plhs[],
const mxArray *mx_fm,
const mxArray *mx_h1,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm,
const mxArray *mx_K)
{
	const int M = mxGetM(mx_fm);	/* # of time samples */
	const int N = mxGetN(mx_fm);	/* # of realizations */
	const int J = *((int *) mxGetData(mx_J));
	const int K = *((int *) mxGetData(mx_K));
	const int L = *((int *) mxGetData(mx_L));

	const double *r_fm = mxGetPr(mx_fm);
	const double *i_fm = mxGetPi(mx_fm);
	const double *p_tm = mxGetPr(mx_tm);
	const double *r_h1 = mxGetPr(mx_h1);
	double *r_ck, *i_ck;
	int nn;

	if (N != 1)
		fprintf(stderr, "Caution: multiple columns?");

	Call(mxIsComplexDouble, (mx_fm))
	Call(mxIsRealDouble, (mx_tm))

	/* J, L, K must be scalar */
	if (!mxIsScalarInt32(mx_J))
		Fail("J must be scalar int32")
	if (!mxIsScalarInt32(mx_K))
		Fail("K must be scalar int32")
	if (!mxIsScalarInt32(mx_L))
		Fail("L must be scalar int32")

	/* check h table size */
	if ((int) mxGetM(mx_h1) != J*L+1) {
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			J, L, (int) mxGetM(mx_h1));
		Fail("h size problem")
	}
	if (mxGetN(mx_h1) != 1)
		Fail("h must be col vector")

	if (M != (int) mxGetM(mx_tm) || 1 != mxGetN(mx_tm))
		Fail("t_m must be Mx1 col vector")

	/* create a new array and set the output pointer to it */
	plhs[0] = mxCreateDoubleMatrix(K, N, mxCOMPLEX);
	r_ck = mxGetPr(plhs[0]);
	i_ck = mxGetPi(plhs[0]);

	/* call the C subroutine N times; once for each realization */
	if (mxIsComplexDouble(mx_h1)) {
		const double *i_h1 = mxGetPi(mx_h1);
		for (nn=0; nn < N; ++nn) {
			interp1_table_complex_adj(r_ck, i_ck, K,
				r_h1, i_h1, J, L, p_tm, M, r_fm, i_fm);
			r_ck += K; i_ck += K;
			r_fm += M; i_fm += M;
		}
	}

	else if (mxIsRealDouble(mx_h1)) {
		for (nn=0; nn < N; ++nn) {
			interp1_table_real_adj(r_ck, i_ck, K,
				r_h1, J, L, p_tm, M, r_fm, i_fm);
			r_ck += K; i_ck += K;
			r_fm += M; i_fm += M;
		}
	}

	else
		Fail("h must be real or complex double (preferably real)")

	return 1;
}


/*
* The gateway routine.
* Usage: ck = function(fm, h_table, J, L, tm, K)
*/
void mexFunction(
const int nlhs, mxArray *plhs[],
const int nrhs, const mxArray *prhs[])
{
	/* check for the proper number of arguments */
	if (nrhs != 6)
		mexFail("6 inputs needed: (f, h, J, L, t, K)")
	if (nlhs > 1)
		mexFail("Less than one output arguments.")

	if (!interp1_table_adj_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]))
		mexFail("interp1_table_adj_mex() failed")

	return;
}
