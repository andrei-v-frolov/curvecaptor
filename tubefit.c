/* $Id: tubefit.c,v 1.1 2001/09/04 03:32:49 frolov Exp $ */

/*
 * Scanner Calibration Reasonably Easy (scarse)
 * Data fitting and approximation routines.
 * 
 * Copyright (C) 2000 Scarse Project
 * Distributed under the terms of GNU Public License.
 * 
 * Maintainer: Andrei Frolov <andrei@phys.ualberta.ca>
 * 
 */

#define SELF "calibrate"

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <fnmatch.h>
#include <math.h>




/**********************************************************************/

/******************* Options and defaults *****************************/

/* Usage */
char *usage_msg[] = {
	"Input/output device calibrator, Version 0.1",
	"Author: Andrei Frolov <andrei@phys.ualberta.ca>",
	"",
	"Usage: " SELF " -d [scanner|display|printer] [...] profile.icm",
	"For now, only scanner calibration is fully operational.",
	"",
	" General options:",
	"  -h		print this message",
	"  -v		verbosity (cumulative)",
	"  --longopt	expand option macro 'longopt' as defined in etc/" SELF ".options",
	"",
	" ICC profile generation:",
	"  -D		dump measured data to stdout, instead of running ipb",
	"  -O options	pass options to ICC profile builder, see ipb help",
	"",
	" Calibration target options:",
	"  -l		list all known target types (or batches if type was given by -t)",
	"  -t type	specify target type (or IT8.7 target layout file to parse)",
	"  -b batch	specify target batch (or IT8.7 target data file to parse)",
	"  -g geometry	size and position of target grid (geometry is like in X)",
	"  -i file	import calibration target data from raster image",
	"  -r file	render calibration target data into raster image and exit",
	NULL
};



/**********************************************************************/

#define TINY 1.0e-8

#define ppow(x,g) ((x > 0.0) ? pow(x, g) : -pow(-x, g))


/******************* Basic error handlers *****************************/

void usage()
{
	int i = 0; char *s;
	
	while ((s = usage_msg[i++]))
		fprintf(stderr, "%s\n", s);
	
	exit(1);
}

void warning(char *fmt, ...)
{
	va_list args;
	
	fprintf(stderr,"%s: ", SELF);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}

void error(char *fmt, ...)
{
	va_list args;
	
	fprintf(stderr, "%s: ", SELF);
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	
	exit(-1);
}

void fatal(char *msg)
{
	error("Fatal error: %s\nTerminating process...", msg);
}


/******************* Foo or die routines ******************************/

/* allocate memory or die */
void *xmalloc(size_t size)
{
	register void *p = malloc(size);
	
	if (!p) fatal("virtual memory exhausted in malloc()");
	return p;
}

/* reallocate memory or die */
void *xrealloc(void *addr, size_t size)
{
	register void *p = realloc(addr, size);
	
	if (!p) fatal("virtual memory exhausted in realloc()");
	return p;
}

/* duplicate string or die */
char *xstrdup(const char *s)
{
	register char *p = strdup(s);
	
	if (!p) fatal("virtual memory exhausted in strdup()");
	return p;
}


/******************* Vectors and matrices *****************************/

/* allocate a vector with subscript range v[0..n] */
double *vector(unsigned long n)
{
	register double *v = (double *)xmalloc((size_t)(n*sizeof(double)));
	
	return v;
}

/* grow a vector to subscript range v[0..n] */
double *grow_vector(double *v, unsigned long n)
{
	v = (double *)xrealloc(v, (size_t)(n*sizeof(double)));
	
	return v;
}

/* allocate an int vector with subscript range v[0..n] */
int *ivector(unsigned long n)
{
	register int *v = (int *)xmalloc((size_t)(n*sizeof(int)));
	
	return v;
}

/* grow an int vector to subscript range v[0..n] */
int *grow_ivector(int *v, unsigned long n)
{
	v = (int *)xrealloc(v, (size_t)(n*sizeof(int)));
	
	return v;
}

/* allocate an unsigned long vector with subscript range v[0..n] */
unsigned long *uvector(unsigned long n)
{
	register unsigned long *v = (unsigned long *)xmalloc((size_t)(n*sizeof(unsigned long)));
	
	return v;
}

/* grow an unsigned long vector to subscript range v[0..n] */
unsigned long *grow_uvector(unsigned long *v, unsigned long n)
{
	v = (unsigned long *)xrealloc(v, (size_t)(n*sizeof(unsigned long)));
	
	return v;
}

/* free a vector allocated with vector() */
void free_vector(void *v)
{
	free((void *)(v));
}


/* copy vector data with subscript range [0..n] */
void vcopy(double src[], double dest[], unsigned long n)
{
	register unsigned long i;
	
	for (i = 0; i < n; i++) dest[i] = src[i];
}

/* apply gamma transformation to a vector */
void vgamma(double src[], double dest[], unsigned long n, double gamma)
{
	register unsigned long i;
	
	for (i = 0; i < n; i++) dest[i] = ppow(src[i], gamma);
}


/* allocate a matrix with subscript range m[0..nr][0..nc] */
double **matrix(unsigned long nr, unsigned long nc)
{
	register unsigned long i;
	register double **m = (double **)xmalloc((size_t)(nr*sizeof(double *)+sizeof(unsigned long)));
	
	*(((unsigned long *)(m))++) = nr;
	
	for (i = 0; i < nr; i++) m[i] = vector(nc);
	
	return m;
}

/* grow a matrix to subscript range m[0..nr][0..nc] */
double **grow_matrix(double **m, unsigned long nr, unsigned long nc)
{
	register unsigned long i;
	unsigned long old_nr = *(--((unsigned long *)(m)));
	
	/* Reallocate row index if necessary */
	if (nr != old_nr)
		m = (double **)xrealloc(m, (size_t)(nr*sizeof(double *)+sizeof(unsigned long)));
	
	*(((unsigned long *)(m))++) = nr;
	
	/* Reallocate rows */
	for (i = 0; i < old_nr; i++) m[i] = grow_vector(m[i], nc);
	for (i = old_nr; i < nr; i++) m[i] = vector(nc);
	
	return m;
}

/* free a matrix allocated by matrix() */
void free_matrix(double **m)
{
	register unsigned long i;
	unsigned long nr = *((unsigned long *)(m)-1);
	
	for (i = 0; i < nr; i++) free_vector(m[i]);
	
	free((void *)((unsigned long *)(m)-1));
}










/**********************************************************************/

/* Linear regression: least square fit of a straight line y=ax+b*/
void linregr(double x[], double y[], int n, double *a, double *b)
{
	int i;
	double xa, ya, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
	
	for (i = 0; i < n; i++) {
		sx += x[i]; sy += y[i];
	}
	xa = sx/n; ya = sy/n;
	
	for (i = 0; i < n; i++) {
		sxx += (x[i]-xa)*x[i];
		sxy += (x[i]-xa)*y[i];
	}
	
	*a = sxy/sxx; *b = ya - (*a)*xa;
}

/* Get coordinate map for traced axis positions */
void axismap(double **X, int nx, double **Y, int ny, double A[2][2], double B[2])
{
	double det, M[2][2], T[2][2];
	
	linregr(X[0], X[1], nx, &(M[0][0]), &(T[0][0]));
	linregr(X[0], X[2], nx, &(M[1][0]), &(T[1][0]));
	linregr(Y[0], Y[1], ny, &(M[0][1]), &(T[0][1]));
	linregr(Y[0], Y[2], ny, &(M[1][1]), &(T[1][1]));
	
	det = M[0][0]*M[1][1]-M[0][1]*M[1][0];
	A[0][0] =  M[1][1]/det;
	A[0][1] = -M[0][1]/det;
	A[1][0] = -M[1][0]/det;
	A[1][1] =  M[0][0]/det;
	
	if (fabs(M[0][0]) >= fabs(M[1][0])) B[0] = T[0][0]; else B[0] = T[0][1];
	if (fabs(M[1][1]) >= fabs(M[0][1])) B[1] = T[1][1]; else B[1] = T[1][0];
}

/* Find parameter index in a template */
int pindex(const char *t, const char *p)
{
	int i = 0;
	char *q = strstr(t, p);
	
	if (!q) return -1;
	while (q-- > t) { if (*q == ' ') i++; }
	
	return i;
}

/* Read curve data */
double **read_curve(FILE *fp, int *pts)
{
	int i;
	double A[2][2], B[2], **m;
	int n[3] = {0, 0, 0}, size[3] = {32, 32, 128};
	double **D[3] = {matrix(3, size[0]), matrix(3, size[1]), matrix(3, size[2])};
	
	char l[13];
	double t[3];
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	
	
	/* Read curve data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#') continue;
		if (*buffer == ';') break;
		
		for (i = 0; i < 3; i++) if (n[i] >= size[i]) {
			D[i] = grow_matrix(D[i], 3, size[i]<<=1);
		}
		
		t[0] = t[1] = t[2] = 0.0;
		if (sscanf(buffer, "%12[^ =]=%lf %lf %lf", &l, &(t[0]), &(t[1]), &(t[2])) < 4)
			error("syntax error in curve data");
		
		i = pindex("Vp Ip Vg", l); if (i == -1)
			error("Unknown curve parameter '%s'", l);
		
		D[i][0][n[i]] = t[0];
		D[i][1][n[i]] = t[1];
		D[i][2][n[i]] = t[2];
		
		n[i]++;
	}
	free(buffer);
	
	/* Translate data */
	axismap(D[0], n[0], D[1], n[1], A, B);
	
	*pts = n[2]; m = matrix(3, *pts);
	for (i = 0; i < *pts; i++) {
		double x = D[2][1][i] - B[0], y = D[2][2][i] - B[1];
		
		m[0][i] = D[2][0][i];
		m[1][i] = A[0][0]*x + A[0][1]*y;
		m[2][i] = A[1][0]*x + A[1][1]*y;
	}
	
	free_matrix(D[0]);
	free_matrix(D[1]);
	free_matrix(D[2]);
	
	return m;
}


















/**********************************************************************/

/*
 * Non-linear model fitting:
 *   chi^2 optimization via simulated annealing
 * 
 * Tube models:
 *   - power-law curve with clipping
 *   - 3x3 matrix transform with clipping
 * 
 */


/******************* Simulated annealing ******************************/

/* n-dimensional amoeba matrix layout:
 *
 *                    |   [0..n-1]   |  [n]  
 *             -------+--------------+-------
 *             [0..n] |simplex coords| f(x)  
 *             -------+--------------+-------
 *  S[i][j] =   [n+1] |best pt coords| f(b)  
 *              [n+2] |vertex sum    |   0   
 *              [n+3] |tmp           |   .   
 *             -------+--------------+-------
 */

/* Restart amoeba around best point */
void restart_amoeba(double **S, int n, double (*func)(double []), double lambda)
{
	int i, j;
	double *best = S[n+1], *sum = S[n+2];
	
	for (i = 0; i <= n; i++)
		for (j = 0; j <= n; j++)
			S[i][j] = best[j];
	
	for (i = 0; i < n; i++) {
		S[i][i] += lambda;
		S[i][n] = (*func)(S[i]);
	}
	S[n][n] = (*func)(S[n]);
	
	for (i = 0; i < n; i++)
		sum[i] = (n+1) * best[j] + lambda;
	sum[n] = 0.0;
}

/* Create new amoeba */
double **new_amoeba(double x[], int n, double (*func)(double []), double lambda)
{
	int i;
	double **S = matrix(n+4, n+1), *best = S[n+1];
	
	for (i = 0; i < n; i++)
		best[i] = x[i];
	best[n] = HUGE;
	
	restart_amoeba(S, n, func, lambda);
	
	return S;
}


/* Try to crawl in given direction */
static double crawl(double **S, int n, double (*func)(double []), double T, int dir, double amount)
{
	int i;
	double y, *x = S[dir], *best = S[n+1], *sum = S[n+2], *try = S[n+3];
	
	/* Give it a try... */
	for (i = 0; i < n; i++)
		try[i] = (1.0-amount) * (sum[i] - x[i])/n + amount * x[i];
	try[n] = (*func)(try);
	
	/* Best move ever? */
	if (try[n] < best[n])
		for (i = 0; i <= n; i++) best[i] = try[i];
	
	y = try[n] - T * (1.0 + fabs(best[n]))/100.0 * log(RAND_MAX/rand());
	
	/* Favourable move? */
	if (y < x[n]) {
		for (i = 0; i < n; i++) {
			sum[i] += try[i] - x[i];
			x[i] = try[i];
		}
		
		x[n] = try[n];
	}
	
	return y;
}

/* Shrink the amoeba around given vertex */
static void shrink(double **S, int n, double (*func)(double []), double T, int dir, double amount)
{
	int i, j;
	double *x = S[dir], *best = S[n+1], *sum = S[n+2];
	
	/* Shrink the amoeba */
	for (i = 0; i <= n; i++) if (i != dir) {
		for (j = 0; j < n; j++)
			S[i][j] = (1.0-amount) * x[j] + amount * S[i][j];
		S[i][n] = (*func)(S[i]);
		
		if (S[i][n] < best[n])
			for (j = 0; j <= n; j++) best[j] = S[i][j];
	}
	
	/* Update vertex sum */
	for (i = 0; i < n; i++) {
		sum[i] = 0.0;
		
		for (j = 0; j <= n; j++)
			sum[i] += S[j][i];
	}
}

/* Minimize the function by simulated annealing */
void anneal(double **S, int n, double (*func)(double []), double T0, int maxsteps, double tol)
{
	int i, k, lo, hi;
	double T, y, ylo, yhi, yhi2;
	
	for (k = 0; k < maxsteps; ) {
		/* Cooling schedule */
		T = T0 * exp(-log(100.0)*k/maxsteps);
		
		/* Rank simplex vertices */
		lo = hi = 0;
		ylo = HUGE; yhi = yhi2 = -HUGE;
		
		for (i = 0; i <= n; i++) {
			y = S[i][n] + T * (1.0 + fabs(S[i][n]))/100.0 * log(RAND_MAX/rand());
			
			if (y < ylo) { lo = i; ylo = y; }
			if (y > yhi) { yhi2 = yhi; hi = i; yhi = y; }
			else if (y > yhi2) { yhi2 = y; }
		}
		
		/* Are we done yet? */
		if (2.0*fabs(S[hi][n]-S[lo][n])/(fabs(S[hi][n])+fabs(S[lo][n])) < tol) break;
		
		/* Make a move: try reflect first */
		y = crawl(S, n, func, T, hi, -1.0); k++;
		
		if (y <= ylo) {
			/* was good, try expanding */
			y = crawl(S, n, func, T, hi, 2.0); k++;
		} else if (y >= yhi2 ) {
			/* no good, try contracting */
			y = crawl(S, n, func, T, hi, 0.5); k++;
			
			/* if that didn't work, try shrinking */
			if (y >= yhi2) { shrink(S, n, func, T, lo, 0.5); k += n; }
		}
	}
}


/******************* Non-linear curve fit *****************************/

/* Tube model structure */
typedef struct {
	char *name;
	int params;
	double (*curve)(double p[], double V[]);
	double *guess;
} model;

/* Curve data is passed as global variable */
static int _pts_;
static double _eps2_;
static double **_data_;
static model *_model_;

/* Minimization criterion is least square */
static double chi2(double p[])
{
	int i, n = _pts_;
	double t, s = 0.0;
	
	for (i = 0; i < n; i++) {
		double I = _data_[2][i];
		double V[] = {_data_[0][i], _data_[1][i]};
		
		t = I - (*(_model_->curve))(p, V); s += t*t;
	}
	
	return s/n/_eps2_;
}

/* Fit parametric curve model to the data using simulated annealing */
double *fit_curve(double **data, int n, double eps, model *m)
{
	int i, D = m->params;
	double *p = vector(D+1), **S;
	
	_pts_ = n;
	_data_ = data;
	_eps2_ = eps*eps;
	_model_ = m;
	
	for(i = 0; i < D; i++) p[i] = m->guess ? m->guess[i] : 0.0;
	
	S = new_amoeba(p, D, chi2, 1.0); anneal(S, D, chi2, 1.0, 1000, 1.0e-6);
	restart_amoeba(S, D, chi2, 0.3); anneal(S, D, chi2, 0.1, 1000, 1.0e-6);
	restart_amoeba(S, D, chi2, 0.1); anneal(S, D, chi2, 0.0, 9999, 1.0e-6);
	
	for (i = 0; i <= D; i++) {
		p[i] = S[D+1][i]; printf("%g\n", p[i]);
	}
	free_matrix(S);
	
	return p;
}


/******************* Vacuum tube models *******************************/

/* uramp() like in Spice 3f4 */
static double uramp(double x)
{
	return (x > 0.0) ? x : 0.0;
}

/* Vacuum diode; Child-Langmuir law */
static double diode(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double K = p[0];
	
	I = K * pow(uramp(Vp), 1.5);
	
	return I;
}

/* Vacuum diode; Child-Langmuir law with contact potential */
static double diode_cp(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double K = p[0], Vc = p[1];
	
	I = K * pow(uramp(Vp+Vc), 1.5);
	
	return I;
}

/* Vacuum diode; Perugini model */
static double init_perugini[] = {0.0, 0.0, 0.0, 1.5};
static double diode_perugini(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double Ka = p[0], Kb = p[1], Vc = p[2], gamma = p[3];
	
	I = (Ka + Kb*Vp) * pow(uramp(Vp+Vc), gamma);
	
	return I;
}

/* Vacuum triode; Child-Langmuir law */
static double triode(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double K = p[0], mu = p[1];
	
	I = K * pow(uramp(mu*Vg + Vp), 1.5);
	
	return I;
}

/* Vacuum triode; Child-Langmuir law with contact potential */
static double triode_cp(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double K = p[0], mu = p[1], Vc = p[2];
	
	I = K * pow(uramp(mu*Vg + Vp + Vc), 1.5);
	
	return I;
}

/* Vacuum triode; Rydel model (4 parameters) */
static double triode_rydel4(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double Ka = p[0], Kb = p[1], mu = p[2], Vc = p[3];
	
	I = (Ka + Kb*Vg) * pow(uramp(mu*Vg + Vp + Vc), 1.5);
	
	return I;
}

/* Vacuum triode; Rydel model (5 parameters) */
static double triode_rydel5(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double Ka = p[0], Kb = p[1], mu = p[2], Vc = p[3], C = p[4];
	
	I = (Ka + Kb*Vg) * pow(uramp(mu*Vg + Vp + Vc), 1.5) * Vp/(Vp+C);
	
	return I;
}

/* Vacuum triode; Rydel model for grid current */
static double triode_rydel_grid(double p[], double V[])
{
	double I;
	double Vg = V[0], Vp = V[1];
	double K = p[0], A = p[1], B = p[2];
	
	I = K * pow((A+Vp)/(B+Vp), 4.0) * pow(uramp(Vg), 1.5);
	
	return I;
}

/* Vacuum triode; Scott model */
static double init_scott[] = {0.0, 1.0, 30.0, 1.5};
static double triode_scott(double p[], double V[])
{
	double I, U;
	double Vg = V[0], Vp = V[1];
	double K = p[0], Kp = p[1], mu = p[2], gamma = p[3];
	
	U = log(1.0 + exp(Kp * (mu*Vg + Vp)))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}

/* Vacuum triode; Koren model */
static double init_koren[] = {0.0, 1.0, 0.0, 30.0, 1.5};
static double triode_koren(double p[], double V[])
{
	double I, U;
	double Vg = V[0], Vp = V[1];
	double K = p[0], Kp = p[1], Kv = p[2], mu = p[3], gamma = p[4];
	
	U = Vp * log(1.0 + exp(Kp + Kp*mu*Vg/sqrt(Kv + Vp*Vp)))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}




/* Tube model index */
model mindex[] = {
	/* Vacuum diode models */
	{"diode; Child-Langmuir law", 1, diode},
	{"diode; Child-Langmuir law with contact potential", 2, diode_cp},
	{"diode; Perugini model", 4, diode_perugini, init_perugini},
	
	/* Vacuum triode models */
	{"triode; Child-Langmuir law", 2, triode},
	{"triode; Child-Langmuir law with contact potential", 3, triode_cp},
	{"triode; Rydel model (4 parameters)", 4, triode_rydel4},
	{"triode; Rydel model (5 parameters)", 5, triode_rydel5},
	{"triode; Rydel model for grid current", 3, triode_rydel_grid},
	{"triode; Scott model", 4, triode_scott, init_scott},
	{"triode; Koren model", 5, triode_koren, init_koren},
};












/**********************************************************************/

/* Main routine */
int main(int argc, char *argv[])
{
	char c;
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv")) != -1)
	switch (c) {
	/* General options */
		case 'h':				/* Help message */
			usage(); break;
		//case 'v':				/* Verbosity */
		//	verbose++; break;
	
	/* Default response */
		default:
			usage();
	}
	
	if (argc != optind) usage();
	
	{
		int i, n;
		double **m = read_curve(stdin, &n);
		
		for (i = 0; i < 10; i++) {
			printf("%s:\n", mindex[i].name);
			fit_curve(m, n, 0.2, mindex+i);
			printf("\n");
		}
	}
	
	return 0;
}
