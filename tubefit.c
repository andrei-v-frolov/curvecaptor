/* $Id: tubefit.c,v 1.3 2001/09/12 03:37:59 frolov Exp $ */

/*
 * Curve Captor - vacuum tube curve capture and model builder tool
 * Numerical backend - data fitting and approximation routines.
 * 
 * Copyright (C) 2001 Andrei Frolov <andrei@phys.ualberta.ca>
 * Distributed under the terms of GNU Public License.
 * 
 */

#define SELF "tubefit"

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>



/******************* Options and defaults *****************************/

/* Usage */
char *usage_msg[] = {
	"Vacuum tube model builder, Version 0.1",
	"Author: Andrei Frolov <andrei@phys.ualberta.ca>",
	"",
	"Usage: " SELF " -[2|3|4|5] [...] < curve.dat",
	"",
	" General options:",
	"  -h		print this message",
	"  -v		verbosity (cumulative)",
	"",
	" Operation mode:",
	"  -[2|3|4|5]	device type (2=diode, 3=triode, etc)",
	"  -f format	specify format for tagged input data",
	"  -d		dump data in plain format for later use",
	"",
	" Model fit options:",
	"  ---",
	NULL
};


/* Options */
static int verbose = 0;

static int                dtype = 3;
static char             *format = NULL;
static int          output_only = 0;



/**********************************************************************/

/*
 * Numerical and utility functions:
 *  error handlers, vectors and matrices, other misc stuff.
 *
 */

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

/*
 * Non-linear model fitting:
 *   chi^2 optimization via simulated annealing
 * 
 * Tube models:
 *   - Diode:  Child-Langmuir, Perugini
 *   - Triode: Child-Langmuir, Rydel, Scott, Koren
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
	double Vp = V[0];
	double K = p[0];
	
	I = K * pow(uramp(Vp), 1.5);
	
	return I;
}

/* Vacuum diode; Child-Langmuir law with contact potential */
static double diode_cp(double p[], double V[])
{
	double I;
	double Vp = V[0];
	double K = p[0], Vc = p[1];
	
	I = K * pow(uramp(Vp+Vc), 1.5);
	
	return I;
}

/* Vacuum diode; Perugini model */
static double init_perugini[] = {0.0, 0.0, 0.0, 1.5};
static double diode_perugini(double p[], double V[])
{
	double I;
	double Vp = V[0];
	double Ka = p[0], Kb = p[1], Vc = p[2], gamma = p[3];
	
	I = (Ka + Kb*Vp) * pow(uramp(Vp+Vc), gamma);
	
	return I;
}


/* Vacuum triode; Child-Langmuir law */
static double triode(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double K = p[0], mu = p[1];
	
	I = K * pow(uramp(mu*Vg + Vp), 1.5);
	
	return I;
}

/* Vacuum triode; Child-Langmuir law with contact potential */
static double triode_cp(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double K = p[0], mu = p[1], Vc = p[2];
	
	I = K * pow(uramp(mu*Vg + Vp + Vc), 1.5);
	
	return I;
}

/* Vacuum triode; Rydel model (4 parameters) */
static double triode_rydel4(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double Ka = p[0], Kb = p[1], mu = p[2], Vc = p[3];
	
	I = (Ka + Kb*Vg) * pow(uramp(mu*Vg + Vp + Vc), 1.5);
	
	return I;
}

/* Vacuum triode; Rydel model (5 parameters) */
static double triode_rydel5(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double Ka = p[0], Kb = p[1], mu = p[2], Vc = p[3], C = p[4];
	
	I = (Ka + Kb*Vg) * pow(uramp(mu*Vg + Vp + Vc), 1.5) * Vp/(Vp+C);
	
	return I;
}

/* Vacuum triode; Rydel model for grid current */
static double triode_rydel_grid(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double K = p[0], A = p[1], B = p[2];
	
	I = K * pow((A+Vp)/(B+Vp), 4.0) * pow(uramp(Vg), 1.5);
	
	return I;
}

/* Vacuum triode; Scott model */
static double init_scott[] = {0.0, 1.0, 30.0, 1.5};
static double triode_scott(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
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
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kp = p[1], Kv = p[2], mu = p[3], gamma = p[4];
	
	U = Vp * log(1.0 + exp(Kp + Kp*mu*Vg/sqrt(Kv + Vp*Vp)))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}


/******************* Vacuum tube model index **************************/

/* Tube model structure */
typedef struct {
	int dtype;		/* 2=diode, 3=triode, etc */
	char *name;		/* Long model name */
	char *macro;		/* Macro implementing Spice model */
	int params;		/* # of model parameters to fit */
	double (*curve)(double p[], double V[]);
	double *p;		/* Parameter values: initial guess/model fit */
} model;

/* Tube model index */
model mindex[] = {
	/* Vacuum diode models */
	{2, "Child-Langmuir law", "diode", 1, diode},
	{2, "Child-Langmuir law with contact potential", "diode_cp", 2, diode_cp},
	{2, "Perugini model", "perugini", 4, diode_perugini, init_perugini},
	
	/* Vacuum triode models */
	{3, "Child-Langmuir law", "triode", 2, triode},
	{3, "Child-Langmuir law with contact potential", "triode_cp", 3, triode_cp},
	{3, "Rydel model (4 parameters)", "rydel4", 4, triode_rydel4},
	{3, "Rydel model (5 parameters)", "rydel5", 5, triode_rydel5},
	{3, "Rydel model for grid current", "rydelg", 3, triode_rydel_grid},
	{3, "Scott model", "scott", 4, triode_scott, init_scott},
	{3, "Koren model", "koren", 5, triode_koren, init_koren},
};


/******************* Vacuum tube model fit ****************************/

/* Curve data is passed as global variable */
static int _pts_;
static double **_data_;
static model *_model_;

/* Minimization criterion is least square */
static double chi2(double p[])
{
	int i, n = _pts_;
	double t, s = 0.0;
	
	for (i = 0; i < n; i++) {
		double I = _data_[3][i];
		double V[] = {_data_[0][i], _data_[1][i], _data_[2][i]};
		
		t = I - (*(_model_->curve))(p, V); s += t*t;
	}
	
	return s/n;
}

/* Fit parametric curve model to the data using simulated annealing */
double *fit_curve(double **data, int n, model *m)
{
	int i, D = m->params;
	double *p = vector(D+1), **S;
	
	_pts_ = n;
	_data_ = data;
	_model_ = m;
	
	for(i = 0; i < D; i++) p[i] = m->p ? m->p[i] : 0.0;
	
	S = new_amoeba(p, D, chi2, 1.0); anneal(S, D, chi2, 1.0, 1000, 1.0e-6);
	restart_amoeba(S, D, chi2, 0.3); anneal(S, D, chi2, 0.1, 1000, 1.0e-6);
	restart_amoeba(S, D, chi2, 0.1); anneal(S, D, chi2, 0.0, 9999, 1.0e-6);
	
	for (i = 0; i <= D; i++) p[i] = S[D+1][i]; m->p = p;
	
	free_matrix(S);
	
	return p;
}

/* Try all appropriate models */
void try_all_models(double **data, int n)
{
	int i, j;
	double *p;
	
	for (i = 0; i < sizeof(mindex)/sizeof(model); i++) {
		model *m = &(mindex[i]);
		int D = m->params;
		
		if (m->dtype != dtype) continue;
		else p = fit_curve(data, n, m);
		
		printf("* %s: mean fit error %g mA\n", m->name, sqrt(p[D]));
		
		printf("*   %s(", m->macro);
		for (j = 0; j < D; j++)
			printf("%.10g%s", p[j], ((j < D-1) ? "," : ""));
		printf(")\n");
	}
}



/****************** Data I/O routines *********************************/

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

/* Read curve data in tagged format (untranslated axis units) */
double **read_tagged_data(FILE *fp, const char *format, int *pts)
{
	int i;
	double A[2][2], B[2], **m;
	
	char *tag[3] = {NULL, NULL, NULL};
	int c[3], n[3] = {0, 0, 0}, size[3] = {32, 32, 128};
	double **D[3] = {matrix(3, size[0]), matrix(3, size[1]), matrix(3, size[2])};
	
	size_t bsize = 128;
	char l[13]; double t[3];
	char *buffer = (char *)xmalloc(bsize);
	
	
	/* Read curve data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#' || *buffer == '\n') continue;
		if (*buffer == ';') break;
		
		for (i = 0; i < 3; i++) if (n[i] >= size[i]) {
			D[i] = grow_matrix(D[i], 3, size[i]<<=1);
		}
		
		t[0] = t[1] = t[2] = 0.0;
		if (sscanf(buffer, " %12[^ =]=%lf %lf %lf", l, &(t[0]), &(t[1]), &(t[2])) < 4)
			error("syntax error in curve data");
		
		if ((i = pindex(format, l)) == -1) continue;
		
		if (!tag[i]) tag[i] = xstrdup(l);
		
		D[i][0][n[i]] = t[0];
		D[i][1][n[i]] = t[1];
		D[i][2][n[i]] = t[2];
		
		n[i]++;
	}
	free(buffer);
	
	/* Translate data */
	*pts = n[2]; m = matrix(4, *pts);
	axismap(D[0], n[0], D[1], n[1], A, B);
	
	for (i = 0; i < 3; i++)
		if ((c[i] = pindex("Vp Vg Vs Ip|Ig|Is", tag[i])) == -1)
			error("Unknown curve parameter '%s'", tag[i]);
	
	for (i = 0; i < *pts; i++) {
		double x = D[2][1][i] - B[0], y = D[2][2][i] - B[1];
		
		m[0][i] = m[1][i] = m[2][i] = m[3][i] = 0.0;
		
		m[c[0]][i] = A[0][0]*x + A[0][1]*y;
		m[c[1]][i] = A[1][0]*x + A[1][1]*y;
		m[c[2]][i] = D[2][0][i];
	}
	
	free_matrix(D[0]);
	free_matrix(D[1]);
	free_matrix(D[2]);
	
	return m;
}

/* Read curve data in plain format */
double **read_data(FILE *fp, int *pts)
{
	int n = 0, size = 128;
	double **m = matrix(4, size);
	
	size_t bsize = 128;
	char *buffer = (char *)xmalloc(bsize);
	
	/* Read curve data */
	while (getline(&buffer, &bsize, fp) != -1) {
		if (*buffer == '#' || *buffer == '\n') continue;
		if (*buffer == ';') break;
		
		if (n >= size) { m = grow_matrix(m, 4, size<<=1); }
		
		m[0][n] = m[1][n] = m[2][n] = m[3][n] = 0.0;
		
		switch (dtype) {
			case 2:
				if (sscanf(buffer, " %lf %lf",
					&(m[0][n]), &(m[3][n])) < 2)
					error("missing parameters in diode data");
				break;
			case 3:
				if (sscanf(buffer, " %lf %lf %lf",
					&(m[1][n]), &(m[0][n]), &(m[3][n])) < 3)
					error("missing parameters in triode data");
				break;
			case 4:
			case 5:
				if (sscanf(buffer, " %lf %lf %lf %lf",
					&(m[2][n]), &(m[1][n]), &(m[0][n]), &(m[3][n])) < 4)
					error("missing parameters in tetrode/pentode data");
				break;
		}
		
		n++;
	}
	free(buffer);
	
	*pts = n; return m;
}

/* Write curve data in plain format */
void write_data(FILE *fp, double **m, int n)
{
	int i;
	
	for (i = 0; i < n; i++) switch (dtype) {
		case 2:
			fprintf(fp, "%12.10g %12.10g\n", m[0][i], m[3][i]); break;
		case 3:
			fprintf(fp, "%12.10g %12.10g %12.10g\n", m[1][i], m[0][i], m[3][i]); break;
		case 4:
		case 5:
			fprintf(fp, "%12.10g %12.10g %12.10g %12.10g\n", m[2][i], m[1][i], m[0][i], m[3][i]); break;
	}
}



/**********************************************************************/

double *fft256(double F[]);

/* Find a root of f(x)=0 bracketed in interval [x1,x2] */
double zbrent(double (*f)(double), double x1, double x2, double eps)
{
	int i;
	double a = x1, b = x2, c = x2, d, e, min1, min2;
	double fa = (*f)(a), fb = (*f)(b), fc = fb, p, q, r, s, tol, xm;
	
	#define ITMAX 100
	#define EPS 3.0e-8
	
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		error("Root must be bracketed in zbrent()");
	
	for (i = 0; i < ITMAX; i++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a; fc = fa; e = d = b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a = b; fa = fb;
			b = c; fb = fc;
			c = a; fc = fa;
		}
		
		tol = 2.0*EPS*fabs(b) + 0.5*eps; xm=0.5*(c-b);
		
		if (fabs(xm) <= tol || fb == 0.0) return b;
		if (fabs(e) >= tol && fabs(fa) > fabs(fb)) {
			s = fb/fa;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			
			if (p > 0.0) q = -q; p = fabs(p);
			
			min1 = 3.0*xm*q - fabs(tol*q); min2 = fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) { e = d; d = p/q; }
			else { d = xm; e = d; }
		} else {
			d = xm; e = d;
		}
		
		a = b; fa = fb;
		if (fabs(d) > tol) b += d;
		else b += (xm > 0.0) ? fabs(tol) : -fabs(tol);
		fb=(*f)(b);
	}
	
	error("Maximum number of iterations exceeded in zbrent()");
	
	#undef ITMAX
	#undef EPS
	
	return 0.0;
}


/* Working point is passed as global variable */
static double _VI_[6]; /* Vp Vg Vs I (Vb 1/R) */


/* Find grid bias for a given working point */
static double _bias_(double Vg)
{
	_VI_[1] = Vg; return (*(_model_->curve))(_model_->p, _VI_) - _VI_[3];
}

double bias(model *m, double Vp, double Ip)
{
	_model_ = m; _VI_[0] = _VI_[2] = Vp; _VI_[3] = Ip;
	
	return zbrent(_bias_, -300, 100, 1.0e-6);
}

/* Find tube output on a given lloadline */
static double _loadline_(double Vp)
{
	_VI_[0] = _VI_[2] = Vp; return (*(_model_->curve))(_model_->p, _VI_) - _VI_[3] + _VI_[5]*(Vp-_VI_[4]);
}

double output(model *m, double Vg, double Vb, double I0, double R)
{
	_model_ = m; _VI_[1] = Vg; _VI_[3] = I0;
	_VI_[4] = Vb; _VI_[5] = (R != 0.0) ? 1.0/R : 0.0;
	
	return zbrent(_loadline_, 0, 1000, 1.0e-6);
}


/* Return the square amplitude of j-th harmonic */
static double harmonic(double *H, int j)
{
	int i = j<<1; double hj = H[i]*H[i] + H[i+1]*H[i+1];
	
	if (j) { i = (256-j)<<1; hj += H[i]*H[i] + H[i+1]*H[i+1]; hj *= 2.0; }
	
	return hj;
}

/* Analyze distortion spectrum for a loadline */
void distortion(model *m, double Vin, double Vb, double I0, double R)
{
	int i;
	double Vg, *Vp = vector(256);
	double Vbias = bias(m, Vb, I0);
	double *H, A, H2, He = 0.0, Ho = 0.0;
	
	for (i = 0; i < 256; i++) {
		Vg = Vbias + Vin*sin(M_PI*i/128.0);
		Vp[i] = output(m, Vg, Vb, I0, R);
	}
	
	H = fft256(Vp);
	A = sqrt(harmonic(H, 1));
	H2 = sqrt(harmonic(H, 2))/A;
	
	for (i = 128; i > 2; i--) {
		if (i%2) Ho += harmonic(H, i);
		else He += harmonic(H, i);
	}
	He = sqrt(He)/A; Ho = sqrt(Ho)/A;
	
	printf("%g %g %g %g    ", A, H2*100.0, He*100.0, Ho*100.0);
	for (i = 0; i < 10; i++) printf("%g ", sqrt(harmonic(H, i)));
	printf("\n");
	
	free_vector(H);
}

/* Tune loadline */
void loadtune(model *m, double Vin, double Vb, double I0)
{
	double R;
	
	for (R = 0.01; R < 30.0; R *= pow(10.0, 0.1)) {
		printf("%g    ", R); distortion(m, Vin, Vb, I0, R);
	}
}



/**********************************************************************/

int main(int argc, char *argv[])
{
	char c;
	int n; double **d;
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv2345f:d")) != -1)
	switch (c) {
	/* General options */
		case 'h':				/* Help message */
			usage(); break;
		case 'v':				/* Verbosity */
			verbose++; break;
	
	/* Operation mode */
		case '2':				/* Diode */
			dtype = 2; break;
		case '3':				/* Triode */
			dtype = 3; break;
		case '4':				/* Tetrode */
			dtype = 4; break;
		case '5':				/* Pentode */
			dtype = 5; break;
		case 'f':				/* Tagged data format */
			format = optarg; break;
		case 'd':				/* Output data only */
			output_only = 1; break;
	
	/* Default response */
		default:
			usage();
	}
	
	if (argc != optind) usage();
	
	/* Read data */
	if (format) d = read_tagged_data(stdin, format, &n);
	else d = read_data(stdin, &n);
	
	if (output_only) { write_data(stdout, d, n); exit(0); }
	
	/* Do model fit */
	try_all_models(d, n);
	
	//printf("%g\n", bias(mindex+9, 150.0, 100.0));
	//printf("%g\n", output(mindex+9, -0.0, 150.0, 30.0, 1.714));
	loadtune(mindex+9, 50.0, 125.0, 100.0);
	return 0;
}
