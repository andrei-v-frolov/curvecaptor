/* $Id: tubefit.c,v 1.7 2001/09/20 04:26:35 frolov Exp $ */

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
	" Operation parameters:",
	"  -[2|3|4|5]	valve type (2=diode, 3=triode, etc)",
	"  -P Pa	rated anode dissipation power",
	"  -L Vp,Ip[,R]	loadline: working point and load resistance",
	"",
	" Model fit options:",
	"  -C flags	apply cuts to tube data to fit specific requirements",
	"		(g = negative grid; p = rated power; l = loadline)",
	"",
	" I/O functions:",
	"  -f format	specify format for tagged input data",
	"  -d		dump data in plain format for later use",
	NULL
};


/* Options */
static int verbose = 0;

static int                vtype = 3;
static double                Pa = 0.0;
static double                V0 = 0.0;
static double                I0 = 0.0;
static double                RL = 0.0;

static char             *format = NULL;
static int          output_only = 0;

static int             grid_cut = 0;
static int            power_cut = 0;
static int         loadline_cut = 0;



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

/* A way to refer to previous model's parameters */
#define PRIOR_MAGIC 0x00505259
#define PRMREF(m,k) ((((PRIOR_MAGIC<<4) + (m))<<4) + k)
#define PRIOR(k) PRMREF(1,k)


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
static double init_diode_cp[] = {PRIOR(0), 0.0};
static double diode_cp(double p[], double V[])
{
	double I;
	double Vp = V[0];
	double K = p[0], Vc = p[1];
	
	I = K * pow(uramp(Vp+Vc), 1.5);
	
	return I;
}

/* Vacuum diode; Perugini model */
static double init_perugini[] = {PRIOR(0), 0.0, PRIOR(1), 1.5};
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
static double init_triode_cp[] = {PRIOR(0), PRIOR(1), 0.0};
static double triode_cp(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double K = p[0], mu = p[1], Vc = p[2];
	
	I = K * pow(uramp(mu*Vg + Vp + Vc), 1.5);
	
	return I;
}

/* Vacuum triode; Rydel model (4 parameters) */
static double init_rydel4[] = {PRIOR(0), 0.0, PRIOR(1), PRIOR(2)};
static double triode_rydel4(double p[], double V[])
{
	double I;
	double Vp = V[0], Vg = V[1];
	double Ka = p[0], Kb = p[1], mu = p[2], Vc = p[3];
	
	I = (Ka + Kb*Vg) * pow(uramp(mu*Vg + Vp + Vc), 1.5);
	
	return I;
}

/* Vacuum triode; Rydel model (5 parameters) */
static double init_rydel5[] = {PRIOR(0), PRIOR(1), PRIOR(2), PRIOR(3), 0.0};
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

/* Vacuum triode; Koren model (4 parameters) */
static double init_koren4[] = {PRMREF(5,0), 5.0, PRMREF(5,1), 1.5};
static double triode_koren4(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kp = p[1], mu = p[2], gamma = p[3];
	
	U = Vp * log(1.0 + exp(Kp + Kp*mu*Vg/Vp))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}

/* Vacuum triode; Koren model (5 parameters) */
static double init_koren5[] = {PRIOR(0), PRIOR(1), 0.0, PRIOR(2), PRIOR(3)};
static double triode_koren5(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kp = p[1], Kv = p[2], mu = p[3], gamma = p[4];
	
	U = Vp * log(1.0 + exp(Kp + Kp*mu*Vg/sqrt(1000.0*Kv + Vp*Vp)))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}


/* Vacuum triode; modified Koren model (6 parameters) */
static double init_koren6[] = {PRMREF(2,0), PRMREF(2,1), 0.0, PRMREF(2,2), 0.0, PRMREF(2,3)};
static double triode_koren6(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kp = p[1], Kg = p[2], mu = p[3], nu = p[4], gamma = p[5];
	
	U = Vp * log(1.0 + exp(Kp + Kp*(mu+nu*Vg/1000.0)*Vg/Vp))/Kp + Kg*Vg;
	I = K * pow(uramp(U), gamma);
	
	return I;
}


/******************* Vacuum tube model index **************************/

/* Tube model structure */
typedef struct {
	int vtype;		/* 2=diode, 3=triode, etc */
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
	{2, "Child-Langmuir law with contact potential", "diode_cp", 2, diode_cp, init_diode_cp},
	{2, "Perugini model", "perugini", 4, diode_perugini, init_perugini},
	
	/* Vacuum triode models */
	{3, "Child-Langmuir law", "triode", 2, triode},
	{3, "Child-Langmuir law with contact potential", "triode_cp", 3, triode_cp, init_triode_cp},
	{3, "Rydel model (4 parameters)", "rydel4", 4, triode_rydel4, init_rydel4},
	{3, "Rydel model (5 parameters)", "rydel5", 5, triode_rydel5, init_rydel5},
	{3, "Rydel model for grid current", "rydelg", 3, triode_rydel_grid},
	{3, "Koren model (4 parameters)", "koren4", 4, triode_koren4, init_koren4},
	{3, "Koren model (5 parameters)", "koren5", 5, triode_koren5, init_koren5},
	{3, "Modified Koren model (6 parameters)", "koren6", 6, triode_koren6, init_koren6},
};


/******************* Vacuum tube model fit ****************************/

/* Curve data is passed as global variable */
static int _pts_;
static double **_data_;
static double *_weight_;
static model *_model_;

/* Minimization criterion is weighted least square */
static double chi2(double p[])
{
	int i, n = _pts_;
	double t, s = 0.0, norm = 0.0;
	
	for (i = 0; i < n; i++) {
		double w = _weight_[i];
		double I = _data_[3][i];
		double V[] = {_data_[0][i], _data_[1][i], _data_[2][i]};
		
		t = I - (*(_model_->curve))(p, V); s += w*t*t; norm += w;
	}
	
	return s/norm;
}

/* Fit parametric curve model to the data using simulated annealing */
double *fit_curve(double **data, int n, model *m)
{
	int i, j, k, D = m->params;
	double *w = vector(n), *p = vector(D+1), **S;
	
	_pts_ = n;
	_data_ = data;
	_weight_ = w;
	_model_ = m;
	
	/* initialize point weights */
	for (i = 0; i < n; i++) {
		w[i] = 1.0;
		
		if (grid_cut) {
			double V = data[1][i];
			
			if (V > 0.0) w[i] = 0.0;
		}
		
		if (power_cut) {
			double V = data[0][i], I = data[3][i];
			double P = V*I/1000.0;
			
			w[i] *= (1.0 - tanh(16.0 * (P - 1.1*Pa)/Pa))/2.0;
		}
		
		if (loadline_cut) {
			double V = data[0][i], I = data[3][i];
			double G = (RL != 0.0) ? 1000.0/RL : 0.0;
			double dI = (I-I0) + G * (V-V0);
			
			w[i] *= exp(-16.0 * dI*dI/(I0*I0));
		}
	}
	
	/* initialize parameter vector */
	for (i = 0; i < D; i++) {
		if (!m->p) { p[i] = 0.0; continue; }
		if (((int)(m->p[i]) >> 8) != PRIOR_MAGIC) { p[i] = m->p[i]; continue; }
		
		j = ((int)(m->p[i]) & 0xF0) >> 4;
		k = (int)(m->p[i]) & 0x0F;
		p[i] = (m-j)->p[k];
	}
	
	/* optimize */
	S = new_amoeba(p, D, chi2, 1.0); anneal(S, D, chi2, 100.0, 8000, 1.0e-6);
	restart_amoeba(S, D, chi2, 0.3); anneal(S, D, chi2,  10.0, 4000, 1.0e-6);
	restart_amoeba(S, D, chi2, 0.1); anneal(S, D, chi2,   0.0, 2000, 1.0e-6);
	
	for (i = 0; i <= D; i++) p[i] = S[D+1][i]; m->p = p;
	
	free_vector(w);
	free_matrix(S);
	
	return p;
}

/* Try all appropriate models and return the best one */
model *best_model(double **data, int n)
{
	int i, j, best = 0;
	double *p, min = HUGE;
	
	for (i = 0; i < sizeof(mindex)/sizeof(model); i++) {
		model *m = &(mindex[i]);
		int D = m->params;
		
		if (m->vtype != vtype) continue;
		else p = fit_curve(data, n, m);
		
		if (p[D] < min) { best = i; min = p[D]; }
		
		printf("* %s: mean fit error %g mA\n", m->name, sqrt(p[D]));
		
		printf("*   %s(", m->macro);
		for (j = 0; j < D; j++)
			printf("%.10g%s", p[j], ((j < D-1) ? "," : ""));
		printf(")\n");
	}
	
	return &(mindex[best]);
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
		
		switch (vtype) {
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
	
	for (i = 0; i < n; i++) switch (vtype) {
		case 2:
			fprintf(fp, "%12.10g %12.10g\n", m[0][i], m[3][i]); break;
		case 3:
			fprintf(fp, "%12.10g %12.10g %12.10g\n", m[1][i], m[0][i], m[3][i]); break;
		case 4:
		case 5:
			fprintf(fp, "%12.10g %12.10g %12.10g %12.10g\n", m[2][i], m[1][i], m[0][i], m[3][i]); break;
	}
}


/* */
void plot_curves(FILE *fp, double **data, int n, model *m, double Vmax, double Imax, double Vgm, double Vgs)
{
	int i; double Vp, Vg, Ip;
	
	#define X(V) (40.0 + 560.0*(V)/Vmax)
	#define Y(I) (460.0 - 440.0*(I)/Imax)
	
	
	/* Auto ranges */
	if ((Vmax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[0][i] > Vmax) Vmax = data[0][i];
		i = 5.0 * pow(10.0, floor(log10(Vmax))-1); Vmax = i * ceil(Vmax/i);
	}
	
	if ((Imax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[3][i] > Imax) Imax = data[3][i];
		i = 5.0 * pow(10.0, floor(log10(Imax))-1); Imax = i * ceil(Imax/i);
	}
	
	if ((Vgm == 0.0) && data) {
		for (i = 0; i < n; i++) {
			if (fabs(data[1][i]) > fabs(Vgm)) Vgm = data[1][i];
			if ((Vgs == 0.0) && (data[1][i] != 0.0)) Vgs = data[1][i];
		}
	}
	
	
	/* Axis grid */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
			X(0),Y(0), X(Vmax),Y(0), X(Vmax),Y(Imax), X(0),Y(Imax));
	
	for (i = 0; i <= 10; i++) {
		fprintf(fp, "line %g %g %g %g -fill black -width 1 -dash .\n",
				X(Vmax*i/10), Y(0), X(Vmax*i/10), Y(Imax));
		fprintf(fp, "text %g %g -anchor n -text %g -fill black\n",
				X(Vmax*i/10), Y(0)+5, Vmax*i/10);
		fprintf(fp, "line %g %g %g %g -fill black -width 1 -dash .\n",
				X(0), Y(Imax*i/10), X(Vmax), Y(Imax*i/10));
		fprintf(fp, "text %g %g -anchor e -text %g -fill black\n",
				X(0)-5, Y(Imax*i/10), Imax*i/10);
	}
	
	/* Data points */
	if (data) for (i = 0; i < n; i++) {
		Vp = data[0][i]; Vg = data[1][i]; Ip = data[3][i];
		
		if (fabs(Vg/Vgs - rint(Vg/Vgs)) < 1.0e-4)
			fprintf(fp, "oval %g %g %g %g -outline black -width 1\n",
					X(Vp)-1, Y(Ip)+1, X(Vp)+2, Y(Ip)-2);
	}
	
	/* Plate curves */
	if (m) for (Vg = 0.0; fabs(Vg) <= fabs(Vgm); Vg += Vgs) {
		fprintf(fp, "line ");
		
		for (Vp = 0.0, Ip = 0.0; (Vp <= Vmax) && (Ip <= Imax); Vp += Vmax/200) {
			double V[] = {Vp, Vg, Vp};
			
			Ip = (*(m->curve))(m->p, V);
			fprintf(fp, "%g %g ", X(Vp), Y(Ip));
		}
		
		fprintf(fp, "-smooth 1 -fill black -width 2\n");
		fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline white -fill white -width 1\n",
				X(0),Y(Imax)-1, X(Vmax),Y(Imax)-1, X(Vmax),0.0, X(0),0.0);
		fprintf(fp, "text %g %g -anchor nw -text %g -fill black\n",
				X(Vp)+2, Y((Ip > Imax) ? Imax : Ip)+2, Vg);
	}
	
	/* Power */
	if (Pa > 0.0) {
		fprintf(fp, "line ");
		
		for (Vp = 1000.0*Pa/Imax; Vp <= Vmax; Vp += Vmax/200)
			fprintf(fp, "%g %g ", X(Vp), Y(1000.0*Pa/Vp));
		
		fprintf(fp, "-smooth 1 -fill red -width 3\n");
		fprintf(fp, "text %g %g -anchor s -text \"%g W\" -fill red\n",
				X(0.95*Vmax), Y(1053.0*Pa/Vmax)-5, Pa);
	}
	
	/* Loadline */
	if (V0 > 0.0) {
		double G = (RL != 0.0) ? 1000.0/RL : 0.0;
		double V1 = 0.0, I1 = I0+G*V0, V2 = Vmax, I2 = I0+G*(V0-Vmax);
		
		if ((I1 > Imax) && (G != 0.0)) { I1 = Imax; V1 = V0 - (Imax-I0)/G; }
		if ((I2 < 0.0) && (G != 0.0)) { I2 = 0.0; V2 = V0 + I0/G; }
		
		fprintf(fp, "line ");
		fprintf(fp, "%g %g %g %g ", X(V1), Y(I1), X(V2), Y(I2));
		fprintf(fp, "-fill blue -width 3\n");
		
		fprintf(fp, "text %g %g -anchor sw -text \"%g R\" -fill blue\n",
				X(V2)+2, Y(I2)-5, RL);
	}
	
	#undef X
	#undef Y
}



/****************** Bias and loadlines ********************************/

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


/* Working point and other stuff is passed in global variables */
static double _Vp_, _Vg_, _I_, _Vb_, _G_, _Vo_;


/* Find grid bias for a given working point */
static double _bias_(double Vg)
{
	double V[] = {_Vp_, Vg, _Vp_}, I = (*(_model_->curve))(_model_->p, V);
	
	return I - _I_;
}

double bias(model *m, double Vp, double Ip)
{
	_model_ = m; _Vp_ = Vp; _I_ = Ip;
	
	return zbrent(_bias_, -1000, 500, 1.0e-12);
}


/* Find tube output on a given loadline */
static double _loadline_(double Vp)
{
	double V[] = {Vp, _Vg_, Vp}, I = (*(_model_->curve))(_model_->p, V);
	
	return (I-_I_) + _G_*(Vp-_Vb_);
}

double output(model *m, double Vg, double Vb, double I0, double R)
{
	_model_ = m; _Vg_ = Vg; _I_ = I0;
	_Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return zbrent(_loadline_, 0, 3000, 1.0e-12);
}


/* Find grid drive for a given output voltage swing */
static double _swing_(double Vi)
{
	double Vg = _Vg_, Vmin, Vmax;
	
	_Vg_ = Vg+Vi; Vmin = zbrent(_loadline_, 0, 3000, 1.0e-12);
	_Vg_ = Vg-Vi; Vmax = zbrent(_loadline_, 0, 3000, 1.0e-12);
	
	_Vg_ = Vg; return fabs(Vmax-Vmin) - 2.0*_Vo_;
}

double drive(model *m, double Vo, double Vb, double I0, double R)
{
	_model_ = m; _Vg_ = bias(m, Vb, I0); _Vo_ = Vo;
	_I_ = I0; _Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return zbrent(_swing_, 0, 1000, 1.0e-12);
}


/****************** Distortion analysis *******************************/

double *fft256(double F[]);

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


/* Tune load for a fixed tube bias */
void loadtune(model *m, double V, double I, double Vo)
{
	double R, Vbias, Vdrive;
	
	for (R = 100; R < 150.0e3; R *= pow(10.0, 0.02)) {
		Vbias = bias(m, V, I);
		Vdrive = drive(m, Vo, V, I, R);
		printf("%g %g %g %g %g    ", R, V, I, Vbias, Vdrive);
		distortion(m, Vdrive, V, I, R);
	}
}

/* Tune working point for a fixed plate power and small signal output */
void signaltune(model *m, double P, double Vo)
{
	double R, V, I, Vbias, Vdrive;
	
	for (R = 100; R < 150.0e3; R *= pow(10.0, 0.02)) {
		V = sqrt(P*R);
		I = sqrt(P/R)*1000.0;
		Vbias = bias(m, V, I);
		Vdrive = drive(m, Vo, V, I, 0.0);
		printf("%g %g %g %g %g    ", R, V, I, Vbias, Vdrive);
		distortion(m, Vdrive, V, I, 0.0);
	}
}

/* Tune working point for a fixed plate and output power */
void powertune(model *m, double P, double kpd)
{
	double R, V, I, Vbias, Vdrive;
	
	for (R = 10; R < 15.0e3; R *= pow(10.0, 0.02)) {
		V = sqrt(P*R);
		I = sqrt(P/R)*1000.0;
		Vbias = bias(m, V, I);
		Vdrive = drive(m, sqrt(kpd)*V, V, I, R);
		printf("%g %g %g %g %g    ", R, V, I, Vbias, Vdrive);
		distortion(m, Vdrive, V, I, R);
	}
}



/**********************************************************************/

int main(int argc, char *argv[])
{
	char c;
	int n; double **d;
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv2345P:L:f:dC:")) != -1)
	switch (c) {
	/* General options */
		case 'h':				/* Help message */
			usage(); break;
		case 'v':				/* Verbosity */
			verbose++; break;
	
	/* Operation parameters */
		case '2':				/* Diode */
			vtype = 2; break;
		case '3':				/* Triode */
			vtype = 3; break;
		case '4':				/* Tetrode */
			vtype = 4; break;
		case '5':				/* Pentode */
			vtype = 5; break;
		case 'P':				/* Rated power */
			Pa = strtod(optarg, NULL); break;
		case 'L':				/* Loadline */
			{
				char *p = optarg, *q;
				
				V0 = strtod(p, &q);
				if (p == q) usage();
				
				if (*(q++) != ',') usage();
				I0 = strtod(q, &p);
				if (p == q) usage();
				
				if (*(p++) == ',') {
					RL = strtod(p, &q);
					if (p == q) usage();
				}
			}
			break;
	
	/* I/O functions */
		case 'f':				/* Tagged data format */
			format = optarg; break;
		case 'd':				/* Output data only */
			output_only = 1; break;
	
	/* Model fit options */
		case 'C':				/* Data cuts */
			while (*optarg) switch (*(optarg++)) {
				case 'g':
					grid_cut = 1; break;
				case 'p':
					power_cut = 1; break;
				case 'l':
					loadline_cut = 1; break;
				default:
					usage();
			}
			break;
	
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
	//best_model(d, n);
	
	//fit_curve(d, n, mindex+3);
	//plot_curves(stdout, d, n, mindex+3, 0, 0, 0, 0);
	plot_curves(stdout, d, n, best_model(d, n), 0, 0, 0, 0);
	
	//printf("%g\n", bias(mindex+9, 125.0, 100.0));
	//printf("%g\n", output(mindex+9, -0.0, 150.0, 30.0, 1.714));
	
	//powertune(mindex+9, 13.0, 0.5);
	//signaltune(mindex+9, 13.0, 1.0);
	//loadtune(mindex+9, 125.0, 100.0, 1.0);
	
	return 0;
}
