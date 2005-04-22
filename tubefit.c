/* $Id: tubefit.c,v 1.16 2005/04/22 06:50:03 afrolov Exp $ */

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

#define WAVE_PTS 64



/******************* Options and defaults *****************************/

/* Usage */
char *usage_msg[] = {
	"Vacuum tube model builder, Version " VERSION,
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
	"  -P Pa 	rated anode dissipation power",
	"  -L Vp,Ip[,R]	loadline: working point and load resistance",
	"  -[I|O] V	specify input/output AC signal amplitude",
	"",
	" Model fit options:",
	"  -C flags	apply cuts to tube data to fit specific requirements",
	"		(g = negative grid; p = rated power; l = loadline)",
	"  -M model	use user-supplied model instead of finding best fitting one",
	"",
	" I/O functions:",
	"  -f format	specify format for tagged input data",
	"  -d		dump data in plain format for later use",
	"  -m		[GUI] list available models and their fits",
	"  -p circuit	[GUI] produce plate curves plot (SE or composite)",
	"		(SE = single ended; CF = cathode follower; PP = push-pull)",
	"  -w		[GUI] do waveform analysis",
	NULL
};


/* Options */
static int verbose = 0;

static int                vtype = 3;
static double                Pa = 0.0;
static double                V0 = 0.0;
static double                I0 = 0.0;
static double                RL = 0.0;
static double               Vin = 0.0;
static double              Vout = 0.0;

static char             *format = NULL;
static int          output_only = 0;
static int          list_models = 0;
static int          plot_curves = 0;
static int             waveform = 0;

static int             grid_cut = 0;
static int            power_cut = 0;
static int         loadline_cut = 0;

static char         *user_model = NULL;



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


/* create an nxn identity matrix */
double **identm(unsigned long n)
{
	register unsigned long i, j;
	register double **m = matrix(n, n);
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) m[i][j] = 0.0; m[i][i] = 1.0;
	}
	
	return m;
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


/******************* Multidimensional minimization ********************/

/* directional minimization along a vector xi in n dimensions */
double vmin(double p[], double xi[], int n, double (*f)(double []), double eps)
{
	double a, b, c, u, v, w, x, fa, fb, fc, fu, fv, fw, fx;
	double q, r, s, t, tol, e = 0.0, d = 0.0;
	
	int i, maxiter = 84; // maximal number of iterations
	
	#define GOLD  1.61803398874989484820458683436563811772030918
	#define CGOLD 0.38196601125010515179541316563436188227969082
	
	#define EVAL(X,F) { double t[n]; for (i = 0; i < n; i++) t[i] = p[i] + (X)*xi[i]; (F) = (*f)(t); }
	#define SWAP(A,B) { double T = (A); (A) = (B); (B) = T; } 
	#define SIGN(A,B) ((B) >= 0.0 ? fabs(A) : -fabs(A))
	
	/* initial bracketing of a minimum */
	a = 0.0; b = 1.0; EVAL(a,fa); EVAL(b,fb);
	if (fb > fa) { SWAP(a,b); SWAP(fa,fb); }
	c = b + GOLD*(b-a); EVAL(c,fc);
	
	while (fb > fc) {
		a = b; b = c; fa = fb; fb = fc;
		c = b + GOLD*(b-a); EVAL(c,fc);
	}
	
	/* Brent's minimization */
	x = w = v = b; fx = fw = fv = fb; if (a > c) SWAP(a,c);
	
	while (maxiter--) {
		b = (a+c)/2.0; tol = eps * sqrt(1.0 + x*x);
		if (fabs(x-b) <= (2.0*tol - (c-a)/2.0)) break;
		
		if (fabs(e) > tol) {
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			s = (x-v)*q-(x-w)*r;
			
			q = 2.0*(q-r); if (q > 0.0) s = -s; else q = -q;
			
			t = e; e = d;
			
			if (fabs(s) >= fabs(0.5*q*t) || s <= q*(a-x) || s >= q*(c-x)) {
				e = (x >= b ? a-x : c-x); d = CGOLD * e;
			} else {
				d = s/q; u = x+d; if (u-a < 2.0*tol || c-u < 2.0*tol) d = SIGN(tol,b-x);
			}
		} else { e = (x >= b ? a-x : c-x); d = CGOLD * e; }
		
		u = (fabs(d) >= tol ? x+d : x+SIGN(tol,d)); EVAL(u,fu);
				
		if (fu <= fx) {
			if (u >= x) a = x; else c = x;
			v = w; w = x; x = u;
			fv = fw; fw = fx; fx = fu;
		} else {
			if (u < x) a = u; else c = u;
			if (fu <= fw || w == x) { v = w; w = u; fv = fw; fw = fu; }
			else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
		}
	}
	
	/* update direction vectors */
	for (i = 0; i < n; i++) { xi[i] *= x; p[i] = p[i] + xi[i]; }
	
	return fx;
}

/* n-dimensional minimization (using Powell's method) */
double mmin(double p[], double **e, int n, double (*f)(double []), double eps)
{
	int i, j, s; double fp, fo, fx, sd, t; double po[n], px[n], xi[n];
	
	int maxiter = 1000; // maximal number of iterations
	
	fp = (*f)(p); for (j = 0; j < n; j++) po[j] = p[j];
	
	while (maxiter--) {
		fo = fp; s = 0; sd = 0.0;
		
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) xi[j] = e[j][i];
			
			t = fp; fp = vmin(p, xi, n, f, eps);
			if (fabs(fp-t) > sd) { sd = fabs(fp-t); s = i; }
		}
		
		if (2.0*fabs(fo-fp) <= eps*eps*(fabs(fo)+fabs(fp))) return fp;
		
		for (j = 0; j < n; j++) {
			px[j] = 2.0*p[j]-po[j];
			xi[j] = p[j]-po[j];
			po[j] = p[j];
		}
		
		fx = (*f)(px);
		
		if (fx < fo && 2.0*(fo-2.0*fp+fx)*(fo-fp-sd)*(fo-fp-sd) < sd*(fo-fx)*(fo-fx)) {
			fp = vmin(p, xi, n, f, eps);
			
			for (j = 0; j < n; j++) {
				e[j][s]=e[j][n-1]; e[j][n-1]=xi[j];
			}
		}
	}
	
	return fp;
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
static double init_koren5[] = {PRIOR(0), PRIOR(1), PRIOR(2), 0.0, PRIOR(3)};
static double triode_koren5(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kp = p[1], mu = p[2], Kv = p[3], gamma = p[4];
	
	U = Vp * log(1.0 + exp(Kp + Kp*mu*Vg/sqrt(1000.0*Kv + Vp*Vp)))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}

/* Vacuum triode; modified Koren model (6 parameters) */
static double init_koren6[] = {PRMREF(2,0), 0.0, PRMREF(2,1), PRMREF(2,2), 0.0, PRMREF(2,3)};
static double triode_koren6(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kc = p[1], Kp = p[2], mu = p[3], nu = p[4], gamma = p[5];
	
	U = Vp * log(1.0 + Kc + exp(Kp + Kp*(mu+nu*Vg/1000.0)*Vg/Vp))/Kp;
	I = K * pow(uramp(U), gamma);
	
	return I;
}

/* Vacuum triode; modified Koren model (8 parameters) */
static double init_koren8[] = {PRIOR(0), PRIOR(1), PRIOR(2), PRIOR(3), PRIOR(4), 0.0, 0.0, PRIOR(5)};
static double triode_koren8(double p[], double V[])
{
	double I, U;
	double Vp = V[0], Vg = V[1];
	double K = p[0], Kc = p[1], Kp = p[2], mu = p[3], nu = p[4], Kv = p[5], Vc = p[6], gamma = p[7];
	
	U = Vp * log(1.0 + Kc + exp(Kp + Kp*(mu+nu*Vg/1000.0)*Vg/sqrt(Kv*Kv+(Vp-Vc)*(Vp-Vc))))/Kp;
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
	{3, "Modified Koren model (8 parameters)", "koren8", 8, triode_koren8, init_koren8},
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

/* Initialize global fit context */
static void init_fit_context(model *m, double **data, int n)
{
	int i;
	double *w = vector(n);
	
	/* point weights */
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
	
	_pts_ = n;
	_data_ = data;
	_weight_ = w;
	_model_ = m;
}

/* Fit parametric curve model to the data using simulated annealing */
double *fit_curve(model *m, double **data, int n)
{
	int i, j, k, D = m->params;
	double *p = vector(D+1), **e = identm(D);
	
	init_fit_context(m, data, n);
	
	/* initialize parameter vector */
	for (i = 0; i < D; i++) {
		if (!m->p) { p[i] = 0.0; continue; }
		if (((int)(m->p[i]) >> 8) != PRIOR_MAGIC) { p[i] = m->p[i]; continue; }
		
		j = ((int)(m->p[i]) & 0xF0) >> 4;
		k = (int)(m->p[i]) & 0x0F;
		p[i] = (m-j)->p[k];
	}
	
	/* optimize */
	p[D] = mmin(p, e, D, chi2, 1.0e-6); m->p = p;
	
	free_vector(_weight_);
	free_matrix(e);
	
	return p;
}


/* Try all appropriate models and return the best one */
model *best_model(FILE *fp, double **data, int n)
{
	double *p, min = HUGE;
	int i, j, k, best = 0, invoke = 0;
	
	for (i = 0, k = 0; i < sizeof(mindex)/sizeof(model); i++) {
		model *m = &(mindex[i]);
		int D = m->params;
		
		if (m->vtype != vtype) continue;
		else p = fit_curve(m, data, n);
		if (p[D] == HUGE) continue;
		
		if (p[D] < min) { min = p[D]; best = i; invoke = k; }
		
		if (verbose) {
			fprintf(stderr, "* %s: mean fit error %g mA\n", m->name, sqrt(p[D]));
			
			fprintf(stderr, "*   %s(", m->macro);
			for (j = 0; j < D; j++)
				fprintf(stderr, "%.10g%s", p[j], ((j < D-1) ? "," : ""));
			fprintf(stderr, ")\n"); fflush(stderr);
		}
		
		if (fp) {
			fprintf(fp, "add command -label \"%s \\[mean fit error %g mA\\]\"", m->name, sqrt(p[D]));
			
			fprintf(fp, " -command { set macro \"%s\"; set mparams \"", m->macro);
			for (j = 0; j < D; j++)
				fprintf(fp, "%.10g%s", p[j], ((j < D-1) ? "," : ""));
			fprintf(fp, "\" }\n");
		}
		
		k++;
	}
	
	if (fp) fprintf(fp, "invoke %i\n", invoke);
	
	return &(mindex[best]);
}

/* Read model specification from a string */
model *read_model(const char *s, double **data, int n)
{
	int i, j, idx;
	char *macro = xstrdup(s), *p = strchr(macro, '('), *q;
	
	if (!p) error("Invalid model specification '%s'", s); *(p++) = 0;
	
	for (i = 0, idx = -1; i < sizeof(mindex)/sizeof(model); i++) {
		model *m = &(mindex[i]);
		int D = m->params;
		
		if (strcmp(macro, m->macro)) continue; else idx = i;
		
		m->p = vector(D+1);
		for (j = 0; j < D; j++) {
			m->p[j] = strtod(p, &q);
			if ((p == q) || (*q && (*q != ')') && (*(q++) != ',')))
				error("Error parsing model parameter '%s'", p);
			p = q;
		}
		
		init_fit_context(m, data, n);
		m->p[D] = chi2(m->p);
		free_vector(_weight_);
		
		break;
	}
	
	if (idx == -1) error("Model '%s' not found", macro);
	
	return &(mindex[idx]);
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



/****************** Numerical root finder *****************************/

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



/****************** Bias and loadlines ********************************/

/* Working point and other stuff is passed in global variables */
static double _Vp_, _Vg_, _I_, _Vb_, _Vbias_, _G_, _Vout_;


/* Find grid bias for a given working point */
static double _bias_(double Vg)
{
	double V[] = {_Vp_, Vg, _Vp_}, I = (*(_model_->curve))(_model_->p, V);
	
	return I - _I_;
}

double bias(model *m, double Vp, double Ip)
{
	_model_ = m; _Vp_ = Vp; _I_ = Ip;
	
	return (_bias_(0.0) > 0.0) ?
		zbrent(_bias_, -Vp, 0.0, 1.0e-12):
		zbrent(_bias_, 0.0, +Vp, 1.0e-12);
}


/****************** Single-ended loadline *****************************/

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
	double Vmin, Vmax;
	
	_Vg_ = _Vbias_+Vi; Vmin = zbrent(_loadline_, 0, 3000, 1.0e-12);
	_Vg_ = _Vbias_-Vi; Vmax = zbrent(_loadline_, 0, 3000, 1.0e-12);
	
	return fabs(Vmax-Vmin) - 2.0*_Vout_;
}

double drive(model *m, double Vout, double Vb, double I0, double R)
{
	_model_ = m; _Vbias_ = bias(m, Vb, I0); _Vout_ = Vout;
	_I_ = I0; _Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return zbrent(_swing_, 0, Vout, 1.0e-12);
}


/****************** Cathode follower loadline *************************/

/* Find tube output on a given loadline */
static double _cf_loadline_(double Vp)
{
	double V[] = {Vp, _Vg_ + (Vp-_Vb_), Vp}, I = (*(_model_->curve))(_model_->p, V);
	
	return (I-_I_) + _G_*(Vp-_Vb_);
}

double cf_output(model *m, double Vg, double Vb, double I0, double R)
{
	_model_ = m; _Vg_ = Vg; _I_ = I0;
	_Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return zbrent(_cf_loadline_, 0, 3000, 1.0e-12);
}


/* Find grid drive for a given output voltage swing */
static double _cf_swing_(double Vi)
{
	double Vmin, Vmax;
	
	_Vg_ = _Vbias_+Vi; Vmin = zbrent(_cf_loadline_, 0, 3000, 1.0e-12);
	_Vg_ = _Vbias_-Vi; Vmax = zbrent(_cf_loadline_, 0, 3000, 1.0e-12);
	
	return fabs(Vmax-Vmin) - 2.0*_Vout_;
}

double cf_drive(model *m, double Vout, double Vb, double I0, double R)
{
	_model_ = m; _Vbias_ = bias(m, Vb, I0); _Vout_ = Vout;
	_I_ = I0; _Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return zbrent(_cf_swing_, 0, 2.0*Vout, 1.0e-12);
}


/****************** Push-pull loadline ********************************/

/* Find tube output on a given push-pull loadline */
static double _pp_loadline_(double Vp)
{
	double V1[] = {_Vb_+Vp, _Vbias_+_Vg_, _Vb_+Vp},
	       V2[] = {_Vb_-Vp, _Vbias_-_Vg_, _Vb_-Vp},
	       Ip1 = (_Vb_+Vp > 0) ? (*(_model_->curve))(_model_->p, V1) : 0.0,
	       Ip2 = (_Vb_-Vp > 0) ? (*(_model_->curve))(_model_->p, V2) : 0.0;
	
	return (Ip2-Ip1) - 2.0*Vp*_G_;
}

double pp_output(model *m, double Vg, double Vbias, double Vb, double R)
{
	_model_ = m; _Vg_ = Vg; _Vbias_ = Vbias;
	_Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return zbrent(_pp_loadline_, -Vb, Vb, 1.0e-12);
}


/* Find grid drive for a given push-pull output voltage swing */
static double _pp_loadline2_(double Vg)
{
	double V1[] = {_Vb_+_Vp_, _Vbias_+Vg, _Vb_+_Vp_},
	       V2[] = {_Vb_-_Vp_, _Vbias_-Vg, _Vb_-_Vp_},
	       Ip1 = (_Vb_+_Vp_ > 0) ? (*(_model_->curve))(_model_->p, V1) : 0.0,
	       Ip2 = (_Vb_-_Vp_ > 0) ? (*(_model_->curve))(_model_->p, V2) : 0.0;
	
	return (Ip2-Ip1) - 2.0*_Vp_*_G_;
}

double pp_drive(model *m, double Vp, double Vbias, double Vb, double R)
{
	_model_ = m; _Vp_ = Vp; _Vbias_ = Vbias;
	_Vb_ = Vb; _G_ = (R != 0.0) ? 1000.0/R : 0.0;
	
	return (Vp > 0.0) ?
		zbrent(_pp_loadline2_, -Vp, 0.0, 1.0e-12):
		zbrent(_pp_loadline2_, 0.0, +Vp, 1.0e-12);
}



/****************** Distortion analysis *******************************/

#ifdef FFTW

#include <rfftw.h>

/* Calculate power spectrum and harmonic distortion */
static void distortion(double *f, double *S, double *D, int N)
{
	int k; double F[N], N2 = N*N;
	rfftw_plan p = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
	
	rfftw_one(p, f, F);
	
	/* Power spectrum */
	S[0] = F[0]*F[0]/N2;				/* DC component */
	for (k = 1; k < (N+1)/2; k++)			/* k-th harmonic */
		S[k] = 4.0*(F[k]*F[k] + F[N-k]*F[N-k])/N2;
	if (N % 2 == 0) S[N/2] = 4.0*F[N/2]*F[N/2]/N2;	/* Nyquist freq */
	
	/* Harmonic distortion */
	D[0] = D[1] = 0.0;
	for (k = N/2; k > 1; k--) {
		D[k] = sqrt(S[k]/S[1]);
		if (k > 2) D[k % 2] += S[k];
	}
	D[0] = sqrt(D[0]/S[1]);
	D[1] = sqrt(D[1]/S[1]);
	
	rfftw_destroy_plan(p);
}
#endif /* FFTW */



/****************** Curves and waveforms ******************************/

/* Plot SE plate curves with Tk toolkit */
void se_plate_curves(FILE *fp, model *m, double **data, int n, double Vmax, double Imax, double Vgm, double Vgs)
{
	int i; double Vp, Vg, Ip, t;
	
	/* Canvas layout: 720x576       */
	/* Margins: l=40,r=40,t=40,b=20 */
	#define X(V) (40.0 + 640.0*(V)/Vmax)
	#define Y(I) (556.0 - 516.0*(I)/Imax)
	
	
	/* Auto ranges */
	if ((Vmax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[0][i] > Vmax) Vmax = data[0][i];
		t = 5.0 * pow(10.0, floor(log10(Vmax))-1); Vmax = t * ceil(Vmax/t);
	}
	
	if ((Imax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[3][i] > Imax) Imax = data[3][i];
		t = 5.0 * pow(10.0, floor(log10(Imax))-1); Imax = t * ceil(Imax/t);
	}
	
	if ((Vgm == 0.0) && data) {
		for (i = 0; i < n; i++) {
			if (fabs(data[1][i]) > fabs(Vgm)) Vgm = data[1][i];
			if ((Vgs == 0.0) && (data[1][i] != 0.0)) Vgs = data[1][i];
		}
	}
	
	if (Vgs == 0.0) { Vgs = -1.0; }
	
	
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
					X(Vp)-2, Y(Ip)+2, X(Vp)+2, Y(Ip)-2);
	}
	
	/* Plate curves */
	if (m) for (Vg = 0.0; fabs(Vg) <= fabs(Vgm); Vg += Vgs) {
		fprintf(fp, "line ");
		
		for (Vp = 0.0, Ip = 0.0; (Vp <= Vmax) && (Ip <= Imax); Vp += Vmax/200) {
			double V[] = {Vp, Vg, Vp};
			
			Ip = (*(m->curve))(m->p, V);
			fprintf(fp, "%g %g ", X(Vp), Y(Ip));
		}
		
		fprintf(fp, "%g %g -smooth 1 -fill black -width 2\n", X(Vp), Y(Ip));
		fprintf(fp, "text %g %g -anchor nw -text %g -fill black\n",
				X(Vp)+2, Y((Ip > Imax) ? Imax : Ip)+2, Vg);
	}
	
	/* Power */
	if (Pa > 0.0) {
		fprintf(fp, "line ");
		
		for (Vp = 1000.0*Pa/Imax; Vp <= Vmax; Vp += Vmax/200)
			fprintf(fp, "%g %g ", X(Vp), Y(1000.0*Pa/Vp));
		
		fprintf(fp, "-smooth 1 -fill red -width 3\n");
		fprintf(fp, "text %g %g -anchor s -text \"%g W\" -font {helvetica 24 bold} -fill red\n",
				X(0.93*Vmax), Y(1075.0*Pa/Vmax)-5, Pa);
	}
	
	/* Loadline */
	if (V0 > 0.0) {
		double Vbias = bias(m, V0, I0);
		double G = (RL != 0.0) ? 1000.0/RL : 0.0;
		double V1 = 0.0, I1 = I0+G*V0, V2 = Vmax, I2 = I0+G*(V0-Vmax);
		
		if ((I1 > Imax) && (G != 0.0)) { I1 = Imax; V1 = V0 - (Imax-I0)/G; }
		if ((I2 < 0.0) && (G != 0.0)) { I2 = 0.0; V2 = V0 + I0/G; }
		
		fprintf(fp, "line %g %g %g %g -fill blue -width 3\n",
				X(V1), Y(I1), X(V2), Y(I2));
		fprintf(fp, "text %g %g -anchor sw -text \"%g R\" -font {helvetica 14 bold} -fill blue\n",
				X(V2)+5, Y(I2)-5, RL);
		fprintf(fp, "oval %g %g %g %g -outline blue -width 2\n",
				X(V0)-5, Y(I0)+5, X(V0)+5, Y(I0)-5);
		fprintf(fp, "text %g %g -anchor ne -text %.3g -fill blue\n",
				X(V0)-5, Y(I0)+5, Vbias);
		
		if (Vout > 0.0) { Vin = drive(m, Vout, V0, I0, RL); }
		
		if (Vin > 0.0) {
			V1 = output(m, Vbias+Vin, V0, I0, RL); I1 = I0+G*(V0-V1);
			V2 = output(m, Vbias-Vin, V0, I0, RL); I2 = I0+G*(V0-V2);
			
			fprintf(fp, "line %g %g %g %g -fill cyan -width 3 -arrow both\n",
					X(V1), Y(I1), X(V2), Y(I2));
		}
		
		if (Vin > 0.0 && waveform) {
			int N = WAVE_PTS;
			double *V = vector(N), *S = vector(N/2+1), *D = vector(N/2+1);
			
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V0), Y(0), X(V0), Y(Imax));
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V1), Y(0), X(V1), Y(Imax));
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V2), Y(0), X(V2), Y(Imax));
			
			fprintf(fp, "line ");
			for (i = 0; i < N; i++) {
				V[i] = output(m, Vbias + Vin*sin(2.0*M_PI*i/N), V0, I0, RL);
				fprintf(fp, "%g %g ", X(V[i]), Y(I0 + 0.2*Imax*(i-N/2)/N));
			}
			fprintf(fp, "-smooth 1 -fill green -width 2\n");
			
			#ifdef FFTW
			distortion(V, S, D, N);
			
			fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
					X(0.5*Vmax),Y(0.9*Imax), X(0.97*Vmax),Y(0.9*Imax), X(0.97*Vmax),Y(0.5*Imax), X(0.5*Vmax),Y(0.5*Imax));
			fprintf(fp, "text %g %g -anchor n -text \""
					"THD: %.3g%% (%.3g%% second, %.3g%% odd, %.3g%% even)"
					"\" -fill black -font {helvetica 8}\n",
					X(0.735*Vmax), Y(0.9*Imax)+3,
					100.0*sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]), 100.0*D[2], 100.0*D[1], 100.0*D[0]);
			
			fprintf(fp, "text %g %g -anchor ne -text \"", X(0.95*Vmax), Y(0.85*Imax));
			for (i = 0; i < 10; i++)
				fprintf(fp, "%2i %11.6fV %11.6f%% %8.2fdb\\n", i,
					sqrt(S[i]), 100.0*sqrt(S[i]/S[1]), 10.0*log10(S[i]/S[1]));
			fprintf(fp, "\" -fill black -font {courier 6}\n");
			
			for (i = 2; i <= 16; i++) {
				double hdb = 20*log10(D[i]), hf = -120.0;
				double x = 0.5 + 0.03*(i-1), y = 0.53 - ((hdb > hf) ? 0.35*(hdb-hf)/hf : 0.0);
				
				fprintf(fp, "line %g %g %g %g -fill yellow -width 14\n",
						X(x*Vmax), Y(0.53*Imax), X(x*Vmax), Y(y*Imax));
				fprintf(fp, "text %g %g -anchor n -justify center -text \"%i\" -fill black -font {helvetica 7}\n",
						X(x*Vmax), Y(0.53*Imax)+2, i);
				fprintf(fp, "text %g %g -anchor s -justify center -text \"%.3g\" -fill black -font {helvetica 6}\n",
						X(x*Vmax), Y(y*Imax)-2, -hdb);
			}
			#endif /* FFTW */
			
			free_vector(V); free_vector(S); free_vector(D);
		}
		
		/* Loadline summary */
		fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
				X(0.5*Vmax),Y(0.96*Imax), X(0.97*Vmax),Y(0.96*Imax), X(0.97*Vmax),Y(0.9*Imax), X(0.5*Vmax),Y(0.9*Imax));
		fprintf(fp, "text %g %g -anchor w -text \""
				"Working point: %.3g V, %.3g mA\\n"
				"Bias: %.3g V (Rk = %.4g)"
				"\" -fill black -font {helvetica 8}\n",
				X(0.5*Vmax)+5, Y(0.93*Imax),
				V0, I0, Vbias, fabs(Vbias/I0*1000.0));
		
		/* Signal summary */
		if (Vin > 0.0) {
			fprintf(fp, "text %g %g -anchor e -text \""
					"Input: %.3g V\\n"
					"Output: %.3g V (%.4g W)"
					"\" -fill black -font {helvetica 8}\n",
					X(0.97*Vmax)-5, Y(0.93*Imax),
					Vin, fabs(V2-V1)/2.0, (RL != 0.0) ? (V2-V1)*(V2-V1)/8.0/RL : 0.0);
		}
	}
	
	/* Annotation */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
			X(0),Y(Imax), X(Vmax),Y(Imax), X(Vmax),10.0, X(0),10.0);
	fprintf(fp, "text %g %g -anchor w -text \"%s: mean fit error %g mA\\n%s(",
			X(0)+3, (Y(Imax)+12.0)/2.0, m->name, sqrt(m->p[m->params]), m->macro);
	for (i = 0; i < m->params; i++)
		fprintf(fp, "%.10g%s", m->p[i], ((i < m->params-1) ? "," : ""));
	fprintf(fp, ")\" -fill black\n");
	
	
	#undef X
	#undef Y
}


/* Plot cathode follower plate curves with Tk toolkit */
void cf_plate_curves(FILE *fp, model *m, double **data, int n, double Vmax, double Imax, double Vgm, double Vgs)
{
	int i; double Vp, Vg, Ip, t;
	
	/* Canvas layout: 720x576       */
	/* Margins: l=40,r=40,t=40,b=20 */
	#define X(V) (40.0 + 640.0*(V)/Vmax)
	#define Y(I) (556.0 - 516.0*(I)/Imax)
	
	
	/* Auto ranges */
	if ((Vmax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[0][i] > Vmax) Vmax = data[0][i];
		t = 5.0 * pow(10.0, floor(log10(Vmax))-1); Vmax = t * ceil(Vmax/t);
	}
	
	if ((Imax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[3][i] > Imax) Imax = data[3][i];
		t = 5.0 * pow(10.0, floor(log10(Imax))-1); Imax = t * ceil(Imax/t);
	}
	
	if ((Vgm == 0.0) && data) {
		for (i = 0; i < n; i++) {
			if (fabs(data[1][i]) > fabs(Vgm)) Vgm = data[1][i];
			if ((Vgs == 0.0) && (data[1][i] != 0.0)) Vgs = data[1][i];
		}
	}
	
	if (Vgs != 0.0) { Vgs *= ceil(fabs(Vmax/10.0/Vgs)); } else { Vgs = Vmax/10.0; }
	
	
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
	
	/* Plate curves */
	if (m) for (Vg = V0; Vg >= V0-Vmax+Vgm; Vg -= fabs(Vgs)) {
		fprintf(fp, "line ");
		
		for (Vp = 0.0, Ip = 0.0; (Vp <= Vmax) && (Ip <= Imax); Vp += Vmax/200) {
			double V[] = {Vp, Vg + (Vp-V0), Vp};
			
			Ip = (*(m->curve))(m->p, V);
			fprintf(fp, "%g %g ", X(Vp), Y(Ip));
		}
		
		fprintf(fp, "%g %g -smooth 1 -fill black -width 2\n", X(Vp), Y(Ip));
		fprintf(fp, "text %g %g -anchor nw -text %g -fill black\n",
				X(Vp)+2, Y((Ip > Imax) ? Imax : Ip)+2, Vg);
	}
	
	/* Class A2 limit */
	if (m) {
		double Va = 0.0;
		
		fprintf(fp, "line ");
		
		for (Vp = 0.0, Ip = 0.0; (Vp <= Vmax) && (Ip <= Imax); Vp += Vmax/200) {
			double V[] = {Vp, 0.0, Vp};
			
			Ip = (*(m->curve))(m->p, V);
			fprintf(fp, "%g %g ", X(Vp), Y(Ip));
			if (!Va && Ip > Imax/2.0) Va = Vp;
		}
		
		fprintf(fp, "%g %g -smooth 1 -fill yellow -width 3\n", X(Vp), Y(Ip));
		fprintf(fp, "text %g %g -anchor se -text \"A2\" -font {helvetica 36 bold} -fill yellow\n",
				X(Va)-10, Y(Imax/2.0));
	}
	
	/* Power */
	if (Pa > 0.0) {
		fprintf(fp, "line ");
		
		for (Vp = 1000.0*Pa/Imax; Vp <= Vmax; Vp += Vmax/200)
			fprintf(fp, "%g %g ", X(Vp), Y(1000.0*Pa/Vp));
		
		fprintf(fp, "-smooth 1 -fill red -width 3\n");
		fprintf(fp, "text %g %g -anchor s -text \"%g W\" -font {helvetica 24 bold} -fill red\n",
				X(0.93*Vmax), Y(1075.0*Pa/Vmax)-5, Pa);
	}
	
	/* Loadline */
	if (V0 > 0.0) {
		double Vbias = bias(m, V0, I0);
		double G = (RL != 0.0) ? 1000.0/RL : 0.0;
		double V1 = 0.0, I1 = I0+G*V0, V2 = Vmax, I2 = I0+G*(V0-Vmax);
		
		if ((I1 > Imax) && (G != 0.0)) { I1 = Imax; V1 = V0 - (Imax-I0)/G; }
		if ((I2 < 0.0) && (G != 0.0)) { I2 = 0.0; V2 = V0 + I0/G; }
		
		fprintf(fp, "line %g %g %g %g -fill blue -width 3\n",
				X(V1), Y(I1), X(V2), Y(I2));
		fprintf(fp, "text %g %g -anchor sw -text \"%g R\" -font {helvetica 14 bold} -fill blue\n",
				X(V2)+5, Y(I2)-5, RL);
		fprintf(fp, "oval %g %g %g %g -outline blue -width 2\n",
				X(V0)-5, Y(I0)+5, X(V0)+5, Y(I0)-5);
		fprintf(fp, "text %g %g -anchor ne -text %.3g -fill blue\n",
				X(V0)-5, Y(I0)+5, Vbias);
		
		if (Vout > 0.0) { Vin = cf_drive(m, Vout, V0, I0, RL); }
		
		if (Vin > 0.0) {
			V1 = cf_output(m, Vbias+Vin, V0, I0, RL); I1 = I0+G*(V0-V1);
			V2 = cf_output(m, Vbias-Vin, V0, I0, RL); I2 = I0+G*(V0-V2);
			
			fprintf(fp, "line %g %g %g %g -fill cyan -width 3 -arrow both\n",
					X(V1), Y(I1), X(V2), Y(I2));
		}
		
		if (Vin > 0.0 && waveform) {
			int N = WAVE_PTS;
			double *V = vector(N), *S = vector(N/2+1), *D = vector(N/2+1);
			
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V0), Y(0), X(V0), Y(Imax));
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V1), Y(0), X(V1), Y(Imax));
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V2), Y(0), X(V2), Y(Imax));
			
			fprintf(fp, "line ");
			for (i = 0; i < N; i++) {
				V[i] = cf_output(m, Vbias + Vin*sin(2.0*M_PI*i/N), V0, I0, RL);
				fprintf(fp, "%g %g ", X(V[i]), Y(I0 + 0.2*Imax*(i-N/2)/N));
			}
			fprintf(fp, "-smooth 1 -fill green -width 2\n");
			
			#ifdef FFTW
			distortion(V, S, D, N);
			
			fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
					X(0.5*Vmax),Y(0.9*Imax), X(0.97*Vmax),Y(0.9*Imax), X(0.97*Vmax),Y(0.5*Imax), X(0.5*Vmax),Y(0.5*Imax));
			fprintf(fp, "text %g %g -anchor n -text \""
					"THD: %.3g%% (%.3g%% second, %.3g%% odd, %.3g%% even)"
					"\" -fill black -font {helvetica 8}\n",
					X(0.735*Vmax), Y(0.9*Imax)+3,
					100.0*sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]), 100.0*D[2], 100.0*D[1], 100.0*D[0]);
			
			fprintf(fp, "text %g %g -anchor ne -text \"", X(0.95*Vmax), Y(0.85*Imax));
			for (i = 0; i < 10; i++)
				fprintf(fp, "%2i %11.6fV %11.6f%% %8.2fdb\\n", i,
					sqrt(S[i]), 100.0*sqrt(S[i]/S[1]), 10.0*log10(S[i]/S[1]));
			fprintf(fp, "\" -fill black -font {courier 6}\n");
			
			for (i = 2; i <= 16; i++) {
				double hdb = 20*log10(D[i]), hf = -120.0;
				double x = 0.5 + 0.03*(i-1), y = 0.53 - ((hdb > hf) ? 0.35*(hdb-hf)/hf : 0.0);
				
				fprintf(fp, "line %g %g %g %g -fill yellow -width 14\n",
						X(x*Vmax), Y(0.53*Imax), X(x*Vmax), Y(y*Imax));
				fprintf(fp, "text %g %g -anchor n -justify center -text \"%i\" -fill black -font {helvetica 7}\n",
						X(x*Vmax), Y(0.53*Imax)+2, i);
				fprintf(fp, "text %g %g -anchor s -justify center -text \"%.3g\" -fill black -font {helvetica 6}\n",
						X(x*Vmax), Y(y*Imax)-2, -hdb);
			}
			#endif /* FFTW */
			
			free_vector(V); free_vector(S); free_vector(D);
		}
		
		/* Loadline summary */
		fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
				X(0.5*Vmax),Y(0.96*Imax), X(0.97*Vmax),Y(0.96*Imax), X(0.97*Vmax),Y(0.9*Imax), X(0.5*Vmax),Y(0.9*Imax));
		fprintf(fp, "text %g %g -anchor w -text \""
				"Working point: %.3g V, %.3g mA\\n"
				"Bias: %.3g V (Rk = %.4g)"
				"\" -fill black -font {helvetica 8}\n",
				X(0.5*Vmax)+5, Y(0.93*Imax),
				V0, I0, Vbias, fabs(Vbias/I0*1000.0));
		
		/* Signal summary */
		if (Vin > 0.0) {
			fprintf(fp, "text %g %g -anchor e -text \""
					"Input: %.3g V\\n"
					"Output: %.3g V (%.4g W)"
					"\" -fill black -font {helvetica 8}\n",
					X(0.97*Vmax)-5, Y(0.93*Imax),
					Vin, fabs(V2-V1)/2.0, (RL != 0.0) ? (V2-V1)*(V2-V1)/8.0/RL : 0.0);
		}
	}
	
	/* Annotation */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
			X(0),Y(Imax), X(Vmax),Y(Imax), X(Vmax),10.0, X(0),10.0);
	fprintf(fp, "text %g %g -anchor w -text \"%s: mean fit error %g mA\\n%s(",
			X(0)+3, (Y(Imax)+12.0)/2.0, m->name, sqrt(m->p[m->params]), m->macro);
	for (i = 0; i < m->params; i++)
		fprintf(fp, "%.10g%s", m->p[i], ((i < m->params-1) ? "," : ""));
	fprintf(fp, ")\" -fill black\n");
	
	
	#undef X
	#undef Y
}


/* Plot composite PP plate curves with Tk toolkit */
void pp_plate_curves(FILE *fp, model *m, double **data, int n, double Vmax, double Imax, double Vgm, double Vgs)
{
	int i; double Vp, Vg, Ip, t;
	double Vbias = bias(m, V0, I0);
	
	/* Canvas layout: 720x576       */
	/* Margins: l=40,r=40,t=40,b=20 */
	#define X(V) (360.0 + 320.0*(V)/Vmax)
	#define Y(I) (298.0 - 258.0*(I)/Imax)
	
	
	/* Auto ranges */
	if ((Vmax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[0][i] > Vmax) Vmax = data[0][i]; Vmax -= V0;
		t = 5.0 * pow(10.0, floor(log10(Vmax))-1); Vmax = t * ceil(Vmax/t);
	}
	
	if ((Imax == 0.0) && data) {
		for (i = 0; i < n; i++) if (data[3][i] > Imax) Imax = data[3][i]; Imax -= I0;
		t = 5.0 * pow(10.0, floor(log10(Imax))-1); Imax = t * ceil(Imax/t);
	}
	
	if ((Vgm == 0.0) && data) {
		for (i = 0; i < n; i++) {
			if (fabs(data[1][i]) > fabs(Vgm)) Vgm = data[1][i];
			if ((Vgs == 0.0) && (data[1][i] != 0.0)) Vgs = data[1][i];
		}
		Vgm -= Vbias;
	}
	
	if (Vgs == 0.0) { Vgs = -1.0; }
	
	
	/* Axis grid */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
			X(-Vmax),Y(-Imax), X(Vmax),Y(-Imax), X(Vmax),Y(Imax), X(-Vmax),Y(Imax));
	
	for (i = 0; i <= 10; i++) {
		fprintf(fp, "line %g %g %g %g -fill black -width 1 -dash .\n",
				X(Vmax*(i-5)/5), Y(-Imax), X(Vmax*(i-5)/5), Y(Imax));
		fprintf(fp, "line %g %g %g %g -fill black -width 1 -dash .\n",
				X(-Vmax), Y(Imax*(i-5)/5), X(Vmax), Y(Imax*(i-5)/5));
		fprintf(fp, "text %g %g -anchor e -text %g -fill black\n",
				X(-Vmax)-5, Y(Imax*(i-5)/5), Imax*(i-5)/5);
	}
	
	/* Plate curves */
	if (m) for (Vg = - floor(fabs(Vgm/Vgs)) * Vgs; fabs(Vg) <= fabs(Vgm); Vg += Vgs) {
		int class, cc = -1; double A;
		static const char *color[] = {"black", "grey", "orange", "yellow"};
		
		fprintf(fp, "line ");
		
		for (Vp = -Vmax; Vp <= Vmax; Vp += Vmax/400) {
			double V1[] = {V0+Vp, Vbias+Vg, V0+Vp}, Ip1 = (V0+Vp > 0) ? (*(m->curve))(m->p, V1) : 0.0,
			       V2[] = {V0-Vp, Vbias-Vg, V0-Vp}, Ip2 = (V0-Vp > 0) ? (*(m->curve))(m->p, V2) : 0.0;
			
			Ip = Ip1 - Ip2; A = (Ip1 < Ip2) ? Ip1/Ip2 : Ip2/Ip1;
			
			if (Ip < -1.07*Imax) continue;
			if (Ip > +1.07*Imax) break;
			
			class = ((A < 0.03) ? 0x02 : 0x00) |
                                ((fabs(Vg) > -Vbias) ? 0x01 : 0x00);
			
			if ((cc != -1) && (class != cc))
				fprintf(fp, "%g %g -smooth 1 -fill %s -width 2\nline ", X(Vp), Y(Ip), color[cc]);
			
			fprintf(fp, "%g %g ", X(Vp), Y(Ip)); cc = class;
		}
		
		fprintf(fp, "%g %g -smooth 1 -fill %s -width 2\n", X(Vp), Y(Ip), color[cc]);
		fprintf(fp, "text %g %g -anchor nw -text %+g -fill black\n",
				X(Vp)+2, Y((Ip > Imax) ? Imax : Ip)+2, Vg);
	}
	
	/* Power */
	if (Pa > 0.0) {
		fprintf(fp, "line ");
		
		for (Vp = V0-1000.0*Pa/(I0+Imax); Vp >= -Vmax; Vp -= Vmax/200) {
			Ip = I0 - 1000.0*Pa/(V0-Vp);
			fprintf(fp, "%g %g ", X(Vp), Y(Ip));
		}
		
		fprintf(fp, "-smooth 1 -fill red -width 3\n");
		
		fprintf(fp, "line ");
		
		for (Vp = 1000.0*Pa/(I0+Imax)-V0; Vp <= Vmax; Vp += Vmax/200) {
			Ip = 1000.0*Pa/(Vp+V0) - I0;
			fprintf(fp, "%g %g ", X(Vp), Y(Ip));
		}
		
		fprintf(fp, "-smooth 1 -fill red -width 3\n");
		
		fprintf(fp, "text %g %g -anchor s -text \"%g W\" -fill red\n",
				X(0.95*Vmax), Y(Ip)-10, Pa);
		fprintf(fp, "text %g %g -anchor n -text \"%g W\" -fill red\n",
				X(-0.95*Vmax), Y(-Ip)+10, Pa);
	}
	
	/* Loadline */
	if (V0 > 0.0) {
		double G = (RL != 0.0) ? 1000.0/RL : 0.0;
		double V1 = Vmax, I1 = 2.0*V1*G;
		
		if (I1 > Imax) { I1 = Imax; V1 = Imax*RL/2000.0; }
		
		fprintf(fp, "line %g %g %g %g -fill blue -width 3\n",
				X(-V1), Y(I1), X(V1), Y(-I1));
		fprintf(fp, "text %g %g -anchor sw -text \"%g R\" -fill blue\n",
				X(V1)+2, Y(-I1)-5, RL);
		fprintf(fp, "oval %g %g %g %g -outline blue -width 2\n",
				X(0)-5, Y(0)+5, X(0)+5, Y(0)-5);
		fprintf(fp, "text %g %g -anchor ne -text %.3g -fill blue\n",
				X(0)-5, Y(0)+5, Vbias);
		
		
		if (Vout > 0.0) { Vin = fabs(pp_drive(m, Vout, Vbias, V0, RL)); }
		
		if (Vin > 0.0) {
			V1 = pp_output(m, Vin, Vbias, V0, RL); I1 = 2.0*V1*G;
			
			fprintf(fp, "line %g %g %g %g -fill cyan -width 3 -arrow both\n",
					X(-V1), Y(I1), X(V1), Y(-I1));
		}
		
		
		if (Vin > 0.0 && waveform) {
			int N = 2*WAVE_PTS;
			double *V = vector(N), *S = vector(N/2+1), *D = vector(N/2+1);
			
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(0), Y(-Imax), X(0), Y(Imax));
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(-V1), Y(-Imax), X(-V1), Y(Imax));
			fprintf(fp, "line %g %g %g %g -fill green -width 1 -dash -\n",
					X(V1), Y(-Imax), X(V1), Y(Imax));
			
			fprintf(fp, "line ");
			for (i = 0; i < N; i++) {
				V[i] = pp_output(m, Vin*sin(2.0*M_PI*i/N), Vbias, V0, RL);
				fprintf(fp, "%g %g ", X(V[i]), Y(-0.4*Imax + 0.4*Imax*(i-N/2)/N));
			}
			fprintf(fp, "-smooth 1 -fill green -width 2\n");
			
			#ifdef FFTW
			distortion(V, S, D, N);
			
			fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
					X(0),Y(0.8*Imax), X(0.94*Vmax),Y(0.8*Imax), X(0.94*Vmax),Y(0), X(0),Y(0));
			fprintf(fp, "text %g %g -anchor n -text \""
					"THD: %.3g%%"
					"\" -fill black -font {helvetica 8}\n",
					X(0.47*Vmax), Y(0.8*Imax)+3,
					100.0*sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]));
			
			fprintf(fp, "text %g %g -anchor ne -text \"", X(0.9*Vmax), Y(0.7*Imax));
			for (i = 1; i < 10; i += 2)
				fprintf(fp, "%2i %11.6fV %11.6f%% %8.2fdb\\n", i,
					sqrt(S[i]), 100.0*sqrt(S[i]/S[1]), 10.0*log10(S[i]/S[1]));
			fprintf(fp, "\" -fill black -font {courier 6}\n");
			
			for (i = 3; i < 32; i += 2) {
				double hdb = 20*log10(D[i]), hf = -120.0;
				double x = 0.0 + 0.03*(i-1), y = 0.06 - ((hdb > hf) ? 0.70*(hdb-hf)/hf : 0.0);
				
				fprintf(fp, "line %g %g %g %g -fill yellow -width 14\n",
						X(x*Vmax), Y(0.06*Imax), X(x*Vmax), Y(y*Imax));
				fprintf(fp, "text %g %g -anchor n -justify center -text \"%i\" -fill black -font {helvetica 7}\n",
						X(x*Vmax), Y(0.06*Imax)+2, i);
				fprintf(fp, "text %g %g -anchor s -justify center -text \"%.3g\" -fill black -font {helvetica 6}\n",
						X(x*Vmax), Y(y*Imax)-2, -hdb);
			}
			#endif /* FFTW */
			
			free_vector(V); free_vector(S); free_vector(D);
		}
		
		/* Loadline summary */
		fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
				X(0),Y(0.92*Imax), X(0.94*Vmax),Y(0.92*Imax), X(0.94*Vmax),Y(0.8*Imax), X(0),Y(0.8*Imax));
		fprintf(fp, "text %g %g -anchor w -text \""
				"Working point: %.3g V, %.3g mA\\n"
				"Bias: %.3g V (Rk = %.4g)"
				"\" -fill black -font {helvetica 8}\n",
				X(0)+5, Y(0.86*Imax),
				V0, I0, Vbias, fabs(Vbias/I0*1000.0));
		
		/* Signal summary */
		if (Vin > 0.0) {
			fprintf(fp, "text %g %g -anchor e -text \""
					"Input: %.3g V\\n"
					"Output: %.3g V (%.4g W)"
					"\" -fill black -font {helvetica 8}\n",
					X(0.94*Vmax)-5, Y(0.86*Imax),
					Vin, fabs(V1), (RL != 0.0) ? 2.0*V1*V1/RL : 0.0);
		}
	}
	
	/* V axis annotation (cover-up) */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline white -fill white -width 1\n",
			X(-Vmax),Y(-Imax)+1, X(Vmax),Y(-Imax)+1, X(Vmax),Y(-1.1*Imax), X(-Vmax),Y(-1.1*Imax));
	for (i = 0; i <= 10; i++)
		fprintf(fp, "text %g %g -anchor n -text %g -fill black\n",
				X(Vmax*(i-5)/5), Y(-Imax)+5, Vmax*(i-5)/5);
	
	/* Push-Pull legend */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
			X(0),Y(-0.92*Imax), X(-0.94*Vmax),Y(-0.92*Imax), X(-0.94*Vmax),Y(-0.8*Imax), X(0),Y(-0.8*Imax));
	fprintf(fp, "text %g %g -anchor c -text \""
			"Push-Pull Composite"
			"\" -fill black -font {helvetica 18 bold}\n",
			X(-0.47*Vmax), Y(-0.86*Imax),
			V0, I0, Vbias, fabs(Vbias/I0*1000.0));
	
	/* Annotation */
	fprintf(fp, "polygon %g %g %g %g %g %g %g %g -outline black -fill white -width 1\n",
			X(-Vmax),Y(Imax), X(Vmax),Y(Imax), X(Vmax),10.0, X(-Vmax),10.0);
	fprintf(fp, "text %g %g -anchor w -text \"%s: mean fit error %g mA\\n%s(",
			X(-Vmax)+3, (Y(Imax)+12.0)/2.0, m->name, sqrt(m->p[m->params]), m->macro);
	for (i = 0; i < m->params; i++)
		fprintf(fp, "%.10g%s", m->p[i], ((i < m->params-1) ? "," : ""));
	fprintf(fp, ")\" -fill black\n");
	
	
	#undef X
	#undef Y
}


/* Circuit index */
static struct {char *name; void (*curves)(FILE *fp, model *m, double **data, int n, double Vmax, double Imax, double Vgm, double Vgs); } circuit[] = {
	{"NONE", NULL},
	{"SE",   se_plate_curves},
	{"CF",   cf_plate_curves},
	{"PP",   pp_plate_curves},
	NULL
};



/**********************************************************************/

int main(int argc, char *argv[])
{
	char c;
	int n; double **d; model *m;
	
	/* Parse options */
	while ((c = getopt(argc, argv, "hv2345P:L:I:O:f:dmp:wC:M:")) != -1)
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
		case 'I':				/* Drive voltage */
			Vin = strtod(optarg, NULL); break;
		case 'O':				/* Swing voltage */
			Vout = strtod(optarg, NULL); break;
	
	/* I/O functions */
		case 'f':				/* Tagged data format */
			format = optarg; break;
		case 'd':				/* Output data only */
			output_only = 1; break;
		case 'm':				/* Models list */
			list_models = 1; break;
		case 'p':				/* Plot curves */
			{
				int i;
				
				for (i = 0; circuit[i].name; i++)
					if (!strcmp(optarg, circuit[i].name)) plot_curves = i;
				
				if (!plot_curves) usage();
			}
			break;
		case 'w':				/* Waveform analysis */
			waveform = 1; break;
	
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
		case 'M':				/* User model */
			user_model = optarg; break;
	
	/* Default response */
		default:
			usage();
	}
	
	if (argc != optind) usage();
	if (verbose) fprintf(stderr, "%s\n", usage_msg[0]);
	
	
	/* Read data */
	if (format) d = read_tagged_data(stdin, format, &n);
	else d = read_data(stdin, &n);
	
	if (output_only) { write_data(stdout, d, n); exit(0); }
	
	/* Build model */
	if (user_model) { m = read_model(user_model, d, n); }
	else { m = best_model((list_models ? stdout : NULL), d, n); }
	
	/* Produce requested output */
	if (plot_curves) (*(circuit[plot_curves].curves))(stdout, m, d, n, 0, 0, 0, 0);
	
	return 0;
}
