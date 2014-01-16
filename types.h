/****************************************************************
 *
 * types.h: potfit types definition file
 *
 ****************************************************************
 *
 * Copyright 2002-2013
 *	Institute for Theoretical and Applied Physics
 *	University of Stuttgart, D-70550 Stuttgart, Germany
 *	http://potfit.sourceforge.net/
 *
 ****************************************************************
 *
 *   This file is part of potfit.
 *
 *   potfit is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   potfit is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

typedef enum Param_T {
  PARAM_STR,
  PARAM_INT,
  PARAM_DOUBLE
} param_t;

typedef struct {
  double x;
  double y;
  double z;
} vector;

/* This is the order of VASP for stresses */
typedef struct {
  double xx;
  double yy;
  double zz;
  double xy;
  double yz;
  double zx;
} sym_tens;

typedef struct {
  /* neighbor properties */
  int   type;			/* type of neighboring atom */
  int   nr;			/* number of neighboring atom */
  double r;			/* r */
  double r2;			/* r^2 */
  double inv_r;			/* 1/r */
  vector dist;			/* real distance vector */
  vector dist_r;		/* distance vector divided by r = normalized distance vector */

  /* data to access the spline tables at the correct position */
  int   slot[SLOTS];		/* the slot, belonging to the neighbor distance */
  double shift[SLOTS];		/* how far into the slot we have to go, in [0..1] */
  double step[SLOTS];		/* step size */
  int   col[SLOTS];		/* coloumn of interaction for this neighbor */

#ifdef ADP
  sym_tens sqrdist;		/* real squared distance */
  double u_val, u_grad;		/* value and gradient of u(r) */
  double w_val, w_grad;		/* value and gradient of w(r) */
#endif

#ifdef COULOMB
  double fnval_el;		/* stores tail of electrostatic potential */
  double grad_el;		/* stores tail of first derivative of electrostatic potential */
  double ggrad_el;		/* stores tail of second derivative of electrostatic potential */
#endif

#ifdef THREEBODY
  double f;			/* value of the cutoff function f_c */
  double df;			/* derivative of the cutoff function f_c */
  int   ijk_start;		/* index of the first entry for angular part */
#endif

#ifdef MEAM
  double drho;
#endif

#ifdef TERSOFF
  vector dzeta;
#endif
} neigh_t;

#ifdef THREEBODY
typedef struct {
  double cos;
#ifdef MEAM
  int   slot;
  double shift;
  double step;
  double g;
  double dg;
#endif
} angle_t;
#endif

#ifdef STIWEB
/* pointers to access Stillinger-Weber parameters directly */
typedef struct {
  int   init;
  double **A;
  double **B;
  double **p;
  double **q;
  double **delta;
  double **a1;
  double **gamma;
  double **a2;
  double ****lambda;
} sw_t;
#endif

#ifdef TERSOFF
#ifndef TERSOFFMOD
/* pointers to access Tersoff parameters directly */
typedef struct {
  int   init;
  double **A;
  double **B;
  double **lambda;
  double **mu;
  double **gamma;
  double **n;
  double **c;
  double **d;
  double **h;
  double **S;
  double **R;
  double **chi;
  double **omega;
  double *c2;
  double *d2;
  double one;
} tersoff_t;
#else
/* pointers to access Tersoff parameters directly */
typedef struct {
  int   init;
  double **A;
  double **B;
  double **lambda;
  double **mu;
  double **eta;
  double **delta;
  double **alpha;
  double **beta;
  double **c1;
  double **c2;
  double **c3;
  double **c4;
  double **c5;
  double **h;
  double **R1;
  double **R2;
} tersoff_t;
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */

typedef struct {
  int   type;
  int   num_neigh;
  vector pos;
  vector force;
  double absforce;
  int   conf;			/* Which configuration... */

#ifdef CONTRIB
  int   contrib;		/* Does this atom contribute to the error sum? */
#endif

#if defined EAM || defined ADP || defined MEAM
  double rho;			/* embedding electron density */
  double gradF;			/* gradient of embedding fn. */
#ifdef TBEAM
  double rho_s;			/* embedding electron density */
  double gradF_s;		/* gradient of embedding fn. */
#endif				/* TBEAM */
#endif				/* EAM || ADP || MEAM */

#ifdef ADP
  vector mu;
  sym_tens lambda;
  double nu;
#endif

#ifdef DIPOLE
  vector E_stat;		/* static field-contribution */
  vector p_sr;			/* short-range dipole moment */
  vector E_ind;			/* induced field-contribution */
  vector p_ind;			/* induced dipole moment */
  vector E_old;			/* stored old induced field */
  vector E_tot;			/* temporary induced field */
#endif

#ifdef THREEBODY
  int   num_angles;
#ifdef MEAM
  double rho_eam;		/* Store EAM density */
#endif
#endif

  neigh_t *neigh;		/* dynamic array for neighbors */
#ifdef THREEBODY
  angle_t *angle_part;		/* dynamic array for angular neighbors */
#endif
} atom_t;

typedef struct {
  int   len;			/* total length of the table */
  int   idxlen;			/* number of changeable potential values */
  int   ncols;			/* number of columns */
  double *begin;		/* first value in the table */
  double *end;			/* last value in the table */
  double *step;			/* table increment */
  double *invstep;		/* inverse of increment */
  int  *first;			/* index of first entry */
  int  *last;			/* index of last entry */
  double *xcoord;		/* the x-coordinates of sampling points */
  double *table;		/* the actual data */
  double *d2tab;		/* second derivatives of table data for spline int */
  int  *idx;			/* indirect indexing */
} pot_table_t;

#ifdef APOT
/* function pointer for analytic potential evaluation */
typedef void (*fvalue_pointer) (double, double *, double *);

typedef struct {
  /* potentials */
  int   number;			/* number of analytic potentials */
  int   invar_pots;		/* number of invariant analytic potentials */
  int  *idxpot;			/* indirect index for potentials */
  char **names;			/* name of analytic potentials */
  double *begin;		/* starting position of potential */
  double *end;			/* end position of potential = cutoff radius */
  int  *n_par;			/* number of parameters for analytic potential */

  /* parameters */
  int   total_par;		/* total number of parameters for all potentials */
#ifdef COULOMB
  int   total_ne_par;		/* total number of non-electrostatic parameters */
#endif
  int  *idxparam;		/* indirect index for potential parameters */
  int **invar_par;		/* array of invariant parameters */
  char ***param_name;		/* name of parameters */
  double **pmin;		/* minimum values for parameters */
  double **values;		/* parameter values for analytic potentials */
  double **pmax;		/* maximum values for parameters */

  /* global parameters */
  int   globals;		/* number of global parameters */
  int  *n_glob;			/* number of global parameter usage */
  int ***global_idx;		/* index of global parameters */

#ifdef PAIR
  double *chempot;		/* chemical potentials */
#endif

#ifdef COULOMB
  double *ratio;		/* stoichiometric ratio */
  double *charge;		/* charges */
  double last_charge;		/* last charge determined on the basis of charge neutrality */
  double *dp_kappa;		/* parameter kappa */
  int   sw_kappa;		/* switch for kappa-optimization */
#endif
#ifdef DIPOLE
  double *dp_alpha;		/* polarisability */
  double *dp_b;			/* parameter for short-range-dipole-moment */
  double *dp_c;			/* parameter for short-range-dipole-moment */
#endif

#ifdef STIWEB
  sw_t  sw;
#endif

#ifdef TERSOFF
  tersoff_t tersoff;
#endif

  fvalue_pointer *fvalue;	/* function pointers for analytic potentials */
} apot_table_t;

typedef struct {
  char **name;			/* identifier of the potential */
  int  *n_par;			/* number of parameters */
  fvalue_pointer *fvalue;	/* function pointer */
} function_table_t;

#endif /* APOT */
