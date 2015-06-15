/****************************************************************
 *
 * types.h: potfit types definition file
 *
 ****************************************************************
 *
 * Copyright 2002-2014
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

/****************************************************************
 *
 *  plain old vector
 *
 ****************************************************************/

typedef struct {
  double x;
  double y;
  double z;
} vector;

/****************************************************************
 *
 *  symmetric tensor in the same order as used by VASP
 *
 ****************************************************************/

typedef struct {
  double xx;
  double yy;
  double zz;
  double xy;
  double yz;
  double zx;
} sym_tens;

/****************************************************************
 *
 *  neighbor table (each atom has one for each neighbor)
 *
 ****************************************************************/

typedef struct {
  /* neighbor properties */
  int type;      /* type of neighboring atom */
  int nr;        /* number of neighboring atom */
  double r;      /* r */
  double r2;     /* r^2 */
  double inv_r;  /* 1/r */
  vector dist;   /* real distance vector */
  vector dist_r; /* distance vector divided by r = normalized distance vector */

  /* data to access the spline tables at the correct position */
  int slot[SLOTS];     /* the slot, belonging to the neighbor distance */
  double shift[SLOTS]; /* how far into the slot we have to go, in [0..1] */
  double step[SLOTS];  /* step size */
  int col[SLOTS];      /* coloumn of interaction for this neighbor */

#if defined(ADP)
  sym_tens sqrdist;     /* real squared distance */
  double u_val, u_grad; /* value and gradient of u(r) */
  double w_val, w_grad; /* value and gradient of w(r) */
#endif

#if defined(COULOMB)
  double fnval_el;      /* stores tail of electrostatic potential */
  double grad_el;       /* stores tail of first derivative of electrostatic potential */
  double ggrad_el;      /* stores tail of second derivative of electrostatic potential */
#endif

#if defined(THREEBODY)
  double f;             /* value of the cutoff function f_c */
  double df;            /* derivative of the cutoff function f_c */
  int ijk_start;        /* index of the first entry for angular part */
#endif

#if defined(MEAM)
  double drho;
#endif

#if defined(TERSOFF)
  vector dzeta;
#endif
} neigh_t;

/****************************************************************
 *
 *  angular neighbor table (each atom has one for each triple of neighbors)
 *
 ****************************************************************/

#if defined(THREEBODY)
typedef struct {
  double cos;
#if defined(MEAM)
  int slot;
  double shift;
  double step;
  double g;
  double dg;
#endif // MEAM
} angle_t;
#endif // THREEBODY

/****************************************************************
 *
 *  pointers to access Stillinger-Weber parameters directly
 *
 ****************************************************************/

#if defined(STIWEB)
typedef struct {
  int init;
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

/****************************************************************
 *
 *  pointers to access Tersoff/TersoffMOD parameters directly
 *
 ****************************************************************/

#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
typedef struct {
  int init;
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
typedef struct {
  int init;
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
#endif // !TERSOFFMOD
#endif // TERSOFF

/****************************************************************
 *
 *  atom data
 *
 ****************************************************************/

typedef struct {
  int type;
  int num_neigh;
  vector pos;
  vector force;
  double absforce;
  int conf; /* Which configuration... */

#if defined(CONTRIB)
  int contrib; /* Does this atom contribute to the error sum? */
#endif

#if defined(EAM) || defined(ADP) || defined(MEAM)
  double rho;   /* embedding electron density */
  double gradF; /* gradient of embedding fn. */
#if defined(TBEAM)
  double rho_s;   /* embedding electron density */
  double gradF_s; /* gradient of embedding fn. */
#endif // TBEAM
#endif // EAM || ADP || MEAM

#if defined(ADP)
  vector mu;
  sym_tens lambda;
  double nu;
#endif

#if defined(DIPOLE)
  vector E_stat; /* static field-contribution */
  vector p_sr;   /* short-range dipole moment */
  vector E_ind;  /* induced field-contribution */
  vector p_ind;  /* induced dipole moment */
  vector E_old;  /* stored old induced field */
  vector E_tot;  /* temporary induced field */
#endif

#if defined(THREEBODY)
  int num_angles;
#if defined(MEAM)
  double rho_eam; /* Store EAM density */
#endif // MEAM
#endif // THREEBODY

  neigh_t *neigh; /* dynamic array for neighbors */
#if defined(THREEBODY)
  angle_t *angle_part; /* dynamic array for angular neighbors */
#endif
} atom_t;

/****************************************************************
 *
 *  potential table: holds tabulated potential data
 *
 ****************************************************************/

typedef struct {
  int len;         /* total length of the table */
  int idxlen;      /* number of changeable potential values */
  int ncols;       /* number of columns */
  double *begin;   /* first value in the table */
  double *end;     /* last value in the table */
  double *step;    /* table increment */
  double *invstep; /* inverse of increment */
  int *first;      /* index of first entry */
  int *last;       /* index of last entry */
  double *xcoord;  /* the x-coordinates of sampling points */
  double *table;   /* the actual data */
  double *d2tab;   /* second derivatives of table data for spline int */
  int *idx;        /* indirect indexing */
} pot_table_t;

#if defined(APOT)
/* function pointer for analytic potential evaluation */
typedef void (*fvalue_pointer)(double, double *, double *);

/****************************************************************
 *
 *  potential table: holds tabulated potential data
 *
 ****************************************************************/

typedef struct {
  /* potentials */
  int number;     /* number of analytic potentials */
  int invar_pots; /* number of invariant analytic potentials */
  int *idxpot;    /* indirect index for potentials */
  char **names;   /* name of analytic potentials */
  double *begin;  /* starting position of potential */
  double *end;    /* end position of potential = cutoff radius */
  int *n_par;     /* number of parameters for analytic potential */

  /* parameters */
  int total_par; /* total number of parameters for all potentials */
#if defined(COULOMB)
  int total_ne_par; /* total number of non-electrostatic parameters */
#endif
  int *idxparam;      /* indirect index for potential parameters */
  int **invar_par;    /* array of invariant parameters */
  char ***param_name; /* name of parameters */
  double **pmin;      /* minimum values for parameters */
  double **values;    /* parameter values for analytic potentials */
  double **pmax;      /* maximum values for parameters */

  /* global parameters */
  int globals;       /* number of global parameters */
  int *n_glob;       /* number of global parameter usage */
  int ***global_idx; /* index of global parameters */

#if defined(PAIR)
  double *chempot; /* chemical potentials */
#endif

#if defined(COULOMB)
  double *ratio;      /* stoichiometric ratio */
  double *charge;     /* charges */
  double last_charge; /* last charge determined on the basis of charge neutrality */
  double *dp_kappa;   /* parameter kappa */
  int sw_kappa;       /* switch for kappa-optimization */
#endif
#if defined(DIPOLE)
  double *dp_alpha; /* polarisability */
  double *dp_b;     /* parameter for short-range-dipole-moment */
  double *dp_c;     /* parameter for short-range-dipole-moment */
#endif

#if defined(STIWEB)
  sw_t sw;
#endif

#if defined(TERSOFF)
  tersoff_t tersoff;
#endif

  fvalue_pointer *fvalue; /* function pointers for analytic potentials */
} apot_table_t;






















typedef struct {
  /* system variables */
  int myid;     /* Which process am I? (0 if serial) */
  int num_cpus; /* How many cpus are there */

  /* general MPI variables */
  int firstatom;
  int firstconf;
  int myatoms;
  int myconf;

#ifdef MPI
  /* MPI datatypes */
  MPI_Datatype MPI_ATOM;
  MPI_Datatype MPI_NEIGH;
#ifdef THREEBODY
  MPI_Datatype MPI_ANGL;
#endif
  MPI_Datatype MPI_STENS;
  MPI_Datatype MPI_VECTOR;

  /* variables needed for atom distribution with MPI */
  int *atom_dist;
  int *atom_len;
  int *conf_dist;
  int *conf_len;
#endif /* MPI */
} potfit_mpi_config;

typedef struct {
   atom_t *atoms;           /* atoms array */
   atom_t *conf_atoms;      /* Atoms in configuration */
   char **elements;         /* element names from vasp2force */
   int have_elements;
   int **na_type;           /* number of atoms per type */
   int *cnfstart;           /* Nr. of first atom in config */
   int *conf_uf;
  #ifdef STRESS
   int *conf_us;
  #endif /* STRESS */
   int *inconf;             /* Nr. of atoms in each config */
   int *useforce;           /* Should we use forces */
  #ifdef STRESS
   int *usestress;          /* Should we use stresses */
  #endif /* STRESS */
   int natoms;      /* number of atoms */
   int nconf;       /* number of configurations */
  #ifdef CONTRIB
   int have_contrib_box;    /* do we have a box of contrib. atoms? */
   int n_spheres;   /* number of spheres of contrib. atoms */
   double *r_spheres;       /* radii of the spheres of contrib. atoms */
  #endif /* CONTRIB */

   double *coheng;          /* Cohesive energy for each config */
   double *conf_vol;
   double *conf_weight;     /* weight of configuration */
   double *force_0;         /* the forces we aim at */
   double *rcut;
   double *rmin;
   double *volume;          /* volume of cell */
   double rcutmin;      /* minimum of all cutoff values */
   double rcutmax;        /* maximum of all cutoff values */
  #ifdef STRESS
   sym_tens *conf_stress;
   sym_tens *stress;        /* Stresses in each config */
  #endif /* STRESS */
   vector box_x, box_y, box_z;
  #ifdef CONTRIB
   vector cbox_o;           /* origin of box of contrib. atoms */
   vector cbox_a, cbox_b, cbox_c;   /* box vectors for box of contrib. atoms */
   vector *sphere_centers;  /* centers of the spheres of contrib. atoms */
  #endif /* CONTRIB */
   vector tbox_x, tbox_y, tbox_z;
} potfit_configurations;

typedef struct {
  /* general settings (from parameter file) */
   int imdpotsteps;      /* resolution of IMD potential */
   int ntypes;     /* number of atom types */
   int opt;         /* optimization flag */
   int rng_seed;        /* seed for RNG */
   int usemaxch;    /* use maximal changes file */
   int write_output_files;
   int write_lammps_files;
   int write_pair;
   int writeimd;
   int write_lammps;        /* write output also in LAMMPS format */
  #ifdef EVO
   double evo_threshold;
  #else /* EVO */
   char* anneal_temp;
  #endif /* EVO */
   double eweight;
   double sweight;
   double extend; /* how far should one extend imd pot */
  #ifdef APOT
   int compnodes;   /* how many additional composition nodes */
   int enable_cp;   /* switch chemical potential on/off */
   double apot_punish_value;
   double plotmin;        /* minimum for plotfile */
  #endif /* APOT */
  double global_cell_scale;      /* global scaling parameter */
} potfit_parameters;

typedef struct {
  char* config;        /* file with atom configuration */
  char* distfile;      /* file for distributions */
  char* endpot;        /* file for end potential */
  char* flagfile;      /* break if file exists */
  char* imdpot;        /* file for IMD potential */
  char* maxchfile;     /* file with maximal changes */
  char* output_prefix; /* prefix for all output files */
  char* output_lammps; /* lammps output files */
  char* plotfile;      /* file for plotting */
  char* plotpointfile; /* write points for plotting */
  char* startpot;      /* file with start potential */
  char* tempfile;      /* backup potential file */
} potfit_filenames;

typedef struct {
  /* memory management */
  char **pointer_names;
  int num_pointers;
  void **pointers;
  double *u_address;
} potfit_memory;

typedef struct {
  /* potential variables */
  int *gradient;           /* Gradient of potential fns.  */
  int *invar_pot;
  int format;     /* format of potential table */
  int have_grad;   /* Is gradient specified?  */
  int have_invar;  /* Are invariant pots specified?  */
  #ifdef APOT
  int *smooth_pot;
  int cp_start;    /* cp in opt_pot.table */
  int global_idx;  /* index for global parameters in opt_pot table */
  int global_pot;  /* number of "potential" for global parameters */
  int have_globals;        /* do we have global parameters? */
  double *calc_list;       /* list of current potential in the calc table */
  double *compnodelist;    /* list of the composition nodes */
  #endif /* APOT */

  /* potential tables */
  pot_table_t opt_pot;     /* potential in the internal representation used for minimisation */
  pot_table_t calc_pot;    /* the potential table used for force calculations */
  #ifdef APOT
  apot_table_t apot_table; /* potential in analytic form */
  #endif /* APOT */
} potfit_potentials;

typedef struct {
  /* optimization variables */
   int fcalls;
   int mdim;
   int ndim;
   int ndimtot;
   int paircol;     /* How many columns for pair potential */
   double d_eps;

  /* pointers for force-vector */
   int energy_p;    /* pointer to energies */
  #ifdef STRESS
   int stress_p;    /* pointer to stresses */
  #endif /* STRESS */
  #if defined EAM || defined ADP || defined MEAM
   int dummy_p;     /* pointer to dummy constraints */
   int limit_p;     /* pointer to limiting constraints */
  #endif /* EAM || ADP || MEAM */
  #ifdef APOT
   int punish_par_p;        /* pointer to parameter punishment contraints */
   int punish_pot_p;        /* pointer to potential punishment constraints */
  #endif /* APOT */
} potfit_calculation;

typedef struct {
  /* misc. stuff - has to belong somewhere */
   int *idx;
   int init_done;
   int plot;        /* plot output flag */
  #if defined EAM || defined ADP || defined MEAM
   double *lambda;          /* embedding energy slope... */
  #endif
   double *maxchange;       /* Maximal permissible change */
  /* variables needed for electrostatic options */
  #ifdef COULOMB
   double dp_eps; /* ??? in eV A */
   double dp_cut;        /* cutoff-radius for long-range interactions */
  #endif /* COULOMB */
  #ifdef DIPOLE
   double dp_tol;      /* dipole iteration precision */
   double dp_mix */
  #endif /* DIPOLE */
} potfit_unknown;

#endif /* APOT */
