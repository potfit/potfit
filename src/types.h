/****************************************************************
 *
 * types.h: potfit types definition file
 *
 ****************************************************************
 *
 * Copyright 2002-2017 - the potfit development team
 *
 * https://www.potfit.net/
 *
 ****************************************************************
 *
 * This file is part of potfit.
 *
 * potfit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * potfit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with potfit; if not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************/

#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

#define POTFIT_SUCCESS 0
#define POTFIT_ERROR_MPI_CLEAN_EXIT 1
#define POTFIT_ERROR 2

typedef enum {
  POTENTIAL_FORMAT_ANALYTIC = 0,
  POTENTIAL_FORMAT_TABULATED_EQ_DIST = 3,
  POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST = 4,
  POTENTIAL_FORMAT_KIM = 5,
  POTENTIAL_FORMAT_UNKNOWN = 999
} POTENTIAL_FORMAT;

// plain old vector

typedef struct {
  double x;
  double y;
  double z;
} vector;

// symmetric tensor in the same order as used by VASP

typedef struct {
  double xx;
  double yy;
  double zz;
  double xy;
  double yz;
  double zx;
} sym_tens;

//  neighbor table (each atom has one for each neighbor)

typedef struct {
  // neighbor properties
  int type;       // type of neighboring atom
  int nr;         // number of neighboring atom
  double r;       // r
  double r2;      // r^2
  double inv_r;   // 1/r
  vector dist;    // real distance vector
  vector dist_r;  // distance vector divided by r = normalized distance vector

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
  double fnval_el; /* stores tail of electrostatic potential */
  double
      grad_el; /* stores tail of first derivative of electrostatic potential */
  double ggrad_el; /* stores tail of second derivative of electrostatic
                      potential */
#endif

#if defined(THREEBODY)
  double f;      /* value of the cutoff function f_c */
  double df;     /* derivative of the cutoff function f_c */
  int ijk_start; /* index of the first entry for angular part */
#endif

#if defined(MEAM)
  double drho;
#endif

#if defined(TERSOFF)
  vector dzeta;
#endif
} neigh_t;

// angular neighbor table (each atom has one for each triple of neighbors)

#if defined(THREEBODY)
typedef struct {
  double cos;
#if defined(MEAM)
  int slot;
  double shift;
  double step;
  double g;
  double dg;
#endif  // MEAM
} angle_t;
#endif  // THREEBODY

// pointers to access Stillinger-Weber parameters directly

#if defined(STIWEB)
typedef struct {
  int init;
  double** A;
  double** B;
  double** p;
  double** q;
  double** delta;
  double** a1;
  double** gamma;
  double** a2;
  double**** lambda;
} sw_t;
#endif

//  pointers to access Tersoff/TersoffMOD parameters directly

#if defined(TERSOFF)
typedef struct {
  int init;
  double** A;
  double** B;
  double** lambda;
  double** mu;
#if defined(TERSOFFMOD)
  double** eta;
  double** delta;
  double** alpha;
  double** beta;
  double** c1;
  double** c2;
  double** c3;
  double** c4;
  double** c5;
  double** h;
  double** R1;
  double** R2;
#else
  double** gamma;
  double** n;
  double** c;
  double** d;
  double** h;
  double** S;
  double** R;
  double** chi;
  double** omega;
  double* c2;
  double* d2;
  double one;
#endif  // TERSOFFMOD
} tersoff_t;
#endif  // TERSOFF

// atom data table

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
#endif            // TBEAM
#endif            // EAM || ADP || MEAM

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
#endif            // MEAM
#endif            // THREEBODY

  neigh_t* neigh; /* dynamic array for neighbors */
#if defined(THREEBODY)
  angle_t* angle_part; /* dynamic array for angular neighbors */
#endif
} atom_t;

// potential table: holds tabulated potential data

typedef struct {
  int len;         /* total length of the table */
  int idxlen;      /* number of changeable potential values */
  int ncols;       /* number of columns */
  double* begin;   /* first value in the table */
  double* end;     /* last value in the table */
  double* step;    /* table increment */
  double* invstep; /* inverse of increment */
  int* first;      /* index of first entry */
  int* last;       /* index of last entry */
  double* xcoord;  /* the x-coordinates of sampling points */
  double* table;   /* the actual data */
  double* d2tab;   /* second derivatives of table data for spline int */
  int* idx;        /* indirect indexing */
} pot_table_t;

#if defined(APOT)
// function pointer for analytic potential evaluation
typedef void (*fvalue_pointer)(const double, const double*, double*);

// potential table: holds analytic potential data

typedef struct {
  /* potentials */
  int number;     /* number of analytic potentials */
  int invar_pots; /* number of invariant analytic potentials */
  int* idxpot;    /* indirect index for potentials */
  char** names;   /* name of analytic potentials */
  double* begin;  /* starting position of potential */
  double* end;    /* end position of potential = cutoff radius */
  int* n_par;     /* number of parameters for analytic potential */

  /* parameters */
  int total_par; /* total number of parameters for all potentials */
#if defined(COULOMB)
  int total_ne_par; /* total number of non-electrostatic parameters */
#endif
  int* idxparam;      /* indirect index for potential parameters */
  int** invar_par;    /* array of invariant parameters */
  char*** param_name; /* name of parameters */
  double** pmin;      /* minimum values for parameters */
  double** values;    /* parameter values for analytic potentials */
  double** pmax;      /* maximum values for parameters */

  /* global parameters */
  int globals;       /* number of global parameters */
  int* n_glob;       /* number of global parameter usage */
  int*** global_idx; /* index of global parameters */

#if defined(PAIR)
  double* chempot; /* chemical potentials */
#endif

#if defined(COULOMB)
  double* ratio;      /* stoichiometric ratio */
  double* charge;     /* charges */
  double last_charge; /* last charge determined on the basis of charge
                         neutrality */
  double* dp_kappa;   /* parameter kappa */
  int sw_kappa;       /* switch for kappa-optimization */
#endif
#if defined(DIPOLE)
  double* dp_alpha; /* polarisability */
  double* dp_b;     /* parameter for short-range-dipole-moment */
  double* dp_c;     /* parameter for short-range-dipole-moment */
#endif

#if defined(STIWEB)
  sw_t sw;
#endif

#if defined(TERSOFF)
  tersoff_t tersoff;
#endif

  fvalue_pointer* fvalue; /* function pointers for analytic potentials */
} apot_table_t;

#endif  // APOT

// potfit_calculation table: holds data for force calculation

typedef struct {
  int fcalls;  /* number of force calculations */
  int mdim;    /* total number of entries in force vector */
  int ndim;    /* number of free optimization parameters in force vector */
  int ndimtot; /* total number of optimization parameters in force vector */
  int paircol; /* How many columns for pair potential ( ntypes*(ntypes+1)/2 ) */

  double d_eps;   /* abortion criterion for powell_lsq */
  double* force;  //

  struct {
    double* vecu_bracket;
    double* vecu_brent;
    double* f_vec3;
  } linmin;

  /* pointers for force-vector */
  int energy_p; /* offset of energies in force vector */
#if defined(STRESS)
  int stress_p; /* offset of stresses in force vector */
#endif          // STRESS
#if defined(EAM) || defined(ADP) || defined(MEAM)
  int dummy_p; /* pointer to dummy constraints */
  int limit_p; /* pointer to limiting constraints */
#endif         // EAM || ADP || MEAM
#if defined(APOT)
  int punish_par_p; /* pointer to parameter punishment contraints */
  int punish_pot_p; /* pointer to potential punishment constraints */
#else
  double* maxchange;  // Maximal permissible change
#if defined EAM || defined ADP || defined MEAM
  double* lambda;  // embedding energy slope...
#endif
#endif  // APOT
} potfit_calculation;

// potfit_configurations table: holds reference configurations data

typedef struct {
  int natoms; /* total number of atoms */
  int nconf;  /* total number of configurations */

  atom_t* atoms;      /* atoms array */
  atom_t* conf_atoms; /* Atoms in configuration */

  char const** elements; /* element names from configuration files */

  int** na_type; /* number of atoms per atom type */

  int* cnfstart; /* index of first atom in each config */
  int* inconf;   /* number of atoms in each config */
  int* conf_uf;  /* local array of "use forces in config X" */
  int* useforce; /* global array of "use forces in config X" */

  double* coheng;      /* cohesive energy for each config */
  double* conf_vol;    /* local volume of configuration cell */
  double* volume;      /* global volume of configuration cell */
  double* conf_weight; /* weight of each configuration */
  double* force_0;     /* ab-initio reference forces */

  double* rcut; /* cutoff radius for each atom type */
  double* rmin; /* minimum distance for each atom type */

  double rcutmin; /* minimum of all cutoff values */
  double rcutmax; /* maximum of all cutoff values */

#if defined(STRESS)
  int* conf_us;          /* local array of "use stresses in config X" */
  int* usestress;        /* global array of "use stresses in config X" */
  sym_tens* conf_stress; /* local stress of each configuration */
  sym_tens* stress;      /* global stress of each configuration */
#endif                   // STRESS

// variables needed for electrostatic options
#if defined(COULOMB)
  double dp_cut;  // cutoff-radius for long-range interactions
#endif            // COULOMB
#if defined(DIPOLE)
  double dp_tol;  // dipole iteration precision
  double dp_mix;  // ???
#endif            // DIPOLE
} potfit_configurations;

// potfit_filenames: holds all kinds of filenames

typedef struct {
  const char* config;        /* file with atom configuration */
  const char* distfile;      /* file for distributions */
  const char* endpot;        /* file for end potential */
  const char* flagfile;      /* break if file exists */
  const char* imdpot;        /* file for IMD potential */
  const char* maxchfile;     /* file with maximal changes */
  const char* output_prefix; /* prefix for all output files */
  const char* output_lammps; /* lammps output files */
  const char* plotfile;      /* file for plotting */
  const char* plotpointfile; /* write points for plotting */
  const char* startpot;      /* file with start potential */
  const char* tempfile;      /* backup potential file */
} potfit_filenames;

// potfit_mpi_config: holds information needed for MPI calculation

typedef struct {
  int init_done;  // potfit setup completed
  int myid;       // index of current process
  int num_cpus;   // total numer of processes

  int firstatom; /* index of first atom for this process */
  int firstconf; /* index of first configuration for this process */
  int myatoms;   /* number of atoms for this process */
  int myconf;    /* number of configurations for this process */

#if defined(MPI)
  int* atom_dist; /* atom distribution for each process (starting index) */
  int* atom_len;  /* atom distribution for each process (number of atoms) */
  int* conf_dist; /* config distribution for each process (starting index) */
  int* conf_len;  /* config distribution for each process (number of configs) */

  /* MPI datatypes */
  MPI_Datatype MPI_ATOM;
  MPI_Datatype MPI_NEIGH;
#if defined(THREEBODY)
  MPI_Datatype MPI_ANGL;
#endif  // THREEBODY
  MPI_Datatype MPI_STENS;
  MPI_Datatype MPI_VECTOR;
#endif  // MPI
} potfit_mpi_config;

// potfit_parameters: holds information from parameter file

typedef struct {
  int imdpotsteps; /* resolution of IMD potential */
  int ntypes;      /* number of atom types */
  int opt;         /* optimization flag */
  int rng_seed;    /* seed for RNG */
  int usemaxch;    /* use maximal changes file */

  int plot;  // plot output flag

  int write_output_files;
  int write_lammps_files;
  int write_pair;
  int writeimd;
  int write_lammps; /* write output also in LAMMPS format */

#if defined(EVO)
  double evo_threshold;
#else
  const char* anneal_temp;
#endif  // EVO
  double eweight;
  double sweight;
  double extend; /* how far should one extend imd pot */
#if defined(APOT)
  int compnodes; /* how many additional composition nodes */
  int enable_cp; /* switch chemical potential on/off */
  double apot_punish_value;
  double plotmin;           /* minimum for plotfile */
#endif                      // APOT
  double global_cell_scale; /* global scaling parameter */
} potfit_parameters;

// potfit_potentials: holds information from potential file

typedef struct {
  const char* interaction_name;
  int* gradient; /* Gradient of potential fns.  */
  int* invar_pot;
  POTENTIAL_FORMAT format_type; /* format of potential table */
  int have_invar;               /* Are invariant pots specified?  */
#if defined(APOT)
  int* smooth_pot;
  int cp_start;         /* cp in opt_pot.table */
  int global_idx;       /* index for global parameters in opt_pot table */
  int global_pot;       /* number of "potential" for global parameters */
  int have_globals;     /* do we have global parameters? */
  double* calc_list;    /* list of current potential in the calc table */
  double* compnodelist; /* list of the composition nodes */
#endif                  // APOT

  /* potential tables */
  pot_table_t opt_pot;  /* potential in the internal representation used for
  minimisation */
  pot_table_t calc_pot; /* the potential table used for force calculations */
#if defined(APOT)
  apot_table_t apot_table; /* potential in analytic form */
#endif                     // APOT
} potfit_potentials;

#if defined(KIM)

/* sturct of all free parameters of data type `double' */
typedef struct {
  char** name;
  double** value;       /* the pointers to parameters */
  int* rank;
  int** shape;
  int Nparam;           /* number of free parameters */
  double** nestedvalue; /* nest all values(pointers) of optimizable parameters into this long list */
  int Nnestedvalue;     /* the length of nestedvalue */
} FreeParamType;

/* varialbes */
typedef struct {
  void** pkimObj;                    /* pointers to kim objects */
  FreeParamType* FreeParamAllConfig; /* free param struct for all configurations */
  char kim_model_name[255]; /* kim model name (read in from input )*/
  char NBC_method[64];      /* neighbor list and boundary conditions */
  int kim_model_has_forces;
  int kim_model_has_virial;
  int is_half_neighbors;    /* using half neighbor list? 1 = half, 0 = full */
  int kim_model_has_energy; /* flag, whether KIM model has the routine to
                             * compute energy, forces, and virial */
  int num_opt_param;        /* number of optimizalbe params( read in from input) */
  char** name_opt_param;    /* optimizable parameter names (read in from input) */
  int* size_opt_param;      /* size of each parameter (number of values each
                             * parameter name represents) */
  double* box_side_len;     /* box_side_len is used to enable MI_OPBC in KIM. */
} potfit_kim;

#endif // KIM

#endif  // TYPES_H_INCLUDED
