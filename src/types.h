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

typedef struct
{
  double x;
  double y;
  double z;
} vector;

/****************************************************************
 *
 *  symmetric tensor in the same order as used by VASP
 *
 ****************************************************************/

typedef struct
{
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

typedef struct
{
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
  double fnval_el; /* stores tail of electrostatic potential */
  double grad_el;  /* stores tail of first derivative of electrostatic potential */
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

/****************************************************************
 *
 *  angular neighbor table (each atom has one for each triple of neighbors)
 *
 ****************************************************************/

#if defined(THREEBODY)
typedef struct
{
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

/****************************************************************
 *
 *  pointers to access Stillinger-Weber parameters directly
 *
 ****************************************************************/

#if defined(STIWEB)
typedef struct
{
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

/****************************************************************
 *
 *  pointers to access Tersoff/TersoffMOD parameters directly
 *
 ****************************************************************/

#if defined(TERSOFF)
typedef struct
{
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

/****************************************************************
 *
 *  atom data
 *
 ****************************************************************/

typedef struct
{
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

/****************************************************************
 *
 *  potential table: holds tabulated potential data
 *
 ****************************************************************/

typedef struct
{
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
/* function pointer for analytic potential evaluation */
typedef void (*fvalue_pointer)(double, double*, double*);

/****************************************************************
 *
 *  potential table: holds analytic potential data
 *
 ****************************************************************/
// Addition by AI 17/7/2015
#ifdef LMP
typedef struct { // lattice_t
  double xx;
  double xy;
  double xz;
  double yx;
  double yy;
  double yz;
  double zx;
  double zy;
  double zz;
} lattice_t;
#endif
// End of addition

typedef struct
{
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

//<<<<<<< HEAD
// Added by AI 17/7/2015
#if defined (LMP)
  int rf_sr0; // Num of Params for 0-Body Interaction (General)
  int rf_sr1; // Num of Params for 1-Body Interaction (Atom)
  int rf_sr2; // Num of Params for 2-Body Interaction (Bond)
  int rf_sr3; // Num of Params for 3-Body Interaction (Angle)
  int rf_sr4; // Num of Params for 4-Body Interaction (Torsion)
  int rf_srO; // Num of Params for 2-Body Interaction (Off-Diagonal)
  int rf_srH; // Num of Params for 3-Body Interaction (Hydrogen-Bond)

  int rf_comb0; // Num of Combination for 0-Body Interaction (General)
  int rf_comb1; // Num of Combination for 1-Body Interaction (Atom)
  int rf_comb2; // Num of Combination for 2-Body Interaction (Bond)
  int rf_comb3; // Num of Combination for 3-Body Interaction (Angle)
  int rf_comb4; // Num of Combination for 4-Body Interaction (Torsion)
  int rf_combO; // Num of Combination for 2-Body Interaction (Off-Diagonal)
  int rf_combH; // Num of Combination for 3-Body Interaction (Hydrogen-Bond)

// 0-Body (General)
  double *vpar[39];

// 1-Body (Atom)
  double *rat;     // Correction for Overcoordination
  double *aval;    // Overcoordination Energy
  double *amas;    // Mass (Not Used)
  double *rvdw;    // vdW Energy
  double *eps;     // vdW Energy
  double *gam;     // Coulomb Energy/ Charge Distribution
  double *rapt;    // Determine Bond Order
  double *stlp;    // Determine Nr of Lone Pairs
  double *alf;     // vdW Energy
  double *vop;     // vdW Energy
  double *valf;    // Valency Angle Energy
  double *valp1;   // Undercoodination
  double *valp2;   // (Not Used)
  double *chi;     // Charge Distribution
  double *eta;     // Charge Distribution 
  double *vnphb;   // Hydrogen Bond
  double *vnq;     // Determine Bond Order
  double *vlp1;    // Lone Pair Energy
  double *vincr;   // (Not Used)
  double *bo131;   // Correction for Overcoodination
  double *bo132;   // Correction for Overcoodination
  double *bo133;   // Correction for Overcoodination
  double *sigqeq;  // (Not Used)
  double *def;     // (Not Used)
  double *vovun;   // Over/Undercoordination Energy
  double *vval1;   // Valency Angle Energy
  double *vrom;    // (Not Used)
  double *vval3;   // Correction for Overcoodination
  double *vval4;   // Valency Angle Energy
  double *rcore2;  // (Not Used)
  double *ecore2;  // (Not Used)
  double *acore2;  // (Not Used)

// 2-Body (Bond)
  double *de1;     // Bond Energy/ Overcoordination Energy
  double *de2;     // Bond Energy
  double *de3;     // Bond Energy
  double *psi;     // Bond Energy
  double *pdo;     // Determine Bond Order
  double *v13cor;  // Correction for Overcoodination
  double *popi;    // Determine Bond Order
  double *vover;   // Overcoordination Energy
  double *psp;     // Bond Energy
  double *pdp;     // Determine Bond Order
  double *ptp;     // Determine Bond Order
  double *bom;     // (Not Used)
  double *bop1;    // Determine Bond Order
  double *bop2;    // Determine Bond Order
  double *ovc;     // Correction for Overcoodination
  double *vuncor;  // (Not Used)
// Off-Diagonal
  double *deodmh;  // vdW Energy
  double *rodmh;   // vdW Energy
  double *godmh;   // vdW Energy
  double *rsig;    // Determine Bond Order
  double *rpi;     // Determine Bond Order
  double *rpi2;    // Determine Bond Order
// Hydrogen Bond
  double *rhb;     // Hydrogen Bond Energy
  double *dehb;    // Hydrogen Bond Energy
  double *vhb1;    // Hydrogen Bond Energy
  double *vhb2;    // Hydrogen Bond Energy

// 3-Body (Angle)
  double *th0;     // Valency Angle Energy
  double *vka;     // Valency Angle Energy
  double *vka3;    // Valency Angle Energy
  double *vka8;    // Valency Angle Conjugation Energy
  double *vkac;    // Valency Angle Energy
  double *vkap;    // Valency Angle Energy
  double *vval2;   // Valency Angle Energy

// 4-Body (Torsion)
  double *v1;      // Torsion Angle Energy
  double *v2;      // Torsion Angle Energy
  double *v3;      // Torsion Angle Energy
  double *v4;      // Torsion Angle Energy
  double *vconj;   // Conjugation Energy
  double *v2bo;    // (Not Used)
  double *v3bo;    // (Not Used)
#endif
// End of Addition

  fvalue_pointer *fvalue; /* function pointers for analytic potentials */
} apot_table_t;

#endif /* APOT */

/****************************************************************
 *
 *  potfit_calculation: holds data for force calculation
 *
 ****************************************************************/

typedef struct
{
  int fcalls;  /* number of force calculations */
  int mdim;    /* total number of entries in force vector */
  int ndim;    /* number of free optimization parameters in force vector */
  int ndimtot; /* total number of optimization parameters in force vector */
  int paircol; /* How many columns for pair potential ( ntypes*(ntypes+1)/2 ) */

  double d_eps;   /* ??? */
  double* force;  //

  struct
  {
    double* vecu_bracket;
    double* vecu_brent;
    double* f_vec3;
  } linmin;

  /* pointers for force-vector */
  int energy_p; /* offset of energies in force vector */
#if defined(STRESS)
  int stress_p; /* offset of stresses in force vector */
#endif          /* STRESS */
#if defined(EAM) || defined(ADP) || defined(MEAM)
  int dummy_p; /* pointer to dummy constraints */
  int limit_p; /* pointer to limiting constraints */
#endif         /* EAM || ADP || MEAM */
#if defined(APOT)
  int punish_par_p; /* pointer to parameter punishment contraints */
  int punish_pot_p; /* pointer to potential punishment constraints */
#endif              /* APOT */
} potfit_calculation;

/****************************************************************
 *
 *  potfit_configurations: holds reference configurations data
 *
 ****************************************************************/

typedef struct
{
  int natoms; /* total number of atoms */
  int nconf;  /* total number of configurations */

  atom_t* atoms;      /* atoms array */
  atom_t* conf_atoms; /* Atoms in configuration */

  char** elements; /* element names from configuration files */

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
#endif                   /* STRESS */

// Addition by AI 17/7/2015
#ifdef LMP
 lattice_t *lattice;
#endif /* LMP */
// End of addition
} potfit_configurations;

/****************************************************************
 *
 *  potfit_filenames: holds all kinds of filenames
 *
 ****************************************************************/

typedef struct
{
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

/****************************************************************
 *
 *  potfit_mpi_config: holds information needed for MPI calculation
 *
 ****************************************************************/

typedef struct
{
  int myid;     /* index of current process */
  int num_cpus; /* total numer of processes */

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
#endif /* MPI */
} potfit_mpi_config;

/****************************************************************
 *
 *  potfit_parameters: holds information from parameter file
 *
 ****************************************************************/

typedef struct
{
  int imdpotsteps; /* resolution of IMD potential */
  int ntypes;      /* number of atom types */
  int opt;         /* optimization flag */
  int rng_seed;    /* seed for RNG */
  int usemaxch;    /* use maximal changes file */

  int write_output_files;
  int write_lammps_files;
  int write_pair;
  int writeimd;
  int write_lammps; /* write output also in LAMMPS format */

#if defined(EVO)
  double evo_threshold;
#else
  char* anneal_temp;
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

/****************************************************************
 *
 *  potfit_potentials: holds information from potential file
 *
 ****************************************************************/

typedef struct
{
  int* gradient; /* Gradient of potential fns.  */
  int* invar_pot;
  int format;     /* format of potential table */
  int have_invar; /* Are invariant pots specified?  */
#if defined(APOT)
  int* smooth_pot;
  int cp_start;         /* cp in opt_pot.table */
  int global_idx;       /* index for global parameters in opt_pot table */
  int global_pot;       /* number of "potential" for global parameters */
  int have_globals;     /* do we have global parameters? */
  double* calc_list;    /* list of current potential in the calc table */
  double* compnodelist; /* list of the composition nodes */
#endif                  /* APOT */

  /* potential tables */
  pot_table_t opt_pot;  /* potential in the internal representation used for
  minimisation */
  pot_table_t calc_pot; /* the potential table used for force calculations */
#ifdef APOT
  apot_table_t apot_table; /* potential in analytic form */
#endif /* APOT */
// Addition by AI 17/7/2015
#ifdef LMP
void* lammpsObj;
#endif /* LMP */
// End of addition
} potfit_potentials;

/****************************************************************
 *
 *  potfit_unknown: strange parameters, need to be merge into any existing
 *struct
 *
 ****************************************************************/

typedef struct
{
  /* misc. stuff - has to belong somewhere */
  int* idx;
  int init_done;
  int plot; /* plot output flag */
#if defined EAM || defined ADP || defined MEAM
  double* lambda; /* embedding energy slope... */
#endif
  double* maxchange; /* Maximal permissible change */
/* variables needed for electrostatic options */
#ifdef COULOMB
  double dp_eps; /* ??? in eV A */
  double dp_cut; /* cutoff-radius for long-range interactions */
#endif           /* COULOMB */
#ifdef DIPOLE
  double dp_tol; /* dipole iteration precision */
  double dp_mix; /* */
#endif           /* DIPOLE */
  char* interaction_name;
} potfit_unknown;

/****************************************************************
 *
 *  ???
 *
 ****************************************************************/

typedef struct
{
  double d_eps; /* abortion criterion for powell_lsq */
} potfit_optimization;
