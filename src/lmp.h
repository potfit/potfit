#include "string.h"
/*#include "/home/kubo/src/mpi.h"
#include "/home/kubo/src/lammps.h"
#include "/home/kubo/src/input.h"
#include "/home/kubo/src/compute.h"
#include "/home/kubo/src/modify.h"
#include "/home/kubo/src/update.h"
#include "/home/kubo/src/atom.h"
#include "/home/kubo/src/force.h"
#include "/home/kubo/src/pair.h"
#include "/home/kubo/src/domain.h"
#include "/home/kubo/src/bond.h"
#include "/home/kubo/src/angle.h"
#include "/home/kubo/src/dihedral.h"
#include "/home/kubo/src/improper.h"
#include "/home/kubo/src/kspace.h"*/

//#include "lammps.h"
//#include "input.h"
//#include "string.h"
//#include "compute.h"
//#include "modify.h"
//#include "update.h"
//#include "atom.h"
//#include "force.h"
//#include "pair.h"
//#include "domain.h"
//#include "bond.h"
//#include "angle.h"
//#include "dihedral.h"
//#include "improper.h"
//#include "kspace.h"

#ifdef __cplusplus
  extern "C" {
#endif

double lammps(int, int, int, double *,double *,double *,int *,double *,double*, void*);
void open_lammps(void **);
void close_lammps(void *);
double get_energy(void *);
void get_force(void *,double *);
void get_coord(void *,double *);
void read_dummy(void *,char*);
void input_coord(void *,double *,int *,int);
void input_tilt(void *,double *);
void input_box(void *,double *);
void input_num_and_type(void *,int,int);
void lammps_input(void *,char *);
void make_dummy(int,int);
void remove_dummy();
int get_numtotal(void *);

#ifdef __cplusplus
}
#endif
