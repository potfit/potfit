#include "stdio.h"
#include "string.h"
#include "mpi.h"
#include "lmp.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "pair_reax_c.h"
#include "universe.h"
#include "modify.h"
#include "compute.h"
#include "time.h"
#include "unistd.h"

using namespace	LAMMPS_NS;

void make_dummy(int me, int numatoms,int numtypes)
{
    char str[256];
    sprintf(str, "data.dummy.%d", me);

    //if (access(str, F_OK) != -1)
//	return;

    FILE *dummy = fopen(str,"w");
    fprintf(dummy,"\nMasses\n");
    fprintf(dummy,"\n");
    for(int i=1;i<=numtypes;i++) fprintf(dummy,"%d 1.0\n",i);
    fprintf(dummy,"\n");
    fprintf(dummy,"Atoms\n");
    fprintf(dummy,"\n");
    for(int i=1;i<=numatoms;i++) fprintf(dummy,"%d 1 0 0.0 0.0 0.0\n",i);
    fprintf(dummy,"\n");
    fclose(dummy);
}

void read_dummy(void *ptr,char *str)
{
	LAMMPS *lammps = (LAMMPS *)ptr;
	lammps->domain->box_exist = 0;
	lammps->input->one(str);
}

void open_lammps(void **ptr)
{
	//MPI_Init(0,NULL);
	LAMMPS *lammps = new LAMMPS(0,NULL,MPI_COMM_SELF);
	*ptr = (LAMMPS *)lammps;
}

void lammps_input(void *ptr,char *str)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	lammps ->input->one(str);
}

void close_lammps(void *ptr)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	delete lammps;
	//MPI_Finalize();
}

void input_num_and_type(void *ptr,int numatoms,int numtypes)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	lammps->atom->natoms = numatoms;
	lammps->atom->ntypes = numtypes;
}

void input_box(void *ptr,double *boxsizehi)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	double *boxsizelo = (double *)malloc(sizeof(double) * 3);
	for(int i=0;i<3;i++)
	{
		boxsizelo[i] = 0.0;
		lammps->domain->boxlo[i] = boxsizelo[i];
		lammps->domain->boxhi[i] = boxsizehi[i];
	}
	free(boxsizelo);boxsizelo = NULL;
}

void input_tilt(void *ptr,double *tilt)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	lammps->domain->triclinic = 1;
	lammps->domain->xy = tilt[0];
	lammps->domain->xz = tilt[1];
	lammps->domain->yz = tilt[2];
}

void input_coord(void *ptr,double *coords,int *type,int numatoms)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	double **c = lammps->atom->x;
	int *types = lammps->atom->type;
	int offset;
	for(int i=0;i<numatoms;i++)
	{
		offset = 3*i;
		c[i][0] = coords[offset+0];
		c[i][1] = coords[offset+1];
		c[i][2] = coords[offset+2];
/*		c[i][0] = coords[i][0];
		c[i][1] = coords[i][1];
		c[i][2] = coords[i][2];  */
		types[i] = type[i];
	}
}

void get_coord(void *ptr,double *ccopy)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	double **x = lammps->atom->x;
	int *tag = lammps->atom->tag;
	int nlocal = lammps->atom->nlocal;
	int id,offset;
	for(int i=0;i<nlocal;i++)
	{
	    id = tag[i];
	    offset = 3*(id-1);
    	ccopy[offset+0] = x[id-1][0];
    	ccopy[offset+1] = x[id-1][1];
    	ccopy[offset+2] = x[id-1][2];
	}
}

void get_force(void *ptr,double *f_list)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	double **f = lammps->atom->f;
	int *tag = lammps->atom->tag;
	int nlocal = lammps->atom->nlocal;
	int id,offset;
	for(int i=0;i<nlocal;i++)
	{
    	id = tag[i];
    	offset = 3*(id-1);
    	f_list[offset+0] = f[id-1][0] * 4.336e-2;
    	f_list[offset+1] = f[id-1][1] * 4.336e-2;
    	f_list[offset+2] = f[id-1][2] * 4.336e-2;
/*	    f_list[id-1][0] = f[id-1][0];
	    f_list[id-1][1] = f[id-1][1];
	    f_list[id-1][2] = f[id-1][2];  */
	}

/*	FILE *fp = fopen("force.txt","w");
	for(int i=0;i<nlocal;i++)
	{
		id = tag[i];
		offset = 3*(id-1);
		fprintf(fp,"%d\t%lf\t%lf\t%lf\n",id,f_list[offset+0],f_list[offset+1],f_list[offset+2]);
	} */
}

double get_energy(void *ptr)
{
    LAMMPS *lammps = (LAMMPS *) ptr;

    int nlocal = lammps->atom->nlocal;
    int npair = nlocal;
    int nbond = nlocal;
    int ntotal = nlocal;
    if(lammps->force->newton) npair += lammps->atom->nghost;
    if(lammps->force->newton_bond) nbond += lammps->atom->nghost;
    if(lammps->force->newton) ntotal += lammps->atom->nghost;
    double *e_list = (double *)malloc(sizeof(double) * ntotal);

    for(int i=0;i<ntotal;i++) e_list[i] = 0.0;
    if(lammps->force->pair)
    {
        double *eatomp = lammps->force->pair->eatom;
        for(int i=0;i<npair;i++) e_list[i] += eatomp[i]; //pari to local
    }
    if(lammps->force->bond)
    {
        double *eatomb = lammps->force->bond->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatomb[i]; //bond to local
    }
    if(lammps->force->angle)
    {
        double *eatoma = lammps->force->angle->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatoma[i];
    }
    if(lammps->force->dihedral)
    {
        double *eatomd = lammps->force->dihedral->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatomd[i];
    }
    if(lammps->force->improper)
    {
        double *eatomi = lammps->force->improper->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatomi[i];
    }
    if(lammps->force->kspace)
    {
        double *eatomk = lammps->force->kspace->eatom;
        for(int i=0;i<nlocal;i++) e_list[i] = eatomk[i];
    }
	double sum = 0.0;
	for(int i=0;i<ntotal;i++) sum += e_list[i];

	free(e_list);e_list = NULL;
	
	return sum;// * 4.336e-2;

/*	int *tag = lammps->atom->tag;
	int id;
    FILE *fp = fopen("energy.txt","w");
    for(int i=0;i<nlocal;i++)
     {
        id = tag[i];
        fprintf(fp,"%d\t%lf\n",id,e_list[id]);
      } */
}

int get_numtotal(void *ptr)
{
	LAMMPS *lammps = (LAMMPS *) ptr;
	int numlocal = lammps->atom->nlocal;
	int numtotal = numlocal;
	if(lammps->force->newton) numtotal += lammps->atom->nghost;
	return numtotal;
}

void remove_dummy(int me)
{
    char str[256];
    sprintf(str, "data.dummy.%d", me);
    //remove(str);
}

/* ----------------------------------------------------------------------
 *    main program to drive LAMMPS
 *    units:distance = Angstroms,energy = eV,force = eV/Angstrome
 *    ------------------------------------------------------------------------- */

double lammps(int me, int numatoms,int numtypes,double *boxsize,double *tilt,double *coords,int *type,double *f_list,double *param, void* ptr)
{
	LAMMPS* lammps = (LAMMPS* )ptr;
	//return 0;
	//printf("OK!!!\n");
	double kcal_to_ev = 4.3393122e-2;
        clock_t t1, t2; float diff;
	char str[256];
	//sprintf(str,"read_data data.dummy.%d", me);
	char input1[][256]={"units real","atom_style charge"};
	char pair_style[256] = "pair_style reax/c NULL";
	//char variables[16][256] = {"compute reax all pair reax/c", "variable eb equal c_reax[1]", "variable ea equal c_reax[2]","variable elp equal c_reax[3]","variable emol equal c_reax[4]","variable ev equal c_reax[5]","variable epen equal c_reax[6]","variable ecoa equal c_reax[7]","variable ehb equal c_reax[8]","variable et equal c_reax[9]","variable eco equal c_reax[10]","variable ew equal c_reax[11]","variable ep equal c_reax[12]","variable efi equal c_reax[13]","variable eqeq equal c_reax[14]","thermo_style custom step temp epair v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq"};
	//char pair_style[256] = "pair_style reax 10 0 1 1.0e-10";
	
	char input2[][256]={"compute pe all pe/atom", "neighbor 2.0 bin","neigh_modify every 10 delay 0 check no","fix 1 all nve", "fix 2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c","thermo_style custom etotal","timestep 1.0"};
	//char input2[][256]={"compute pe all pe/atom", "neighbor 2.0 bin","neigh_modify every 10 delay 0 check no","fix 1 all nve", "fix 2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c","timestep 1.0"};
	
	//char input2[][256]={"pair_coeff * * ffield.reax 1 2 3","compute pe all pe/atom", "neighbor 2.0 bin","neigh_modify every 10 delay 0 check no","fix 1 all nve", "thermo_style custom etotal","timestep 1.0"};	
        char run[256] = {"run 0"};

    char pair_coeff[256],numt[256],space[256];
    int narg = 3 + numtypes;

    char **ch = (char**)malloc(sizeof(char*) * narg);

    for(int i=0;i<narg;i++) ch[i] = (char*)malloc(sizeof(char) * 256);

// Kubo 20130222 ---------------
//    ch[0] = "*";ch[1] = "*";ch[2] = "ffield.reax";
    strcpy(ch[0],"*");
    strcpy(ch[1],"*");
    strcpy(ch[2],"ffield.reax");
// End -------------------------

    for(int i=3;i<narg;i++)
    {
        sprintf(numt,"%d",i-2);
        strcpy(ch[i],numt);
    }

	//make_dummy(me,numatoms,numtypes);

	// input units and atom_style
	for(int i=0;i<2;i++) 
	lammps->input->one(input1[i]);
	
	// input number of atoms and atom type
        lammps->atom->natoms = numatoms;
        lammps->atom->ntypes = numtypes;

	// input box size
        for(int i=0;i<3;i++)
        {
                lammps->domain->boxlo[i] = 0.0;
                lammps->domain->boxhi[i] = boxsize[i];
        }
	
	// input tilt angle if necessary
	if(tilt[0] != 0 || tilt[1] != 0 || tilt[2] != 0)
	{
	        lammps->domain->triclinic = 1;
	      	lammps->domain->xy = tilt[0];
	        lammps->domain->xz = tilt[1];
	        lammps->domain->yz = tilt[2];	
	}
	
	lammps->input->one("processors 1 1 1 map xyz");
        lammps->domain->box_exist = 1;
        //lammps->input->one(str);
        //lammps->input->one("lattice sc 0.8");
	for (int i = 0; i < numatoms; i++)
	{
		lammps->input->one("create_atoms 1 single 0 0 0 units box");
	}
	//lammps->input->one(str);
	
	
        double **c = lammps->atom->x;
        int *types = lammps->atom->type;
        int offset;
        for(int i=0;i<numatoms;i++)
        {
                offset = 3*i;
                c[i][0] = coords[offset+0];
                c[i][1] = coords[offset+1];
                c[i][2] = coords[offset+2];
                types[i] = type[i];
        }
	
	// input pair_style and pair_coeff
        lammps->input->one(pair_style);
	lammps->force->pair->coeff2(narg,ch,param);

	// input other commands
	for(int i=0; i<7;i++) 
	        lammps ->input->one(input2[i]);

	//for(int i=0; i<16;i++)
	//	lammps->input->one(variables[i]);

	// run
        lammps ->input->one(run);
  
	// get force
	for(int i=0;i<3*numatoms;i++) f_list[i] = 0.0;
        double **f = lammps->atom->f;
        int *tag = lammps->atom->tag;
        int nlocal = lammps->atom->nlocal;
        int id;
        for(int i=0;i<nlocal;i++)
        {
	        id = tag[i];
	        offset = 3*(id-1);
	        f_list[offset+0] = f[i][0] * kcal_to_ev;//4.336e-2;
	        f_list[offset+1] = f[i][1] * kcal_to_ev;//4.336e-2;
	        f_list[offset+2] = f[i][2] * kcal_to_ev;//4.336e-2;
        }

	// get energy
    int npair = nlocal;
    int nbond = nlocal;
    int ntotal = nlocal;
    if(lammps->force->newton) npair += lammps->atom->nghost;
    if(lammps->force->newton_bond) nbond += lammps->atom->nghost;
    if(lammps->force->newton) ntotal += lammps->atom->nghost;
    double *e_list = (double *)malloc(sizeof(double) * ntotal);

    for(int i=0;i<ntotal;i++) e_list[i] = 0.0;
    if(lammps->force->pair)
    {
        double *eatomp = lammps->force->pair->eatom;
        for(int i=0;i<npair;i++) e_list[i] = e_list[i] += eatomp[i];
    }
    if(lammps->force->bond)
    {
        double *eatomb = lammps->force->bond->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatomb[i];
    }
    if(lammps->force->angle)
    {
        double *eatoma = lammps->force->angle->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatoma[i];
    }
    if(lammps->force->dihedral)
    {
        double *eatomd = lammps->force->dihedral->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatomd[i];
    }
    if(lammps->force->improper)
    {
        double *eatomi = lammps->force->improper->eatom;
        for(int i=0;i<nbond;i++) e_list[i] += eatomi[i];
    }
    if(lammps->force->kspace)
    {
        double *eatomk = lammps->force->kspace->eatom;
        for(int i=0;i<nlocal;i++) e_list[i] = eatomk[i];
    }
        double sum = 0.0;
        for(int i=0;i<ntotal;i++) sum += e_list[i];
	
        free(e_list);e_list = NULL; 

 	lammps->input->one("clear");

        for(int i=0;i<narg;i++){free(ch[i]);ch[i] = NULL;} // Kubo 20130222
	
	free(ch);ch = NULL;
	
	return sum * kcal_to_ev;//4.336e-2;
}
