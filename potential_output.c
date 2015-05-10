/****************************************************************
 *
 * potential_output.c: Routines for writing the potential table
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

#include "potfit.h"

#include "elements.h"
#include "functions.h"
#include "potential.h"
#include "splines.h"
#include "utils.h"

#define NPLOT 1000

#ifdef COULOMB

/****************************************************************
 *
 * write coulomb-potential
 *
 ****************************************************************/

void write_coulomb_table()
{
  apot_table_t *apt = &apot_table;
  if (ntypes == 2) {
    int   i, j;
    double value, c1;
    FILE *outfile;
    char *filename1 = "Coulomb_00";
    char *filename2 = "Coulomb_01";
    char *filename3 = "Coulomb_11";

    c1 = -apt->charge[0] / 2;

    outfile = fopen(filename1, "a");
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < atoms[i].num_neigh; j++) {
	value = apt->charge[0] * apt->charge[0] * atoms[i].neigh[j].fnval_el;
	fprintf(outfile, "%f\t%f\n", atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);

    outfile = fopen(filename2, "a");
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < atoms[i].num_neigh; j++) {
	value = apt->charge[0] * c1 * atoms[i].neigh[j].fnval_el;
	fprintf(outfile, "%f\t%f\n", atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);

    outfile = fopen(filename3, "a");
    for (i = 0; i < natoms; i++) {
      for (j = 0; j < atoms[i].num_neigh; j++) {
	value = c1 * c1 * atoms[i].neigh[j].fnval_el;
	fprintf(outfile, "%f\t%f\n", atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);
  } else {
    printf("Coulomb-outfiles are only available in case of two atom types.");
  }

  return;
}

#endif /* COULOMB */

#ifdef APOT

/****************************************************************
 *
 *  write potential table (format 0)
 *
 ****************************************************************/

void write_pot_table0(apot_table_t *apt, char *filename)
{
  int   i, j;
  FILE *outfile;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 0 %d", apt->number);
  fprintf(outfile, "\n#T %s", interaction_name);

  if (have_elements) {
    /* write elements */
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);

    /* write order of the interactions */
    fprintf(outfile, "\n##");
    /* pair potentials */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);

#if defined EAM || defined ADP || defined MEAM
    /* transfer functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif /* EAM || ADP || MEAM */

#ifdef ADP
    /* dipole terms */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
    /* quadrupole terms */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#endif /* ADP */

#ifdef MEAM
    /* f terms */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
    /* g terms */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif /* MEAM */

#ifdef STIWEB
    /* stiweb_3 terms */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
    /* lambda terms */
    int   k = 0;
    for (i = 0; i < ntypes; i++)
      for (j = 0; j < ntypes; j++)
	for (k = j; k < ntypes; k++)
	  fprintf(outfile, " %s-%s-%s", elements[i], elements[j], elements[k]);
#endif /* STIWEB */

#ifdef TERSOFF
    /* mixing terms */
    for (i = 0; i < ntypes; i++)
      for (j = i + 1; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#endif /* TERSOFF */

  }

  /* write invariant switch for individual potentials */
  if (have_invar) {
    fprintf(outfile, "\n#I");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }

  /* end tag */
  fprintf(outfile, "\n#E\n\n");

#ifdef PAIR
  /* write chemical potentials if enabled */
  if (enable_cp) {
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, "cp_%s %.10f %.2f %.2f\n", elements[i], apt->chempot[i],
	apt->pmin[apt->number][i], apt->pmax[apt->number][i]);
#ifdef CN
    if (compnodes > 0)
      fprintf(outfile, "cn %d\n", compnodes);
    for (j = 0; j < compnodes; j++)
      fprintf(outfile, "%.2f %.10f %.2f %.2f\n", compnodelist[j],
	apt->chempot[ntypes + j], apt->pmin[apt->number][ntypes + j], apt->pmax[apt->number][ntypes + j]);
#endif /* CN */
    fprintf(outfile, "\n");
  }
#endif /* PAIR */

#ifdef COULOMB
  fprintf(outfile, "elstat\n");
  for (i = 0; i < ntypes - 1; i++)
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number][i],
      apt->charge[i], apt->pmin[apt->number][i], apt->pmax[apt->number][i]);
  fprintf(outfile, "charge_%s\t %f\n", elements[ntypes - 1], apt->last_charge);
  fprintf(outfile, "%s\t\t %f\t %f\t %f\n", apt->param_name[apt->number + 1][0],
    apt->dp_kappa[0], apt->pmin[apt->number + 1][0], apt->pmax[apt->number + 1][0]);

#ifdef DIPOLE
  for (i = 0; i < ntypes; i++)
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 2][i],
      apt->dp_alpha[i], apt->pmin[apt->number + 2][i], apt->pmax[apt->number + 2][i]);
  for (i = 0; i < paircol; i++) {
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 3][i],
      apt->dp_b[i], apt->pmin[apt->number + 3][i], apt->pmax[apt->number + 3][i]);
  }
  for (i = 0; i < paircol; i++) {
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 4][i],
      apt->dp_c[i], apt->pmin[apt->number + 4][i], apt->pmax[apt->number + 4][i]);
  }
#endif /* DIPOLE */

  fprintf(outfile, "\n");
#endif /* COULOMB */

  /* write globals section */
  if (have_globals) {
    fprintf(outfile, "global\t%d\n", apt->globals);
    for (i = 0; i < apt->globals; i++)
      fprintf(outfile, "%s\t%18.8f\t%12.4f\t%12.4f\n",
	apt->param_name[global_pot][i], apt->values[global_pot][i], apt->pmin[global_pot][i],
	apt->pmax[global_pot][i]);
    fprintf(outfile, "\n");
  }

  /* write data */
  for (i = 0; i < apt->number; i++) {
    if (smooth_pot[i]) {
      fprintf(outfile, "type\t%s_sc\n", apt->names[i]);
    } else {
      fprintf(outfile, "type\t%s\n", apt->names[i]);
    }
    fprintf(outfile, "cutoff\t%f\n", apot_table.end[i]);
    fprintf(outfile, "# rmin\t%f\n", apt->begin[i]);
    for (j = 0; j < apt->n_par[i]; j++) {
      if (apt->param_name[i][j][strlen(apt->param_name[i][j]) - 1] != '!') {
	fprintf(outfile, "%s\t%18.8f\t%12.4f\t%12.4f\n", apt->param_name[i][j],
	  apt->values[i][j], apt->pmin[i][j], apt->pmax[i][j]);
      } else {
	fprintf(outfile, "%s\n", apt->param_name[i][j]);
      }
    }
    if (i != (apt->number - 1))
      fprintf(outfile, "\n");
  }
  fclose(outfile);

  return;
}

#else /* APOT */

/****************************************************************
 *
 *  write potential table (format 3)
 *
 ****************************************************************/

void write_pot_table3(pot_table_t *pt, char *filename)
{
  FILE *outfile = NULL, *outfile2 = NULL;
  int   i, j, flag = 0;
  double r;

  if (*plotpointfile != '\0')
    flag = 1;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* if needed: open file for plotpoints */
  if (flag) {
    outfile2 = fopen(plotpointfile, "w");
    if (NULL == outfile)
      error(1, "Could not open file %s\n", filename);
  }

  /* write header */
  fprintf(outfile, "#F 3 %d", pt->ncols);
  fprintf(outfile, "\n#T %s", interaction_name);
  if (have_elements) {
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    fprintf(outfile, "\n##");
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#if defined EAM || defined MEAM
    /* transfer functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#ifdef TBEAM
    /* transfer functions s-band */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions s-band */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif /* TBEAM */
#endif /* EAM || MEAM */
#ifdef MEAM
    /* pre-anglpart */
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
    /* angl part */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif /* MEAM */
  }
  if (have_invar) {
    fprintf(outfile, "\n#I");
    for (i = 0; i < pt->ncols; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }
  fprintf(outfile, "\n#G");
  for (i = 0; i < pt->ncols; i++)
    fprintf(outfile, " %d", gradient[i]);
  fprintf(outfile, "\n#E\n");

  /* write info block */
  for (i = 0; i < pt->ncols; i++) {
    fprintf(outfile, "%.16e %.16e %d\n", pt->begin[i], pt->end[i], pt->last[i] - pt->first[i] + 1);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < pt->ncols; i++) {
    r = pt->begin[i];
    /* write gradient */
    fprintf(outfile, "%.16e %.16e\n", pt->table[pt->first[i] - 2], pt->table[pt->first[i] - 1]);
    for (j = pt->first[i]; j <= pt->last[i]; j++) {
      fprintf(outfile, "%.16e\n", pt->table[j]);
      if (flag)
	fprintf(outfile2, "%.6e %.6e %d\n", r, pt->table[j], j);
      r += pt->step[i];
    }
    fprintf(outfile, "\n");
    if (flag)
      fprintf(outfile2, "\n\n");
  }
  fclose(outfile);
  if (flag)
    fclose(outfile2);
}

/****************************************************************
 *
 *  write potential table (format 4)
 *
 ****************************************************************/

void write_pot_table4(pot_table_t *pt, char *filename)
{
  FILE *outfile = NULL, *outfile2 = NULL;
  int   i, j, flag = 0;

  if (*plotpointfile != '\0')
    flag = 1;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* if needed: open file for plotpoints */
  if (flag) {
    outfile2 = fopen(plotpointfile, "w");
    if (NULL == outfile)
      error(1, "Could not open file %s\n", filename);
  }

  /* write header */
  fprintf(outfile, "#F 4 %d", pt->ncols);
  fprintf(outfile, "\n#T %s", interaction_name);
  if (have_elements) {
    fprintf(outfile, "\n#C");
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    fprintf(outfile, "\n##");
    for (i = 0; i < ntypes; i++)
      for (j = i; j < ntypes; j++)
	fprintf(outfile, " %s-%s", elements[i], elements[j]);
#ifdef EAM
    /* transfer functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#ifdef TBEAM
    /* transfer functions s-band */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
    /* embedding functions s-band */
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %s", elements[i]);
#endif /* TBEAM */
#endif /* EAM */
  }
  if (have_invar) {
    fprintf(outfile, "\n#I");
    for (i = 0; i < pt->ncols; i++)
      fprintf(outfile, " %d", invar_pot[i]);
  }
  fprintf(outfile, "\n#G");
  for (i = 0; i < pt->ncols; i++)
    fprintf(outfile, " %d", gradient[i]);
  fprintf(outfile, "\n#E\n");

  /* write info block */
  for (i = 0; i < pt->ncols; i++) {
    fprintf(outfile, "%d\n", pt->last[i] - pt->first[i] + 1);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < pt->ncols; i++) {
    fprintf(outfile, "%.16e %.16e\n", pt->table[pt->first[i] - 2], pt->table[pt->first[i] - 1]);
    for (j = pt->first[i]; j <= pt->last[i]; j++) {
      fprintf(outfile, "%.16e %.16e\n", pt->xcoord[j], pt->table[j]);
      if (flag)
	fprintf(outfile2, "%.6e %.6e %d\n", pt->xcoord[j], pt->table[j], j);
    }
    fprintf(outfile, "\n");
    if (flag)
      fprintf(outfile2, "\n\n");
  }
  fclose(outfile);
  if (flag)
    fclose(outfile2);
}

#endif /* APOT */

/****************************************************************
 *
 *  write potential table for IMD (format 2)
 *
 ****************************************************************/

#if defined STIWEB || defined TERSOFF

void write_pot_table_imd(pot_table_t *pt, char *prefix)
{
  int   i;
  char  filename[255];
  FILE *outfile;

  /* access pt to shut up compiler warnings */
  i = pt->len;

#ifdef STIWEB
  sprintf(filename, "%s.sw.pot", prefix);
#endif /* STIWEB */
#ifdef TERSOFF
#ifndef TERSOFFMOD
  sprintf(filename, "%s.tersoff.pot", prefix);
#else
  sprintf(filename, "%s.tersoffmod.pot", prefix);
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */
  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  fprintf(outfile, "# IMD potential file written by %s\n#\n", POTFIT_VERSION);
#ifdef STIWEB
  fprintf(outfile, "# This is a Stillinger-Weber potential. Compile IMD with the 'stiweb' option.\n#\n");
#endif /* STIWEB */
#ifdef TERSOFF
#ifndef TERSOFFMOD
  fprintf(outfile, "# This is a Tersoff potential. Compile IMD with the 'tersoff2' option.\n#\n");
#else
  fprintf(outfile, "# This is a modified Tersoff potential. Compile IMD with the 'tersoffmod2' option.\n#\n");
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */

  /* write warning header */
  fprintf(outfile, "# WARNING:\n");
  fprintf(outfile, "# DO NOT USE THIS FILE AS A POTENTIAL FILE IN IMD !!!!!\n");
  fprintf(outfile, "# COPY THE CONTENTS OF THIS FILE INTO THE IMD PARAMETER FILE\n\n");

#ifdef STIWEB
  write_imd_data_pair(outfile, "stiweb_a\t", 0, 0);	/* A_ij */
  write_imd_data_pair(outfile, "stiweb_b\t", 0, 1);	/* B_ij */
  write_imd_data_pair(outfile, "stiweb_p\t", 0, 2);	/* p_ij */
  write_imd_data_pair(outfile, "stiweb_q\t", 0, 3);	/* q_ij */
  write_imd_data_pair(outfile, "stiweb_de\t", 0, 4);	/* delta_ij */
  write_imd_data_pair(outfile, "stiweb_a1\t", 0, 5);	/* a1_ij */
  write_imd_data_pair(outfile, "stiweb_ga\t", paircol, 0);	/* gamma_ij */
  write_imd_data_pair(outfile, "stiweb_a2\t", paircol, 1);	/* a2_ij */

  /* lambda_ijk */
  /* strange format (e.g. binary): 000, 100, 001, 101, 011, 111 */
  int   j;
  fprintf(outfile, "stiweb_la\t");
  for (j = 0; j < paircol; j++)
    for (i = 0; i < ntypes; i++)
      fprintf(outfile, " %g", apot_table.values[apot_table.number - 1][i * paircol + j]);
  fprintf(outfile, "\n");
#endif /* STIWEB */

#ifdef TERSOFF
#ifndef TERSOFFMOD
  write_imd_data_pair(outfile, "ters_a\t\t", 0, 0);	/* A_ij */
  write_imd_data_pair(outfile, "ters_b\t\t", 0, 1);	/* B_ij */
  write_imd_data_pair(outfile, "ters_la\t\t", 0, 2);	/* lambda_ij */
  write_imd_data_pair(outfile, "ters_mu\t\t", 0, 3);	/* mu_ij */
  write_imd_data_pair(outfile, "ters_ga\t\t", 0, 4);	/* gamma_ij */
  write_imd_data_pair(outfile, "ters_n\t\t", 0, 5);	/* n_ij */
  write_imd_data_pair(outfile, "ters_c\t\t", 0, 6);	/* c_ij */
  write_imd_data_pair(outfile, "ters_d\t\t", 0, 7);	/* d_ij */
  write_imd_data_pair(outfile, "ters_h\t\t", 0, 8);	/* h_ij */
  write_imd_data_pair(outfile, "ters_r_cut\t\t", 0, 9);	/* r_cut = S_ij */
  write_imd_data_pair(outfile, "ters_r0\t\t", 0, 10);	/* r0 = R_ij */

  /* chi and omega only for mixing potentials */
  if (ntypes > 1) {
    /* chi_ij, chi_ii = 1.0 */
    fprintf(outfile, "ters_chi\t\t");
    for (i = 0; i < ntypes * (ntypes - 1) / 2.0; i++)
      fprintf(outfile, " %g", apot_table.values[paircol + i][0]);
    fprintf(outfile, "\n");

    /* omega_ij, chi_ii = 1.0 */
    fprintf(outfile, "ters_om\t\t");
    for (i = 0; i < ntypes * (ntypes - 1) / 2.0; i++)
      fprintf(outfile, " %g", apot_table.values[paircol + i][1]);
    fprintf(outfile, "\n");
  }
#else
  write_imd_data_pair(outfile, "ters_a\t\t", 0, 0);	/* A_ij */
  write_imd_data_pair(outfile, "ters_b\t\t", 0, 1);	/* B_ij */
  write_imd_data_pair(outfile, "ters_la\t\t", 0, 2);	/* lambda_ij */
  write_imd_data_pair(outfile, "ters_mu\t\t", 0, 3);	/* mu_ij */
  write_imd_data_pair(outfile, "ters_eta\t\t", 0, 4);	/* eta_ij */
  write_imd_data_pair(outfile, "ters_delta\t\t", 0, 5);	/* delta_ij */
  write_imd_data_pair(outfile, "ters_alpha\t\t", 0, 6);	/* alpha_ij */
  write_imd_data_pair(outfile, "ters_beta\t\t", 0, 7);	/* beta_ij */
  write_imd_data_pair(outfile, "ters_c1\t\t", 0, 8);	/* c1_ij */
  write_imd_data_pair(outfile, "ters_c2\t\t", 0, 9);	/* c2_ij */
  write_imd_data_pair(outfile, "ters_c3\t\t", 0, 10);	/* c3_ij */
  write_imd_data_pair(outfile, "ters_c4\t\t", 0, 11);	/* c4_ij */
  write_imd_data_pair(outfile, "ters_c5\t\t", 0, 12);	/* c5_ij */
  write_imd_data_pair(outfile, "ters_h\t\t", 0, 13);	/* h_ij */
  write_imd_data_pair(outfile, "ters_r0\t\t", 0, 14);	/* R1_ij */
  write_imd_data_pair(outfile, "ters_r_cut\t\t", 0, 15);	/* R2_ij */
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */

  fclose(outfile);
  printf("Parameters for IMD potential written to\t%s\n", filename);
}

#else

void write_pot_table_imd(pot_table_t *pt, char *prefix)
{
  int   i, j, k, m, m2, col1, col2;
  double r2;
#if defined(APOT)
  double temp = 0.0;
#if defined (EAM) || defined(ADP) || defined(MEAM)
  double root = 0.0;
#endif
#endif
#if !defined(APOT) && (defined (EAM) || defined(ADP) || defined(MEAM))
  double temp = 0.0;
  double temp2;
  double root;
#endif /* !APOT */
  double *r2begin, *r2end, *r2step;
  FILE *outfile;
  char  filename[255];

  /* allocate memory */
  r2begin = (double *)malloc(ntypes * ntypes * sizeof(double));
  r2end = (double *)malloc(ntypes * ntypes * sizeof(double));
  r2step = (double *)malloc(ntypes * ntypes * sizeof(double));
  if ((r2begin == NULL) || (r2end == NULL) || (r2step == NULL))
    error(1, "Cannot allocate memory in  write_pot_table_imd");

  /* pair potential part (over r^2) */
  sprintf(filename, "%s_phi.imd.pt", prefix);
  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = dsquare((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = dsquare(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif /* APOT */
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2], r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
#ifdef MEAM
	fprintf(outfile, "%.16e\n", splint_ne_lin(pt, pt->table, col1, sqrt(r2)));
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* MEAM */
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD pair potential written to \t\t%s\n", filename);

#if defined EAM || defined ADP || defined MEAM
  /* write transfer function (over r^2) */
  sprintf(filename, "%s_rho.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      col1 = (ntypes * (ntypes + 1)) / 2 + j;
      col2 = i * ntypes + j;
#ifdef APOT
      r2begin[col2] = dsquare((plotmin == 0 ? 0.1 : plotmin));
#else
      /* Extrapolation possible  */
      r2begin[col2] = dsquare(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif /* APOT */
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2], r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      col1 = (ntypes * (ntypes + 1)) / 2 + j;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne_lin(pt, pt->table, col1, sqrt(r2)));
#endif
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD transfer function written to \t%s\n", filename);

  /* write embedding function (over r) */
  sprintf(filename, "%s_F.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    col1 = (ntypes * (ntypes + 3)) / 2 + i;
#ifdef APOT
    r2begin[i] = 0;
    r2end[i] = pt->end[col1];
#else
    /* pad with zeroes */
    r2begin[i] = pt->begin[col1] - extend * pt->step[col1];
    /* extrapolation */
    r2end[i] = pt->end[col1] + extend * pt->step[col1];
#endif /* APOT */
    r2step[i] = (r2end[i] - r2begin[i]) / imdpotsteps;
    fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[i], r2end[i], r2step[i]);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < ntypes; i++) {
    r2 = r2begin[i];
    col1 = (ntypes * (ntypes + 3)) / 2 + i;
    root = (pt->begin[col1] > 0.0) ? pt->table[pt->first[col1]] / sqrt(pt->begin[col1]) : 0.0;
    root += (pt->end[col1] < 0.0) ? pt->table[pt->last[col1]] / sqrt(-pt->end[col1]) : 0.0;
    for (k = 0; k <= imdpotsteps; k++) {
#ifdef APOT
      apot_table.fvalue[col1] (r2, apot_table.values[col1], &temp);
#else
      temp = splint_ne_lin(pt, pt->table, col1, r2);
      temp2 = r2 - pt->end[col1];
      temp += (temp2 > 0.0) ? 5e2 * (temp2 * temp2 * temp2) : 0.0;
#endif /* APOT */
      fprintf(outfile, "%.16e\n", temp);
      r2 += r2step[i];
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  printf("IMD embedding function written to \t%s\n", filename);
#endif /* EAM || ADP */

#ifdef ADP
  /* write dipole function (over r^2) */
  sprintf(filename, "%s_upot.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = dsquare((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = dsquare(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif /* APOT */
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2], r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD dipole potential written to \t%s\n", filename);

  /* write quadrupole function (over r^2) */
  sprintf(filename, "%s_wpot.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += 2 * paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = dsquare((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = dsquare(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif /* APOT */
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2], r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col1 += 2 * paircol + 2 * ntypes;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD quadrupole potential written to \t%s\n", filename);
#endif /* APOT */

#ifdef MEAM
  /* write f_r2 for MEAM */
  sprintf(filename, "%s_f_meam.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes * ntypes);

  /* write info block */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = (ntypes * (ntypes + 5)) / 2;
      col1 += i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      /* Extrapolation possible  */
#ifdef APOT
      r2begin[col2] = dsquare((plotmin == 0 ? 0.1 : plotmin));
#else
      r2begin[col2] = dsquare(MAX(pt->begin[col1] - extend * pt->step[col1], 0));
#endif /* APOT */
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2], r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (i = 0; i < ntypes; i++) {
    m += i;
    m2 = 0;
    for (j = 0; j < ntypes; j++) {
      m2 += j;
      col1 = (ntypes * (ntypes + 5)) / 2;
      col1 += i < j ? i * ntypes + j - m : j * ntypes + i - m2;
      col2 = i * ntypes + j;
      r2 = r2begin[col2];
      for (k = 0; k < imdpotsteps; k++) {
#ifdef APOT
	apot_table.fvalue[col1] (sqrt(r2), apot_table.values[col1], &temp);
	temp =
	  smooth_pot[col1] ? temp * cutoff(sqrt(r2), apot_table.end[col1],
	  apot_table.values[col1][apot_table.n_par[col1] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne_lin(pt, pt->table, col1, sqrt(r2)));
#endif /* APOT */
	r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }

  fclose(outfile);
  printf("IMD MEAM f potential data written to\t%s\n", filename);

  /* write g(cos) for MEAM */
  sprintf(filename, "%s_g_meam.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", ntypes);

  /* write info block */
  for (i = 0; i < ntypes; i++) {
    col1 = (ntypes * (ntypes + 3)) + i;
    /* from -1 to +1 */
    r2begin[i] = pt->begin[col1];
    r2end[i] = pt->end[col1];
    r2step[i] = (r2end[i] - r2begin[i]) / imdpotsteps;
    fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[i], r2end[i], r2step[i]);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (i = 0; i < ntypes; i++) {
    r2 = r2begin[i];
    col1 = (ntypes * (ntypes + 3)) + i;
    for (k = 0; k <= imdpotsteps; k++) {
      fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, r2));
      r2 += r2step[i];
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  printf("IMD MEAM g potential data written to\t%s\n", filename);
#endif /* MEAM */

  free(r2begin);
  free(r2end);
  free(r2step);

  /* write endpot for IMD with electrostatics */
#if defined COULOMB && defined APOT
  apot_table_t *apt = &apot_table;
  int   ncols;
  ncols = ntypes * (ntypes + 1) / 2;
  sprintf(filename, "%s_charges.imd", prefix);

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  fprintf(outfile, "charge\t\t");
  for (i = 0; i < ntypes - 1; i++)
    fprintf(outfile, "%f\t", apt->charge[i]);
  fprintf(outfile, "%f\n", apt->last_charge);

  if ((strcmp(apt->names[0], "ms") == 0) && (strcmp(apt->names[1], "ms") == 0)
    && (strcmp(apt->names[2], "ms") == 0)) {
    fprintf(outfile, "ms_D\t\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][0]);
    fprintf(outfile, "\nms_gamma\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][1]);
    fprintf(outfile, "\nms_r0\t\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][2]);
    fprintf(outfile, "\n");
  } else if ((strcmp(apt->names[0], "buck") == 0)
    && (strcmp(apt->names[1], "buck") == 0)
    && (strcmp(apt->names[2], "buck") == 0)) {
    fprintf(outfile, "buck_a\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][0]);
    fprintf(outfile, "\nbuck_sigma\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][1]);
    fprintf(outfile, "\nbuck_c\t");
    for (i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][2]);
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "\new_rcut\t\t%f\n", dp_cut);
  fprintf(outfile, "ew_kappa\t\t%f\n", apt->dp_kappa[0]);
  fprintf(outfile, "r_cut\t\t");
  for (i = 0; i < ncols; i++)
    fprintf(outfile, "%f\t", apot_table.end[0]);
  fprintf(outfile, "\n\n");

#ifdef DIPOLE
  fprintf(outfile, "dp_alpha\t");
  for (i = 0; i < ntypes; i++)
    fprintf(outfile, "%f\t", apt->dp_alpha[i]);
  fprintf(outfile, "\ndp_b\t\t");
  for (i = 0; i < ncols; i++)
    fprintf(outfile, "%f\t", apt->dp_b[i]);
  fprintf(outfile, "\ndp_c\t\t");
  for (i = 0; i < ncols; i++)
    fprintf(outfile, "%f\t", apt->dp_c[i]);
  fprintf(outfile, "\n");
#endif /* DIPOLE */

  fclose(outfile);
  printf("Electrostatic table for IMD written to \t%s\n", filename);
#endif /* COULOMB && APOT */
}

#endif /* STIWEB */

/****************************************************************
 *
 *  write plot version of potential table
 *
 ****************************************************************/

void write_plotpot_pair(pot_table_t *pt, char *filename)
{
  FILE *outfile;
  int   i, j;
#ifndef APOT
  int   k = 0, l;
#else
  double h = 0.0;
  double temp = 0.0;
#endif /* APOT */
  double r, r_step;

  /* access pt to shut up compiler warnings */
  i = pt->len;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write data */
#ifndef APOT
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = pt->begin[k];
      r_step = (pt->end[k] - r) / (NPLOT - 1);
      for (l = 0; l < NPLOT - 1; l++) {
	fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r));
	r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
#if defined EAM || defined ADP || defined MEAM
  for (i = paircol; i < paircol + ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = paircol + ntypes; i < paircol + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
#endif
#ifdef MEAM
  for (i = paircol; i < paircol + ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = paircol + ntypes; i < paircol + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
  for (i = paircol + 2 * ntypes; i < 2 * paircol + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = 2 * paircol + 2 * ntypes; i < 2 * paircol + 3 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
#endif /* MEAM */
#ifdef ADP
  for (i = paircol + 2 * ntypes; i < 2 * (paircol + ntypes); i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - r) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
    k++;
  }
  for (i = 2 * (paircol + ntypes); i < 3 * paircol + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - r) / (NPLOT - 1);
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
    k++;
  }
#endif /* ADP */
#else /* APOT */
  for (i = 0; i < apot_table.number; i++) {
    r = (plotmin == 0 ? 0.1 : plotmin);
    r_step = (apot_table.end[i] - r) / (NPLOT - 1);
    h = apot_table.values[i][apot_table.n_par[i] - 1];
    for (j = 0; j < NPLOT; j++) {
      apot_table.fvalue[i] (r, apot_table.values[i], &temp);
      temp = smooth_pot[i] ? temp * cutoff(r, apot_table.end[i], h) : temp;
      if (isnan(temp))
	temp = 10e30;
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    if (i != (apot_table.number - 1))
      fprintf(outfile, "\n\n");
  }
#endif /* APOT */

  fclose(outfile);
  printf("Potential plotting data written to \t%s\n", filename);
}

/****************************************************************
 *
 *  write potential table for LAMMPS
 *
 ****************************************************************/

#if defined STIWEB || defined TERSOFF

void write_pot_table_lammps(pot_table_t *pt)
{
  FILE *outfile;
  char  filename[255];
  int   i, j, k;
  int   col_j;
#ifdef STIWEB
  double epsilon, sigma, scale_p, scale_q;
#endif /* STIWEB */

  if (!have_elements) {
    warning("There are no elements listed in you configuration file.\n");
    warning("A LAMMPS potential cannot be written without them.\n");
    return;
  }

  /* access pt to shut up compiler warnings */
  i = pt->len;

#ifdef STIWEB
  /* check if final potential is LAMMPS compliant (a1==a2) */
  if (apot_table.values[0][5] != apot_table.values[paircol][1]) {
    warning("Your potential is not supported by LAMMPS.\n");
    warning("Please ensure that the values a1 and a2 are the same,\n");
    warning("  either by using global parameters or by fixing them.\n");

    return;
  }
#endif /* STIWEB */

  /* open file */
  if (strcmp(output_prefix, "") != 0)
    strcpy(filename, output_prefix);
  else
    strcpy(filename, endpot);
#ifdef STIWEB
  sprintf(filename, "%s.lammps.sw", filename);
#endif /* STIWEB */
#ifdef TERSOFF
#ifndef TERSOFFMOD
  sprintf(filename, "%s.lammps.tersoff", filename);
#else
  sprintf(filename, "%s.lammps.tersoffmod", filename);
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* initialize periodic table */
  init_elements();
#ifdef STIWEB
  update_stiweb_pointers(opt_pot.table);
#endif /* STIWEB */
#ifdef TERSOFF
  update_tersoff_pointers(opt_pot.table);
#endif /* TERSOFF */

  /* write header data */
#ifdef STIWEB
  fprintf(outfile, "# Stillinger-Weber parameters fitted with %s\n", POTFIT_VERSION);
  fprintf(outfile, "# these entries are in LAMMPS \"metal\" units:\n");
  fprintf(outfile, "#   epsilon = eV; sigma = Angstroms\n");
  fprintf(outfile, "#   other quantities are unitless\n\n");
  fprintf(outfile, "# format of a single entry (one or more lines):\n");
  fprintf(outfile, "#   element_1 element_2 element_3\n");
  fprintf(outfile, "#   epsilon sigma a lambda gamma costheta0 A B p q tol\n\n");
#endif /* STIWEB */
#ifdef TERSOFF
#ifndef TERSOFFMOD
  fprintf(outfile, "# Tersoff parameters fitted with %s\n", POTFIT_VERSION);
  fprintf(outfile, "# these entries are in LAMMPS \"metal\" units:\n");
  fprintf(outfile, "#   A,B = eV; lambda1,2,3 = 1/Angstroms; R,D = Angstroms\n");
  fprintf(outfile, "#   other quantities are unitless\n\n");
  fprintf(outfile, "# format of a single entry (one or more lines):\n");
  fprintf(outfile, "#   element_1 element_2 element_3\n");
  fprintf(outfile, "#   m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A\n\n");
#else
  fprintf(outfile, "# modified Tersoff parameters fitted with %s\n", POTFIT_VERSION);
  fprintf(outfile, "# these entries are in LAMMPS \"metal\" units:\n");
  fprintf(outfile, "#   A,B = eV; lambda1,2,3 = 1/Angstroms; R,D = Angstroms\n");
  fprintf(outfile, "#   other quantities are unitless\n\n");
  fprintf(outfile, "# format of a single entry (one or more lines):\n");
  fprintf(outfile, "#   element_1 element_2 element_3\n");
  fprintf(outfile, "#   beta, alpha, h, eta, beta_ters, lambda2, B, R, D, lambda1, A, n\n#   c1, c2, c3, c4, c5\n\n");
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */

  /* write data, one line per triple of elements */
  for (i = 0; i < ntypes; i++) {
    for (j = 0; j < ntypes; j++) {
      for (k = 0; k < ntypes; k++) {
	col_j = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2) : j * ntypes + i - ((j * (j + 1)) / 2);
	fprintf(outfile, "%s %s %s", elements[i], elements[j], elements[k]);
#ifdef STIWEB
	/* calculate scaling factors */
	/* energy scaling factor */
	epsilon = 1.0;
	/* distance scaling factor */
	sigma = apot_table.values[col_j][4];
	power_1(&scale_p, &sigma, &apot_table.values[col_j][2]);
	power_1(&scale_q, &sigma, &apot_table.values[col_j][3]);
	scale_p = 1.0 / (epsilon * scale_p);
	scale_q = 1.0 / (epsilon * scale_q);
	fprintf(outfile, " %g", epsilon);	/* epsilon */
	fprintf(outfile, " %g", sigma);	/* sigma */
	fprintf(outfile, " %g", apot_table.values[col_j][5] / sigma);	/* a */
	fprintf(outfile, " %g", *(apot_table.sw.lambda[i][j][k]) / epsilon);	/* lambda */
	fprintf(outfile, " %g", apot_table.values[paircol + col_j][0] / sigma);	/* gamma */
	fprintf(outfile, " %g", -1.0 / 3.0);	/* costheta0 */
	fprintf(outfile, " %g", apot_table.values[col_j][1] * scale_q);	/* A */
	fprintf(outfile, " %g", apot_table.values[col_j][0] * scale_p / apot_table.values[col_j][1] / scale_q);	/* B */
	fprintf(outfile, " %g", apot_table.values[col_j][2]);	/* p */
	fprintf(outfile, " %g", apot_table.values[col_j][3]);	/* q */
	fprintf(outfile, " %g", 0.0);	/* tol */
#endif /* STIWEB */
#ifdef TERSOFF
#ifndef TERSOFFMOD
	/* potfit calculates the Tersoff_2 potential model */
	/* lammps requires the parameters in Tersoff_1 format */
	/* m has no effect, set it to 1.0 */
	fprintf(outfile, " %g", 1.0);
	fprintf(outfile, " %g", *(apot_table.tersoff.omega[col_j]));	/* gamma_ijk = omega_ik */
	fprintf(outfile, " %g", 0.0);	/* lambda3 = 0 */
	fprintf(outfile, " %g", apot_table.values[col_j][6]);	/* c */
	fprintf(outfile, " %g", apot_table.values[col_j][7]);	/* d */
	fprintf(outfile, " %g", apot_table.values[col_j][8]);	/* costheta0 = h */
	fprintf(outfile, " %g", apot_table.values[col_j][5]);	/* n */
	fprintf(outfile, " %g", apot_table.values[col_j][4]);	/* beta */
	fprintf(outfile, " %g", apot_table.values[col_j][3]);	/* lambda2 */
	fprintf(outfile, " %g", apot_table.values[col_j][1] * *(apot_table.tersoff.chi[col_j]));	/* B (includes chi) */
	fprintf(outfile, " %g", 0.5 * (apot_table.values[col_j][9] + apot_table.values[col_j][10]));	/* R */
	fprintf(outfile, " %g", 0.5 * (apot_table.values[col_j][9] - apot_table.values[col_j][10]));	/* D */
	fprintf(outfile, " %g", apot_table.values[col_j][2]);	/* lambda1 */
	fprintf(outfile, " %g", apot_table.values[col_j][0]);	/* A */
#else /* !TERSOFFMOD */
	/* order of parameters for a LAMMPS potential: */
	/* beta, alpha, h, eta, beta_ters (dummy), lambda2, B, R, D, lambda1, A, n, c1,..,c5 */
	fprintf(outfile, " %g", apot_table.values[col_j][7]);	/* beta */
	fprintf(outfile, " %g", apot_table.values[col_j][6]);	/* alpha */
	fprintf(outfile, " %g", apot_table.values[col_j][13]);	/* h */
	fprintf(outfile, " %g", apot_table.values[col_j][4]);	/* eta */
	fprintf(outfile, " %g", 1.0);	/* beta_ters (dummy) */
	fprintf(outfile, " %g", apot_table.values[col_j][3]);	/* lambda2 = mu */
	fprintf(outfile, " %g", apot_table.values[col_j][1]);	/* B */
	fprintf(outfile, " %g", 0.5 * (apot_table.values[col_j][14] + apot_table.values[col_j][15]));	/* R = 0.5 * (R1+R2) */
	fprintf(outfile, " %g", 0.5 * (apot_table.values[col_j][15] - apot_table.values[col_j][14]));	/* D = 0.5 * (R2-R1) */
	fprintf(outfile, " %g", apot_table.values[col_j][2]);	/* lambda1 = lambda */
	fprintf(outfile, " %g", apot_table.values[col_j][0]);	/* A */
	fprintf(outfile, " %g", 1.0 / (2.0 * apot_table.values[col_j][5]));	/* n = 1/(2*delta) */
	fprintf(outfile, " %g", apot_table.values[col_j][8]);	/* c1 */
	fprintf(outfile, " %g", apot_table.values[col_j][9]);	/* c2 */
	fprintf(outfile, " %g", apot_table.values[col_j][10]);	/* c3 */
	fprintf(outfile, " %g", apot_table.values[col_j][11]);	/* c4 */
	fprintf(outfile, " %g", apot_table.values[col_j][12]);	/* c5 */
#endif /* !TERSOFFMOD */
#endif /* TERSOFF */
	fprintf(outfile, "\n");
      }
    }
  }

  fclose(outfile);
  printf("Potential in LAMMPS format written to \t%s\n", filename);
}

#else

 /* in DYNAMO multi-element setfl format */
void write_pot_table_lammps(pot_table_t *pt)
{
#ifdef COULOMB
  printf("Potential in LAMMPS format is not available for coulomb interactions.\n");
  return;
#elif defined DIPOLE
  printf("Potential in LAMMPS format is not available for dipole interactions.\n");
  return;
#else /* COULOMB */

  FILE *outfile;
  char  filename[255];
  int   i, j;
  int   k = 0, l;
  double drho, dr, r, temp;

  /* open file */
  if (strcmp(output_prefix, "") != 0)
    sprintf(filename, "%s.lammps.%s", output_prefix, interaction_name);
  else
    sprintf(filename, "%s.lammps.%s", endpot, interaction_name);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* initialize periodic table */
  init_elements();

  /* write LAMMPS formated header */
  /* lines 1,2,3 = comments (ignored) */
  fprintf(outfile, "LAMMPS %s potential generated by potfit\n", interaction_name);
  fprintf(outfile, "potfit version %s (%s, %s)\n", POTFIT_VERSION, __DATE__, __TIME__);
  fprintf(outfile, "-----\n");
  /* line 4: Nelements Element1 Element2 ... ElementN */
  fprintf(outfile, "%d", ntypes);
  for (i = 0; i < ntypes; i++)
    fprintf(outfile, " %s", elements[i]);
  fprintf(outfile, "\n");

  temp = 999.9;
  for (i=0;i<ntypes;i++) {
    k = paircol + ntypes + i;
#ifdef APOT
    temp = MIN(temp, apot_table.end[k]);
#else
	temp = MIN(temp, opt_pot.end[k]);
#endif
  }

  drho = temp / (imdpotsteps - 1);
  dr = rcutmin / (imdpotsteps - 1);

  /* line 5: Nrho, drho, Nr, dr, cutoff */
  fprintf(outfile, "%d %f %d %f %f\n", imdpotsteps, drho, imdpotsteps, dr, rcutmin);

  /* one block for every atom type */
#if defined EAM || defined ADP
  for (i = 0; i < ntypes; i++) {
    /* atomic number, mass, lattice constant, lattice type */
    fprintf(outfile, "%3d %f 0 ???\n", ele_number_from_name(elements[i]), ele_mass_from_name(elements[i]));
    r = 0.0;
    /* embedding function F(n) */
    k = paircol + ntypes + i;
    for (j = 0; j < imdpotsteps; j++) {
#ifdef APOT
	apot_table.fvalue[k] (r, apot_table.values[k], &temp);
	temp =
	  smooth_pot[k] ? temp * cutoff(r, apot_table.end[k],
	  apot_table.values[k][apot_table.n_par[k] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, k, r));
#endif /* APOT */
      r += drho;
    }
    r = 0.0;
    k = paircol + i;
    /* transfer function rho(r) */
    for (j = 0; j < imdpotsteps; j++) {
#ifdef APOT
	apot_table.fvalue[k] (r, apot_table.values[k], &temp);
	temp =
	  smooth_pot[k] ? temp * cutoff(r, apot_table.end[k],
	  apot_table.values[k][apot_table.n_par[k] - 1]) : temp;
	fprintf(outfile, "%.16e\n", temp);
#else
	fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, k, r));
#endif /* APOT */
      r += dr;
    }
  }
#endif /* EAM || ADP */

  /* pair potentials */
  for (i = 0; i < ntypes; i++)
    for (j = 0; j <= i; j++) {
      r = 0.0;
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2) : j * ntypes + i - ((j * (j + 1)) / 2);
      for (l = 0; l < imdpotsteps; l++) {
#ifdef APOT
	apot_table.fvalue[k] (r, apot_table.values[k], &temp);
	temp =
	  smooth_pot[k] ? temp * cutoff(r, apot_table.end[k],
	  apot_table.values[k][apot_table.n_par[k] - 1]) : temp;
	fprintf(outfile, "%.16e\n", r * temp);
#else
	fprintf(outfile, "%.16e\n", r * splint_ne(pt, pt->table, k, r));
#endif /* APOT */
	r += dr;
      }
    }

#ifdef ADP
  /* dipole distortion */
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = 0.0;
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2) : j * ntypes + i - ((j * (j + 1)) / 2);
      k += paircol + 2 * ntypes;
      for (l = 0; l < imdpotsteps; l++) {
	fprintf(outfile, "%e\n", splint_ne(pt, pt->table, k, r));
	r += dr;
      }
    }
  /* quadrupole distortion */
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = 0.0;
      k = (i <= j) ? i * ntypes + j - ((i * (i + 1)) / 2) : j * ntypes + i - ((j * (j + 1)) / 2);
      k += 2 * (paircol + ntypes);
      for (l = 0; l < imdpotsteps; l++) {
	fprintf(outfile, "%e\n", splint_ne(pt, pt->table, k, r));
	r += dr;
      }
    }
#endif /* ADP */

  fclose(outfile);
  printf("Potential in LAMMPS format written to \t%s\n", filename);

  return;

#endif /* COULOMB */
}

#endif /* STIWEB || TERSOFF */

/****************************************************************
 *
 *  write alternate plot version of potential table
 *  (same intervals for all pair and transfer functions)
 *
 ****************************************************************/

void write_altplot_pair(pot_table_t *pt, char *filename)
{
  int   i, j, k, l;
  double r, rmin = 100.0, rmax = 0.0, r_step;
#if defined EAM || defined MEAM
  double temp;
#endif /* EAM || MEAM */
  FILE *outfile;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write data */
  k = 0;
  for (i = 0; i < ntypes; i++) {
    for (j = i; j < ntypes; j++) {
      rmin = MIN(rmin, pt->begin[k]);
      rmax = MAX(rmax, pt->end[k]);
      k++;
    }
    rmin = MIN(rmin, pt->begin[paircol + i]);
    rmax = MAX(rmax, pt->end[paircol + i]);
  }
  k = 0;
  r_step = (rmax - rmin) / (NPLOT - 1);
  for (i = 0; i < ntypes; i++)
    for (j = i; j < ntypes; j++) {
      r = rmin;
      for (l = 0; l < NPLOT - 1; l++) {
	fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r));
	r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
#if defined EAM || defined MEAM
  j = k;
  for (i = j; i < j + ntypes; i++) {
    r = rmin;
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, r <= pt->end[i] ? splint_ne(pt, pt->table, i, r) : 0);
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = j + ntypes; i < j + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT; l++) {
      temp = splint_ne(pt, pt->table, i, r);
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
#endif /* EAM || MEAM */
#ifdef MEAM
  j = k;
  for (i = j; i < j + ntypes; i++) {
    r = rmin;
    for (l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, r <= pt->end[i] ? splint_ne(pt, pt->table, i, r) : 0);
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (i = j + ntypes; i < j + 2 * ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (l = 0; l < NPLOT; l++) {
      temp = splint_ne(pt, pt->table, i, r);
      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
#endif /* MEAM */
  fclose(outfile);
  printf("Potential plotting data written to %s\n", filename);
}

#ifdef PDIST

/****************************************************************
 *
 * write_pairdist(pot_table_t *pt, char *filename)
 *    - write distribution function of function access
 *
 ****************************************************************/

void write_pairdist(pot_table_t *pt, char *filename)
{
  int  *freq;			/* frequency... */
  int   h, i, j, typ1, typ2, col;
#if defined(EAM)
  int k = 0;
  int l = 0;
#endif /* EAM */
  double rr;
  atom_t *atom;
  neigh_t *neigh;
  FILE *outfile;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* initialize distribution vector */
  freq = (int *)malloc(ndimtot * sizeof(int));
  for (i = 0; i < ndimtot; i++)
    freq[i] = 0;

  for (h = firstconf; h < firstconf + myconf; h++) {
    for (i = 0; i < inconf[h]; i++) {
      atom = atoms + i + cnfstart[h];
      typ1 = atom->type;

      /* pair potentials */
      for (j = 0; j < atom->num_neigh; j++) {
	neigh = atom->neigh + j;
	typ2 = neigh->type;
	col = (typ1 <= typ2) ? typ1 * ntypes + typ2 - ((typ1 * (typ1 + 1)) / 2)
	  : typ2 * ntypes + typ1 - ((typ2 * (typ2 + 1)) / 2);
	/* this has already been calculated */
	if (neigh->r < pt->end[col])
	  freq[neigh->slot[0]]++;
#ifdef EAM
	/* transfer function */
	col = paircol + typ2;
	if (neigh->r < pt->end[col])
	  freq[neigh->slot[1]]++;
#endif /* EAM */
      }
#ifdef EAM
      /* embedding function - get index first */
      col = paircol + ntypes + typ1;
      if (format == 3) {
	rr = atom->rho - pt->begin[col];
#ifdef RESCALE
	if (rr < 0.0)
	  error(1, "short distance");
	j = (int)(rr * pt->invstep[col]) + pt->first[col];
#else
	if (rr < 0.0)
	  rr = 0.0;		/* extrapolation */
	j = MIN((int)(rr * pt->invstep[col]) + pt->first[col], pt->last[col]);
#endif /* RESCALE */
      } else {			/* format ==4 */
	rr = atom->rho;
	k = pt->first[col];
	l = pt->last[col];
	while (l - k > 1) {
	  j = (k + l) >> 1;
	  if (pt->xcoord[j] > rr)
	    l = j;
	  else
	    k = j;
	}
	j = k;
      }
      freq[j]++;
#endif /* EAM */
    }
  }
  /* finished calculating data - write it to output file */
  j = 0;
  for (col = 0; col < pt->ncols; col++) {
    for (i = pt->first[col]; i < pt->last[col]; i++) {
      rr = 0.5 * (pt->xcoord[i] + pt->xcoord[i + 1]);
      fprintf(outfile, "%f %d\n", rr, freq[i]);
    }
    fprintf(outfile, "\n\n");
  }
  fclose(outfile);
  printf("Distribution data written to\t\t%s\n", filename);
}

#endif /* PDIST */

#ifdef APOT

void write_imd_data_pair(FILE *outfile, char *string, int offset, int y)
{
  int   i = 0;

  fprintf(outfile, "%s", string);
  for (i = 0; i < paircol; i++)
    fprintf(outfile, " %g", apot_table.values[offset + i][y]);
  fprintf(outfile, "\n");
}

#endif /* APOT */
