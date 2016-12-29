/****************************************************************
 *
 * potential_output.c: Routines for writing the potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2016 - the potfit development team
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

#include "potfit.h"

#include "elements.h"
#include "functions.h"
#include "memory.h"
#include "potential_output.h"
#include "splines.h"
#include "utils.h"

#define NPLOT 1000

/****************************************************************
 *
 *  write potential table for IMD (format 2)
 *
 ****************************************************************/

#if defined STIWEB || defined TERSOFF

void write_pot_table_imd(char const* prefix)
{
  char filename[255];

  FILE* outfile = NULL;

#if defined(STIWEB)
  sprintf(filename, "%s.sw.pot", prefix);
#endif  // STIWEB

#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
  sprintf(filename, "%s.tersoff.pot", prefix);
#else
  sprintf(filename, "%s.tersoffmod.pot", prefix);
#endif  // !TERSOFFMOD
#endif  // TERSOFF

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  fprintf(outfile, "# IMD potential file written by %s\n#\n", POTFIT_VERSION);
#if defined(STIWEB)
  fprintf(outfile,
          "# This is a Stillinger-Weber potential. Compile IMD with the "
          "'stiweb' option.\n#\n");
#endif  // STIWEB
#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
  fprintf(outfile,
          "# This is a Tersoff potential. Compile IMD with the 'tersoff2' "
          "option.\n#\n");
#else
  fprintf(outfile,
          "# This is a modified Tersoff potential. Compile IMD with the "
          "'tersoffmod2' option.\n#\n");
#endif  // !TERSOFFMOD
#endif  // TERSOFF

  /* write warning header */
  fprintf(outfile, "# WARNING:\n");
  fprintf(outfile, "# DO NOT USE THIS FILE AS A POTENTIAL FILE IN IMD !!!!!\n");
  fprintf(outfile,
          "# COPY THE CONTENTS OF THIS FILE INTO THE IMD PARAMETER FILE\n\n");

#if defined(STIWEB)
  write_imd_data_pair(outfile, "stiweb_a\t", 0, 0);               /* A_ij */
  write_imd_data_pair(outfile, "stiweb_b\t", 0, 1);               /* B_ij */
  write_imd_data_pair(outfile, "stiweb_p\t", 0, 2);               /* p_ij */
  write_imd_data_pair(outfile, "stiweb_q\t", 0, 3);               /* q_ij */
  write_imd_data_pair(outfile, "stiweb_de\t", 0, 4);              /* delta_ij */
  write_imd_data_pair(outfile, "stiweb_a1\t", 0, 5);              /* a1_ij */
  write_imd_data_pair(outfile, "stiweb_ga\t", g_calc.paircol, 0); /* gamma_ij */
  write_imd_data_pair(outfile, "stiweb_a2\t", g_calc.paircol, 1); /* a2_ij */

  /* lambda_ijk */
  /* strange format (e.g. binary): 000, 100, 001, 101, 011, 111 */
  fprintf(outfile, "stiweb_la\t");
  for (int j = 0; j < g_calc.paircol; j++)
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(outfile, " %g",
              g_pot.apot_table
                  .values[g_pot.apot_table.number - 1][i * g_calc.paircol + j]);
  fprintf(outfile, "\n");
#endif  // STIWEB

#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
  write_imd_data_pair(outfile, "ters_a\t\t", 0, 0);     /* A_ij */
  write_imd_data_pair(outfile, "ters_b\t\t", 0, 1);     /* B_ij */
  write_imd_data_pair(outfile, "ters_la\t\t", 0, 2);    /* lambda_ij */
  write_imd_data_pair(outfile, "ters_mu\t\t", 0, 3);    /* mu_ij */
  write_imd_data_pair(outfile, "ters_ga\t\t", 0, 4);    /* gamma_ij */
  write_imd_data_pair(outfile, "ters_n\t\t", 0, 5);     /* n_ij */
  write_imd_data_pair(outfile, "ters_c\t\t", 0, 6);     /* c_ij */
  write_imd_data_pair(outfile, "ters_d\t\t", 0, 7);     /* d_ij */
  write_imd_data_pair(outfile, "ters_h\t\t", 0, 8);     /* h_ij */
  write_imd_data_pair(outfile, "ters_r_cut\t\t", 0, 9); /* r_cut = S_ij */
  write_imd_data_pair(outfile, "ters_r0\t\t", 0, 10);   /* r0 = R_ij */

  /* chi and omega only for mixing potentials */
  if (g_param.ntypes > 1) {
    /* chi_ij, chi_ii = 1.0 */
    fprintf(outfile, "ters_chi\t\t");
    for (int i = 0; i < g_param.ntypes * (g_param.ntypes - 1) / 2.0; i++)
      fprintf(outfile, " %g", g_pot.apot_table.values[g_calc.paircol + i][0]);
    fprintf(outfile, "\n");

    /* omega_ij, chi_ii = 1.0 */
    fprintf(outfile, "ters_om\t\t");
    for (int i = 0; i < g_param.ntypes * (g_param.ntypes - 1) / 2.0; i++)
      fprintf(outfile, " %g", g_pot.apot_table.values[g_calc.paircol + i][1]);
    fprintf(outfile, "\n");
  }
#else
  write_imd_data_pair(outfile, "ters_a\t\t", 0, 0);      /* A_ij */
  write_imd_data_pair(outfile, "ters_b\t\t", 0, 1);      /* B_ij */
  write_imd_data_pair(outfile, "ters_la\t\t", 0, 2);     /* lambda_ij */
  write_imd_data_pair(outfile, "ters_mu\t\t", 0, 3);     /* mu_ij */
  write_imd_data_pair(outfile, "ters_eta\t\t", 0, 4);    /* eta_ij */
  write_imd_data_pair(outfile, "ters_delta\t\t", 0, 5);  /* delta_ij */
  write_imd_data_pair(outfile, "ters_alpha\t\t", 0, 6);  /* alpha_ij */
  write_imd_data_pair(outfile, "ters_beta\t\t", 0, 7);   /* beta_ij */
  write_imd_data_pair(outfile, "ters_c1\t\t", 0, 8);     /* c1_ij */
  write_imd_data_pair(outfile, "ters_c2\t\t", 0, 9);     /* c2_ij */
  write_imd_data_pair(outfile, "ters_c3\t\t", 0, 10);    /* c3_ij */
  write_imd_data_pair(outfile, "ters_c4\t\t", 0, 11);    /* c4_ij */
  write_imd_data_pair(outfile, "ters_c5\t\t", 0, 12);    /* c5_ij */
  write_imd_data_pair(outfile, "ters_h\t\t", 0, 13);     /* h_ij */
  write_imd_data_pair(outfile, "ters_r0\t\t", 0, 14);    /* R1_ij */
  write_imd_data_pair(outfile, "ters_r_cut\t\t", 0, 15); /* R2_ij */
#endif  // !TERSOFFMOD
#endif  // TERSOFF

  fclose(outfile);
  printf("Parameters for IMD potential written to\t%s\n", filename);
}

#else

void write_pot_table_imd(char const* prefix)
{
  char filename[255];
  pot_table_t* pt = &g_pot.opt_pot;

  /* allocate memory */
  double r2begin[g_param.ntypes * g_param.ntypes];
  double r2end[g_param.ntypes * g_param.ntypes];
  double r2step[g_param.ntypes * g_param.ntypes];

  /* pair potential part (over r^2) */
  sprintf(filename, "%s_phi.imd.pt", prefix);

  /* open file */
  FILE* outfile = fopen(filename, "w");

  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes * g_param.ntypes);

  /* write info block */
  int m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 =
          i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      int col2 = i * g_param.ntypes + j;
/* Extrapolation possible  */
#if defined(APOT)
      r2begin[col2] = dsquare((g_param.plotmin == 0 ? 0.1 : g_param.plotmin));
#else
      r2begin[col2] =
          dsquare(MAX(pt->begin[col1] - g_param.extend * pt->step[col1], 0));
#endif  // APOT
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / g_param.imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
              r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 =
          i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      int col2 = i * g_param.ntypes + j;
      double r2 = r2begin[col2];
      for (int k = 0; k < g_param.imdpotsteps; k++) {
#if defined(APOT)
        double temp = 0;
        g_pot.apot_table.fvalue[col1](sqrt(r2), g_pot.apot_table.values[col1],
                                      &temp);
        temp =
            g_pot.smooth_pot[col1]
                ? temp *
                      apot_cutoff(
                          sqrt(r2), g_pot.apot_table.end[col1],
                          g_pot.apot_table
                              .values[col1][g_pot.apot_table.n_par[col1] - 1])
                : temp;
        fprintf(outfile, "%.16e\n", temp);
#else
#if defined(MEAM)
        fprintf(outfile, "%.16e\n",
                splint_ne_lin(pt, pt->table, col1, sqrt(r2)));
#else
        fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif  // MEAM
#endif  // APOT
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
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes * g_param.ntypes);

  /* write info block */
  for (int i = 0; i < g_param.ntypes; i++) {
    for (int j = 0; j < g_param.ntypes; j++) {
      int col1 = (g_param.ntypes * (g_param.ntypes + 1)) / 2 + j;
      int col2 = i * g_param.ntypes + j;
#if defined(APOT)
      r2begin[col2] = dsquare((g_param.plotmin == 0 ? 0.1 : g_param.plotmin));
#else
      /* Extrapolation possible  */
      r2begin[col2] =
          dsquare(MAX(pt->begin[col1] - g_param.extend * pt->step[col1], 0));
#endif  // APOT
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / g_param.imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
              r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  for (int i = 0; i < g_param.ntypes; i++) {
    for (int j = 0; j < g_param.ntypes; j++) {
      int col1 = (g_param.ntypes * (g_param.ntypes + 1)) / 2 + j;
      int col2 = i * g_param.ntypes + j;
      double r2 = r2begin[col2];
      for (int k = 0; k < g_param.imdpotsteps; k++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[col1])(sqrt(r2),
                                         g_pot.apot_table.values[col1], &temp);
        temp =
            g_pot.smooth_pot[col1]
                ? temp *
                      apot_cutoff(
                          sqrt(r2), g_pot.apot_table.end[col1],
                          g_pot.apot_table
                              .values[col1][g_pot.apot_table.n_par[col1] - 1])
                : temp;
        fprintf(outfile, "%.16e\n", temp);
#else
        fprintf(outfile, "%.16e\n",
                splint_ne_lin(pt, pt->table, col1, sqrt(r2)));
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
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes);

  /* write info block */
  for (int i = 0; i < g_param.ntypes; i++) {
    int col1 = (g_param.ntypes * (g_param.ntypes + 3)) / 2 + i;
#if defined(APOT)
    r2begin[i] = 0;
    r2end[i] = pt->end[col1];
#else
    /* pad with zeroes */
    r2begin[i] = pt->begin[col1] - g_param.extend * pt->step[col1];
    /* extrapolation */
    r2end[i] = pt->end[col1] + g_param.extend * pt->step[col1];
#endif  // APOT
    r2step[i] = (r2end[i] - r2begin[i]) / g_param.imdpotsteps;
    fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[i], r2end[i], r2step[i]);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (int i = 0; i < g_param.ntypes; i++) {
    double r2 = r2begin[i];
    int col1 = (g_param.ntypes * (g_param.ntypes + 3)) / 2 + i;
    double root = (pt->begin[col1] > 0.0)
                      ? pt->table[pt->first[col1]] / sqrt(pt->begin[col1])
                      : 0.0;
    root += (pt->end[col1] < 0.0)
                ? pt->table[pt->last[col1]] / sqrt(-pt->end[col1])
                : 0.0;
    double temp = 0.0;
    for (int k = 0; k <= g_param.imdpotsteps; k++) {
#if defined(APOT)
      (*g_pot.apot_table.fvalue[col1])(r2, g_pot.apot_table.values[col1],
                                       &temp);
#else
      temp = splint_ne_lin(pt, pt->table, col1, r2);
      double temp2 = r2 - pt->end[col1];
      temp += (temp2 > 0.0) ? 5e2 * (temp2 * temp2 * temp2) : 0.0;
#endif  // APOT
      fprintf(outfile, "%.16e\n", temp);
      r2 += r2step[i];
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  printf("IMD embedding function written to \t%s\n", filename);
#endif  // EAM || ADP

#if defined(ADP)
  /* write dipole function (over r^2) */
  sprintf(filename, "%s_upot.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes * g_param.ntypes);

  /* write info block */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 =
          i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      col1 += g_calc.paircol + 2 * g_param.ntypes;
      int col2 = i * g_param.ntypes + j;
/* Extrapolation possible  */
#if defined(APOT)
      r2begin[col2] = dsquare((g_param.plotmin == 0 ? 0.1 : g_param.plotmin));
#else
      r2begin[col2] =
          dsquare(MAX(pt->begin[col1] - g_param.extend * pt->step[col1], 0));
#endif  // APOT
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / g_param.imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
              r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 =
          i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      col1 += g_calc.paircol + 2 * g_param.ntypes;
      int col2 = i * g_param.ntypes + j;
      double r2 = r2begin[col2];
      for (int k = 0; k < g_param.imdpotsteps; k++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[col1])(sqrt(r2),
                                         g_pot.apot_table.values[col1], &temp);
        temp =
            g_pot.smooth_pot[col1]
                ? temp *
                      apot_cutoff(
                          sqrt(r2), g_pot.apot_table.end[col1],
                          g_pot.apot_table
                              .values[col1][g_pot.apot_table.n_par[col1] - 1])
                : temp;
        fprintf(outfile, "%.16e\n", temp);
#else
        fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif  // APOT
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
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes * g_param.ntypes);

  /* write info block */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 =
          i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      col1 += 2 * g_calc.paircol + 2 * g_param.ntypes;
      int col2 = i * g_param.ntypes + j;
/* Extrapolation possible  */
#if defined(APOT)
      r2begin[col2] = dsquare((g_param.plotmin == 0 ? 0.1 : g_param.plotmin));
#else
      r2begin[col2] =
          dsquare(MAX(pt->begin[col1] - g_param.extend * pt->step[col1], 0));
#endif  // APOT
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / g_param.imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
              r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 =
          i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      col1 += 2 * g_calc.paircol + 2 * g_param.ntypes;
      int col2 = i * g_param.ntypes + j;
      double r2 = r2begin[col2];
      for (int k = 0; k < g_param.imdpotsteps; k++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[col1])(sqrt(r2),
                                         g_pot.apot_table.values[col1], &temp);
        temp =
            g_pot.smooth_pot[col1]
                ? temp *
                      apot_cutoff(
                          sqrt(r2), g_pot.apot_table.end[col1],
                          g_pot.apot_table
                              .values[col1][g_pot.apot_table.n_par[col1] - 1])
                : temp;
        fprintf(outfile, "%.16e\n", temp);
#else
        fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, sqrt(r2)));
#endif  // APOT
        r2 += r2step[col2];
      }
      fprintf(outfile, "%.16e\n", 0.0);
      fprintf(outfile, "\n");
    }
  }
  fclose(outfile);
  printf("IMD quadrupole potential written to \t%s\n", filename);
#endif  // APOT

#if defined(MEAM)
  /* write f_r2 for MEAM */
  sprintf(filename, "%s_f_meam.imd.pt", prefix);
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes * g_param.ntypes);

  /* write info block */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 = (g_param.ntypes * (g_param.ntypes + 5)) / 2;
      col1 += i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      int col2 = i * g_param.ntypes + j;
/* Extrapolation possible  */
#if defined(APOT)
      r2begin[col2] = dsquare((g_param.plotmin == 0 ? 0.1 : g_param.plotmin));
#else
      r2begin[col2] =
          dsquare(MAX(pt->begin[col1] - g_param.extend * pt->step[col1], 0));
#endif  // APOT
      r2end[col2] = dsquare(pt->end[col1]);
      r2step[col2] = (r2end[col2] - r2begin[col2]) / g_param.imdpotsteps;
      fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[col2], r2end[col2],
              r2step[col2]);
    }
  }
  fprintf(outfile, "\n");

  /* write data */
  m = 0;
  for (int i = 0; i < g_param.ntypes; i++) {
    m += i;
    int m2 = 0;
    for (int j = 0; j < g_param.ntypes; j++) {
      m2 += j;
      int col1 = (g_param.ntypes * (g_param.ntypes + 5)) / 2;
      col1 += i < j ? i * g_param.ntypes + j - m : j * g_param.ntypes + i - m2;
      int col2 = i * g_param.ntypes + j;
      double r2 = r2begin[col2];
      for (int k = 0; k < g_param.imdpotsteps; k++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[col1])(sqrt(r2),
                                         g_pot.apot_table.values[col1], &temp);
        temp =
            g_pot.smooth_pot[col1]
                ? temp *
                      apot_cutoff(
                          sqrt(r2), g_pot.apot_table.end[col1],
                          g_pot.apot_table
                              .values[col1][g_pot.apot_table.n_par[col1] - 1])
                : temp;
        fprintf(outfile, "%.16e\n", temp);
#else
        fprintf(outfile, "%.16e\n",
                splint_ne_lin(pt, pt->table, col1, sqrt(r2)));
#endif  // APOT
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
  fprintf(outfile, "#F 2 %d\n#E\n", g_param.ntypes);

  /* write info block */
  for (int i = 0; i < g_param.ntypes; i++) {
    int col1 = (g_param.ntypes * (g_param.ntypes + 3)) + i;
    /* from -1 to +1 */
    r2begin[i] = pt->begin[col1];
    r2end[i] = pt->end[col1];
    r2step[i] = (r2end[i] - r2begin[i]) / g_param.imdpotsteps;
    fprintf(outfile, "%.16e %.16e %.16e\n", r2begin[i], r2end[i], r2step[i]);
  }
  fprintf(outfile, "\n");

  /* write data */
  for (int i = 0; i < g_param.ntypes; i++) {
    double r2 = r2begin[i];
    int col1 = (g_param.ntypes * (g_param.ntypes + 3)) + i;
    for (int k = 0; k <= g_param.imdpotsteps; k++) {
      fprintf(outfile, "%.16e\n", splint_ne(pt, pt->table, col1, r2));
      r2 += r2step[i];
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
  printf("IMD MEAM g potential data written to\t%s\n", filename);
#endif  // MEAM

/* write endpot for IMD with electrostatics */
#if defined(COULOMB) && defined(APOT)
  apot_table_t* apt = &g_pot.apot_table;

  int ncols = g_param.ntypes * (g_param.ntypes + 1) / 2;

  sprintf(filename, "%s_charges.imd", prefix);

  /* open file */
  outfile = fopen(filename, "w");

  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  fprintf(outfile, "charge\t\t");

  for (int i = 0; i < g_param.ntypes - 1; i++)
    fprintf(outfile, "%f\t", apt->charge[i]);

  fprintf(outfile, "%f\n", apt->last_charge);

  if ((strcmp(apt->names[0], "ms") == 0) &&
      (strcmp(apt->names[1], "ms") == 0) &&
      (strcmp(apt->names[2], "ms") == 0)) {
    fprintf(outfile, "ms_D\t\t");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][0]);
    fprintf(outfile, "\nms_gamma\t");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][1]);
    fprintf(outfile, "\nms_r0\t\t");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][2]);
    fprintf(outfile, "\n");
  } else if ((strcmp(apt->names[0], "buck") == 0) &&
             (strcmp(apt->names[1], "buck") == 0) &&
             (strcmp(apt->names[2], "buck") == 0)) {
    fprintf(outfile, "buck_a\t");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][0]);
    fprintf(outfile, "\nbuck_sigma\t");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][1]);
    fprintf(outfile, "\nbuck_c\t");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, "%f\t", apt->values[i][2]);
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "\new_rcut\t\t%f\n", g_config.dp_cut);
  fprintf(outfile, "ew_kappa\t\t%f\n", apt->dp_kappa[0]);
  fprintf(outfile, "r_cut\t\t");
  for (int i = 0; i < ncols; i++)
    fprintf(outfile, "%f\t", g_pot.apot_table.end[0]);
  fprintf(outfile, "\n\n");

#if defined(DIPOLE)
  fprintf(outfile, "dp_alpha\t");
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, "%f\t", apt->dp_alpha[i]);
  fprintf(outfile, "\ndp_b\t\t");
  for (int i = 0; i < ncols; i++)
    fprintf(outfile, "%f\t", apt->dp_b[i]);
  fprintf(outfile, "\ndp_c\t\t");
  for (int i = 0; i < ncols; i++)
    fprintf(outfile, "%f\t", apt->dp_c[i]);
  fprintf(outfile, "\n");
#endif  // DIPOLE

  fclose(outfile);
  printf("Electrostatic table for IMD written to \t%s\n", filename);
#endif  // COULOMB && APOT
}

#endif  // STIWEB

#if defined(APOT)

void write_imd_data_pair(FILE* outfile, char* string, int offset, int y)
{
  int i = 0;

  fprintf(outfile, "%s", string);
  for (i = 0; i < g_calc.paircol; i++)
    fprintf(outfile, " %g", g_pot.apot_table.values[offset + i][y]);
  fprintf(outfile, "\n");
}

#endif  // APOT
