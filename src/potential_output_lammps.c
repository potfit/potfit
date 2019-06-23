/****************************************************************
 *
 * potential_output.c: Routines for writing the potential table
 *
 ****************************************************************
 *
 * Copyright 2002-2018 - the potfit development team
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
#include "force.h"
#include "functions.h"
#include "potential_output.h"
#include "splines.h"
#include "utils.h"

#define NPLOT 1000

#if defined STIWEB || defined TERSOFF

void write_pot_table_lammps()
{
  FILE* outfile;
  char filename[255];
  int i, j, k;
  int col_j;
#if defined(STIWEB)
  double epsilon, sigma, scale_p, scale_q;
#endif  // STIWEB

  if (g_config.elements == NULL) {
    warning("There are no elements listed in you configuration file.\n");
    warning("A LAMMPS potential cannot be written without them.\n");
    return;
  }

#if defined(STIWEB)
  /* check if final potential is LAMMPS compliant (a1==a2) */
  if (g_pot.apot_table.values[0][5] !=
      g_pot.apot_table.values[g_calc.paircol][1]) {
    warning("Your potential is not supported by LAMMPS.\n");
    warning("Please ensure that the values a1 and a2 are the same,\n");
    warning("  either by using global parameters or by fixing them.\n");

    return;
  }
#endif  // STIWEB

  /* open file */
  if (strcmp(g_files.output_prefix, "") != 0)
    strcpy(filename, g_files.output_prefix);
  else
    strcpy(filename, g_files.endpot);
#if defined(STIWEB)
  strcat(filename, ".lammps.sw");
#endif  // STIWEB
#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
  strcat(filename, ".lammps.tersoff");
#else
  strcat(filename, ".lammps.tersoffmod");
#endif  // !TERSOFFMOD
#endif  // TERSOFF
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

#if defined(STIWEB)
  update_stiweb_pointers(g_pot.opt_pot.table);
#endif  // STIWEB
#if defined(TERSOFF)
  update_tersoff_pointers(g_pot.opt_pot.table);
#endif  // TERSOFF

/* write header data */
#if defined(STIWEB)
  fprintf(outfile, "# Stillinger-Weber parameters fitted with %s\n",
          POTFIT_VERSION);
  fprintf(outfile, "# these entries are in LAMMPS \"metal\" units:\n");
  fprintf(outfile, "#   epsilon = eV; sigma = Angstroms\n");
  fprintf(outfile, "#   other quantities are unitless\n\n");
  fprintf(outfile, "# format of a single entry (one or more lines):\n");
  fprintf(outfile, "#   element_1 element_2 element_3\n");
  fprintf(outfile,
          "#   epsilon sigma a lambda gamma costheta0 A B p q tol\n\n");
#endif  // STIWEB
#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
  fprintf(outfile, "# Tersoff parameters fitted with %s\n", POTFIT_VERSION);
  fprintf(outfile, "# these entries are in LAMMPS \"metal\" units:\n");
  fprintf(outfile,
          "#   A,B = eV; lambda1,2,3 = 1/Angstroms; R,D = Angstroms\n");
  fprintf(outfile, "#   other quantities are unitless\n\n");
  fprintf(outfile, "# format of a single entry (one or more lines):\n");
  fprintf(outfile, "#   element_1 element_2 element_3\n");
  fprintf(outfile,
          "#   m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, "
          "lambda1, A\n\n");
#else
  fprintf(outfile, "# modified Tersoff parameters fitted with %s\n",
          POTFIT_VERSION);
  fprintf(outfile, "# these entries are in LAMMPS \"metal\" units:\n");
  fprintf(outfile,
          "#   A,B = eV; lambda1,2,3 = 1/Angstroms; R,D = Angstroms\n");
  fprintf(outfile, "#   other quantities are unitless\n\n");
  fprintf(outfile, "# format of a single entry (one or more lines):\n");
  fprintf(outfile, "#   element_1 element_2 element_3\n");
  fprintf(outfile,
          "#   beta, alpha, h, eta, beta_ters, lambda2, B, R, D, lambda1, A, "
          "n\n#   c1, c2, c3, c4, c5\n\n");
#endif  // !TERSOFFMOD
#endif  // TERSOFF

  /* write data, one line per triple of elements */
  for (i = 0; i < g_param.ntypes; i++) {
    for (j = 0; j < g_param.ntypes; j++) {
      for (k = 0; k < g_param.ntypes; k++) {
        col_j = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                         : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
        fprintf(outfile, "%s %s %s", g_config.elements[i], g_config.elements[j],
                g_config.elements[k]);
#if defined(STIWEB)
        /* calculate scaling factors */
        /* energy scaling factor */
        epsilon = 1.0;
        /* distance scaling factor */
        sigma = g_pot.apot_table.values[col_j][4];
        power_1(&scale_p, &sigma, &g_pot.apot_table.values[col_j][2]);
        power_1(&scale_q, &sigma, &g_pot.apot_table.values[col_j][3]);
        scale_p = 1.0 / (epsilon * scale_p);
        scale_q = 1.0 / (epsilon * scale_q);
        fprintf(outfile, " %g", epsilon); /* epsilon */
        fprintf(outfile, " %g", sigma);   /* sigma */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][5] / sigma); /* a */
        fprintf(outfile, " %g",
                *(g_pot.apot_table.sw.lambda[i][j][k]) / epsilon); /* lambda */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[g_calc.paircol + col_j][0] /
                    sigma);                  /* gamma */
        fprintf(outfile, " %g", -1.0 / 3.0); /* costheta0 */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][1] * scale_q); /* A */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][0] * scale_p /
                                    g_pot.apot_table.values[col_j][1] /
                                    scale_q);                       /* B */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][2]); /* p */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][3]); /* q */
        fprintf(outfile, " %g", 0.0);                               /* tol */
#endif                                                              // STIWEB
#if defined(TERSOFF)
#if !defined(TERSOFFMOD)
        /* potfit calculates the Tersoff_2 potential model */
        /* lammps requires the parameters in Tersoff_1 format */
        /* m has no effect, set it to 1.0 */
        fprintf(outfile, " %g", 1.0);
        fprintf(outfile, " %g",
                *(g_pot.apot_table.tersoff
                      .omega[col_j])); /* gamma_ijk = omega_ik */
        fprintf(outfile, " %g", 0.0);  /* lambda3 = 0 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][6]); /* c */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][7]); /* d */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][8]); /* costheta0 = h */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][5]); /* n */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][4]); /* beta */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][3]); /* lambda2 */
        fprintf(
            outfile, " %g",
            g_pot.apot_table.values[col_j][1] *
                *(g_pot.apot_table.tersoff.chi[col_j])); /* B (includes chi) */
        fprintf(outfile, " %g",
                0.5 * (g_pot.apot_table.values[col_j][9] +
                       g_pot.apot_table.values[col_j][10])); /* R */
        fprintf(outfile, " %g",
                0.5 * (g_pot.apot_table.values[col_j][9] -
                       g_pot.apot_table.values[col_j][10])); /* D */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][2]); /* lambda1 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][0]); /* A */
#else                                                               /* !TERSOFFMOD */
        /* order of parameters for a LAMMPS potential: */
        /* beta, alpha, h, eta, beta_ters (dummy), lambda2, B, R, D, lambda1, A,
         * n, c1,..,c5 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][7]);  /* beta */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][6]);  /* alpha */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][13]); /* h */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][4]);  /* eta */
        fprintf(outfile, " %g", 1.0); /* beta_ters (dummy) */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][3]); /* lambda2 = mu */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][1]); /* B */
        fprintf(
            outfile, " %g",
            0.5 * (g_pot.apot_table.values[col_j][14] +
                   g_pot.apot_table.values[col_j][15])); /* R = 0.5 * (R1+R2) */
        fprintf(
            outfile, " %g",
            0.5 * (g_pot.apot_table.values[col_j][15] -
                   g_pot.apot_table.values[col_j][14])); /* D = 0.5 * (R2-R1) */
        fprintf(outfile, " %g",
                g_pot.apot_table.values[col_j][2]); /* lambda1 = lambda */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][0]); /* A */
        fprintf(
            outfile, " %g",
            1.0 / (2.0 *
                   g_pot.apot_table.values[col_j][5])); /* n = 1/(2*delta) */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][8]);  /* c1 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][9]);  /* c2 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][10]); /* c3 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][11]); /* c4 */
        fprintf(outfile, " %g", g_pot.apot_table.values[col_j][12]); /* c5 */
#endif  // !TERSOFFMOD
#endif  // TERSOFF
        fprintf(outfile, "\n");
      }
    }
  }

  fclose(outfile);
  printf("Potential in LAMMPS format written to \t%s\n", filename);
}

#else

#if defined(PAIR)

void write_lammps_table_pair()
{
  char filename[255];

  /* open file */
  if (strcmp(g_files.output_prefix, "") != 0)
    sprintf(filename, "%s.lammps.%s", g_files.output_prefix,
            g_pot.interaction_name);
  else
    sprintf(filename, "%s.lammps.%s", g_files.endpot, g_pot.interaction_name);

  FILE* outfile = fopen(filename, "w");

  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write LAMMPS formated header */
  /* comments (ignored) */
  fprintf(outfile, "# LAMMPS %s potential generated by potfit\n",
          g_pot.interaction_name);
  fprintf(outfile, "# potfit version %s (%s, %s)\n\n", POTFIT_VERSION, __DATE__,
          __TIME__);

  int idx = 0;
  for (int i = 0; i < g_param.ntypes; ++i) {
    for (int j = i; j < g_param.ntypes; ++j, ++idx) {
#if defined(APOT)
      fprintf(outfile, "%s_", g_pot.apot_table.names[idx]);
#else // APOT
      fprintf(outfile, "POTENTIAL_");
#endif // APOT
      fprintf(outfile, "%s_%s\n", g_config.elements[i], g_config.elements[j]);
      fprintf(outfile, "N %d\n\n", g_param.imdpotsteps);

      // write pot data here
#if defined(APOT)
      double r = g_param.plotmin == 0 ? g_pot.calc_pot.begin[idx] : g_param.plotmin;
#else
      double r = g_pot.calc_pot.begin[idx];
#endif
      double dr = (g_pot.calc_pot.end[idx] - r) / (g_param.imdpotsteps - 1);

      for (int k = 0; k < g_param.imdpotsteps; ++k, r += dr) {
        double fval = 0.0;
        double grad = 0.0;

#if defined(APOT)
        double a = 0.0;
        double b = 0.0;
        (*g_pot.apot_table.fvalue[idx])(r + dr, g_pot.apot_table.values[idx], &a);
        (*g_pot.apot_table.fvalue[idx])(r, g_pot.apot_table.values[idx], &fval);
        (*g_pot.apot_table.fvalue[idx])(r - dr, g_pot.apot_table.values[idx], &b);

        if (fval == 0.0 && (isnan(fval) || isinf(fval)))
          g_pot.apot_table.fvalue[idx](r + 0.001, g_pot.apot_table.values[idx], &fval);

        if (g_pot.smooth_pot[idx]) {
          a *= apot_cutoff(r + dr, g_pot.apot_table.end[idx], g_pot.apot_table.values[idx][g_pot.apot_table.n_par[idx] - 1]);
          fval *= apot_cutoff(r, g_pot.apot_table.end[idx], g_pot.apot_table.values[idx][g_pot.apot_table.n_par[idx] - 1]);
          b *= apot_cutoff(r - dr, g_pot.apot_table.end[idx], g_pot.apot_table.values[idx][g_pot.apot_table.n_par[idx] - 1]);
        }

        grad = - (a - b) / (2.0 * dr);
#else
        fval = g_splint_comb(&g_pot.calc_pot, g_pot.calc_pot.table, idx, r, &grad);
#endif  // APOT
        fprintf(outfile, "%d %f %.16e %.16e\n", k + 1, r, fval, grad);
      }
      fprintf(outfile, "\n");
    }
  }

  fclose(outfile);
  printf("Potential in LAMMPS format written to \t%s\n", filename);
}

#endif // PAIR

#if defined(EAM) || defined(ADP)

void write_lammps_table_eam_adp()
{
  char filename[255];

  /* open file */
  if (strcmp(g_files.output_prefix, "") != 0)
    sprintf(filename, "%s.lammps.%s", g_files.output_prefix,
            g_pot.interaction_name);
  else
    sprintf(filename, "%s.lammps.%s", g_files.endpot, g_pot.interaction_name);

  FILE* outfile = fopen(filename, "w");

  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write LAMMPS formated header */
  /* lines 1,2,3 = comments (ignored) */
  fprintf(outfile, "LAMMPS %s potential generated by potfit\n",
          g_pot.interaction_name);
  fprintf(outfile, "potfit version %s (%s, %s)\n", POTFIT_VERSION, __DATE__,
          __TIME__);
  fprintf(outfile, "-----\n");
  /* line 4: Nelements Element1 Element2 ... ElementN */
  fprintf(outfile, "%d", g_param.ntypes);
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);
  fprintf(outfile, "\n");

  double temp = 999.9;

  for (int i = 0; i < g_param.ntypes; i++) {
    int k = g_calc.paircol + g_param.ntypes + i;
#if defined(APOT)
    temp = MIN(temp, g_pot.apot_table.end[k]);
#else
    temp = MIN(temp, g_pot.opt_pot.end[k]);
#endif
  }

  double drho = temp / (g_param.imdpotsteps - 1);
  double dr = g_config.rcutmin / (g_param.imdpotsteps - 1);

  /* line 5: Nrho, drho, Nr, dr, cutoff */
  fprintf(outfile, "%d %f %d %f %f\n", g_param.imdpotsteps, drho,
          g_param.imdpotsteps, dr, g_config.rcutmin);

/* one block for every atom type */
#if defined EAM || defined ADP
  for (int i = 0; i < g_param.ntypes; i++) {
    /* atomic number, mass, lattice constant, lattice type */
    fprintf(outfile, "%3d %f 0 ???\n",
            ele_number_from_name(g_config.elements[i]),
            ele_mass_from_name(g_config.elements[i]));
    double r = 0.0;
    /* embedding function F(n) */
    int k = g_calc.paircol + g_param.ntypes + i;
    for (int j = 0; j < g_param.imdpotsteps; j++) {
#if defined(APOT)
      (*g_pot.apot_table.fvalue[k])(r, g_pot.apot_table.values[k], &temp);
      if (r == 0.0 && (isnan(temp) || isinf(temp)))
        g_pot.apot_table.fvalue[k](r + 0.001, g_pot.apot_table.values[k],
                                   &temp);
      if (g_pot.smooth_pot[k]) {
        temp *= apot_cutoff(
            r, g_pot.apot_table.end[k],
            g_pot.apot_table.values[k][g_pot.apot_table.n_par[k] - 1]);
      }
      fprintf(outfile, "%.16e\n", temp);
#else
      fprintf(outfile, "%.16e\n",
              g_splint(&g_pot.calc_pot, g_pot.calc_pot.table, k, r));
#endif  // APOT
      r += drho;
    }
    r = 0.0;
    k = g_calc.paircol + i;
    /* transfer function rho(r) */
    for (int j = 0; j < g_param.imdpotsteps; j++) {
#if defined(APOT)
      (*g_pot.apot_table.fvalue[k])(r, g_pot.apot_table.values[k], &temp);
      if (r == 0.0 && (isnan(temp) || isinf(temp)))
        g_pot.apot_table.fvalue[k](r + 0.001, g_pot.apot_table.values[k],
                                   &temp);
      if (g_pot.smooth_pot[k]) {
        temp *= apot_cutoff(
            r, g_pot.apot_table.end[k],
            g_pot.apot_table.values[k][g_pot.apot_table.n_par[k] - 1]);
      }
      fprintf(outfile, "%.16e\n", temp);
#else
      fprintf(outfile, "%.16e\n",
              g_splint(&g_pot.calc_pot, g_pot.calc_pot.table, k, r));
#endif  // APOT
      r += dr;
    }
  }
#endif  // EAM || ADP

  /* pair potentials */
  for (int i = 0; i < g_param.ntypes; i++) {
    for (int j = 0; j <= i; j++) {
      double r = 0.0;
      int k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                       : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      for (int l = 0; l < g_param.imdpotsteps; l++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[k])(r, g_pot.apot_table.values[k], &temp);
        if (r == 0.0 && (isnan(temp) || isinf(temp)))
          g_pot.apot_table.fvalue[k](r + 0.001, g_pot.apot_table.values[k],
                                     &temp);
        if (g_pot.smooth_pot[k]) {
          temp *= apot_cutoff(
              r, g_pot.apot_table.end[k],
              g_pot.apot_table.values[k][g_pot.apot_table.n_par[k] - 1]);
        }
        fprintf(outfile, "%.16e\n", r * temp);
#else
        fprintf(outfile, "%.16e\n",
                r * g_splint(&g_pot.calc_pot, g_pot.calc_pot.table, k, r));
#endif  // APOT
        r += dr;
      }
    }
  }

#if defined(ADP)
  /* dipole distortion */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++) {
      double r = 0.0;
      int k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                       : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      k += g_calc.paircol + 2 * g_param.ntypes;
      for (int l = 0; l < g_param.imdpotsteps; l++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[k])(r, g_pot.apot_table.values[k], &temp);
        if (r == 0.0 && (isnan(temp) || isinf(temp)))
          g_pot.apot_table.fvalue[k](r + 0.001, g_pot.apot_table.values[k],
                                     &temp);
        if (g_pot.smooth_pot[k]) {
          temp *= apot_cutoff(
              r, g_pot.apot_table.end[k],
              g_pot.apot_table.values[k][g_pot.apot_table.n_par[k] - 1]);
        }
        fprintf(outfile, "%.16e\n", r * temp);
#else
        fprintf(outfile, "%e\n",
                g_splint(&g_pot.calc_pot, g_pot.calc_pot.table, k, r));
#endif // APOT
        r += dr;
      }
    }
  /* quadrupole distortion */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++) {
      double r = 0.0;
      int k = (i <= j) ? i * g_param.ntypes + j - ((i * (i + 1)) / 2)
                       : j * g_param.ntypes + i - ((j * (j + 1)) / 2);
      k += 2 * (g_calc.paircol + g_param.ntypes);
      for (int l = 0; l < g_param.imdpotsteps; l++) {
#if defined(APOT)
        double temp = 0.0;
        (*g_pot.apot_table.fvalue[k])(r, g_pot.apot_table.values[k], &temp);
        if (r == 0.0 && (isnan(temp) || isinf(temp)))
          g_pot.apot_table.fvalue[k](r + 0.001, g_pot.apot_table.values[k],
                                     &temp);
        if (g_pot.smooth_pot[k]) {
          temp *= apot_cutoff(
              r, g_pot.apot_table.end[k],
              g_pot.apot_table.values[k][g_pot.apot_table.n_par[k] - 1]);
        }
        fprintf(outfile, "%.16e\n", r * temp);
#else
        fprintf(outfile, "%e\n",
                g_splint(&g_pot.calc_pot, g_pot.calc_pot.table, k, r));
#endif // APOT
        r += dr;
      }
    }
#endif  // ADP

  fclose(outfile);
  printf("Potential in LAMMPS format written to \t%s\n", filename);
}

#endif // EAM || ADP

/* in DYNAMO multi-element setfl format */
void write_pot_table_lammps()
{
#if defined(PAIR)
  write_lammps_table_pair();
#elif defined(EAM) || defined(ADP)
  write_lammps_table_eam_adp();
#endif

  printf("Writing a LAMMPS potential for %s is not supported\n", g_pot.interaction_name);
}

#endif  // STIWEB || TERSOFF
