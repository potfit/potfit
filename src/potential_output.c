/****************************************************************
 *
 * potential_output.c: Routines for writing the potential table
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

#include "potfit.h"

#include "functions.h"
#include "memory.h"
#include "potential_output.h"
#include "splines.h"

#define NPLOT 1000

void write_pot_table0(char const* filename);
void write_pot_table3(char const* filename);
void write_pot_table4(char const* filename);
void write_pot_table5(char const* filename);

/****************************************************************
 *
 * write output potential
 *
 ****************************************************************/

void write_pot_table_potfit(char const* filename)
{
  switch (g_pot.format_type) {
    case POTENTIAL_FORMAT_UNKNOWN:
      error(1, "Unknown potential format detected! (%s:%d)\n", __FILE__,
            __LINE__);
    case POTENTIAL_FORMAT_ANALYTIC:
      write_pot_table0(filename);
      break;
    case POTENTIAL_FORMAT_TABULATED_EQ_DIST:
      write_pot_table3(filename);
    case POTENTIAL_FORMAT_TABULATED_NON_EQ_DIST:
      write_pot_table4(filename);
      break;
    case POTENTIAL_FORMAT_KIM:
      write_pot_table5(filename);
      break;
  }
}

/****************************************************************
 *
 *  write potential table (format 0)
 *
 ****************************************************************/

void write_pot_table0(char const* filename)
{
#if defined(APOT)
  apot_table_t* apt = &g_pot.apot_table;

  FILE* outfile = NULL;

  /* open file */
  outfile = fopen(filename, "w");
  if (NULL == outfile)
    error(1, "Could not open file %s\n", filename);

  /* write header */
  fprintf(outfile, "#F 0 %d", apt->number);
  fprintf(outfile, "\n#T %s", g_pot.interaction_name);

  /* write elements */
  fprintf(outfile, "\n#C");
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);

  /* write order of the interactions */
  fprintf(outfile, "\n##");
  /* pair potentials */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);

#if defined(ANG)
  /* f terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);
  /* g terms */
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);
#endif  // ANG

#if defined(EAM) || defined(ADP) || defined(MEAM)
  /* transfer functions */
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);
  /* embedding functions */
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);
#endif  // EAM || ADP || MEAM

#if defined(ADP)
  /* dipole terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);
  /* quadrupole terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);
#endif  // ADP

#if defined(MEAM)
  /* f terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);
  /* g terms */
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);
#endif  // MEAM

#if defined(STIWEB)
  /* stiweb_3 terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);
  /* lambda terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = 0; j < g_param.ntypes; j++)
      for (int k = j; k < g_param.ntypes; k++)
        fprintf(outfile, " %s-%s-%s", g_config.elements[i],
                g_config.elements[j], g_config.elements[k]);
#endif  // STIWEB

#if defined(TERSOFF)
  /* mixing terms */
  for (int i = 0; i < g_param.ntypes; i++)
    for (int j = i + 1; j < g_param.ntypes; j++)
      fprintf(outfile, " %s-%s", g_config.elements[i], g_config.elements[j]);
#endif  // TERSOFF

  /* write invariant switch for individual potentials */
  if (g_pot.have_invar) {
    fprintf(outfile, "\n#I");
    for (int i = 0; i < apt->number; i++)
      fprintf(outfile, " %d", g_pot.invar_pot[i]);
  }

  /* end tag */
  fprintf(outfile, "\n#E\n\n");

#if defined(PAIR)
  /* write chemical potentials if enabled */
  if (g_param.enable_cp) {
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(outfile, "cp_%s %.10f %.2f %.2f\n", g_config.elements[i],
              apt->chempot[i], apt->pmin[apt->number][i],
              apt->pmax[apt->number][i]);
#if defined(CN)
    if (compnodes > 0)
      fprintf(outfile, "cn %d\n", compnodes);
    for (j = 0; j < compnodes; j++)
      fprintf(outfile, "%.2f %.10f %.2f %.2f\n", compnodelist[j],
              apt->chempot[ntypes + j], apt->pmin[apt->number][ntypes + j],
              apt->pmax[apt->number][ntypes + j]);
#endif  // CN
    fprintf(outfile, "\n");
  }
#endif  // PAIR

#if defined(COULOMB)
  fprintf(outfile, "elstat\n");
  for (int i = 0; i < g_param.ntypes - 1; i++)
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number][i],
            apt->charge[i], apt->pmin[apt->number][i],
            apt->pmax[apt->number][i]);
  fprintf(outfile, "charge_%s\t %f\n", g_config.elements[g_param.ntypes - 1],
          apt->last_charge);
  fprintf(outfile, "%s\t\t %f\t %f\t %f\n", apt->param_name[apt->number + 1][0],
          apt->dp_kappa[0], apt->pmin[apt->number + 1][0],
          apt->pmax[apt->number + 1][0]);

#if defined(DIPOLE)
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 2][i],
            apt->dp_alpha[i], apt->pmin[apt->number + 2][i],
            apt->pmax[apt->number + 2][i]);
  for (int i = 0; i < g_calc.paircol; i++) {
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 3][i],
            apt->dp_b[i], apt->pmin[apt->number + 3][i],
            apt->pmax[apt->number + 3][i]);
  }
  for (int i = 0; i < g_calc.paircol; i++) {
    fprintf(outfile, "%s\t %f\t %f\t %f\n", apt->param_name[apt->number + 4][i],
            apt->dp_c[i], apt->pmin[apt->number + 4][i],
            apt->pmax[apt->number + 4][i]);
  }
#endif  // DIPOLE

  fprintf(outfile, "\n");
#endif  // COULOMB

  /* write globals section */
  if (g_pot.have_globals) {
    fprintf(outfile, "global\t%d\n", apt->globals);
    for (int i = 0; i < apt->globals; i++)
      fprintf(outfile, "%s\t%18.8f\t%12.4f\t%12.4f\n",
              apt->param_name[g_pot.global_pot][i],
              apt->values[g_pot.global_pot][i], apt->pmin[g_pot.global_pot][i],
              apt->pmax[g_pot.global_pot][i]);
    fprintf(outfile, "\n");
  }

  /* write data */
  for (int i = 0; i < apt->number; i++) {
    if (g_pot.smooth_pot[i]) {
      fprintf(outfile, "type\t%s_sc\n", apt->names[i]);
    } else {
      fprintf(outfile, "type\t%s\n", apt->names[i]);
    }

    fprintf(outfile, "cutoff\t%f\n", apt->end[i]);
    fprintf(outfile, "# rmin\t%f\n", apt->begin[i]);

    for (int j = 0; j < apt->n_par[i]; j++) {
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
#endif  // APOT
}

/****************************************************************
 *
 *  write potential table (format 3)
 *
 ****************************************************************/

void write_pot_table3(char const* filename)
{
  FILE* pfile_table = NULL;
  FILE* pfile_plot = NULL;
  int plot_flag = 0;
  pot_table_t* pt = &g_pot.opt_pot;

  if (g_files.plotpointfile != NULL)
    plot_flag = 1;

  /* open file */
  pfile_table = fopen(filename, "w");
  if (pfile_table == NULL)
    error(1, "Could not open file %s\n", filename);

  /* if needed: open file for plotpoints */
  if (plot_flag) {
    pfile_plot = fopen(g_files.plotpointfile, "w");
    if (pfile_table == NULL)
      error(1, "Could not open file %s\n", filename);
  }

  /* write header */
  fprintf(pfile_table, "#F 3 %d", pt->ncols);
  fprintf(pfile_table, "\n#T %s", g_pot.interaction_name);

  if (g_config.elements != NULL) {
    fprintf(pfile_table, "\n#C");
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
    fprintf(pfile_table, "\n##");

    for (int i = 0; i < g_param.ntypes; i++)
      for (int j = i; j < g_param.ntypes; j++)
        fprintf(pfile_table, " %s-%s", g_config.elements[i],
                g_config.elements[j]);
#if defined(EAM) || defined(MEAM)
    /* transfer functions */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
    /* embedding functions */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
#if defined(TBEAM)
    /* transfer functions s-band */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
    /* embedding functions s-band */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
#endif  // TBEAM
#endif  // EAM || MEAM

#if defined(MEAM)
    /* pre-anglpart */
    for (int i = 0; i < g_param.ntypes; i++)
      for (int j = i; j < g_param.ntypes; j++)
        fprintf(pfile_table, " %s-%s", g_config.elements[i],
                g_config.elements[j]);
    /* angl part */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
#endif  // MEAM
  }

  if (g_pot.have_invar) {
    fprintf(pfile_table, "\n#I");
    for (int i = 0; i < pt->ncols; i++)
      fprintf(pfile_table, " %d", g_pot.invar_pot[i]);
  }

  fprintf(pfile_table, "\n#G");

  for (int i = 0; i < pt->ncols; i++)
    fprintf(pfile_table, " %d", g_pot.gradient[i]);

  fprintf(pfile_table, "\n#E\n");

  /* write info block */
  for (int i = 0; i < pt->ncols; i++)
    fprintf(pfile_table, "%.16e %.16e %d\n", pt->begin[i], pt->end[i],
            pt->last[i] - pt->first[i] + 1);
  fprintf(pfile_table, "\n");

  /* write data */
  for (int i = 0; i < pt->ncols; i++) {
    double r = pt->begin[i];
    /* write gradient */
    fprintf(pfile_table, "%.16e %.16e\n", pt->table[pt->first[i] - 2],
            pt->table[pt->first[i] - 1]);

    for (int j = pt->first[i]; j <= pt->last[i]; j++) {
      fprintf(pfile_table, "%.16e\n", pt->table[j]);
      if (plot_flag)
        fprintf(pfile_plot, "%.6e %.6e %d\n", r, pt->table[j], j);
      r += pt->step[i];
    }
    fprintf(pfile_table, "\n");

    if (plot_flag)
      fprintf(pfile_plot, "\n\n");
  }

  fclose(pfile_table);

  if (plot_flag)
    fclose(pfile_plot);
}

/****************************************************************
 *
 *  write potential table (format 4)
 *
 ****************************************************************/

void write_pot_table4(char const* filename)
{
  FILE* pfile_table = NULL;
  FILE* pfile_plot = NULL;
  int plot_flag = 0;
  pot_table_t* pt = &g_pot.opt_pot;

  if (g_files.plotpointfile != NULL)
    plot_flag = 1;

  /* open file */
  pfile_table = fopen(filename, "w");
  if (pfile_table == NULL)
    error(1, "Could not open file %s\n", filename);

  /* if needed: open file for plotpoints */
  if (plot_flag) {
    pfile_plot = fopen(g_files.plotpointfile, "w");
    if (pfile_plot == NULL)
      error(1, "Could not open file %s\n", filename);
  }

  /* write header */
  fprintf(pfile_table, "#F 4 %d", pt->ncols);
  fprintf(pfile_table, "\n#T %s", g_pot.interaction_name);

  if (g_config.elements != NULL) {
    fprintf(pfile_table, "\n#C");

    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);

    fprintf(pfile_table, "\n##");

    for (int i = 0; i < g_param.ntypes; i++)
      for (int j = i; j < g_param.ntypes; j++)
        fprintf(pfile_table, " %s-%s", g_config.elements[i],
                g_config.elements[j]);

#if defined(EAM)
    /* transfer functions */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);

    /* embedding functions */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);

#if defined(TBEAM)
    /* transfer functions s-band */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);

    /* embedding functions s-band */
    for (int i = 0; i < g_param.ntypes; i++)
      fprintf(pfile_table, " %s", g_config.elements[i]);
#endif  // TBEAM
#endif  // EAM
  }

  if (g_pot.have_invar) {
    fprintf(pfile_table, "\n#I");
    for (int i = 0; i < pt->ncols; i++)
      fprintf(pfile_table, " %d", g_pot.invar_pot[i]);
  }

  fprintf(pfile_table, "\n#G");

  for (int i = 0; i < pt->ncols; i++)
    fprintf(pfile_table, " %d", g_pot.gradient[i]);

  fprintf(pfile_table, "\n#E\n");

  /* write info block */
  for (int i = 0; i < pt->ncols; i++)
    fprintf(pfile_table, "%d\n", pt->last[i] - pt->first[i] + 1);

  fprintf(pfile_table, "\n");

  /* write data */
  for (int i = 0; i < pt->ncols; i++) {
    fprintf(pfile_table, "%.16e %.16e\n", pt->table[pt->first[i] - 2],
            pt->table[pt->first[i] - 1]);

    for (int j = pt->first[i]; j <= pt->last[i]; j++) {
      fprintf(pfile_table, "%.16e %.16e\n", pt->xcoord[j], pt->table[j]);

      if (plot_flag)
        fprintf(pfile_plot, "%.6e %.6e %d\n", pt->xcoord[j], pt->table[j], j);
    }
    fprintf(pfile_table, "\n");

    if (plot_flag)
      fprintf(pfile_plot, "\n\n");
  }
  fclose(pfile_table);

  if (plot_flag)
    fclose(pfile_plot);
}

/****************************************************************
 *
 *  write potential table (format 5)
 *
 ****************************************************************/

void write_pot_table5(char const* filename)
{
#if defined(KIM)
  pot_table_t* pt = &g_pot.opt_pot;

  // open file
  FILE* outfile = fopen(filename, "w");
  if (outfile == NULL)
    error(1, "Could not open file %s\n", filename);

  // write header
  fprintf(outfile, "#F 5 1");
  fprintf(outfile, "\n#C");
  for (int i = 0; i < g_param.ntypes; i++)
    fprintf(outfile, " %s", g_config.elements[i]);
  fprintf(outfile, "\n#E\n\n");

  // write KIM Model name
  fprintf(outfile, "# KIM Model name\n");
  fprintf(outfile, "model\t%s\n\n", g_kim.model_name);

  // write cutoff
  fprintf(outfile, "# cutoff distance\n");
  fprintf(outfile, "cutoff\t%f\n\n", g_config.rcutmax);

  // number of opt params
  fprintf(outfile, "# the number of optimizable parameters that will be listed below\n");
  fprintf(outfile, "# use 'kim_opt_param list' to get a list of supported parameters for the current model\n");
  fprintf(outfile, "kim_opt_param\t%d\n\n", g_kim.num_opt_param);

  // write data
  int k = 0;
  fprintf(outfile, "# parameters\n");
  for (int i = 0; i < g_kim.num_opt_param; i++) {
    fprintf(outfile, "%s\n", g_kim.freeparams.name[g_kim.idx_opt_param[i]]);
    for (int j = 0; j < g_kim.size_opt_param[i]; j++) {
      fprintf(outfile, "  %f\t", pt->table[k]);
      fprintf(outfile, "%f\t%f\n", g_pot.apot_table.pmin[0][k], g_pot.apot_table.pmax[0][k]);
      k++;
    }
  }

  fclose(outfile);
#endif
}

/****************************************************************
  write plot version of potential table
****************************************************************/

void write_plotpot_pair(pot_table_t* pt, const char* filename)
{
  double r = 0.0;
  double r_step = 0.0;

  // open file
  FILE* outfile = fopen(filename, "w");
  if (outfile == NULL)
    error(1, "Could not open file %s for writing\n", filename);

#if !defined(APOT)
  int k = 0;

  // write pair data
  for (int i = 0; i < g_param.ntypes; i++) {
    for (int j = i; j < g_param.ntypes; j++) {
      r = pt->begin[k];
      r_step = (pt->end[k] - r) / (NPLOT - 1);
      for (int l = 0; l < NPLOT - 1; l++) {
        fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, k, r));
        r += r_step;
      }
      fprintf(outfile, "%e %e\n\n\n", r, 0.0);
      k++;
    }
  }

#if defined(EAM) || defined(ADP) || defined(MEAM)
  for (int i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (int l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (int i = g_calc.paircol + g_param.ntypes;
       i < g_calc.paircol + 2 * g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (int l = 0; l < NPLOT; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
#endif  // EAM || ADP || MEAM

#if defined(MEAM)
  for (int i = g_calc.paircol; i < g_calc.paircol + g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (int l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (int i = g_calc.paircol + g_param.ntypes;
       i < g_calc.paircol + 2 * g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (int l = 0; l < NPLOT; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "\n\n\n");
  }
  for (int i = g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * g_calc.paircol + 2 * g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (int l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (int i = 2 * g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * g_calc.paircol + 3 * g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - pt->begin[i]) / (NPLOT - 1);
    for (int l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
#endif  // MEAM

#if defined(ADP)
  for (int i = g_calc.paircol + 2 * g_param.ntypes;
       i < 2 * (g_calc.paircol + g_param.ntypes); i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - r) / (NPLOT - 1);
    for (int l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
  for (int i = 2 * (g_calc.paircol + g_param.ntypes);
       i < 3 * g_calc.paircol + 2 * g_param.ntypes; i++) {
    r = pt->begin[i];
    r_step = (pt->end[i] - r) / (NPLOT - 1);
    for (int l = 0; l < NPLOT - 1; l++) {
      fprintf(outfile, "%e %e\n", r, splint_ne(pt, pt->table, i, r));
      r += r_step;
    }
    fprintf(outfile, "%e %e\n\n\n", r, 0.0);
  }
#endif  // ADP

#else   // APOT

  for (int i = 0; i < g_pot.apot_table.number; i++) {
    r = (g_param.plotmin == 0 ? 0.1 : g_param.plotmin);
    r_step = (g_pot.apot_table.end[i] - r) / (NPLOT - 1);
    double h = g_pot.apot_table.values[i][g_pot.apot_table.n_par[i] - 1];
    for (int j = 0; j < NPLOT; j++) {
      double temp = 0.0;

      (*g_pot.apot_table.fvalue[i])(r, g_pot.apot_table.values[i], &temp);

      if (g_pot.smooth_pot[i])
        temp *= apot_cutoff(r, g_pot.apot_table.end[i], h);

      if (isnan(temp))
        temp = 10e30;

      fprintf(outfile, "%e %e\n", r, temp);
      r += r_step;
    }
    if (i != (g_pot.apot_table.number - 1))
      fprintf(outfile, "\n\n");
  }
#endif  // APOT

  fclose(outfile);
  printf("Potential plotting data written to \t%s\n", filename);
}

#if defined(PDIST)

/****************************************************************
  write_pairdist(pot_table_t *pt, char *filename)
    - write distribution function of function access
****************************************************************/

void write_pairdist(pot_table_t* pt, const char* filename)
{
  // open file
  FILE* outfile = fopen(filename, "w");
  if (outfile == NULL)
    error(1, "Could not open file %s\n", filename);

  // initialize distribution vector
  int freq[g_calc.ndimtot];
  memset(freq, 0, sizeof(freq));

  for (int h = g_mpi.firstconf; h < g_mpi.firstconf + g_mpi.myconf; h++) {
    for (int i = 0; i < g_config.inconf[h]; i++) {
      atom_t* atom = g_config.atoms + i + g_config.cnfstart[h];
      int typ1 = atom->type;

      // pair potentials
      for (int j = 0; j < atom->num_neigh; j++) {
        neigh_t* neigh = atom->neigh + j;
        int typ2 = neigh->type;
        int col = (typ1 <= typ2)
                  ? typ1 * g_param.ntypes + typ2 - ((typ1 * (typ1 + 1)) / 2)
                  : typ2 * g_param.ntypes + typ1 - ((typ2 * (typ2 + 1)) / 2);
        // this has already been calculated
        if (neigh->r < pt->end[col])
          freq[neigh->slot[0]]++;
#if defined(EAM)
        // transfer function
        col = g_calc.paircol + typ2;
        if (neigh->r < pt->end[col])
          freq[neigh->slot[1]]++;
#endif  // EAM
      }
#if defined(EAM)
      // embedding function - get index first
      int col = g_calc.paircol + g_param.ntypes + typ1;
      int j = 0;
      if (g_pot.format_type == POTENTIAL_FORMAT_TABULATED_EQ_DIST) {
        double rr = atom->rho - pt->begin[col];
#if defined(RESCALE)
        if (rr < 0.0)
          error(1, "short distance\n");
        j = (int)(rr * pt->invstep[col]) + pt->first[col];
#else
        if (rr < 0.0)
          rr = 0.0; /* extrapolation */
        j = MIN((int)(rr * pt->invstep[col]) + pt->first[col], pt->last[col]);
#endif         // RESCALE
      } else { /* format ==4 */
        int k = pt->first[col];
        int l = pt->last[col];
        while (l - k > 1) {
          j = (k + l) >> 1;
          if (pt->xcoord[j] > atom->rho)
            l = j;
          else
            k = j;
        }
        j = k;
      }
      freq[j]++;
#endif  // EAM
    }
  }

  // finished calculating data - write it to output file
  for (int col = 0; col < pt->ncols; col++) {
    for (int i = pt->first[col]; i < pt->last[col]; i++) {
      double rr = 0.5 * (pt->xcoord[i] + pt->xcoord[i + 1]);
      fprintf(outfile, "%f %d\n", rr, freq[i]);
    }
    fprintf(outfile, "\n\n");
  }
  fclose(outfile);
  printf("Distribution data written to\t\t%s\n", filename);
}

#endif  // PDIST

#if defined(COULOMB)

/****************************************************************
 *
 * write coulomb-potential
 *
 ****************************************************************/

void write_coulomb_table()
{
  apot_table_t* apt = &g_pot.apot_table;

  if (g_param.ntypes == 2) {
    double value, c1;
    FILE* outfile;
    char* filename1 = "Coulomb_00";
    char* filename2 = "Coulomb_01";
    char* filename3 = "Coulomb_11";

    c1 = -apt->charge[0] / 2;

    outfile = fopen(filename1, "a");

    for (int i = 0; i < g_config.natoms; i++) {
      for (int j = 0; j < g_config.atoms[i].num_neigh; j++) {
        value = apt->charge[0] * apt->charge[0] *
                g_config.atoms[i].neigh[j].fnval_el;
        fprintf(outfile, "%f\t%f\n", g_config.atoms[i].neigh[j].r, value);
      }
    }

    fclose(outfile);

    outfile = fopen(filename2, "a");

    for (int i = 0; i < g_config.natoms; i++) {
      for (int j = 0; j < g_config.atoms[i].num_neigh; j++) {
        value = apt->charge[0] * c1 * g_config.atoms[i].neigh[j].fnval_el;
        fprintf(outfile, "%f\t%f\n", g_config.atoms[i].neigh[j].r, value);
      }
    }

    fclose(outfile);

    outfile = fopen(filename3, "a");

    for (int i = 0; i < g_config.natoms; i++) {
      for (int j = 0; j < g_config.atoms[i].num_neigh; j++) {
        value = c1 * c1 * g_config.atoms[i].neigh[j].fnval_el;
        fprintf(outfile, "%f\t%f\n", g_config.atoms[i].neigh[j].r, value);
      }
    }
    fclose(outfile);
  } else {
    printf("Coulomb-outfiles are only available in case of two atom types.");
  }
}

#endif  // COULOMB
