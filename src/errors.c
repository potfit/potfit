/****************************************************************
 *
 * errors.c: Error reports for the individual components
 *
 ****************************************************************
 *
 * Copyright 2002-2015
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

#include "utils.h"

void write_errors(double* force, double tot)
{
  int i, j;
  char file[255];
  FILE* outfile;
  double sqr;
  double rms[3];
  char* component[6];

/* set spline density corrections to 0 */
#if defined EAM || defined ADP || defined MEAM
  double* totdens = (double*)malloc(g_param.ntypes * sizeof(double));
  reg_for_free(totdens, "totdens");
#endif /* EAM || ADP || MEAM */

#if defined EAM || defined ADP || defined MEAM
#ifndef MPI
  /* Not much sense in printing rho when not communicated... */
  if (g_param.write_output_files)
  {
    strcpy(file, g_files.output_prefix);
    strcat(file, ".rho_loc");
    outfile = fopen(file, "w");
    if (NULL == outfile)
      error(1, "Could not open file %s\n", file);
  }
  else
  {
    outfile = stdout;
    printf("Local electron density rho\n");
  }
  for (i = 0; i < g_param.ntypes; i++)
  {
    totdens[i] = 0.0;
  }
  fprintf(outfile, "#    atomtype\trho\n");
#ifdef MEAM
  fprintf(outfile, "#    atomtype\trho\trho_eam\trho_meam\n");
#endif /* MEAM */
  for (i = 0; i < g_config.natoms; i++)
  {
#if defined EAM || defined ADP
    fprintf(outfile, "%d\t%d\t%f\n", i, g_config.atoms[i].type, g_config.atoms[i].rho);
#elif defined MEAM
    fprintf(outfile, "%d\t%d\t%f\t%f\t%f\n", i, g_config.atoms[i].type,
            g_config.atoms[i].rho, g_config.atoms[i].rho_eam,
            g_config.atoms[i].rho - g_config.atoms[i].rho_eam);
#endif /* EAM || ADP */
    totdens[g_config.atoms[i].type] += g_config.atoms[i].rho;
  }
  fprintf(outfile, "\n");
  for (i = 0; i < g_param.ntypes; i++)
  {
    totdens[i] /= (double)g_config.na_type[g_config.nconf][i];
    fprintf(outfile, "Average local electron density at atom sites type %d: %f\n", i,
            totdens[i]);
  }
  if (g_param.write_output_files)
  {
    printf("Local electron density data written to \t%s\n", file);
    fclose(outfile);
  }
#endif /* !MPI */
#endif /* EAM || ADP || MEAM */

/* prepare for error calculations */
#ifdef CONTRIB
  int contrib_atoms = 0;
#endif /* CONTRIB */
  double f_sum = 0.0;
  double e_sum = 0.0;
  double s_sum = 0.0;

  /* write force deviations */
  if (g_param.write_output_files)
  {
    strcpy(file, g_files.output_prefix);
    strcat(file, ".force");
    outfile = fopen(file, "w");
    if (NULL == outfile)
      error(1, "Could not open file %s\n", file);
  }
  else
  {
    outfile = stdout;
    printf("Forces:\n");
  }

  for (i = 0; i < 6; i++)
  {
    component[i] = (char*)malloc(3 * sizeof(char));
    if (NULL == component[i])
      error(1, "Could not allocate memory for component strings");
  }

  strcpy(component[0], "x");
  strcpy(component[1], "y");
  strcpy(component[2], "z");

  for (i = 0; i < 3 * g_config.natoms; i++)
  {
#ifdef CONTRIB
    if (0 == g_config.atoms[i / 3].contrib)
      sqr = 0.0;
    else
#endif /* CONTRIB */
      sqr = g_config.conf_weight[g_config.atoms[i / 3].conf] * dsquare(force[i]);
    f_sum += sqr;
#ifdef FWEIGHT
    if (i > 2 && i % 3 == 0 &&
        g_config.atoms[i / 3].conf != g_config.atoms[i / 3 - 1].conf)
      fprintf(outfile, "\n\n");
    if (i == 0)
      fprintf(outfile, "#conf:atom\ttype\tdf^2\t\tf\t\tf0\t\tdf/f0\t\t|f|\n");
    fprintf(
        outfile, "%3d:%6d:%s\t%4s\t%20.18f\t%11.6f\t%11.6f\t%14.8f\t%14.8f\n",
        g_config.atoms[i / 3].conf, i / 3, component[i % 3],
        g_config.elements[g_config.atoms[i / 3].type], sqr,
        force[i] * (FORCE_EPS + g_config.atoms[i / 3].absforce) + g_config.force_0[i],
        g_config.force_0[i],
        (force[i] * (FORCE_EPS + g_config.atoms[i / 3].absforce)) / g_config.force_0[i],
        g_config.atoms[i / 3].absforce);
#else
    if (i > 2 && i % 3 == 0 &&
        g_config.atoms[i / 3].conf != g_config.atoms[i / 3 - 1].conf)
      fprintf(outfile, "\n\n");
    if (i == 0)
      fprintf(outfile, "#conf:atom\ttype\tdf^2\t\tf\t\tf0\t\tdf/f0\n");
    fprintf(outfile, "%3d:%6d:%s\t%4s\t%e\t%e\t%e\t%e\n", g_config.atoms[i / 3].conf,
            i / 3, component[i % 3], g_config.elements[g_config.atoms[i / 3].type], sqr,
            force[i] + g_config.force_0[i], g_config.force_0[i],
            force[i] / g_config.force_0[i]);
#endif /* FWEIGHT */
  }
  if (g_param.write_output_files)
  {
    printf("Force data written to \t\t\t%s\n", file);
    fclose(outfile);
  }

  /* write energy deviations */
  if (g_param.eweight != 0)
  {
    if (g_param.write_output_files)
    {
      strcpy(file, g_files.output_prefix);
      strcat(file, ".energy");
      outfile = fopen(file, "w");
      if (NULL == outfile)
        error(1, "Could not open file %s\n", file);
    }
    else
    {
      outfile = stdout;
      printf("Cohesive Energies\n");
    }

    if (g_param.write_output_files)
    {
      fprintf(outfile, "# global energy weight w is %f\n", g_param.eweight);
      fprintf(outfile, "# nr.\tconf_w\tw*de^2\t\te\t\te0\t\t|e-e0|\t\te-e0\t\tde/e0\n");
    }
    else
    {
      fprintf(outfile, "energy weight is %f\n", g_param.eweight);
      fprintf(outfile, "conf\tconf_w\t(w*de)^2\te\t\te0\t\tde/e0\n");
    }

    for (i = 0; i < g_config.nconf; i++)
    {
      sqr =
          g_config.conf_weight[i] * g_param.eweight * dsquare(force[g_calc.energy_p + i]);
      e_sum += sqr;
      if (g_param.write_output_files)
      {
        fprintf(outfile, "%3d\t%6.2f\t%10.6f\t%13.10f\t%13.10f\t%f\t%f\t%f\n", i,
                g_config.conf_weight[i], sqr / g_config.conf_weight[i],
                force[g_calc.energy_p + i] + g_config.force_0[g_calc.energy_p + i],
                g_config.force_0[g_calc.energy_p + i], fabs(force[g_calc.energy_p + i]),
                force[g_calc.energy_p + i],
                force[g_calc.energy_p + i] / g_config.force_0[g_calc.energy_p + i]);
      }
      else
        fprintf(outfile, "%d\t%.4f\t%f\t%f\t%f\t%f\n", i, g_config.conf_weight[i], sqr,
                force[g_calc.energy_p + i] + g_config.force_0[g_calc.energy_p + i],
                g_config.force_0[g_calc.energy_p + i],
                force[g_calc.energy_p + i] / g_config.force_0[g_calc.energy_p + i]);
    }
    if (g_param.write_output_files)
    {
      printf("Energy data written to \t\t\t%s\n", file);
      fclose(outfile);
    }
  }
  else
  {
    printf("Energy data not written (energy weight was 0).\n");
  }

#if defined(STRESS)
  /* write stress deviations */
  if (g_param.sweight != 0)
  {
    if (g_param.write_output_files)
    {
      strcpy(file, g_files.output_prefix);
      strcat(file, ".stress");
      outfile = fopen(file, "w");
      if (NULL == outfile)
        error(1, "Could not open file %s\n", file);
      fprintf(outfile, "# global stress weight w is %f\n", g_param.sweight);
    }
    else
    {
      outfile = stdout;
      fprintf(outfile, "Stresses on unit cell\n");
    }
    strcpy(component[0], "xx");
    strcpy(component[1], "yy");
    strcpy(component[2], "zz");
    strcpy(component[3], "xy");
    strcpy(component[4], "yz");
    strcpy(component[5], "zx");

    fprintf(outfile, "#\tconf_w\tw*ds^2\t\ts\t\ts0\t\tds/s0\n");

    for (i = g_calc.stress_p; i < g_calc.stress_p + 6 * g_config.nconf; i++)
    {
      sqr = g_config.conf_weight[(i - g_calc.stress_p) / 6] * g_param.sweight *
            dsquare(force[i]);
      s_sum += sqr;
      fprintf(outfile, "%3d-%s\t%6.2f\t%14.8f\t%12.10f\t%12.10f\t%14.8f\n",
              (i - g_calc.stress_p) / 6, component[(i - g_calc.stress_p) % 6],
              g_config.conf_weight[(i - g_calc.stress_p) / 6], sqr,
              force[i] + g_config.force_0[i], g_config.force_0[i],
              force[i] / g_config.force_0[i]);
    }
    if (g_param.write_output_files)
    {
      printf("Stress data written to \t\t\t%s\n", file);
      fclose(outfile);
    }
  }
  else
  {
    printf("Stress data not written (stress weight was 0).\n");
  }
#endif /* STRESS */

#if (defined EAM || defined ADP || defined MEAM) && !defined NOPUNISH
  /* write EAM punishments */
  if (g_param.write_output_files)
  {
    strcpy(file, g_files.output_prefix);
    strcat(file, ".punish");
    outfile = fopen(file, "w");
    if (NULL == outfile)
      error(1, "Could not open file %s\n", file);
    fprintf(outfile, "Limiting constraints\n");
    fprintf(outfile, "#conf\tp^2\t\tpunishment\n");
  }
  else
  {
    outfile = stdout;
    printf("Punishment Constraints\n");
  }
  for (i = g_calc.limit_p; i < g_calc.dummy_p; i++)
  {
    sqr = dsquare(force[i]);
    if (g_param.write_output_files)
      fprintf(outfile, "%d\t%f\t%f\n", i - g_calc.limit_p, sqr,
              force[i] + g_config.force_0[i]);
    else
      fprintf(outfile, "%d %f %f %f %f\n", i - g_calc.limit_p, sqr,
              force[i] + g_config.force_0[i], g_config.force_0[i],
              force[i] / g_config.force_0[i]);
  }
  if (g_param.write_output_files)
  {
    fprintf(outfile, "\nDummy Constraints\n");
    fprintf(outfile, "element\tU^2\t\tU'^2\t\tU\t\tU'\n");
    for (i = g_calc.dummy_p; i < g_calc.dummy_p + g_param.ntypes; i++)
    {
#ifndef RESCALE
      sqr = dsquare(force[i]);
      fprintf(outfile, "%s\t%f\t%f\t%f\t%g\n", g_config.elements[i - g_calc.dummy_p], 0.0,
              sqr, 0.0, force[i]);
#else
      sqr = dsquare(force[i + g_param.ntypes]);
      fprintf(outfile, "%s\t%f\t%f\t%f\t%f\n", g_config.elements[i - g_calc.dummy_p], sqr,
              dsquare(force[i]), force[i + g_param.ntypes], force[i]);
#endif /* !RESCALE */
    }
#ifndef RESCALE
    fprintf(outfile, "\nNORESCALE: <n>!=1\n");
    fprintf(outfile, "<n>=%f\n",
            force[g_calc.dummy_p + g_param.ntypes] / DUMMY_WEIGHT + 1);
    fprintf(outfile, "Additional punishment of %f added.\n",
            dsquare(force[g_calc.dummy_p + g_param.ntypes]));
#endif /* !RESCALE */
    printf("Punishment constraints data written to \t%s\n", file);
    fclose(outfile);
  }
  else
  {
    fprintf(outfile, "Dummy Constraints\n");
    for (i = g_calc.dummy_p; i < g_calc.dummy_p + g_param.ntypes; i++)
    {
      sqr = dsquare(force[i]);
      fprintf(outfile, "%s\t%f\t%f\n", g_config.elements[i - g_calc.dummy_p], sqr,
              force[i]);
    }
  }
#endif /* (EAM || ADP || MEAM) && !NOPUNISH */

  /* final error report */
  printf("\n###### error report ######\n");
  if (g_param.write_output_files)
  {
    strcpy(file, g_files.output_prefix);
    strcat(file, ".error");
    outfile = fopen(file, "w");
    if (NULL == outfile)
      error(1, "Could not open file %s\n", file);
#ifndef STRESS
    fprintf(outfile, "total error sum %f, count %d (%d forces, %d energies)\n", tot,
            g_calc.mdim - 6 * g_config.nconf, 3 * g_config.natoms, g_config.nconf);
#else
    fprintf(outfile,
            "total error sum %f, count %d (%d forces, %d energies, %d stresses)\n", tot,
            g_calc.mdim, 3 * g_config.natoms, g_config.nconf, 6 * g_config.nconf);
#endif /* !STRESS */
  }
#ifndef STRESS
  printf("total error sum %f, count %d (%d forces, %d energies)\n", tot,
         g_calc.mdim - 6 * g_config.nconf, 3 * g_config.natoms, g_config.nconf);
#else
  printf("total error sum %f, count %d (%d forces, %d energies, %d stresses)\n", tot,
         g_calc.mdim, 3 * g_config.natoms, g_config.nconf, 6 * g_config.nconf);
#endif /* !STRESS */

  /* calculate the rms errors for forces, energies, stress */
  rms[0] = 0.0; /* rms rms for forces */
  rms[1] = 0.0; /* energies */
  rms[2] = 0.0; /* stresses */

  /* forces */
  for (i = 0; i < g_config.natoms; i++)
  {
#ifdef CONTRIB
    if (atoms[i].contrib)
    {
      contrib_atoms++;
#endif /* CONTRIB */
#if defined(FWEIGHT)
      rms[0] += dsquare(force[3 * i + 0] * (FORCE_EPS + g_config.atoms[i].absforce));
      rms[0] += dsquare(force[3 * i + 1] * (FORCE_EPS + g_config.atoms[i].absforce));
      rms[0] += dsquare(force[3 * i + 2] * (FORCE_EPS + g_config.atoms[i].absforce));
#else
    rms[0] += dsquare(force[3 * i + 0]);
    rms[0] += dsquare(force[3 * i + 1]);
    rms[0] += dsquare(force[3 * i + 2]);
#endif
#ifdef CONTRIB
    }
#endif /* CONTRIB */
  }
#ifdef CONTRIB
  rms[0] = sqrt(rms[0] / (3 * contrib_atoms));
#else
  rms[0] = sqrt(rms[0] / (3 * g_config.natoms));
#endif /* CONTRIB */

  /* energies */
  if (g_param.eweight != 0)
  {
    for (i = 0; i < g_config.nconf; i++)
      rms[1] += dsquare(force[3 * g_config.natoms + i]);
    if (isnan(rms[1]))
      rms[1] = 0.0;
    rms[1] = sqrt(rms[1] / g_config.nconf);
  }

  /* stresses */
  if (g_param.sweight != 0)
  {
    for (i = 0; i < g_config.nconf; i++)
      for (j = 0; j < 6; j++)
        rms[2] += dsquare(force[3 * g_config.natoms + g_config.nconf + 6 * i + j]);
    if (isnan(rms[2]))
      rms[2] = 0.0;
    rms[2] = sqrt(rms[2] / (6 * g_config.nconf));
  }

  if (g_param.write_output_files)
  {
    fprintf(outfile, "sum of force-errors  = %e\t\t( %.3f%% - av: %f)\n", f_sum,
            f_sum / tot * 100, f_sum / (3 * g_config.natoms));
    if (g_param.eweight != 0)
      fprintf(outfile, "sum of energy-errors = %e\t\t( %.3f%% )\n", e_sum,
              e_sum / tot * 100);
#ifdef STRESS
    if (g_param.sweight != 0)
      fprintf(outfile, "sum of stress-errors = %e\t\t( %.3f%% )\n", s_sum,
              s_sum / tot * 100);
#endif /* STRESS */
    if ((tot - f_sum - e_sum - s_sum) > 0.01 && g_param.opt == 1)
    {
      fprintf(outfile,
              "\n --> Warning <--\nThis sum contains punishments! Please check "
              "your results.\n");
      fprintf(outfile, "Total sum of punishments = %e\t\t( %.3f%% )\n\n",
              tot - f_sum - e_sum - s_sum, (tot - f_sum - e_sum - s_sum) / tot * 100);
    }
    fprintf(outfile, "rms-errors:\n");
    fprintf(outfile, "force \t%e\t(%f meV/A)\n", rms[0], rms[0] * 1000);
    if (g_param.eweight != 0)
      fprintf(outfile, "energy \t%e\t(%f meV)\n", rms[1], rms[1] * 1000);
#ifdef STRESS
    if (g_param.sweight != 0)
      fprintf(outfile, "stress \t%e\t(%f MPa)\n", rms[2], rms[2] / 160.2 * 1000);
#endif /* STRESS */
    fprintf(outfile, "\n");
    fprintf(outfile, "\tforce [meV/A]\tenergy [meV]\tstress [MPa]\terror sum\n");
    fprintf(outfile, "RMS:\t%e\t%e\t%e\t%f\n", rms[0] * 1000, rms[1] * 1000,
            rms[2] / 160.2 * 1000, tot);
  }
  printf("sum of force-errors  = %e\t\t( %.3f%% - av: %e)\n", f_sum, f_sum / tot * 100,
         f_sum / (3 * g_config.natoms));
  if (g_param.eweight != 0)
    printf("sum of energy-errors = %e\t\t( %.3f%% )\n", e_sum, e_sum / tot * 100);
#ifdef STRESS
  if (g_param.sweight != 0)
    printf("sum of stress-errors = %e\t\t( %.3f%% )\n", s_sum, s_sum / tot * 100);
#endif /* STRESS */
  if ((tot - f_sum - e_sum - s_sum) > 0.01 && g_param.opt == 1)
  {
    printf(
        "\n --> Warning <--\nThis sum contains punishments! Check your "
        "results.\n");
    printf("sum of punishments = %e\t\t( %.3f%% )\n\n", tot - f_sum - e_sum - s_sum,
           (tot - f_sum - e_sum - s_sum) / tot * 100);
  }
  printf("rms-errors:\n");
  printf("force \t%e\t(%f meV/A)\n", rms[0], rms[0] * 1000);
  if (g_param.eweight != 0)
    printf("energy \t%e\t(%f meV)\n", rms[1], rms[1] * 1000);
#ifdef STRESS
  if (g_param.sweight != 0)
    printf("stress \t%e\t(%f MPa)\n", rms[2], rms[2] / 160.2 * 1000);
#endif /* STRESS */
  if (g_param.write_output_files)
  {
    fclose(outfile);
  }

  for (i = 0; i < 6; i++)
    free(component[i]);
}
