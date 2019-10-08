/****************************************************************
 *
 * uq.c: Potential Ensemble method for Uncertainty Quantification
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
 *****************************************************************/

#include "potfit.h"

#if defined(MKL)
#include <mkl_lapack.h>
#elif defined(__ACCELERATE__)
#include <Accelerate/Accelerate.h>
#else
#error No math library defined!
#endif  // MKL

#include "errors.h"
#include "force.h"
#include "memory.h"
#include "potential_input.h"
#include "potential_output.h"
#include "random.h"
#include "uq.h"
#include "utils.h"

#define VERY_SMALL 1.E-12
#define VERY_LARGE 1.E12

// Only for analytic potentials at the moment

#if defined(UQ) && defined(APOT)

double single_param_pert_cost(double, int);
void subsection_pert(double*, double*, int, double);
double hess_bracketing(double, int, int);
int calc_hessian(double**, double*);
int calc_h0_eigenvectors(double**, double, double, double**, double*);
int calc_svd(double**, double**, double*);
double generate_mc_sample(double** const, double** const, double, double, double*, int*, FILE*);
int mc_moves(double**,double*, double*, double, FILE*);

/****************************************************************
 *
 *   main potential ensemble routine
 *
 ****************************************************************/

void ensemble_generation(double cost_0)
{
  /* open file */
  FILE* outfile = fopen(g_files.ensemblefile, "w");
  if (outfile == NULL)
    error(1, "Could not open file %s for writing\n", g_files.ensemblefile);

  /* Initialise variables to 0 */
  int pot_attempts = 0;
  double acc_prob = 0.0;

  /* Default input parameters if not specified */
  if (g_param.eig_max == 0) {
    g_param.eig_max = 1;
  }

  /* Calculate the best fit hessian */
  double** hessian = mat_double(g_pot.opt_pot.idxlen, g_pot.opt_pot.idxlen);

  int hess_counter = 1;
  int flag = calc_hessian(hessian, &cost_0);
  while (flag == 1){
    /* Check that we haven't been through this thing 10 times, if so error */
    if (hess_counter == 10){
          error(1,
          "Too many recalculations of the hessian implies the potential is "
          "poorly fit.\nIt is advised to rerun parameter optimisation and use "
          "the true minimum.\n");
    }

    hess_counter += 1;
    flag = calc_hessian(hessian, &cost_0);

  }

  printf("Hessian calulated, finding it's eigenvalues.\n");
  fflush(stdout);

  /* Print eigenvalues and eigenvectors of hessian to ensemblefile */
  fprintf(outfile, "# hessian:\n# ");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      fprintf(outfile, "%-10.8f ", hessian[i][j]);
    }
    fprintf(outfile, "\n# ");
  }

  int m = 0;
  // Initial lower bound for eigenvalues - this range is adjusted until all are found
  double vl = -1.0;
  // Initial upper bound for eigenvalues
  double vu = 10000.0;
  int count = 0;

  /* Allocate memory for the eigenvectors */
  double** v_0 = mat_double(g_pot.opt_pot.idxlen, g_pot.opt_pot.idxlen);
  double eigenvalues[g_pot.opt_pot.idxlen];

  /* Enter loop to find eigenvalues, increasing the range by a factor of ten
     each until all eigenvalues are found. If after 10 iterations the
     eigenvalues have not be found, revert to signular value decomposition */
  while (m < g_pot.opt_pot.idxlen) {
    vl *= 10;
    vu *= 10;

    if (g_param.use_svd == 1) {
      m = calc_svd(hessian, v_0, eigenvalues);
    } else {
      m = calc_h0_eigenvectors(hessian, vl, vu, v_0, eigenvalues);

      if (count > 10) {
        printf("WARNING: NOT CONVERGING! Use singular value decomposition.\n");

        m = calc_svd(hessian, v_0, eigenvalues);
        fflush(stdout);

        /* Check all eigenvalues have been found, otherwise exit */
        if (m != g_pot.opt_pot.idxlen) {
          error(1, "Only %d of %d eigenvalues were found!", m,
                g_pot.opt_pot.idxlen);
        }
      }
      count += 1;
    }
  }

  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    if (eigenvalues[i] < 0) {
      printf(
          "WARNING: Eigenvalue %d is negative = %.4f, has the best fit minimum "
          "been found?!\nThis implies the best fit potential is not at the "
          "miminum!\n",
          i, eigenvalues[i]);
    }
  }

  /* Print eigenvalues and eigenvectors of hessian to ensemblefile */
  fprintf(outfile, "\n# eigenvalues: (l1, l2, ..., lN)\n# ");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    fprintf(outfile, "%-11.8f ", eigenvalues[i]);
  }
  fprintf(outfile, "\n#\n");
  fprintf(outfile, "# eigenvectors: (v1, v2, ..., vN)\n# ");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      fprintf(outfile, "%-10.8f ", v_0[j][i]);  // Rows contain eigenvectors
    }
    fprintf(outfile, "\n# ");
  }

  /* Print initial best fit */
  fprintf(outfile, "\n# id      ");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    fprintf(outfile, "param_%-4d ", i + 1);
  }
  fprintf(outfile, "cost       weight     accepted   attempts   acc_prob\n");
  fflush(outfile);

  /* Write initial cost to file */
  fprintf(outfile, "%-10d", pot_attempts);
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    fprintf(outfile, "%-10.8lf ", g_pot.opt_pot.table[g_pot.opt_pot.idx[i]]);
  }
  fprintf(outfile, "%.8lf ", cost_0);  // cost_temp);

  printf("Beginning MCMC ensemble generation.\n");
  fflush(stdout);

  /* Initialise variables and take first Monte Carlo step */
  int weight = 1;
  double cost_attempt = cost_0;

  /* run until number of moves specified in param file are accepted */
  for (int i = 0; i <= g_param.acc_moves; i++) {
    if (i == g_param.acc_moves) {
      pot_attempts += weight;
      acc_prob =
          (((double)i + 1.0)) /
          (double)
              pot_attempts; /* Add one to include the MC step outside loop */

      /* For the final configuration (i = (g_param.acc_moves - 1) print the
       * remaining weight calculated and then exit */
      fprintf(outfile, "%-10d 1          %-10d %-10.2f\n", weight, pot_attempts,
              acc_prob);
      continue;
    }

    double cost = generate_mc_sample(hessian, v_0, cost_attempt, cost_0,
                                     eigenvalues, &weight, outfile);

    cost_attempt = cost;

    pot_attempts += weight;
    acc_prob =
        (((double)i + 1.0)) /
        (double)pot_attempts; /* Add one to include the MC step outside loop */

    /* Write weight from best cost set - i.e. how many trials before this new
     * parameter set was accepted */
    fprintf(outfile, "%-10d 1          %-10d %-10.2f\n", weight, pot_attempts,
            acc_prob);

    /* Write accepted move to file */
    fprintf(outfile, "%-10d", i + 1);
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      fprintf(outfile, "%-10.8lf ", g_pot.opt_pot.table[g_pot.opt_pot.idx[j]]);
    }
    fprintf(outfile, "%.8lf ", cost);

    char file[255];
    sprintf(file, "%s.ensemble_pot_%d", g_files.output_prefix, i + 1);

    if (cost < cost_0) {
      update_apot_table(g_pot.opt_pot.table);
      write_pot_table_potfit(file);

      printf(
          "WARNING: New best fit parameter set found for potential %s. Old "
          "cost = %.8lf, new cost = %.8lf\n",
          file, cost_0, cost);
    }
    /* Write potential input files for parameter ensemble */
    else if ((g_param.write_ensemble != 0) &&
             ((i + 1) % g_param.write_ensemble == 0)) {
      update_apot_table(g_pot.opt_pot.table);
      write_pot_table_potfit(file);
    }
  }

  fclose(outfile);
  printf("UQ ensemble parameters written to %s\n", g_files.ensemblefile);
}

/****************************************************************
 *
 *    Bracketing function for hessian finite difference
 *    perturbation range - finds perturbation value
 *
 ****************************************************************/

double hess_bracketing(double cost_aim, int index, int dir)
{
  /* SET INITIAL PERT ORDER OF MAGNITUDE */

  double ub = 0.0; /* Upper bound to be set */
  double pert = 0.0;

  if (dir == 1) {
    pert = -0.000001; /* Initial pert guess */
#if defined(DEBUG)
    printf("Calculating negative perturbation direction.\n");
#endif
  } else {
    pert = 0.000001; /* Initial pert guess */
#if defined(DEBUG)
    printf("Calculating positive perturbation direction.\n");
#endif
  }

  double lb = 0.0; /* Lower bound to be set */
  double pert_cost =
      single_param_pert_cost(pert, index); /* Set initial pert cost */
  double lb_cost_r = 0.0;
  double pert_range = 0.0;

  /* If the initial guess is too large reduce the perturbation to set lb */
  while (pert_cost > cost_aim) {
    ub = pert;
    pert /= 10;
    pert_cost = single_param_pert_cost(pert, index);
  }
  lb = pert;
#if defined(DEBUG)
  printf("INITIAL lB IS %g, pert_cost = %.6g, cost_aim = %.6g\n", lb, pert_cost,
         cost_aim);
#endif

  /* Set upper bound order of magnitude */
  while (pert_cost < cost_aim) {
    /* Fill as cost(lb) */
    lb_cost_r = pert_cost;

    lb = pert;
#if defined(DEBUG)
    printf("INITIAL lB IS %g, pert_cost = %.6g, cost_aim = %.6g\n", lb,
           pert_cost, cost_aim);
#endif
    pert *= 10;
    pert_cost = single_param_pert_cost(pert, index);
  }
  ub = pert;

#if defined(DEBUG)
  printf("INITIAL UB IS %g, pert_cost = %.6g, cost_aim = %.6g\n", ub, pert_cost,
         cost_aim);
#endif

  /* Equivalent to cost(ub) - cost(lb) */
  pert_range = pert_cost - lb_cost_r;

#if defined(DEBUG)
  printf("pert range = %g\n", pert_range);
#endif

  /* If the bounds are somehow muddled up, swap them */
  if (fabs(ub) < fabs(lb)) {
    double temp;
    SWAP(lb, ub, temp);
  }

  /* SUBDIVIDE INTERVAL INTO 10 AND EVALUATE SECTIONS */
  /* REPEAT UNTIL RANGE IS WITHIN 5% OF COST_AIM VALUE */
  subsection_pert(&lb, &ub, index, cost_aim);

#if defined(DEBUG)
  printf("FINAL RANGE = [%g,%g]\n", lb, ub);
  fflush(stdout);
#endif

  /* FIT A LINE TO THE RANGE TO FIND THE CORRESPONDING PERT VALUE */
  /* Join lb and ub by a line, use the gradient to calculate */
  /* pert value (i.e. x) corresponding to cost_aim. */

  /* Calculate the costs for the bounds */
  double lb_cost = single_param_pert_cost(lb, index);
  double ub_cost = single_param_pert_cost(ub, index);

  double grad = (ub_cost - lb_cost) / (ub - lb);
  double lb_pert = ((cost_aim - lb_cost) / grad) + lb;
  double ub_pert = ((cost_aim - ub_cost) / grad) + ub;

  /* Take pert as the average of the two values, although they should be the
   * same */
  /* Always return the absolute value of perturbation so as not to confuse
   * finite difference algorithm */
  double perturbation = fabs(lb_pert + ub_pert) / 2.0;

#if defined(DEBUG)
  printf("grad = %g lb_pert = %g ub_pert = %g pert[%d] = %g\n", grad, lb_pert,
         ub_pert, index, perturbation);
#endif

  return perturbation;
}

/****************************************************************
 *
 *    Subdivide perturbations bracketing cost_T
 *    to find smaller range.
 *
 ****************************************************************/

void subsection_pert(double* lb, double* ub, int index, double cost_aim)
{
  /* SUBDIVIDE INTERVAL INTO 10 AND EVALUATE SECTIONS */
  double range[11];
  double cost_range[11];
  double percentage_aim;

  /* Enter ub, lb and their costs as already known */
  range[0] = *lb;
  range[10] = *ub;
  cost_range[0] = single_param_pert_cost(*lb, index);
  cost_range[10] = single_param_pert_cost(*ub, index);

  /* Check the costs bracket the cost_aim as expected, intuitive (righthand,
   * positive) x-direction is pos = 1.*/
  if ((cost_range[0] > cost_aim) || (cost_range[10] < cost_aim)) {
    /* This could be because it's the negative x-direction which is steeper, so
     * check */
    *lb *= (-1);
    *ub *= (-1);
    cost_range[0] = single_param_pert_cost(*lb, index);
    cost_range[10] = single_param_pert_cost(*ub, index);

    /* If changing pert signs fixed it, continue with new values */
    /* However if the same thing happens, something else is wrong - abort. */
    if ((cost_range[0] > cost_aim) || (cost_range[10] < cost_aim)) {
      printf(
          "WARNING: Cost_aim = %g, is not bracketed in range [%g,%g],\n "
          "Cost_lb = %g, cost_ub = %g\n",
          cost_aim, *lb, *ub, cost_range[0], cost_range[10]);
      error(1, "Bracketing of parameter perturbation failed.\n");
    }
  }

#if defined(DEBUG)
  printf("subsection: [%g, %g] = cost[%g,%g], cost_aim = %g\n", range[0],
         range[10], cost_range[0], cost_range[10], cost_aim);
#endif

  /* Fill remaining range with values and costs */
  for (int i = 1; i < 11; i++) {
    range[i] = range[0] + (i * 0.1 * (range[10] - range[0]));
    cost_range[i] = single_param_pert_cost(range[i], index);

    /* When cost becomes greater than cost aim, break */
    if (cost_range[i] >= cost_aim) {
#if defined(DEBUG)
      printf(
          "Range is lb = %g ub = %g, cost_range = [%g, %g], percentage of "
          "cost_aim = %g%%\n",
          range[i - 1], range[i], cost_range[i - 1], cost_range[i],
          (100 * (cost_range[i] - cost_range[i - 1])) / cost_aim);
      fflush(stdout);
#endif

      *lb = range[i - 1];
      *ub = range[i];
      break;
    } else if (i == 10) {
      /* Else the change is in the last section */

#if defined(DEBUG)
      printf(
          "FINAL SECTION: Range is lb = %g ub = %g, cost_range = [%g, %g], "
          "percentage of cost_aim = %g%%\n",
          range[i], range[i + 1], cost_range[i - 1], cost_range[i],
          (100 * (cost_range[i] - cost_range[i - 1])) / cost_aim);
      fflush(stdout);
#endif

      *lb = range[i];
    }
  }

  /* ITERATE TO FIND SECTION BOUNDS TO WITHIN 5% OF COST_AIM */
  percentage_aim = (cost_range[10] - cost_range[0]) / cost_aim;

#if defined(DEBUG)
  printf("percentage_aim = %g\n", percentage_aim);
#endif

  /* If the range is within 5% of the cost_aim value, return and move to
   * checking the negative x-direction */
  if (percentage_aim <= 0.05) {
    return;
  } else {
    /* Is range is still too large, subdivide into 10 again and repeat */
    subsection_pert(lb, ub, index, cost_aim);
  }
}

/****************************************************************
 *
 *    Perturb parameter[index] by perturbation
 *
 ****************************************************************/

double single_param_pert_cost(double pert, int index)
{
  double old_value = g_pot.opt_pot.table[g_pot.opt_pot.idx[index]];

  // Perturb parameter by perturbation
  g_pot.opt_pot.table[g_pot.opt_pot.idx[index]] *= (1.0 + pert);

  // Calculate the lb perturbation cost
  const double cost = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

  // Return to best fit params
  g_pot.opt_pot.table[g_pot.opt_pot.idx[index]] = old_value;

  return cost;
}

/****************************************************************
 *
 *    Calculate the best fit potential hessian
 *
 ****************************************************************/

int calc_hessian(double** hessian, double* cost)
{
  /*  Implementing equation 5.7.10 from Numerical recipes in C

    Create the Hessian of analytic potential parameters
    For N parameters, require:
    diagonal: 2N cost evaluations
    off-diagonal: 2N(N-1) cost evaluations (4 per hessian element)

    Allocate memory to store:
    - the cost evaluations per hessian element
    - the size of each parameter perturbation (i.e. 0.0001*parameter)
    - the final hessian elements */
  double param_perturb_dist[g_pot.opt_pot.idxlen];
  double two_cost = 2.0 * *cost;
  double new_cost_param_values[g_pot.opt_pot.idxlen + 1];



  /* Initialise values for possible better fit found */
  for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
    new_cost_param_values[j] = 0.0;
  }
  new_cost_param_values[g_pot.opt_pot.idxlen] = VERY_LARGE;

  /* FIND PERTURBATION VALUES FOR HESSIAN CURVATURE CALCULATION  */

  double cost_aim =
      *cost + ((2.0 * *cost * g_param.uq_temp) / g_pot.opt_pot.idxlen);
  double pert[g_pot.opt_pot.idxlen];

  /* Find the correct perturbation value for each parameter */
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    /* If user specified perturbation, use this instead. */
    if (g_param.hess_pert < 0) {
      printf("Cost aim = %g\n", cost_aim);
      fflush(stdout);
      /* Take the minimum perturbation of each direction */
      double pert_pos = hess_bracketing(cost_aim, i, 0);
      double pert_neg = hess_bracketing(cost_aim, i, 1);
      printf("PERT VALUES: pos_pert = %g, neg_pert = %g diff = %g\n", pert_pos,
             pert_neg, fabs(pert_pos - pert_neg));
      if (pert_pos <= pert_neg) {
        pert[i] = pert_pos;
        printf("Using positive perturbtion.\n");
      } else {
        pert[i] = pert_neg;
        printf("Using negative perturbtion.\n");
      }

      // Dont allow perturbations more than tha value of the parameter
      if (pert[i] > 1.0) {
        pert[i] = 1.0;
      }

      printf(
          "FINAL PERT VALUE %.8lf for param %d = %g (percentage of param = "
          "%g%%)\n",
          pert[i], i, g_pot.opt_pot.table[g_pot.opt_pot.idx[i]],
          fabs((pert[i] * 100) / g_pot.opt_pot.table[g_pot.opt_pot.idx[i]]));

    } else {
      pert[i] = g_param.hess_pert;
    }
  } /* parameter loop */

  /* Pre-calculate each parameter perturbation for diagonal entries */
  for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
    param_perturb_dist[j] = pert[j] * g_pot.opt_pot.table[g_pot.opt_pot.idx[j]];

    // THIS SHOULD BE HIGHER
    if (g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] == 0) {
      param_perturb_dist[j] = pert[j];
      printf("parameter %d is 0. Using set perturbation of %f.\n", j, pert[j]);
    }
  }

#if defined(DEBUG)
  printf("i i cost_plus cost_minus diff\n");
#endif

  /* For diagonal entries, use (c_(i+1) - 2*cost +
   * c_(i-1))/(param_perturb_dist[i]^2) */
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    double cost_plus;
    double cost_minus;

    g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] += param_perturb_dist[i];
    cost_plus = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

    if ((cost_plus < *cost) &&
        (cost_plus < new_cost_param_values[g_pot.opt_pot.idxlen])) {
      /* If new minima is found, store these values */
      for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
        new_cost_param_values[j] = g_pot.opt_pot.table[g_pot.opt_pot.idx[j]];
      }
      new_cost_param_values[g_pot.opt_pot.idxlen] = cost_plus;
    }

    g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] -= 2 * param_perturb_dist[i];
    cost_minus = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

    if ((cost_minus < *cost) &&
        (cost_minus < new_cost_param_values[g_pot.opt_pot.idxlen])) {
      /* If new minima is found, store these values */
      for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
        new_cost_param_values[j] = g_pot.opt_pot.table[g_pot.opt_pot.idx[j]];
      }
      new_cost_param_values[g_pot.opt_pot.idxlen] = cost_minus;
    }

    /* Reset original param values without perturbation */
    g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] += param_perturb_dist[i];

    /* print a warning if either cost_plus or cost_minus are less than 10^(-12)
     * or a new minima is found */
    if ((cost_plus < VERY_SMALL) || (cost_minus < VERY_SMALL)) {
      printf(
          "WARNING: The change in cost_plus/cost_minus when calculating the "
          "hessian is less than 10^(-12).\n This will affect precision. "
          "Consider changing the scale of cost perturbation. \n");
    }

    hessian[i][i] = cost_plus - two_cost + cost_minus;
    hessian[i][i] /= (param_perturb_dist[i] * param_perturb_dist[i]);

#if defined(DEBUG)
    if ((cost_plus > cost_aim) || (cost_minus > cost_aim)) {
      printf("%d %d %.3f %.3f * %.3f\n", i, i, cost_plus, cost_minus,
             fabs(cost_plus - cost_minus));
    } else {
      printf("%d %d %.3f %.3f   %.3f\n", i, i, cost_plus, cost_minus,
             fabs(cost_plus - cost_minus));
    }
    fflush(stdout);
#endif
  }

  /* For off-diagonal entries:
     Use
     [c_(i+1)(j+1)-c_(i+1)(j-1)-c_(i-1)(j+1)+c_(i-1)(j-1)]/(param_perturb_dist[i]*param_perturb_dist[j]*4)
   */
#if defined(DEBUG)
  printf(
      "i j cost_2plus cost_2minus cost_pm cost_mp cost_2diff "
      "cost_pm_mp_diff\n");
#endif

  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    for (int j = (i + 1); j < g_pot.opt_pot.idxlen; j++) {
      double cost_2plus;
      double cost_2minus;
      double cost_pm;
      double cost_mp;

      /* c_(i+1)(j+1) */
      g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] += param_perturb_dist[i];
      g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] += param_perturb_dist[j];
      cost_2plus = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

      /* c_(i+1)(j-1) */
      g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] -= 2 * param_perturb_dist[j];
      cost_pm = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

      if ((cost_pm < *cost) &&
          (cost_pm < new_cost_param_values[g_pot.opt_pot.idxlen])) {
        /* If new minima is found, store these values */
        for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
          new_cost_param_values[j] = g_pot.opt_pot.table[g_pot.opt_pot.idx[j]];
        }
        new_cost_param_values[g_pot.opt_pot.idxlen] = cost_pm;
      }

      /* c_(i-1)(j+1) */
      g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] -= 2 * param_perturb_dist[i];
      g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] += 2 * param_perturb_dist[j];
      cost_mp = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

      if ((cost_mp < *cost) &&
          (cost_mp < new_cost_param_values[g_pot.opt_pot.idxlen])) {
        /* If new minima is found, store these values */
        for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
          new_cost_param_values[j] = g_pot.opt_pot.table[g_pot.opt_pot.idx[j]];
        }
        new_cost_param_values[g_pot.opt_pot.idxlen] = cost_mp;
      }

      /* c_(i-1)(j-1) */
      g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] -= 2 * param_perturb_dist[j];
      cost_2minus = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

      /* reset to start */
      g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] += param_perturb_dist[i];
      g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] += param_perturb_dist[j];

      /* print a warning if either cost_pm or cost_mp are less than 10^(-12) or
       * a new minima is found */
      if ((cost_pm < VERY_SMALL) || (cost_mp < VERY_SMALL)) {
        printf(
            "WARNING: The change in cost_pm/cost_mp when calculating the "
            "hessian is less than 10^(-12).\nThis will affect precision. "
            "Consider changing the scale of cost perturbation. \n");
      }

      hessian[i][j] = cost_2plus + cost_2minus - cost_pm - cost_mp;
      hessian[i][j] /= (4 * param_perturb_dist[i] * param_perturb_dist[j]);

      hessian[j][i] = hessian[i][j];

#if defined(DEBUG)
      if ((cost_2plus > cost_aim) || (cost_2minus > cost_aim) ||
          (cost_pm > cost_aim) || (cost_mp > cost_aim)) {
        printf("%d %d %.3f %.3f %.3f %.3f * %.3f %.3f\n", i, j, cost_2plus,
               cost_2minus, cost_pm, cost_mp, fabs(cost_2plus - cost_2minus),
               fabs(cost_pm - cost_mp));
      } else {
        printf("%d %d %.3f %.3f %.3f %.3f   %.3f %.3f\n", i, j, cost_2plus,
               cost_2minus, cost_pm, cost_mp, fabs(cost_2plus - cost_2minus),
               fabs(cost_pm - cost_mp));
      }
      fflush(stdout);
#endif
    }
  }

  /* If new cost value is found, return parameters */
  if (new_cost_param_values[g_pot.opt_pot.idxlen] != VERY_LARGE) {
    printf(
        "WARNING: A new cost minimum has been found.\nOriginal cost = %f,\t "
        "New cost = %f.\nCalculation restarting with new best fit potential "
        "values.\n\n",
        *cost, new_cost_param_values[g_pot.opt_pot.idxlen]);

    printf("NEW COST MINIMA VALUES:\n");
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      printf("Param %d = %.8lf\n", j, new_cost_param_values[j]);
    }
    printf("Cost = %f\n", new_cost_param_values[g_pot.opt_pot.idxlen]);
    fflush(stdout);

    /* Move old potential to temp file */
    if (g_files.tempfile && strlen(g_files.tempfile)) {
#if defined(APOT)
      update_apot_table(g_pot.opt_pot.table);
#endif  // APOT
      write_pot_table_potfit(g_files.tempfile);
    }

    /* Set new cost potential parameters as the best fit potential  */
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      g_pot.opt_pot.table[g_pot.opt_pot.idx[j]] = new_cost_param_values[j];
    }

    /* Write out new end potential */
#if defined(APOT)
    update_apot_table(g_pot.opt_pot.table);
#endif  // APOT
    write_pot_table_potfit(g_files.endpot);

    // will not work with MPI
#if defined(PDIST) && !defined(MPI)
    write_pairdist(&g_pot.opt_pot, g_files.distfile);
#endif  // PDIST && !MPI

    /* write the error files for forces, energies, stresses, ... */
    write_errors(g_calc.force, new_cost_param_values[g_pot.opt_pot.idxlen]);
    fflush(stdout);
    /* Update new cost minima value and return flag of 1 to recalculate the Hessian */
    *cost = new_cost_param_values[g_pot.opt_pot.idxlen];
    return 1;
  } /* If a new cost is found */

  return 0;
}

/***********************************************************************
 *
 *    Find hessian eigenvalues and eigenvectors by eigen decomposition
 *
 ***********************************************************************/

int calc_h0_eigenvectors(double** hessian, double lower_bound,
                         double upper_bound, double** v_0, double* w)
{
  char jobz = 'V';   /* Compute eigenvectors and eigenvalues */
  char range = 'A';  //'V'; /* all eigenvalues in the half-open interval (VL,VU]
                     //will be found */
  char uplo = 'U';   /* Upper triangle of A is stored */
  int lda = g_pot.opt_pot
                .idxlen; /* leading dimension of the array A. lda >= max(1,N) */
  const char inp = 'S';
  /* 2*DLAMCH('S');  absolute error tolerance for eigenvalues */
  double abstol = 2 * DLAMCH( &inp);
  int ldz = g_pot.opt_pot.idxlen; /* Dimension of array z */
  int il = 0;
  int iu = 0;

  int m; /* number eigenvalues found */
  int iwork[5 * g_pot.opt_pot.idxlen];
  int lwork = 8 * g_pot.opt_pot.idxlen;
  double work[lwork];
  int ifail[g_pot.opt_pot.idxlen]; /* contains indices of unconverged
                                      eigenvectors if info > 0  */
  int info = 0;

#if defined(MKL)
  dsyevx_(&jobz, &range, &uplo, &g_pot.opt_pot.idxlen, &hessian[0][0], &lda,
          &lower_bound, &upper_bound, &il, &iu, &abstol, &m, w, &v_0[0][0],
          &ldz, work, &lwork, iwork, ifail, &info);
#elif defined(__ACCELERATE__)
  dsyevx_(&jobz, &range, &uplo, &g_pot.opt_pot.idxlen, &hessian[0][0], &lda,
          &lower_bound, &upper_bound, &il, &iu, &abstol, &m, w, &v_0[0][0],
          &ldz, work, &lwork, iwork, ifail, &info);
#endif

  if (info == 0) {
    printf("Eigenvalue calculation completed successfully.\n");
  } else if (info > 0) {
    printf("%d eigenvectors failed to converge.\n", info);
    m = g_pot.opt_pot.idxlen - m;
  } else {
    error(1,
          "Illegal argument supplied to eigenvalue decomposition function.\n");
  }
  fflush(stdout);

  return m;
}

/********************************************************
 *
 *    Find hessian eigenvalues and eigenvectors by SVD
 *
 *********************************************************/

int calc_svd(double** hessian, double** u, double* s)
{
  char jobu = 'A';  /* Compute left singular vectors */
  char jobvt = 'N'; /* Compute right singular vectors */
  int lda = g_pot.opt_pot
                .idxlen; /* leading dimension of the array A. lda >= max(1,N) */

  int lwork = 5 * g_pot.opt_pot.idxlen;
  double work[lwork];
  int info = 0;
  double vt[g_pot.opt_pot.idxlen][g_pot.opt_pot.idxlen];

#if defined(MKL)
  dgesvd_(&jobu, &jobvt, &g_pot.opt_pot.idxlen, &g_pot.opt_pot.idxlen,
          &hessian[0][0], &lda, s, &u[0][0], &lda, &vt[0][0],
          &g_pot.opt_pot.idxlen, work, &lwork, &info);
#elif defined(__ACCELERATE__)
  dgesvd_(&jobu, &jobvt, &g_pot.opt_pot.idxlen, &g_pot.opt_pot.idxlen,
          &hessian[0][0], &lda, s, &u[0][0], &lda, &vt[0][0],
          &g_pot.opt_pot.idxlen, work, &lwork, &info);
#endif

#if defined(DEBUG)
  printf("------------------------------------------------------\n\n");
  printf(
      "Singular values and singular vectors of the best fit hessian for U:\n");
  printf("eigenvalue, eigenvector (x, y ,z)\n");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    printf("%-11.12f ", s[i]);
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      printf("%-10.12f ", u[i][j]);  // Rows contain singularvectors
    }
    printf("\n");
  }

  printf("\n------------------------------------------------------\n\n");

  printf("------------------------------------------------------\n\n");
  printf(
      "Singular values and singular vectors of the best fit hessian for "
      "V^T:\n");
  printf("eigenvalue, eigenvector (x, y ,z)\n");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    printf("%-11.12f ", s[i]);
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      printf("%-10.12f ", vt[j][i]);  // Columns contain singularvectors
    }
    printf("\n");
  }

  printf("\n------------------------------------------------------\n\n");
  fflush(stdout);
#endif

  if (info == 0) {
    return lda;
  } else {

    error(1, "Finding all eigenvalues by singular value decomposition (SVD) "
        "unsuccessful.\n");
    return 0;
  }
}

/********************************************************
 *
 *    Initiate MC steps until one is accepted
 *
 *********************************************************/

double generate_mc_sample(double** const a, double** const v_0,
                          double cost_before, double cost_0, double* w,
                          int* weight, FILE* outfile)
{
  /* If smooth cutoff is enabled, there is an extra parameter (h), which we are
   * adjusting unless it's bounds are chosen to be fixed */

  /* store old parameters incase proposed move isn't accepted */
  double old_params[g_pot.opt_pot.idxlen];
  double* cost_ptr = &cost_before;

  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    old_params[i] = g_pot.opt_pot.table[g_pot.opt_pot.idx[i]];
  }

  int count = 0;
  int mc_decision = 0;

  /* Keep generating trials for this hessian until a move is accepted
     This saves multiple calculations of the same hessian when a move isn't
     accepted */
  while (mc_decision == 0) {
    count++; /* If move not accepted, count current parameters again */

    /* reset parameters to initials params */
    for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
      g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] = old_params[i];
    }

    /* call function continuously until we accept a move for this set of
     * eigenvalues */
    mc_decision = mc_moves(v_0, w, cost_ptr, cost_0, outfile);
  }
  *weight = count;

  /* Return new cost */
  return cost_before;
}

/********************************************************
 *
 *    Take MC step and decide if accepted
 *
 *********************************************************/

int mc_moves(double** v_0, double* w, double* cost_before, double cost_0,
             FILE* outfile)
{
  double lambda[g_pot.opt_pot.idxlen];
  double R;
  double cost_after;
  double delta[g_pot.opt_pot.idxlen];

  R = sqrt(g_param.acc_rescaling);

  /* If eigenvalue is less than eig_pert (defaults to 1), replace it with 1. */
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
#if defined(MIN_STEP)
    if (w[i] > g_param.eig_max) {
      w[i] = g_param.eig_max;
    }
#else /* Use max(lambda,1) */
    /* If negative eigenvalue, set to the minimum of: it's absolute value or the
     * smallest positive eigenvalue */
    if (w[i] < 0) {
      /* Find smallest positive eigenvalue */
      double min_pos_eigenvalue = VERY_LARGE;
      for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
        if ((w[j] > 0) && (w[j] < min_pos_eigenvalue)) {
          min_pos_eigenvalue = w[j];
        }
      }
      /* Set as the smaller of fabs(w[i]) or min_pos_eigenvalue */
      if (fabs(w[i]) <= min_pos_eigenvalue) {
        w[i] = fabs(w[i]);
      } else {
        w[i] = min_pos_eigenvalue;
      }
    } /* end w[i] < 0 */
    /* If eigenvalue is < 1 (including any previously negative), set to 1 */
    if (w[i] < g_param.eig_max) {
      w[i] = g_param.eig_max;
    }
#endif
    double r = R * normdist();
    w[i] = fabs(w[i]);  // Ensured above, leave just incase MIN_STEP used...
    lambda[i] = 1 / sqrt(w[i]);
    lambda[i] *= r;
  }

  /* Matrix multiplication (delta_param[i] = Sum{1}{g_pot.opt_pot.idxlen}
   * [v_0[i][j] * (r[j]/lambda[j])] ) */

  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    delta[i] = 0;
    for (int j = 0; j < g_pot.opt_pot.idxlen; j++) {
      // Rows contain eigenvectors, however we want the
      // ith parameter contribution of each (i.e. columns of v_0)
      delta[i] += v_0[j][i] * lambda[j];
    }
    g_pot.opt_pot.table[g_pot.opt_pot.idx[i]] += delta[i];
  }

  cost_after = calc_forces(g_pot.opt_pot.table, g_calc.force, 0);

  double cost_diff = cost_after - *cost_before;

  /* Accept downhill moves outright */
  if (cost_diff < 0) {
    *cost_before = cost_after;

    return 1;
  }

  /* Monte Carlo step (seeded from srand(time(NULL)) in generate_mc_sample() )
     generate uniform random number [0,1], if less than probability then accept
     change if step accepted, move new cost to cost_before for next cycle */
  double probability = exp(-(g_pot.opt_pot.idxlen * (cost_diff)) /
                           (2.0 * cost_0 * g_param.uq_temp));

  double mc_rand_number = eqdist();

  if (mc_rand_number <= probability) {
    *cost_before = cost_after;

    return 1;
  }

#if defined(DEBUG)
  /* Print out unsuccessful moves */
  fprintf(outfile, "-          ");
  for (int i = 0; i < g_pot.opt_pot.idxlen; i++) {
    fprintf(outfile, "%-10.6g ", g_pot.opt_pot.table[g_pot.opt_pot.idx[i]]);
  }
  fprintf(outfile, "%-10.6g 1          0          -          -\n", cost_after);
#endif

  /* If move not accepted, return 0 */
  return 0;
}

#endif /* UQ && APOT */
