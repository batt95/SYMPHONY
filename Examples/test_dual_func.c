#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

#include "symphony.h"
// feb223
// These paths should be changed
// I still don't figure out how to use makefile to 
// include these files automatically
#include "/Users/feb223/projects/coin/RVF/SYMPHONY/include/sym_types.h"
#include "/Users/feb223/projects/coin/RVF/SYMPHONY/include/sym_master.h"

void generate_rand_array(int m, int l, int u, double *res) {
  for (int i = 0; i < m; i++) {
    res[i] = (rand() % (u - l + 1)) + l;
  }
}

void set_rhs(sym_environment *env, double *rhs, int m){
   for (int i = 0; i < m; i++){
      switch (env->mip->sense[i]){
         case 'E' :
            sym_set_row_upper(env, i, rhs[i]);
            sym_set_row_lower(env, i, rhs[i]);
            break;
         case 'G':
            sym_set_row_lower(env, i, rhs[i]);
            break;
         case 'L':
            sym_set_row_upper(env, i, rhs[i]);
            break;
         default:
            printf("Error\n");
      }
   }
}

int change_lbub_from_disj(double *lb, double *ub, int n, disjunction_desc* disj){
	int i, idx;
	// Fill the lb
	for (i = 0; i < disj->lblen; i++){
		idx = disj->lbvaridx[i];
		lb[idx] = disj->lb[i];
	}
	// Fill the ub
	for (i = 0; i < disj->ublen; i++){
		idx = disj->ubvaridx[i];
		ub[idx] = disj->ub[i];
	}

	return 0;
}

int reset_lbub_from_disj(double *lb, double *ub, 
						double *orig_lb, double *orig_ub, int n, 
						disjunction_desc* disj){
	int i, idx;
	// Revert the lb
	for (i = 0; i < disj->lblen; i++){
		idx = disj->lbvaridx[i];
		lb[idx] = orig_lb[idx];
	}
	// Revert the ub
	for (i = 0; i < disj->ublen; i++){
		idx = disj->ubvaridx[i];
		ub[idx] = orig_ub[idx];
	}
	
	return 0;
}

int main(int argc, char **argv)
{   
   
   int termcode;
   int numTrain = 10;
   int numTests = 10;
   double warmObjVal, coldObjVal, dualFuncObj;

   sym_environment *env_cold = sym_open_environment(); 
   sym_environment *env_warm = sym_open_environment(); 

   sym_parse_command_line(env_warm, argc, argv); 
   // sym_parse_command_line(env_cold, argc, argv); 

   sym_load_problem(env_warm);
   // sym_load_problem(env_cold);
   sym_read_lp(env_cold, "dummy.lp");

   int m = env_warm->mip->m;
   int n = env_warm->mip->n;
   double *rhs = (double*)malloc(sizeof(double) * m);

   sym_set_int_param(env_warm, "verbosity", -2);
   // sym_set_int_param(env_cold, "verbosity", -2);

   // feb223
   // Those are the parameters to be set in order to 
   // keep the branch-and-bound tree valid for RHS changes
   // (E.g. Cuts and reduced cost fixing are not RHS-invariant)
   sym_set_int_param(env_warm, "keep_warm_start", TRUE);
   sym_set_int_param(env_warm, "keep_dual_function_description", TRUE);
   sym_set_int_param(env_warm, "should_use_rel_br", FALSE);
   sym_set_int_param(env_warm, "use_hot_starts", FALSE);
   sym_set_int_param(env_warm, "should_warmstart_node", TRUE);
   sym_set_int_param(env_warm, "sensitivity_analysis", TRUE);
   sym_set_int_param(env_warm, "sensitivity_rhs", true);
   sym_set_int_param(env_warm, "sensitivity_bounds", TRUE);
   sym_set_int_param(env_warm, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env_warm, "do_primal_heuristic", FALSE);
   sym_set_int_param(env_warm, "prep_level", -1);
   sym_set_int_param(env_warm, "tighten_root_bounds", FALSE);
   sym_set_int_param(env_warm, "max_sp_size", 100);
   sym_set_int_param(env_warm, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env_warm, "generate_cgl_cuts", FALSE);
   // sym_set_int_param(env_warm, "max_active_nodes", 1);

   // First solve 
   if ((termcode = sym_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);

   // rhs[0] = 3000;
   // rhs[1] = 3000;

   // set_rhs(env_warm, rhs, m);
   // sym_warm_solve(env_warm);

   // sym_build_dual_func(env_warm);
   // print_dual_function(env_warm);

   // rhs[0] = 2000;
   // rhs[1] = 2000;
   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   sym_get_obj_val(env_warm, &warmObjVal);
   printf("WARM OBJ : %.5f\n", warmObjVal);

   // return 0;

   double *lb = (double*)malloc(sizeof(double) * n);
   double *ub = (double*)malloc(sizeof(double) * n);

   memcpy(lb, env_warm->mip->lb, sizeof(double) * n);
   memcpy(ub, env_warm->mip->ub, sizeof(double) * n);

   dual_func_desc *df = env_warm->warm_start->dual_func;

   const double* elem = df->duals->getElements();
   const int* indices = df->duals->getIndices();
   const double* elem_r = df->rays->getElements();
   const int* indices_r = df->rays->getIndices();
   int first, last, j;

   double bigM = 10000000, epsilon = 1e-4;
   // row data
   int *col_indices = (int*)malloc(sizeof(int) * 10000);
   double *col_elems = (double*)malloc(sizeof(double) * 10000);
   int numelems = 0;
   double intercept = 0;

   // col data
   int *q_indices = (int*)malloc(sizeof(int) * 10000);
   int q_len = 0;
   int *u_indices = (int*)malloc(sizeof(int) * 10000);
   int u_len = 0;
   int *v_indices = (int*)malloc(sizeof(int) * 10000);
   int v_len = 0;
   int col_idx = 1;
   char varname[33];

   // Add beta variables
   for (int i = 0; i < m; i++){
      snprintf(varname, sizeof(varname), "b%d", i);
      sym_add_col(env_cold, 0, NULL, NULL, -bigM, bigM, 0, FALSE, varname);
      col_idx++;
   }

   // Add phi variable
   snprintf(varname, sizeof(varname), "p");
   sym_add_col(env_cold, 0, NULL, NULL, -bigM, bigM, 1, FALSE, varname);
   col_idx++;
   
   for (int i = 0; i < df->num_terms; i++){

      // Add q_t variables
      snprintf(varname, sizeof(varname), "q%d", i);
      sym_add_col(env_cold, 0, NULL, NULL, -bigM, bigM, 0, FALSE, varname);
      q_indices[q_len++] = col_idx++;
      // Add u_t variables
      snprintf(varname, sizeof(varname), "u%d", i);
      sym_add_col(env_cold, 0, NULL, NULL, 0, 1, 0, TRUE, varname);
      u_indices[u_len++] = col_idx++;

      if (df->disj[i].dual_idx >= 0){
         // Add first constraint
         // phi
         col_indices[numelems] = m + 1;
         col_elems[numelems++] = 1;
         // q_t
         col_indices[numelems] = q_indices[q_len - 1];
         col_elems[numelems++] = -1;
         sym_add_row(env_cold, numelems, col_indices, col_elems, 'L', 0, 0);
         
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;

         // Add second constraint
         // phi
         col_indices[numelems] = m + 1;
         col_elems[numelems++] = 1;
         // q_t
         col_indices[numelems] = q_indices[q_len - 1];
         col_elems[numelems++] = -1;
         // u_t
         col_indices[numelems] = u_indices[u_len - 1];
         col_elems[numelems++] = -bigM;
         sym_add_row(env_cold, numelems, col_indices, col_elems, 'G', -bigM, 0);
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;

         // Add third constraint
         // q_t
         col_indices[numelems] = q_indices[q_len - 1];
         col_elems[numelems++] = 1;

         first = df->duals->getVectorFirst(df->disj[i].dual_idx);
		   last = df->duals->getVectorLast(df->disj[i].dual_idx);

         for (j = first; j < last && indices[j] < m; j++){
            col_indices[numelems] = indices[j] + 1;
            col_elems[numelems++] = -elem[j];
         }

         change_lbub_from_disj(lb, ub, n, df->disj + i);
         intercept = 0;
         
         for (; j < last; j++){
            if (elem[j] > epsilon){
               intercept += elem[j] * lb[indices[j] - m];
            } 
            else if (elem[j] < -epsilon){
               intercept += elem[j] * ub[indices[j] - m];
            } 
         }
         
         reset_lbub_from_disj(lb, ub, env_warm->mip->lb, env_warm->mip->ub, n, df->disj + i);

         sym_add_row(env_cold, numelems, col_indices, col_elems, 'G', intercept, 0);
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;
         
      } 
      if (df->disj[i].ray_idx >= 0){
         // Add v_t variables
         snprintf(varname, sizeof(varname), "v%d", i);
         sym_add_col(env_cold, 0, NULL, NULL, 0, 1, 0, TRUE, varname);
         v_indices[v_len++] = col_idx++;
   
         // Add first constraint
         // phi
         col_indices[numelems] = m + 1;
         col_elems[numelems++] = 1;
         // q_t
         col_indices[numelems] = q_indices[q_len - 1];
         col_elems[numelems++] = -1;
         // v_t
         col_indices[numelems] = v_indices[v_len - 1];
         col_elems[numelems++] = -bigM;
         sym_add_row(env_cold, numelems, col_indices, col_elems, 'L', 0, 0);
         
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;

         // Add second constraint
         // phi
         col_indices[numelems] = m + 1;
         col_elems[numelems++] = 1;
         // q_t
         col_indices[numelems] = q_indices[q_len - 1];
         col_elems[numelems++] = -1;
         // u_t
         col_indices[numelems] = u_indices[u_len - 1];
         col_elems[numelems++] = -bigM;
         sym_add_row(env_cold, numelems, col_indices, col_elems, 'G', -bigM, 0);
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;

         // Add u_t + v_t <= 1
         // u_t
         col_indices[numelems] = u_indices[u_len - 1];
         col_elems[numelems++] = 1;
         // v_t
         col_indices[numelems] = v_indices[v_len - 1];
         col_elems[numelems++] = 1;
         sym_add_row(env_cold, numelems, col_indices, col_elems, 'L', 1, 0);
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;

         // Add rays constraints
         // v_t
         col_indices[numelems] = v_indices[v_len - 1];
         col_elems[numelems++] = -bigM;

         first = df->rays->getVectorFirst(df->disj[i].ray_idx);
		   last = df->rays->getVectorLast(df->disj[i].ray_idx);

         for (j = first; j < last && indices_r[j] < m; j++){
            col_indices[numelems] = indices_r[j] + 1;
            col_elems[numelems++] = -elem_r[j];
         }

         change_lbub_from_disj(lb, ub, n, df->disj + i);
         intercept = 0;
         
         for (; j < last; j++){
            if (elem_r[j] > epsilon){
               intercept += elem_r[j] * lb[indices_r[j] - m];
            } 
            else if (elem_r[j] < -epsilon){
               intercept += elem_r[j] * ub[indices_r[j] - m];
            } 
         }

         sym_add_row(env_cold, numelems, col_indices, col_elems, 
                     'L', -intercept, 0);
         sym_add_row(env_cold, numelems, col_indices, col_elems, 
                     'G', -bigM +epsilon -intercept, 0);

         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;

         // Add third constraint
         // q_t
         col_indices[numelems] = q_indices[q_len - 1];
         col_elems[numelems++] = 1;

         first = df->duals->getVectorFirst(0);
		   last = df->duals->getVectorLast(0);

         for (j = first; j < last && indices[j] < m; j++){
            col_indices[numelems] = indices[j] + 1;
            col_elems[numelems++] = -elem[j];
         }

         intercept = 0;
         
         for (; j < last; j++){
            if (elem[j] > epsilon){
               intercept += elem[j] * lb[indices[j] - m];
            } 
            else if (elem[j] < -epsilon){
               intercept += elem[j] * ub[indices[j] - m];
            } 
         }
         
         reset_lbub_from_disj(lb, ub, env_warm->mip->lb, env_warm->mip->ub, 
                              n, df->disj + i);

         sym_add_row(env_cold, numelems, col_indices, col_elems, 
                     'G', intercept, 0);
         // clear
         for (int j = numelems; j >= 0; --j){
            col_indices[j] = 0;
            col_elems[j] = 0;
         }
         numelems = 0;
      }
   }

   // Add fourth constraint
   for (int i = 0; i < u_len; i++){
      col_indices[numelems] = u_indices[i];
      col_elems[numelems++] = 1;
   }
   sym_add_row(env_cold, numelems, col_indices, col_elems, 'E', 1, 0);
   // clear
   for (int j = numelems; j >= 0; --j){
      col_indices[j] = 0;
      col_elems[j] = 0;
   }
   numelems = 0;

   // Fix beta variables
   for (int i = 0; i < m; i++){
      col_indices[numelems] = i + 1;
      col_elems[numelems++] = 1;
      sym_add_row(env_cold, numelems, col_indices, col_elems, 'E', 2000, 0);
      for (int j = numelems; j >= 0; --j){
         col_indices[j] = 0;
         col_elems[j] = 0;
      }
      numelems = 0;
   }

   sym_write_lp(env_cold, "mastertest");
   sym_solve(env_cold);

   rhs[0] = 3000;
   rhs[1] = 3000;

   set_rhs(env_warm, rhs, m);
   sym_solve(env_warm);
   sym_get_obj_val(env_warm, &warmObjVal);
   printf("WARM OBJ : %.5f\n", warmObjVal);
   
   free(rhs);
   sym_close_environment(env_warm);
   sym_close_environment(env_cold);
   return 0;
}  