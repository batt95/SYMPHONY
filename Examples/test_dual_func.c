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

int clean_row_arrays(int *col_indices, double *col_elems,
            int *numelems){
   // clear
   for (int j = *numelems; j >= 0; --j){
      col_indices[j] = 0;
      col_elems[j] = 0;
   }
   *numelems = 0;
}

int add_q_u_v_variables(sym_environment *env,
            disjunction_desc* disj, int i,
            int *q_indices, int *q_len,
            int *u_indices, int *u_len,
            int *v_indices, int *v_len,
            int *col_idx){
   
   char varname[33];
   double bigM = 10000000;

   // Add q_t variables
   snprintf(varname, sizeof(varname), "q%d", i);
   sym_add_col(env, 0, NULL, NULL, -bigM, bigM, 0, FALSE, varname);
   q_indices[(*q_len)++] = (*col_idx)++;

   // Add u_t variables
   snprintf(varname, sizeof(varname), "u%d", i);
   sym_add_col(env, 0, NULL, NULL, 0, 1, 0, TRUE, varname);
   u_indices[(*u_len)++] = (*col_idx)++;

   if (disj->raylen){
      // Add v_t variables
      snprintf(varname, sizeof(varname), "v%d", i);
      sym_add_col(env, 0, NULL, NULL, 0, 1, 0, TRUE, varname);
      v_indices[(*v_len)++] = (*col_idx)++;
   }
   
   return 0;
}

int add_phi_constrs(sym_environment *env,
            disjunction_desc* disj, int i, int m,
            int *q_indices, int q,
            int *u_indices, int u,
            int *v_indices, int v,
            int *col_indices, double *col_elems,
            int *numelems){
   double bigM = 10000000; 

   // Add first constraint
   // phi
   col_indices[*numelems] = m + 1;
   col_elems[(*numelems)++] = 1;
   // q_t
   col_indices[(*numelems)] = q_indices[q];
   col_elems[(*numelems)++] = -1;

   if (disj->raylen){
      // v_t
      col_indices[(*numelems)] = v_indices[v];
      col_elems[(*numelems)++] = -bigM;
   }

   sym_add_row(env, (*numelems), col_indices, col_elems, 'L', 0, 0);
   clean_row_arrays(col_indices, col_elems, numelems);

   // Add second constraint
   // phi
   col_indices[(*numelems)] = m + 1;
   col_elems[(*numelems)++] = 1;
   // q_t
   col_indices[(*numelems)] = q_indices[q];
   col_elems[(*numelems)++] = -1;
   // u_t
   col_indices[(*numelems)] = u_indices[u];
   col_elems[(*numelems)++] = -bigM;

   sym_add_row(env, (*numelems), col_indices, col_elems, 'G', -bigM, 0);
   clean_row_arrays(col_indices, col_elems, numelems);

   return 0;
}

int add_dual_constrs(sym_environment *env,
            dual_func_desc *df,
            disjunction_desc* disj, int i, int m, int n,
            int *q_indices, int q,
            int *u_indices, int u,
            int *v_indices, int v,
            int *col_indices, double *col_elems,
            int *numelems,
            double *lb, double *ub){
   
   double intercept, epsilon = 1e-4;
   int first, last, row, j;

   const double* elem = df->duals->getElements();
   const int* indices = df->duals->getIndices();

   for (int d = 0; d < disj->leaflen; d++){
      // Add dual constraints
      // q_t
      col_indices[(*numelems)] = q_indices[q];
      col_elems[(*numelems)++] = 1;

      row = disj->leaf_idx[d];

      first = df->duals->getVectorFirst(row);
      last = df->duals->getVectorLast(row);

      for (j = first; j < last && indices[j] < m; j++){
         col_indices[(*numelems)] = indices[j] + 1;
         col_elems[(*numelems)++] = -elem[j];
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

      sym_add_row(env, (*numelems), col_indices, col_elems, 'G', intercept, 0);
      clean_row_arrays(col_indices, col_elems, numelems);
   }

   return 0;
}

int add_ray_constrs(sym_environment *env,
            dual_func_desc *df,
            disjunction_desc* disj, int i, int m, int n,
            int *q_indices, int q,
            int *u_indices, int u,
            int *v_indices, int v,
            int *col_indices, double *col_elems,
            int *numelems,
            double *lb, double *ub){

   if (!disj->raylen){
      return 0;
   }           

   double bigM = 10000000;
   double intercept = 0, epsilon = 1e-4;
   int first, last, row, j;

   const double* elem = df->rays->getElements();
   const int* indices = df->rays->getIndices();

   // Add u_t + v_t <= 1
   // u_t
   col_indices[(*numelems)] = u_indices[u];
   col_elems[(*numelems)++] = 1;
   // v_t
   col_indices[(*numelems)] = v_indices[v];
   col_elems[(*numelems)++] = 1;
   sym_add_row(env, (*numelems), col_indices, col_elems, 'L', 1, 0);
   clean_row_arrays(col_indices, col_elems, numelems);

   for (int r = 0; r < disj->raylen; r++){
      // Add rays constraints
      // v_t
      col_indices[(*numelems)] = v_indices[v];
      col_elems[(*numelems)++] = -bigM;

      row = disj->ray_idx[r];

      first = df->rays->getVectorFirst(row);
      last = df->rays->getVectorLast(row);

      for (j = first; j < last && indices[j] < m; j++){
         col_indices[(*numelems)] = indices[j] + 1;
         col_elems[(*numelems)++] = -elem[j];
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

      sym_add_row(env, (*numelems), col_indices, col_elems, 
                  'L', -intercept, 0);
      sym_add_row(env, (*numelems), col_indices, col_elems, 
                  'G', - bigM + epsilon - intercept, 0);

      clean_row_arrays(col_indices, col_elems, numelems);
   }

   return 0;

}

int add_u_constr(sym_environment *env,
            int *u_indices, int u_len,
            int *col_indices, double *col_elems,
            int *numelems){

   // Add fourth constraint
   for (int i = 0; i < u_len; i++){
      col_indices[(*numelems)] = u_indices[i];
      col_elems[(*numelems)++] = 1;
   }
   sym_add_row(env, (*numelems), col_indices, col_elems, 'E', 1, 0);
   clean_row_arrays(col_indices, col_elems, numelems);
}



int main(int argc, char **argv)
{   
   
   int termcode;
   double warmObjVal, coldObjVal, dualFuncObj;

   sym_environment *master, *recourse = sym_open_environment(); 

   // sym_read_lp(recourse, "mastertest.lp");
   // sym_solve(recourse);

   // sym_close_environment(recourse);
   // return 0;

   sym_set_int_param(recourse, "verbosity", -2);
   sym_parse_command_line(recourse, argc, argv); 
   sym_load_problem(recourse);

   sym_set_int_param(recourse, "keep_warm_start", TRUE);
   sym_set_int_param(recourse, "keep_dual_function_description", TRUE);
   sym_set_int_param(recourse, "should_use_rel_br", FALSE);
   sym_set_int_param(recourse, "use_hot_starts", FALSE);
   sym_set_int_param(recourse, "should_warmstart_node", TRUE);
   sym_set_int_param(recourse, "sensitivity_analysis", TRUE);
   sym_set_int_param(recourse, "sensitivity_rhs", true);
   sym_set_int_param(recourse, "sensitivity_bounds", TRUE);
   sym_set_int_param(recourse, "set_obj_upper_lim", FALSE);
   sym_set_int_param(recourse, "do_primal_heuristic", FALSE);
   sym_set_int_param(recourse, "prep_level", -1);
   sym_set_int_param(recourse, "tighten_root_bounds", FALSE);
   sym_set_int_param(recourse, "max_sp_size", 100);
   sym_set_int_param(recourse, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(recourse, "generate_cgl_cuts", FALSE); 

   int m = recourse->mip->m;
   int n = recourse->mip->n;
   int num_rhs = 5;
   double *rhs = (double*)malloc(sizeof(double) * m * num_rhs);
   double *recourse_objs = (double*)malloc(sizeof(double) * num_rhs);
   double *master_objs = (double*)malloc(sizeof(double) * num_rhs);

   rhs[0] = 1500;
   rhs[1] = 1500;

   rhs[2] = 700;
   rhs[3] = 700;

   rhs[4] = 2000;
   rhs[5] = 2000;

   rhs[6] = 5000;
   rhs[7] = 5000;

   rhs[8] = 3500;
   rhs[9] = 3500;

   int *q_indices = NULL, *u_indices = NULL, *v_indices = NULL;
   int q_len = 0, u_len = 0, v_len = 0;
   int col_idx = 1;
   int *col_indices = NULL;
   double *col_elems = NULL;
   int numelems = 0;
   double intercept = 0;
   dual_func_desc *df = NULL;
   disjunction_desc *disj;
   int max_num_cols;
   double bigM = 10000000, epsilon = 1e-4;
   int num_scen = 1;
   char varname[33];

   double *lb = (double*)malloc(sizeof(double) * n);
   double *ub = (double*)malloc(sizeof(double) * n);

   memcpy(lb, recourse->mip->lb, sizeof(double) * n);
   memcpy(ub, recourse->mip->ub, sizeof(double) * n);

   for (int i = 0; i < num_rhs; i++){
      set_rhs(recourse, rhs + m*i, m);
      sym_warm_solve(recourse);
      sym_build_dual_func(recourse);
      print_dual_function(recourse);

      // free memory
      if (col_indices) free(col_indices);
      if (col_elems)   free(col_elems);
      if (q_indices)   free(q_indices);
      if (u_indices)   free(u_indices);
      if (v_indices)   free(v_indices);
      q_len = 0, u_len = 0, v_len = 0;
      col_idx = 1;

      master = sym_open_environment();
      sym_set_int_param(master, "verbosity", -2);
      // This dummy is needed now since we are starting
      // from an empty problem and SYMPHONY mess it up.
      // In GenBenders there will be the base master problem
      sym_read_lp(master, "dummy.lp");

      df = recourse->warm_start->dual_func;   
      max_num_cols = 1 + m + num_scen + 3 * df->num_terms;

      // row data
      col_indices = (int*)malloc(sizeof(int) * max_num_cols);
      col_elems = (double*)malloc(sizeof(double) * max_num_cols);
      numelems = 0;
      intercept = 0;

      // col data
      q_indices = (int*)malloc(sizeof(int) * df->num_terms);
      u_indices = (int*)malloc(sizeof(int) * df->num_terms);
      v_indices = (int*)malloc(sizeof(int) * df->num_terms);

      // Add beta variables
      for (int i = 0; i < m; i++){
         snprintf(varname, sizeof(varname), "b%d", i);
         sym_add_col(master, 0, NULL, NULL, -bigM, bigM, 0, FALSE, varname);
         col_idx++;
      }

      snprintf(varname, sizeof(varname), "phi");
      sym_add_col(master, 0, NULL, NULL, -bigM, bigM, 1, FALSE, varname);
      col_idx++;

      for (int i = 0; i < df->num_terms; i++){
         disj = df->disj + i;
         add_q_u_v_variables(master, disj, i, 
                              q_indices, &q_len,
                              u_indices, &u_len,
                              v_indices, &v_len, 
                              &col_idx);

         // sym_write_lp(master, "mastertest");

         add_phi_constrs(master, disj, i, m,
                              q_indices, q_len - 1,
                              u_indices, u_len - 1,
                              v_indices, v_len - 1,
                              col_indices, col_elems, &numelems);


         // sym_write_lp(master, "mastertest");

         change_lbub_from_disj(lb, ub, n, disj);

         add_dual_constrs(master, df, disj, i, m, n,
                              q_indices, q_len - 1,
                              u_indices, u_len - 1,
                              v_indices, v_len - 1, 
                              col_indices, col_elems, &numelems,
                              lb, ub);

         // sym_write_lp(master, "mastertest");
         
         add_ray_constrs(master, df, disj, i, m, n,
                              q_indices, q_len - 1,
                              u_indices, u_len - 1,
                              v_indices, v_len - 1, 
                              col_indices, col_elems, &numelems,
                              lb, ub);

         // sym_write_lp(master, "mastertest");
         
         reset_lbub_from_disj(lb, ub, recourse->mip->lb, recourse->mip->ub,
                              n, disj);
      }

      add_u_constr(master, u_indices, u_len, col_indices, col_elems, &numelems);
         

      // Fix beta variables
      for (int j = 0; j < m; j++){
         col_indices[numelems] = j + 1;
         col_elems[numelems++] = 1;
         sym_add_row(master, numelems, col_indices, col_elems, 'E', rhs[m*i + j], 0);
         clean_row_arrays(col_indices, col_elems, &numelems);
      }

      sym_write_lp(master, "mastertest");
      sym_set_int_param(master, "verbosity", 1);
      sym_solve(master);
      sym_get_obj_val(recourse, &warmObjVal);
      printf("recourse obj val: %.5f\n", warmObjVal);
      sym_evaluate_dual_function(recourse, rhs, 0, &dualFuncObj);
      printf("dual func obj val: %.5f\n", dualFuncObj);
      recourse_objs[i] = warmObjVal;
      sym_get_obj_val(master, &warmObjVal);
      master_objs[i] = warmObjVal;
      printf("master obj val: %.5f\n", warmObjVal);
      sym_close_environment(master);
   }
   
   free(rhs);
   free(recourse_objs);
   free(col_indices);
   free(col_elems);
   free(q_indices);
   free(u_indices);
   free(v_indices);
   free(lb);
   free(ub);
   sym_close_environment(recourse);
   return 0;
}  