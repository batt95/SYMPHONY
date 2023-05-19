#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

// feb223
// should we use the Osi Interface?
// #include "OsiSymSolverInterface.hpp"

#include "symphony.h"
// feb223
// These paths should be changed
// I still don't figure out how to use makefile to 
// include these files automatically
#include "../include/sym_types.h"
#include "../include/sym_master.h"

double* generate_rand_array(int n, int l, int u) {
  double* res = (double *)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    res[i] = (rand() % (u - l + 1)) + l;
  }
  return res;
}

void set_rhs(sym_environment *env, double *rhs, int m, char sense){
   switch (sense){
      case 'E' :
         for (int i = 0; i < m; i++){
            sym_set_row_upper(env, i, rhs[i]);
            sym_set_row_lower(env, i, rhs[i]);
         }
         break;
      case 'G':
         for (int i = 0; i < m; i++){
            sym_set_row_upper(env, i, rhs[i]);
         }
         break;
      case 'L':
         for (int i = 0; i < m; i++){
            sym_set_row_lower(env, i, rhs[i]);
         }
         break;

      default:
         printf("Error\n");
   }
   
   
}


void set_obj(sym_environment *env, double *coeff, int n){
   for (int i = 0; i < n; i++){
      sym_set_obj_coeff(env, i, coeff[i]);
   }
}

int check_objs(sym_environment *env1, sym_environment *env2){
   double *ObjVal1 = (double *) malloc(sizeof(double));
   double *ObjVal2 = (double *) malloc(sizeof(double));
   sym_get_obj_val(env1, ObjVal1);
   sym_get_obj_val(env2, ObjVal2);
   return fabs(*ObjVal1 - *ObjVal2) < 1e-5;
}

int main(int argc, char **argv)
{   
      
   double inf = 1e30;
   int i = 0, termcode;
   int numRhs = 3;
   int n = 6;
   int m = 1;

   double *rhs = (double *) malloc(numRhs * sizeof(double));
   double *coeff = (double *) malloc(n * sizeof(double));
   double *warmCoeff = (double *) malloc(n * sizeof(double));
   double *coldCoeff = (double *) malloc(n * sizeof(double));
   double *coldObjVal = (double *) malloc(sizeof(double));
   double *warmObjVal = (double *) malloc(sizeof(double));
   double *coldSol = (double *) malloc(n * sizeof(double));
   double *warmSol = (double *) malloc(n * sizeof(double));
   double *duals = (double *) malloc(m * sizeof(double));

   rhs[0] = 11.5;
   rhs[1] = 4;
   rhs[2] = 10;

   int* rhs_ind = (int *) malloc(m * sizeof(int));
   double *lb_rhs = (double *) malloc(sizeof(double));
   
   sym_environment *env_warm = sym_open_environment(); 
   sym_environment *env_cold = sym_open_environment(); 

   // Maybe we need to store the warm_start_desc
   warm_start_desc * ws; 

   sym_parse_command_line(env_warm, argc, argv); 
   sym_parse_command_line(env_cold, argc, argv);   
   sym_load_problem(env_warm);
   sym_load_problem(env_cold);

   sym_set_int_param(env_warm, "verbosity", -2);
   sym_set_int_param(env_cold, "verbosity", -2);
   // feb223
   // Those are the parameters to be set in order to 
   // keep the branch-and-bound tree valid for RHS changes
   // (E.g. Cuts and reduced cost fixing are not RHS-invariant)
   sym_set_int_param(env_warm, "keep_warm_start", TRUE);
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
   sym_set_int_param(env_warm, "max_active_nodes", 1);

   // First solve both 
   if ((termcode = sym_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }
   
   print_tree(env_warm->warm_start->rootnode);
   ws = sym_get_warm_start(env_warm, TRUE);

   sym_get_obj_val(env_warm, warmObjVal);
   printf("WARM OBJ : %.5f\n", *warmObjVal);

   sym_get_col_solution(env_warm, warmSol);
   for (int j = 0; j < n; j++){
      printf("X%d : %.5f\n", j, warmSol[j]);
   }

   print_tree(ws->rootnode);
   
   // sym_write_lp(env_warm, "LP");

   // -----------------------------------------------------
   // RHS Test 
   // -----------------------------------------------------
   for(int i = 0; i < numRhs; i++){
      printf("=================================\n");
      // Set the new RHS
      sym_set_row_upper(env_warm, 0, rhs[i]);
      sym_set_row_lower(env_warm, 0, rhs[i]);

      // Warm solve
      if ((termcode = sym_warm_solve(env_warm)) < 0){
         printf("PROBLEM INFEASIBLE!\n");
      } 

      ws = sym_get_warm_start(env_warm, TRUE);
      sym_get_obj_val(env_warm, warmObjVal);
      printf("WARM OBJ : %.5f\n", *warmObjVal);

      sym_get_col_solution(env_warm, warmSol);
      for (int j = 0; j < n; j++){
         printf("X%d : %.5f\n", j, warmSol[j]);
      }
      print_tree(ws->rootnode);
   
   }
   printf("ENDING\n");

   sym_close_environment(env_warm);
   sym_close_environment(env_cold);

   return 0;
}  