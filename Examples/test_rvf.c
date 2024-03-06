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
#include "../include/sym_types.h"
#include "../include/sym_master.h"

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

double *linspace(double a, double b, int pieces){

   double *res = (double *)malloc(sizeof(double) * pieces);
   double step = fabs(a + b)/((double)pieces);

   res[0] = a;
   res[pieces - 1] = b;

   for (int i = 1; i < pieces - 1; i++){
      res[i] = res[i - 1] + step;
   }

   return res;
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
   sym_set_int_param(env_warm, "max_active_nodes", 1);

   sym_set_int_param(env_cold, "keep_warm_start", TRUE);
   sym_set_int_param(env_cold, "keep_dual_function_description", TRUE);
   sym_set_int_param(env_cold, "should_use_rel_br", FALSE);
   sym_set_int_param(env_cold, "use_hot_starts", FALSE);
   sym_set_int_param(env_cold, "should_warmstart_node", TRUE);
   sym_set_int_param(env_cold, "sensitivity_analysis", TRUE);
   sym_set_int_param(env_cold, "sensitivity_rhs", true);
   sym_set_int_param(env_cold, "sensitivity_bounds", TRUE);
   sym_set_int_param(env_cold, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env_cold, "do_primal_heuristic", FALSE);
   sym_set_int_param(env_cold, "prep_level", -1);
   sym_set_int_param(env_cold, "tighten_root_bounds", FALSE);
   sym_set_int_param(env_cold, "max_sp_size", 100);
   sym_set_int_param(env_cold, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env_cold, "generate_cgl_cuts", FALSE);
   sym_set_int_param(env_cold, "max_active_nodes", 1);

   //----------------------
   // linspace numpy-like 
   //----------------------
   // double a = -16, b = 5, zerotol = 1e-7;
   // int pieces = 2000;
   // double *zeta_lst = linspace(a, b, pieces);
   // double *zeta_lst_1 = linspace(-15, 5, 10);
   // double *rvf_lst  = (double*)malloc(sizeof(double) * pieces);
   // double *df_lst   = (double*)malloc(sizeof(double) * pieces);

   int num_objs = 2;
   double *rhs = (double*)malloc(sizeof(double) * num_objs);
   
   // First solve 
   if ((termcode = sym_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }
   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 0, &dualFuncObj);


   printf("======================================\n");
   rhs[0] = -2747.0;
   rhs[1] = -3132.0;

   set_rhs(env_warm, rhs, num_objs);
    
   if ((termcode = sym_warm_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   if (sym_is_proven_optimal(env_warm)){
      printf("  PROVEN OPTIMAL!\n");
   } else if (sym_is_proven_primal_infeasible(env_warm)){
      printf("  PROVEN INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

   printf("======================================\n");
   rhs[0] = -2747.0;
   rhs[1] = -2740.5;

   set_rhs(env_warm, rhs, num_objs);
    
   if ((termcode = sym_warm_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   if (sym_is_proven_optimal(env_warm)){
      printf("  PROVEN OPTIMAL!\n");
   } else if (sym_is_proven_primal_infeasible(env_warm)){
      printf("  PROVEN INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

   printf("======================================\n");
   rhs[0] = -2747.0;
   rhs[1] = -2349.0;

   set_rhs(env_warm, rhs, num_objs);
    
   if ((termcode = sym_warm_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   if (sym_is_proven_optimal(env_warm)){
      printf("  PROVEN OPTIMAL!\n");
   } else if (sym_is_proven_primal_infeasible(env_warm)){
      printf("  PROVEN INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

   printf("======================================\n");
   rhs[0] = -2747.0;
   rhs[1] = -1957.5;

   set_rhs(env_warm, rhs, num_objs);
    
   if ((termcode = sym_warm_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   if (sym_is_proven_optimal(env_warm)){
      printf("  PROVEN OPTIMAL!\n");
   } else if (sym_is_proven_primal_infeasible(env_warm)){
      printf("  PROVEN INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

   printf("======================================\n");
   rhs[0] = -2747.0;
   rhs[1] = -1566.0;

   set_rhs(env_warm, rhs, num_objs);
    
   if ((termcode = sym_warm_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   if (sym_is_proven_optimal(env_warm)){
      printf("  PROVEN OPTIMAL!\n");
   } else if (sym_is_proven_primal_infeasible(env_warm)){
      printf("  PROVEN INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

   // // check_dual_solutions(env_cold->mip, env_warm->warm_start->dual_func);

   // rhs[0] = -55.5;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;
      
   // // Set the new RHS
   // set_rhs(env_warm, rhs, m);

   // if ((termcode = sym_warm_solve(env_warm)) < 0){
   //    printf("PROBLEM INFEASIBLE!\n");
   // }

   // sym_get_obj_val(env_warm, &warmObjVal);
   // printf("WARM OBJ : %.5f\n", warmObjVal);

   // sym_build_dual_func(env_warm);

   // rhs[0] = 40.0/9.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -55.5;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -73.0/6.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -10;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -73.0/6.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;
      
   // // Set the new RHS
   // set_rhs(env_warm, rhs, m);

   // if ((termcode = sym_warm_solve(env_warm)) < 0){
   //    printf("PROBLEM INFEASIBLE!\n");
   // }

   // sym_get_obj_val(env_warm, &warmObjVal);
   // printf("WARM OBJ : %.5f\n", warmObjVal);

   // sym_build_dual_func(env_warm);

   // rhs[0] = 40.0/9.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -55.5;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -73.0/6.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -10;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

   // rhs[0] = -10;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;
      
   // // Set the new RHS
   // set_rhs(env_warm, rhs, m);

   // if ((termcode = sym_warm_solve(env_warm)) < 0){
   //    printf("PROBLEM INFEASIBLE!\n");
   // }

   // sym_get_obj_val(env_warm, &warmObjVal);
   // printf("WARM OBJ : %.5f\n", warmObjVal);

   // sym_build_dual_func(env_warm);

   // rhs[0] = 40.0/9.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
   // printf("Dual function evaluates %.7f\n", dualFuncObj);

   // rhs[0] = -55.5;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
   // printf("Dual function evaluates %.7f\n", dualFuncObj);

   // rhs[0] = -73.0/6.0;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
   // printf("Dual function evaluates %.7f\n", dualFuncObj);

   // rhs[0] = -10;
   // rhs[1] = 4;
   // rhs[2] = 5;
   // rhs[3] = 5;

   // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
   // printf("Dual function evaluates %.7f\n", dualFuncObj);


   // Memory clean-ups
   sym_close_environment(env_warm);
   sym_close_environment(env_cold);

   // free(zeta_lst);
   // free(rvf_lst);
   // free(df_lst);

   return 0;
}  

