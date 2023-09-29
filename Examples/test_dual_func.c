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

int main(int argc, char **argv)
{   
   
   int termcode;
   int numTests = 10;
   double warmObjVal, coldObjVal, dualFuncObj;

   sym_environment *env_cold = sym_open_environment(); 
   sym_environment *env_warm = sym_open_environment(); 

   sym_parse_command_line(env_warm, argc, argv); 
   sym_parse_command_line(env_cold, argc, argv); 

   sym_load_problem(env_warm);
   sym_load_problem(env_cold);

   int m = env_cold->mip->m;
   double *rhs = (double*)malloc(sizeof(double) * m);

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

   // First solve 
   if ((termcode = sym_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }

   sym_build_dual_func(env_warm);
   
   sym_get_obj_val(env_warm, &warmObjVal);
   printf("WARM OBJ : %.5f\n", warmObjVal);

   // -----------------------------------------------------
   // RHS Test 
   // -----------------------------------------------------
   for(int i = 0; i < numTests; i++){
      printf("=================================\n");
      
      printf("RHS %d\n", i);

      // Generate a random new rhs
      generate_rand_array(m, 0, 2000, rhs);
      
      // Set the new RHS
      set_rhs(env_cold, rhs, m);

      sym_write_lp(env_cold, "peppe.lp");

      // Cold solve
      if ((termcode = sym_solve(env_cold)) < 0){
         printf("PROBLEM INFEASIBLE!\n");
      } 

      sym_get_obj_val(env_cold, &coldObjVal);
      printf("OPT OBJ : %.7f\n", coldObjVal);

      printf("RHS: (%.2f, %.2f)\n", rhs[0], rhs[1]);
      sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
      printf("DUAL FUNC: %.5f\n", dualFuncObj);

      assert(dualFuncObj - coldObjVal < 1e-5 || fabs(dualFuncObj - coldObjVal) < 1e-5);
   
   }
   printf("ENDING\n");

   sym_close_environment(env_warm);
   sym_close_environment(env_cold);
   return 0;
}  