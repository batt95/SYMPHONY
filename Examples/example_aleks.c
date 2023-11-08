#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

#include "symphony.h"
// feb223
// These paths should be changed
#include "../include/sym_types.h"
#include "../include/sym_master.h"

int main(int argc, char **argv)
{   
   int termcode = 0;
   double objVal = 0.0;
   int n = 6;

   double *sol = (double*)malloc(sizeof(double) * n);
   // Create SYMPHONY environment
   sym_environment *env = sym_open_environment(); 

   // Parse command line parameters and instance name
   sym_parse_command_line(env, argc, argv); 

   // Load the problem
   sym_load_problem(env);

   // Reduce outputs
   sym_set_int_param(env, "verbosity", -2);

   // Set any type of parameters
   // These are the set of parameters I use when 
   // I want to work on B&B tree for things that are valid 
   // for all RHS (e.g., like the dual functions)
   // See the documentation for further details
   // on what these parameters do
   sym_set_int_param(env, "keep_warm_start", TRUE);
   sym_set_int_param(env, "keep_dual_function_description", TRUE); // This parameter is my WIP
   sym_set_int_param(env, "should_use_rel_br", FALSE);
   sym_set_int_param(env, "use_hot_starts", FALSE);
   sym_set_int_param(env, "should_warmstart_node", TRUE);
   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "sensitivity_rhs", true);
   sym_set_int_param(env, "sensitivity_bounds", TRUE);
   sym_set_int_param(env, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env, "do_primal_heuristic", FALSE);
   sym_set_int_param(env, "prep_level", -1);
   sym_set_int_param(env, "tighten_root_bounds", FALSE);
   sym_set_int_param(env, "max_sp_size", 100);
   sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env, "generate_cgl_cuts", FALSE);
   sym_set_int_param(env, "max_active_nodes", 1);

   // Solve 
   if ((termcode = sym_solve(env)) < 0){
      printf("PROBLEM INFEASIBLE!\n");
   }
   
   // Access objective value
   sym_get_obj_val(env, &objVal);

   // Access solution vector
   sym_get_col_solution(env, sol);
   printf("=======================\n");
   printf("        SOLUTION       \n");
   printf("=======================\n");
   printf(" - objval: %.5f\n", objVal);
   for (int j = 0; j < n; j++){
      printf(" - x[%d] : %.5f\n", j, sol[j]);
   }
   printf("=======================\n\n");   

   // Access the B&B Tree and get all useful information
   sym_read_tree_info(env);

   // Change RHS of the 0-th constraint ('<')
   sym_set_row_upper(env, 0, 1000);
   // sym_set_row_lower(env, 0, 10);

   // Solve again using sym_warm_solve()
   // calling sym_solve() destroys previous information
   if ((termcode = sym_warm_solve(env)) < 0){
      printf("PROBLEM INFEASIBLE!\n");
   } 

   // Access objective value
   sym_get_obj_val(env, &objVal);

   // Access solution vector
   sym_get_col_solution(env, sol);
   printf("=======================\n");
   printf("        SOLUTION       \n");
   printf("=======================\n");
   printf(" - objval: %.5f\n", objVal);
   for (int j = 0; j < n; j++){
      printf(" - x[%d] : %.5f\n", j, sol[j]);
   }
   printf("=======================\n\n");   

   // Access the modified B&B Tree and get all useful information
   sym_read_tree_info(env);

   // Close the environment
   sym_close_environment(env);

   return 0;
}  