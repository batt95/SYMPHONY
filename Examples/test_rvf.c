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
   sym_set_int_param(env_warm, "use_hot_starts", TRUE);
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
   // sym_set_int_param(env_warm, "max_presolve_iter", 0);
   // sym_set_int_param(env_warm, "limit_strong_branching_time", 0);

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
   int num_zetas = 81;
   double zetas[][2] = {
        {-2747.0, -3132.0}, {-2747.0, -2740.5}, {-2747.0, -2349.0}, {-2747.0, -1957.5}, {-2747.0, -1566.0}, {-2747.0, -1174.5}, {-2747.0, -783.0}, {-2747.0, -391.5}, {-2747.0, -0.0},
        {-2403.625, -3132.0}, {-2403.625, -2740.5}, {-2403.625, -2349.0}, {-2403.625, -1957.5}, {-2403.625, -1566.0}, {-2403.625, -1174.5}, {-2403.625, -783.0}, {-2403.625, -391.5}, {-2403.625, -0.0},
        {-2060.25, -3132.0}, {-2060.25, -2740.5}, {-2060.25, -2349.0}, {-2060.25, -1957.5}, {-2060.25, -1566.0}, {-2060.25, -1174.5}, {-2060.25, -783.0}, {-2060.25, -391.5}, {-2060.25, -0.0},
        {-1716.875, -3132.0}, {-1716.875, -2740.5}, {-1716.875, -2349.0}, {-1716.875, -1957.5}, {-1716.875, -1566.0}, {-1716.875, -1174.5}, {-1716.875, -783.0}, {-1716.875, -391.5}, {-1716.875, -0.0},
        {-1373.5, -3132.0}, {-1373.5, -2740.5}, {-1373.5, -2349.0}, {-1373.5, -1957.5}, {-1373.5, -1566.0}, {-1373.5, -1174.5}, {-1373.5, -783.0}, {-1373.5, -391.5}, {-1373.5, -0.0},
        {-1030.125, -3132.0}, {-1030.125, -2740.5}, {-1030.125, -2349.0}, {-1030.125, -1957.5}, {-1030.125, -1566.0}, {-1030.125, -1174.5}, {-1030.125, -783.0}, {-1030.125, -391.5}, {-1030.125, -0.0},
        {-686.75, -3132.0}, {-686.75, -2740.5}, {-686.75, -2349.0}, {-686.75, -1957.5}, {-686.75, -1566.0}, {-686.75, -1174.5}, {-686.75, -783.0}, {-686.75, -391.5}, {-686.75, -0.0},
        {-343.375, -3132.0}, {-343.375, -2740.5}, {-343.375, -2349.0}, {-343.375, -1957.5}, {-343.375, -1566.0}, {-343.375, -1174.5}, {-343.375, -783.0}, {-343.375, -391.5}, {-343.375, -0.0},
        {0.0, -3132.0}, {0.0, -2740.5}, {0.0, -2349.0}, {0.0, -1957.5}, {0.0, -1566.0}, {0.0, -1174.5}, {0.0, -783.0}, {0.0, -391.5}, {0.0, -0.0}
    };
   double rhs[2] = {0, 0};
   
   // First solve 
   if ((termcode = sym_solve(env_warm)) < 0){
      printf("WARM: PROBLEM INFEASIBLE!\n");
   }
   sym_build_dual_func(env_warm);
   // print_dual_function(env_warm);
   sym_evaluate_dual_function(env_warm, rhs, 0, &dualFuncObj);
   print_dual_function(env_warm);

   if (sym_is_proven_optimal(env_warm)){
      sym_get_obj_val(env_warm, &warmObjVal);
      printf(" DF: %.10f\n", dualFuncObj);
      printf("RVF: %.10f\n", warmObjVal);
      assert((fabs(dualFuncObj - warmObjVal) < 1e-5));
   } else if (sym_is_proven_primal_infeasible(env_warm)) {
      assert(dualFuncObj > 1e19);
      printf("RVF: INFEASIBLE\n");
      printf(" DF: %.10f\n", dualFuncObj);
   }

   printf("======================================\n");

   for (int i = 0; i < num_zetas; i++){
      rhs[0] = zetas[i][0];
      rhs[1] = zetas[i][1];
      printf("======================================\n");
      printf(" RHS: (%.2f, %.2f)\n", rhs[0], rhs[1]);
      printf("======================================\n");
      
      set_rhs(env_warm, rhs, num_objs);
    
      if ((termcode = sym_warm_solve(env_warm)) < 0){
         printf("WARM: PROBLEM INFEASIBLE!\n");
      }

      if (sym_is_proven_optimal(env_warm)){
         printf("  PROVEN OPTIMAL!\n");
         sym_get_obj_val(env_warm, &warmObjVal);
         printf("RVF WARM: %.10f\n", warmObjVal);
      } else if (sym_is_proven_primal_infeasible(env_warm)){
         printf("  PROVEN INFEASIBLE!\n");
      }

      sym_build_dual_func(env_warm);
      // print_dual_function(env_warm);
      sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);
      print_dual_function(env_warm);

      set_rhs(env_cold, rhs, num_objs);
      if ((termcode = sym_solve(env_cold)) < 0){
         printf("WARM: PROBLEM INFEASIBLE!\n");
      }
     
      if (sym_is_proven_optimal(env_cold)){
         sym_get_obj_val(env_cold, &warmObjVal);
         printf(" DF: %.10f\n", dualFuncObj);
         printf("RVF: %.10f\n", warmObjVal);
         if (fabs(dualFuncObj - warmObjVal) >= 1e-5){
            print_dual_function(env_warm);
            sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);
            print_dual_function(env_warm);
         }
         assert((fabs(dualFuncObj - warmObjVal) < 1e-5));
      } else if (sym_is_proven_primal_infeasible(env_cold)) {
         assert(dualFuncObj > 1e19);
         printf("RVF: INFEASIBLE\n");
         printf(" DF: %.10f\n", dualFuncObj);
      }
   }

   printf("LP TIME from DUAL FUNC: %.5f\n", sym_get_lp_time_dual_func(env_warm));

   // Memory clean-ups
   sym_close_environment(env_warm);
   sym_close_environment(env_cold);

   // free(zeta_lst);
   // free(rvf_lst);
   // free(df_lst);

   return 0;
}  



// #include <cstddef>
// #include <cstdio>
// #include <cstdlib>
// #include <cassert>
// #include <math.h>
// #include <time.h>

// #include "symphony.h"
// // feb223
// // These paths should be changed
// // I still don't figure out how to use makefile to 
// // include these files automatically
// #include "../include/sym_types.h"
// #include "../include/sym_master.h"

// void generate_rand_array(int m, int l, int u, double *res) {
//   for (int i = 0; i < m; i++) {
//     res[i] = (rand() % (u - l + 1)) + l;
//   }
// }

// void set_rhs(sym_environment *env, double *rhs, int m){
//    for (int i = 0; i < m; i++){
//       switch (env->mip->sense[i]){
//          case 'E' :
//             sym_set_row_upper(env, i, rhs[i]);
//             sym_set_row_lower(env, i, rhs[i]);
//             break;
//          case 'G':
//             sym_set_row_lower(env, i, rhs[i]);
//             break;
//          case 'L':
//             sym_set_row_upper(env, i, rhs[i]);
//             break;
//          default:
//             printf("Error\n");
//       }
//    }
// }

// double *linspace(double a, double b, int pieces){

//    double *res = (double *)malloc(sizeof(double) * pieces);
//    double step = fabs(a + b)/((double)pieces);

//    res[0] = a;
//    res[pieces - 1] = b;

//    for (int i = 1; i < pieces - 1; i++){
//       res[i] = res[i - 1] + step;
//    }

//    return res;
// }

// int main(int argc, char **argv)
// {   
   
//    int termcode;
//    int numTrain = 10;
//    int numTests = 10;
//    double warmObjVal, coldObjVal, dualFuncObj;

//    sym_environment *env_cold = sym_open_environment(); 
//    sym_environment *env_warm = sym_open_environment(); 

//    sym_parse_command_line(env_warm, argc, argv); 
//    sym_parse_command_line(env_cold, argc, argv); 

//    sym_load_problem(env_warm);
//    sym_load_problem(env_cold);

//    sym_set_int_param(env_warm, "verbosity", -2);
//    sym_set_int_param(env_cold, "verbosity", -2);

//    // feb223
//    // Those are the parameters to be set in order to 
//    // keep the branch-and-bound tree valid for RHS changes
//    // (E.g. Cuts and reduced cost fixing are not RHS-invariant)
//    sym_set_int_param(env_warm, "keep_warm_start", TRUE);
//    sym_set_int_param(env_warm, "keep_dual_function_description", TRUE);
//    sym_set_int_param(env_warm, "should_use_rel_br", FALSE);
//    sym_set_int_param(env_warm, "use_hot_starts", FALSE);
//    sym_set_int_param(env_warm, "should_warmstart_node", TRUE);
//    sym_set_int_param(env_warm, "sensitivity_analysis", TRUE);
//    sym_set_int_param(env_warm, "sensitivity_rhs", true);
//    sym_set_int_param(env_warm, "sensitivity_bounds", TRUE);
//    sym_set_int_param(env_warm, "set_obj_upper_lim", FALSE);
//    sym_set_int_param(env_warm, "do_primal_heuristic", FALSE);
//    sym_set_int_param(env_warm, "prep_level", -1);
//    sym_set_int_param(env_warm, "tighten_root_bounds", FALSE);
//    sym_set_int_param(env_warm, "max_sp_size", 100);
//    sym_set_int_param(env_warm, "do_reduced_cost_fixing", FALSE);
//    sym_set_int_param(env_warm, "generate_cgl_cuts", FALSE);
//    sym_set_int_param(env_warm, "max_active_nodes", 1);

//    sym_set_int_param(env_cold, "keep_warm_start", TRUE);
//    sym_set_int_param(env_cold, "keep_dual_function_description", TRUE);
//    sym_set_int_param(env_cold, "should_use_rel_br", FALSE);
//    sym_set_int_param(env_cold, "use_hot_starts", FALSE);
//    sym_set_int_param(env_cold, "should_warmstart_node", TRUE);
//    sym_set_int_param(env_cold, "sensitivity_analysis", TRUE);
//    sym_set_int_param(env_cold, "sensitivity_rhs", true);
//    sym_set_int_param(env_cold, "sensitivity_bounds", TRUE);
//    sym_set_int_param(env_cold, "set_obj_upper_lim", FALSE);
//    sym_set_int_param(env_cold, "do_primal_heuristic", FALSE);
//    sym_set_int_param(env_cold, "prep_level", -1);
//    sym_set_int_param(env_cold, "tighten_root_bounds", FALSE);
//    sym_set_int_param(env_cold, "max_sp_size", 100);
//    sym_set_int_param(env_cold, "do_reduced_cost_fixing", FALSE);
//    sym_set_int_param(env_cold, "generate_cgl_cuts", FALSE);
//    sym_set_int_param(env_cold, "max_active_nodes", 1);

//    //----------------------
//    // linspace numpy-like 
//    //----------------------
//    // double a = -16, b = 5, zerotol = 1e-7;
//    // int pieces = 2000;
//    // double *zeta_lst = linspace(a, b, pieces);
//    // double *zeta_lst_1 = linspace(-15, 5, 10);
//    // double *rvf_lst  = (double*)malloc(sizeof(double) * pieces);
//    // double *df_lst   = (double*)malloc(sizeof(double) * pieces);

//    int num_objs = 2;
//    double *rhs = (double*)malloc(sizeof(double) * num_objs);
   
//    // First solve 
//    if ((termcode = sym_solve(env_warm)) < 0){
//       printf("WARM: PROBLEM INFEASIBLE!\n");
//    }
//    sym_build_dual_func(env_warm);
//    print_dual_function(env_warm);
//    sym_evaluate_dual_function(env_warm, rhs, 0, &dualFuncObj);


//    printf("======================================\n");
//    rhs[0] = -2747.0;
//    rhs[1] = -3132.0;

//    set_rhs(env_warm, rhs, num_objs);
    
//    if ((termcode = sym_warm_solve(env_warm)) < 0){
//       printf("WARM: PROBLEM INFEASIBLE!\n");
//    }

//    if (sym_is_proven_optimal(env_warm)){
//       printf("  PROVEN OPTIMAL!\n");
//    } else if (sym_is_proven_primal_infeasible(env_warm)){
//       printf("  PROVEN INFEASIBLE!\n");
//    }

//    sym_build_dual_func(env_warm);
//    print_dual_function(env_warm);
//    sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

//    printf("======================================\n");
//    rhs[0] = -2747.0;
//    rhs[1] = -2740.5;

//    set_rhs(env_warm, rhs, num_objs);
    
//    if ((termcode = sym_warm_solve(env_warm)) < 0){
//       printf("WARM: PROBLEM INFEASIBLE!\n");
//    }

//    if (sym_is_proven_optimal(env_warm)){
//       printf("  PROVEN OPTIMAL!\n");
//    } else if (sym_is_proven_primal_infeasible(env_warm)){
//       printf("  PROVEN INFEASIBLE!\n");
//    }

//    sym_build_dual_func(env_warm);
//    print_dual_function(env_warm);
//    sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

//    printf("======================================\n");
//    rhs[0] = -2747.0;
//    rhs[1] = -2349.0;

//    set_rhs(env_warm, rhs, num_objs);
    
//    if ((termcode = sym_warm_solve(env_warm)) < 0){
//       printf("WARM: PROBLEM INFEASIBLE!\n");
//    }

//    if (sym_is_proven_optimal(env_warm)){
//       printf("  PROVEN OPTIMAL!\n");
//    } else if (sym_is_proven_primal_infeasible(env_warm)){
//       printf("  PROVEN INFEASIBLE!\n");
//    }

//    sym_build_dual_func(env_warm);
//    print_dual_function(env_warm);
//    sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

//    printf("======================================\n");
//    rhs[0] = -2747.0;
//    rhs[1] = -1957.5;

//    set_rhs(env_warm, rhs, num_objs);
    
//    if ((termcode = sym_warm_solve(env_warm)) < 0){
//       printf("WARM: PROBLEM INFEASIBLE!\n");
//    }

//    if (sym_is_proven_optimal(env_warm)){
//       printf("  PROVEN OPTIMAL!\n");
//    } else if (sym_is_proven_primal_infeasible(env_warm)){
//       printf("  PROVEN INFEASIBLE!\n");
//    }

//    sym_build_dual_func(env_warm);
//    print_dual_function(env_warm);
//    sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

//    printf("======================================\n");
//    rhs[0] = -2747.0;
//    rhs[1] = -1566.0;

//    set_rhs(env_warm, rhs, num_objs);
    
//    if ((termcode = sym_warm_solve(env_warm)) < 0){
//       printf("WARM: PROBLEM INFEASIBLE!\n");
//    }

//    if (sym_is_proven_optimal(env_warm)){
//       printf("  PROVEN OPTIMAL!\n");
//    } else if (sym_is_proven_primal_infeasible(env_warm)){
//       printf("  PROVEN INFEASIBLE!\n");
//    }

//    sym_build_dual_func(env_warm);
//    print_dual_function(env_warm);
//    sym_evaluate_dual_function(env_warm, rhs, 2, &dualFuncObj);

//    // // check_dual_solutions(env_cold->mip, env_warm->warm_start->dual_func);

//    // rhs[0] = -55.5;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;
      
//    // // Set the new RHS
//    // set_rhs(env_warm, rhs, m);

//    // if ((termcode = sym_warm_solve(env_warm)) < 0){
//    //    printf("PROBLEM INFEASIBLE!\n");
//    // }

//    // sym_get_obj_val(env_warm, &warmObjVal);
//    // printf("WARM OBJ : %.5f\n", warmObjVal);

//    // sym_build_dual_func(env_warm);

//    // rhs[0] = 40.0/9.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -55.5;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -73.0/6.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -10;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -73.0/6.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;
      
//    // // Set the new RHS
//    // set_rhs(env_warm, rhs, m);

//    // if ((termcode = sym_warm_solve(env_warm)) < 0){
//    //    printf("PROBLEM INFEASIBLE!\n");
//    // }

//    // sym_get_obj_val(env_warm, &warmObjVal);
//    // printf("WARM OBJ : %.5f\n", warmObjVal);

//    // sym_build_dual_func(env_warm);

//    // rhs[0] = 40.0/9.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -55.5;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -73.0/6.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -10;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);

//    // rhs[0] = -10;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;
      
//    // // Set the new RHS
//    // set_rhs(env_warm, rhs, m);

//    // if ((termcode = sym_warm_solve(env_warm)) < 0){
//    //    printf("PROBLEM INFEASIBLE!\n");
//    // }

//    // sym_get_obj_val(env_warm, &warmObjVal);
//    // printf("WARM OBJ : %.5f\n", warmObjVal);

//    // sym_build_dual_func(env_warm);

//    // rhs[0] = 40.0/9.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
//    // printf("Dual function evaluates %.7f\n", dualFuncObj);

//    // rhs[0] = -55.5;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
//    // printf("Dual function evaluates %.7f\n", dualFuncObj);

//    // rhs[0] = -73.0/6.0;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
//    // printf("Dual function evaluates %.7f\n", dualFuncObj);

//    // rhs[0] = -10;
//    // rhs[1] = 4;
//    // rhs[2] = 5;
//    // rhs[3] = 5;

//    // sym_evaluate_dual_function(env_warm, rhs, m, &dualFuncObj);
//    // printf("Dual function evaluates %.7f\n", dualFuncObj);


//    // Memory clean-ups
//    sym_close_environment(env_warm);
//    sym_close_environment(env_cold);

//    // free(zeta_lst);
//    // free(rvf_lst);
//    // free(df_lst);

//    return 0;
// }  

