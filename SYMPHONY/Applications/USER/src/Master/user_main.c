/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/*===========================================================================*/

#define CALL_FUNCTION(f) \
if ((termcode = f) < 0){                                                    \
   printf("Error detected: termcode = %i\n", termcode);                     \
   printf("Exiting...\n\n");                                                \
   exit(termcode);                                                          \
}

/*===========================================================================*\
   This file contains the main() for the master process.

   Note that, if you want to use the OSI SYMPHONY interface, you should set the
   USE_OSI_INTERFACE flag and define the COINROOT path in the SYMPHONY 
   Makefile. Otherwise, the C callable library functions will be used by 
   default. See below for the usage.
\*===========================================================================*/

#if defined(USE_OSI_INTERFACE)

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{
   OsiSymSolverInterface si;

   /* Parse the command line */
   si.parseCommandLine(argc, argv);
   
   /* Read in the problem */
   si.loadProblem();

   /* Find a priori problem bounds */
   si.findInitialBounds();

   /* Solve the problem */
   si.branchAndBound();
   
   return(0);
}

#else

#include "symphony.h"
#include "sym_master.h"
#include <stdlib.h>

#include "user.h"

#define USER_FUNC_ERROR -1

int user_test(sym_environment *env);

int main(int argc, char **argv)
{

   int termcode;

   /* Create a SYMPHONY environment */
   sym_environment *env = sym_open_environment();

   /* Print version info */
   sym_version();
   
   if (!env){
      printf("Error initializing environement\n");
      exit(0);
   }

   /* Create the data structure for storing the problem instance.*/
   user_problem *prob = (user_problem *)calloc(1, sizeof(user_problem));
   prob->mip = (MIPdesc *)calloc(1, sizeof(MIPdesc));
   
   CALL_FUNCTION( sym_set_user_data(env, (void *)prob) );

   CALL_FUNCTION( sym_parse_command_line(env, argc, argv) );

   CALL_FUNCTION( user_read_data(prob, prob->par.infile) );
   
   CALL_FUNCTION( user_load_problem(env, prob) );

   CALL_FUNCTION( sym_solve(env) );
     
   // Added by Suresh for debugging
   /*
   printf("\n\n");
   double obj_val = 0.0;
   int n = env->best_sol.xlength;
   int i = 0;
   for (i = 0; i < n; i++) {
      if (env->best_sol.xind[i] < prob->origvar_num) {
         obj_val += env->best_sol.xval[i]*prob->origobj_coeffs[env->best_sol.xind[i]];
      } else {
         continue;
      }
   }
   printf("\n\nORIGINAL OBJECTIVE FUNCTION VALUE = %06f\n\n", obj_val);
   */


   CALL_FUNCTION( sym_close_environment(env) );
   
   return(0);
   
}

/*===========================================================================*\
\*===========================================================================*/

int user_read_data(user_problem *prob, char *infile)
{
   int j, counter;
   CoinMpsIO mps;

   mps.messageHandler()->setLogLevel(0);
   
   mps.setInfinity(mps.getInfinity()); // TODO: What exactly is this doing here?

   if (mps.readMps(infile,"")){
      return(USER_FUNC_ERROR);
   }
   
   prob->mip->m  = mps.getNumRows();
   prob->mip->n  = mps.getNumCols();
   prob->mip->obj_sense = 1.0; // TODO: Remove this 'minimize' assumption.
   prob->infty = mps.getInfinity();
   
   prob->mip->obj    = (double *) malloc(DSIZE * prob->mip->n);
   prob->mip->rhs    = (double *) malloc(DSIZE * prob->mip->m);
   prob->mip->sense  = (char *)   malloc(CSIZE * prob->mip->m);
   prob->mip->rngval = (double *) malloc(DSIZE * prob->mip->m);
   prob->mip->ub     = (double *) malloc(DSIZE * prob->mip->n);
   prob->mip->lb     = (double *) malloc(DSIZE * prob->mip->n);
   prob->mip->colname = (char **) malloc(sizeof(char *) * prob->mip->n);  
   prob->mip->rowname = (char **) malloc(sizeof(char *) * prob->mip->m);  
   
   memcpy(prob->mip->obj, const_cast <double *> (mps.getObjCoefficients()),
	  DSIZE * prob->mip->n);
   memcpy(prob->mip->rhs, const_cast <double *> (mps.getRightHandSide()),
	  DSIZE * prob->mip->m);
   memcpy(prob->mip->sense, const_cast <char *> (mps.getRowSense()),
	  CSIZE * prob->mip->m);
   memcpy(prob->mip->rngval, const_cast <double *> (mps.getRowRange()),
	  DSIZE * prob->mip->m);
   memcpy(prob->mip->ub, const_cast <double *> (mps.getColUpper()),
	  DSIZE * prob->mip->n);
   memcpy(prob->mip->lb, const_cast <double *> (mps.getColLower()),
	  DSIZE * prob->mip->n);

   // Save names
   for (j = 0; j < prob->mip->n; j++){
      prob->mip->colname[j] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      strncpy(prob->mip->colname[j], mps.columnName(j), MAX_NAME_SIZE);
      prob->mip->colname[j][MAX_NAME_SIZE-1] = 0;
   }

   for (j = 0; j < prob->mip->m; j++){
      prob->mip->rowname[j] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      strncpy(prob->mip->rowname[j], mps.rowName(j), MAX_NAME_SIZE);
      prob->mip->rowname[j][MAX_NAME_SIZE-1] = 0;
   }

   // Added by Suresh for debugging
   /*
   prob->origvar_num = prob->mip->n;
   prob->origobj_coeffs = (double *) malloc(DSIZE * prob->mip->n);
   memcpy(prob->origobj_coeffs, const_cast <double *> (mps.getObjCoefficients()),
	  DSIZE * prob->mip->n);
     */

   //user defined matind, matval, matbeg--fill as column ordered
   const CoinPackedMatrix * matrixByCol= mps.getMatrixByCol();
   
   prob->mip->matbeg = (int *) malloc(ISIZE * (prob->mip->n + 1));
   memcpy(prob->mip->matbeg, const_cast<int *>(matrixByCol->getVectorStarts()),
	  ISIZE * (prob->mip->n + 1));
   
   prob->mip->matval = (double *) malloc(DSIZE*prob->mip->matbeg[prob->mip->n]);
   prob->mip->matind = (int *)    malloc(ISIZE*prob->mip->matbeg[prob->mip->n]);
   
   memcpy(prob->mip->matval, const_cast<double *> (matrixByCol->getElements()),
	  DSIZE * prob->mip->matbeg[prob->mip->n]);
   memcpy(prob->mip->matind, const_cast<int *> (matrixByCol->getIndices()), 
	  ISIZE * prob->mip->matbeg[prob->mip->n]);

   //user defined matind_row, matval_row, matbeg_row--fill as row ordered
   const CoinPackedMatrix * matrixByRow= mps.getMatrixByRow();
   
   prob->matbeg_row = (int *) malloc(ISIZE * (prob->mip->m + 1));
   memcpy(prob->matbeg_row, const_cast<int *>(matrixByRow->getVectorStarts()),
	  ISIZE * (prob->mip->m + 1));
   
   prob->matval_row = (double *) malloc(DSIZE*prob->matbeg_row[prob->mip->m]);
   prob->matind_row = (int *)    malloc(ISIZE*prob->matbeg_row[prob->mip->m]);
   
   memcpy(prob->matval_row, const_cast<double *> (matrixByRow->getElements()),
	  DSIZE * prob->matbeg_row[prob->mip->m]);
   memcpy(prob->matind_row, const_cast<int *> (matrixByRow->getIndices()), 
	  ISIZE * prob->matbeg_row[prob->mip->m]);
  
   //count number of different type constraints.
   prob->con_sense_e = 0;
   prob->con_sense_l = 0;
   prob->con_sense_g = 0;
   prob->ubinfty = 0;
   prob->lbinfty = 0;
   prob->infubind     = (int *) calloc(prob->mip->n, ISIZE);
   prob->inflbind     = (int *) calloc(prob->mip->n, ISIZE);
   prob->infubsofar   = (int *) malloc(ISIZE * prob->mip->n);
   prob->inflbsofar   = (int *) malloc(ISIZE * prob->mip->n);
   for (j = 0; j < prob->mip->m; j++) {
      if (prob->mip->sense[j] == 'E') {
         prob->con_sense_e++;
      } else if (prob->mip->sense[j] == 'L') {
         prob->con_sense_l++;
      } else if (prob->mip->sense[j] == 'G') {
         prob->con_sense_g++;
      } else {
         printf("\nNOOOOOOO!! ERROR!! ERROR!!\n");
      }
   }

   // count number of infinity UBs and -infinity LBs, and also their indicators.
   for (j = 0; j < prob->mip->n; j++) {
      prob->infubsofar[j] = prob->ubinfty;
      if (prob->mip->ub[j] >= prob->infty) {
         prob->ubinfty++;
         prob->infubind[j] = 1;
      }
      prob->inflbsofar[j] = prob->lbinfty;
      if (prob->mip->lb[j] <= -prob->infty) {
         prob->lbinfty++;
         prob->inflbind[j] = 1;
      }
   }

   prob->tempub     = (double *) malloc(DSIZE * (prob->mip->n - prob->ubinfty));
   prob->templb     = (double *) malloc(DSIZE * (prob->mip->n - prob->lbinfty));

   // Fill temporary UB arrays
   counter = 0;
   for (j = 0; j < prob->mip->n; j++) {
      if (prob->mip->ub[j] >= prob->infty) {
         continue;
      } else {
         prob->tempub[counter] = prob->mip->ub[j];
         counter++;
      }
   }

   // Fill temporary LB arrays
   counter = 0;
   for (j = 0; j < prob->mip->n; j++) {
      if (prob->mip->lb[j] <= -prob->infty) {
         continue;
      } else {
         prob->templb[counter] = prob->mip->lb[j];
         counter++;
      }
   }

   /* Suresh: debug code START */
   /*
   printf("\nORIG_OBJ\n");
   for (j = 0; j < prob->mip->n; j++) {
      printf("%s\tx%d\t%2.4e\n", prob->mip->colname[j], j, prob->mip->obj[j]);
   }
   printf("\nORIG_LB\n");
   for (j = 0; j < prob->mip->n; j++) {
      printf("%s\tx%d\t%2.4e\n", prob->mip->colname[j], j, prob->mip->lb[j]);
   }
   printf("\nORIG_UB\n");
   for (j = 0; j < prob->mip->n; j++) {
      printf("%s\tx%d\t%2.4e\n", prob->mip->colname[j], j, prob->mip->ub[j]);
   }
   printf("\nORIG_SENSE\n");
   for (j = 0; j < prob->mip->m; j++) {
      printf("%s\tc%d\t%c\n", prob->mip->rowname[j], j, prob->mip->sense[j]);
   }
   printf("\nORIG_RHS\n");
   for (j = 0; j < prob->mip->m; j++) {
      printf("%s\tc%d\t%2.4e\n", prob->mip->rowname[j], j, prob->mip->rhs[j]);
   }
   */
   /*
   int i;
   printf("\nORIG_COEFF\n");
   for (i = 0; i < prob->mip->n; i++) {
      if ((prob->mip->matbeg[i+1] - prob->mip->matbeg[i]) > 0) {
         for (j = prob->mip->matbeg[i]; j < prob->mip->matbeg[i+1]; j++) {
            printf("x%d\tc%d\t%2.4e\n", i, prob->mip->matind[j], prob->mip->matval[j]);
         }
      }
   }
   */
   /* Suresh: debug code END */

   return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*\
\*===========================================================================*/

int user_load_problem(sym_environment *env, user_problem *prob){
   
   int i, j, index, index1, n, m, nz, nz_index = 0, *column_starts, *matrix_indices;
   double *matrix_values, *lb, *ub, *obj, *rhs, *rngval;
   char *sense, *is_int, obj_sense = 1.0;
   
   // number of upper level variables
   int num_upperlevel_vars = 1;
   // TODO: Suresh: following two names are misleading. Fix them later.
   // number of nonzeros in lower level constraint coeffs corresponding to upper level variables
   int nz_upperlevel = prob->mip->matbeg[num_upperlevel_vars];
   // number of nonzeros in lower level constraint coeffs corresponding to lower level variables
   int nz_lowerlevel = prob->mip->matbeg[prob->mip->n] - nz_upperlevel;
   // number of lower level finite upper bound constraints on variables
   int num_lowerlevel_ubcons = (prob->mip->n - num_upperlevel_vars - (prob->ubinfty - prob->infubsofar[num_upperlevel_vars]));
   // number of lower level dual variables on finite upper bound constraints
   int num_lowerlevel_dual_ubvars = num_lowerlevel_ubcons;
   // number of lower level finite lower bound constrains on variables
   int num_lowerlevel_lbcons = (prob->mip->n - num_upperlevel_vars - (prob->lbinfty - prob->inflbsofar[num_upperlevel_vars]));
   // number of lower level dual variables on finite lower bound constraints
   int num_lowerlevel_dual_lbvars = num_lowerlevel_lbcons;
   // number of nonzeros per row in lower level part of coefficient matrix
   int *nz_lowerlevel_row;

   /* set up the inital LP data */
   n = prob->mip->n + prob->mip->m + num_lowerlevel_dual_ubvars + num_lowerlevel_dual_lbvars;
   m = prob->mip->m + num_lowerlevel_ubcons + num_lowerlevel_lbcons + (prob->mip->n - num_upperlevel_vars);
   nz = prob->mip->matbeg[prob->mip->n] + 2*num_lowerlevel_ubcons + 2*num_lowerlevel_lbcons + nz_lowerlevel;
   prob->colnum = n;
   prob->rownum = m;

   /* Allocate the arrays */
   column_starts  = (int *) malloc((n + 1) * ISIZE);
   matrix_indices = (int *) malloc((nz) * ISIZE);
   matrix_values  = (double *) malloc((nz) * DSIZE);
   obj            = (double *) malloc(n * DSIZE);
   lb             = (double *) malloc(n * DSIZE);
   ub             = (double *) malloc(n * DSIZE);
   rhs            = (double *) malloc(m * DSIZE);
   sense          = (char *) malloc(m * CSIZE);
   rngval         = (double *) calloc(m, DSIZE); /* TODO:Correct the assumption that this is zero always */
   is_int         = (char *) malloc(n * CSIZE);
   
   /* Fill out the appropriate data structures */
   if (prob->mip->obj_sense > 0.0) {
      for (i = 0; i < n; i++) {
         if (i < prob->mip->n) {
            obj[i] = prob->mip->obj[i];
         } else {
            obj[i] = 0.0;
         }
      }
   } else {
      for (i = 0; i < n; i++) {
         if (i < prob->mip->n) {
            obj[i] = -prob->mip->obj[i];
         } else {
            obj[i] = 0.0;
         }
      }
   }
   env->mip->colname = (char **) malloc(sizeof(char *) * n);  
   env->mip->rowname = (char **) malloc(sizeof(char *) * m);  

   /* The original upper level variables */
   for (i = 0; i < num_upperlevel_vars; i++) {
      is_int[i] = FALSE;
      ub[i] = prob->mip->ub[i];
      lb[i] = prob->mip->lb[i];
      env->mip->colname[i] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      strncpy(env->mip->colname[i], prob->mip->colname[i], MAX_NAME_SIZE);
      env->mip->colname[i][MAX_NAME_SIZE-1] = 0;
   }
   /* The original lower level variables */
   for (; i < prob->mip->n; i++) {
      is_int[i] = FALSE;
      ub[i] = prob->infty;
      lb[i] = -prob->infty;
      env->mip->colname[i] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      strncpy(env->mip->colname[i], prob->mip->colname[i], MAX_NAME_SIZE);
      env->mip->colname[i][MAX_NAME_SIZE-1] = 0;
   }
   /* The duals on original lower level constraints */
   for (;i < prob->mip->n + prob->mip->m; i++){
      is_int[i] = FALSE;
      if (prob->mip->sense[i - prob->mip->n] == 'L') {
         ub[i] = 0;
         lb[i] = -prob->infty;
      } else if (prob->mip->sense[i - prob->mip->n] == 'G') {
         ub[i] = prob->infty;
         lb[i] = 0;
      } else {
         ub[i] = prob->infty;
         lb[i] = -prob->infty;
      }
      env->mip->colname[i] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      env->mip->colname[i][0] = 'D';
      env->mip->colname[i][1] = '_';
      strncpy(env->mip->colname[i]+2, prob->mip->rowname[i - prob->mip->n],
            MAX_NAME_SIZE-2);
      env->mip->colname[i][MAX_NAME_SIZE-1] = 0;
   }
   /* The duals on lower level variable upper bound constraints */
   for (j = num_upperlevel_vars; j < prob->mip->n; j++) {
      if (prob->mip->ub[j] >= prob->infty) {
         continue;
      }
      index = prob->mip->n + prob->mip->m + ((j - num_upperlevel_vars) - (prob->infubsofar[j] - prob->infubsofar[num_upperlevel_vars]));
      is_int[index] = FALSE;
      ub[index] = 0.0;
      lb[index] = -prob->infty;
      env->mip->colname[index] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      env->mip->colname[index][0] = 'U';
      env->mip->colname[index][1] = '_';
      strncpy(env->mip->colname[index]+2, prob->mip->colname[j],
            MAX_NAME_SIZE-2);
      env->mip->colname[index][MAX_NAME_SIZE-1] = 0;
      i++;
   }

   /* The duals on lower level variable lower bound constraints */
   for (j = num_upperlevel_vars; j < prob->mip->n; j++) {
      if (prob->mip->lb[j] <= -prob->infty) {
         continue;
      }
      index = prob->mip->n + prob->mip->m + num_lowerlevel_dual_ubvars + ((j - num_upperlevel_vars) - (prob->inflbsofar[j] - prob->inflbsofar[num_upperlevel_vars]));
      is_int[index] = FALSE;
      ub[index] = prob->infty;
      lb[index] = 0.0;
      env->mip->colname[index] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      env->mip->colname[index][0] = 'L';
      env->mip->colname[index][1] = '_';
      strncpy(env->mip->colname[index]+2, prob->mip->colname[j],
            MAX_NAME_SIZE-2);
      env->mip->colname[index][MAX_NAME_SIZE-1] = 0;
      i++;
   }

   /* The original lower level constraints */
   for (i = 0; i < prob->mip->m; i++) {
      sense[i] = prob->mip->sense[i];
      rhs[i] = prob->mip->rhs[i];
      env->mip->rowname[i] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      strncpy(env->mip->rowname[i], prob->mip->rowname[i], MAX_NAME_SIZE);
      env->mip->rowname[i][MAX_NAME_SIZE-1] = 0;
   }
   /* The lower level upper bound constraints on variables */
   for (j = num_upperlevel_vars; j < prob->mip->n; j++) {
      if (prob->mip->ub[j] >= prob->infty) {
         continue;
      }
      index = prob->mip->m + ((j - num_upperlevel_vars) - (prob->infubsofar[j] - prob->infubsofar[num_upperlevel_vars]));
      sense[index] = 'L';
      rhs[index] = prob->mip->ub[j];
      env->mip->rowname[index] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      env->mip->rowname[index][0] = 'U';
      env->mip->rowname[index][1] = '_';
      strncpy(env->mip->rowname[index]+2, prob->mip->colname[j],
            MAX_NAME_SIZE-2);
      env->mip->rowname[index][MAX_NAME_SIZE-1] = 0;
      i++;
   }

   /* The lower level lower bound constraints on variables */
   for (j = num_upperlevel_vars; j < prob->mip->n; j++) {
      if (prob->mip->lb[j] <= -prob->infty) {
         continue;
      }
      index = prob->mip->m + num_lowerlevel_ubcons + ((j - num_upperlevel_vars) - (prob->inflbsofar[j] - prob->inflbsofar[num_upperlevel_vars]));
      sense[index] = 'G';
      rhs[index] = prob->mip->lb[j];
      env->mip->rowname[index] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      env->mip->rowname[index][0] = 'L';
      env->mip->rowname[index][1] = '_';
      strncpy(env->mip->rowname[index]+2, prob->mip->colname[j],
            MAX_NAME_SIZE-2);
      env->mip->rowname[index][MAX_NAME_SIZE-1] = 0;
      i++;
   }

   /* The lower level dual constraints */
   for (j = num_upperlevel_vars; j < prob->mip->n; j++) {
      index = prob->mip->m + num_lowerlevel_ubcons + num_lowerlevel_lbcons + ((j - num_upperlevel_vars));
      sense[index] = 'E';
      // Suresh: zero'd rhs for debugging to check if entire code works properly
//      rhs[index] = 0;
      rhs[index] = -prob->mip->obj[j];
      env->mip->rowname[index] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
      env->mip->rowname[index][0] = 'E';
      env->mip->rowname[index][1] = '_';
      strncpy(env->mip->rowname[index]+2, prob->mip->colname[j],
	      MAX_NAME_SIZE-2);
      env->mip->rowname[index][MAX_NAME_SIZE-1] = 0;
      i++;
   }

   /* column sparse format entries corresponding to original upper level variables */
   column_starts[0] = 0;
   for (i = 0; i < num_upperlevel_vars; i++) {
      column_starts[i+1] = column_starts[i] + prob->mip->matbeg[i+1] - prob->mip->matbeg[i];
      if (prob->mip->matbeg[i + 1] - prob->mip->matbeg[i] > 0) {
         for (j = (prob->mip->matbeg[i]); j < (prob->mip->matbeg[i + 1]); j++) {
            matrix_values[nz_index] = prob->mip->matval[j];
            matrix_indices[nz_index] = prob->mip->matind[j];
            nz_index++;
         }
      }
   }
   /* column sparse format entries corresponding to original lower level variables */
   for (; i < prob->mip->n; i++) {
      column_starts[i+1] = column_starts[i] + prob->mip->matbeg[i+1] - prob->mip->matbeg[i] + 2 - prob->infubind[i] - prob->inflbind[i];
      if (prob->mip->matbeg[i + 1] - prob->mip->matbeg[i] > 0) {
         for (j = (prob->mip->matbeg[i]); j < (prob->mip->matbeg[i + 1]); j++) {
            matrix_values[nz_index] = prob->mip->matval[j];
            matrix_indices[nz_index] = prob->mip->matind[j];
            nz_index++;
         }
      }
      if (!prob->infubind[i]) {
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + ((i - num_upperlevel_vars) - (prob->infubsofar[i] - prob->infubsofar[num_upperlevel_vars]));
         nz_index++;
      }
      if (!prob->inflbind[i]) {
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + num_lowerlevel_ubcons + ((i - num_upperlevel_vars) - (prob->inflbsofar[i] - prob->inflbsofar[num_upperlevel_vars]));
         nz_index++;
      }
   }
   /* column sparse format entries corresponding to dual variables on original lower level constraints */
   // At first, find number of nonzeros per row in the lower level coefficient part of matrix
   nz_lowerlevel_row = (int *) malloc((nz_lowerlevel) * ISIZE);
   for (i = 0; i < prob->mip->m; i++) {
      nz_lowerlevel_row[i] = prob->matbeg_row[i+1] - prob->matbeg_row[i];
      if ((prob->matbeg_row[i+1] - prob->matbeg_row[i]) > 0) {
         for (j = prob->matbeg_row[i]; j < prob->matbeg_row[i+1]; j++) {
            if (prob->matind_row[j] < num_upperlevel_vars) {
               nz_lowerlevel_row[i]--;
            }
         }
      }
   }
   // Now update column sparse format entries as required
   index = 0;
   for (i = 0; i < prob->mip->m; i++) {
      column_starts[prob->mip->n + 1 + index] = column_starts[prob->mip->n + index] + nz_lowerlevel_row[i];
      index++;
      for (j = prob->matbeg_row[i]; j < prob->matbeg_row[i+1]; j++) {
         if (prob->matind_row[j] >= num_upperlevel_vars) {
            matrix_values[nz_index] = prob->matval_row[j];
            matrix_indices[nz_index] = prob->mip->m + num_lowerlevel_ubcons + num_lowerlevel_lbcons + (prob->matind_row[j] - num_upperlevel_vars);
            nz_index++;
         }
      }
   }
   /* column sparse format entries corresponding to dual variables on lower level primal UB constraints */
   for (i = num_upperlevel_vars; i < prob->mip->n; i++) {
      if (!prob->infubind[i]) {
         column_starts[prob->mip->n + 1 + index] = column_starts[prob->mip->n + index] + 1;
         index++;
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + num_lowerlevel_ubcons + num_lowerlevel_lbcons + (i - num_upperlevel_vars);
         nz_index++;
      }
   }
   /* column sparse format entries corresponding to dual variables on lower level primal LB constraints */
   for (i = num_upperlevel_vars; i < prob->mip->n; i++) {
      if (!prob->inflbind[i]) {
         column_starts[prob->mip->n + 1 + index] = column_starts[prob->mip->n + index] + 1;
         index++;
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + num_lowerlevel_ubcons + num_lowerlevel_lbcons + (i - num_upperlevel_vars);
         nz_index++;
      }
   }

   /* Assign memory and fill complementarity constraint indices for variables */
   prob->ccind          =   (int *) malloc(n * ISIZE);
   index = 0;
   index1 = 0;
   for (i = 0; i < n;) {
      if (i < prob->mip->n) {
         prob->ccind[i] = -1;
         i++;
      } else if (i < prob->mip->n + prob->mip->m) {
         if (sense[i - prob->mip->n] == 'E') {
            prob->ccind[i] = -1;
            i++;
            index1++;
         } else {
            prob->ccind[i] = index1;
            i++;
            index1++;
         }
      } else {
         prob->ccind[i] = index1;
         i++;
         index1++;
      }
   }

   /* Load the problem to SYMPHONY */   
   sym_explicit_load_problem(env, n, m, column_starts, matrix_indices,
			     matrix_values, lb, ub, is_int, obj, 0, sense, rhs,
			     rngval, true);

   /* Suresh: debug code START */
   /*
   printf("\nNEW_OBJ\n");
   for (j = 0; j < n; j++) {
      printf("%s\tx%d\t%f\n", env->mip->colname[j], j, obj[j]);
   }
   */
   /*
   printf("\nNEW_LB\n");
   for (j = 0; j < n; j++) {
      printf("%s\tx%d\t%2.4e\n", env->mip->colname[j], j, lb[j]);
   }
   printf("\nNEW_UB\n");
   for (j = 0; j < n; j++) {
      printf("%s\tx%d\t%2.4e\n", env->mip->colname[j], j, ub[j]);
   }
   printf("\nNEW_SENSE\n");
   for (j = 0; j < m; j++) {
      printf("%s\tc%d\t%c\n", env->mip->rowname[j], j, sense[j]);
   }
   */
   /*
   printf("\nNEW_RHS\n");
   for (j = 0; j < m; j++) {
      printf("%s\tc%d\t%f\n", env->mip->rowname[j], j, rhs[j]);
   }
   */
   /*
   printf("\nNEW_COEFF\n");
   for (i = 0; i < n; i++) {
      if ((column_starts[i+1] - column_starts[i]) > 0) {
         for (j = column_starts[i]; j < column_starts[i+1]; j++) {
            printf("x%d\tc%d\t%f\n", i, matrix_indices[j], matrix_values[j]);
         }
      }
   }
   */
   /*
   printf("\nNEW_ROW_NAMES\n");
   for(i = 0; i < m; i++) {
      printf("%s\tc%d\n", env->mip->rowname[i], i);
   }
   printf("\nNEW_COL_NAMES\n");
   for(i = 0; i < n; i++) {
      printf("%s\tx%d\n", env->mip->colname[i], i);
   }
   */
   /*
   printf("\nNEW_CCIND\n");
   for (i = 0; i < n; i++) {
      printf("%s\tx%d\t%d\n", env->mip->colname[i], i, prob->ccind[i]);
   }
   */
   /* Suresh: debug code END */

   /* Change prob->mip values to final problem values */
   prob->mip->n = n;
   prob->mip->m = m;
   prob->mip->nz = nz;
   prob->mip->obj_sense = obj_sense;

   prob->mip->obj    = (double *) realloc(prob->mip->obj, DSIZE * prob->mip->n);
   prob->mip->rhs    = (double *) realloc(prob->mip->rhs, DSIZE * prob->mip->m);
   prob->mip->sense  = (char *)   realloc(prob->mip->sense, CSIZE * prob->mip->m);
   prob->mip->rngval = (double *) realloc(prob->mip->rngval, DSIZE * prob->mip->m);
   prob->mip->ub     = (double *) realloc(prob->mip->ub, DSIZE * prob->mip->n);
   prob->mip->lb     = (double *) realloc(prob->mip->lb, DSIZE * prob->mip->n);
   prob->mip->is_int = (char *)   malloc(CSIZE * prob->mip->n);
   /* Default values for vvind, vvnum, feasible and rowact */
   prob->feasible    = USER__DO_NOT_BRANCH;
   prob->vvind       = (int *)    calloc(prob->mip->n, ISIZE);
   prob->vvnum       = 0;
   prob->rowact      = (double *) calloc(prob->mip->m, DSIZE);
   
   memcpy(prob->mip->obj, obj, DSIZE * prob->mip->n);
   memcpy(prob->mip->rhs, rhs, DSIZE * prob->mip->m);
   memcpy(prob->mip->sense, sense, CSIZE * prob->mip->m);
   memset(prob->mip->rngval, 0, DSIZE * m);                     // TODO: Fix this assumption.
   memcpy(prob->mip->ub, ub, DSIZE * prob->mip->n);
   memcpy(prob->mip->lb, lb, DSIZE * prob->mip->n);
   memcpy(prob->mip->is_int, is_int, CSIZE * prob->mip->n);

   prob->mip->matbeg = (int *) realloc(prob->mip->matbeg, ISIZE * (prob->mip->n + 1));
   memcpy(prob->mip->matbeg, column_starts, ISIZE * (prob->mip->n + 1));
   
   prob->mip->matval = (double *) realloc(prob->mip->matval, DSIZE*prob->mip->matbeg[prob->mip->n]);
   prob->mip->matind = (int *)    realloc(prob->mip->matind, ISIZE*prob->mip->matbeg[prob->mip->n]);
   
   memcpy(prob->mip->matval, matrix_values, DSIZE * prob->mip->matbeg[prob->mip->n]);
   memcpy(prob->mip->matind, matrix_indices, ISIZE * prob->mip->matbeg[prob->mip->n]);

   /* TODO: Delete mat*_row vectors here? */
   FREE(column_starts);
   FREE(matrix_indices);
   FREE(matrix_values);
   FREE(lb);
   FREE(ub);
   FREE(obj);
   FREE(sense);
   FREE(rhs);
   FREE(rngval);
   FREE(is_int);
   FREE(nz_lowerlevel_row);
   /* TODO: Is it good to free tempub, templb, infubind, inflbind, ifubsofar, inflbsofar here? */
   FREE(prob->tempub);
   FREE(prob->templb);
   FREE(prob->inflbsofar);
   FREE(prob->infubsofar);
   FREE(prob->inflbind);
   FREE(prob->infubind);
   FREE(prob->matval_row);
   FREE(prob->matind_row);
   FREE(prob->matbeg_row);

   return (FUNCTION_TERMINATED_NORMALLY);
}

#endif
