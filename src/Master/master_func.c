/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* This file was developed by Menal Guzelsoy for the SYMPHONY OSI interface. */
/*                                                                           */
/* (c) Copyright 2000-2022 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __PVM__
#include <pvmtev.h>
#endif

#include "sym_messages.h"
#include "sym_master.h"
#include "sym_master_u.h"
#include "sym_macros.h"
#include "sym_pack_cut.h"
#include "sym_pack_array.h"
#include "sym_lp_solver.h"
#include "sym_primal_heuristics.h"
// #include "sym_lp.h"
#include "sym_tm.h"

// #define DEBUG_DUAL_FUNC
// #define CHECK_DUAL_FUNC

/*===========================================================================*/
/*===========================================================================*/

int resolve_node(sym_environment *env, bc_node *node)
{
	node_desc *desc = &node->desc;
	LPdata *lp_data = (LPdata *)calloc(1, sizeof(LPdata));
	lp_sol *best_sol = &(env->best_sol);
	branch_desc *bpath;
	branch_obj *bobj;
	bc_node **path, *n;
	int level = node->bc_level;

	int *matbeg = 0, *matind = 0, size, nzcnt = 0, total_nzcnt = 0;
	int return_value, iterd = 0;
	double *matval = 0, colsol;
	int i, j = 0, cnt = 0, *xind = 0;
	double *xval = 0, *rhs = 0, lpetol = 9.9999999999999995e-07;
	char *sense = 0;
	cut_data *cut;

	/*------------------------------------------------------------------------*\
	 * Now go through the branching stuff
	\*----------------------------------------------------------------------- */

	path = (bc_node **)malloc((2 * (level + 1) + BB_BUNCH) * sizeof(bc_node *));

	for (i = level, n = node; i >= 0; n = n->parent, i--)
		path[i] = n;

	/*------------------------------------------------------------------------*\
	 * Read in the basis.
	 * This is cplex style. sorry about it... Still, it
	 * might be ok if {VAR,SLACK}_{B,LB,UB} are properly defined
	\*----------------------------------------------------------------------- */

	node_desc *new_desc;

	int varexp_ind = 0, cutexp_ind = 0, nfexp_ind = 0;
	int bv_ind = 0, br_ind = 0, ev_ind = 0, er_ind = 0;
	array_desc extravar = {EXPLICIT_LIST, 0, 0, NULL};
	array_desc extrarow = {EXPLICIT_LIST, 0, 0, NULL};
	array_desc not_fixed = {EXPLICIT_LIST, 0, 0, NULL};
	basis_desc basis;

	int *list, *stat;

	char deal_with_nf = (desc->nf_status == NF_CHECK_AFTER_LAST ||
						 desc->nf_status == NF_CHECK_UNTIL_LAST);

	memset((char *)(&basis), 0, sizeof(basis_desc));

	/*------------------------------------------------------------------------*\
	 * First go up in the search tree to the root and record for every field
	 * the node where the first explicit description occurs.
	\*------------------------------------------------------------------------*/

	if (desc->uind.type == NO_DATA_STORED)
	{
		varexp_ind = -1;
	}
	else
	{
		for (i = level, n = node; !varexp_ind && i > 0; n = n->parent, i--)
			if (n->desc.uind.type == EXPLICIT_LIST)
				varexp_ind = i;
	}
	if (desc->cutind.type == NO_DATA_STORED)
	{
		cutexp_ind = -1;
	}
	else
	{
		for (i = level, n = node; !cutexp_ind && i > 0; n = n->parent, i--)
			if (n->desc.cutind.type == EXPLICIT_LIST)
				cutexp_ind = i;
	}
	if (deal_with_nf)
	{
		for (i = level, n = node; !nfexp_ind && i > 0; n = n->parent, i--)
			if (n->desc.not_fixed.type == EXPLICIT_LIST)
				nfexp_ind = i;
	}
	if ((basis.basis_exists = desc->basis.basis_exists) == TRUE)
	{
		for (i = level, n = node; !bv_ind && i > 0; n = n->parent, i--)
			if (n->desc.basis.basevars.type == EXPLICIT_LIST)
				bv_ind = i;
		for (i = level, n = node; !br_ind && i > 0; n = n->parent, i--)
			if (n->desc.basis.baserows.type == EXPLICIT_LIST)
				br_ind = i;
		for (i = level, n = node; !ev_ind && i > 0; n = n->parent, i--)
			if (n->desc.basis.extravars.type == EXPLICIT_LIST)
				ev_ind = i;
		for (i = level, n = node; !er_ind && i > 0; n = n->parent, i--)
			if (n->desc.basis.extrarows.type == EXPLICIT_LIST)
				er_ind = i;
	}
	else
	{
		ev_ind = er_ind = level;
	}

	/* An upper estimate on the total length of arrays */
	if (varexp_ind >= 0)
	{
		extravar.size = (n = path[varexp_ind])->desc.uind.size;
		for (i = varexp_ind + 1; i <= level; i++)
			extravar.size += path[i]->desc.uind.added;
	}
	if (cutexp_ind >= 0)
	{
		extrarow.size = (n = path[cutexp_ind])->desc.cutind.size;
		for (i = cutexp_ind + 1; i <= level; i++)
			extrarow.size += path[i]->desc.cutind.added;
	}
	if (deal_with_nf && nfexp_ind >= 0)
	{
		not_fixed.size = (n = path[nfexp_ind])->desc.not_fixed.size;
		for (i = nfexp_ind + 1; i <= level; i++)
			not_fixed.size += path[i]->desc.not_fixed.added;
	}
	else
	{
		not_fixed.size = 0;
	}

	/* If the LP function is compiled into the tree manager as a single
	   executable, then we allocate new memory for these arrays since
	   these arrays will be used directly instead of being passed
	   through PVM. */
	if (extravar.size)
	{
		extravar.list = (int *)malloc(extravar.size * ISIZE);
		if (basis.basis_exists)
			basis.extravars.stat = (int *)malloc(extravar.size * ISIZE);
	}
	if (extrarow.size)
	{
		extrarow.list = (int *)malloc(extrarow.size * ISIZE);
		if (basis.basis_exists)
			basis.extrarows.stat = (int *)malloc(extrarow.size * ISIZE);
	}
	if (not_fixed.size)
		not_fixed.list = (int *)malloc(not_fixed.size * ISIZE);
	if (env->base->varnum && basis.basis_exists)
		basis.basevars.stat = (int *)malloc(env->base->varnum * ISIZE);
	if (env->base->cutnum && basis.basis_exists)
		basis.baserows.stat = (int *)malloc(env->base->cutnum * ISIZE);

	/* The extra variables (uind) and the corresponding basis part */
	if (varexp_ind >= 0)
	{
		extravar.size = (n = path[varexp_ind])->desc.uind.size;
		if (extravar.size > 0)
			memcpy(extravar.list, n->desc.uind.list, ISIZE * extravar.size);
		for (i = varexp_ind + 1; i <= ev_ind; i++)
			modify_list(&extravar, &path[i]->desc.uind);
		if (basis.basis_exists)
		{
			/* at this point i == ev_ind */

			if (path[ev_ind]->desc.basis.extravars.size > 0)
				memcpy(basis.extravars.stat,
					   path[ev_ind]->desc.basis.extravars.stat,
					   path[ev_ind]->desc.basis.extravars.size * ISIZE);
			for (i = ev_ind + 1; i <= level; i++)
			{
				modify_list_and_stat(&extravar, basis.extravars.stat,
									 &path[i]->desc.uind,
									 &path[i]->desc.basis.extravars);
			}
			/* Although we send an explicit list, the type is sent over to show
			   to the LP process how extravars are stored in TM */
			basis.extravars.type = node->desc.basis.extravars.type;
			basis.extravars.size = extravar.size;
			basis.extravars.list = NULL;
			/* Now extravar.list/extravar.size and basis.extravars are OK */
			/* Fix basis.basevars */
			basis.basevars.type = EXPLICIT_LIST;
			basis.basevars.size = path[bv_ind]->desc.basis.basevars.size;
			basis.basevars.list = NULL;
			if (basis.basevars.size > 0)
			{
				memcpy(basis.basevars.stat,
					   path[bv_ind]->desc.basis.basevars.stat,
					   basis.basevars.size * ISIZE);
				for (i = bv_ind + 1; i <= level; i++)
				{
					list = path[i]->desc.basis.basevars.list;
					stat = path[i]->desc.basis.basevars.stat;
					for (j = path[i]->desc.basis.basevars.size - 1; j >= 0; j--)
						basis.basevars.stat[list[j]] = stat[j];
				}
			}
		}
	}

	/* Now take care of cutind and the corresponding basis part */
	if (cutexp_ind >= 0)
	{
		extrarow.size = (n = path[cutexp_ind])->desc.cutind.size;
		if (extrarow.size > 0)
			memcpy(extrarow.list, n->desc.cutind.list, ISIZE * extrarow.size);
		for (i = cutexp_ind + 1; i <= er_ind; i++)
			modify_list(&extrarow, &path[i]->desc.cutind);
		if (basis.basis_exists)
		{
			/* at this point i == er_ind */
			if (path[er_ind]->desc.basis.extrarows.size > 0)
				memcpy(basis.extrarows.stat,
					   path[er_ind]->desc.basis.extrarows.stat,
					   path[er_ind]->desc.basis.extrarows.size * ISIZE);
			for (i = er_ind + 1; i <= level; i++)
			{
				modify_list_and_stat(&extrarow, basis.extrarows.stat,
									 &path[i]->desc.cutind,
									 &path[i]->desc.basis.extrarows);
			}
			/* Same trick as above */
			basis.extrarows.type = node->desc.basis.extrarows.type;
			basis.extrarows.size = extrarow.size;
			basis.extrarows.list = NULL;
			/* Now extrarow.list/extrarow.size and basis.extrarows are OK */
			/* Fix basis.baserows */
			basis.baserows.type = EXPLICIT_LIST;
			basis.baserows.size = path[br_ind]->desc.basis.baserows.size;
			basis.baserows.list = NULL;
			if (basis.baserows.size > 0)
			{
				memcpy(basis.baserows.stat,
					   path[br_ind]->desc.basis.baserows.stat,
					   basis.baserows.size * ISIZE);
				for (i = br_ind + 1; i <= level; i++)
				{
					list = path[i]->desc.basis.baserows.list;
					stat = path[i]->desc.basis.baserows.stat;
					for (j = path[i]->desc.basis.baserows.size - 1; j >= 0; j--)
						basis.baserows.stat[list[j]] = stat[j];
				}
			}
		}
	}

	/* Finally the not fixed ones */
	if (deal_with_nf)
	{
		not_fixed.size = (n = path[nfexp_ind])->desc.not_fixed.size;
		if (not_fixed.size > 0)
			memcpy(not_fixed.list, n->desc.not_fixed.list, ISIZE * not_fixed.size);
		for (i = nfexp_ind + 1; i <= level; i++)
			modify_list(&not_fixed, &path[i]->desc.not_fixed);
	}

	new_desc = (node_desc *)calloc(1, sizeof(node_desc));

	new_desc->nf_status = desc->nf_status;
	new_desc->basis = basis;
	if (deal_with_nf)
		new_desc->not_fixed = not_fixed;
	new_desc->uind = extravar;
	new_desc->cutind = extrarow;

#if 0
   /* The cuts themselves */
   if (extrarow.size > 0){
      new_desc->cuts = (cut_data **)
	 malloc(extrarow.size*sizeof(cut_data *));
      for (i = 0; i < extrarow.size; i++){
	 new_desc->cuts[i] = tm->cuts[extrarow.list[i]];
      }
   }
#endif

	/* User defined description */
	new_desc->desc_size = desc->desc_size;
	if (new_desc->desc_size > 0)
		memcpy((char *)new_desc->desc, (char *)desc->desc, new_desc->desc_size);

	/*------------------------------------------------------------------------*\
	 * Load the lp problem (load_lp is an lp solver dependent routine).
	\*----------------------------------------------------------------------- */

	//   lp_data->mip->m += new_desc->cutind.size;

	lp_data->m = env->mip->m;
	lp_data->n = env->mip->n;

	open_lp_solver(lp_data);
	// This is a bit of a hack to avoid making a full copy of env->mip
	// Should be calling load_lp_prob(lp_data, 0, 0);
	// Perhaps make MIPdesc an argument of load_lp_prob??

	MIPdesc *mip = env->mip;
	lp_data->si->loadProblem(lp_data->n, lp_data->m,
							 mip->matbeg, mip->matind,
							 mip->matval, mip->lb,
							 mip->ub, mip->obj,
							 mip->sense, mip->rhs,
							 mip->rngval);

	/*------------------------------------------------------------------------*\
	 * Now go through the branching stuff
	\*----------------------------------------------------------------------- */

	bpath = (branch_desc *)malloc((2 * (level + 1) + BB_BUNCH) * sizeof(branch_desc));

	for (i = 0; i < level; i++, bpath++)
	{
		for (j = path[i]->bobj.child_num - 1; j >= 0; j--)
			if (path[i]->children[j] == path[i + 1])
				break;
		bobj = &path[i]->bobj;
		bpath->type = bobj->type;
		bpath->name = bobj->name;
		bpath->sense = bobj->sense[j];
		bpath->rhs = bobj->rhs[j];
		bpath->range = bobj->range[j];
		bpath->branch = bobj->branch[j];
	}

	bpath = bpath - level;
	if (level)
	{
		for (i = 0; i < level; i++, bpath++)
		{
			if (bpath->type == BRANCHING_VARIABLE)
			{
				j = bpath->name;
				switch (bpath->sense)
				{
				case 'E':
					change_lbub(lp_data, j, bpath->rhs, bpath->rhs);
					break;
				case 'L':
					change_ub(lp_data, j, bpath->rhs);
					break;
				case 'G':
					change_lb(lp_data, j, bpath->rhs);
					break;
				case 'R':
					change_lbub(lp_data, j, bpath->rhs, bpath->rhs + bpath->range);
					break;
				}
			}
			else
			{ /* BRANCHING_CUT */
				j = bpath->name;
				change_row(lp_data, j, bpath->sense, bpath->rhs, bpath->range);
			}
		}
	}
	bpath = bpath - level;
	/*------------------------------------------------------------------------*\
	 * Add cuts here
	 * FIXME! ASSUMING ALL THE CUTS ARE EXPLICIT ROW!
	\*----------------------------------------------------------------------- */
	desc = new_desc;

	if (desc->cutind.size > 0)
	{
		lp_data->tmp.c = (char *)calloc(lp_data->m, CSIZE);
		lp_data->tmp.d = (double *)calloc(lp_data->m, DSIZE);
		size = desc->cutind.size;
		sense = (char *)malloc(size * CSIZE);
		rhs = (double *)malloc(size * DSIZE);
		matbeg = (int *)calloc(size + 1, ISIZE);
		matbeg[0] = 0;

		for (i = 0; i < env->warm_start->cut_num; i++)
		{
			for (j = 0; j < desc->cutind.size; j++)
			{
				if (i == desc->cutind.list[j])
				{
					cut = env->warm_start->cuts[i];
					nzcnt = ((int *)(cut->coef))[0];
					total_nzcnt += nzcnt;
					sense[j] = cut->sense;
					rhs[j] = cut->rhs;
					matbeg[j + 1] = matbeg[j] + nzcnt;
				}
			}
		}

		matind = (int *)malloc(total_nzcnt * ISIZE);
		matval = (double *)malloc(total_nzcnt * DSIZE);

		for (i = 0; i < env->warm_start->cut_num; i++)
		{
			for (j = 0; j < desc->cutind.size; j++)
			{
				if (i == desc->cutind.list[j])
				{
					cut = env->warm_start->cuts[i];
					nzcnt = matbeg[j + 1] - matbeg[j];
					memcpy(matind + matbeg[j], (int *)(cut->coef + ISIZE),
						   ISIZE * nzcnt);
					memcpy(matval + matbeg[j],
						   (double *)(cut->coef + (1 + nzcnt) * ISIZE),
						   DSIZE * nzcnt);
				}
			}
		}
		add_rows(lp_data, size, total_nzcnt, rhs, sense, matbeg, matind, matval);
	}

	/*----------------------------------------------------------------------- */
	/* Load The Basis */
	/*----------------------------------------------------------------------- */

	if (desc->basis.basis_exists == TRUE)
	{
		int *rstat, *cstat;
		if (desc->basis.extravars.size == 0)
		{
			cstat = desc->basis.basevars.stat;
		}
		else if (desc->basis.basevars.size == 0)
		{
			cstat = desc->basis.extravars.stat;
		}
		else
		{							 /* neither is zero */
			cstat = lp_data->tmp.i1; /* n */
			memcpy(cstat,
				   desc->basis.basevars.stat, desc->basis.basevars.size * ISIZE);
			memcpy(cstat + desc->basis.basevars.size,
				   desc->basis.extravars.stat, desc->basis.extravars.size * ISIZE);
		}
		if (desc->basis.extrarows.size == 0)
		{
			rstat = desc->basis.baserows.stat;
		}
		else if (desc->basis.baserows.size == 0)
		{
			rstat = desc->basis.extrarows.stat;
		}
		else
		{							 /* neither is zero */
			rstat = lp_data->tmp.i2; /* m */
			memcpy(rstat,
				   desc->basis.baserows.stat, desc->basis.baserows.size * ISIZE);
			memcpy(rstat + desc->basis.baserows.size,
				   desc->basis.extrarows.stat, desc->basis.extrarows.size * ISIZE);
		}
		load_basis(lp_data, cstat, rstat);
	}

	size_lp_arrays(lp_data, FALSE, TRUE, lp_data->m, lp_data->n, nzcnt);

	// This is a hack. If lp_data->slacks is allocated, we crash
	FREE(lp_data->slacks);

	return_value = dual_simplex(lp_data, &iterd);

	if (return_value == LP_D_UNBOUNDED || return_value == LP_ABANDONED ||
		return_value == LP_D_INFEASIBLE)
	{
		node->feasibility_status = INFEASIBLE_PRUNED;
		node->node_status = NODE_STATUS__PRUNED;
#ifdef SENSITIVITY_ANALYSIS
		if (lp_data->has_ray)
		{
			if (!node->rays)
			{
				node->rays = (double *)malloc((lp_data->m + lp_data->n) * DSIZE);
			}
			memcpy(node->rays, lp_data->raysol, (lp_data->m + lp_data->n) * DSIZE);
		}
#endif
	}

	if (return_value == LP_OPTIMAL || return_value == LP_D_OBJLIM ||
		return_value == LP_D_ITLIM)
	{
		// lp_data->x = (double *)malloc(DSIZE*lp_data->n);
		get_x(lp_data);
		for (i = lp_data->n - 1; i >= 0; i--)
		{
			colsol = lp_data->x[i];
			if (env->mip->is_int[i])
			{
				if (colsol - floor(colsol) > lp_data->lpetol &&
					ceil(colsol) - colsol > lp_data->lpetol)
				{
					break;
				}
			}
		}
		if (i < 0)
		{
			node->node_status = NODE_STATUS__PRUNED;
			node->feasibility_status = FEASIBLE_PRUNED;

			FREE(node->sol_ind);
			FREE(node->sol);

			node->sol_ind = xind = (int *)malloc(ISIZE * lp_data->n);
			node->sol = xval = (double *)malloc(DSIZE * lp_data->n);

			for (i = 0; i < lp_data->n; i++)
			{
				if (lp_data->x[i] > lpetol || lp_data->x[i] < -lpetol)
				{
					xind[cnt] = i;
					xval[cnt++] = lp_data->x[i];
				}
			}
			node->sol_size = cnt;

			if ((best_sol->has_sol &&
				 lp_data->objval < best_sol->objval) ||
				!best_sol->has_sol)
			{

				best_sol->has_sol = TRUE;
				env->has_ub = TRUE;

				double gran_round = floor(env->par.lp_par.granularity + 0.5);
				if (env->par.lp_par.granularity > lp_data->lpetol &&
					fabs(gran_round - env->par.lp_par.granularity) < lp_data->lpetol)
				{
					/* we have granularity. symphony now uses granularity to set ub on
					 * lp-solver using granularity. so we round the solution to the
					 * nearest integer so that this tighter ub does not cut off other
					 * good solutions.
					 */
					best_sol->objval = floor(lp_data->objval + 0.5);
				}
				else
				{
					best_sol->objval = lp_data->objval;
				}
				env->ub = best_sol->objval;
				env->has_ub = TRUE;

				FREE(best_sol->xind);
				FREE(best_sol->xval);

				best_sol->xind = (int *)malloc(ISIZE * cnt);
				best_sol->xval = (double *)malloc(DSIZE * cnt);
				memcpy(best_sol->xind, xind, ISIZE * cnt);
				memcpy(best_sol->xval, xval, DSIZE * cnt);
				best_sol->xlength = cnt;
				best_sol->xlevel = node->bc_level;
				best_sol->xindex = node->bc_index;
				best_sol->lpetol = lpetol;
			}
			else
			{
				node->sol_size = 0;
				FREE(node->sol_ind);
				FREE(node->sol);
				node->node_status = NODE_STATUS__CANDIDATE;
			}
		}
#ifdef SENSITIVITY_ANALYSIS
		if (lp_data->dualsol)
		{
			if (!node->duals)
			{
				node->duals = (double *)malloc(lp_data->m * DSIZE);
			}
			if (!node->dj)
			{
				node->dj = (double *)malloc(lp_data->n * DSIZE);
			}
			memcpy(node->duals, lp_data->dualsol, lp_data->m * DSIZE);
			memcpy(node->dj, lp_data->dj, lp_data->n * DSIZE);
		}
#endif
	}

	for (i = 0; i < 2 * (level + 1) + BB_BUNCH; i++)
	{
		path[i] = 0;
	}

	FREE(path);
	FREE(bpath);
	FREE(sense);
	FREE(rhs);
	FREE(matind);
	FREE(matbeg);
	FREE(matval);
	/* FIXME- for now just copy free_node_desc here. Decide where to carry
	   resolve_node() */
	//   free_node_desc(&desc);
	if (desc)
	{
		node_desc *n = desc;
		FREE(n->cutind.list);
		FREE(n->uind.list);
		if (n->nf_status == NF_CHECK_AFTER_LAST ||
			n->nf_status == NF_CHECK_UNTIL_LAST)
			FREE(n->not_fixed.list);
		if (n->basis.basis_exists)
		{
			FREE(n->basis.basevars.list);
			FREE(n->basis.basevars.stat);
			FREE(n->basis.extravars.list);
			FREE(n->basis.extravars.stat);
			FREE(n->basis.baserows.list);
			FREE(n->basis.baserows.stat);
			FREE(n->basis.extrarows.list);
			FREE(n->basis.extrarows.stat);
		}
		if (n->desc_size > 0)
			FREE(n->desc);
		FREE(desc);
	}
	node->lower_bound = lp_data->objval;
	free_lp_arrays(lp_data);
	close_lp_solver(lp_data);
	FREE(lp_data);

	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int update_tree_bound(sym_environment *env, bc_node *root, int *cut_num,
					  int *cuts_ind, char *cru_vars, int change_type)
{

	int i, resolve = 0;
	char deletable = TRUE;

	if (root)
	{

		check_trim_tree(env, root, cut_num, cuts_ind, change_type);

		if (root->node_status == NODE_STATUS__PRUNED ||
			root->node_status == NODE_STATUS__TIME_LIMIT ||
			root->node_status == NODE_STATUS__ITERATION_LIMIT ||
			root->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
			root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
		{
			if (change_type == OBJ_COEFF_CHANGED || change_type == RHS_CHANGED ||
				change_type == COL_BOUNDS_CHANGED || change_type == COLS_ADDED)
			{
				if (root->feasibility_status == FEASIBLE_PRUNED ||
					root->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
					root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
				{

					check_better_solution(env, root, FALSE, change_type);
				}

				if (root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION &&
					root->bobj.child_num > 0)
				{
					for (i = 0; i < root->bobj.child_num; i++)
					{
						if (!update_tree_bound(env, root->children[i], cut_num, cuts_ind, cru_vars, change_type))
						{
							deletable = FALSE;
						}
					}
					if (change_type == COL_BOUNDS_CHANGED && root->bobj.child_num > 0)
					{
						update_branching_decisions(env, root, change_type);
					}
				}
				else
				{
					if (root->node_status == NODE_STATUS__WSPRUNED)
						root->node_status = NODE_STATUS__PRUNED;
					else
						root->node_status = NODE_STATUS__WARM_STARTED;
					if (resolve == 0)
					{
						root->lower_bound = -MAXDOUBLE;
					}
					else
					{
						resolve_node(env, root);
					}
				}
				if (resolve == 0)
				{
					root->feasibility_status = 0; // or? ROOT_NODE;
				}
			}
		}
		else
		{
			if (root->bobj.child_num > 0)
			{
				if (env->mip->var_type_modified)
				{
					if (!env->mip->is_int[root->children[0]->bobj.name])
					{
						for (i = 0; i < root->bobj.child_num; i++)
						{
							ws_free_subtree(env, root->children[i], change_type, TRUE, TRUE);
						}
						root->bobj.child_num = 0;
						root->node_status = NODE_STATUS__WARM_STARTED;
						if (resolve == 0)
						{
							root->lower_bound = -MAXDOUBLE;
						}
						else
						{
							resolve_node(env, root);
						}
					}
				}
				else
				{
					if (change_type == COL_BOUNDS_CHANGED && root->bobj.child_num > 0)
					{
						update_branching_decisions(env, root, change_type);
					}
					for (i = 0; i < root->bobj.child_num; i++)
					{
						if (!update_tree_bound(env, root->children[i], cut_num, cuts_ind, cru_vars, change_type))
						{
							deletable = FALSE;
						}
					}
				}
			}
			else
			{
				if (root->node_status == NODE_STATUS__WSPRUNED)
					root->node_status = NODE_STATUS__PRUNED;
				else
					root->node_status = NODE_STATUS__WARM_STARTED;
				if (resolve == 0)
				{
					root->lower_bound = -MAXDOUBLE;
				}
				else
				{
					resolve_node(env, root);
				}
			}
		}

		/* should be before resolve!!!*/
		if (change_type == COLS_ADDED)
		{
			update_node_desc(env, root, change_type);
		}

		if (env->warm_start->trim_tree == ON_CRU_VARS)
		{
			if (deletable && root->bobj.child_num)
			{
				for (i = 0; i < root->bobj.child_num; i++)
				{
					ws_free_subtree(env, root->children[i], change_type, FALSE, TRUE);
				}
				root->node_status = NODE_STATUS__WARM_STARTED;
				if (resolve == 0)
				{
					root->lower_bound = -MAXDOUBLE;
				}
				else
				{
					resolve_node(env, root);
				}
				root->bobj.child_num = 0;
				if (root->bc_level)
				{
					if (cru_vars[root->parent->bobj.name])
					{
						deletable = FALSE;
					}
				}
			}
			if (!deletable && root->bobj.child_num)
			{
				for (i = 0; i < root->bobj.child_num; i++)
				{
					register_cuts(root->children[i], cut_num, cuts_ind);
				}
			}
			if (root->bc_level)
			{
				if (cru_vars[root->parent->bobj.name])
				{
					deletable = FALSE;
				}
			}
		}
	}
	return deletable;
}

/*===========================================================================*/
/*===========================================================================*/

void register_cuts(bc_node *root, int *cut_num, int *cuts_ind)
{

	int i, r_cnum = root->desc.cutind.size;
	int *c_list = root->desc.cutind.list;
	int c_ind;
	if (r_cnum > 0)
	{
		for (i = 0; i < r_cnum; i++)
		{
			c_ind = c_list[i];
			if (cuts_ind[c_ind] < 0)
			{
				cuts_ind[c_ind] = c_list[i] = *cut_num;
				(*cut_num)++;
			}
			else
			{
				c_list[i] = cuts_ind[c_ind];
			}
		}
	}
}

/*===========================================================================*/
/*===========================================================================*/

void update_node_desc(sym_environment *env, bc_node *root, int change_type)
{

	node_desc *desc = &(root->desc);
	int i, tmp_size, colnum = env->mip->n;
	int added_colnum = env->mip->new_col_num;
	int old_colnum = colnum - added_colnum;

	if (change_type == COLS_ADDED && added_colnum > 0)
	{
		if (desc->uind.type == EXPLICIT_LIST)
		{
			tmp_size = desc->uind.size;
			if (tmp_size < colnum)
			{
				desc->uind.list = (int *)realloc((int *)desc->uind.list,
												 (tmp_size + added_colnum) * ISIZE);
				for (i = 0; i < added_colnum; i++)
				{
					desc->uind.list[tmp_size + i] = old_colnum + i;
				}
				desc->uind.size += added_colnum;
			}
		}
		if (desc->basis.basis_exists &&
			desc->basis.extravars.type == EXPLICIT_LIST)
		{
			tmp_size = desc->basis.extravars.size;
			if (tmp_size < colnum)
			{
				/*  if(desc->basis.extravars.list){ */
				/* 	       desc->basis.extravars.list =  */
				/* 		  (int*)realloc((int*)desc->basis.extravars.list,  */
				/* 				(tmp_size + added_colnum) *ISIZE); */
				/* 	       for(i = 0; i < added_colnum; i++){ */
				/* 		  desc->basis.extravars.list[tmp_size + i] = old_colnum + i; */
				/* 	       } */
				/* 	    } */
				desc->basis.extravars.stat =
					(int *)realloc((int *)desc->basis.extravars.stat,
								   (tmp_size + added_colnum) * ISIZE);
				for (i = 0; i < added_colnum; i++)
				{
					desc->basis.extravars.stat[tmp_size + i] = VAR_AT_LB;
				}
				desc->basis.extravars.size += added_colnum;
			}
		}
	}
}

/*===========================================================================*/
/*===========================================================================*/
void update_branching_decisions(sym_environment *env, bc_node *root, int change_type)
{
	int i, j;
	int type = 0;
	int resolve = 0;
	int deleted = 0;

	if (change_type == COL_BOUNDS_CHANGED && root->bobj.child_num > 0)
	{
		for (i = 0; i < root->bobj.child_num; i++)
		{
			if (root->bobj.type == BRANCHING_VARIABLE)
			{
				switch (root->bobj.sense[i])
				{
				case 'E':
					printf("error1-update_warm_start_tree\n");
					exit(0);
				case 'L':
					if (type == 0)
					{
						if (root->bobj.rhs[i] < env->mip->lb[root->bobj.name])
						{
							root->bobj.rhs[i] = ceil(env->mip->lb[root->bobj.name]);
						}
						else if (root->bobj.rhs[i] > env->mip->ub[root->bobj.name])
						{
							root->bobj.rhs[i] = ceil(env->mip->ub[root->bobj.name]) - 1;
						}
					}
					else
					{
						if (root->bobj.rhs[i] < env->mip->lb[root->bobj.name])
						{
							root->children[i]->node_status = NODE_STATUS__WSPRUNED;
							root->children[i]->feasibility_status = INFEASIBLE_PRUNED;
							for (j = 0; j < root->children[i]->bobj.child_num; j++)
							{
								ws_free_subtree(env, root->children[i]->children[j], change_type, FALSE, TRUE);
							}
							root->children[i]->bobj.child_num = 0;
						}
						else if (root->bobj.rhs[i] > env->mip->ub[root->bobj.name])
						{
							root->bobj.rhs[i] = floor(env->mip->ub[root->bobj.name]);
						}
					}
					break;
				case 'G':

					if (type == 0)
					{
						if (root->bobj.rhs[i] > env->mip->ub[root->bobj.name])
						{
							root->bobj.rhs[i] = floor(env->mip->ub[root->bobj.name]);
						}
						else if (root->bobj.rhs[i] < env->mip->lb[root->bobj.name])
						{
							root->bobj.rhs[i] = floor(env->mip->lb[root->bobj.name]) + 1;
						}
					}
					else
					{
						if (root->bobj.rhs[i] > env->mip->ub[root->bobj.name])
						{
							root->children[i]->node_status = NODE_STATUS__WSPRUNED;
							root->children[i]->feasibility_status = INFEASIBLE_PRUNED;
							for (j = 0; j < root->children[i]->bobj.child_num; j++)
							{
								ws_free_subtree(env, root->children[i]->children[j], change_type, FALSE, TRUE);
							}
							root->children[i]->bobj.child_num = 0;
						}
						else if (root->bobj.rhs[i] < env->mip->lb[root->bobj.name])
						{
							root->bobj.rhs[i] = ceil(env->mip->lb[root->bobj.name]);
						}
					}

					break;
				case 'R':
					printf("error2-update_warm_start_tree\n");
					exit(0);
				}
			}
			else
			{
				printf("error3-update_warm_start_tree\n");
				exit(0);
			}
		}
		root->bobj.child_num -= deleted;
		if (root->bobj.child_num <= 0)
		{
			if (root->node_status != NODE_STATUS__WSPRUNED)
				root->node_status = NODE_STATUS__WARM_STARTED;
			if (resolve == 0)
			{
				root->lower_bound = MAXDOUBLE;
			}
			else
			{
				resolve_node(env, root);
			}
		}
	}
}

/*===========================================================================*/
/*===========================================================================*/
void check_trim_tree(sym_environment *env, bc_node *root, int *cut_num,
					 int *cuts_ind, int change_type)
{
	int i;
	char trim_type = env->warm_start->trim_tree;
	int level = env->warm_start->trim_tree_level;
	int index = env->warm_start->trim_tree_index;
	problem_stat *stat = &(env->warm_start->stat);
	int trim_subtree = FALSE;

	if (trim_type)
	{
		register_cuts(root, cut_num, cuts_ind);
	}

	if (trim_type == TRIM_INDEX)
	{
		if (root->bobj.child_num > 0)
		{
			for (i = 0; i < root->bobj.child_num; i++)
			{
				if (root->children[i]->bc_index <= index)
				{
					break;
				}
			}
			if (i >= root->bobj.child_num)
			{
				trim_subtree = TRUE;
			}
		}
	}
	else if (trim_type == TRIM_LEVEL)
	{
		if (root->bc_level >= level)
		{
			trim_subtree = TRUE;
		}
	}

	if (trim_subtree && root->bobj.child_num)
	{
		for (i = 0; i < root->bobj.child_num; i++)
		{
			ws_free_subtree(env, root->children[i], change_type, TRUE, FALSE);
		}
		root->bobj.child_num = 0;
	}
	else
	{
		for (i = 0; i < root->bobj.child_num; i++)
		{
			root->children[i]->bc_index = stat->tree_size++;
			stat->created++;
		}
	}

	if (root->node_status == NODE_BRANCHED_ON && root->bobj.child_num)
	{
		stat->analyzed++;
	}
}
/*===========================================================================*/
/*===========================================================================*/

void cut_ws_tree_index(sym_environment *env, bc_node *root, int index,
					   problem_stat *stat, int change_type)
{

	int i;

	if (root)
	{
		if (root->node_status != NODE_STATUS__CANDIDATE)
		{
			stat->analyzed++;
		}
		if (root->bobj.child_num > 0)
		{
			for (i = 0; i < root->bobj.child_num; i++)
			{
				if (root->children[i]->bc_index <= index)
				{
					break;
				}
			}
			if (i < root->bobj.child_num)
			{

				for (i = 0; i < root->bobj.child_num; i++)
				{
					root->children[i]->bc_index = ++stat->tree_size;
					stat->created++;
				}
				for (i = root->bobj.child_num - 1; i >= 0; i--)
				{
					cut_ws_tree_index(env, root->children[i], index, stat,
									  change_type);
				}
			}
			else
			{
				for (i = root->bobj.child_num - 1; i >= 0; i--)
					ws_free_subtree(env, root->children[i], change_type, TRUE, FALSE);
				root->bobj.child_num = 0;

				if (root->node_status == NODE_STATUS__BRANCHED_ON)
				{
					root->node_status = NODE_STATUS__WARM_STARTED;
				}
			}
		}
	}
}

/*===========================================================================*/
/*===========================================================================*/

void cut_ws_tree_level(sym_environment *env, bc_node *root, int level,
					   problem_stat *stat, int change_type)
{

	int i;
	if (root)
	{
		if (root->node_status != NODE_STATUS__CANDIDATE)
		{
			stat->analyzed++;
		}

		if (root->bc_level < level)
		{
			for (i = 0; i < root->bobj.child_num; i++)
			{
				root->children[i]->bc_index = ++stat->tree_size;
				stat->created++;
			}
			for (i = root->bobj.child_num - 1; i >= 0; i--)
			{
				cut_ws_tree_level(env, root->children[i], level, stat, change_type);
			}
		}
		if (root->bc_level == level)
		{
			for (i = root->bobj.child_num - 1; i >= 0; i--)
				ws_free_subtree(env, root->children[i], change_type, TRUE, FALSE);
			root->bobj.child_num = 0;

			if (root->node_status == NODE_STATUS__BRANCHED_ON)
			{
				root->node_status = NODE_STATUS__WARM_STARTED;
			}
		}
	}
}

/*===========================================================================*/
/*===========================================================================*/

void ws_free_subtree(sym_environment *env, bc_node *root, int change_type, int check_solution,
					 int update_stats)
{

	int i;

	if (root == NULL)
		return;

	if (check_solution)
	{
		if (root->feasibility_status == FEASIBLE_PRUNED ||
			root->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
			root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
		{
			check_better_solution(env, root, TRUE, change_type);
		}
	}

	for (i = root->bobj.child_num - 1; i >= 0; i--)
	{
		ws_free_subtree(env, root->children[i], change_type, check_solution, update_stats);
	}

	if (update_stats)
	{
		env->warm_start->stat.analyzed--;
		env->warm_start->stat.created--;
		env->warm_start->stat.tree_size--;
	}

	free_tree_node(root);
}

/*===========================================================================*/
/*===========================================================================*/
void check_better_solution(sym_environment *env, bc_node *root, int delete_node,
						   int change_type)
{

	int i, j, new_solution = FALSE;
	double upper_bound = 0.0, lpetol = 9.9999999999999995e-07, valuesi = 0;
	double obj[2] = {0.0, 0.0}, gamma, tau, objval;
	int *matbeg, *matind;
	double *rowact = NULL, *matval, *colsol = NULL;
	int feasible = TRUE;

#ifdef USE_SYM_APPLICATION
	/* for now, just assume it is infeasible */

#else

	MIPdesc *mip = env->mip;
	lp_sol *best_sol = &(env->best_sol);

	if (env->mip->var_type_modified == TRUE)
	{
		for (i = root->sol_size - 1; i >= 0; i--)
		{
			if (mip->is_int[root->sol_ind[i]])
			{
				valuesi = root->sol[i];
				if (valuesi - floor(valuesi) > lpetol &&
					ceil(valuesi) - valuesi > lpetol)
				{
					break;
				}
			}
		}
		feasible = i < 0 ? TRUE : FALSE;
	}

	if (feasible)
	{
		if (change_type == OBJ_COEFF_CHANGED || change_type == COLS_ADDED)
		{
			if (!env->par.mc_warm_start)
			{
				for (i = 0; i < root->sol_size; i++)
				{
					upper_bound += mip->obj[root->sol_ind[i]] * root->sol[i];
				}
				if ((env->has_ub &&
					 upper_bound < env->ub) ||
					!env->has_ub)
				{

					if (!env->has_ub)
					{
						env->has_ub = TRUE;
						best_sol->has_sol = TRUE;
					}

					env->ub = upper_bound;
					best_sol->objval = upper_bound;
					new_solution = true;
				}
			}
			else
			{

				gamma = env->par.lp_par.mc_gamma;
				tau = env->par.lp_par.mc_tau;

				for (i = 0; i < root->sol_size; i++)
				{
					obj[0] += mip->obj1[root->sol_ind[i]] * root->sol[i];
					obj[1] += mip->obj2[root->sol_ind[i]] * root->sol[i];
				}

				if (gamma != -1.0)
				{
					objval = gamma * obj[0] + tau * obj[1];
				}
				else
				{
					objval = tau * obj[1];
				}
				if ((env->has_mc_ub && (objval < env->mc_ub - lpetol -
													 MAX(0, MIN(gamma, tau)))) ||
					!env->has_mc_ub)
				{

					if (!env->has_mc_ub)
					{
						env->has_mc_ub = TRUE;
						env->has_ub = TRUE;
						best_sol->has_sol = TRUE;
					}

					env->mc_ub = env->ub = best_sol->objval = objval;
					env->obj[0] = obj[0];
					env->obj[1] = obj[1];
					new_solution = TRUE;
				}
			}
		}
		else if (change_type == RHS_CHANGED)
		{

			colsol = (double *)calloc(env->mip->n, DSIZE);
			rowact = (double *)calloc(env->mip->m, DSIZE);
			matbeg = env->mip->matbeg;
			matval = env->mip->matval;
			matind = env->mip->matind;

			for (i = 0; i < root->sol_size; i++)
			{
				upper_bound += mip->obj[root->sol_ind[i]] * root->sol[i];
				colsol[root->sol_ind[i]] = root->sol[i];
			}

			for (i = 0; i < env->mip->n; i++)
			{
				for (j = matbeg[i]; j < matbeg[i + 1]; j++)
				{
					rowact[matind[j]] += matval[j] * colsol[i];
				}
			}

			for (i = 0; i < env->mip->m; i++)
			{
				switch (env->mip->sense[i])
				{
				case 'L':
					if (rowact[i] > env->mip->rhs[i] + lpetol)
						feasible = FALSE;
					break;
				case 'G':
					if (rowact[i] < env->mip->rhs[i] - lpetol)
						feasible = FALSE;
					break;
				case 'E':
					if (!((rowact[i] > env->mip->rhs[i] - lpetol) &&
						  (rowact[i] < env->mip->rhs[i] + lpetol)))
						feasible = FALSE;
					break;
				case 'R':
					if (rowact[i] > env->mip->rhs[i] + lpetol ||
						rowact[i] < env->mip->rhs[i] - env->mip->rngval[i] - lpetol)
						feasible = FALSE;
					break;
				case 'N':
				default:
					break;
				}

				if (!feasible)
					break;
			}

			if (feasible)
			{
				if ((env->has_ub &&
					 upper_bound < env->ub) ||
					!env->has_ub)
				{

					if (!env->has_ub)
					{
						env->has_ub = TRUE;
						best_sol->has_sol = TRUE;
					}

					env->ub = upper_bound;
					best_sol->objval = upper_bound;
					new_solution = true;
				}
			}
		}
		else if (change_type == COL_BOUNDS_CHANGED)
		{

			colsol = (double *)calloc(env->mip->n, DSIZE);

			for (i = 0; i < root->sol_size; i++)
			{
				upper_bound += mip->obj[root->sol_ind[i]] * root->sol[i];
				colsol[root->sol_ind[i]] = root->sol[i];
			}

			for (i = 0; i < env->mip->n; i++)
			{
				if ((colsol[i] < env->mip->lb[i] - lpetol) ||
					(colsol[i] > env->mip->ub[i] + lpetol))
					feasible = FALSE;
			}

			if (feasible)
			{
				if ((env->has_ub &&
					 upper_bound < env->ub) ||
					!env->has_ub)
				{

					if (!env->has_ub)
					{
						env->has_ub = TRUE;
						best_sol->has_sol = TRUE;
					}

					env->ub = upper_bound;
					best_sol->objval = upper_bound;
					new_solution = true;
				}
			}
		}
	}

	if (new_solution)
	{

		best_sol->xlevel = root->bc_level;
		best_sol->xindex = root->bc_index;
		best_sol->xlength = root->sol_size;
		best_sol->lpetol = lpetol;

		FREE(best_sol->xind);
		FREE(best_sol->xval);

		//      if(delete_node){
		best_sol->xind = root->sol_ind;
		best_sol->xval = root->sol;
		root->sol_ind = 0;
		root->sol = 0;
		//      } else {
		//	 best_sol->xind = (int *)malloc(ISIZE*root->sol_size);
		//	 best_sol->xval = (double *)malloc(DSIZE*root->sol_size);
		//	 memcpy(best_sol->xind, root->sol_ind, ISIZE*root->sol_size);
		//	 memcpy(best_sol->xval, root->sol, DSIZE* root->sol_size);
		//}
	}
	FREE(rowact);
	FREE(colsol);
	FREE(root->sol);
	FREE(root->sol_ind);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int copy_node(warm_start_desc *ws, bc_node *n_to, bc_node *n_from)
{

	if (!n_to || !n_from)
	{
		printf("copy_node(): Empty node_structure(s)!\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
		;
	}

	n_to->bc_index = n_from->bc_index;
	n_to->bc_level = n_from->bc_level;

	n_to->lp = n_from->lp;
	n_to->cg = n_from->cg;
	n_to->cp = n_from->cp;

	n_to->lower_bound = n_from->lower_bound;
	n_to->opt_estimate = n_from->opt_estimate;
	n_to->node_status = n_from->node_status;
	n_to->feasibility_status = n_from->feasibility_status;
	n_to->sol_size = n_from->sol_size;

	if (n_from->feasibility_status == FEASIBLE_PRUNED ||
		n_from->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
		n_from->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
	{
		if (n_from->sol)
		{
			n_to->sol = (double *)malloc(n_from->sol_size * DSIZE);
			n_to->sol_ind = (int *)malloc(n_from->sol_size * ISIZE);
			memcpy(n_to->sol, n_from->sol, n_from->sol_size * DSIZE);
			memcpy(n_to->sol_ind, n_from->sol_ind, n_from->sol_size * ISIZE);
		}
	}

#ifdef SENSITIVITY_ANALYSIS
	// FIXME: This is a bit fragile. Assumes variable list doesn't
	// change, etc.
	if (n_from->duals)
	{
		// Why desc.uind.size (number of variables)?
		// duals should be of size == m
		n_to->duals = (double *)malloc(ws->rootnode->desc.uind.size * DSIZE);
		memcpy(n_to->duals, n_from->duals,
			   ws->rootnode->desc.uind.size * DSIZE);
	}

	// feb223
	// FIXME: intcpt can be inconsistent when send_node_desc()
	// is called
	n_to->intcpt = n_from->intcpt;

	if (n_from->dj)
	{
		n_to->dj = (double *)malloc(ws->rootnode->desc.uind.size * DSIZE);
		memcpy(n_to->dj, n_from->dj,
			   ws->rootnode->desc.uind.size * DSIZE);
	}

	if (n_from->rays)
	{
		n_to->rays = (double *)malloc(n_from->desc.uind.size * DSIZE);
		memcpy(n_to->rays, n_from->rays,
			   n_from->desc.uind.size * DSIZE);
	}
#endif

#ifdef TRACE_PATH
	n_to->optimal_path = n_from->optimal_path;
#endif

	n_to->bobj = n_from->bobj;

#if defined(COMPILING_FOR_LP) || defined(SYM_COMPILE_IN_LP)

	// FIXME, Do we need this while writing to file

#if 0
   if (n_from->bobj.row){   
      n_to->bobj.row = (waiting_row*) malloc(sizeof(waiting_row));
      memcpy(n_to->bobj.row, n_from->bobj.row, sizeof(waiting_row));
      
      n_to->bobj.row->matind = (int*)malloc(sizeof(int)*n_to->bobj.row->nzcnt);
      n_to->bobj.row->matval = 
	 (double*)malloc(sizeof(double)*n_to->bobj.row->nzcnt);
      
      memcpy(n_to->bobj.row->matind, n_from->bobj.row->matind, 
	     ISIZE*n_to->bobj.row->nzcnt);
      memcpy(n_to->bobj.row->matval, n_from->bobj.row->matval, 
	     DSIZE*n_to->bobj.row->nzcnt);
      
      n_to->bobj.row->cut = (cut_data*)malloc( sizeof(cut_data));
      memcpy(n_to->bobj.row->cut, n_from->bobj.row->cut, sizeof(cut_data));
      
      n_to->bobj.row->cut->coef = 
	 (char*)malloc(sizeof(char)*n_to->bobj.row->cut->size);
      memcpy(n_to->bobj.row->cut->coef, n_from->bobj.row->cut->coef, 
	     CSIZE*n_to->bobj.row->cut->size);   
   }
#endif

#endif

#ifndef MAX_CHILDREN_NUM
	n_to->bobj.sense = (char *)malloc(n_to->bobj.child_num * CSIZE);
	n_to->bobj.rhs = (double *)malloc(n_to->bobj.child_num * DSIZE);
	n_to->bobj.range = (double *)malloc(n_to->bobj.child_num * DSIZE);
	n_to->bobj.branch = (int *)malloc(n_to->bobj.child_num * ISIZE);

	n_to->bobj.objval = (double *)malloc(n_to->bobj.child_num * DSIZE);
	n_to->bobj.termcode = (int *)malloc(n_to->bobj.child_num * ISIZE);
	n_to->bobj.iterd = (int *)malloc(n_to->bobj.child_num * ISIZE);
	n_to->bobj.feasible = (int *)malloc(n_to->bobj.child_num * ISIZE);

#endif

	memcpy(n_to->bobj.sense, n_from->bobj.sense,
		   n_to->bobj.child_num * CSIZE);
	memcpy(n_to->bobj.rhs, n_from->bobj.rhs,
		   n_to->bobj.child_num * DSIZE);
	memcpy(n_to->bobj.range, n_from->bobj.range,
		   n_to->bobj.child_num * DSIZE);
	memcpy(n_to->bobj.branch, n_from->bobj.branch,
		   n_to->bobj.child_num * ISIZE);

#if 0
   /* We shouldn't need to copy this. It's only used in the LP */
   memcpy(n_to->bobj.objval, n_from->bobj.objval, 
	  n_to->bobj.child_num*DSIZE); 
   memcpy(n_to->bobj.termcode, n_from->bobj.termcode, 
	  n_to->bobj.child_num*ISIZE); 
   memcpy(n_to->bobj.iterd, n_from->bobj.iterd, 
	  n_to->bobj.child_num*ISIZE); 
   memcpy(n_to->bobj.feasible, n_from->bobj.feasible, 
	  n_to->bobj.child_num*ISIZE);
#endif

	n_to->desc = n_from->desc;

	if (n_to->desc.uind.size)
	{
		n_to->desc.uind.list = (int *)malloc(n_to->desc.uind.size * ISIZE);
		memcpy(n_to->desc.uind.list, n_from->desc.uind.list,
			   n_to->desc.uind.size * ISIZE);
	}
	if (n_to->desc.basis.basis_exists)
	{
		if (n_to->desc.basis.basevars.size)
		{
			n_to->desc.basis.basevars.stat =
				(int *)malloc(n_to->desc.basis.basevars.size * ISIZE);
			memcpy(n_to->desc.basis.basevars.stat,
				   n_from->desc.basis.basevars.stat,
				   n_to->desc.basis.basevars.size * ISIZE);
			if (n_to->desc.basis.basevars.type == WRT_PARENT)
			{
				n_to->desc.basis.basevars.list =
					(int *)malloc(n_to->desc.basis.basevars.size * ISIZE);
				memcpy(n_to->desc.basis.basevars.list,
					   n_from->desc.basis.basevars.list,
					   n_to->desc.basis.basevars.size * ISIZE);
			}
		}

		if (n_to->desc.basis.extravars.size)
		{
			n_to->desc.basis.extravars.stat =
				(int *)malloc(n_to->desc.basis.extravars.size * ISIZE);
			memcpy(n_to->desc.basis.extravars.stat,
				   n_from->desc.basis.extravars.stat,
				   n_to->desc.basis.extravars.size * ISIZE);
			if (n_to->desc.basis.extravars.type == WRT_PARENT)
			{
				n_to->desc.basis.extravars.list =
					(int *)malloc(n_to->desc.basis.extravars.size * ISIZE);
				memcpy(n_to->desc.basis.extravars.list,
					   n_from->desc.basis.extravars.list,
					   n_to->desc.basis.extravars.size * ISIZE);
			}
		}

		if (n_to->desc.basis.baserows.size)
		{
			n_to->desc.basis.baserows.stat =
				(int *)malloc(n_to->desc.basis.baserows.size * ISIZE);
			memcpy(n_to->desc.basis.baserows.stat,
				   n_from->desc.basis.baserows.stat,
				   n_to->desc.basis.baserows.size * ISIZE);
			if (n_to->desc.basis.baserows.type == WRT_PARENT)
			{
				n_to->desc.basis.baserows.list =
					(int *)malloc(n_to->desc.basis.baserows.size * ISIZE);
				memcpy(n_to->desc.basis.baserows.list,
					   n_from->desc.basis.baserows.list,
					   n_to->desc.basis.baserows.size * ISIZE);
			}
		}

		if (n_to->desc.basis.extrarows.size)
		{
			n_to->desc.basis.extrarows.stat =
				(int *)malloc(n_to->desc.basis.extrarows.size * ISIZE);
			memcpy(n_to->desc.basis.extrarows.stat,
				   n_from->desc.basis.extrarows.stat,
				   n_to->desc.basis.extrarows.size * ISIZE);
			if (n_to->desc.basis.extrarows.type == WRT_PARENT)
			{
				n_to->desc.basis.extrarows.list =
					(int *)malloc(n_to->desc.basis.extrarows.size * ISIZE);
				memcpy(n_to->desc.basis.extrarows.list,
					   n_from->desc.basis.extrarows.list,
					   n_to->desc.basis.extrarows.size * ISIZE);
			}
		}
	}
	if (n_to->desc.not_fixed.size)
	{
		n_to->desc.not_fixed.list =
			(int *)malloc(n_to->desc.not_fixed.size * ISIZE);
		memcpy(n_to->desc.not_fixed.list, n_from->desc.not_fixed.list,
			   n_to->desc.not_fixed.size * ISIZE);
	}

	if (n_to->desc.cutind.size)
	{
		n_to->desc.cutind.list = (int *)malloc(n_to->desc.cutind.size * ISIZE);
		memcpy(n_to->desc.cutind.list, n_from->desc.cutind.list,
			   n_to->desc.cutind.size * ISIZE);
	}

	if (n_to->desc.desc_size)
	{
		n_to->desc.desc = (char *)malloc(n_to->desc.desc_size * CSIZE);
		memcpy(n_to->desc.desc, n_from->desc.desc,
			   n_to->desc.desc_size * CSIZE);
	}

	if (n_to->desc.bnd_change)
	{
		n_to->desc.bnd_change = (bounds_change_desc *)
			calloc(sizeof(bounds_change_desc), 1);
		if (n_from->desc.bnd_change->num_changes)
		{
			n_to->desc.bnd_change->index = (int *)malloc(ISIZE * n_from->desc.bnd_change->num_changes);
			n_to->desc.bnd_change->lbub = (char *)malloc(CSIZE * n_from->desc.bnd_change->num_changes);
			n_to->desc.bnd_change->value = (double *)malloc(DSIZE * n_from->desc.bnd_change->num_changes);
			memcpy(n_to->desc.bnd_change->index,
				   n_from->desc.bnd_change->index,
				   ISIZE * n_from->desc.bnd_change->num_changes);
			memcpy(n_to->desc.bnd_change->lbub,
				   n_from->desc.bnd_change->lbub,
				   CSIZE * n_from->desc.bnd_change->num_changes);
			memcpy(n_to->desc.bnd_change->value,
				   n_from->desc.bnd_change->value,
				   DSIZE * n_from->desc.bnd_change->num_changes);
		}
		n_to->desc.bnd_change->num_changes =
			n_from->desc.bnd_change->num_changes;
	}

	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int copy_tree(warm_start_desc *ws,
			  bc_node *root_to, bc_node *root_from)
{
	int i, childNum;

	if (!root_to || !root_from)
	{
		printf("copy_tree(): Empty root node(s)!\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	if (root_from)
	{
		copy_node(ws, root_to, root_from);
		childNum = root_to->bobj.child_num;
		if (childNum)
		{
			root_to->children = (bc_node **)calloc(sizeof(bc_node *), childNum);
			for (i = 0; i < childNum; i++)
			{
				root_to->children[i] = (bc_node *)calloc(1, sizeof(bc_node));
				root_to->children[i]->parent = root_to;
				copy_tree(ws, root_to->children[i], root_from->children[i]);
			}
		}
	}
	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int write_node(bc_node *node, FILE *f)
{
	int i;

	if (!node)
	{
		printf("write_node(): Empty node!\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	fprintf(f, "\n\n");

	fprintf(f, " NODE_INDEX      : %i\n", node->bc_index);
	fprintf(f, " NODE_LEVEL      : %i\n", node->bc_level);
	fprintf(f, " LOWER_BOUND     : %.4f\n", node->lower_bound);
	fprintf(f, " NODE_STATUS     : %i\n", (int)node->node_status);
	fprintf(f, " FEAS_STATUS     : %i\n", (int)node->feasibility_status);
	fprintf(f, " NODE_LP         : %i\n", node->lp);
	fprintf(f, " NODE_CG         : %i\n", node->cg);
	fprintf(f, " NODE_CP         : %i\n", node->cp);
	fprintf(f, " OPT_ESTIMATE    : %.4f\n", node->opt_estimate);

#ifdef TRACE_PATH
	fprintf(f, " OPTIMAL_PATH    : %c\n", node->optimal_path);
#endif
	if (node->parent)
	{
		fprintf(f, " PARENT_INDEX    : %i\n", node->parent->bc_index);
	}
	else
	{
		fprintf(f, " PARENT_INDEX    : -1\n");
	}
	fprintf(f, " CHILDREN        : %i %i %i\n", (int)node->bobj.type,
			node->bobj.name, node->bobj.child_num);
	for (i = 0; i < node->bobj.child_num; i++)
	{
		fprintf(f, "%i %c %f %f %i\n", node->children[i]->bc_index,
				node->bobj.sense[i], node->bobj.rhs[i],
				node->bobj.range[i], node->bobj.branch[i]);
	}

	fprintf(f, " NODE_DESCRIPTION                 : %i\n", node->desc.nf_status);
	fprintf(f, " USER_INDICES(Type,Size,Added)    : %i %i %i\n",
			(int)node->desc.uind.type, node->desc.uind.size,
			node->desc.uind.added);

	for (i = 0; i < node->desc.uind.size; i++)
	{
		fprintf(f, " %i", node->desc.uind.list[i]);
	}
	fprintf(f, "\n");

	fprintf(f, " NOT_FIXED(Type,Size,Added)   : %i %i %i\n",
			(int)node->desc.not_fixed.type, node->desc.not_fixed.size,
			node->desc.not_fixed.added);

	for (i = 0; i < node->desc.not_fixed.size; i++)
	{
		fprintf(f, " %i", node->desc.not_fixed.list[i]);
	}

	fprintf(f, "\n");

	fprintf(f, " CUT_INDICES(Type,Size,Added)   : %i %i %i\n",
			(int)node->desc.cutind.type, node->desc.cutind.size,
			node->desc.cutind.added);

	for (i = 0; i < node->desc.cutind.size; i++)
	{
		fprintf(f, " %i", node->desc.cutind.list[i]);
	}
	fprintf(f, "\n");

	fprintf(f, " BASIS          : %i\n", (int)node->desc.basis.basis_exists);
	fprintf(f, " BASE_VARIABLES : %i %i\n", (int)node->desc.basis.basevars.type,
			node->desc.basis.basevars.size);
	if (node->desc.basis.basevars.type == WRT_PARENT)
	{
		for (i = 0; i < node->desc.basis.basevars.size; i++)
		{
			fprintf(f, " %i %i", node->desc.basis.basevars.list[i],
					node->desc.basis.basevars.stat[i]);
		}
	}
	else
	{
		for (i = 0; i < node->desc.basis.basevars.size; i++)
		{
			fprintf(f, " %i", node->desc.basis.basevars.stat[i]);
		}
	}
	fprintf(f, "\n");

	fprintf(f, " EXTRA_VARIABLES : %i %i\n",
			(int)node->desc.basis.extravars.type,
			node->desc.basis.extravars.size);
	if (node->desc.basis.extravars.type == WRT_PARENT)
	{
		for (i = 0; i < node->desc.basis.extravars.size; i++)
		{
			fprintf(f, " %i %i", node->desc.basis.extravars.list[i],
					node->desc.basis.extravars.stat[i]);
		}
	}
	else
	{
		for (i = 0; i < node->desc.basis.extravars.size; i++)
		{
			fprintf(f, " %i", node->desc.basis.extravars.stat[i]);
		}
	}

	fprintf(f, "\n");
	fprintf(f, " BASE_ROWS      : %i %i\n",
			(int)node->desc.basis.baserows.type,
			node->desc.basis.baserows.size);
	if (node->desc.basis.baserows.type == WRT_PARENT)
	{
		for (i = 0; i < node->desc.basis.baserows.size; i++)
		{
			fprintf(f, " %i %i", node->desc.basis.baserows.list[i],
					node->desc.basis.baserows.stat[i]);
		}
	}
	else
	{
		for (i = 0; i < node->desc.basis.baserows.size; i++)
		{
			fprintf(f, " %i", node->desc.basis.baserows.stat[i]);
		}
	}
	fprintf(f, "\n");

	fprintf(f, " EXTRA_ROWS       : %i %i\n",
			(int)node->desc.basis.extrarows.type,
			node->desc.basis.extrarows.size);
	if (node->desc.basis.extrarows.type == WRT_PARENT)
	{
		for (i = 0; i < node->desc.basis.extrarows.size; i++)
		{
			fprintf(f, " %i %i", node->desc.basis.extrarows.list[i],
					node->desc.basis.extrarows.stat[i]);
		}
	}
	else
	{
		for (i = 0; i < node->desc.basis.extrarows.size; i++)
		{
			fprintf(f, " %i", node->desc.basis.extrarows.stat[i]);
		}
	}

	fprintf(f, "\n");
	fprintf(f, " USER_DESC_SIZE_&_ELEMENTS       : %i\n",
			node->desc.desc_size);
	for (i = 0; i < node->desc.desc_size; i++)
	{
		fprintf(f, " %i", (int)node->desc.desc[i]);
	}

	fprintf(f, "\n");

// #ifdef SENSITIVITY_ANALYSIS
// 	fprintf(f, "================================\n");
// 	fprintf(f, " DUALS    : ");
// 	if (node->duals){
// 		for (int i = 0; i < 1; i++){
// 			fprintf(f, "%.3f ", node->duals[i]);
// 		}
// 		fprintf(f, "\n");
// 	}
// 	if (node->dj){
// 		fprintf(f, " DJs    : ");
// 		for (int i = 0; i < 4; i++){
// 			fprintf(f, "%.3f ", node->dj[i]);
// 		}
// 		fprintf(f, "\n");
// 	}
// 	if (node->rays){
// 		fprintf(f, " RAY    : ");
// 		for (int i = 0; i < 4 + 1; i++){
// 			fprintf(f, "%.3f ", node->rays[i]);
// 		}
// 		fprintf(f, "\n");
// 	}
// 	fprintf(f, " BASE INDEX    : ");
// 	for (int i = 0; i < node->basis_len; i++){
// 		fprintf(f, "%d ", node->basis_idx[i]);
// 	}
// 	fprintf(f, "\n");
// 	fprintf(f, "================================\n");
// #endif

	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int write_tree(bc_node *root, FILE *f)
{
	int i;
	if (!root)
	{
		printf("write_tree(): Empty root node!\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	write_node(root, f);

	for (i = 0; i < root->bobj.child_num; i++)
	{
		write_tree(root->children[i], f);
	}

	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int read_node(bc_node *node, FILE *f)
{
	char str[80];
	int i = 0, num = 0, ch = 0;

	if (!node || !f)
	{
		printf("read_node(): Empty node or unable to read from file!\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	fscanf(f, "%s %s %i", str, str, &node->bc_index);
	fscanf(f, "%s %s %i", str, str, &node->bc_level);
	fscanf(f, "%s %s %lf", str, str, &node->lower_bound);
	fscanf(f, "%s %s %i", str, str, &ch);
	node->node_status = (char)ch;
	fscanf(f, "%s %s %i", str, str, &node->lp);
	fscanf(f, "%s %s %i", str, str, &node->cg);
	fscanf(f, "%s %s %i", str, str, &node->cp);
	fscanf(f, "%s %s %lf", str, str, &node->opt_estimate);

#ifdef TRACE_PATH
	fscanf(f, "%s %s %c", str, str, &node->optimal_path);
#endif
	fscanf(f, "%s %s %i", str, str, &num);
	fscanf(f, "%s %s %i %i %i", str, str, &ch, &node->bobj.name,
		   &node->bobj.child_num);
	node->bobj.type = (char)ch;
	if (node->bobj.child_num)
	{
#ifndef MAX_CHILDREN_NUM
		node->bobj.sense = (char *)malloc(node->bobj.child_num * CSIZE);
		node->bobj.rhs = (double *)malloc(node->bobj.child_num * DSIZE);
		node->bobj.range = (double *)malloc(node->bobj.child_num * DSIZE);
		node->bobj.branch = (int *)malloc(node->bobj.child_num * ISIZE);
#endif
		for (i = 0; i < node->bobj.child_num; i++)
		{
			fscanf(f, "%i %c %lf %lf %i", &num, &node->bobj.sense[i],
				   &node->bobj.rhs[i], &node->bobj.range[i],
				   &node->bobj.branch[i]);
		}
	}

	fscanf(f, "%s %s %i", str, str, &node->desc.nf_status);
	fscanf(f, "%s %s %i %i %i", str, str, &ch, &node->desc.uind.size,
		   &node->desc.uind.added);
	node->desc.uind.type = (char)ch;

	if (node->desc.uind.size)
	{
		node->desc.uind.list = (int *)malloc(node->desc.uind.size * ISIZE);
		for (i = 0; i < node->desc.uind.size; i++)
		{
			fscanf(f, "%i", &node->desc.uind.list[i]);
		}
	}

	fscanf(f, "%s %s %i %i %i", str, str, &ch, &node->desc.not_fixed.size,
		   &node->desc.not_fixed.added);
	node->desc.not_fixed.type = (char)ch;

	if (node->desc.not_fixed.size)
	{
		node->desc.not_fixed.list =
			(int *)malloc(node->desc.not_fixed.size * ISIZE);
		for (i = 0; i < node->desc.not_fixed.size; i++)
		{
			fscanf(f, "%i", &node->desc.not_fixed.list[i]);
		}
	}

	fscanf(f, "%s %s %i %i %i", str, str, &ch, &node->desc.cutind.size,
		   &node->desc.cutind.added);
	node->desc.cutind.type = (char)ch;
	if (node->desc.cutind.size)
	{
		node->desc.cutind.list = (int *)malloc(node->desc.cutind.size * ISIZE);
		for (i = 0; i < node->desc.cutind.size; i++)
		{
			fscanf(f, "%i", &node->desc.cutind.list[i]);
		}
	}

	fscanf(f, "%s %s %i", str, str, &ch);
	node->desc.basis.basis_exists = (char)ch;
	fscanf(f, "%s %s %i %i", str, str, &ch, &node->desc.basis.basevars.size);
	node->desc.basis.basevars.type = (char)ch;
	if (node->desc.basis.basevars.size)
	{
		node->desc.basis.basevars.stat =
			(int *)malloc(node->desc.basis.basevars.size * ISIZE);
		if (node->desc.basis.basevars.type == WRT_PARENT)
		{
			node->desc.basis.basevars.list =
				(int *)malloc(node->desc.basis.basevars.size * ISIZE);
			for (i = 0; i < node->desc.basis.basevars.size; i++)
			{
				fscanf(f, "%i %i", &node->desc.basis.basevars.list[i],
					   &node->desc.basis.basevars.stat[i]);
			}
		}
		else
		{
			for (i = 0; i < node->desc.basis.basevars.size; i++)
				fscanf(f, "%i", &node->desc.basis.basevars.stat[i]);
		}
	}

	fscanf(f, "%s %s %i %i", str, str, &ch, &node->desc.basis.extravars.size);
	node->desc.basis.extravars.type = (char)ch;
	if (node->desc.basis.extravars.size)
	{
		node->desc.basis.extravars.stat =
			(int *)malloc(node->desc.basis.extravars.size * ISIZE);
		if (node->desc.basis.extravars.type == WRT_PARENT)
		{
			node->desc.basis.extravars.list =
				(int *)malloc(node->desc.basis.extravars.size * ISIZE);
			for (i = 0; i < node->desc.basis.extravars.size; i++)
			{
				fscanf(f, "%i %i", &node->desc.basis.extravars.list[i],
					   &node->desc.basis.extravars.stat[i]);
			}
		}
		else
		{
			for (i = 0; i < node->desc.basis.extravars.size; i++)
				fscanf(f, "%i", &node->desc.basis.extravars.stat[i]);
		}
	}

	fscanf(f, "%s %s %i %i", str, str, &ch, &node->desc.basis.baserows.size);
	node->desc.basis.baserows.type = (char)ch;
	if (node->desc.basis.baserows.size)
	{
		node->desc.basis.baserows.stat =
			(int *)malloc(node->desc.basis.baserows.size * ISIZE);
		if (node->desc.basis.baserows.type == WRT_PARENT)
		{
			node->desc.basis.baserows.list =
				(int *)malloc(node->desc.basis.baserows.size * ISIZE);
			for (i = 0; i < node->desc.basis.baserows.size; i++)
				fscanf(f, "%i %i", &node->desc.basis.baserows.list[i],
					   &node->desc.basis.baserows.stat[i]);
		}
		else
		{
			for (i = 0; i < node->desc.basis.baserows.size; i++)
				fscanf(f, "%i", &node->desc.basis.baserows.stat[i]);
		}
	}

	fscanf(f, "%s %s %i %i", str, str, &ch, &node->desc.basis.extrarows.size);
	node->desc.basis.extrarows.type = (char)ch;
	if (node->desc.basis.extrarows.size)
	{
		node->desc.basis.extrarows.stat =
			(int *)malloc(node->desc.basis.extrarows.size * ISIZE);
		if (node->desc.basis.extrarows.type == WRT_PARENT)
		{
			node->desc.basis.extrarows.list =
				(int *)malloc(node->desc.basis.extrarows.size * ISIZE);
			for (i = 0; i < node->desc.basis.extrarows.size; i++)
				fscanf(f, "%i %i", &node->desc.basis.extrarows.list[i],
					   &node->desc.basis.extrarows.stat[i]);
		}
		else
		{
			for (i = 0; i < node->desc.basis.extrarows.size; i++)
				fscanf(f, "%i", &node->desc.basis.extrarows.stat[i]);
		}
	}

	fscanf(f, "%s %s %i", str, str, &node->desc.desc_size);
	if (node->desc.desc_size)
	{
		node->desc.desc =
			(char *)malloc(node->desc.desc_size * CSIZE);
		for (i = 0; i < node->desc.desc_size; i++)
		{
			fscanf(f, "%i", &ch);
			node->desc.desc[i] = (char)ch;
		}
	}

	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int read_tree(bc_node *root, FILE *f)
{
	if (!root || !f)
	{
		printf("read_tree(): Empty node or unable to write!\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	read_node(root, f);

	int i, childNum = root->bobj.child_num;

	if (childNum != 0)
	{
		root->children = (bc_node **)malloc(sizeof(bc_node *) * childNum);
		for (i = 0; i < childNum; i++)
		{
			root->children[i] = (bc_node *)calloc(1, sizeof(bc_node));
			root->children[i]->parent = root;
			read_tree(root->children[i], f);
		}
	}

	return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc *create_copy_warm_start(warm_start_desc *ws)
{
	int i, num = 0, allocated_cut_num = 0;
	warm_start_desc *ws_copy;

	if (!ws)
	{
		printf("create_copy_warm_start():");
		printf("The warm start description is empty!\n");
		return (NULL);
	}

	ws_copy = (warm_start_desc *)calloc(1, sizeof(warm_start_desc));
	memcpy(ws_copy, ws, sizeof(warm_start_desc));
	num = ws_copy->cut_num;
	allocated_cut_num = ws_copy->allocated_cut_num;
	ws_copy->cuts = (cut_data **)
		calloc(allocated_cut_num, sizeof(cut_data *));
	for (i = 0; i < num; i++)
	{
		ws_copy->cuts[i] = (cut_data *)calloc(1, sizeof(cut_data));
		memcpy(ws_copy->cuts[i], ws->cuts[i], sizeof(cut_data));
		ws_copy->cuts[i]->coef =
			(char *)calloc(ws_copy->cuts[i]->size, CSIZE);
		memcpy(ws_copy->cuts[i]->coef, ws->cuts[i]->coef,
			   CSIZE * ws_copy->cuts[i]->size);
	}
	ws_copy->rootnode = (bc_node *)calloc(1, sizeof(bc_node));
	copy_tree(ws, ws_copy->rootnode, ws->rootnode);

	if (ws->best_sol.xlength)
	{
		ws_copy->best_sol.xind = (int *)malloc(ISIZE * ws->best_sol.xlength);
		ws_copy->best_sol.xval = (double *)malloc(DSIZE * ws->best_sol.xlength);

		memcpy(ws_copy->best_sol.xind, ws->best_sol.xind,
			   ISIZE * ws->best_sol.xlength);
		memcpy(ws_copy->best_sol.xval, ws->best_sol.xval,
			   DSIZE * ws->best_sol.xlength);
	}

	return (ws_copy);
}

/*===========================================================================*/
/*===========================================================================*/

MIPdesc *create_copy_mip_desc(MIPdesc *mip)
{
	MIPdesc *mip_copy;
	int i;

	if (mip)
	{
		mip_copy = (MIPdesc *)calloc(1, sizeof(MIPdesc));
		memcpy(mip_copy, mip, sizeof(MIPdesc));

		if (mip->n)
		{
			mip_copy->obj = (double *)malloc(DSIZE * mip_copy->n);

			mip_copy->ub = (double *)malloc(DSIZE * mip_copy->n);
			mip_copy->lb = (double *)malloc(DSIZE * mip_copy->n);
			mip_copy->is_int = (char *)malloc(CSIZE * mip_copy->n);
			mip_copy->matbeg = (int *)malloc(ISIZE * (mip_copy->n + 1));

			memcpy(mip_copy->obj, mip->obj, DSIZE * mip_copy->n);
			memcpy(mip_copy->ub, mip->ub, DSIZE * mip_copy->n);
			memcpy(mip_copy->lb, mip->lb, DSIZE * mip_copy->n);
			memcpy(mip_copy->is_int, mip->is_int, CSIZE * mip_copy->n);
			memcpy(mip_copy->matbeg, mip->matbeg, ISIZE * (mip_copy->n + 1));

			if (mip->obj1)
			{
				mip_copy->obj1 = (double *)malloc(DSIZE * mip_copy->n);
				memcpy(mip_copy->obj1, mip->obj1, DSIZE * mip_copy->n);
			}
			if (mip->obj2)
			{
				mip_copy->obj2 = (double *)malloc(DSIZE * mip_copy->n);
				memcpy(mip_copy->obj2, mip->obj2, DSIZE * mip_copy->n);
			}
		}

		if (mip->m)
		{
			mip_copy->rhs = (double *)malloc(DSIZE * mip_copy->m);
			mip_copy->sense = (char *)malloc(CSIZE * mip_copy->m);
			mip_copy->rngval = (double *)malloc(DSIZE * mip_copy->m);

			memcpy(mip_copy->rhs, mip->rhs, DSIZE * mip_copy->m);
			memcpy(mip_copy->sense, mip->sense, CSIZE * mip_copy->m);
			memcpy(mip_copy->rngval, mip->rngval, DSIZE * mip_copy->m);
		}

		if (mip->nz)
		{
			mip_copy->matval = (double *)malloc(DSIZE * mip_copy->nz);
			mip_copy->matind = (int *)malloc(ISIZE * mip_copy->nz);
			memcpy(mip_copy->matval, mip->matval, DSIZE * mip_copy->nz);
			memcpy(mip_copy->matind, mip->matind, ISIZE * mip_copy->nz);
		}

		/* will not be used anywhere other than presolve where they are
	   always initialized */
		mip_copy->mip_inf = 0;
		mip->cru_vars = 0;
		mip->orig_sense = 0;
		mip->orig_ind = 0;
		if (mip->row_matbeg)
		{
			mip_copy->row_matbeg = (int *)malloc(ISIZE * (mip_copy->m + 1));
			mip_copy->row_matind = (int *)malloc(ISIZE * mip_copy->nz);
			mip_copy->row_matval = (double *)malloc(DSIZE * mip_copy->nz);
			mip_copy->row_lengths = (int *)malloc(ISIZE * mip_copy->m);
			mip_copy->col_lengths = (int *)malloc(ISIZE * mip_copy->n);

			memcpy(mip_copy->row_matbeg, mip->row_matbeg,
				   ISIZE * (mip_copy->m + 1));
			memcpy(mip_copy->row_matind, mip->row_matind, ISIZE * mip_copy->nz);
			memcpy(mip_copy->row_matval, mip->row_matval, DSIZE * mip_copy->nz);
			memcpy(mip_copy->row_lengths, mip->row_lengths, ISIZE * mip_copy->m);
			memcpy(mip_copy->col_lengths, mip->col_lengths, ISIZE * mip_copy->n);
		}

		if (mip->colname)
		{
			mip_copy->colname = (char **)calloc(sizeof(char *), mip_copy->n);

			for (i = 0; i < mip_copy->n; i++)
			{
				/* FIXME! Resctricting col_name to 30 chars! */
				if (mip->colname[i])
				{
					mip_copy->colname[i] = (char *)malloc(CSIZE * MAX_NAME_SIZE);
					strncpy(mip_copy->colname[i], mip->colname[i], MAX_NAME_SIZE);
					mip_copy->colname[i][MAX_NAME_SIZE - 1] = 0;
				}
			}
		}

		if (mip->fixed_n)
		{
			memcpy(mip_copy->fixed_ind, mip->fixed_ind, ISIZE * mip->fixed_n);
			memcpy(mip_copy->fixed_val, mip->fixed_val, DSIZE * mip->fixed_n);
		}
	}
	else
	{
		printf("create_copy_mip_desc():");
		printf("Trying to copy an empty mip desc!\n");
		return (NULL);
	}

	return (mip_copy);
}

/*===========================================================================*/
/*===========================================================================*/

sym_environment *create_copy_environment(sym_environment *env)
{
	int i, j, num;
	sym_environment *env_copy;
	params *par;
	lp_sol *sol;
#if defined(SYM_COMPILE_IN_TM) && defined(SYM_COMPILE_IN_CP)
	cp_cut_data *cp_cut;
	cut_data *cut;
#endif
	sp_desc *sp;

	if (!env)
	{
		printf("create_copy_sym_environment(): The given problem is empty!\n");
		printf("Unable to copy.\n");
		return (NULL);
	}
	env_copy = (sym_environment *)calloc(1, sizeof(sym_environment));
	memcpy(env_copy, env, sizeof(sym_environment));

	//   initialize_u(env_copy);

	/* Note that, if some modifications on the user function have been done
	   after initialization, it will not be reflected here, since SYMPHONY
	   doesn't know anoything about the user structure! For a temporary
	   solution, the user pointer will be directed to the original user
	   structure! So, be careful from now on that, further modifications
	   on the user structure of either the original or the clone env will
	   affect the both! */

	env_copy->user = env->user;

	/*========================================================================*/
	/*   copy params */

	par = &(env_copy->par);

	if (par->tm_par.lp_mach_num)
		par->tm_par.lp_machs =
			(char **)malloc(sizeof(char *) * par->tm_par.lp_mach_num);
	if (par->tm_par.cg_mach_num)
		par->tm_par.cg_machs =
			(char **)malloc(sizeof(char *) * par->tm_par.cg_mach_num);
	if (par->tm_par.cp_mach_num)
		par->tm_par.cp_machs =
			(char **)malloc(sizeof(char *) * par->tm_par.cp_mach_num);
	for (i = 0; i < par->tm_par.lp_mach_num; i++)
	{
		par->tm_par.lp_machs[i] =
			(char *)malloc(CSIZE * (MACH_NAME_LENGTH + 1));
		memcpy(par->tm_par.lp_machs[i], env->par.tm_par.lp_machs[i],
			   CSIZE * (MACH_NAME_LENGTH + 1));
	}

	for (i = 0; i < par->tm_par.cg_mach_num; i++)
	{
		par->tm_par.cg_machs[i] =
			(char *)malloc(CSIZE * (MACH_NAME_LENGTH + 1));
		memcpy(par->tm_par.cg_machs[i], env->par.tm_par.cg_machs[i],
			   CSIZE * (MACH_NAME_LENGTH + 1));
	}

	for (i = 0; i < par->tm_par.cp_mach_num; i++)
	{
		par->tm_par.cp_machs[i] =
			(char *)malloc(CSIZE * (MACH_NAME_LENGTH + 1));
		memcpy(par->tm_par.cp_machs[i], env->par.tm_par.cp_machs[i],
			   CSIZE * (MACH_NAME_LENGTH + 1));
	}

	/*========================================================================*/

	/* copy lp_sol */

	sol = &(env_copy->best_sol);
	if (sol->xlength)
	{
		sol->xind = (int *)malloc(ISIZE * sol->xlength);
		sol->xval = (double *)malloc(DSIZE * sol->xlength);
		memcpy(sol->xind, env->best_sol.xind, ISIZE * sol->xlength);
		memcpy(sol->xval, env->best_sol.xval, DSIZE * sol->xlength);
	}

	/*========================================================================*/

	/* Copy solution pool */

	if (env->sp)
	{
		sp = env_copy->sp = (sp_desc *)malloc(sizeof(sp_desc));
		memcpy(sp, env->sp, sizeof(sp_desc));
		sp->solutions =
			(sp_solution **)calloc(sp->max_solutions, sizeof(sp_solution *));
		for (i = 0; i < sp->num_solutions; i++)
		{
			sp->solutions[i] = (sp_solution *)malloc(sizeof(sp_solution));
			sp->solutions[i]->xlength = env->sp->solutions[i]->xlength;
			if (sp->solutions[i]->xlength)
			{
				sp->solutions[i]->xind =
					(int *)malloc(ISIZE * sp->solutions[i]->xlength);
				sp->solutions[i]->xval =
					(double *)malloc(DSIZE * sp->solutions[i]->xlength);
				memcpy(sp->solutions[i]->xind, env->sp->solutions[i]->xind,
					   ISIZE * sp->solutions[i]->xlength);
				memcpy(sp->solutions[i]->xval, env->sp->solutions[i]->xval,
					   DSIZE * sp->solutions[i]->xlength);
			}
		}
	}

	/*========================================================================*/

	/* copy mip */
	if (env->mip)
	{
		// free_mip_desc(env_copy->mip);
		if (env->prep_mip)
		{
			env_copy->prep_mip = create_copy_mip_desc(env->prep_mip);
			if (env->orig_mip)
			{
				env_copy->orig_mip = create_copy_mip_desc(env->orig_mip);
				env_copy->mip = env_copy->orig_mip;
			}
			else
			{
				env_copy->mip = create_copy_mip_desc(env->mip);
			}
		}
		else
		{
			env_copy->mip = create_copy_mip_desc(env->mip);
			env_copy->prep_mip = env_copy->orig_mip = 0;
		}
	}

	/*========================================================================*/

	/* copy base_desc */

	if (env->base)
	{
		env_copy->base = (base_desc *)calloc(1, sizeof(base_desc));
		memcpy(env_copy->base, env->base, sizeof(base_desc));

		if (env->base->varnum)
		{
			env_copy->base->userind = (int *)malloc(ISIZE * env->base->varnum);
			memcpy(env_copy->base->userind, env->base->userind,
				   ISIZE * env->base->varnum);
		}
	}

	/*========================================================================*/

	/* copy root_desc */

	if (env->rootdesc)
	{
		node_desc *desc = env_copy->rootdesc =
			(node_desc *)calloc(1, sizeof(node_desc));
		memcpy(desc, env->rootdesc, sizeof(node_desc));

		if (desc->uind.size)
		{
			desc->uind.list = (int *)malloc(desc->uind.size * ISIZE);
			memcpy(desc->uind.list, env->rootdesc->uind.list,
				   desc->uind.size * ISIZE);
		}

		if (desc->not_fixed.size)
		{
			desc->not_fixed.list =
				(int *)malloc(desc->not_fixed.size * ISIZE);
			memcpy(desc->not_fixed.list, env->rootdesc->not_fixed.list,
				   desc->not_fixed.size * ISIZE);
		}

		if (desc->cutind.size)
		{
			desc->cutind.list = (int *)malloc(desc->cutind.size * ISIZE);
			memcpy(desc->cutind.list, env->rootdesc->cutind.list,
				   desc->cutind.size * ISIZE);
		}

		if (desc->desc_size)
		{
			desc->desc = (char *)malloc(desc->desc_size * CSIZE);
			memcpy(desc->desc, env->rootdesc->desc,
				   desc->desc_size * CSIZE);
		}

		if (desc->bnd_change)
		{
			desc->bnd_change = (bounds_change_desc *)
				calloc(sizeof(bounds_change_desc), 1);
			if (env->rootdesc->bnd_change->num_changes)
			{
				desc->bnd_change->index = (int *)malloc(ISIZE * env->rootdesc->bnd_change->num_changes);
				desc->bnd_change->lbub = (char *)malloc(CSIZE * env->rootdesc->bnd_change->num_changes);
				desc->bnd_change->value = (double *)malloc(DSIZE * env->rootdesc->bnd_change->num_changes);
				memcpy(desc->bnd_change->index,
					   env->rootdesc->bnd_change->index,
					   ISIZE * env->rootdesc->bnd_change->num_changes);
				memcpy(desc->bnd_change->lbub,
					   env->rootdesc->bnd_change->lbub,
					   CSIZE * env->rootdesc->bnd_change->num_changes);
				memcpy(desc->bnd_change->value,
					   env->rootdesc->bnd_change->value,
					   DSIZE * env->rootdesc->bnd_change->num_changes);
			}
			desc->bnd_change->num_changes =
				env->rootdesc->bnd_change->num_changes;
		}
	}

	/*========================================================================*/
	/* jump the tm */

	/*========================================================================*/
	/* copy the warm start */

	if (env->warm_start)
	{
		env_copy->warm_start = create_copy_warm_start(env->warm_start);
	}

	/*========================================================================*/
	/*copy the cut pool */

#if defined(SYM_COMPILE_IN_TM) && defined(SYM_COMPILE_IN_CP)
	if (env_copy->par.tm_par.max_cp_num > 1)
	{
		env_copy->cp =
			(cut_pool **)malloc(env_copy->par.tm_par.max_cp_num *
								sizeof(cut_pool *));
		for (i = 0; i < env_copy->par.tm_par.max_cp_num; i++)
		{
			env_copy->cp[i] = (cut_pool *)calloc(1, sizeof(cut_pool));
			env_copy->cp[i]->par = env_copy->par.cp_par;
#ifdef USE_SYM_APPLICATION
			user_send_cp_data(env_copy->user, &env_copy->cp[i]->user);
#else
			env_copy->cp[i]->user = env_copy->user;
#endif
		}
		num = env_copy->par.tm_par.max_cp_num;
	}
	else
	{
		num = 0;
	}

	if (num)
	{
		for (i = 0; i < num; i++)
		{
			memcpy(env_copy->cp[i], env->cp[i], sizeof(cut_pool));
			env_copy->cp[i]->cuts =
				(cp_cut_data **)malloc(env_copy->cp[i]->allocated_cut_num *
									   sizeof(cp_cut_data *));
			for (j = 0; j < env_copy->cp[i]->cut_num; j++)
			{
				env_copy->cp[i]->cuts[j] =
					(cp_cut_data *)calloc(1, sizeof(cp_cut_data));
				cp_cut = env_copy->cp[i]->cuts[j];
				memcpy(cp_cut, env->cp[i]->cuts[j], sizeof(cp_cut_data));
				cp_cut->cut.coef = (char *)malloc(cp_cut->cut.size * CSIZE);
				memcpy(cp_cut->cut.coef, env->cp[i]->cuts[j]->cut.coef,
					   cp_cut->cut.size * CSIZE);
			}

			sol = &(env_copy->cp[i]->cur_sol);

			sol->xind = (int *)malloc(ISIZE * sol->max_sol_length);
			sol->xval = (double *)malloc(DSIZE * sol->max_sol_length);
			memcpy(sol->xind, env->cp[i]->cur_sol.xind, ISIZE * sol->max_sol_length);
			memcpy(sol->xval, env->cp[i]->cur_sol.xval, DSIZE * sol->max_sol_length);

			num = env_copy->cp[i]->cuts_to_add_num;
			if (num)
			{
				env_copy->cp[i]->cuts_to_add =
					(cut_data **)malloc(num * sizeof(cut_data *));
				for (j = 0; j < num; j++)
				{
					env_copy->cp[i]->cuts_to_add[j] =
						(cut_data *)calloc(1, sizeof(cut_data));
					cut = env_copy->cp[i]->cuts_to_add[j];
					memcpy(cut, env_copy->cp[i]->cuts_to_add[j], sizeof(cut_data));
					cut->coef = (char *)malloc(cut->size * CSIZE);
					memcpy(cut->coef, env->cp[i]->cuts_to_add[j]->coef,
						   cut->size * CSIZE);
				}
			}
		}
	}
#endif

	return (env_copy);
}
/*===========================================================================*/
/*===========================================================================*/
// Anahita
void get_dual_pruned PROTO((bc_node * root, MIPdesc *mip,
							double **dual_pieces, int *cur_piece_no,
							int MAX_ALLOWABLE_NUM_PIECES))
{
#ifdef SENSITIVITY_ANALYSIS
	int i, j;
	bc_node *child;

	if (root)
	{
		// this is an inermediate node, do nothing immediately,
		// dive to children

		for (i = 0; i < root->bobj.child_num; i++)
		{

			child = root->children[i];

			// temp
			bool leaf_node = 1;
			if (child->children)
			{
				leaf_node = 0;
			}
			printf("level %d, child->node_status is %d and child->feasibility_status is %d.\n", child->bc_index, child->node_status, child->feasibility_status);

			if (child->node_status == NODE_STATUS__ITERATION_LIMIT)
			{
				printf("This child is pruned because of iteration limit");
			}
			// end temp
			if (child->node_status == NODE_STATUS__PRUNED || leaf_node)
			{
				if (*cur_piece_no > MAX_ALLOWABLE_NUM_PIECES)
				{
					printf("MAX_ALLOWABLE_NUM_PIECES reached...terminating..\n");
					exit(1);
				}
				if (child->feasibility_status == FEASIBLE_PRUNED ||
					child->feasibility_status == OVER_UB_PRUNED ||
					child->feasibility_status == PRUNED_HAS_CAN_SOLUTION || leaf_node)
				{

					// if pruned node is int feasible or it is
					// over UB pruned, get the dual info

					// allocate memory for dual_pieces
					dual_pieces[*cur_piece_no] = (double *)malloc((1 + mip->m) * sizeof(double));

					// write dual info ***CHECK***
					dual_pieces[*cur_piece_no][0] = child->intcpt; // buildCpx
					// dual_pieces[*cur_piece_no] [0] = child->lower_bound;//buildCpx-3
					for (j = 0; j < mip->m; j++)
					{
						dual_pieces[*cur_piece_no][j + 1] = child->duals[j];
					}
					// increment
					(*cur_piece_no)++;
				}

				else if (child->feasibility_status == INFEASIBLE_PRUNED)
				{
					printf("Infeasible node...Examining rays/duals\n");
					// sleep(2);//check if the node was pruned by dual bound exceed

					if (child->duals && child->rays)
					{
						printf("HUH! Both duals and rays\n");
					}
					if (child->duals && !child->rays)
					{
						printf("HUH! Just duals\n");
						exit(1);
					}
					if (!child->duals && child->rays)
					{
						printf("HUH! Just rays\n");
						exit(1);
					}

					// Combine the solutions

					double Lambda = 1;

					dual_pieces[*cur_piece_no] = (double *)malloc((1 + mip->m) * sizeof(double));

					dual_pieces[*cur_piece_no][0] = Lambda * child->intcpt;
					for (j = 0; j < mip->m; j++)
					{
						dual_pieces[*cur_piece_no][j + 1] = 99999 * Lambda * child->rays[j];
					}
					// increment
					(*cur_piece_no)++;

					/* //write dual+ray info */
					/* /\* here child-duals has the dual info */
					/* 	  and child-lowerbound is the dual obj-->first we get the */
					/* 	  intercept of the duals. Later we will work on rays *\/ */
					/* child->lower_bound; */

					/* /\* */
					/* 	  child-incpt has the reduced cost info that is the variable bounds*their */
					/* 	  duals, with rays *\/ */
					/* dual_pieces[*cur_piece_no] [0] = child->lower_bound; */
					/* for (j = 0; j < mip->m; j++){ */
					/* 	  dual_pieces[*cur_piece_no][j+1] = child->duals[j] + Lambda* */
					/* 	     child->rays[i]; */
					/* } */

				} // end of infeasible nodes
				else
				{ // status unknown
					printf("get_dual_pruned(): Unknown error!\n");
					exit(1);
				}
			} // if not pruned, recurse
			else
			{
				get_dual_pruned(child, mip, dual_pieces, cur_piece_no,
								MAX_ALLOWABLE_NUM_PIECES);
			}
		} // child loop
	}	  // if root end

#else
	printf("get_dual_pruned():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
#endif
}

/*===========================================================================*/
/*===========================================================================*/
double get_lb_for_new_rhs(bc_node *node, MIPdesc *mip, branch_desc *bpath,
						  int rhs_cnt, int *new_rhs_ind, double *new_rhs_val,
						  int lb_cnt, int *new_lb_ind, double *new_lb_val,
						  int ub_cnt, int *new_ub_ind, double *new_ub_val)
{
#ifdef SENSITIVITY_ANALYSIS
	int i, j, k, retval;
	double min = -SYM_INFINITY;
	bc_node *child;
	double lb = 0, objval = -SYM_INFINITY;
	int level = node->bc_level;
	branch_obj *bobj;

	if (!node)
	{
		printf("Warning: NULL pointer in get_lb_for_new_rhs()\n");
		return (-SYM_INFINITY);
	}

#if 0
   bc_node **path, *n;
   branch_desc *bpath ;

   path = (bc_node **) malloc((2*(level+1)+BB_BUNCH)*sizeof(bc_node *));
   bpath = (branch_desc *) malloc 
      ((2*(level+1)+BB_BUNCH)*sizeof(branch_desc));

   for (i = level, n = node; i >= 0; n = n->parent, i--)
      path[i] = n;

   for (i = 0; i < level; i++){
      for (j = path[i]->bobj.child_num - 1; j >= 0; j--)
	 if (path[i]->children[j] == path[i+1])
	    break;
      bobj = &path[i]->bobj;
      bpath[i].type = bobj->type;
      bpath[i].name = bobj->name;
      bpath[i].sense = bobj->sense[j];
      bpath[i].rhs = bobj->rhs[j];
      bpath[i].range = bobj->range[j];
      bpath[i].branch = bobj->branch[j];
   }
#else
	if (level > 0)
	{
		for (j = node->parent->bobj.child_num - 1; j >= 0; j--)
			if (node->parent->children[j] == node)
				break;

		bobj = &(node->parent->bobj);
		bpath[level - 1].type = bobj->type;
		bpath[level - 1].name = bobj->name;
		bpath[level - 1].sense = bobj->sense[j];
		bpath[level - 1].rhs = bobj->rhs[j];
		bpath[level - 1].range = bobj->range[j];
		bpath[level - 1].branch = bobj->branch[j];
	}
#endif

	// Start with the previous bound
	if (node->feasibility_status == FEASIBLE_PRUNED ||
		node->feasibility_status == OVER_UB_PRUNED ||
		node->feasibility_status == NODE_BRANCHED_ON)
	{
		node->C_LP = node->lower_bound;
#ifdef CHECK_DUAL_SOLUTION
		for (i = 0; i < mip->m; i++)
		{
			lb += node->duals[i] * mip->rhs[i];
		}
		for (i = 0; i < mip->n; i++)
		{
			if (node->dj[i] >= 0)
			{
				lb += node->dj[i] * mip->lb[i];
			}
			else
			{
				lb += node->dj[i] * mip->ub[i];
			}
		}
		// This is just to check the lower bound value
		for (i = 0; i < level; i++)
		{
			if (bpath[i]->type == BRANCHING_VARIABLE)
			{
				switch (bpath[i]->sense)
				{
				case 'E':
					if (bpath[i]->rhs < mip->ub[j])
					{
						lb += node->dj[j] * (bpath[i]->rhs - mip->ub[j]);
					}
					else
					{
						lb += node->dj[j] * (bpath[i]->rhs - mip->lb[j]);
					}
					break;
				case 'L':
					if (node->dj[j] <= 0)
					{
						lb += node->dj[j] * (bpath[i]->rhs - mip->ub[j]);
					}
					break;
				case 'G':
					if (node->dj[j] >= 0)
					{
						lb += node->dj[j] * (bpath[i]->rhs - mip->lb[j]);
					}
					break;
				case 'R':
					printf("Warning: Ranged constraints not handled!\n");
					exit(1);
					break;
				}
			}
			else
			{ /* BRANCHING_CUT */
				printf("Warning: Branching cuts not handled!\n");
				exit(1);
			}
		}
		if (fabs(node->lower_bound - lb) >= .5)
		{
			printf("Error!\n");
		}
#endif

		// Modify
		char sense;
		for (i = 0; i < lb_cnt; i++)
		{
			if (node->dj[new_lb_ind[i]] <= 0)
			{
				continue;
			}
			for (j = 0; j < level; j++)
			{
				if (bpath[j].type == BRANCHING_VARIABLE)
				{
					k = bpath[j].name; // assuming no extra vars!
					sense = bpath[j].sense;
					if (new_lb_ind[i] == k)
					{
						if (new_lb_val[i] > bpath[j].rhs)
						{
							if (sense == 'G')
							{
								node->C_LP +=
									node->dj[k] * (new_lb_val[i] - bpath[j].rhs);
								break;
							}
							else
							{
								// New lower bound makes this node infeasible
								return (node->C_LP = node->B_IP = SYM_INFINITY);
							}
						}
						// At this point, we don't need to continue the loop
						break;
					}
				}
			}
			// This means we did not branch on new_ub_ind[i], so we adjust as usual
			if (j == level)
			{
				node->C_LP += node->dj[new_lb_ind[i]] *
							  (new_lb_val[i] - mip->lb[new_lb_ind[i]]);
			}
		}
		for (i = 0; i < ub_cnt; i++)
		{
			if (node->dj[new_ub_ind[i]] >= 0)
			{
				continue;
			}
			for (j = 0; j < level; j++)
			{
				if (bpath[j].type == BRANCHING_VARIABLE)
				{
					k = bpath[j].name; // assuming no extra vars!
					sense = bpath[j].sense;
					if (new_ub_ind[i] == k)
					{
						if (new_ub_val[i] < bpath[j].rhs)
						{
							if (sense == 'L')
							{
								node->C_LP +=
									node->dj[k] * (new_ub_val[i] - bpath[j].rhs);
								break;
							}
							else
							{
								// New lower bound makes this node infeasible
								return (node->C_LP = node->B_IP = SYM_INFINITY);
							}
						}
						// At this point, we don't need to continue the loop
						break;
					}
				}
			}
			if (j == level)
			{
				// This means we did not branch on new_ub_ind[i], so we adjust as usual
				node->C_LP += node->dj[new_ub_ind[i]] *
							  (new_ub_val[i] - mip->ub[new_ub_ind[i]]);
			}
		}
		for (i = 0; i < rhs_cnt; i++)
		{
			node->C_LP += node->duals[new_rhs_ind[i]] *
						  (new_rhs_val[i] - mip->rhs[new_rhs_ind[i]]);
		}
		node->B_IP = node->C_LP;
	}
	else if (node->feasibility_status == INFEASIBLE_PRUNED)
	{

		node->B_IP = SYM_INFINITY;
		retval = check_feasibility_new_rhs(node, mip, bpath,
										   rhs_cnt,
										   new_rhs_ind, new_rhs_val,
										   lb_cnt,
										   new_lb_ind, new_lb_val,
										   ub_cnt,
										   new_ub_ind, new_ub_val,
										   &objval);
		if (retval == LP_OPTIMAL || retval == LP_D_OBJLIM ||
			retval == LP_D_ITLIM)
		{
			node->B_IP = objval;
		}
		else
		{
			node->B_IP = SYM_INFINITY;
		}
	}
	else
	{
		printf("get_lb_for_new_rhs(): Unknown feasiblility status!\n");
		exit(1);
	}

	if (node->bobj.child_num > 0)
	{
		min = SYM_INFINITY;
	}
	for (i = 0; i < node->bobj.child_num; i++)
	{

		child = node->children[i];

		child->B_IP = get_lb_for_new_rhs(child, mip, bpath,
										 rhs_cnt, new_rhs_ind, new_rhs_val,
										 lb_cnt, new_lb_ind, new_lb_val,
										 ub_cnt, new_ub_ind, new_ub_val);

		if (child->B_IP < min)
		{
			min = child->B_IP;
		}
	}

#ifdef CHECK_DUAL_SOLUTION
	if (node->bobj.child_num == 0)
	{
		printf("%.3f %.3f %.3f\n", node->B_IP, node->lower_bound, lb);
	}
#endif

	return (node->B_IP > min ? node->B_IP : min);

#else

	printf("get_lb_for_new_rhs():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
	return (-SYM_INFINITY);

#endif
}

/*===========================================================================*/
/*===========================================================================*/

double get_ub_for_new_obj(bc_node *root, MIPdesc *mip, int cnt,
						  int *ind, double *val)
{
#ifdef SENSITIVITY_ANALYSIS
	int i, j, n;
	double inf = SYM_INFINITY;
	double min = inf;
	bc_node *child;
	double valuesi = 0.0, lpetol = 9.9999999999999995e-07;
	double objval = inf;
	double objval_pr = inf;
	double *sol = NULL;
	if (root)
	{

		for (n = 0; n < root->bobj.child_num; n++)
		{

			child = root->children[n];

			if (child->node_status == NODE_STATUS__PRUNED ||
				child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
			{

				if (child->feasibility_status == FEASIBLE_PRUNED ||
					child->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
					child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
				{

					/* too much memory to keep the OVER_UB_PRUNED solutions -
				   if sensitivity-analysis is set,
				   try the feasibility test and if succeeded then keep them*/

					if (sol)
						FREE(sol);
					sol = (double *)calloc(mip->n, DSIZE);

					objval = 0.0;
					for (j = 0; j < child->sol_size; j++)
					{
						objval += mip->obj[child->sol_ind[j]] * child->sol[j];
						sol[child->sol_ind[j]] = child->sol[j];
					}

					for (j = 0; j < cnt; j++)
					{
						objval += sol[ind[j]] *
								  (val[j] - mip->obj[ind[j]]);
					}

					if (child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
					{
						objval_pr = get_ub_for_new_obj(child, mip, cnt, ind, val);
						objval = objval < objval_pr ? objval : objval_pr;
					}
				}
			}
			else
			{
				objval = get_ub_for_new_obj(child, mip, cnt, ind, val);
			}

			if (objval < min)
			{
				min = objval;
			}
		}
	}

	if (sol)
		FREE(sol);

	return (min);
#else
	printf("get_ub_for_new_obj():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
	return (SYM_INFINITY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

/* send in a row oriented mip description! */
double get_ub_for_new_rhs(bc_node *root, MIPdesc *mip, int cnt,
						  int *ind, double *val)
{
#ifdef SENSITIVITY_ANALYSIS
	int i, j, k, n;
	double inf = SYM_INFINITY;
	double min = inf;
	bc_node *child;
	double valuesi = 0.0, lpetol = 9.9999999999999995e-07;
	double objval = inf, row_val = 0.0;
	double objval_pr;
	int feasible;
	int nonzeros;
	double *sol = NULL;

	int *matbeg = mip->matbeg, *matind = mip->matind;
	double *matval = mip->matval;

	if (root)
	{

		for (n = 0; n < root->bobj.child_num; n++)
		{

			child = root->children[n];

			if (child->node_status == NODE_STATUS__PRUNED ||
				child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
			{

				if (child->feasibility_status == FEASIBLE_PRUNED ||
					child->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
					child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
				{

					if (sol)
						FREE(sol);
					sol = (double *)calloc(mip->n, DSIZE);
					objval = 0.0;

					for (j = 0; j < child->sol_size; j++)
					{
						objval += mip->obj[child->sol_ind[j]] * child->sol[j];
						sol[child->sol_ind[j]] = child->sol[j];
					}

					/* see whether it is feasible for the new rhs! */
					feasible = TRUE;

					for (i = 0; i < cnt; i++)
					{
						row_val = 0.0;
						for (j = matbeg[ind[i]]; j < matbeg[ind[i] + 1]; j++)
						{
							row_val += matval[j] * sol[matind[j]];
						}
						switch (mip->sense[ind[i]])
						{
						case 'L':
							if (row_val > val[i] + lpetol)
							{
								feasible = FALSE;
							}
							break;
						case 'G':
							if (row_val < val[i] - lpetol)
							{
								feasible = FALSE;
							}
							break;
						case 'E':
							if (!((row_val > val[i] - lpetol) &&
								  (row_val < val[i] + lpetol)))
							{
								feasible = FALSE;
							}
							break;
						case 'R':
							if (row_val > val[i] + lpetol || row_val < val[i] - mip->rngval[i] - lpetol)
							{
								feasible = FALSE;
							}
							break;
						case 'N':
							break;
						}
						if (!feasible)
						{
							break;
						}
					}

					if (!feasible)
					{
						objval = inf;
					}

					if (child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION)
					{
						objval_pr = get_ub_for_new_rhs(child, mip, cnt, ind, val);
						objval = objval < objval_pr ? objval : objval_pr;
					}
				}
			}
			else
			{
				objval = get_ub_for_new_rhs(child, mip, cnt, ind, val);
			}

			if (objval < min)
			{
				min = objval;
			}
		}
	}

	return (min);
#else
	printf("get_ub_for_new_rhs():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
	return (SYM_INFINITY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

double check_feasibility_new_rhs(bc_node *node, MIPdesc *mip, branch_desc *bpath,
								 int rhs_cnt,
								 int *new_rhs_ind, double *new_rhs_val,
								 int lb_cnt,
								 int *new_lb_ind, double *new_lb_val,
								 int ub_cnt,
								 int *new_ub_ind, double *new_ub_val,
								 double *objval)
{
#ifdef SENSITIVITY_ANALYSIS
	int i, j;
	int level = node->bc_level;
	branch_obj *bobj;
	int retval, iterd;
	LPdata *lp_data;
	double lb, ub;

	double *old_obj = mip->obj;
	//   mip->obj = (double*)calloc(mip->n, DSIZE);

	lp_data = (LPdata *)calloc(1, sizeof(LPdata));
	lp_data->mip = mip;
	lp_data->n = mip->n;
	lp_data->m = mip->m;

	// load_the problem

	open_lp_solver(lp_data);
	load_lp_prob(lp_data, 0, 0);

	// change the rhs!

	// FIXME! change_rhs needs lp_data->tmp.c and lp_data->tmp.d???

	lp_data->tmp.c = (char *)calloc(mip->m, CSIZE);
	lp_data->tmp.d = (double *)calloc(mip->m, DSIZE);

	change_rhs(lp_data, rhs_cnt, new_rhs_ind, new_rhs_val);
	for (i = 0; i < lb_cnt; i++)
	{
		change_lbub(lp_data, new_lb_ind[i],
					new_lb_val[i], mip->lb[new_lb_ind[i]]);
	}
	for (i = 0; i < ub_cnt; i++)
	{
		change_lbub(lp_data, new_ub_ind[i],
					mip->lb[new_ub_ind[i]], new_ub_val[i]);
	}

	// add the branching changes

	for (i = 0; i < level; i++)
	{
		if (bpath[i].type == BRANCHING_VARIABLE)
		{
			j = bpath[i].name; // assuming no extra vars!
			switch (bpath[i].sense)
			{
			case 'E':
				change_lbub(lp_data, j, bpath[i].rhs, bpath[i].rhs);
				break;
			case 'L':
				get_ub(lp_data, j, &ub);
				if (bpath[i].rhs < ub)
				{
					change_ub(lp_data, j, bpath[i].rhs);
				}
				break;
			case 'G':
				get_lb(lp_data, j, &lb);
				if (bpath[i].rhs > lb)
				{
					change_lb(lp_data, j, bpath[i].rhs);
				}
				break;
			case 'R':
				change_lbub(lp_data, j, bpath[i].rhs,
							bpath[i].rhs + bpath[i].range);
				break;
			}
		}
		else
		{ /* BRANCHING_CUT */
			j = bpath[i].name;
			change_row(lp_data, j, bpath[i].sense,
					   bpath[i].rhs, bpath[i].range);
		}
	}

	// see whether it is feasible!

	size_lp_arrays(lp_data, FALSE, TRUE, mip->m, mip->n, mip->nz);
	// get_slacks call crashes here, but we can turn it off by deallocating
	FREE(lp_data->slacks);
	if ((retval = dual_simplex(lp_data, &iterd)) == LP_OPTIMAL)
	{
		*objval = lp_data->objval;
	}
	close_lp_solver(lp_data);

	lp_data->mip = NULL;
	FREE(lp_data);

#if 0
   FREE(bpath);

   for(i=0; i<(level+1); i++){
      if(path[i])
	 path[i]=NULL;
   }

   FREE(path);
#endif

	return (retval);
#else
	printf("check_feasibility_new_rhs():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
	return (0);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

char check_solution(sym_environment *env, lp_sol *sol, double *colsol)
{
	int i, j, nz = 0, *matBeg, *matInd;
	double value, *rowAct = NULL, *matVal;
	char feasible;
	double lpetol = 9.9999999999999995e-07;
	double *tmp_sol;
	double *obj = env->mip->obj;

	if (colsol)
	{
		tmp_sol = colsol;
	}
	else
	{
		tmp_sol = (double *)calloc(env->mip->n, DSIZE);
		for (i = 0; i < sol->xlength; i++)
		{
			assert(sol->xind[i] < env->mip->n);
			tmp_sol[sol->xind[i]] = sol->xval[i];
		}
	}

	sol->objval = 0;
	/* step 1. check for bounds and integrality */
	for (i = env->mip->n - 1; i >= 0; i--)
	{
		sol->objval += tmp_sol[i] * obj[i];
		if (tmp_sol[i] < env->mip->lb[i] - lpetol ||
			tmp_sol[i] > env->mip->ub[i] + lpetol)
			break;
		if (!env->mip->is_int[i])
			continue; /* Not an integer variable */
		value = tmp_sol[i];
		if (tmp_sol[i] > env->mip->lb[i] && tmp_sol[i] < env->mip->ub[i] && tmp_sol[i] - floor(tmp_sol[i]) > lpetol &&
			ceil(tmp_sol[i]) - tmp_sol[i] > lpetol)
		{
			break;
		}
	}

	feasible = i < 0 ? true : false;

	/* step 2. check for the constraint matrix */

	if (feasible)
	{
		rowAct = (double *)calloc(env->mip->m, DSIZE);
		matBeg = env->mip->matbeg;
		matVal = env->mip->matval;
		matInd = env->mip->matind;

		for (i = 0; i < env->mip->n; i++)
		{
			for (j = matBeg[i]; j < matBeg[i + 1]; j++)
			{
				rowAct[matInd[j]] += matVal[j] * tmp_sol[i];
			}
		}

		for (i = 0; i < env->mip->m; i++)
		{
			switch (env->mip->sense[i])
			{
			case 'L':
				if (rowAct[i] > env->mip->rhs[i] + lpetol)
					feasible = FALSE;
				break;
			case 'G':
				if (rowAct[i] < env->mip->rhs[i] - lpetol)
					feasible = FALSE;
				break;
			case 'E':
				if (!((rowAct[i] > env->mip->rhs[i] - lpetol) &&
					  (rowAct[i] < env->mip->rhs[i] + lpetol)))
					feasible = FALSE;
				break;
			case 'R':
				if (rowAct[i] > env->mip->rhs[i] + lpetol ||
					rowAct[i] < env->mip->rhs[i] - env->mip->rngval[i] - lpetol)
					feasible = FALSE;
				break;
			case 'N':
			default:
				break;
			}

			if (!feasible)
				break;
		}
		FREE(rowAct);
	}

	if (!colsol)
	{
		FREE(tmp_sol);
	}

	return (feasible);
}

/*===========================================================================*/
/*===========================================================================*/

int trim_warm_tree(sym_environment *env, bc_node *n)
{
	int i, not_pruned = 0;

	/* There isn't anything to do if this is a leaf. */
	if (n->bobj.child_num == 0)
		return (0);

	/* There isn't anything to do if all children are pruned, and we are
	   better off to go down if only one is not pruned. */
	for (i = n->bobj.child_num - 1; i >= 0; i--)
		if (n->children[i]->node_status != NODE_STATUS__PRUNED)
			if (++not_pruned > 1)
				break;
	if (not_pruned == 0)
		return (0);
	if (not_pruned == 1)
	{
		for (i = n->bobj.child_num - 1; i >= 0; i--)
			if (n->children[i]->node_status != NODE_STATUS__PRUNED)
			{
				trim_warm_tree(env, n->children[i]);
				break;
			}
		return (0);
	}

	/* So there are at least two not pruned. */
	for (i = n->bobj.child_num - 1; i >= 0; i--)
		if (n->children[i]->lower_bound + env->par.tm_par.granularity <
			env->warm_start->ub)
			break;

	/* if all children have high objval */
	if (i < 0)
	{
		/* get rid of the children */
		for (i = n->bobj.child_num - 1; i >= 0; i--)
			free_subtree(n->children[i]);
		/* free the children description */
		FREE(n->children);
		n->bobj.child_num = 0;
#ifndef MAX_CHILDREN_NUM
		FREE(n->bobj.sense);
		FREE(n->bobj.rhs);
		FREE(n->bobj.range);
		FREE(n->bobj.branch);
#endif
	}
	else
	{
		/* try to trim every child */
		for (i = n->bobj.child_num - 1; i >= 0; i--)
			trim_warm_tree(env, n->children[i]);
	}
	return (0);
}

/*===========================================================================*/

void free_disjunction(dual_func_desc *dual_func){

	for(int j = dual_func->num_terms - 1; j >= 0; --j){
		FREE(dual_func->disj[j].lbvaridx);
		FREE(dual_func->disj[j].lb);
		FREE(dual_func->disj[j].ubvaridx);
		FREE(dual_func->disj[j].ub);
	}
		
	FREE(dual_func->disj);
	FREE(dual_func->feas_stati);
}

void free_master(sym_environment *env)
{
	int i;
	MIPdesc *tmp;

	FREE(env->best_sol.xind);
	FREE(env->best_sol.xval);

	if ((tmp = env->mip))
	{
		free_mip_desc(env->mip);
		FREE(env->mip);
	}

	if (env->prep_mip && env->prep_mip != tmp)
	{
		free_mip_desc(env->prep_mip);
		FREE(env->prep_mip);
	}
	else
	{ // We made a copy, so don't free it again
		env->prep_mip = NULL;
	}

	if (env->rootdesc)
	{
		FREE(env->rootdesc->desc);
		FREE(env->rootdesc->uind.list);
		FREE(env->rootdesc->not_fixed.list);
		FREE(env->rootdesc->cutind.list);
		FREE(env->rootdesc);
	}

	if (env->base)
	{
		FREE(env->base->userind);
		FREE(env->base);
	}

#ifdef SYM_COMPILE_IN_TM
	if (env->warm_start)
	{
		free_subtree(env->warm_start->rootnode);
		if (env->warm_start->best_sol.has_sol)
		{
			FREE(env->warm_start->best_sol.xind);
			FREE(env->warm_start->best_sol.xval);
		}
		if (env->warm_start->cuts)
		{
			for (i = env->warm_start->cut_num - 1; i >= 0; i--)
			{
				if (env->warm_start->cuts[i])
				{
					FREE(env->warm_start->cuts[i]->coef);
				}
				FREE(env->warm_start->cuts[i]);
			}
		}

		FREE(env->warm_start->cuts);
#ifdef SENSITIVITY_ANALYSIS 
		if(env->warm_start->dual_func){

			if (env->warm_start->dual_func->duals)
				delete env->warm_start->dual_func->duals;


			if (env->warm_start->dual_func->rays){
				for (i = 0; i < env->warm_start->dual_func->num_rays; i++){
					FREE(env->warm_start->dual_func->rays[i]);
				}
				FREE(env->warm_start->dual_func->rays);
			}
				
			free_disjunction(env->warm_start->dual_func);

			dual_hash *tmp, *curr;
			HASH_ITER(hh, env->warm_start->dual_func->hashtb, curr, tmp) {
				HASH_DEL(env->warm_start->dual_func->hashtb, curr);  
				FREE(curr->basis_idx);
				FREE(curr);             
			}
			FREE(env->warm_start->dual_func->hashtb);
			if (env->warm_start->dual_func->si){
				delete env->warm_start->dual_func->si;
			} 
			FREE(env->warm_start->dual_func);
		}
#endif
		FREE(env->warm_start);
	}
#ifdef SYM_COMPILE_IN_CP
	if (env->cp)
	{
		for (i = 0; i < env->par.tm_par.max_cp_num; i++)
		{
			env->cp[i]->msgtag = YOU_CAN_DIE;
			cp_close(env->cp[i]);
		}
		FREE(env->cp);
	}
#endif
#ifdef SYM_COMPILE_IN_LP
	if (env->sp)
	{
		sp_free_sp(env->sp);
		FREE(env->sp);
	}
#endif
#endif

	return;
}

/*===========================================================================*/
/*===========================================================================*/

// feb223

int get_num_leaf_nodes(bc_node *node, int *num_rays)
{
	if (node == NULL)
	{
		return 0;
	}
	else if (!node->bobj.child_num)
	{
		// Count how many rays are there
		if (
			(node->rays) 
			// && 
			// (node->feasibility_status == INFEASIBLE_PRUNED) ||
			// (node->feasibility_status == OVER_UB_PRUNED)
			){
			(*num_rays)++;
		}
		return 1;
	}
	else
	{
		int num_leaf_subtrees = 0;

		if ((node->rays) 
			// && 
			// (node->feasibility_status == INFEASIBLE_PRUNED) ||
			// (node->feasibility_status == OVER_UB_PRUNED)
			){
			(*num_rays)++;
		}

		for (int i = 0; i < node->bobj.child_num; i++)
		{
			num_leaf_subtrees += get_num_leaf_nodes(node->children[i], num_rays);
		}
		return num_leaf_subtrees;
	}
	printf("Warning: Error from get_num_leaf_nodes()\n");
	return (FUNCTION_TERMINATED_ABNORMALLY);
}

/*===========================================================================*/

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

void printDisjunction(disjunction_desc disj){
	printf("LBs: ");
    for(int i = 0; i < disj.lblen; i++){
        printf("%d : %.2f, ", disj.lbvaridx[i], disj.lb[i]);
    }
    printf("\n");
	printf("UBs: ");
    for(int i = 0; i < disj.ublen; i++){
        printf("%d : %.2f, ", disj.ubvaridx[i], disj.ub[i]);
    }
	printf("\n");
	printf("--------------------\n");
}

void print_dual_function(warm_start_desc *ws)
{
	int k = 0;
	int tot_duals = 0;

	if (!ws || !ws->dual_func){
		printf("WARNING: No dual function in the current warm start!\n");
		return;
	}

	printf("====================================\n");
	printf("  DUAL FUNCTION INFO\n");
	printf("====================================\n");
	if (ws->dual_func->policy == DUALS_LEAF_ONLY)
		printf("Policy: DUALS_LEAF_ONLY\n");
	else if (ws->dual_func->policy == DUALS_SAVE_ALL)
		printf("Policy: DUALS_SAVE_ALL\n");
	printf("Granularity: %.10f\n", ws->dual_func->granularity);
	printf("Num rays: %d\n", ws->dual_func->num_rays);
	printf("Num dual sol: %d\n", ws->dual_func->num_pieces);
	printf("Num terms: %d\n", ws->dual_func->num_terms);
	printf("CPU Time spent in LPs: %.5f\n", ws->dual_func->lp_cpu_time);
	printf("====================================\n");

	// printf("==========================\n");
	// printf("RAYS\n");
	// printf("==========================\n");
	// if (ws->dual_func->num_rays == 0){
	// 	printf("No rays!");
	// 	printf("\n------------------------\n");
	// } else {
	// 	for (int i = 0; i < ws->dual_func->num_rays; i++){
	// 		for (int j = 0; j < ws->m + ws->n; j++){
	// 			printf("%.10f ", ws->dual_func->rays[i][j]);
	// 		}
	// 		printf("\n------------------------\n");
	// 	}
	// }

	// CoinPackedMatrix *duals = ws->dual_func->duals;
	// printf("===============================\n");
	// printf("DUAL SOLUTIONS + REDUCED COSTS \n");
	// printf("===============================\n");

    // const double* elem = duals->getElements();
    // const int* indices = duals->getIndices();

	// int major = duals->getMajorDim();
	// int minor = duals->getMinorDim();
	// int num_elem = duals->getNumElements();
	// int peppe = duals->getVectorLast(0);
	// int ord = duals->isColOrdered();
	// if (major == 0){
	// 	printf("No Dual Solutions!");
	// 	printf("\n------------------------\n");
	// } else {
	// 	for (int i = 0; i < major; i++){
	// 		const CoinBigIndex last = duals->getVectorLast(i);
	// 		for (CoinBigIndex j = duals->getVectorFirst(i); j < last; ++j){
	// 			// printf("%d, %d - %d \n", j, last, indices[j]);
	// 			printf("%d : %.10f, ", indices[j] < ws->m ? indices[j] : indices[j] /*- ws->m*/, elem[j]);
	// 		}
	// 		printf("\n------------------------\n");
   	// 	}
	// }
    
	
	// printf("==========================\n");
	// printf("DISJUNCTION\n");
	// printf("==========================\n");
	// if (ws->dual_func->num_terms == 0){
	// 	printf("Only root node in the B&B Tree!\n");
	// 	printf("Feasibility Status: ");
	// 		switch (ws->rootnode->feasibility_status)
	// 		{
	// 		case ROOT_NODE:
	// 			printf("ROOT_NODE\n");
	// 			break;
			
	// 		case INFEASIBLE_PRUNED:
	// 			printf("INFEASIBLE_PRUNED\n");
	// 			break;
			
	// 		case FEASIBLE_PRUNED:
	// 			printf("FEASIBLE_PRUNED\n");
	// 			break;

	// 		case OVER_UB_PRUNED:
	// 			printf("OVER_UB_PRUNED\n");
	// 			break;
			
	// 		default:
	// 			printf("UNKNOWN\n");
	// 			break;
	// 		}
	// } else {
	// 	for (int i = 0; i < ws->dual_func->num_terms; i++)
	// 	{
	// 		printf("Feasibility Status: ");
	// 		switch (ws->dual_func->feas_stati[i])
	// 		{
	// 		case ROOT_NODE:
	// 			printf("ROOT_NODE\n");
	// 			break;
			
	// 		case INFEASIBLE_PRUNED:
	// 			printf("INFEASIBLE_PRUNED\n");
	// 			break;
			
	// 		case FEASIBLE_PRUNED:
	// 			printf("FEASIBLE_PRUNED\n");
	// 			break;

	// 		case OVER_UB_PRUNED:
	// 			printf("OVER_UB_PRUNED\n");
	// 			break;
			
	// 		default:
	// 			printf("UNKNOWN\n");
	// 			break;
	// 		}
	// 		// printDisjunction(ws->dual_func->disj[i]);
	// 	}
	// }
	

	// printf("==========================\n");
	// printf("BASIS\n");
	// printf("==========================\n");
	// for (dual_hash *s = ws->dual_func->hashtb; s != NULL; s = (dual_hash *)(s->hh.next)) {
    //     for (int i = 0; i < s->len; i++){
	// 		printf("%d ", s->basis_idx[i]);
	// 	}
	// 	printf("\n");
    // }
}

/*===========================================================================*/

int add_dual_to_table(dual_hash **hashtb, dual_hash *toAdd){
    dual_hash *s;
	int is_added = 0;
	int len = sizeof(int) * toAdd->len;
	// printf("Basis toAdd: ");
	// for (int i = 0; i < toAdd->len; i++){
	// 	printf("%d ", toAdd->basis_idx[i]);
	// }
	// printf("\n--------\n");
    HASH_FIND(hh, *hashtb, toAdd->basis_idx, len, s);
    if (s == NULL){
        s = toAdd;
		// TODO: FREE this
        HASH_ADD_KEYPTR(hh, *hashtb, s->basis_idx, len, s);
		is_added = 1;
    } else {
        is_added = 0;
    }
	// FREE(s);
	return is_added;
}

// int count_leaf = 0;
// int optimal_leaf_found = 0;

#ifdef CHECK_DUAL_FUNC
	int count_leaf = 0;
#endif

void collect_duals(sym_environment *env, bc_node *node, MIPdesc *mip,
				  branch_desc *bpath, int* curr_term, int* curr_piece, int* curr_ray, double** rays,
				  int *duals_index_row, int *duals_index_col, double *duals_val,
				  int *nnz_duals, int *max_poss_nnz, int *curr_leaf){

	if (node == NULL)
	{
		printf("Warning: NULL pointer in get_leaf_node_data()\n");
		return;
	}
	const CoinPackedMatrix *A;
	// if (!node->parent) count_leaf = 0;

	warm_start_desc *ws = env->warm_start;

	int i, j, k, l;
	int is_added = 0;
	int varidx, nnz = 0;
	int lbchange = 0, ubchange = 0; 
	int level = node->bc_level, child_num = node->bobj.child_num;

	double     dualobj = 0;
	double ray_times_b = 0;
	double     zerotol = 1e-9;

	double *lb = NULL;
	double *ub = NULL;

	branch_obj *bobj;
	dual_hash *dual;
	disjunction_desc disj;
	
	// if level>0, save the branching constraint
	if (level > 0)
	{
		for (j = 0; j < node->parent->bobj.child_num; j++)
			if (node->parent->children[j] == node)
				break;

		bobj = &(node->parent->bobj);
		
		// Save branching info
		if (bobj->type == BRANCHING_VARIABLE){
			bobj = &(node->parent->bobj);
			bpath[level - 1].type = bobj->type;
			bpath[level - 1].name = bobj->name;
			bpath[level - 1].sense = bobj->sense[j];
			bpath[level - 1].rhs = bobj->rhs[j];
			bpath[level - 1].range = bobj->range[j];
			bpath[level - 1].branch = bobj->branch[j];
		} else {
			// BRANCHING_CUT
			printf("Branching type not handled!\n");
			exit(1);
		}
	}

	// TODO: deal with ITERATION_LIMIT and TIME_LIMIT
	// check feasibility status of this node
	if (
		// (!is_lower_bound_inherited) && 
	   (node->feasibility_status == ROOT_NODE ||
		node->feasibility_status == FEASIBLE_PRUNED ||
		node->feasibility_status == OVER_UB_PRUNED ||
		node->feasibility_status == NODE_BRANCHED_ON )
		){
		if (node->duals && node->dj && (!node->rays) &&
		((ws->dual_func->policy == DUALS_SAVE_ALL) ||
		 (ws->dual_func->policy == DUALS_LEAF_ONLY && !child_num)))
		{
			dual = NULL;
			if (node->basis_idx && (node->basis_len > 0)){
				dual = (dual_hash *)malloc(sizeof(dual_hash));
				dual->basis_idx = (int *)malloc(ISIZE * node->basis_len);
				memcpy(dual->basis_idx, node->basis_idx, ISIZE * node->basis_len);
				dual->len = node->basis_len;
				// COMMENT THIS when using the hash table
				// if (dual){
				// 	FREE(dual->basis_idx);
				// 	FREE(dual);
				// }
			} 
			// try to add this dual into the hashtable
			// TODO: are the basis_idx correctly set??? Are they unique???
			if ((dual && (is_added = add_dual_to_table(&(ws->dual_func->hashtb), dual)))){
				// successfully added, collect reduced costs
				// find nnzs and fill reduced costs related structures
				for (i = 0; i < ws->m; i++)
				{
					(*max_poss_nnz)++;
					if (fabs(node->duals[i]) >= zerotol)
					{
						duals_index_row[*nnz_duals] = *curr_piece;
						duals_index_col[*nnz_duals] = i;
						duals_val[*nnz_duals] = node->duals[i];
						(*nnz_duals)++;
					}
				}
				
				for (i = 0; i < ws->n; i++)
				{
					(*max_poss_nnz)++;
					if (fabs(node->dj[i]) >= zerotol)
					{
						duals_index_row[*nnz_duals] = *curr_piece;
						duals_index_col[*nnz_duals] = ws->m + i;
						duals_val[*nnz_duals] = node->dj[i];
						(*nnz_duals)++;
					}
				}
				(*curr_piece)++;

			} else {
				// This dual is a duplicate, we can FREE it
				if (dual){
					FREE(dual->basis_idx);
					FREE(dual);
				}
				
			}
		}
	}

	// If this is a child, we have a term of the disjunction
	if (!child_num)
	{
		ws->dual_func->feas_stati[(*curr_leaf)++] = node->feasibility_status; 
		// printf(" - %d Fesibility status of this leaf: %d\n", count_leaf, node->feasibility_status);
		// printf("   %d Lower bound at this leaf: %.5f\n", count_leaf, node->lower_bound);
		if (level > 0)
		{
			disj.lb = NULL;
			disj.lbvaridx = NULL;
			disj.ub = NULL;
			disj.ubvaridx = NULL;
			disj.lblen = disj.ublen = 0;   
			lb = (double *)malloc(DSIZE * ws->n);
			ub = (double *)malloc(DSIZE * ws->n);
			memcpy(lb, mip->lb, DSIZE * ws->n);
			memcpy(ub, mip->ub, DSIZE * ws->n);
			for (int i = 0; i < level; i++){
				varidx = bpath[i].name;
				switch (bpath[i].sense)
				{
				case 'E':
					lb[varidx] = bpath[i].rhs;
					ub[varidx] = bpath[i].rhs;
					break;
				case 'L':
					// if (bpath[i].rhs < mip->ub[varidx])
						ub[varidx] = bpath[i].rhs;
					break;
				case 'G':
					// if (bpath[i].rhs > mip->lb[varidx])
						lb[varidx] = bpath[i].rhs;
					break;
				case 'R':
					printf("Warning: Ranged constraints not handled!\n");
					exit(1);
					break;
				}
			}

			// count how many lb/ub changed from the original mip
			for (int i = 0; i < ws->n; i++){
				if (lb[i] > mip->lb[i])
					lbchange++;
				if (ub[i] < mip->ub[i])
					ubchange++;	
			}

			disj.lblen = lbchange;
			disj.ublen = ubchange;
		
			if (lbchange > 0){
				disj.lbvaridx = (int *)malloc(sizeof(int) * lbchange);
				disj.lb = (double *)malloc(sizeof(double) * lbchange);
			}

			if (ubchange > 0){
				disj.ubvaridx = (int *)malloc(sizeof(int) * ubchange);
				disj.ub = (double *)malloc(sizeof(double) * ubchange);
			}
			
			// fill
			l = k = 0;
			for (int i = 0; i < ws->n; i++){
				if (lb[i] > mip->lb[i]){
					disj.lbvaridx[l] = i;
					disj.lb[l] = lb[i];
					l++;
				}
				if (ub[i] < mip->ub[i]){
					disj.ubvaridx[k] = i;
					disj.ub[k] = ub[i];
					k++;
				}
			}

			ws->dual_func->disj[*curr_term] = disj;

			(*curr_term)++;

			// Must collect the ray
			if ((node->feasibility_status == INFEASIBLE_PRUNED) ||
				((node->feasibility_status == OVER_UB_PRUNED) && (node->rays))){
				if ((node->feasibility_status == INFEASIBLE_PRUNED) && !node->rays){
					// This should never happen now, but if it does, resolve the LP
					// and get the ray
#ifdef CHECK_DUAL_FUNC
					printf(" WARNING in collect_duals():\n");
					printf(" Infeasible node without rays! Solving LP...\n");
#endif
					clock_t start, end;
					double cpu_time_used;
					start = clock();
					if (!ws->dual_func->si){
						ws->dual_func->si = new OsiXSolverInterface();
						ws->dual_func->si->setupForRepeatedUse();
						ws->dual_func->si->setHintParam(OsiDoReducePrint);
						ws->dual_func->si->messageHandler()->setLogLevel(0);
						ws->dual_func->si->setCleanupScaling(1);
					}
					ws->dual_func->si->loadProblem(ws->n, ws->m,
							mip->matbeg, mip->matind,
							mip->matval, lb,
							ub, mip->obj,
							mip->sense, mip->rhs,
							mip->rngval);

					ws->dual_func->si->initialSolve();

					if (ws->dual_func->si->isProvenPrimalInfeasible()){
						std::vector<double *> vRays = ws->dual_func->si->getDualRays(1, false);
						// check that there is at least one ray
						int raysReturned = static_cast<unsigned int>(vRays.size());
						assert(raysReturned == 1);	
						if (vRays[0])
						{
							double *ray = vRays[0];
							int i;
							bool is_ray_empty = TRUE;

							const double *inverseRowScale = ws->dual_func->si->getModelPtr()->inverseRowScale();
							const double *rowScale = ws->dual_func->si->getModelPtr()->rowScale();
							double *rayA = (double *)malloc(sizeof(double) * ws->n);

							A = ws->dual_func->si->getMatrixByCol();
							
							double norm = 0;

							for (i = 0; i < ws->m; i++)
							{
								if (fabs(ray[i]) > 1e-5){
									is_ray_empty = FALSE;
									break;
								}
							} 

							if (!is_ray_empty)
							{	
								for (i = 0; i < ws->m; i++)
								{
									norm += ray[i] * ray[i];
								}
								norm = sqrt(norm);

								if (norm > 1e-7){
									for (i = 0; i < ws->m; i++)
									{
										ray[i] /= norm;
									}
								}
									
								A->transposeTimes(ray, rayA);
								node->rays = (double *)malloc(sizeof(double) * (ws->n + ws->m));
								memcpy(node->rays, ray, (ws->m) * DSIZE);
								memcpy(node->rays + ws->m, rayA, (ws->n) * DSIZE);
							} else {
#ifdef CHECK_DUAL_FUNC
								printf(" WARNING in collect_duals():\n");
								printf("  Ray is all zeros!\n");
#endif
							}
							
							delete[] vRays[0];
							FREE(rayA);
						}
						else
						{
#ifdef CHECK_DUAL_FUNC
							printf(" WARNING in collect_duals():\n");
							printf("  Still no ray after LP solve!\n");
#endif
						}

					} else {
#ifdef CHECK_DUAL_FUNC
						printf(" WARNING in collect_duals():\n");
						printf("  This node should not be infeasible!\n");
#endif
					}	
					end = clock();
					cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
					ws->dual_func->lp_cpu_time += cpu_time_used;
				}
				if (node->rays){
					bool is_duplicate = false;
					bool is_all_zero = true;
					for (int i = 0; i < ws->m; i++){
						ray_times_b -= node->rays[i] * mip->rhs[i];
					}

					for (int i = ws->m; i < ws->m + ws->n; i++){
						if (node->rays[i] > zerotol){
							ray_times_b += node->rays[i] * lb[i - ws->m];
						} 
						else if (node->rays[i] < -zerotol){
							ray_times_b += node->rays[i] * ub[i - ws->m];
						}
					}

					assert(ray_times_b > 1e-5);
					// Check that the ray is not all zeros
					for (i = 0; i < ws->m + ws->n; i++)
					{
						if (fabs(node->rays[i]) > 1e-5){
							is_all_zero = false;
							break;
						}
					}

					if (!is_all_zero){
						// Check duplicates among current rays
						for (i = 0; i < (*curr_ray) && (!is_duplicate); i++){
							for (j = 0; j < ws->m; j++){
								if (fabs(node->rays[j] - rays[i][j]) > 1e-7){
									// These rays are different, no need to check more
									break;
								}
								if (j == ws->m - 1)
									// We already have this ray
									is_duplicate = true;
							}
						}

						// Check duplicates among rays in dual_func_desc
						for (i = 0; i < ws->dual_func->num_rays && (!is_duplicate); i++){
							for (j = 0; j < ws->m; j++){
								if (fabs(node->rays[j] - ws->dual_func->rays[i][j]) > 1e-7){
									// These rays are different, no need to check more
									break;
								}
								if (j == ws->m - 1)
									// We already have this ray
									is_duplicate = true;
							}
						}

						if (!is_duplicate){
							// Allocate space for curr ray
							rays[*curr_ray] = (double *)malloc(DSIZE * (ws->m + ws->n));
							memcpy(rays[*curr_ray], node->rays, DSIZE * (ws->m + ws->n));
							(*curr_ray)++;
						}
					}
				}
			}

#ifdef CHECK_DUAL_FUNC
			// Let's do some checks
			// Dual Obj Value of this node 
			if (!lb)
				lb = (double *)malloc(DSIZE * ws->n);

			if (!ub)
				ub = (double *)malloc(DSIZE * ws->n);

			memcpy(lb, mip->lb, DSIZE * ws->n);
			memcpy(ub, mip->ub, DSIZE * ws->n);

			change_lbub_from_disj(lb, ub, ws->n, ws->dual_func->disj + ((*curr_term) - 1));
			
			dualobj = 0;
			for (int i = 0; i < ws->m; i++){
				dualobj += mip->rhs[i] * node->duals[i];
			}
			for (int i = 0; i < ws->n; i++){
				if (node->dj[i] > 1e-5)
					dualobj += lb[i] * node->dj[i];
				else if (node->dj[i] < -1e-5)
					dualobj += ub[i] * node->dj[i];
			}
			// if (count_leaf == 85)
			// 	if (node->rays)
			// 		printf("Has ray! %.5f\n", ray_times_b);
			printf(" - %d Fesibility status of this leaf: %d\n", count_leaf, node->feasibility_status);
			printf("   %d Lower bound at this leaf: %.5f\n", count_leaf, node->lower_bound);
			printf("   %d Dual Obj Val at this leaf: %.5f\n", count_leaf++, dualobj);
			reset_lbub_from_disj(lb, ub, mip->lb, mip->ub, ws->n, ws->dual_func->disj + ((*curr_term) - 1));
			if (node->feasibility_status == INFEASIBLE_PRUNED){
				// There must be a ray checking its primal infeasibility
				if (node->rays){
					if (ray_times_b < 1e-5){
						printf(" WARNING in collect_duals():\n");
						printf("  INFEASIBLE_PRUNED: Farkas proof at this node doesn't work!\n");
					}
				} else {
					printf(" WARNING in collect_duals():\n");
					printf("  No ray at this node!\n");
				}
			} else if (node->feasibility_status == OVER_UB_PRUNED){
				// Either there is a ray proving unboundedness...
				if (node->rays){
					if (ray_times_b < 1e-5){
						printf(" WARNING in collect_duals():\n");
						printf("  OVER_UB_PRUNED: Farkas proof at this node doesn't work!\n");
					}
				} else {
					// ... or the dual obj value is over the UB
					if (dualobj + ws->dual_func->granularity - ws->ub < 1e-5){
						printf(" WARNING in collect_duals():\n");
						printf( "bc node %d :", node->bc_index);
						printf("  OVER_UB_PRUNED: dual obj val not over UB!\n");
						printf("  dualobj = %.5f < %.5f = UB\n", dualobj, ws->ub);
					}
				}	
			} 
			

#endif
			FREE(lb);
			FREE(ub);
		}
	} else {
		// if child_num > 0, then do recursion on child nodes
		for (j = 0; j < child_num; j++)
			collect_duals(env, node->children[j], mip, bpath, curr_term, curr_piece, curr_ray, rays,
						duals_index_row, duals_index_col, duals_val, nnz_duals, max_poss_nnz, curr_leaf);

	}
}

int build_dual_func(sym_environment *env)
{

#ifdef SENSITIVITY_ANALYSIS
#ifdef CHECK_DUAL_FUNC
	count_leaf = 0;
#endif
	warm_start_desc *ws = env->warm_start;
	MIPdesc *mip = env->mip;
	ws->m = mip->m;
	ws->n = mip->n;

	if (!ws || !env)
	{
		printf("Warning: NULL pointer in build_dual_func()\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	int num_rays = 0;
	int num_leaf = get_num_leaf_nodes(ws->rootnode, &num_rays);

	// printf("There are %d rays in the current tree!\n", num_rays);

 	if (!ws->dual_func)
	{
		ws->dual_func = (dual_func_desc *)calloc(1, sizeof(dual_func_desc));
		ws->dual_func->hashtb = NULL;
		ws->dual_func->duals  = NULL;
		ws->dual_func->num_pieces = 0;
		ws->dual_func->policy = DUALS_SAVE_ALL;
		// ws->dual_func->policy = DUALS_LEAF_ONLY;
		ws->dual_func->granularity = env->par.tm_par.granularity;
		ws->dual_func->rays = NULL;
		ws->dual_func->num_rays = 0;
		ws->dual_func->disj = NULL;
		ws->dual_func->feas_stati = NULL;
		ws->dual_func->num_terms = 0;
	} else {
		// Delete the previous disjunction, there will be a new one
		free_disjunction(ws->dual_func);
	}
	
	ws->dual_func->num_terms = num_leaf;
	// ws->dual_func->num_rays = num_rays;

	// Allocate the space for the expected number of duals
	// based on the policy
	int num_pieces;
	if ((ws->dual_func->policy == DUALS_SAVE_ALL))
	{
		num_pieces = ws->stat.tree_size;
	}
	else if (ws->dual_func->policy == DUALS_LEAF_ONLY)
	{
		num_pieces = num_leaf;
	}
	else
	{
		printf("build_dual_func():\n");
		printf("Policy not supported.\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	// Allocate memory for disjunction
	ws->dual_func->disj = (disjunction_desc*)malloc(sizeof(disjunction_desc) * num_leaf);
	ws->dual_func->feas_stati = (int*)malloc(sizeof(int) * num_leaf);

	// Allocate memory for rays to be collected
	double **rays = (double**)malloc(sizeof(double*) * num_rays);

	// Duals and reduced costs
	int nnz_duals = 0;
	int *duals_index_row = (int *)malloc(ISIZE * num_pieces * (ws->m + ws->n));
	int *duals_index_col = (int *)malloc(ISIZE * num_pieces * (ws->m + ws->n));
	double *duals_val = (double *)malloc(DSIZE * num_pieces * (ws->m + ws->n));

	// branching path
	branch_desc *bpath = (branch_desc *)malloc(sizeof(branch_desc) * (ws->stat.max_depth));
	int len = 0;

	int curr_term = 0;  // current disjunction term
	int curr_piece = 0; // current dual piece
	int curr_ray = 0;   // current ray
	int tot_piece = 0; // DEBUG: to check how sparse it is

	int curr_leaf = 0; // Debug: check the stati of leaves

	collect_duals(env, ws->rootnode, mip, bpath, &curr_term, &curr_piece, &curr_ray, rays,
				 duals_index_row, duals_index_col, duals_val, &nnz_duals, &tot_piece, &curr_leaf);

	// Update the actual number of disjunction terms
	// Infeasible leaf lead to infeasible disjunction term
	ws->dual_func->num_terms = curr_term;
#ifdef CHECK_DUAL_FUNC
	if (num_leaf != curr_term){
		printf(" WARNING in build_dual_func():\n");
		printf(" Not all leaves have been collected! num_leaf = %d, dual_func->num_terms %d\n", 
		num_leaf, ws->dual_func->num_terms);
	}
#endif
	// if (nnz_duals)
	// 	printf("DEBUG: matrix sparsity %.3f\n", ((double)(nnz_duals)/(double)(tot_piece)));
	// else
	// 	printf("DEBUG: No new dual pieces\n");

	// If at least one new dual have been found, construct the Dual Function
	if (curr_piece > 0)
	{
		CoinPackedMatrix *df = new CoinPackedMatrix(false, duals_index_row,
									duals_index_col, duals_val, nnz_duals);

		// Fix the number of columns to m + n
		df->setDimensions(-1, ws->m + ws->n);

		if (!ws->dual_func->duals)
		{
			// New Dual Function
			ws->dual_func->duals = df;
		}
		else
		{
			// Already existing Dual Function
			// We append to the previous, duplicates have been already checked
			ws->dual_func->duals->bottomAppendPackedMatrix(*df);
			delete df;
		}

		// There should be a limit on the number of dual solutions we can collect
		// ws->dual_func->num_pieces = ws->dual_func->duals->getMajorDim();
		ws->dual_func->num_pieces += curr_piece;
#ifdef CHECK_DUAL_FUNC
		if (ws->dual_func->num_pieces !=  ws->dual_func->duals->getMajorDim()){
			printf(" WARNING in build_dual_func():\n");
			printf(" Some pieces are lost in CoinPackedMatrix! dual_func->num_pieces = %d, duals->getMajorDim() = %d\n", 
			ws->dual_func->num_pieces, ws->dual_func->duals->getMajorDim());
		}	
		int len_ht = 0;
		dual_hash *s;
		for (s = ws->dual_func->hashtb; s != NULL; s = (dual_hash*)s->hh.next){
			len_ht++;
		}

		if (len_ht != ws->dual_func->num_pieces){
			printf(" WARNING in build_dual_func():\n");
			printf(" Some pieces are lost in HashTable! dual_func->num_pieces = %d, len_hashtable = %d\n", 
			ws->dual_func->num_pieces, len_ht);
		}	
#endif
	}

	if (curr_ray > 0){
		// Allocate space for new rays
		ws->dual_func->rays = (double **)realloc(ws->dual_func->rays, 
								sizeof(double*) * (ws->dual_func->num_rays + curr_ray));
		
		if (ws->dual_func->rays == NULL){
			printf("build_dual_func():\n");
			printf("realloc() error.\n");
		} else {

			// Move new rays into dual_func_desc
			for (int i = 0; i < curr_ray; i++){
				ws->dual_func->rays[ws->dual_func->num_rays + i] = 
					(double *)malloc(DSIZE * (ws->m + ws->n));
				memcpy(ws->dual_func->rays[ws->dual_func->num_rays + i], 
						rays[i], (ws->m + ws->n) * DSIZE);
			}

			ws->dual_func->num_rays += curr_ray;
		}
	}

	for (int i = 0; i < curr_ray; i++){
		FREE(rays[i]);
	}
	FREE(rays);
	FREE(bpath);
	FREE(duals_index_row);
	FREE(duals_index_col);
	FREE(duals_val);

	// printf("==========================\n");
	// printf("RAYS\n");
	// printf("==========================\n");
	// for (int i = 0; i < ws->dual_func->num_rays; i++){
	// 	for (int j = 0; j < ws->m + ws->n; j++){
	// 		printf("%.10f ", ws->dual_func->rays[i][j]);
	// 	}
	// 	printf("\n------------------------\n");
	// }

	// print_dual_function(ws);
	return (FUNCTION_TERMINATED_NORMALLY);
#else
	printf("build_dual_func():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
	return (FUNCTION_TERMINATED_ABNORMALLY);

#endif
}

int evaluate_dual_function(warm_start_desc *ws, MIPdesc *mip, 
						double *new_rhs, int size_new_rhs, 
						double *dual_bound){
#ifdef SENSITIVITY_ANALYSIS
	if (!ws)
	{
		printf("Warning: NULL pointer in evaluate_dual_func()\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	if (!ws->dual_func)
	{
		printf("Warning: None or empty dual function in evaluate_dual_func()\n");
		return (FUNCTION_TERMINATED_ABNORMALLY);
	}

	if (ws->dual_func->num_pieces == 0){
		*dual_bound = -SYM_INFINITY;
		return (FUNCTION_TERMINATED_NORMALLY);
	}

#ifdef DEBUG_DUAL_FUNC
	// For sanity check
	double sanity_objVal = 0;
	OsiXSolverInterface *si;
	if (!ws->dual_func->si){
		ws->dual_func->si = new OsiXSolverInterface();
		ws->dual_func->si->setupForRepeatedUse();
		ws->dual_func->si->setHintParam(OsiDoReducePrint);
		ws->dual_func->si->messageHandler()->setLogLevel(0);
		ws->dual_func->si->setCleanupScaling(1);
	} 
	si = ws->dual_func->si;
#endif

	// Build the new rhs to be of the correct size, ws->m
	double *rhs = (double *)malloc(ws->m * sizeof(double));

	memcpy(rhs, new_rhs, CoinMin(size_new_rhs, ws->m) * sizeof(double));
	if (size_new_rhs < ws->m){
		for (int i = size_new_rhs; i < ws->m; i++){
			rhs[i] = mip->rhs[i];
		}
	}

	int i, u, l, idx, curr_is_infty;
	CoinPackedMatrix *duals = ws->dual_func->duals;
	disjunction_desc *disj;
	CoinBigIndex first = 0, last = 0, j = 0;

	double zerotol = 1e-9;
	double granularity = ws->dual_func->granularity;
	double global_best_bound = SYM_INFINITY;
	double local_best_bound;
	double *ray;
	double ray_times_b = 0;
	int    is_ray_infty = 0;
	bool *is_term_feas = NULL;
	double rhs_times_pi_plus_dj = 0;
	int is_infty = 0;
	
	// Allocate space for rhs * pi (i.e. dual solution)
	double *rhs_times_pi = (double *)calloc(ws->dual_func->num_pieces, DSIZE);
	int *dj_start = (int *)malloc(ws->dual_func->num_pieces * ISIZE);
	double *lb = (double *)malloc(ws->n * sizeof(double));
	double *ub = (double *)malloc(ws->n * sizeof(double)); 

	// Root lb/ub
	memcpy(lb, mip->lb, DSIZE * ws->n);
	memcpy(ub, mip->ub, DSIZE * ws->n);

	const double* elem = duals->getElements();
    const int* indices = duals->getIndices();

	// Start by computing rhs * pi 
	for (i = 0; i < ws->dual_func->num_pieces; i++){
		// if (i == 173){
		// 	printf("here!");
		// }
		first = duals->getVectorFirst(i);
		last = duals->getVectorLast(i);
		for (j = first; j < last && indices[j] < ws->m; ++j){
			rhs_times_pi[i] += elem[j] * rhs[indices[j]];
		}
		// the last value of j is where the djs starts for the curr dual piece
		if (j < last){
			dj_start[i] = j;	
		} else {
			// this should only happen when all djs are zero
			dj_start[i] = last;
		}
	}

	if (ws->dual_func->num_terms == 0){
		// The current b&b tree consists just of the root node
		local_best_bound = -SYM_INFINITY;
		// Now add the lb/ub * dj from the original MIP
		for (i = 0; i < ws->dual_func->num_pieces; i++){
			rhs_times_pi_plus_dj = rhs_times_pi[i];
			is_infty = 0;
			last = duals->getVectorLast(i);
			for (j = dj_start[i]; j < last; j++)
			{
				idx = indices[j] - ws->m;
				if (elem[j] >= zerotol) // dj > 0 ?
				{
					if (lb[idx] > -SYM_INFINITY){
						rhs_times_pi_plus_dj += elem[j] * lb[idx];
					} else {
						is_infty -= elem[j];
					}
				}
				else if (elem[j] <= -zerotol) // dj < 0 ?
				{
					if (ub[idx] < SYM_INFINITY){
						rhs_times_pi_plus_dj += elem[j] * ub[idx];
					} else {
						is_infty -= elem[j];	
					}
				}
			}

			if ((is_infty == 0) && (rhs_times_pi_plus_dj > local_best_bound)){
				local_best_bound = rhs_times_pi_plus_dj;
			}
		}
		
		global_best_bound = local_best_bound;
		goto TERM_EVAL_DUAL_FUNC;
	}

	// =========== feb223 ====================================
	// Use the dual rays in the Dual Function description to check
	// if they yield a Farkas proof of primal infeasibility 
	// for some term of the disjunction.
	// 
	// Given the LP:
	// min 		cx
	// s.t.		Ax <= b,
	//		l <= x <= u,			
	// and let r be a dual ray, then the Farkas proof goes as
	// 
	// 				rAz - rb > 0, where
	// 
	// z[j]  =	l[j], if rA[j] > 0,
	// 	     =	u[j], if rA[j] < 0 
	// yield the minimum violation.
	// 
	// CLP returns by default (r, -rA) but r is for the scaled 
	// LP it is internally solving and it may not prove infeasibility
	// for the unscaled problem. This has already taken care
	// of in get_dual_ray().
	// =========================================================

	is_term_feas = (bool *)malloc(ws->dual_func->num_terms * sizeof(bool));
	CoinFillN(is_term_feas, ws->dual_func->num_terms, true);
	
	if (ws->dual_func->num_rays > 0){
		for (int t = 0; t < ws->dual_func->num_terms; t++){

			disj = ws->dual_func->disj + t;
			// Install lb/ub
			change_lbub_from_disj(lb, ub, ws->n, disj);

			for (int r = 0; (r < ws->dual_func->num_rays) && (is_term_feas[t]); r++){
				ray_times_b = 0;
				ray = ws->dual_func->rays[r];
				for (int i = 0; i < ws->m; i++){
					ray_times_b -= ray[i] * rhs[i];
				}

				for (int i = ws->m; i < ws->m + ws->n; i++){
					if (ray[i] > zerotol){
						ray_times_b += ray[i] * lb[i - ws->m];
					} 
					else if (ray[i] < -zerotol){
						ray_times_b += ray[i] * ub[i - ws->m];
					}
				}

				if (ray_times_b > 1e-5){
					// This term of the disjunction is infeasible for this rhs
					// and should not be considered
					is_term_feas[t] = false;
				}
			}

			// Reset lb/ub
			reset_lbub_from_disj(lb, ub, mip->lb, mip->ub, ws->n, disj);
		}
	}

#ifdef DEBUG_DUAL_FUNC
	// This is a sanity check if the rays prove primal infeasibility
	if (ws->dual_func->num_rays > 0){
		clock_t start, end;
		double cpu_time_used;
		start = clock();
		for (int t = 0; t < ws->dual_func->num_terms; t++){

			disj = ws->dual_func->disj + t;
			change_lbub_from_disj(lb, ub, ws->n, disj);
		
			si->loadProblem(ws->n, ws->m,
								mip->matbeg, mip->matind,
								mip->matval, lb,
								ub, mip->obj,
								mip->sense, rhs,
								mip->rngval);
			
			si->initialSolve();

			if (si->isProvenPrimalInfeasible()){
				if (is_term_feas[t]){
					// This case could happen when the leaf is OVER_UB_PRUNED
					// and CLP stopped before proving dual unboundedness. 
					// This is OK, as long as we have the appropriate dual solution

					// printf("Warning: FarkasProof not working!\n");
					// is_term_feas[t] = FALSE;
					// si->writeLp("proof_not_working", "lp");
					// si->writeMps("proof_not_working");
					// print_dual_function(ws);
					// assert(!is_term_feas[t]);
				}
			} else if (si->isProvenOptimal()){
				if (!is_term_feas[t]){
					printf("Warning: FarkasProof not working!\n");
					// If this happens, it is problematic!!
					// Farkas Proof is giving a false positive

					// is_term_feas[t] = TRUE;
					// si->writeLp("proof_not_working", "lp");
					// si->writeMps("proof_not_working");
					// print_dual_function(ws);
					// assert(is_term_feas[t]);
				}
			}

			reset_lbub_from_disj(lb, ub, mip->lb, mip->ub, ws->n, disj);
		}
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		ws->dual_func->lp_cpu_time += cpu_time_used;
	}

#endif
	
	// Now adjust the lb/ub * dj from the disjunction terms
	for (int t = 0; t < ws->dual_func->num_terms; t++){
		
		if (!is_term_feas[t])
			continue;

		disj = ws->dual_func->disj + t;
		local_best_bound = -SYM_INFINITY;

		// Sanity check on dual optimal value
		disj = ws->dual_func->disj + t;

		change_lbub_from_disj(lb, ub, ws->n, disj);

#ifdef DEBUG_DUAL_FUNC
		si->loadProblem(ws->n, ws->m,
							mip->matbeg, mip->matind,
							mip->matval, lb,
							ub, mip->obj,
							mip->sense, rhs,
							mip->rngval);
		
		si->initialSolve();
		if (si->isProvenOptimal()){
			sanity_objVal = si->getObjValue();
		} else {
			// Should never happen!
			sanity_objVal = SYM_INFINITY;
		}
#endif
		for (i = 0; i < ws->dual_func->num_pieces; i++){
			// if (i == 173 && t == 22){
			// 	printf("here!");
			// }

			rhs_times_pi_plus_dj = rhs_times_pi[i];
			is_infty = 0;
			j = dj_start[i];
			last = duals->getVectorLast(i);

			for (j = dj_start[i]; j < last; j++)
			{
				idx = indices[j] - ws->m;
				if (elem[j] >= zerotol) // dj > 0 ?
				{
					if (lb[idx] > -SYM_INFINITY){
						rhs_times_pi_plus_dj += elem[j] * lb[idx];
					} else {
						is_infty -= elem[j];
					}
				}
				else if (elem[j] <= -zerotol) // dj < 0 ?
				{
					if (ub[idx] < SYM_INFINITY){
						rhs_times_pi_plus_dj += elem[j] * ub[idx];
					} else {
						is_infty -= elem[j];	
					}
				}
			}
			
			if ((is_infty == 0) &&
				(rhs_times_pi_plus_dj > local_best_bound)){
				local_best_bound = rhs_times_pi_plus_dj;
			} 
			// if (t == 146){
			// 	printf("local_best_bound: %.5f\n", local_best_bound);
			// }
		}

#ifdef DEBUG_DUAL_FUNC
		// Equality with sanity_objVal may not hold since
		// CLP might have stopped the solution of LPs prematurely 
		// due to the curr UB in the Tree but this is not a problem.
		// Hence we check that local_best_bound 
		// is a valid dual bound to sanity_objVal
		assert((local_best_bound - sanity_objVal) < 0.001);
#endif	
		
		if (fabs(local_best_bound - global_best_bound) > 1e-5 && // Not equals
		    local_best_bound < global_best_bound){ 
			global_best_bound = local_best_bound;
		}

		// if (t <= 146)
		// 	printf("t = %d, global_best_bound: %.5f\n", t, global_best_bound);

		reset_lbub_from_disj(lb, ub, mip->lb, mip->ub, ws->n, disj);
	}
	
TERM_EVAL_DUAL_FUNC:

	// Assuming minimization here ???
	// TODO: global_best_bound may be fractional
	// if (granularity > 1e-5 && // granularity is 1.0 or 2.0, ...
	//     fabs(floor(global_best_bound + 0.5) - global_best_bound) > 1e-5){
	// 	global_best_bound = round(global_best_bound / granularity) * granularity;
	// }
	*dual_bound = global_best_bound;

	FREE(rhs_times_pi);
	FREE(dj_start);
	FREE(lb);
	FREE(ub);
	FREE(rhs);
	if (is_term_feas);	
		FREE(is_term_feas);


	return (FUNCTION_TERMINATED_NORMALLY);
#else
	printf("evaluate_dual_func():\n");
	printf("Sensitivity analysis features are not enabled.\n");
	printf("Please rebuild SYMPHONY with these features enabled\n");
	return (FUNCTION_TERMINATED_ABNORMALLY);

#endif
}

