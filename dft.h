/*
    source: https://github.com/mrquincle/fbm/
*/

#include "util/misc.h"
#include <math.h>

#define Max_trees 9		/* Maximum number of diffusion trees in model */

#define chld(n,i) ((n).child[0]<(n).child[1] ? (n).child[i] \
                                             : (n).child[1-i])

#define npts(n,i) ((n).child[0]<(n).child[1] ? (n).n_points[i] \
                                             : (n).n_points[1-i])

#define totpts(n) ((n).n_points[0] + (n).n_points[1])




typedef struct
{
  int child[2];			/* The two children, -ve for non-terminals,
				   +ve for terminal nodes (ie, data points) */

  int n_points[2];		/* Number of data points ultimately reached
				   via each child. */
                   

} dft_tree_node;


typedef struct
{
  int root;		/* Parent indexes, offset for use with +/- indexes */
  double *divt;		/* Divergence times, offset for indexing from 1 */
  dft_tree_node *nodes;	/* Nodes in tree, indexed from 1 */

} dft_state[Max_trees];



/* TABLES OF PRECOMPUTED VALUES. */

static double *log_factorial;	/* Logs of factorial from 0 to N_train */
static double *sum_reciprocals;	/* Table of sum from i=1 to n of 1/i, for
				   n from 0 to N_train */

/* INITIALIZE THE TABLES OF PRECOMPUTED VALUES.*/
static void initialize_tables (int N_train)
{
  int n;

  if (log_factorial==0) 
  { 
    log_factorial = (double*) chk_alloc (N_train+1, sizeof *log_factorial); // chk_alloc is defined in misc.c

    log_factorial[0] = 0;

    for (n = 1; n<=N_train; n++)
    { log_factorial[n] = log_factorial[n-1] + log((double)n);
    }
  }

  if (sum_reciprocals==0)
  { 
    sum_reciprocals = (double*) chk_alloc (N_train+1, sizeof *sum_reciprocals);

    sum_reciprocals[0] = 0;
    for (n = 1; n<=N_train; n++)
    { sum_reciprocals[n] = sum_reciprocals[n-1] + 1.0/n;
    }
  }
}


double dft_div 
( double c1,		/* Hyperparameters, including divergence func.*/
  int dt,			/* Index of tree (from 0) */
  double t			/* Time argument, non-negative, <= 1 */
)
{ if (t>1) abort();
  if (t==1) 
  { return c1 ? 1 : 0;
  }
  return c1 / (1-t);
}


double dft_log_prob_div
( double c1,		/* Hyperparameters for diffusion tree model */
  int dt,			    /* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  int root,			    /* Root of sub-tree to look at, zero for all */
  double div0			/* Divergence time of branch to root */
)
{
  double lp, div1;
  int n;

  

  if (root==0)
  { 
      root = st[dt].root;
  }

  if (root>0) return 0; // root > 0 means --> terminal node

  n = totpts(st[dt].nodes[-root]);
  if (n<2) abort();

  div1 = st[dt].divt[-root]; // dft_state[tree_index].divt[-root] what is divt?
  if (div1<div0) abort();
   
  // note shuning: log prob of no divergence;
  if (div1>div0)
  { lp = - sum_reciprocals[n-1] * (dft_cdiv(c1,dt,div1)-dft_cdiv(c1,dt,div0)); // dft_cdiv is defined inside dft-div.c
  }
  
  // note shuning: log prob of divergence;
  lp += log (dft_div (c1, dt, div1));

  // lp = log prob of divergence - log prob of not divergence ()
  
  // note Shuning: recursive function for child nodes..... think in details;
  lp += dft_log_prob_div (c1, dt, st, chld(st[dt].nodes[-root],0), div1);
  lp += dft_log_prob_div (c1, dt, st, chld(st[dt].nodes[-root],1), div1);

  return lp;
}

double dft_cdiv 
( double c1,
  int dt,			/* Index of tree (from 0) */
  double t			/* Time argument, non-negative, <= 1 */
)
{ 
    if (t>1) abort();
    return c1 * log(1-t);
}