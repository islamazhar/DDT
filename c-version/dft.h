/*
    source: https://github.com/mrquincle/fbm/
*/

#include "util/misc.h"
#include <math.h>
#include <iostream>

using namespace std;

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
  int *parents;
  double *divt;		/* Divergence times, offset for indexing from 1 */
  dft_tree_node *nodes;	/* Nodes in tree, indexed from 1 */

} dft_state[Max_trees];


double dft_cdiv 
( double c1,
  int dt,			/* Index of tree (from 0) */
  double t			/* Time argument, non-negative, <= 1 */
)
{ 
    if (t>1) abort();
    return c1 * log(1-t);
}

/* TABLES OF PRECOMPUTED VALUES. */

double *log_factorial = 0;	/* Logs of factorial from 0 to N_train */
double *sum_reciprocals = 0;	/* Table of sum from i=1 to n of 1/i, for
				   n from 0 to N_train */

/* INITIALIZE THE TABLES OF PRECOMPUTED VALUES.*/
static void initialize_tables (int N_train)
{
  int n;

  if (log_factorial==0) 
  { 
    log_factorial = (double*) malloc ((N_train+1)*sizeof(double)); // chk_alloc is defined in misc.c

    log_factorial[0] = 0;

    for (n = 1; n<=N_train; n++)
    { log_factorial[n] = log_factorial[n-1] + log((double)n);
    }
  }

  if (sum_reciprocals==0)
  { 
    sum_reciprocals = (double*) malloc ((N_train+1)*sizeof(double));

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
  // cout << t << " = " << c1  << endl;
  return c1 / (1-t);
}


double dft_log_prob_div
( double c1,
  int dt,			    /* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  int root,			    /* Root of sub-tree to look at, zero for all */
  double div0,			/* Divergence time of branch to root */
  int n_cases
)
{
  double lp=0, div1;
  int n;

  if (root == 0) { // find the root of the current tree. The root of the current tree has parent = 0
        for(int i=1;i<n_cases;i++) {
            if(st[dt].parents[i] == 0) {
                root = i;
                break;
            }
        }
  }

  if (root == -1 ) return 0; // root == -1 means --> terminal node

  n = totpts(st[dt].nodes[root]);
  cout << n << " " << root << endl;
  if (n<2) abort(); // unnecessary

  div1 = st[dt].divt[root];
  if (div1<div0) abort();
   
  div1 = st[dt].divt[root];
  if (div1<div0) abort();
   

  if (div1>div0)
  { lp = - sum_reciprocals[n-1] * (dft_cdiv(1,dt,div1)-dft_cdiv(1,dt,div0));
  }
  

  lp += log (dft_div (c1, dt, div1));

 // cout << chld(st[dt].nodes[cur_node],0) << " " << chld(st[dt].nodes[cur_node],1) << endl;

  lp += dft_log_prob_div (c1, dt, st, chld(st[dt].nodes[root],0), div1, n_cases);
  lp += dft_log_prob_div (c1, dt, st, chld(st[dt].nodes[root],1), div1, n_cases);

  return lp;
}
int dft_conv_tree
(
        dft_state dft,		/* Array of nodes to store all dft trees*/
        int dt,             /* The id of tree */
        int cur_node,		/* the id of the current node in the tree*/
        int n_cases         /* Number of cases in tree */
)
{
    if (cur_node == 0) { // find the root of the current tree. The root of the current tree has parent = 0
        for(int i=1;i<n_cases;i++) {
            if(dft[dt].parents[i] == 0) {
                cur_node = i;
                break;
            }
        }
    }

    if (cur_node == -1) return 1; // leaf node

    int l = chld(dft[dt].nodes[cur_node],0);
    int r = chld(dft[dt].nodes[cur_node],1);

    dft[dt].nodes[cur_node].n_points[0] = dft_conv_tree(dft, dt, l, n_cases);
    dft[dt].nodes[cur_node].n_points[1] = dft_conv_tree(dft, dt, r, n_cases);

    return totpts(dft[dt].nodes[cur_node]);
}



