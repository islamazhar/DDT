#include "dft.h"

double dft_log_prob_div
( dft_hypers *hyp,		/* Hyperparameters for diffusion tree model */
  int dt,			    /* Index of tree (from 0) */
  dft_state st,			/* Pointers to state */
  int root,			    /* Root of sub-tree to look at, zero for all */
  double div0			/* Divergence time of branch to root */
)
{
  double lp, div1;
  int n;

  initialize_tables();

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
  { lp = - sum_reciprocals[n-1] * (dft_cdiv(hyp,dt,div1)-dft_cdiv(hyp,dt,div0)); // dft_cdiv is defined inside dft-div.c
  }
  
  // note shuning: log prob of divergence;
  lp += log (dft_div (hyp, dt, div1));

  // lp = log prob of divergence - log prob of not divergence ()
  
  // note Shuning: recursive function for child nodes..... think in details;
  lp += dft_log_prob_div (hyp, dt, st, chld(st[dt].nodes[-root],0), div1);
  lp += dft_log_prob_div (hyp, dt, st, chld(st[dt].nodes[-root],1), div1);

  return lp;
}

double dft_cdiv 
( dft_hypers *h,		/* Hyperparameters, including divergence func.*/
  int dt,			/* Index of tree (from 0) */
  double t			/* Time argument, non-negative, <= 1 */
)
{ if (t>1) abort();
  return h->c0[dt] * t - (h->c1[dt]>0 ? h->c1[dt] * log(1-t) : 0)
                       + (h->c2[dt]>0 ? h->c2[dt] * (1/(1-t) - 1) : 0);
}