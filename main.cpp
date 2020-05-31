#include <stdio.h>
#include <iostream>

#include "dft.h"


using namespace std; 


int main
(
    int argc,
    char **argv
)
{
   
  int num_trees;
  int num_nodes;
  int node_ID;
  int right_child_ID;
  int left_child_ID;
  int div_time;
  int laten_value;
  int is_root;
  int latent_value;

  dft_state st;

  FILE *input = fopen("tree.txt", "r");
  
  if(input == 0) {
      cout << " can not open the input file" << endl;
      exit(1);
  }
  
  fscanf(input, "%d %d", &num_trees, &num_nodes);
  printf("The number of trees = %d %d\n",num_trees);
  
  for(int dt=1; dt <= num_trees; dt++) {
      
      // allocating space for divt and nodes array
      st[dt].divt = (double*) chk_alloc(num_nodes, sizeof(double));
      st[dt].nodes = (dft_tree_node*) chk_alloc(num_nodes, sizeof(dft_tree_node));
      
      // taking input for the num_nodes
      for(int num_node = 0 ;num_node <=num_nodes; num_node++) {
        fscanf(input, "%d %d %d %d %d", &node_ID, &is_root, &right_child_ID, &left_child_ID, &div_time, &latent_value); // latent value not needed right now
        
        st[dt].divt[node_ID] = div_time;
        st[dt].nodes->child[0] = right_child_ID;
        st[dt].nodes->child[0] = left_child_ID;

        if (is_root) {
            st[dt].root = node_ID;
        }
      }
      
    initialize_tables(num_nodes);
    dft_log_prob_div(1, dt, st, 0, 0);
  }
  

  fclose(input);
}