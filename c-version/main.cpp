#include <stdio.h>
#include <iostream>

#include "dft.h"
#include "util/misc.h"


using namespace std; 


int main
(
    int argc,
    char **argv
)
{
  /* place holder vars */
  int num_trees = 1;
  int num_nodes;
  int node_ID;
  int right_child_ID;
  int left_child_ID;
  int parent_id;
  double div_time;


  

  dft_state st;

  FILE *fp = fopen("tree.txt", "r");
  
  if(fp == 0) {
      cout << " can not open the input file" << endl;
      exit(1);
  }
  
  fscanf(fp, "%d",&num_nodes);
  printf("The number of nodes in the tree = %d\n",num_nodes);
  
  for(int dt=1; dt <= num_trees; dt++) {
      
      // allocating space for divt and nodes array
      st[dt].divt = (double*) malloc((num_nodes )*sizeof(double));
      st[dt].parents = (int*) malloc((num_nodes )*sizeof(int));
      st[dt].nodes = (dft_tree_node*) malloc((num_nodes)*sizeof(dft_tree_node));
      
      // taking input for the num_nodes
      for(int num_node = 0; num_node <num_nodes - 1; num_node++) {
        
        fscanf(fp, "%d %d %lf %d %d", &node_ID, &parent_id, &div_time, &right_child_ID, &left_child_ID); // , &latent_value); // latent value not needed right now
        
        //cout << node_ID << " " << is_root << " " << right_child_ID << " "<< " " << div_time << " " << endl;
        
        st[dt].divt[node_ID] = div_time;
        st[dt].parents[node_ID] = parent_id;
        st[dt].nodes[node_ID].child[0] = right_child_ID;
        st[dt].nodes[node_ID].child[1] = left_child_ID;
      }
      


      initialize_tables(num_nodes);
      cout << "Done initialization of the table" << endl;

      int sz = dft_conv_tree(st, dt, 0, num_nodes);
      cout << " Done with converstion " << sz << endl;

      double likelihood = dft_log_prob_div(1, dt, st, 0, 0, num_nodes);
      cout << "Likeihood is = " << likelihood << endl;
  }
  
  //fclose(fp);
  return 0;
}