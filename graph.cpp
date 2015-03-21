#include "graph.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

graph graph_new(int bidders, int items)
{
  graph g = (graph)malloc(sizeof(Graph));
  g->dbidder = (int *)calloc(bidders,sizeof(int));
  g->ditem = (int *)calloc(items,sizeof(int));
  g->bidders = bidders;
  g->items = items;
  g->edges = 0;
  g->b_adj = (int **)malloc(bidders * sizeof(int *));
  g->i_adj = (int **)malloc(items * sizeof(int *));
  // g->adj = (int **)malloc(bidders * sizeof(int *));
  g->adj = (double **)malloc(bidders * sizeof(double *));
  for(int i = 0; i < bidders; i++) {
    // g->adj[i] = (int *)calloc(items, sizeof(int));
    g->adj[i] = (double *)calloc(items, sizeof(double));
    g->b_adj[i] = (int *)calloc(items, sizeof(int));
  }
  for(int j = 0; j < items; j++)
   g->i_adj[j] = (int *)calloc(bidders, sizeof(int));
  return g;
}

// void graph_add_edge(graph g, int bidder, int item, int value)
void graph_add_edge(graph g, int bidder, int item, double value)
{
  g->b_adj[bidder][g->dbidder[bidder]] = item;
  g->i_adj[item][g->ditem[item]] = bidder;
  (g->dbidder[bidder])++;
  (g->ditem[item])++;
  g->edges = g->edges + 1;
  g->adj[bidder][item] = value;
}

void graph_print(graph g)
{
  cout << "Bidders: " << g->bidders << endl;
  cout << "Items: " << g->items << endl;
  for(int i = 0; i < g->bidders; i++)
  {
    cout << i << ": ";
    for(int j = 0; j < g->items; j++)
    {
      if(g->adj[i][j])
        cout << j << "(" << g->adj[i][j] << ") ";
    }
    cout << endl;
  }
}

void graph_save(graph g, char *filename) {
  ofstream file;
  file.open(filename);
  file << g->bidders << " " << g->items << " " << g->edges << endl;
  for(int i = 0; i < g->bidders; i++) {
    for(int j = 0; j < g->items; j++) {
      if(g->adj[i][j])
        file << i << " " << j << " " << g->adj[i][j] << endl;
    }
  }
  file.close();
}
