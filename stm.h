#ifndef STM_H_
#define STM_H_

#include "graph.h"
#include <ilcplex/ilocplex.h>
#include <vector>
using namespace std;

void stm_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true);
void hlms_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true, bool improve_heuristic = false, string heuristic_file = "");
void mst_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true);
void loose_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true);
#endif /* STM_H_ */
