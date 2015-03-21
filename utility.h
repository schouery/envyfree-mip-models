#ifndef UTILITY_H_
#define UTILITY_H_

#include "graph.h"
#include <ilcplex/ilocplex.h>
#include <vector>
using namespace std;

void utility_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true, bool use_presolve = true);
void assignment_vars(graph g, IloModel model, IloNumVarArray x, int **columns);
IloRangeArray assignment_ineq(graph g, IloModel model, IloNumVarArray x, int **columns);

#endif /* UTILITY_H_ */
