#ifndef STM_H_
#define STM_H_

#include "graph.h"
#include <ilcplex/ilocplex.h>
#include <vector>
using namespace std;

void stm_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true, bool use_presolve = true, bool stronger=false, bool restrict_prices=true, bool mst=false);

#endif /* STM_H_ */
