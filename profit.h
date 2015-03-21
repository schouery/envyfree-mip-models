#ifndef PROFIT_H_
#define PROFIT_H_

#include "graph.h"
#include <vector>
using namespace std;

void profit_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer = true, bool use_presolve = true);

#endif /* PROFIT_H_ */
