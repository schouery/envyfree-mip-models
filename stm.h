#ifndef STM_H_
#define STM_H_

#include "graph.h"
#include <ilcplex/ilocplex.h>

void stm_solve(graph g, bool integer = true, bool use_presolve = true, bool stronger=false, bool restrict_prices=true);

#endif /* STM_H_ */
