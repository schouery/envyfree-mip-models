#include "profit.h"
#include "utility.h"
#include "utils.h"
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>
ILOSTLBEGIN

IloRangeArray profit_constraints(graph g, IloModel model, IloNumVarArray x, IloNumVarArray p, IloNumVarArray z, int **columns) {
  string s;
  IloRangeArray c(model.getEnv());
  c.add(assignment_ineq(g, model, x, columns));
  for(int i = 0; i < g->bidders; i++) {
  	for(int edge = 0; edge < g->dbidder[i]; edge++) {
  		int k = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      for(int other_edge = 0; other_edge < g->dbidder[i]; other_edge++) {
        int j = g->b_adj[i][other_edge];
        e += g->adj[i][j]*x[columns[i][j]];
      }
      e += - z[i] + p[k];
      c.add(e >= g->adj[i][k]);
      s = "b_" + itos(i) + "," + itos(k);
      c[c.getSize()-1].setName(s.c_str());
  	}
  }
  for(int i = 0; i < g->bidders; i++) {
    IloNumExpr e(model.getEnv());
    e += z[i];
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      e -= g->adj[i][j]*x[columns[i][j]];
    }    
    c.add(e <= 0);
    s = "c_" + itos(i);
    c[c.getSize()-1].setName(s.c_str());
  }
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      e += z[i] - p[j] - boundp(j)*x[columns[i][j]];
      c.add(e >= - boundp(j));
      s = "e_" + itos(i) + "," + itos(j);
      c[c.getSize()-1].setName(s.c_str());
    }
  }
  return c;
}

void profit_create_vars(graph g, IloModel model, IloNumVarArray x, IloNumVarArray p, IloNumVarArray z, int **columns) {
  IloEnv env = model.getEnv();
  string s;
  assignment_vars(g, model, x, columns);  
  for(int j = 0; j < g->items; j++) {
    s = "p_" + itos(j);
    p.add(IloNumVar(env, 0.0, boundp(j), ILOFLOAT, s.c_str()));
    // p.add(IloNumVar(env, lowerp(j), boundp(j), ILOFLOAT, s.c_str()));
  }
  for(int i = 0; i < g->bidders; i++) {
    s = "z_" + itos(i);
    z.add(IloNumVar(env, 0.0, boundu(i), ILOFLOAT, s.c_str()));
  }
}

void profit_create_objective(graph g, IloModel model, IloNumVarArray x, IloNumVarArray z, int **columns) {
  IloNumExpr e(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    e += z[i];
  }
  model.add(IloMaximize(model.getEnv(), e));
}

void profit_load(graph g, IloCplex cplex, IloModel model, IloNumVarArray x, IloNumVarArray p, IloNumVarArray z, int **columns, vector<int>& allocation, vector<double>& pricing) {
  IloNumVarArray startVar(model.getEnv());
  IloNumArray startVal(model.getEnv());
  for(int j = 0; j < g->items; j++) {
    if(boundp(j) > 0){
      startVar.add(p[j]);
      startVal.add(pricing[j]);
    }
  }
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      startVar.add(x[columns[i][j]]);
      startVal.add(allocation[i] == j ? 1 : 0);
    }
  }
  for(int i = 0; i < g->bidders; i++) {
    if(boundu(i) > 0){
      startVar.add(z[i]);
      startVal.add(allocation[i] != -1 ? pricing[allocation[i]] : 0);
    }
  }
  cplex.addMIPStart(startVar, startVal);
}

void profit_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer, bool use_presolve) {
  int **columns = (int **)malloc(g->bidders * sizeof(int *));
  for(int i = 0; i < g->bidders; i++)
    columns[i] = (int *)calloc(g->items, sizeof(int));
  IloEnv env;
  try {
    if(getVerbosity() != CplexOutput)
      env.setOut(env.getNullStream());
    IloModel model(env);
    IloNumVarArray x(env);
    IloNumVarArray p(env);
    IloNumVarArray z(env);
    IloCplex cplex(model);
    profit_create_vars(g, model, x, p, z, columns);
    profit_create_objective(g, model, x, z, columns);
    model.add(profit_constraints(g, model, x, p, z, columns));
    config_cplex(cplex); 
    if(!integer) {
      model.add(IloConversion(env, x, ILOFLOAT));
    } else {
      profit_load(g, cplex, model, x, p, z, columns, allocation, pricing);
    }
    if(!use_presolve) {
      cplex.setParam(IloCplex::PreInd, 0);
    }
    clock_start();
    if (!cplex.solve()) {
      failed_print(g);
    } else {
      if(integer)
        solution_print(cplex, env, g);  
      else
        relax_print(cplex, env, use_presolve);
    }
  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
  }
  catch (...) {
    cerr << "Unknown exception caught" << endl;
  }
}
