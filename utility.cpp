#include "utility.h"
#include "utils.h"
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>
ILOSTLBEGIN

void assignment_vars(graph g, IloModel model, IloNumVarArray x, int **columns) {
  int k = 0;
  IloEnv env = model.getEnv();
  string s;
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      columns[i][j] = k;
      s = "x_" + itos(i) + "," + itos(j);
      x.add(IloNumVar(env, 0.0, 1.0, ILOBOOL, s.c_str()));
      k++;
    }
  }
}

IloRangeArray assignment_ineq(graph g, IloModel model, IloNumVarArray x, int **columns) {
  string s;
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    IloNumExpr e(model.getEnv());
    for(int k = 0; k < g->dbidder[i]; k++) {
      int j = g->b_adj[i][k];
      e += x[columns[i][j]];
    }
    c.add(e <= 1);
    s = "a_" + itos(i);
    c[c.getSize()-1].setName(s.c_str());
  }
  return c;
}

IloRangeArray utility_u_ineq(graph g, IloModel model, IloNumVarArray x, IloNumVarArray p, IloNumVarArray u, int M, int **columns) {
  string s;
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    // ni = no item
    // IloNumExpr ni(model.getEnv());
    // ni += u[i];
    for(int k = 0; k < g->dbidder[i]; k++) {
      int j = g->b_adj[i][k];
      // pu = positive utility
      IloNumExpr pu(model.getEnv());
      pu += u[i];
      pu += p[j];
      c.add(pu >= g->adj[i][j]);
      s = "e_" + itos(i) + "," + itos(j);
      c[c.getSize()-1].setName(s.c_str());
      // ru = right utility
      // IloNumExpr ru(model.getEnv());
      // ru += u[i];
      // ru += ((boundp(j) + boundu(i)) - g->adj[i][j])*x[columns[i][j]];
      // ru += p[j];
      // c.add(ru <= (boundp(j) + boundu(i)));
      // s = "d_" + itos(i) + "," + itos(j);
      // c[c.getSize()-1].setName(s.c_str());
      IloNumExpr ru(model.getEnv());
      ru += u[i];
      ru += p[j];
      for(int e = 0; e < g->dbidder[i]; e++) {
        int jl = g->b_adj[i][e];
        ru -= g->adj[i][jl]*x[columns[i][jl]];
      }
      ru += boundp(j)*x[columns[i][j]];
      c.add(ru <= boundp(j));
      s = "d_" + itos(i) + "," + itos(j);
      c[c.getSize()-1].setName(s.c_str());
      // ni
      // ni -= boundu(i)*x[columns[i][j]];
    }
    // c.add(ni <= 0);c
  }
  return c;
}

void utility_create_vars(graph g, IloModel model, IloNumVarArray x, IloNumVarArray p, IloNumVarArray u, int **columns) {
  IloEnv env = model.getEnv();
  string s;
  assignment_vars(g, model, x, columns);  
  for(int j = 0; j < g->items; j++) {
    s = "p_" + itos(j);
    p.add(IloNumVar(env, 0.0, boundp(j), ILOFLOAT, s.c_str()));
    // p.add(IloNumVar(env, lowerp(j), boundp(j), ILOFLOAT, s.c_str()));
  }
  for(int i = 0; i < g->bidders; i++) {
    s = "u_" + itos(i);
    u.add(IloNumVar(env, 0.0, boundu(i), ILOFLOAT, s.c_str()));
  }
}

// We fix the price of an item in 0 if nobody want it,
// This is to prevent problems with Cplex
IloRangeArray utility_fix(graph g, IloModel model, IloNumVarArray p) {
  IloRangeArray c(model.getEnv());
  string s;
  for(int j = 0; j < g->items; j++) {
    if(g->ditem[j] == 0) {
      IloNumExpr e(model.getEnv());
      e += p[j];
      c.add(e <= 0);
      s = "fi_" + itos(j);
      c[c.getSize()-1].setName(s.c_str());
    }
  }
  return c;
}

void utility_create_objective(graph g, IloModel model, IloNumVarArray x, IloNumVarArray u, int **columns) {
  IloNumExpr e(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int k = 0; k < g->dbidder[i]; k++) {
      int j = g->b_adj[i][k];
      e += g->adj[i][j] * x[columns[i][j]];
    }
    e -= u[i];
  }
  model.add(IloMaximize(model.getEnv(), e));
}

void utility_solution_print(graph g, IloEnv env, IloCplex cplex, IloNumVarArray x, IloNumVarArray p, IloNumVarArray u, int **columns) {
  IloNumArray vals(env);
  env.out() << "Solution status = " << cplex.getStatus() << endl;
  env.out() << "Solution value  = " << cplex.getObjValue() << endl;
  env.out() << "Solution  = " << endl;
  for(int i = 0; i < g->bidders; i++) {
    bool found = false;
    for(int j = 0; !found && j < g->items; j++) {
      if(g->adj[i][j] && cplex.getValue(x[columns[i][j]]) > 0.5) {
        env.out() << i << " --> " << j;
        env.out() << " pj: " << cplex.getValue(p[j]);
        env.out() << " vij: " << g->adj[i][j];
        env.out() << " ui: " << cplex.getValue(u[i]);
        env.out() << " dif: " << g->adj[i][j] - cplex.getValue(p[j]);
        env.out() <<  endl;
        found = true;
      }
    }
    if(!found) {
      env.out() << i << " not matched " ;
      env.out() << " ui: " << cplex.getValue(u[i]);
      env.out() << " dif: 0";
      env.out() <<  endl;
    }
  }
}

IloRangeArray utility_constraints(graph g, IloModel model, IloNumVarArray x, IloNumVarArray p, IloNumVarArray u, int **columns) {
  IloRangeArray c(model.getEnv());
  c.add(assignment_ineq(g, model, x, columns));
  c.add(utility_u_ineq(g, model, x, p, u, 0, columns));
  // c.add(utility_fix(g, model, p));
  return c;
}

void utility_solve(graph g, bool integer, bool use_presolve) {
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
    IloNumVarArray u(env);
    IloCplex cplex(model);
    utility_create_vars(g, model, x, p, u, columns);
    utility_create_objective(g, model, x, u, columns);
    model.add(utility_constraints(g, model, x, p, u, columns));
    config_cplex(cplex); 
    if(!integer) {
      model.add(IloConversion(env, x, ILOFLOAT));
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
