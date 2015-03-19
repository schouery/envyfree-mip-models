#include "utils.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

void stm_create_vars(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int **columns) {
  int k = 0;
  IloEnv env = model.getEnv();
  string s;
	for(int j = 0; j < g->items; j++) {
    s = "pi_" + itos(j);
    pi.add(IloNumVar(env, 0.0, boundp(j), ILOFLOAT, s.c_str()));
	}
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      columns[i][j] = k;
      s = "theta_" + itos(i) + "," + itos(j);
      theta.add(IloNumVar(env, 0.0, 1.0, ILOBOOL, s.c_str()));
			s = "p_" + itos(i) + "," + itos(j);
			// p.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, s.c_str()));
      p.add(IloNumVar(env, 0.0, g->adj[i][j], ILOFLOAT, s.c_str()));
      k++;
    }
  }	
}

void stm_create_objective(graph g, IloModel model, IloNumVarArray p, int **columns) {
  IloNumExpr e(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      e += p[columns[i][j]];
		}
	}
  model.add(IloMaximize(model.getEnv(), e));
}

IloRangeArray stm_stronger_c1(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns) {
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int k = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      for(int el = 0; el < g->dbidder[i]; el++) {
        int j = g->b_adj[i][el];
        e += g->adj[i][j]*theta[columns[i][j]] - p[columns[i][j]];
      }
      e += pi[k];
      c.add(e >= g->adj[i][k]);
    }
  }
  return c;
}

IloRangeArray stm_c1(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns){
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int k = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      for(int el = 0; el < g->dbidder[i]; el++) {
        int j = g->b_adj[i][el];
        if(j != k) {
          e += g->adj[i][j]*theta[columns[i][j]] - p[columns[i][j]];
          e -= g->adj[i][k]*theta[columns[i][j]];
        }
      }
      e += pi[k];
      c.add(e >= 0);
    }
  }
  return c;  
}

IloRangeArray stm_c2(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns){
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      e += g->adj[i][j]*theta[columns[i][j]] - p[columns[i][j]];
      c.add(e >= 0);
    }
  }  
  return c;
}

IloRangeArray stm_c3(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns){
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    IloNumExpr e(model.getEnv());
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      e += theta[columns[i][j]];
    }
    c.add(e <= 1);
  }
  return c;
}

IloRangeArray stm_c4(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns) {
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      e += p[columns[i][j]] - pi[j];
      c.add(e <= 0);
    }
  }
  return c;
}

IloRangeArray stm_c5(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns) {
  IloRangeArray c(model.getEnv());
  double * Rj = (double *)calloc(g->items, sizeof(double));
  for(int j = 0; j < g->items; j++) {
    for(int i = 0; i < g->bidders; i++) {
      if(Rj[j] < g->adj[i][j])
        Rj[j] = g->adj[i][j];
    }
  }
  for(int i = 0; i < g->bidders; i++) {
    for(int edge = 0; edge < g->dbidder[i]; edge++) {
      int j = g->b_adj[i][edge];
      IloNumExpr e(model.getEnv());
      e += p[columns[i][j]] - pi[j] - Rj[j]*theta[columns[i][j]];
      c.add(e >= -Rj[j]);
    }
  }
  return c;
}

IloRangeArray stm_constraints(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns, bool stronger=false, bool restrict_prices = true) {
  IloRangeArray c(model.getEnv());
  if(stronger)
    c.add(stm_stronger_c1(g,model,theta,p,pi,M,columns));
  else
    c.add(stm_c1(g,model,theta,p,pi,M,columns));
  c.add(stm_c2(g,model,theta,p,pi,M,columns));
  c.add(stm_c3(g,model,theta,p,pi,M,columns));
  if(restrict_prices)
    c.add(stm_c4(g,model,theta,p,pi,M,columns));
  c.add(stm_c5(g,model,theta,p,pi,M,columns));
  return c;
}

void stm_solve(graph g, bool integer, bool use_presolve, bool stronger, bool restrict_prices) {
  int ** columns = (int **)malloc(g->bidders * sizeof(int *));
  for(int i = 0; i < g->bidders; i++)
    columns[i] = (int *)calloc(g->items, sizeof(int));
  int M = bigM(g);
  IloEnv env;
  try {
    if(getVerbosity() != CplexOutput)
      env.setOut(env.getNullStream());
    IloModel model(env);
		IloNumVarArray theta(env);
		IloNumVarArray p(env);
		IloNumVarArray pi(env);
    IloCplex cplex(model);
    stm_create_vars(g, model, theta, p, pi, columns);
    stm_create_objective(g, model, p, columns);
    model.add(stm_constraints(g, model, theta, p, pi, M, columns, stronger, restrict_prices));
    config_cplex(cplex);
    // cplex.exportModel("model_stm.lp");
    if(!integer) {
      model.add(IloConversion(env, theta, ILOFLOAT));
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
