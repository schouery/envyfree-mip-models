#include "utils.h"
#include <ilcplex/ilocplex.h>
#include <vector>
#include <fstream>
ILOSTLBEGIN

using namespace std;

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

void mst_create_vars(graph g, IloModel model, IloNumVarArray v, int ***columns_for_v) {
  IloEnv env = model.getEnv();
  int count = 0;
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      for(int el = 0; el < g->dbidder[i]; el++) {
        int k = g->b_adj[i][el];
        columns_for_v[i][j][k] = count;
        string s = "v_" + itos(i) + "," + itos(j) + "," + itos(k);
        v.add(IloNumVar(env, 0.0, boundp(j), ILOFLOAT, s.c_str()));
        count++;  
      }
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

IloRangeArray hlms_c1(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns) {
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

IloRangeArray mst_c1(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, IloNumVarArray v, int **columns, int ***columns_for_v) {
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      IloNumExpr exp(model.getEnv());
      exp = pi[j];
      for(int el = 0; el < g->dbidder[i]; el++) {
        int k = g->b_adj[i][el];
        exp -= v[columns_for_v[i][j][k]];
      }
      c.add(exp >= 0);
    }
  }
  return c;
}

IloRangeArray mst_c2(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, IloNumVarArray v, int **columns, int ***columns_for_v) {
  IloRangeArray c(model.getEnv());
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      for(int el = 0; el < g->dbidder[i]; el++) {
        int k = g->b_adj[i][el];
        IloNumExpr e(model.getEnv());
        e = v[columns_for_v[i][j][k]] - (g->adj[i][j] - g->adj[i][k])*theta[columns[i][k]] - p[columns[i][k]];
        c.add(e >= 0);
      }      
    }
  }
  return c;
}

IloRangeArray stm_constraints(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int M, int **columns, bool hlms=false, bool restrict_prices = true) {
  IloRangeArray c(model.getEnv());
  if(hlms)
    c.add(hlms_c1(g,model,theta,p,pi,M,columns));
  else
    c.add(stm_c1(g,model,theta,p,pi,M,columns));
  c.add(stm_c2(g,model,theta,p,pi,M,columns));
  c.add(stm_c3(g,model,theta,p,pi,M,columns));
  if(restrict_prices)
    c.add(stm_c4(g,model,theta,p,pi,M,columns));
  c.add(stm_c5(g,model,theta,p,pi,M,columns));
  return c;
}

IloRangeArray mst_constraints(graph g, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, IloNumVarArray v, int M, int **columns, int ***columns_for_v) {
  IloRangeArray c(model.getEnv());
  c.add(mst_c1(g, model, theta, p, pi, v, columns, columns_for_v));
  c.add(mst_c2(g, model, theta, p, pi, v, columns, columns_for_v));
  c.add(stm_c2(g,model,theta,p,pi,M,columns));
  c.add(stm_c3(g,model,theta,p,pi,M,columns));
  c.add(stm_c4(g,model,theta,p,pi,M,columns));
  c.add(stm_c5(g,model,theta,p,pi,M,columns));
  return c;
}

void stm_load(graph g, IloCplex cplex, IloModel model, IloNumVarArray theta, IloNumVarArray p, IloNumVarArray pi, int **columns, vector<int>& allocation, vector<double>& pricing) {
  IloNumVarArray startVar(model.getEnv());
  IloNumArray startVal(model.getEnv());
  for(int j = 0; j < g->items; j++) {
    if(boundp(j) > 0){
      startVar.add(pi[j]);
      startVal.add(pricing[j]);
    }
  }
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      startVar.add(theta[columns[i][j]]);
      startVal.add(allocation[i] == j ? 1 : 0);
    }
  }
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      if(g->adj[i][j] > 0){
        startVar.add(p[columns[i][j]]);
        startVal.add(allocation[i] != -1 ? pricing[allocation[i]] : 0);
      }
    }
  }
  cplex.addMIPStart(startVar, startVal);
}

void stm_family_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer, bool mst, bool hlms, bool loose, bool improve_heuristic = false, string heuristic_file = "") {
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
    IloNumVarArray v(env);
    int ***columns_for_v;
    stm_create_vars(g, model, theta, p, pi, columns);
    stm_create_objective(g, model, p, columns);
    if(mst) {
      columns_for_v = (int ***)malloc(g->bidders * sizeof(int **));
      for (int i = 0; i < g->bidders; ++i) {
        columns_for_v[i] = (int **)malloc(g->items * sizeof(int *));
        for (int j = 0; j < g->items; ++j)
          columns_for_v[i][j] = (int *)malloc(g->items * sizeof(int));
      }
      mst_create_vars(g, model, v, columns_for_v);  
      model.add(mst_constraints(g, model, theta, p, pi, v, M, columns, columns_for_v));
    } else {
      model.add(stm_constraints(g, model, theta, p, pi, M, columns, hlms, loose));  
    }
    //  {
    if (improve_heuristic) {
      config_cplex(cplex, false);
      ifstream file("cplex-params", ios::in);
      if(file) {
        file.close();
        cplex.readParam("cplex-params");
      }
    } else {
      config_cplex(cplex);  
    }
    // cplex.exportModel("model.lp");
    if(!integer) {
      model.add(IloConversion(env, theta, ILOFLOAT));
    } else {
      stm_load(g, cplex, model, theta, p, pi, columns, allocation, pricing);
    }
    clock_start();
    if (!cplex.solve()) {
      failed_print(g);
    } else if(integer) {        
        if(improve_heuristic) {
          ofstream file(heuristic_file.c_str(), ios::out);
          for(int i = 0; i < g->bidders; i++) {
            bool found = false;
            for(int e = 0; e < g->dbidder[i]; e++) {
              int j = g->b_adj[i][e];
              if(cplex.getValue(theta[columns[i][j]]) > 0) {
                file << i << " " << j << " " << cplex.getValue(pi[j]) << endl;
                found = true;
              }
            }
            if(!found) file << i << " -1 0.0" << endl;
          }
          file.close();
        } else {
          solution_print(cplex, env, g);  
        }
    } else
      relax_print(cplex, env);
  }
  catch (IloException& e) {
    cerr << "Concert exception caught: " << e << endl;
  }
  catch (...) {
    cerr << "Unknown exception caught" << endl;
  }
}

void hlms_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer, bool improve_heuristic, string heuristic_file) {
  stm_family_solve(g, allocation, pricing, integer, false, true, false, improve_heuristic, heuristic_file);
}
void mst_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer) {
  stm_family_solve(g, allocation, pricing, integer, true, false, false); 
}
void loose_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer) {
  stm_family_solve(g, allocation, pricing, integer, false, false, true); 
}

void stm_solve(graph g, vector<int>& allocation, vector<double>& pricing, bool integer) {
  stm_family_solve(g, allocation, pricing, integer, false, false, false); 
}
