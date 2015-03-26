#include "utils.h"
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

int util_timelimit = 0;
int util_nodelimit = 0;
Verbose util_verbosity = MinimalOutput;
string log_filename;
ofstream output;

string itos(int number) {
  stringstream out;
  out << number;
  return out.str();
}

string ftos(double number) {
  stringstream out;
  out << number;
  return out.str();
}

int bigM(graph g) {
  int max = 0;
  for(int i = 0; i < g->bidders; i++) {
    for(int c = 0; c < g->items; c++) {
      if(max < g->adj[i][c])
        max = g->adj[i][c];
    }
  }  
  return max + 1;
}

void config_cplex(IloCplex cplex, bool log) {
  if(log && (util_verbosity == MinimalOutput || util_verbosity == UserOutput)) {
    output.open(log_filename.c_str());
    cplex.setOut(output);
  } 
  cplex.setParam(IloCplex::WorkMem, WORKMEM);
  cplex.setParam(IloCplex::Threads, THREADS);
  if(util_timelimit > 0)
    cplex.setParam(IloCplex::TiLim, util_timelimit);
  if(util_nodelimit)
    cplex.setParam(IloCplex::NodeLim, util_nodelimit);
}

void setTimeLimit(int timelimit) {
  util_timelimit = timelimit;
}

void setNodeLimit(int nodelimit) {
  util_nodelimit = nodelimit;
}

void setVerbosity(Verbose verbosity) {
  util_verbosity = verbosity;
}

Verbose getVerbosity() {
  return util_verbosity;
}

inline double gap(double integer, double relaxed) {
  return ((relaxed - integer)/integer);
}

void failed_print(graph g) {
  cout << "  bidders: " << g->bidders << endl;
  cout << "  items: " << g->items << endl;
  cout << "  time: " << clock_time() << endl;
  cout << "  status: " << "Failed" << endl;
}

double clock_time() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1000000.0;
}

double clock_start_time;

void clock_start() {
  clock_start_time  = clock_time();
}

double clock_current_time() {
  return clock_time() - clock_start_time;
}

double lower_power(double value, double power) {
  if(power == 1.0 || value == 0.0)
    return value;
  double p;
  for(p = 1.0; p*power <= value; p *= power);
  return p;
}

double upper_power(double value, double power) {
  if(power == 1.0 || value == 0.0)
    return value;
  double p;
  for(p = 1.0; p < value; p *= power);
  return p;
}

double *bp, *bu, lp;

double boundp(int j){
  return bp[j];
}

double boundu(int i) {
  return bu[i];
}

void init_bounds(graph g, double epsilon) {
  bp = (double *)calloc(g->items, sizeof(double));
  bu = (double *)calloc(g->bidders, sizeof(double));
  for(int i = 0; i < g->bidders; i++) {
    for(int e = 0; e < g->dbidder[i]; e++) {
      int j = g->b_adj[i][e];
      if(bp[j] < g->adj[i][j])
        bp[j] = g->adj[i][j];
      if(bu[i] < g->adj[i][j])
        bu[i] = g->adj[i][j];
    }
  }
}

void solution_print(IloCplex cplex, IloEnv env, graph g) {
  cout << "  bidders: " << g->bidders << endl;
  cout << "  items: " << g->items << endl;
  cout << "  edges: " << g->edges << endl;
  cout << "  time: " << clock_current_time() << endl;
  cout << "  status: " << cplex.getStatus() << endl;
  cout << "  value: " << cplex.getObjValue() << endl;
  cout << "  upper_bound: " << cplex.getBestObjValue() << endl;
  cout << "  gap: " << gap(cplex.getObjValue(), cplex.getBestObjValue()) << endl;
  cout << "  binaries: " << cplex.getNbinVars() << endl;
  cout << "  columns: " << cplex.getNcols() << endl;
  cout << "  iterations: " << cplex.getNiterations() << endl;
  cout << "  nodes: " << cplex.getNnodes() << endl;
  cout << "  nodes_left: " << cplex.getNnodesLeft() << endl;
  cout << "  rows: " << cplex.getNrows() << endl;
}

void relax_print(IloCplex cplex, IloEnv env) {
  cout << "  relax_time: "  << clock_current_time() << endl;
  cout << "  relax_value: " << cplex.getObjValue() << endl;
}

void define_log(int argc, char **argv) {
  log_filename = argv[argc-1];
  log_filename += "_opts";
  for(int i = 1; i < argc-1; i++) {
    log_filename += argv[i];
  }
  log_filename += ".log";
}

void save_log() {
  output.close();
}
