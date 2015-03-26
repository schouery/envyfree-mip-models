#ifndef UTILS_H_
#define UTILS_H_
// #define DEBUG(x) cout << x  
#define EPSILON 0.0001

#ifndef THREADS
#define THREADS 1
#endif
#ifndef WORKMEM
#define WORKMEM 1024 // 4 GB
#endif

#define DEFAULT_FORMULATION 3

#include <string>
#include "graph.h"

enum Verbose{MinimalOutput, UserOutput, CplexOutput};

using namespace std;

string itos(int number);
string ftos(double number);
int bigM(graph g);

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

void config_cplex(IloCplex cple);
void setTimeLimit(int timelimit);
void setNodeLimit(int nodelimit);
void setVerbosity(Verbose verbosity);
Verbose getVerbosity();
void solution_print(IloCplex cplex, IloEnv env, graph g);
void relax_print(IloCplex cplex, IloEnv env);
double clock_time();

double lower_power(double value, double power);

double upper_power(double value, double power);

void failed_print(graph g);

double boundp(int j);

double boundu(int i);

void init_bounds(graph g, double epsilon);

void clock_start();
  
double clock_current_time();

void define_log(int argc, char **argv);

void save_log();



#endif /* UTILS_H_ */
