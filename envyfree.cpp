#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "utility.h"
#include "profit.h"
#include "utils.h"
#include "stm.h"
#define PROFIT_FORMULATION 1
#define STM_FORMULATION 2
#define UTILITY_FORMULATION 3
#define HLMS_FORMULATION 4
#define LOOSE_FORMULATION 5
#define MST_FORMULATION 6
#define IMPROVE_HEURISTIC 7

using namespace std;
char *main_filename = NULL;
double epsilon = 0.0;
bool approx = false;
bool integer_prices = false;
bool greedy = false;

graph read(char *filename) {
  ifstream file;
  file.open(filename);
  int bidders, items, edges;
  file >> bidders >> items >> edges;
  graph g = graph_new(bidders, items);  
  for(int e = 0; e < edges; e++) {
    int b, i;
    double vij;
    file >> b >> i >> vij;
    graph_add_edge(g, b, i, vij);
  }
  file.close();
  return g;
}

void usage(char *argv) {
  cout << argv << " [-U|-S|-H|-P|-L|-M|-I] [-t time] [-n nodes] [-v level] [-h] [-e epsilon] filename" << endl;
  cout << "-U: " << "Utility Formulation" << endl;
  cout << "-S: " << "STM Formulation" << endl;
  cout << "-H: " << "HLMS Formulation" << endl;
  cout << "-L: " << "Loose Formulation" << endl;
  cout << "-P: " << "Profit Formulation" << endl;
  cout << "-M: " << "MST Formulation" << endl;
  cout << "-I: " << "Improves current heuristic" << endl;
  cout << "-t time: " << "Set time as a time limit" << endl;
  cout << "-n nodes: " << "Set nodes as a node limit" << endl;
  cout << "-v level: " << "Set the verbosity level" << endl;
  cout << "\t 0 Minimal Output" << endl;
  cout << "\t 1 User Output" << endl;
  cout << "\t 2 Cplex Output" << endl;  
}

int setup(int argc, char **argv) {
  int c, formulation = DEFAULT_FORMULATION; 
  opterr = 0; 
  while ((c = getopt (argc, argv, "USHPLMIt:n:v:h")) != -1) {
    switch (c) {
      case 'U':
        formulation = UTILITY_FORMULATION;
        break;
      case 'P':
        formulation = PROFIT_FORMULATION;
        break;
      case 'S':
        formulation = STM_FORMULATION;
        break;
      case 'H':
        formulation = HLMS_FORMULATION;
        break;
      case 'L':
        formulation = LOOSE_FORMULATION;
        break;
      case 'M':
        formulation = MST_FORMULATION;
        break;
      case 'I':
        formulation = IMPROVE_HEURISTIC;
        break;
      case 't':
        setTimeLimit(atoi(optarg));
        break;
      case 'n':
        setNodeLimit(atoi(optarg));
        break;
      case 'v':
        setVerbosity((Verbose) atoi(optarg));
        break;
      case 'h':
        return -1;
      case '?':
        if (optopt == 't' || optopt == 'v')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return -1;
      default:
        return -1;
    }
  }
  for (int index = optind; index < argc; index++)
    main_filename = argv[index];
  if(!main_filename)
    return -1;
  return formulation;
}

void read_heuristic(string heuristic_file, graph g, vector<int>& allocation, vector<double>& pricing) {
  ifstream file(heuristic_file.c_str(), ios::in);
  for (int i = 0; i < g->items; ++i)
    pricing.push_back(boundp(i));
  for (int i = 0; i < g->bidders; ++i)
    allocation.push_back(-1);
  for (int i = 0; file && i < g->bidders; ++i) {
    int bidder, item;
    double value;
    file >> bidder >> item >> value;
    allocation[bidder] = item;
    if(item >= 0)
      pricing[item] = value;
  }
  file.close();
}

int main(int argc, char **argv) {
  int option = setup(argc, argv);
  define_log(argc, argv);
  graph g;
  vector<int> allocation;
  vector<double> pricing;
  string heuristic_file(main_filename);
  heuristic_file.append("-heuristic");
  if(option >= 0) {
    g = read(main_filename);
    init_bounds(g, epsilon);
    read_heuristic(heuristic_file, g, allocation, pricing);
  }
  switch(option) {
    case UTILITY_FORMULATION:
      utility_solve(g, allocation, pricing);
      // utility_solve(g, allocation, pricing, false);
      // utility_solve(g, allocation, pricing, false, false);
      break;
    case PROFIT_FORMULATION:
      profit_solve(g, allocation, pricing);
      // profit_solve(g, allocation, pricing, false);
      // profit_solve(g, allocation, pricing, false, false);
      break;
    case STM_FORMULATION:
      stm_solve(g, allocation, pricing);
      // stm_solve(g, allocation, pricing, false);
      // stm_solve(g, allocation, pricing, false, false);
      break;
    case HLMS_FORMULATION:
      hlms_solve(g, allocation, pricing);
      // stm_solve(g, allocation, pricing, false, true, true);
      // stm_solve(g, allocation, pricing, false, false, true);
      break;
    case LOOSE_FORMULATION:
      loose_solve(g, allocation, pricing);
      // stm_solve(g, allocation, pricing, false, true, true, false);
      // stm_solve(g, allocation, pricing, false, false, true, false);
      break;
    case MST_FORMULATION:
      mst_solve(g, allocation, pricing);
      break;
    case IMPROVE_HEURISTIC:
      hlms_solve(g, allocation, pricing, true, true, heuristic_file);
      break;
    default:
    usage(argv[0]);
  }
  save_log();
  return 0;
}
