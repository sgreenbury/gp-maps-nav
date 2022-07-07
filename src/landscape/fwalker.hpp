#ifndef _fwalker_h_included_
#define _fwalker_h_included_

#include <math.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#include "base/utilities.hpp"

using namespace std;

class FW {
  std::set<int> u, v, p;
  std::vector<int> positions_;
  std::vector<int> genos;
  std::vector<double> fitness;
  std::vector<std::vector<int>> phenotype_transitions;

  int timing;
  clock_t start, middle, end;
  std::ofstream times;

  int verbose_;

  int K;
  int L;
  int D;
  int sparse;

  int neutral_mutations_;

  int run;
  int swaps;
  int all;
  int N_p;
  // int phenotype;
  int nd_ref;
  double max_fit;
  int max_fit_p;
  int threshold;

  int targets;
  int target;
  int sources;
  int source;

  int print_component;

  int *sequence;

  int uphill_count;
  int downhill_count;

  int stop_at_fittest_;
  int record_phenotype_transitions_;
  int reseed_fitness_;

  int remap_phenotypes_;

  std::string fitness_assignment;

  std::string file_name_;

  std::ofstream stats;
  std::ofstream fp1;

  std::string geno_fname_;
  std::string outpath_;

  void Init();
  void Clear();

  void AssignFitness();
  void AssignFitness(int max);
  void SetFitness();
  void SingleRun(int geno);

  void Transfer(std::set<int> &from, std::set<int> &to, int x);
  void Transfer(std::unordered_set<int> &from, std::unordered_set<int> &to,
                int x);
  void fStartComponent(int p, int n);
  void fAddEdge(int x, int y);
  void fEndComponent();

 public:
  FW(int BASE, int GENOME_LENGTH, int ND, int TARGETS, int SOURCES, int ALL,
     int SWAPS, int DIMENSION, int NEUTRAL_MUTATIONS, int RUN,
     int THRESHOLD = 10000000, int STOP_AT_FITTEST = 1,
     int RECORD_PHENOTYPE_TRANSITIONS = 0, int RESEED_FITNESS = 0,
     int REMAP_PHENOTYPES = 0, std::string FITNESS_ASSIGNMENT = "target",
     std::string GENO_FNAME = "geno_list0.txt", std::string file_name = "",
     std::string OUTPATH = "", int VERBOSE = 0);

  ~FW();
  void Run();
};

#endif
