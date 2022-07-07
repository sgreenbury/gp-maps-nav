#ifndef _fwalker_frna_h_included_
#define _fwalker_frna_h_included_

#include <math.h>

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "rna_fold.hpp"

// using namespace std;

double CompareStringsLevDist(const std::string &s1, const std::string &s2);

class FW_frna {
  //  std::unordered_set <int> u,v,p;
  std::unordered_set<std::string> u, v;
  std::vector<std::string> u2;
  std::vector<std::string> phenos_start;
  std::vector<std::string> genos_start;
  std::map<std::string, std::string> genos;
  std::unordered_set<std::string> p;
  std::unordered_set<std::string> pheno_set;

  std::map<std::string, double> fitness;

  std::vector<double> neighbour_fitnesses;

  std::vector<char> alphabet;

  int TIMING;
  clock_t start, middle, end;
  std::ofstream times;

  int K;
  int L;
  int sparse;

  // int seed;
  int all;
  int N_p;
  int phenotype;
  int nd_ref;
  double max_fit;
  int max_fit_p;
  int threshold;
  double fit_prob;

  int source;
  int target;

  int targets;
  int target_index;
  int sources;
  int source_index;

  int shape_level = 0;
  std::string trivial_shape;

  int population_size;

  int load_genos_;
  std::string fitness_measure;

  int print_component;

  int neutral_mutations_;

  int record_diffs_only;

  int debug;

  std::string geno_fname;
  std::string pheno_fname;

  std::string file_name_;
  std::string outpath_;

  std::ofstream stats;
  std::ofstream fp1;

  //  System Test;

  void SetAlphabet();
  void Init();
  void Clear();

  double CompareStrings(std::string s1, std::string s2);
  double SetUpFitness(std::string g_p);

  void SingleRun(std::string geno);
  void SingleRunDFS(std::string geno);

  void SingleRunMono(std::string geno);

  void Transfer(std::set<int> &from, std::set<int> &to, int x);
  void Transfer(std::unordered_set<int> &from, std::unordered_set<int> &to,
                int x);
  void Transfer(std::set<std::string> &from, std::set<std::string> &to,
                std::string x);
  void Transfer(std::unordered_set<std::string> &from,
                std::unordered_set<std::string> &to, std::string x);
  void fStartComponent(int p, int n);
  void fAddEdge(std::string x, std::string y);
  void fEndComponent();

  std::string FindGenotype(std::string phys_seq_string);

 public:
  FW_frna(int BASE, int GENOME_LENGTH, int ND, int TESTS, int SAMPLES, int ALL,
          int NEUTRAL_MUTATIONS = 1, int THRESHOLD = 2000000,
          std::string OUTPATH = "./", std::string FILE_NAME = "0",
          std::string GENO_FNAME = "cs_l20.txt", int LOAD_GENOS = 1,
          std::string PHENO_FNAME = "ps_l20.txt");

  FW_frna(int BASE, int GENOME_LENGTH, int ND, int TESTS, int SAMPLES, int ALL,
          int NEUTRAL_MUTATIONS = 1, int THRESHOLD = 2000000,
          int POPULATION_SIZE = 100, int SHAPE_LEVEL = 0,
          std::string FITNESS_MEASURE = "random", std::string OUTPATH = "./",
          std::string FILE_NAME = "0", std::string GENO_FNAME = "cs_l20.txt",
          int LOAD_GENOS = 0, std::string PHENO_FNAME = "ps_l20.txt",
          int RECORD_DIFFS_ONLY = 0, int DEBUG = 0);

  ~FW_frna();
  void Run();
  void RunDFS();
  void RunMono();
};

#endif
