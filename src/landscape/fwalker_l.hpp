#ifndef _fwalker_l_h_included_
#define _fwalker_l_h_included_

#include <math.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "fast_assemble_14.hpp"
#include "pheno_14.hpp"
#include "rna_fold.hpp"

class FW_l {
  std::unordered_set<std::string> u, v;
  std::map<std::string, int> phenos;
  std::vector<Phenotype> phenos_vector;
  std::map<std::string, int> genos;
  std::unordered_set<int> p;

  std::vector<double> fitness;

  std::vector<char> alphabet;

  int timing;
  clock_t start, middle, end;
  std::ofstream times;

  std::string fitness_assignment;

  std::string pheno_fname;

  int K;
  int L;
  int sparse;

  int Nt;
  int assembly_tests;

  int all;
  int N_p;
  int phenotype;
  int nd_ref;
  double max_fit;
  int max_fit_p;
  int threshold;

  int neutral_mutations_;

  int verbose;

  int targets;
  int target;
  int sources;
  int source;

  int print_component;

  std::string outpath_;
  std::string file_name_;

  std::string gp_map;

  std::ofstream stats;
  std::ofstream fp1;

  void SetAlphabet();
  void Init();
  void Init(int m_file_ref, int nd_file_ref);
  void Clear();
  int Lookup(Phenotype P);

  void IntSeq2String(int *&int_seq, std::string &s);
  void String2IntSeq(std::string &s, int *&int_seq);

  void AssignFitness();
  void AssignFitness(int max);
  void SetFitness();

  void SingleRun(std::string geno);

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

 public:
  FW_l(int BASE, int GENOME_LENGTH, int ND, int TESTS, int SAMPLES, int ALL,
       int NEUTRAL_MUTATIONS = 1, int THRESHOLD = 5000000,
       int ASSEMBLY_TESTS = 20, std::string FITNESS_ASSIGNMENT = "target",
       std::string OUTPATH = "./", std::string FILE_NAME = "0",
       std::string PHENO_FNAME = "pheno_list0.txt", std::string GP_MAP = "rna",
       int VERBOSE = 0);
  ~FW_l();
  void Run();
};

#endif
