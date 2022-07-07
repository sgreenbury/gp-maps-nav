#ifndef _cwalker_h_included_
#define _cwalker_h_included_

#include <math.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include "base/utilities.hpp"

class nc {
  std::set<int> w, u, v, s, sp;
  std::vector<int> genotype_list;
  std::vector<int> genotype_component_list;
  std::vector<int> phenotype_component_list;
  std::vector<int> component_volumes;
  std::vector<int> component_surfaces;
  std::vector<long> component_transitions;
  std::vector<long> phenotype_transitions;
  std::vector<std::vector<int> > component_matrix_transitions;
  std::vector<long> component_ind_transitions;

  int K;
  int L;
  int sparse;
  int nd_ref;
  int nd_include;
  int threshold_nd;

  int N_p;
  int phenotype;
  int S;
  double r_c;
  double e_wagner_c;
  double e_cowperthwaite_c;
  double e_wagner_c_nond;
  double e_cowperthwaite_c_nond;

  int swaps;

  int verbose_;
  int debug_;

  // Paths
  std::string outpath_;
  std::string file_name_;
  std::string geno_fname_;

  // ofstream pointer
  std::ofstream fp1;

  int f_p;
  double r_p;
  double e_wagner_p;
  double e_cowperthwaite_p;
  double e_wagner_p_nond;
  double e_cowperthwaite_p_nond;

  int write_components;
  int estimate_component_stats;

  void Init(int K, int L);
  void Clear();

  void Transfer(std::set<int> &from, std::set<int> &to, int x);
  void fStartComponent(int p, int n);
  void fAddEdge(int x, int y);
  void fEndComponent();

  void CalculateComponentTransitions();
  void CalculateComponentStatistics();
  void CalculatePhenotypeStatistics();

  void WriteCompStatsHeader(std::ofstream &file);
  void WritePhenoStatsHeader(std::ofstream &file);
  void fPrintComponentStatistics(std::ofstream &file);
  void fPrintPhenotypeStatistics(std::ofstream &file);

 public:
  nc(int base, int genome_length, int ND_REF, int ND_INCLUDE, int ND_THRESHOLD,
     int SWAPS, int WRITE_COMPONENTS, int ESTIMATE_COMPONENT_STATS,
     std::string GENO_FNAME, std::string FILE_NAME, std::string OUTPATH,
     int VERBOSE = 0, int DEBUG = 0);

  ~nc();
  void Run();
};

#endif
