#ifndef _evo_rna_h_included_
#define _evo_rna_h_included_

#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "absl/strings/str_format.h"

#include "base/utilities.hpp"
#include "individual.hpp"
#include "landscape/rna_fold.hpp"

class Evo_RNA {
 public:
  int id;
  int test;
  int generations;
  int generation;
  double threshold;

  int complete;
  int decrease;
  double previous_max;
  double best_max;
  int N;
  double mu;
  Individual source;
  int K;
  int L;
  Individual dummy;

  int verbose;
  int vverbose = 0;
  int track_mrca;

  int neutral_mutations_ = 1;

  int debug;

  std::string fitness_measure;

  std::string outpath;
  std::string file_name;

  std::string statsfile;
  std::string mrcafile;
  std::string summaryfile;

  std::ofstream file;
  std::ofstream mrca_fp;
  std::ofstream summary_fp;

  int max_fit_more_than_threshold;
  double maj_max_fit;
  double previous_maj_max_fit;

  int default_precision = 3;
  int fitness_precision = 12;

  std::string target;

  std::unordered_set<std::string> pheno_set;
  std::unordered_map<std::string, double> phenotype_fitnesses;
  std::unordered_map<std::string, int> phenotype_counter;

  std::vector<Individual> current_population;
  std::vector<Individual> next_population;
  std::vector<double> fitnesses;

  std::vector<std::vector<int>> history;
  std::vector<std::vector<Individual>> history_genos;
  // int previous_mrca_generation;
  // int mrca_generation;
  std::pair<int, int> current_mrca;
  std::pair<int, int> previous_mrca;

  std::pair<int, int> CalculateMRCA(int generation,
                                    std::vector<std::vector<int>>& history);
  std::deque<std::pair<int, int>> AncestryPath(std::pair<int, int>& start,
                                               std::pair<int, int>& end);
  void AssignMRCA(std::pair<int, int>& mrca);
  void WriteStatsMRCA(std::pair<int, int>& mrca);
  void WriteStatsMRCAHeader();
  void WriteSummaryFileHeader();

  double CompareStrings(std::string s1, std::string s2);
  void CalculateFitness(Individual& ind);

  void Init(int ID, int TEST, int POPULATION_SIZE, int GENERATIONS,
            Individual SOURCE, std::string TARGET,
            std::set<std::string>& PHENO_SET, double MU,
            std::string FITNESS_MEASURE = "random",
            double FITNESS_THRESHOLD = 0.5, int VERBOSE = 0, int TRACK_MRCA = 0,
            int NEUTRAL_MUTATIONS = 1, int DEBUG = 0,
            const std::string OUTPATH = "./",
            const std::string FILE_NAME = "0");
  void WriteStatsFileHeader();
  void InitFirstPopulation();
  void PopulationG2P();
  void CalculateFitnesses();
  void MakeFitnessVector();
  int PickMember();
  void PopulateNextGeneration();
  void AssignNextGeneration2Current();

  void OneGeneration();

  void WriteStats();
  void WriteSummary();
  void Run();
};

#endif
