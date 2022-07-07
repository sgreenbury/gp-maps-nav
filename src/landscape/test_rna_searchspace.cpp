#include <fstream>
#include <iostream>
#include <vector>
// #include <filesystem> // std=c++17 only, currently disabled
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
//#include "absl/flags/internal/usage.h"
#include <map>
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "fwalker_frna.hpp"
#include "rna_fold.hpp"

ABSL_FLAG(int, seed, 1, "Random seed");
ABSL_FLAG(unsigned int, N, 1000, "Number of tests");
ABSL_FLAG(unsigned int, level, 1, "Level of output shape structure");

std::vector<int> LoadGenos(std::string geno_fname_) {
  // Load GP map, remap phenotypes such that form ordered list
  std::ifstream in;
  in.open(geno_fname_);
  std::string line;
  int max = 0;
  std::vector<int> genos;
  while (std::getline(in, line)) {
    int phenotype = String2Number(line);
    if (phenotype > max) max = phenotype;
    genos.push_back(phenotype);
  }
  in.close();
  return genos;
}

std::map<std::string, int> LoadPhenos(std::string pheno_fname_) {
  // Load GP map, remap phenotypes such that form ordered list
  std::ifstream in;
  in.open(pheno_fname_);
  std::string line;
  int max = 0;
  std::map<std::string, int> phenos;
  int count = 0;
  while (std::getline(in, line)) {
    int phenotype = String2Number(line);
    phenos[line] = count;
    count++;
  }
  in.close();
  return phenos;
}

int main(int argc, char *argv[]) {
  std::string message = "Test 'rna_fold.cpp'";
  absl::SetProgramUsageMessage(message);
  absl::ParseCommandLine(argc, argv);

  clock_t start, end;
  start = clock();
  srand(absl::GetFlag(FLAGS_seed));
  rand();

  // Test folds
  std::string n_g;  //(absl::GetFlag(FLAGS_seq));

  int level = absl::GetFlag(FLAGS_level);
  int L = 12;
  int K = 4;
  std::string n_p;
  std::string n_s;
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  // Load genotypes
  std::vector<int> genos = LoadGenos("gp_maps/RNA_12/geno_list0.txt");

  // Load phenotypes
  std::map<std::string, int> phenos =
      LoadPhenos("gp_maps/RNA_12/pheno_list0.txt");

  int *seq;
  seq = new int[L];

  std::vector<char> alpha = std::vector<char>(4, '\0');
  // RNA12   : {'A': 0, 'U': 1, 'C': 2, 'G': 3}
  alpha[0] = 'A';
  alpha[1] = 'U';
  alpha[2] = 'C';
  alpha[3] = 'G';

  start = clock();
  for (unsigned int i = 0; i < (int)pow((double)K, L); i++) {
    // Get base sequence
    int_to_basek(i, seq, L, K);

    if (i % 100000 == 0) {
      std::cout << "Geno : " << i << std::endl;
    }

    // Convert to RNA chars
    for (int j = 0; j < L; j++) {
      char_seq[j] = alpha[seq[j]];
    }
    char_seq[L] = '\0';

    n_g.assign(char_seq);

    Fold(char_seq, phys_seq);
    n_p.assign(phys_seq);
    int n_p_lookup = phenos[n_p];
    int g_orig = genos[i];
    int g_new = n_p_lookup;
    if (g_new != g_orig) {
      std::cout << i << "    " << n_g << "    " << n_p << "    " << g_new
                << "    " << g_orig << std::endl;
    }
  }

  free(char_seq);
  free(phys_seq);

  return 0;
}
