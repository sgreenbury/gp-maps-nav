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
#include "fast_assemble_14.hpp"
#include "fwalker_frna.hpp"
#include "pheno_14.hpp"
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

std::vector<Phenotype> LoadPhenos(std::string pheno_fname) {
  std::vector<Phenotype> phenos;
  // Load the phenotypes
  int N_p = 0;
  while (1) {
    Phenotype P;
    std::string buffer = absl::StrFormat("%s/pheno%d.txt", pheno_fname, N_p);
    std::ifstream in(buffer);
    // If cannot open the file, then finished loading phenotypes
    if (!in.good()) {
      break;
    } else {
      in.close();
    }
    P.fLoad(buffer.c_str());
    phenos.push_back(P);
    std::cout << N_p << std::endl;
    phenos[N_p].PrintPic();
    N_p++;
  }
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
  std::string n_g;

  int L = 8;
  int K = 8;

  int Nt = 2;
  int assembly_tests = 1000;

  // Load phenotypes
  std::vector<Phenotype> phenos = LoadPhenos("gp_maps/s_2_8/phenotypes");

  // Load genotypes
  std::vector<int> genos = LoadGenos("gp_maps/s_2_8/geno_list0.txt");

  // Swap labels for mirror image GP map on disk
  for (int i = 0; i < genos.size(); i++) {
    if (genos[i] == 7) {
      genos[i] = 8;
      continue;
    }
    if (genos[i] == 8) genos[i] = 7;
  }

  int *seq;
  seq = new int[L];

  start = clock();
  for (unsigned int i = 0; i < (int)pow((double)K, L); i++) {
    // Get base sequence
    int_to_basek(i, seq, L, K);

    if (i % 100000 == 0) {
      std::cout << "Geno : " << i << std::endl;
    }

    // Get phenotype from assembly
    Phenotype A = Assemble(seq, Nt, assembly_tests);

    int n_p = -1;
    for (unsigned int i = 0; i < phenos.size(); i++) {
      if (phenos[i].MatchRots(A) == 1) {
        n_p = i;
        break;
      }
    }

    int g_orig = genos[i];
    int g_new = n_p;
    if (genos[i] != n_p) {
      std::cout << "New:" << std::endl;
      A.PrintPic();
      std::cout << "Old:" << std::endl;
      phenos[g_orig].PrintPic();

      std::cout << i << "   ";
      for (int k = 0; k < L; k++) {
        std::cout << seq[k];
        if (k < L - 1) std::cout << ",";
      }
      std::cout << "\t";
      std::cout << g_new << "    " << g_orig << std::endl;
      std::cin.get();
    }
  }

  return 0;
}
