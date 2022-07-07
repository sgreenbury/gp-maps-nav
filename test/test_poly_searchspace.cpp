#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "landscape/fast_assemble_14.hpp"
#include "landscape/fwalker_frna.hpp"
#include "landscape/pheno_14.hpp"
#include "landscape/rna_fold.hpp"

#include "gtest/gtest.h"

std::vector<int> LoadGenos(std::string geno_fname) {
  // Load GP map, remap phenotypes such that form ordered list
  std::ifstream in;
  in.open(geno_fname);
  std::string line;
  std::vector<int> genos;
  while (std::getline(in, line)) {
    int phenotype = String2Number(line);
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

TEST(Poly, Searchspace) {
  // Test S_2_8 polyomino searchspace on disk
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
  for (unsigned int i = 0; i < genos.size(); i++) {
    if (genos[i] == 7) {
      genos[i] = 8;
      continue;
    }
    if (genos[i] == 8) genos[i] = 7;
  }

  int *seq;
  seq = new int[L];

  int mismatch = 0;
  int S = (int)pow((double)K, L);

  // Just test 200,000 genotypes to avoid long test time but all phenotypes
  // covered
  // int TESTS = S;
  int TESTS = 200000;
  for (int idx = 0; idx < TESTS; idx++) {
    int i = rand() % S;

    // Get base sequence
    int_to_basek(i, seq, L, K);

    if (idx % 100000 == 0) {
      std::cout << "Geno : " << i << std::endl;
    }

    // Get phenotype from assembly
    Phenotype A = Assemble(seq, Nt, assembly_tests);

    int n_p = -1;
    for (unsigned int pidx = 0; pidx < phenos.size(); pidx++) {
      if (phenos[pidx].MatchRots(A) == 1) {
        n_p = pidx;
        break;
      }
    }

    int g_orig = genos[i];
    int g_new = n_p;
    if (genos[i] != n_p) {
      mismatch++;

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
    }
  }

  EXPECT_EQ(mismatch, 0);
}
