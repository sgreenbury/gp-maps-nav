#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
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

std::map<std::string, int> LoadPhenos(std::string pheno_fname) {
  // Load GP map, remap phenotypes such that form ordered list
  std::ifstream in;
  in.open(pheno_fname);
  std::string line;
  std::map<std::string, int> phenos;
  int count = 0;
  while (std::getline(in, line)) {
    phenos[line] = count;
    count++;
  }
  in.close();
  return phenos;
}

TEST(RNA, Searchspace) {
  // Test full RNA 12 searchspace on disk
  std::string n_g;

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

  int mismatch = 0;
  int S = (int)pow((double)K, L);
  for (int i = 0; i < S; i++) {
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
      mismatch++;
      std::cout << i << "    " << n_g << "    " << n_p << "    " << g_new
                << "    " << g_orig << std::endl;
    }
  }

  free(char_seq);
  free(phys_seq);

  EXPECT_EQ(mismatch, 0);
}
