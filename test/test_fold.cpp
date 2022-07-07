#include <fstream>
#include <iostream>
// #include <filesystem> // std=c++17 only, currently disabled
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
//#include "absl/flags/internal/usage.h"
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "gtest/gtest.h"
#include "landscape/fwalker_frna.hpp"
#include "landscape/rna_fold.hpp"

TEST(Fold, DotBracketFold) {
  std::string n_g(
      "UCUUCAUGAGGUACUUCACAUCAUAAAUGACAGGAAAAAACAAUCGAAGGAUCUCAAAGA");
  std::string n_p;
  std::string n_s;
  int L = n_g.size();

  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  std::strcpy(char_seq, n_g.c_str());
  char_seq[L] = '\0';
  Fold(char_seq, phys_seq);
  n_p.assign(phys_seq);
  free(char_seq);
  free(phys_seq);
  EXPECT_EQ(n_p,
            "((((..((((((.((((...(((....)))..(.......)....)))).))))))))))");
}

TEST(Fold, ShapeLevel1) {
  std::string n_g(
      "UCUUCAUGAGGUACUUCACAUCAUAAAUGACAGGAAAAAACAAUCGAAGGAUCUCAAAGA");
  std::string n_p;
  std::string n_s;
  int L = n_g.size();

  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  std::strcpy(char_seq, n_g.c_str());
  char_seq[L] = '\0';
  Fold(char_seq, phys_seq);
  n_p.assign(phys_seq);
  n_s = shape(n_p, 1);
  std::cout << n_s << std::endl;
  std::cout << "[_[_[_[]_[]_]_]]" << std::endl;
  free(char_seq);
  free(phys_seq);
  EXPECT_EQ(n_p,
            "((((..((((((.((((...(((....)))..(.......)....)))).))))))))))");
  EXPECT_EQ(n_s, "[_[_[_[]_[]_]_]]");
}
