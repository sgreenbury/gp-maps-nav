#include <fstream>
#include <iostream>
// #include <filesystem> // std=c++17 only, currently disabled
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
//#include "absl/flags/internal/usage.h"
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "fwalker_frna.hpp"
#include "rna_fold.hpp"

ABSL_FLAG(int, seed, 1, "Random seed");
ABSL_FLAG(unsigned int, N, 1000, "Number of tests");
ABSL_FLAG(unsigned int, level, 1, "Level of output shape structure");
ABSL_FLAG(std::string, seq,
          "UCUUCAUGAGGUACUUCACAUCAUAAAUGACAGGAAAAAACAAUCGAAGGAUCUCAAAGA",
          "Example sequence");

int main(int argc, char *argv[]) {
  std::string message = "Test 'rna_fold.cpp'";
  absl::SetProgramUsageMessage(message);
  absl::ParseCommandLine(argc, argv);

  clock_t start, end;
  start = clock();
  srand(absl::GetFlag(FLAGS_seed));
  rand();

  // Test lev distance
  std::string s1 = "tes";
  std::string s2 = "test";
  std::string s3 = "aaaaa";
  std::string s4 = "bbbbb";
  double t;
  t = CompareStringsLevDist(s1, s2);
  std::cout << s1 << " " << s2 << " " << t << "\n";
  t = CompareStringsLevDist(s1, s1);
  std::cout << s1 << " " << s1 << " " << t << "\n";
  t = CompareStringsLevDist(s3, s4);
  std::cout << s3 << " " << s4 << " " << t << "\n";
  t = CompareStringsLevDist(s1, s4);
  std::cout << s3 << " " << s4 << " " << t << "\n";

  s1 = "[_[_[[[[[_[_[_[]_]]_]_e]_]_]_]_]_]_";
  s2 = "_[_[_[_[[]_]_]_]_]_";
  t = CompareStringsLevDist(s1, s2);
  std::cout << s1 << " " << s2 << " " << t << "\n";

  // Test folds
  std::string n_g(absl::GetFlag(FLAGS_seq));

  int level = absl::GetFlag(FLAGS_level);
  int L = n_g.size();
  std::string n_p;
  std::string n_s;
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  unsigned int N = absl::GetFlag(FLAGS_N);
  std::cout << "Sequence: " << n_g << ", length = " << L << std::endl;
  start = clock();
  for (unsigned int i = 0; i < N; i++) {
    std::strcpy(char_seq, n_g.c_str());
    char_seq[L] = '\0';
    Fold(char_seq, phys_seq);
    n_p.assign(phys_seq);
    if (i == N - 1) std::cout << "Fold    : " << n_p << std::endl;
  }

  // Finish program timer
  end = clock();
  std::cout << "Fold time = " << (double)(end - start) / CLOCKS_PER_SEC
            << " seconds" << std::endl;

  start = clock();
  for (unsigned int i = 0; i < N; i++) {
    std::strcpy(char_seq, n_g.c_str());
    char_seq[L] = '\0';
    Fold(char_seq, phys_seq);
    n_p.assign(phys_seq);

    n_s = shape(n_p, level);
    if (i == N - 1) std::cout << "Shape    :" << n_s << std::endl;
  }
  end = clock();
  std::cout << "Fold and shape time = "
            << (double)(end - start) / CLOCKS_PER_SEC << " seconds"
            << std::endl;

  free(char_seq);
  free(phys_seq);

  return 0;
}
