#include <fstream>
#include <iostream>
// #include <filesystem> // std=c++17 only, currently disabled
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
//#include "absl/flags/internal/usage.h"
#include <regex>
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "rna_fold.hpp"

ABSL_FLAG(int, seed, 1, "Seed");
ABSL_FLAG(std::string, path, "", "Path for sequences.");
ABSL_FLAG(std::string, outpath, "", "Path for outputs.");
ABSL_FLAG(int, begin, 0, "Starting number.");
ABSL_FLAG(int, end, 1000000000, "Finish number.");

int main(int argc, char *argv[]) {
  std::string message = "Test 'rna_fold.cpp'";
  absl::SetProgramUsageMessage(message);
  absl::ParseCommandLine(argc, argv);

  clock_t start, end;
  start = clock();

  srand(absl::GetFlag(FLAGS_seed));
  rand();

  unsigned int L = 1;

  std::string g;
  std::string p;
  std::string s;

  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  start = clock();

  std::ifstream infile(absl::GetFlag(FLAGS_path));
  // std::string outname = std::regex_replace(absl::GetFlag(FLAGS_path),
  // 					   std::regex(".csv"), "_all.csv");
  std::string outname = absl::GetFlag(FLAGS_outpath);
  std::ofstream out(outname);
  int counter = 0;
  while (std::getline(infile, g)) {
    if ((counter < absl::GetFlag(FLAGS_begin)) |
        (counter >= absl::GetFlag(FLAGS_end))) {
      counter++;
      continue;
    }
    if (counter % 1000 == 0) {
      std::cout << counter << "\t" << g << std::endl;
    }
    if (g.size() != L) {
      L = g.size();
      free(char_seq);
      free(phys_seq);
      char_seq = (char *)malloc(sizeof(char) * (L + 1));
      char_seq[L] = '\0';
      phys_seq = (char *)malloc(sizeof(char) * (L + 1));
      phys_seq[L] = '\0';
    }

    std::strcpy(char_seq, g.c_str());
    char_seq[L] = '\0';
    Fold(char_seq, phys_seq);
    p.assign(phys_seq);
    out << g << ",";
    for (unsigned int i = 0; i <= 5; i++) {
      out << shape(p, i);
      if (i == 5) {
        out << std::endl;
      } else {
        out << ",";
      }
    }
    counter++;
  }
  end = clock();
  std::cout << "All seqs folded = " << (double)(end - start) / CLOCKS_PER_SEC
            << " seconds" << std::endl;

  free(char_seq);
  free(phys_seq);

  return 0;
}
