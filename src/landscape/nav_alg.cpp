#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "fwalker_l.hpp"

ABSL_FLAG(int, seed, 1, "Random seed");
ABSL_FLAG(bool, all, false, "Test all phenotypes as targets or a random set");
ABSL_FLAG(int, sources, 100, "Number of sources per target to test");
ABSL_FLAG(int, targets, 100, "Number of targets to test");
ABSL_FLAG(int, nd_ref, 0, "Index of UND or del phenotype");
ABSL_FLAG(int, base, 4, "Base K of the genotype");
ABSL_FLAG(int, length, 12, "Length L of the genotype");
ABSL_FLAG(int, swaps, 0,
          "Number of genotype swaps to perform to decorrelate GP map");
ABSL_FLAG(bool, neutral_mutations, true,
          "Allow neutral mutations during the search");
ABSL_FLAG(int, threshold, 10000000,
          "Computational threshold for size of u + v");
ABSL_FLAG(
    int, assembly_tests, 20,
    "Number of times to repeat assembly of polyomino to check for determinism");
ABSL_FLAG(bool, verbose, false, "Print progress");
ABSL_FLAG(std::string, fitness_assignment, "target",
          "\"source\", \"target\" or \"fixed\" of fitness assignment to "
          "phenotypes");
ABSL_FLAG(std::string, file_name, "0",
          "Tag to identify set of files together in outpath");
ABSL_FLAG(std::string, outpath, "./", "Path for saving outputs");
ABSL_FLAG(std::string, pheno_fname, "pheno_list0.txt",
          "Input phenotype file path");
ABSL_FLAG(std::string, gp_map, "rna",
          "The type of GP map, options: \"rna\" or \"polyomino\"");

void WriteFlags(std::string configfile) {
  std::ofstream config(configfile.c_str());
  config << "---" << std::endl;
  config << "base: " << absl::GetFlag(FLAGS_base) << std::endl;
  config << "length: " << absl::GetFlag(FLAGS_length) << std::endl;
  config << "nd_ref: " << absl::GetFlag(FLAGS_nd_ref) << std::endl;
  config << "targets: " << absl::GetFlag(FLAGS_targets) << std::endl;
  config << "sources: " << absl::GetFlag(FLAGS_sources) << std::endl;
  config << "all: " << absl::GetFlag(FLAGS_all) << std::endl;
  config << "swaps: " << absl::GetFlag(FLAGS_swaps) << std::endl;
  config << "neutral_mutations: " << absl::GetFlag(FLAGS_neutral_mutations)
         << std::endl;
  config << "threshold: " << absl::GetFlag(FLAGS_threshold) << std::endl;
  config << "assembly_tests: " << absl::GetFlag(FLAGS_assembly_tests)
         << std::endl;
  config << "fitness_assignment: " << absl::GetFlag(FLAGS_fitness_assignment)
         << std::endl;
  config << "pheno_fname: " << absl::GetFlag(FLAGS_pheno_fname) << std::endl;
  config << "file_name: " << absl::GetFlag(FLAGS_file_name) << std::endl;
  config << "outpath: " << absl::GetFlag(FLAGS_outpath) << std::endl;
  config << "gp_map: " << absl::GetFlag(FLAGS_gp_map) << std::endl;
  config << "verbose: " << absl::GetFlag(FLAGS_verbose) << std::endl;
  config << "seed: " << absl::GetFlag(FLAGS_seed) << std::endl;
  config.close();
}

int main(int argc, char *argv[]) {
  std::string message =
      "Fitness landscape navigability for algorithmic folding of RNA or "
      "polyomino assembly.";
  absl::SetProgramUsageMessage(message);
  absl::ParseCommandLine(argc, argv);

  // Set-up time points
  clock_t start, end;
  start = clock();
  srand(absl::GetFlag(FLAGS_seed));
  rand();

  // Initialize outfile and config file
  std::string outfile =
      absl::StrFormat("%sout_%s.txt", absl::GetFlag(FLAGS_outpath),
                      absl::GetFlag(FLAGS_file_name));
  std::string configfile =
      absl::StrFormat("%sconfig_%s.yaml", absl::GetFlag(FLAGS_outpath),
                      absl::GetFlag(FLAGS_file_name));

  // Write flags to file
  WriteFlags(configfile);

  // Open outfile
  std::ofstream out0(outfile.c_str());

  // Make FW object
  FW_l fw_l(absl::GetFlag(FLAGS_base), absl::GetFlag(FLAGS_length),
            absl::GetFlag(FLAGS_nd_ref), absl::GetFlag(FLAGS_targets),
            absl::GetFlag(FLAGS_sources), absl::GetFlag(FLAGS_all),
            absl::GetFlag(FLAGS_neutral_mutations),
            absl::GetFlag(FLAGS_threshold), absl::GetFlag(FLAGS_assembly_tests),
            absl::GetFlag(FLAGS_fitness_assignment),
            absl::GetFlag(FLAGS_outpath), absl::GetFlag(FLAGS_file_name),
            absl::GetFlag(FLAGS_pheno_fname), absl::GetFlag(FLAGS_gp_map),
            absl::GetFlag(FLAGS_verbose));

  // Perform searches
  fw_l.Run();

  // Finish program timer
  end = clock();
  out0 << "Program time = " << (double)(end - start) / CLOCKS_PER_SEC
       << " seconds" << endl;
  out0.close();
  return 0;
}
