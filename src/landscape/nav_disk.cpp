#include <fstream>
#include <iostream>
// #include <filesystem> // std=c++17 only, currently disabled
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
//#include "absl/flags/internal/usage.h"
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "fwalker.hpp"

ABSL_FLAG(int, seed, 1, "Random seed");
ABSL_FLAG(int, all, 0, "Test all phenotypes as targets or a random set");
ABSL_FLAG(int, sources, 100, "Number of sources per target to test");
ABSL_FLAG(int, targets, 100, "Number of targets to test");
ABSL_FLAG(int, nd_ref, 0, "Index of UND or del phenotype");
ABSL_FLAG(int, base, 4, "Base K of the genotype");
ABSL_FLAG(int, length, 12, "Length L of the genotype");
ABSL_FLAG(int, swaps, 0,
          "Number of genotype swaps to perform to decorrelate GP map");
ABSL_FLAG(int, dimension, 12, "Dimension D of the search");
ABSL_FLAG(bool, neutral_mutations, true,
          "Allow neutral mutations during the search");
ABSL_FLAG(int, run, 0, "Index for parallel runs");
ABSL_FLAG(int, threshold, 10000000,
          "Computational threshold for size of u + v");
ABSL_FLAG(bool, stop_at_fittest, true,
          "Stop the estimation when fittest phenotype found");
ABSL_FLAG(bool, record_phenotype_transitions, false,
          "Record transitions between phenotypes and fitnesses during search");
ABSL_FLAG(bool, reseed_fitness, false,
          "Fix fitness values for a given source and target");
ABSL_FLAG(bool, remap_phenotypes, false,
          "Remap phenotypes labels to contiguous integers");
ABSL_FLAG(bool, verbose, false, "Print progress");
ABSL_FLAG(std::string, fitness_assignment, "target",
          "\"source\", \"target\" or \"fixed\" of fitness assignment to "
          "phenotypes");
ABSL_FLAG(std::string, geno_fname, "geno_list0.txt",
          "Input genotype file path");
ABSL_FLAG(std::string, file_name, "0",
          "Tag to identify set of files together in outpath");
ABSL_FLAG(std::string, outpath, "./", "Path for saving outputs");

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
  config << "dimension: " << absl::GetFlag(FLAGS_dimension) << std::endl;
  config << "neutral_mutations: " << absl::GetFlag(FLAGS_neutral_mutations)
         << std::endl;
  config << "run: " << absl::GetFlag(FLAGS_run) << std::endl;
  config << "threshold: " << absl::GetFlag(FLAGS_threshold) << std::endl;
  config << "stop_at_fittest: " << absl::GetFlag(FLAGS_stop_at_fittest)
         << std::endl;
  config << "record_phenotype_transitions: "
         << absl::GetFlag(FLAGS_record_phenotype_transitions) << std::endl;
  config << "reseed_fitness: " << absl::GetFlag(FLAGS_reseed_fitness)
         << std::endl;
  config << "remap_phenotypes: " << absl::GetFlag(FLAGS_remap_phenotypes)
         << std::endl;
  config << "fitness_assignment: " << absl::GetFlag(FLAGS_fitness_assignment)
         << std::endl;
  config << "geno_fname: " << absl::GetFlag(FLAGS_geno_fname) << std::endl;
  config << "file_name: " << absl::GetFlag(FLAGS_file_name) << std::endl;
  config << "outpath: " << absl::GetFlag(FLAGS_outpath) << std::endl;
  config << "verbose: " << absl::GetFlag(FLAGS_verbose) << std::endl;
  config << "seed: " << absl::GetFlag(FLAGS_seed) << std::endl;
  config.close();
}

int main(int argc, char *argv[]) {
  std::string message =
      "Fitness landscape navigability for variable properties: correlations "
      "and dimensionality.";
  absl::SetProgramUsageMessage(message);
  absl::ParseCommandLine(argc, argv);

  clock_t start, end;
  start = clock();
  srand(absl::GetFlag(FLAGS_seed));
  rand();

  // Make directories, currently disabled, requires c++17
  // std::filesystem::create_directories(absl::GetFlag(FLAGS_outpath));

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

  // Constuct fitness landscape walker class
  FW fw(absl::GetFlag(FLAGS_base), absl::GetFlag(FLAGS_length),
        absl::GetFlag(FLAGS_nd_ref), absl::GetFlag(FLAGS_targets),
        absl::GetFlag(FLAGS_sources), absl::GetFlag(FLAGS_all),
        absl::GetFlag(FLAGS_swaps), absl::GetFlag(FLAGS_dimension),
        absl::GetFlag(FLAGS_neutral_mutations), absl::GetFlag(FLAGS_run),
        absl::GetFlag(FLAGS_threshold), absl::GetFlag(FLAGS_stop_at_fittest),
        absl::GetFlag(FLAGS_record_phenotype_transitions),
        absl::GetFlag(FLAGS_reseed_fitness),
        absl::GetFlag(FLAGS_remap_phenotypes),
        absl::GetFlag(FLAGS_fitness_assignment),
        absl::GetFlag(FLAGS_geno_fname), absl::GetFlag(FLAGS_file_name),
        absl::GetFlag(FLAGS_outpath), absl::GetFlag(FLAGS_verbose));

  fw.Run();

  // Finish program timer
  end = clock();
  out0 << "Program time = " << (double)(end - start) / CLOCKS_PER_SEC
       << " seconds" << endl;
  out0.close();
  return 0;
}
