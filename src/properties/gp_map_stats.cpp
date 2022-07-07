#include <fstream>
#include <iostream>
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
#include "absl/strings/str_format.h"
#include "base/utilities.hpp"
#include "cwalker.hpp"

ABSL_FLAG(int, seed, 1, "Random seed.");
ABSL_FLAG(int, nd_ref, 0, "Index of UND or del phenotype.");
ABSL_FLAG(bool, nd_include, false, "Include statistics for ND phenotype.");
ABSL_FLAG(int, nd_threshold, 1,
          "Maximum number of ground state phenotypes permitted.");
ABSL_FLAG(int, base, 4, "Base K of the genotype.");
ABSL_FLAG(int, length, 12, "Length L of the genotype.");
ABSL_FLAG(int, swaps, 0,
          "Number of genotype swaps to perform to decorrelate GP map.");
ABSL_FLAG(bool, write_components, false,
          "Write the genotypes that belong to each component to file.");
ABSL_FLAG(bool, estimate_component_stats, false,
          "Write statistics at the compoent level to file.");
ABSL_FLAG(bool, record_phenotype_transitions, false,
          "Record transitions between phenotypes and fitnesses during search.");
ABSL_FLAG(bool, verbose, false, "Print progress.");
ABSL_FLAG(bool, debug, false, "Run in debug mode.");
ABSL_FLAG(std::string, geno_fname, "geno_list0.txt",
          "Input genotype file path.");
ABSL_FLAG(std::string, file_name, "0",
          "Tag to identify set of files together in outpath.");
ABSL_FLAG(std::string, outpath, "./", "Path for saving outputs.");

void WriteFlags(std::string configfile) {
  std::ofstream config(configfile.c_str());
  config << "---" << std::endl;
  config << "base: " << absl::GetFlag(FLAGS_base) << std::endl;
  config << "length: " << absl::GetFlag(FLAGS_length) << std::endl;
  config << "nd_ref: " << absl::GetFlag(FLAGS_nd_ref) << std::endl;
  config << "nd_include: " << absl::GetFlag(FLAGS_nd_include) << std::endl;
  config << "nd_threshold: " << absl::GetFlag(FLAGS_nd_threshold) << std::endl;
  config << "swaps: " << absl::GetFlag(FLAGS_swaps) << std::endl;
  config << "write_components: " << absl::GetFlag(FLAGS_write_components)
         << std::endl;
  config << "estimate_component_stats: "
         << absl::GetFlag(FLAGS_estimate_component_stats) << std::endl;
  config << "geno_fname: " << absl::GetFlag(FLAGS_geno_fname) << std::endl;
  config << "file_name: " << absl::GetFlag(FLAGS_file_name) << std::endl;
  config << "outpath: " << absl::GetFlag(FLAGS_outpath) << std::endl;
  config << "seed: " << absl::GetFlag(FLAGS_seed) << std::endl;
  config << "verbose: " << absl::GetFlag(FLAGS_verbose) << std::endl;
  config << "debug: " << absl::GetFlag(FLAGS_verbose) << std::endl;
  config.close();
}

int main(int argc, char *argv[]) {
  std::string message =
      "Measurement of GP map properties at the component and phenotype level "
      "from a geno_list.";

  absl::SetProgramUsageMessage(message);
  absl::ParseCommandLine(argc, argv);

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

  // Check nd_included if component level stats
  if (absl::GetFlag(FLAGS_estimate_component_stats)) {
    assert(absl::GetFlag(FLAGS_nd_include));
  }

  // Make neutral componenent object
  nc NC(absl::GetFlag(FLAGS_base), absl::GetFlag(FLAGS_length),
        absl::GetFlag(FLAGS_nd_ref), absl::GetFlag(FLAGS_nd_include),
        absl::GetFlag(FLAGS_nd_threshold), absl::GetFlag(FLAGS_swaps),
        absl::GetFlag(FLAGS_write_components),
        absl::GetFlag(FLAGS_estimate_component_stats),
        absl::GetFlag(FLAGS_geno_fname), absl::GetFlag(FLAGS_file_name),
        absl::GetFlag(FLAGS_outpath), absl::GetFlag(FLAGS_verbose),
        absl::GetFlag(FLAGS_debug));

  NC.Run();

  // Finish program timer
  end = clock();
  out0 << "Program time = " << (double)(end - start) / CLOCKS_PER_SEC
       << " seconds" << std::endl;
  out0.close();
  return 0;
}
