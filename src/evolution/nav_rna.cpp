#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/usage.h"
#include "absl/strings/str_format.h"
#include "evo_rna.hpp"
#include "individual.hpp"
#include "landscape/rna_fold.hpp"

ABSL_FLAG(int, seed, 1, "Random seed");
ABSL_FLAG(int, sources, 50, "Number of sources per target to test");
ABSL_FLAG(int, targets, 20, "Number of targets to test");
ABSL_FLAG(int, generations, 10000, "Number of generations for evolution");
ABSL_FLAG(int, base, 4, "Base K of the genotype");
ABSL_FLAG(int, length, 20, "Length L of the genotype");
ABSL_FLAG(int, population_size, 100, "Size of population");
ABSL_FLAG(double, mu, 0.05, "Point mutation rate per site");
ABSL_FLAG(double, fitness_threshold, 0.5,
          "Proportion of population required for fitness record");
ABSL_FLAG(std::string, fitness_measure, "random",
          "Fitness measure used: 'random' or 'hamming'");
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
ABSL_FLAG(bool, debug, false, "Run in debug mode");
ABSL_FLAG(bool, track_mrca, false,
          "Track Most Recent Common Ancestor during evolution");
ABSL_FLAG(std::string, source_fitness, "zero",
          "Scheme to set source fitness: 'zero' or 'random'");
ABSL_FLAG(std::string, pheno_fname, "",
          "Input phenotype file path for setting sources and targets");
ABSL_FLAG(std::string, file_name, "0",
          "Tag to identify set of files together in outpath");
ABSL_FLAG(std::string, outpath, "./", "Path for saving outputs");

void WriteFlags(std::string configfile) {
  std::ofstream config(configfile.c_str());
  config << "---" << std::endl;
  config << "base: " << absl::GetFlag(FLAGS_base) << std::endl;
  config << "length: " << absl::GetFlag(FLAGS_length) << std::endl;
  // config << "nd_ref: " << absl::GetFlag(FLAGS_nd_ref) << std::endl;
  config << "targets: " << absl::GetFlag(FLAGS_targets) << std::endl;
  config << "sources: " << absl::GetFlag(FLAGS_sources) << std::endl;
  config << "generations: " << absl::GetFlag(FLAGS_generations) << std::endl;
  config << "population_size: " << absl::GetFlag(FLAGS_population_size)
         << std::endl;
  config << "fitness_threshold: " << absl::GetFlag(FLAGS_fitness_threshold)
         << std::endl;
  config << "fitness_measure: " << absl::GetFlag(FLAGS_fitness_measure)
         << std::endl;
  config << "mu: " << absl::GetFlag(FLAGS_mu) << std::endl;
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
  config << "source_fitness: " << absl::GetFlag(FLAGS_source_fitness)
         << std::endl;
  config << "pheno_fname: " << absl::GetFlag(FLAGS_pheno_fname) << std::endl;
  config << "file_name: " << absl::GetFlag(FLAGS_file_name) << std::endl;
  config << "outpath: " << absl::GetFlag(FLAGS_outpath) << std::endl;
  config << "verbose: " << absl::GetFlag(FLAGS_verbose) << std::endl;
  config << "track_mrca: " << absl::GetFlag(FLAGS_track_mrca) << std::endl;
  config << "seed: " << absl::GetFlag(FLAGS_seed) << std::endl;
  config.close();
}

char RandomBase() {
  int el = rand() % 4;
  char ch = 'b';
  switch (el) {
    case 0:
      ch = 'G';
      break;
    case 1:
      ch = 'C';
      break;
    case 2:
      ch = 'A';
      break;
    case 3:
      ch = 'U';
      break;
  }
  return ch;
}

int main(int argc, char *argv[]) {
  // Program usage message setting
  std::string message =
      "Fitness landscape navigability under evolutionary dynamics.";
  absl::SetProgramUsageMessage(message);

  // Stores all non-flag arguments of command line in cl
  std::vector<char *> cl = absl::ParseCommandLine(argc, argv);

  // Set-up clock and seeding
  clock_t start, end;
  start = clock();
  srand(absl::GetFlag(FLAGS_seed));
  rand();

  // Constants
  std::string outfile =
      absl::StrFormat("%sout_%s.txt", absl::GetFlag(FLAGS_outpath),
                      absl::GetFlag(FLAGS_file_name));
  std::string configfile =
      absl::StrFormat("%sconfig_%s.yaml", absl::GetFlag(FLAGS_outpath),
                      absl::GetFlag(FLAGS_file_name));

  // Open outfile
  std::ofstream out0(outfile.c_str());

  // Write flags to file
  WriteFlags(configfile);

  // Initialize Vienna 1.85, not needed with Vienna>=2
  initialize_fold(absl::GetFlag(FLAGS_length));

  // Load phenotypes for sources and targets
  std::string pheno_fname =
      absl::StrFormat("%s", absl::GetFlag(FLAGS_pheno_fname));
  std::ifstream in(pheno_fname.c_str());
  std::string line;
  std::vector<string> phenos;
  while (std::getline(in, line)) {
    out0 << line << std::endl;
    phenos.push_back(line);
  }
  in.close();

  std::set<std::string> pheno_set;
  for (int i = 0; i < phenos.size(); i++) pheno_set.insert(phenos[i]);

  // Seed new random number gen:
  // Loop over number of targets
  for (unsigned int i = 0; i < absl::GetFlag(FLAGS_targets); i++) {
    std::string target;
    do {
      std::set<std::string>::const_iterator it(pheno_set.begin());
      int el = rand() % pheno_set.size();
      std::advance(it, el);
      target = *it;
    } while (target == std::string(absl::GetFlag(FLAGS_length), '.'));

    // Loop over sources per target
    for (unsigned int j = 0; j < absl::GetFlag(FLAGS_sources); j++) {
      // Write down seed for random number generator
      // Write down the target phenotype.
      Individual ind(absl::GetFlag(FLAGS_base), absl::GetFlag(FLAGS_length));

      // Recopy
      do {
        ind.RandomGenotype();
        ind.G2P();
      } while (ind.phenotype == target);

      // Set the initial fitness for the source based upon source fitness flag
      if (absl::GetFlag(FLAGS_source_fitness) == "random") {
        ind.fitness = (double)rand() / RAND_MAX;
      } else if (absl::GetFlag(FLAGS_source_fitness) == "zero") {
        ind.fitness = 0.0;
      } else {
        ind.fitness = 0.0;
      }

      // Construct the evolution simulation object
      Evo_RNA evo;

      // Inialize the evoltionary run
      evo.Init(i, j, absl::GetFlag(FLAGS_population_size),
               absl::GetFlag(FLAGS_generations), ind, target, pheno_set,
               absl::GetFlag(FLAGS_mu), absl::GetFlag(FLAGS_fitness_measure),
               absl::GetFlag(FLAGS_fitness_threshold),
               absl::GetFlag(FLAGS_verbose), absl::GetFlag(FLAGS_track_mrca),
               absl::GetFlag(FLAGS_neutral_mutations),
               absl::GetFlag(FLAGS_debug), absl::GetFlag(FLAGS_outpath),
               absl::GetFlag(FLAGS_file_name));

      // Perform the evolutionary run
      evo.Run();
    }
  }

  // Finish program timer
  end = clock();
  out0 << "Program time = " << (double)(end - start) / CLOCKS_PER_SEC
       << " seconds" << endl;
  out0.close();

  return 0;
}
