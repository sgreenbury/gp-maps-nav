#include "evo_rna.hpp"

void Evo_RNA::Init(int ID, int TEST, int POPULATION_SIZE, int GENERATIONS,
                   Individual SOURCE, std::string TARGET,
                   std::set<std::string>& PHENO_SET, double MU,
                   std::string FITNESS_MEASURE, double FITNESS_THRESHOLD,
                   int VERBOSE, int TRACK_MRCA, int NEUTRAL_MUTATIONS,
                   int DEBUG, const std::string OUTPATH,
                   const std::string FILE_NAME) {
  // Initialize basic class atttributes
  id = ID;
  test = TEST;
  source = SOURCE;
  target = TARGET;
  K = source.K;
  L = source.L;
  mu = MU;
  N = POPULATION_SIZE;
  current_population.resize(N, source);
  next_population.resize(N, source);
  generations = GENERATIONS;
  phenotype_fitnesses[target] = 1;
  phenotype_fitnesses[source.phenotype] = source.fitness;

  // Initialize dummy individual
  dummy = Individual(K, L);

  // Clear and copy PHENO_SET to pheno_set
  pheno_set.clear();
  for (auto& it : PHENO_SET) {
    pheno_set.insert(it);
  }

  // Choice around whether to make trivial structure unfit.
  phenotype_fitnesses[std::string(L, '.')] = 0.0;

  // Stats file name
  outpath = OUTPATH;
  file_name = FILE_NAME;

  statsfile = absl::StrFormat("%sstats_%s.txt", outpath, file_name);
  mrcafile = absl::StrFormat("%smrca_%s.txt", outpath, file_name);
  summaryfile = absl::StrFormat("%ssummary_%s.txt", outpath, file_name);

  fitnesses.resize(N, 0.0);
  complete = 0;
  previous_max = 0.0;
  decrease = 0;
  best_max = 0.0;

  // Threshold choice for valley decision
  threshold = FITNESS_THRESHOLD * N;

  // Random or hamming fitness
  fitness_measure = FITNESS_MEASURE;

  max_fit_more_than_threshold = 0;

  maj_max_fit = 0.0;
  previous_maj_max_fit = 0.0;

  // Verbose and track_mrca flags
  verbose = VERBOSE;
  track_mrca = TRACK_MRCA;

  // Flag for allowing setting neutral mutant fitness to zero
  neutral_mutations_ = NEUTRAL_MUTATIONS;

  // Debug mode
  debug = DEBUG;

  if (verbose) {
    std::cout << "Source: " << source.phenotype << std::endl;
    std::cout << "Target: " << target << std::endl;
  }

  // Write header of the stats file for convenience
  if (id == 0 and test == 0) {
    WriteStatsFileHeader();
    WriteSummaryFileHeader();
    if (track_mrca) {
      WriteStatsMRCAHeader();
    }
  }

  // Initialize history and MRCA objects
  if (track_mrca) {
    // Make history of size generations (rows) x population (columns)
    history.clear();
    history.resize(generations + 1, std::vector<int>(N, -1));

    // Init the vectors to track Individuals
    history_genos.clear();
    history_genos.resize(generations + 1,
                         std::vector<Individual>(N, Individual(K, L)));

    // Set first generation parents as 0
    for (int i = 0; i < N; i++) {
      history[0][i] = 0;
      history_genos[0][i].AssignGenotype(source.g_seq);
      CalculateFitness(history_genos[0][i]);
    }
    // Set to -2 so that -1 != -2 and we record first MRCA (seed genotype)
    current_mrca = std::pair<int, int>(-2, 0);
    previous_mrca = std::pair<int, int>(-2, 0);
  } else {
    current_mrca = std::pair<int, int>(-1, 0);
    previous_mrca = std::pair<int, int>(-1, 0);
  }

  // Initialize first population
  InitFirstPopulation();
}

void Evo_RNA::WriteStatsFileHeader() {
  file.open(statsfile.c_str(), std::ofstream::out);
  file << "Target"
       << "\t"
       << "Source"
       << "\t"
       << "Generation"
       << "\t"
       << "Average_fitness"
       << "\t"
       << "Max_fitness"
       << "\t"
       << "Proportion_max_fit"
       << "\t"
       << "More_than_threshold_max_fit"
       << "\t"
       << "Majority_max_fitness"
       << "\t"
       << "Previous_majority_max_fitness"
       << "\t"
       << "Decrease"
       << "\t"
       << "Best_max_fitness"
       << "\t"
       << "Phenos_found"
       << "\t"
       << "most_common_phenotype"
       << "\t"
       << "most_common_phenotype_count"
       << "\t"
       << "most_common_phenotype_fitness"
       << "\t"
       << "most_common_phenotype_redundant"
       << "\t"
       << "max_fit_phenotype"
       << "\t"
       << "max_fit_phenotype_count"
       << "\t"
       << "max_fit_phenotype_fitness"
       << "\t"
       << "max_fit_phenotype_redundant" << std::endl;
  file.close();
}

void Evo_RNA::WriteSummaryFileHeader() {
  summary_fp.open(summaryfile.c_str(), std::ofstream::out);
  summary_fp << "Target"
             << "\t"
             << "Source"
             << "\t"
             << "source_genotype"
             << "\t"
             << "source_phenotype"
             << "\t"
             << "target_phenotype" << std::endl;
  summary_fp.close();
}

void Evo_RNA::InitFirstPopulation() {
  for (unsigned int i = 0; i < N; i++) {
    current_population[i].AssignGenotype(source.genotype);
  }
}

void Evo_RNA::PopulationG2P() {
  for (unsigned int i = 0; i < N; i++) {
    current_population[i].G2P();
  }
}

double Evo_RNA::CompareStrings(std::string s1, std::string s2) {
  int sames = 0;
  for (unsigned int i = 0; i < s1.size(); i++) {
    if (s1[i] == s2[i]) sames++;
  }
  return (double)sames / s1.size();
}

void Evo_RNA::CalculateFitness(Individual& ind) {
  // Random
  if (fitness_measure == "random") {
    std::unordered_map<std::string, double>::iterator it =
        phenotype_fitnesses.find(ind.phenotype);
    if (it != phenotype_fitnesses.end()) {
      ind.fitness = it->second;
    } else {
      double f = (double)rand() / RAND_MAX;
      phenotype_fitnesses[ind.phenotype] = f;
      ind.fitness = f;
    }
  } else {
    // Hamming
    if (fitness_measure == "hamming") {
      if (ind.phenotype != std::string(L, '.')) {
        ind.fitness = CompareStrings(ind.phenotype, target);
        phenotype_fitnesses[ind.phenotype] = ind.fitness;
        return;
      }
      ind.fitness = 0.0;
    } else {
      if (fitness_measure == "random_db_only") {
        // Random but only for those that occur in pheno_fname/pheno_set
        std::unordered_map<std::string, double>::iterator it =
            phenotype_fitnesses.find(ind.phenotype);
        if (it != phenotype_fitnesses.end()) {
          ind.fitness = it->second;
        } else {
          double f = 0.0;
          // Only if exists in pheno_set loaded earlier from pheno_fname
          if (pheno_set.find(ind.phenotype) != pheno_set.end()) {
            f = (double)rand() / RAND_MAX;
          }
          phenotype_fitnesses[ind.phenotype] = f;
          ind.fitness = f;
        }
      } else {
        if (fitness_measure == "hamming_db_only") {
          if (ind.phenotype != std::string(L, '.')) {
            if (pheno_set.find(ind.phenotype) != pheno_set.end()) {
              ind.fitness = CompareStrings(ind.phenotype, target);
            } else {
              ind.fitness = 0.0;
            }
            phenotype_fitnesses[ind.phenotype] = ind.fitness;
            return;
          }
          ind.fitness = 0.0;
        } else {
          if (fitness_measure == "random_hamming") {
            std::unordered_map<std::string, double>::iterator it =
                phenotype_fitnesses.find(ind.phenotype);
            if (it != phenotype_fitnesses.end()) {
              ind.fitness = it->second;
            } else {
              // Get a random genotype to use instead of ind for the phenotype
              std::string trivial_structure = std::string(L, '.');
              do {
                dummy.RandomGenotype();
                dummy.G2P();
              } while (dummy.phenotype == trivial_structure or
                       dummy.phenotype == target);

              double f = CompareStrings(dummy.phenotype, target);
              phenotype_fitnesses[ind.phenotype] = f;
              ind.fitness = f;
            }
            phenotype_fitnesses[ind.phenotype] = ind.fitness;
            return;
          }
        }
      }
    }
  }
}

void Evo_RNA::CalculateFitnesses() {
  for (unsigned int i = 0; i < current_population.size(); i++) {
    CalculateFitness(current_population[i]);
    // If neutral mutations are not allowed, then set any cases where they have
    // occurred to have zero fitness
    if (neutral_mutations_ != 1) {
      if ((current_population[i].phenotype ==
           current_population[i].phenotype_prev) and
          (current_population[i].genotype !=
           current_population[i].genotype_prev)) {
        // Set the fitness
        current_population[i].fitness = 0.;
      }
    }
  }
}

void Evo_RNA::MakeFitnessVector() {
  double total = 0.0;
  for (int i = 0; i < N; i++) {
    total += current_population[i].fitness;
  }
  if (total < 0.00000001) {
    for (int i = 0; i < N; i++) {
      fitnesses[i] = (double)(i + 1) / N;
    }
  } else {
    for (int i = 0; i < N; i++) {
      if (i == 0) {
        fitnesses[i] = current_population[i].fitness / total;
        continue;
      }
      fitnesses[i] = fitnesses[i - 1] + current_population[i].fitness / total;
    }
  }
  fitnesses[N - 1] = 1.0;
}

int Evo_RNA::PickMember() {
  double r = (double)rand() / RAND_MAX;
  double lower = 0.0;
  double upper = 0.0;
  for (unsigned int i = 0; i < N; i++) {
    upper = fitnesses[i];
    if (r >= lower and r <= upper) {
      return i;
    }
    lower = fitnesses[i];
  }
  throw std::invalid_argument("Failed to pick a member of population.");
  return -1;
}

void Evo_RNA::PopulateNextGeneration() {
  int j;
  for (int i = 0; i < N; i++) {
    j = PickMember();

    if (track_mrca) {
      // Store the parent index of the next member:
      // E.g. member i at generatation t+1 has parent j from generation t
      history[generation + 1][i] = j;

      // Assign genotype to the history
      // TODO: convert this to a map storage dtype:
      // std::map <std::pair <int, int>, Individual> history_ind;
      history_genos[generation + 1]
                   [i].AssignGenotype(current_population[j].g_seq);

      history_genos[generation + 1][i].fitness = current_population[j].fitness;
    }

    // Store genotype in next_population Individual vector
    next_population[i].AssignGenotype(current_population[j].g_seq);

    // NEW: For checking against previous geno/phenotypes
    if (debug or neutral_mutations_ != 1) {
      current_population[j].genotype_prev = current_population[j].genotype;
      current_population[j].phenotype_prev = current_population[j].phenotype;
      next_population[i] = current_population[j];
    }

    // Point mutations
    for (unsigned int pos = 0; pos < next_population[i].genotype.size();
         pos++) {
      double r_mu = (double)rand() / RAND_MAX;
      if (r_mu < mu) {
        next_population[i].MutateBase(pos);
      }
    }
    next_population[i].G2P();
  }
}

std::pair<int, int> Evo_RNA::CalculateMRCA(
    int generation, std::vector<std::vector<int>>& history) {
  // Create two containers for performing
  std::unordered_set<int> previous_parents;
  std::unordered_set<int> current_parents;

  // Add each member of most current generation
  for (int i = 0; i < N; i++) {
    current_parents.insert(i);
  }

  // Loop back through all generations
  for (int i = generation; i >= 0; i--) {
    // Clear set that stores the unique parents of generation "i"
    previous_parents.clear();
    // Loop over current parents and add their parents
    for (auto& it : current_parents) {
      previous_parents.insert(history[i][it]);
    }
    // If exactly one parent left, MRCA is found
    if (previous_parents.size() == 1) {
      // MRCA found
      std::pair<int, int> mrca_pair =
          std::pair<int, int>(i - 1, *previous_parents.begin());
      if (vverbose) {
        std::cout << "Generation " << generation << ": " << mrca_pair.first
                  << "-" << mrca_pair.second << std::endl;
      }
      return std::pair<int, int>(i - 1, *previous_parents.begin());
    } else {
      // First generation, MRCA is original seed genotype
      if (generation == 0) {
        // MRCA is original seed genotype
        return std::pair<int, int>(-1, 0);
      }
      // MRCA not found yet, clear current parents, add previous parents
      else {
        // Clear current and copy previous parents to current
        current_parents.clear();
        current_parents.insert(previous_parents.begin(),
                               previous_parents.end());
      }
    }
  }
}

std::deque<std::pair<int, int>> Evo_RNA::AncestryPath(
    std::pair<int, int>& start, std::pair<int, int>& end) {
  // Store the ancestors in a deque as get path backwards so need push_front()
  std::deque<std::pair<int, int>> ancestry_path;

  // If the end mrca is before the beginning, return end only
  if (end.first < 0) {
    ancestry_path.push_front(end);
    return ancestry_path;
  }

  // Begin at end
  std::pair<int, int> parent = end;

  // Loop through parents until start
  for (int i = end.first; i > start.first; i--) {
    ancestry_path.push_front(parent);
    parent.second = history[parent.first][parent.second];
    parent.first = i - 1;
  }

  // Return deque of path
  return ancestry_path;
}

void Evo_RNA::AssignNextGeneration2Current() {
  for (int i = 0; i < N; i++) {
    current_population[i].AssignGenotype(next_population[i].g_seq);

    // NEW: For checking changes to genotype
    if (debug or neutral_mutations_ == 0) {
      current_population[i] = next_population[i];
    }
  }
}

void Evo_RNA::AssignMRCA(std::pair<int, int>& mrca) {
  previous_mrca = current_mrca;
  current_mrca = mrca;
}

void Evo_RNA::WriteStatsMRCAHeader() {
  // Write header for the MRCA file, use main stats file as template
  mrca_fp.open(mrcafile.c_str(), std::ofstream::out);
  mrca_fp << "Target"
          << "\t"
          << "Source"
          << "\t"
          << "Generation"
          << "\t"
          << "mrca_generation"
          << "\t"
          << "mrca_member"
          << "\t"
          << "mrca_genotype"
          << "\t"
          << "mrca_phenotype"
          << "\t"
          << "mrca_fitness" << std::endl;
  mrca_fp.close();
}
void Evo_RNA::WriteStatsMRCA(std::pair<int, int>& mrca) {
  // Write out stats for MRCA
  mrca_fp.precision(default_precision);

  if (verbose) {
    std::cout << "MRCA @ generation: " << absl::StrFormat("%4.0f", mrca.first)
              << ", mem: " << absl::StrFormat("%4.0f", mrca.second)
              << ", geno: "
              << history_genos[mrca.first + 1][mrca.second].genotype
              << ", fit: "
              << absl::StrFormat(
                     "%.6f", history_genos[mrca.first + 1][mrca.second].fitness)
              << std::endl;
  }

  // Get ancestry path of the previous_mrca to current_mrca
  std::deque<std::pair<int, int>> ancestry_path =
      AncestryPath(previous_mrca, current_mrca);

  if (vverbose) {
    std::cout << "Previou MRCA: " << previous_mrca.first << "-"
              << previous_mrca.second << std::endl;
    std::cout << "Current MRCA: " << current_mrca.first
              << "::::::" << current_mrca.second << std::endl;
    for (unsigned int i = 0; i < ancestry_path.size(); i++) {
      std::cout << "Ancestor: " << ancestry_path[i].first << "\t"
                << ancestry_path[i].second << std::endl;
    }
  }
  for (unsigned int i = 0; i < ancestry_path.size(); i++) {
    // Write the individual to mrca file
    std::pair<int, int> ind_on_ancestry_path = ancestry_path[i];
    mrca_fp << id << "\t" << test << "\t";

    // If the last element, write the first generation that the ancestor is
    // the MRCA
    if (i == ancestry_path.size() - 1)
      mrca_fp << generation << "\t";
    else
      mrca_fp << "NaN"
              << "\t";

    mrca_fp << ind_on_ancestry_path.first << "\t" << ind_on_ancestry_path.second
            << "\t"
            << history_genos[ind_on_ancestry_path.first + 1]
                            [ind_on_ancestry_path.second]
                                .genotype
            << "\t"
            << history_genos[ind_on_ancestry_path.first + 1]
                            [ind_on_ancestry_path.second]
                                .phenotype
            << "\t";
    mrca_fp << std::setprecision(fitness_precision)
            << history_genos[ind_on_ancestry_path.first + 1]
                            [ind_on_ancestry_path.second]
                                .fitness
            << std::setprecision(default_precision) << std::endl;
  }
}

void Evo_RNA::OneGeneration() {
  // Probably not needed as any copy to of G_1 to G_2 has a G2P call.
  PopulationG2P();
  CalculateFitnesses();
  WriteStats();
  if (generation % ((int)((double)generations / 5)) == 0 and
      (vverbose or debug)) {
    std::cout << generation << std::endl;
    for (int i = 0; i < N; i++) {
      std::cout << current_population[i].phenotype << "\t"
                << current_population[i].fitness << std::endl;
    }
  }

  MakeFitnessVector();

  // Check the G, P, fitnesses are correct
  if (debug) {
    std::cout << "Generation: " << generation << "\n";
    std::cout << "Current genotype\tPrevious genotype\tCurrent "
                 "phenotype\tPrevious phenotype\n";
    for (unsigned int i = 0; i < N; i++) {
      int same_geno = (current_population[i].genotype ==
                       current_population[i].genotype_prev);
      int same_pheno = (current_population[i].phenotype ==
                        current_population[i].phenotype_prev);
      std::cout << current_population[i].genotype << "\t"
                << current_population[i].genotype_prev << "\t"
                << "Same: " << same_geno << "\t|\t"
                << current_population[i].phenotype << "\t"
                << current_population[i].phenotype_prev << "\t"
                << "Same: " << same_pheno << "\t"
                << "Fitness: " << current_population[i].fitness << "\t"
                << "CDF: " << fitnesses[i] << "\t";
      if ((same_geno == 0) & (same_pheno == 1)) {
        std::cout << "NEUTRAL MUTANT!\n";
      } else {
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
    std::cin.get();
  }

  PopulateNextGeneration();

  if (track_mrca) {
    // Get MRCA
    std::pair<int, int> mrca = CalculateMRCA(generation, history);
    // Assign MRCA gneration
    AssignMRCA(mrca);
    // Write stats if different to previous generation
    if (current_mrca.first != previous_mrca.first) {
      WriteStatsMRCA(mrca);
    }
  }

  AssignNextGeneration2Current();
}

void Evo_RNA::WriteSummary() {
  summary_fp.precision(default_precision);
  summary_fp << id << "\t" << test << "\t" << source.genotype << "\t"
             << source.phenotype << "\t" << target << std::endl;
}

void Evo_RNA::WriteStats() {
  // Should make these into Mean and SE functions for simplicity
  // evo_run - generation - source - source phenotype - target - av_fit -
  // max_fit - proportion with max_fit - above threshold? - decrease - best max
  // - phenotypes
  double av_fit = 0.0;
  double max_fit = 0.0;

  // Calculate the average and max fitness
  for (int i = 0; i < N; i++) {
    av_fit += current_population[i].fitness / N;
    if (current_population[i].fitness > max_fit)
      max_fit = current_population[i].fitness;
  }

  // Count the number of max fit
  int count = 0;
  for (int i = 0; i < N; i++) {
    if (current_population[i].fitness >= max_fit - 0.00000001) count++;
  }

  // Max fit occupy more than half?
  max_fit_more_than_threshold = (count > threshold) ? 1 : 0;

  if (max_fit_more_than_threshold == 1) {
    previous_maj_max_fit = maj_max_fit;
    maj_max_fit = max_fit;
  }

  // Get the largest maximum fitness
  best_max = std::max(best_max, max_fit);

  // Has the majority max fitness decreased?
  if (maj_max_fit < previous_maj_max_fit)
    decrease = 1;
  else
    decrease = 0;
  previous_max = max_fit;

  // Count the number of each phenotype
  phenotype_counter.clear();
  for (int i = 0; i < N; i++) {
    std::unordered_map<std::string, int>::iterator it(
        phenotype_counter.find(current_population[i].phenotype));
    if (it != phenotype_counter.end()) {
      it->second++;
    } else {
      phenotype_counter[current_population[i].phenotype] = 1;
    }
  }

  // Record details about most common phenotype
  int most_common_phenotype_count = -1;
  int redundant_most_common_phenotype = 0;
  double most_common_phenotype_fitness = 0.0;
  std::string most_common_phenotype = "";
  for (auto& it : phenotype_counter) {
    if (it.second > most_common_phenotype_count) {
      most_common_phenotype = it.first;
      redundant_most_common_phenotype = 0;
      most_common_phenotype_count = it.second;
      most_common_phenotype_fitness =
          phenotype_fitnesses[most_common_phenotype];
    } else {
      if (it.second == most_common_phenotype_count) {
        redundant_most_common_phenotype = 1;
      }
    }
  }

  // Calculate the average and max fitness
  double max_fit_phenotype_fitness = -1.0;
  int max_fit_phenotype_count = 0;
  std::string max_fit_phenotype = "";
  int redundant_max_fit_phenotype = 0;
  for (int i = 0; i < N; i++) {
    if (current_population[i].fitness > max_fit_phenotype_fitness) {
      max_fit_phenotype_fitness = current_population[i].fitness;
      max_fit_phenotype = current_population[i].phenotype;
    }
  }

  // Count the number of phenotypes with max fitness and record if more than a
  // single max fit phenotype
  for (int i = 0; i < N; i++) {
    if (current_population[i].fitness >= max_fit - 0.00000001) {
      max_fit_phenotype_count++;
      if (current_population[i].phenotype != max_fit_phenotype) {
        redundant_max_fit_phenotype = 1;
      }
    };
  }

  // Write results to file, longer precision for fitness to ensure accurate max
  // is recorded when at target
  file.precision(default_precision);
  file << id << "\t" << test << "\t" << generation << "\t";

  // Write fitness information
  file << std::setprecision(fitness_precision) << av_fit << "\t" << max_fit
       << "\t" << std::setprecision(default_precision) << (double)count / N
       << "\t" << max_fit_more_than_threshold << "\t";

  if (max_fit_more_than_threshold) {
    file << std::setprecision(fitness_precision) << maj_max_fit << "\t"
         << std::setprecision(default_precision);
  } else {
    file << "NaN"
         << "\t";
  }
  file << std::setprecision(fitness_precision) << previous_maj_max_fit << "\t"
       << std::setprecision(default_precision) << decrease << "\t"
       << std::setprecision(fitness_precision) << best_max << "\t"
       << std::setprecision(default_precision) << phenotype_fitnesses.size()
       << "\t";

  // Extra information to understand population better
  file << most_common_phenotype << "\t" << most_common_phenotype_count << "\t"
       << std::setprecision(fitness_precision) << most_common_phenotype_fitness
       << "\t" << std::setprecision(default_precision)
       << redundant_most_common_phenotype << "\t" << max_fit_phenotype << "\t"
       << max_fit_phenotype_count << "\t"
       << std::setprecision(fitness_precision) << max_fit_phenotype_fitness
       << "\t" << std::setprecision(default_precision)
       << redundant_max_fit_phenotype << std::endl;
}

void Evo_RNA::Run() {
  // Capture stuff to do with current generation...
  generation = 0;

  // Open stats file
  file.open(statsfile.c_str(), std::ofstream::out | std::ofstream::app);
  summary_fp.open(summaryfile.c_str(), std::ofstream::out | std::ofstream::app);

  // Open MRCA file if tracking
  if (track_mrca) {
    mrca_fp.open(mrcafile.c_str(), std::ofstream::out | std::ofstream::app);
  }

  // Write to summay file
  WriteSummary();

  // Loop over generations
  for (generation = 0; generation < generations; generation++) {
    OneGeneration();
  }

  // Close files
  file.close();
  summary_fp.close();

  if (track_mrca) {
    mrca_fp.close();
  }
}
