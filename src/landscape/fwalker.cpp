#include "fwalker.hpp"

FW::FW(int BASE, int GENOME_LENGTH, int ND, int TARGETS, int SOURCES, int ALL,
       int SWAPS, int DIMENSION, int NEUTRAL_MUTATIONS, int RUN, int THRESHOLD,
       int STOP_AT_FITTEST, int RECORD_PHENOTYPE_TRANSITIONS,
       int RESEED_FITNESS, int REMAP_PHENOTYPES, std::string FITNESS_ASSIGNMENT,
       std::string GENO_FNAME, std::string FILE_NAME, std::string OUTPATH,
       int VERBOSE) {
  print_component = 0;
  sparse = 1;
  sequence = 0;

  // Run index
  run = RUN;

  // Record times
  timing = 1;

  // Sequence and GP map properties
  L = GENOME_LENGTH;
  K = BASE;
  nd_ref = ND;

  // Amount of decorrelation and dimensionality allowed
  swaps = SWAPS;
  D = DIMENSION;
  neutral_mutations_ = NEUTRAL_MUTATIONS;

  // Navigability estimation settings
  all = ALL;
  sources = SOURCES;
  targets = TARGETS;
  threshold = THRESHOLD;
  stop_at_fittest_ = STOP_AT_FITTEST;
  reseed_fitness_ = RESEED_FITNESS;
  record_phenotype_transitions_ = RECORD_PHENOTYPE_TRANSITIONS;

  // Fitness assignment: fixed, target or source
  fitness_assignment = FITNESS_ASSIGNMENT;

  // File names
  geno_fname_ = GENO_FNAME;
  file_name_ = FILE_NAME;
  outpath_ = OUTPATH;

  // Remap phenotypes
  remap_phenotypes_ = REMAP_PHENOTYPES;

  // Verbosity for printing
  verbose_ = VERBOSE;

  // Initialize
  Init();
}

FW::~FW() { Clear(); }

void FW::Init() {
  // Counters for ruggedness estimation
  uphill_count = 0;
  downhill_count = 0;

  sequence = new int[L];

  int S = (int)pow((double)K, L);
  genos.reserve(S);

  // Make positions_ for randomising the order of mutations
  positions_.resize(L, 0);
  for (unsigned int i = 0; i < positions_.size(); i++) positions_[i] = i;

  // Load GP map, remap phenotypes such that form ordered list
  std::ifstream in;
  in.open(geno_fname_);
  std::string line;
  int max = 0;
  if (remap_phenotypes_) {
    std::set<int> totals;
    std::map<int, int> dict;
    // int count = 0;
    while (std::getline(in, line)) {
      int phenotype = String2Number(line);
      if (phenotype > max) max = phenotype;
      genos.push_back(phenotype);
      if (totals.find(phenotype) == totals.end()) {
        dict[phenotype] = totals.size();
        totals.insert(phenotype);
      }
    }
    in.close();
    totals.clear();

    for (unsigned int i = 0; i < genos.size(); i++) {
      genos[i] = dict[genos[i]];
    }
    int nd_ref = dict[nd_ref];
    if (verbose_) std::cout << "nd: " << nd_ref << endl;

    N_p = dict.size();
    dict.clear();
  } else {
    while (std::getline(in, line)) {
      int phenotype = String2Number(line);
      if (phenotype > max) max = phenotype;
      genos.push_back(phenotype);
    }
    in.close();
    N_p = max + 1;
  }

  // Make the fitness vector for the N_p phenotypes
  fitness.resize(N_p, -1);

  // Perform swaps of genotypes of phenotypes
  int scount = 0;
  int temp = -1;
  while (scount < swaps) {
    int a = rand() % S;
    int b = rand() % S;
    if (genos[a] != genos[b]) {
      temp = genos[a];
      genos[a] = genos[b];
      genos[b] = temp;
      scount++;
    }
  }

  // File name format:
  // <outpath>/config_<file_name>.txt : this is done in nav_disk.cpp
  // <outpath>/out_<file_name>.txt
  // <outpath>/stats_<file_name>.txt
  // <outpath>/times_<file_name>.txt
  // <outpath>/components/component_<file_name>_t-<target>_s-<source>.txt
  // <outpath>/phenotype_transitions/phenotype_transitions_<file_name>_t-<target>_s-<source>.txt
  // <outpath>/phenotype_transitions/fitness_<file_name>_t-<target>_s-<source>.txt
  char buffer[50];
  sprintf(buffer, "%sstats_%s.txt", outpath_.c_str(), file_name_.c_str());
  stats.open(buffer);
  if (timing == 1) {
    sprintf(buffer, "%stimes_%s.txt", outpath_.c_str(), file_name_.c_str());
    times.open(buffer);
  }
}

void FW::Clear() {
  if (sequence != 0) delete[] sequence;
}

void FW::fStartComponent(int t, int s) {
  char *String;
  String = new char[50];
  sprintf(String, "%scomponents/component_%s_t-%d_s-%d.txt", outpath_.c_str(),
          file_name_.c_str(), t, s);
  fp1.open(String);
  delete[] String;
}

void FW::fAddEdge(int x, int y) { fp1 << x << "\t" << y << "\n"; }

void FW::fEndComponent() { fp1.close(); }

void FW::Transfer(std::set<int> &from, std::set<int> &to, int x) {
  from.erase(x);
  to.insert(x);
}

void FW::Transfer(std::unordered_set<int> &from, std::unordered_set<int> &to,
                  int x) {
  from.erase(x);
  to.insert(x);
}

void FW::AssignFitness() {
  max_fit = 0;
  max_fit_p = -1;
  for (unsigned int i = 0; i < fitness.size(); i++) {
    if (i != nd_ref)
      fitness[i] = (double)rand() / RAND_MAX;
    else
      fitness[i] = 0.0;

    if (fitness[i] > max_fit) {
      max_fit = fitness[i];
      max_fit_p = i;
    }
  }
}

void FW::AssignFitness(int max) {
  max_fit = 0;
  max_fit_p = -1;
  for (unsigned int i = 0; i < fitness.size(); i++) {
    if (i != nd_ref)
      if (i == max)
        fitness[i] = 1.0;
      else
        fitness[i] = (double)rand() / RAND_MAX;
    else
      fitness[i] = 0.0;

    if (fitness[i] > max_fit) {
      max_fit = fitness[i];
      max_fit_p = i;
    }
  }
}

void FW::SetFitness() {
  if (all == 1)
    if (target != nd_ref)
      AssignFitness(target);
    else
      return;
  else
    AssignFitness();
}

void FW::SingleRun(int geno) {
  // Clear each of the sets used for set information
  u.clear();
  v.clear();
  p.clear();

  uphill_count = 0;
  downhill_count = 0;

  if (record_phenotype_transitions_) {
    phenotype_transitions.clear();
    phenotype_transitions.resize(N_p, std::vector<int>(N_p, 0));
  }

  // Shuffle the positions_ vector
  std::random_shuffle(positions_.begin(), positions_.end());
  // std::shuffle(positions_.begin(), positions_.end());

  if (source == -1 and target == -1)
    print_component = 1;
  else
    print_component = 0;

  if (print_component == 1) fStartComponent(target, source);

  // Insert source genotype
  u.insert(geno);

  // If u is empty, then no more to search OR if u or v are too big OR fittest
  // found, exit loop
  while (u.size() > 0) {
    if ((u.size() + v.size()) > threshold or
        (p.find(max_fit_p) != p.end() and stop_at_fittest_ == 1))
      break;

    // Take first element of u and attempt all 1-mutants.
    int g = *u.begin();
    int g_p = genos[g];

    int_to_basek(g, sequence, L, K);

    // Test every 1-mutant of genotype "trial"
    int all_downhill = 1;
    for (int j = 0; j < D; j++) {
      int position = positions_[j];
      int save = sequence[position];
      for (int k = 0; k < K; k++) {
        if (save == k) continue;
        sequence[position] = k;

        int n = basek_to_int(sequence, L, K);
        int n_p = genos[n];

        // std::cout << "neighbour: " << n << endl;
        // std::cout << "neighbour p: " << n_p << endl;

        // If fitness n_p is < fitness of g_p or if neutral mutations not
        // allowed and phenotype is the same, continue

        if ((fitness[n_p] < fitness[g_p]) or
            ((neutral_mutations_ != 1 and n_p == g_p))) {
          continue;
        }

        p.insert(n_p);

        if (print_component == 1) fAddEdge(g, n);

        // Added above before v set check
        all_downhill = 0;

        // If already in v, no need to consider
        if (v.find(n) != v.end()) continue;

        // If not already checked then put in u.
        u.insert(n);

        // Count the transtion
        if (record_phenotype_transitions_) {
          phenotype_transitions[g_p][n_p]++;
        }
      }
      // Make jth postion the original base again.
      sequence[position] = save;
    }
    // Move g from u to v
    Transfer(u, v, g);

    // Record whether local optima and not local optima
    uphill_count += 1 - all_downhill;
    downhill_count += all_downhill;
  }
  if (print_component == 1) fEndComponent();

  // Print file
  stats << target << "\t" << source << "\t" << max_fit_p << "\t" << max_fit
        << "\t" << geno << "\t" << genos[geno] << "\t" << p.size() << "\t"
        << (p.find(max_fit_p) != p.end()) << "\t" << u.size() << "\t"
        << v.size() << "\t" << uphill_count << "\t" << downhill_count
        << std::endl;

  // Record phenotype transitions
  if (record_phenotype_transitions_) {
    std::ofstream pts;
    char *pts_file;
    pts_file = new char[100];

    // Phenotype transitions for file_name_ at target and source
    sprintf(pts_file,
            "%sphenotype_transitions/"
            "phenotype_transitions_%s_t-%d_s-%d.txt",
            outpath_.c_str(), file_name_.c_str(), target, source);
    pts.open(pts_file);
    delete[] pts_file;

    for (unsigned i = 0; i < phenotype_transitions.size(); i++) {
      for (unsigned j = 0; j < phenotype_transitions[i].size(); j++) {
        if (j < phenotype_transitions[i].size() - 1) {
          pts << phenotype_transitions[i][j] << "\t";
        } else {
          pts << phenotype_transitions[i][j] << "\n";
        }
      }
    }

    pts.close();

    // Fitness values for file_name_ at target and source
    pts_file = new char[100];
    sprintf(pts_file,
            "%sphenotype_transitions/"
            "fitness_%s_t-%d_s-%d.txt",
            outpath_.c_str(), file_name_.c_str(), target, source);
    pts.open(pts_file);
    delete[] pts_file;

    for (unsigned i = 0; i < fitness.size(); i++) {
      pts << fitness[i] << "\n";
    }

    pts.close();
  }
}

void FW::Run() {
  stats << "Target"
        << "\t"
        << "Source"
        << "\t"
        << "Fit_p"
        << "\t"
        << "Max_Fit"
        << "\t"
        << "Geno"
        << "\t"
        << "Geno_p"
        << "\t"
        << "Phenos_found"
        << "\t"
        << "Fittest_found?"
        << "\t"
        << "u_size"
        << "\t"
        << "v_size"
        << "\t"
        << "uphill_count"
        << "\t"
        << "downhill_count" << std::endl;

  // If fitness landscape "fixed" for all source-target pairs, set fitness
  if (fitness_assignment == "fixed") {
    SetFitness();
  }

  // If small number of phenotypes test all as max individually.
  if (all == 1) targets = N_p;
  target = 0;
  for (target = 0; target < targets; target++) {
    // This makes fitness same for every target independent of other factors
    if (reseed_fitness_) {
      srand(target + 1);
      rand();
    }

    // If target is nd_ref and testing all phenotypes
    if (target == nd_ref and all == 1) continue;

    // If assignment mechanism for every target, set fitness
    if (fitness_assignment == "target") {
      SetFitness();
    }

    if (verbose_) {
      std::cout << "Fitness:" << endl;
      for (unsigned int i = 0; i < fitness.size(); i++) {
        std::cout << i << "\t" << fitness[i] << endl;
      }
    }

    // Timing
    source = 0;
    for (source = 0; source < sources; source++) {
      if (timing == 1) start = clock();

      int seed_for_geno = target * sources + source + 1;

      // This makes fitness same for every test independent of other factors
      if (reseed_fitness_) {
        srand(seed_for_geno);
        rand();
      }
      int geno_p;
      do
        geno_p = rand() % N_p;
      while (geno_p == max_fit_p);

      // If assignment mechanism for every source, set fitness
      if (fitness_assignment == "target") {
        // This makes fitness same for every test independent of other factors
        if (reseed_fitness_) {
          srand(seed_for_geno);
          rand();
        }

        SetFitness();
      }

      // This makes fitness same for every test independent of other factors
      if (reseed_fitness_) {
        srand(seed_for_geno);
        rand();
      }
      int geno;
      do
        geno = rand() % (int)pow((double)K, L);
      while (genos[geno] != geno_p);

      if (timing == 1) middle = clock();

      SingleRun(geno);

      if (timing == 1) end = clock();

      if (timing == 1)
        times << target << "." << source << "\t"
              << "Geno time: " << (double)(middle - start) / CLOCKS_PER_SEC
              << "\t"
              << "Search time: " << (double)(end - middle) / CLOCKS_PER_SEC
              << "\t"
              << "Source time: " << (double)(end - start) / CLOCKS_PER_SEC
              << endl;
    }
  }
  stats.close();
  if (timing == 1) times.close();
}
