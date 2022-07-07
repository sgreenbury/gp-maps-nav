#include "fwalker_frna.hpp"
#include "rna_fold.hpp"

double CompareStringsLevDist(const std::string &s1, const std::string &s2) {
  int max_dist = std::max(s1.size(), s2.size());
  int lev_dist = GenLevDist(s1, s2);
  return 1. - (double)lev_dist / max_dist;
}

FW_frna::FW_frna(int BASE, int GENOME_LENGTH, int ND, int TARGETS, int SOURCES,
                 int ALL, int NEUTRAL_MUTATIONS, int THRESHOLD,
                 std::string OUTPATH, std::string FILE_NAME,
                 std::string GENO_FNAME, int LOAD_GENOS,
                 std::string PHENO_FNAME) {
  print_component = 0;
  sparse = 1;

  TIMING = 1;
  all = ALL;
  L = GENOME_LENGTH;
  K = BASE;
  nd_ref = ND;
  sources = SOURCES;
  targets = TARGETS;
  fit_prob = 0.5;

  threshold = THRESHOLD;

  neutral_mutations_ = NEUTRAL_MUTATIONS;

  outpath_ = OUTPATH;
  file_name_ = FILE_NAME;

  geno_fname = GENO_FNAME;
  pheno_fname = PHENO_FNAME;

  load_genos_ = LOAD_GENOS;

  debug = 0;

  Init();
}

FW_frna::FW_frna(int BASE, int GENOME_LENGTH, int ND, int TARGETS, int SOURCES,
                 int ALL, int NEUTRAL_MUTATIONS, int THRESHOLD,
                 int POPULATION_SIZE, int SHAPE_LEVEL,
                 std::string FITNESS_MEASURE, std::string OUTPATH,
                 std::string FILE_NAME, std::string GENO_FNAME, int LOAD_GENOS,
                 std::string PHENO_FNAME, int RECORD_DIFFS_ONLY, int DEBUG) {
  print_component = 0;
  sparse = 1;

  TIMING = 1;
  all = ALL;
  L = GENOME_LENGTH;
  K = BASE;
  nd_ref = ND;
  sources = SOURCES;
  targets = TARGETS;
  fit_prob = 0.5;

  shape_level = SHAPE_LEVEL;

  threshold = THRESHOLD;

  neutral_mutations_ = NEUTRAL_MUTATIONS;

  outpath_ = OUTPATH;
  file_name_ = FILE_NAME;

  geno_fname = GENO_FNAME;
  pheno_fname = PHENO_FNAME;

  population_size = POPULATION_SIZE;
  fitness_measure = FITNESS_MEASURE;

  load_genos_ = LOAD_GENOS;
  record_diffs_only = RECORD_DIFFS_ONLY;

  debug = DEBUG;

  Init();
}

FW_frna::~FW_frna() { Clear(); }

void FW_frna::SetAlphabet() {
  alphabet.resize(K);
  if (K == 4) {
    alphabet[0] = 'C';
    alphabet[1] = 'G';
    alphabet[2] = 'A';
    alphabet[3] = 'U';
  }
}

void FW_frna::Init() {
  SetAlphabet();

  std::ifstream in_g, in_p;
  in_g.open(geno_fname.c_str());
  in_p.open(pheno_fname.c_str());
  std::string line_p;
  std::string line_g;
  p.clear();
  pheno_set.clear();
  std::string nd(L, '.');
  while (std::getline(in_p, line_p)) {
    if (p.find(line_p) == p.end() and line_p != nd) {
      phenos_start.push_back(line_p);
      p.insert(line_p);
      pheno_set.insert(line_p);

      if (load_genos_) {
        std::getline(in_g, line_g);
        genos_start.push_back(line_g);
      }
    } else {
      if (load_genos_) {
        std::getline(in_g, line_g);
      }
    }
  }
  in_p.close();
  in_g.close();
  p.clear();
  N_p = phenos_start.size();

  std::cout << "Phenotypes: " << N_p << std::endl;

  trivial_shape = shape(nd, shape_level);

  std::string statsfile =
      absl::StrFormat("%sstats_%s.txt", outpath_, file_name_);
  std::string timesfile =
      absl::StrFormat("%stimes_%s.txt", outpath_, file_name_);
  stats.open(statsfile.c_str());
  if (TIMING == 1) {
    times.open(timesfile.c_str());
  }
}

void FW_frna::Clear() {}

void FW_frna::fStartComponent(int p, int n) {
  char *String;
  String = new char[50];
  sprintf(String, "components/component%d-%d.txt", p, n);
  fp1.open(String);
  delete[] String;
}

void FW_frna::fAddEdge(std::string x, std::string y) {
  fp1 << x << "\t" << y << "\n";
}

void FW_frna::fEndComponent() { fp1.close(); }

void FW_frna::Transfer(std::set<int> &from, std::set<int> &to, int x) {
  from.erase(x);
  to.insert(x);
}

void FW_frna::Transfer(std::unordered_set<int> &from,
                       std::unordered_set<int> &to, int x) {
  from.erase(x);
  to.insert(x);
}

void FW_frna::Transfer(std::set<std::string> &from, std::set<std::string> &to,
                       std::string x) {
  from.erase(x);
  to.insert(x);
}

void FW_frna::Transfer(std::unordered_set<std::string> &from,
                       std::unordered_set<std::string> &to, std::string x) {
  from.erase(x);
  to.insert(x);
}

double FW_frna::CompareStrings(std::string s1, std::string s2) {
  int sames = 0;
  for (unsigned int i = 0; i < s1.size(); i++) {
    if (s1[i] == s2[i]) sames++;
  }
  return (double)sames / s1.size();
}

double FW_frna::SetUpFitness(std::string g_p) {
  double f_g = -1;
  if (fitness.find(g_p) != fitness.end()) {
    f_g = fitness[g_p];
  } else {
    if (fitness_measure == "random") {
      f_g = (double)rand() / RAND_MAX;
      fitness[g_p] = f_g;
    } else if (fitness_measure == "hamming") {
      if (g_p != std::string(L, '.')) {
        f_g = CompareStrings(g_p, phenos_start[target]);
        fitness[g_p] = f_g;
      }
    } else if (fitness_measure == "lev") {
      if (g_p != trivial_shape) {
        f_g = CompareStringsLevDist(g_p,
                                    shape(phenos_start[target], shape_level));
        fitness[g_p] = f_g;
      } else {
        f_g = 0.;
        fitness[g_p] = f_g;
      }
    } else if (fitness_measure == "random_db_only") {
      if (pheno_set.find(g_p) != pheno_set.end()) {
        f_g = (double)rand() / RAND_MAX;
      } else {
        f_g = 0.0;
      }

      fitness[g_p] = f_g;
    } else if (fitness_measure == "hamming_db_only") {
      if (pheno_set.find(g_p) != pheno_set.end()) {
        if (g_p != std::string(L, '.')) {
          f_g = CompareStrings(g_p, phenos_start[target]);
          fitness[g_p] = f_g;
        }
      } else {
        f_g = 0.0;
      }

      fitness[g_p] = f_g;
    }
  }
  return f_g;
}

std::string FW_frna::FindGenotype(std::string phys_seq_string) {
  // Initialize the required parts
  double found;
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  // Copy target structure to phys_seq
  std::strcpy(phys_seq, phys_seq_string.c_str());

  // Try to find a sequence, uses random char_seq
  do {
    for (int i = 0; i < L; i++) {
      int ri = rand() % K;
      char_seq[i] = alphabet[ri];
    }
    found = InverseFold(char_seq, phys_seq);
  } while (found != 0);

  std::string char_seq_string(L, '.');
  char_seq_string.assign(char_seq);

  // Free memory
  free(char_seq);
  free(phys_seq);

  return char_seq_string;
}

void FW_frna::SingleRun(std::string geno) {
  // Clear each of the sets used for set information
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  u.clear();
  v.clear();
  p.clear();
  genos.clear();
  fitness.clear();

  fitness[std::string(L, '.')] = 0.0;
  fitness[phenos_start[target]] = 1.0;

  fitness[phenos_start[source]] = (double)rand() / RAND_MAX;

  double f_g = 0.0;
  double f_n = 0.0;

  std::string g;
  g.resize(L);
  std::string g_p;
  g_p.resize(L);
  std::string n;
  n.resize(L);
  std::string n_p;
  n_p.resize(L);

  u.insert(geno);

  std::strcpy(char_seq, geno.c_str());
  char_seq[L] = '\0';
  Fold(char_seq, phys_seq);
  n_p.assign(phys_seq);
  genos[geno] = std::string(phys_seq);

  // If u is empty, then run aground.
  // If u or v are too big, then
  if (source_index == -1 and target_index == -1)
    print_component = 1;
  else
    print_component = 0;

  if (print_component == 1) fStartComponent(target_index, source_index);

  while (u.size() > 0) {
    if ((u.size() + v.size()) > threshold or
        p.find(phenos_start[target]) != p.end())
      break;

    // Take first element of u and attempt all 1-mutants.

    // Copy beginning string of u to g.
    g = *u.begin();
    g_p = genos[g];

    if (fitness.find(g_p) != fitness.end()) {
      f_g = fitness[g_p];
    } else {
      f_g = (double)rand() / RAND_MAX;
      fitness[g_p] = f_g;
    }

    n = g;

    // Target every 1-mutant of genotype "trial"
    for (int j = 0; j < L; j++) {
      char save = n[j];
      for (int k = 0; k < K; k++) {
        if (save == alphabet[k]) continue;
        n[j] = alphabet[k];
        if (genos.find(n) != genos.end()) {
          n_p = genos[n];
        } else {
          // Convert string into phenotype using GP program.
          std::strcpy(char_seq, n.c_str());
          char_seq[L] = '\0';
          Fold(char_seq, phys_seq);
          n_p.assign(phys_seq);
          genos[n] = n_p;
        }

        if (fitness.find(n_p) != fitness.end()) {
          f_n = fitness[n_p];
        } else {
          f_n = (double)rand() / RAND_MAX;
          fitness[n_p] = f_n;
        }

        // If fitness n_p is < fitness of g_p or if neutral mutations not
        // allowed and phenotype is the same, continue
        if ((f_n < f_g) or ((neutral_mutations_ != 1) and n_p == g_p)) {
          continue;
        }

        p.insert(n_p);

        if (print_component == 1) fAddEdge(g, n);

        // If already in v, no need to consider
        if (v.find(n) != v.end()) continue;

        // If not already checked then put in u.
        u.insert(n);
      }
      // Make jth postion the original base again.
      n[j] = save;
    }
    Transfer(u, v, g);
  }
  free(char_seq);
  free(phys_seq);

  if (print_component == 1) fEndComponent();

  // Print file
  stats << target_index << "\t" << source_index << "\t" << phenos_start[source]
        << "\t" << fitness[phenos_start[source]] << "\t" << phenos_start[target]
        << "\t" << fitness[phenos_start[target]] << "\t" << p.size() << "\t"
        << (p.find(phenos_start[target]) != p.end()) << "\t" << u.size() << "\t"
        << v.size() << std::endl;
}

void FW_frna::SingleRunDFS(std::string geno) {
  // Clear each of the sets used for set information
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  u.clear();
  v.clear();
  p.clear();
  u2.clear();
  genos.clear();
  fitness.clear();

  fitness[std::string(L, '.')] = 0.0;
  fitness[phenos_start[target]] = 1.0;
  SetUpFitness(phenos_start[source]);

  double f_g = 0.0;
  double f_n = 0.0;

  std::string g;
  g.resize(L);
  std::string g_p;
  g_p.resize(L);
  std::string n;
  n.resize(L);
  std::string n_p;
  n_p.resize(L);

  u.insert(geno);
  u2.push_back(geno);

  std::strcpy(char_seq, geno.c_str());
  char_seq[L] = '\0';
  Fold(char_seq, phys_seq);
  n_p.assign(phys_seq);
  genos[geno] = std::string(phys_seq);

  // If u is empty, then run aground.
  // If u or v are too big, then
  if (source_index == -1 and target_index == -1)
    print_component = 1;
  else
    print_component = 0;

  if (print_component == 1) fStartComponent(target_index, source_index);

  // Make two new vectors for order searching and shuffle every it
  std::vector<int> Ls(L);
  for (unsigned int j = 0; j < Ls.size(); j++) Ls[j] = j;
  std::vector<int> Ks(K);
  for (unsigned int j = 0; j < Ks.size(); j++) Ks[j] = j;

  while (u2.size() > 0) {
    int CHECKED = 1;
    if ((u.size() + v.size()) > threshold or
        p.find(phenos_start[target]) != p.end())
      break;

    // Shuffle Ls and Ks
    std::random_shuffle(Ls.begin(), Ls.end());
    std::random_shuffle(Ks.begin(), Ks.end());

    g = u2.back();
    if (genos.find(g) != genos.end())
      g_p = genos[g];
    else {
      std::strcpy(char_seq, g.c_str());
      char_seq[L] = '\0';
      Fold(char_seq, phys_seq);
      g_p.assign(phys_seq);
      if (genos.size() < threshold) genos[g] = g_p;
    }

    // Use fitness function to get fitness
    f_g = SetUpFitness(g_p);

    n = g;

    // Test every 1-mutant of genotype "trial"
    std::string aneut = "";
    for (unsigned int j = 0; j < Ls.size(); j++) {
      char save = n[Ls[j]];
      for (unsigned int k = 0; k < Ks.size(); k++) {
        if (save == alphabet[Ks[k]]) continue;

        n[Ls[j]] = alphabet[Ks[k]];
        if (genos.find(n) != genos.end()) {
          n_p = genos[n];
        } else {
          // Convert string into phenotype using GP program.
          std::strcpy(char_seq, n.c_str());
          char_seq[L] = '\0';
          Fold(char_seq, phys_seq);
          n_p.assign(phys_seq);
          if (genos.size() < threshold) genos[n] = n_p;
        }

        // Use fitness function to get fitness
        f_n = SetUpFitness(n_p);

        if (f_n > f_g and 0) {
          std::cout << "Ks[k]      : " << Ks[k] << std::endl;
          std::cout << "Ls[j]      : " << Ls[j] << std::endl;
          std::cout << "g          : " << g << std::endl;
          std::cout << "neighbour  : " << n << std::endl;
          std::cout << "neighbour p: " << n_p << std::endl;
          std::cout << "fitness p  : " << f_n << std::endl;
          std::cout << "u length   : " << u.size() << std::endl;
          std::cout << "v length   : " << v.size() << std::endl;
          std::cout << "genos size : " << genos.size() << std::endl;
          std::cout << "p length   : " << p.size() << std::endl;
        }
        // If fitness n_p is < g_p, then add phenotype to p but don't test
        // genotype.
        if ((f_n < f_g) or ((neutral_mutations_ != 1) and n_p == g_p)) {
          continue;
        }

        // If only single genotype left in u, write down a neutral option
        if (u.size() == 1)
          if (g_p == n_p and v.find(n) == v.end() and u.find(n) == u.end())
            aneut = n;

        p.insert(n_p);

        if (print_component == 1) fAddEdge(g, n);

        // If already in v or u, no need to consider
        if (v.find(n) != v.end() or u.find(n) != u.end()) continue;

        // If not already checked and non-neutral, then put in u.
        if (g_p != n_p) {
          // Add to back of u2 where it will be picked up in next iteration
          u2.push_back(n);
          u.insert(n);
          // Set CHECKED to 0, exit the loop to begin next genotype
          CHECKED = 0;
          break;
        }
      }
      if (CHECKED == 0) break;

      // Make jth postion the original base again.
      n[Ls[j]] = save;
    }
    if (CHECKED == 1) {
      v.insert(g);
      u.erase(g);
      // Assert check that "g" is really at the back
      assert(u2.back() == g);
      u2.pop_back();

      if (u.size() == 0)
        if (aneut.size() != 0) {
          u2.push_back(aneut);
          u.insert(aneut);
          CHECKED = 0;
        }
    }
  }
  free(char_seq);
  free(phys_seq);

  if (print_component == 1) fEndComponent();

  // Print file
  stats << target_index << "\t" << source_index << "\t" << phenos_start[source]
        << "\t" << fitness[phenos_start[source]] << "\t" << phenos_start[target]
        << "\t" << fitness[phenos_start[target]] << "\t" << p.size() << "\t"
        << (p.find(phenos_start[target]) != p.end()) << "\t" << u2.size()
        << "\t" << v.size() << std::endl;
}

void FW_frna::Run() {
  stats << "Target_index"
        << "\t"
        << "Source_index"
        << "\t"
        << "Source"
        << "\t"
        << "Fitness_source"
        << "\t"
        << "Target"
        << "\t"
        << "Fitness_target"
        << "\t"
        << "Phenos_found"
        << "\t"
        << "Fittest_found?"
        << "\t"
        << "u_size"
        << "\t"
        << "v_size" << std::endl;

  target_index = 0;
  for (target_index = 0; target_index < targets; target_index++) {
    std::cout << targets << std::endl;
    // Choose which fRNA phenotype is target and which is source.
    target = rand() % phenos_start.size();

    // Timing
    std::string geno;
    source_index = 0;
    for (source_index = 0; source_index < sources; source_index++) {
      std::cout << target_index << "." << source_index << std::endl;
      if (TIMING == 1) start = clock();

      do {
        source = rand() % phenos_start.size();
      } while (source == target);

      // If load_genos_, then get source from vector
      if (load_genos_) {
        geno = genos_start[source];
      }
      // Else find a genotype with inverse_fold
      else {
        geno = FindGenotype(phenos_start[source]);
      }

      if (TIMING == 1) middle = clock();

      SingleRun(geno);

      if (TIMING == 1) end = clock();

      if (TIMING == 1)
        times << target_index << "." << source_index << "\t"
              << "Geno time: " << (double)(middle - start) / CLOCKS_PER_SEC
              << "\t"
              << "Search time: " << (double)(end - middle) / CLOCKS_PER_SEC
              << "\t"
              << "Source time: " << (double)(end - start) / CLOCKS_PER_SEC
              << std::endl;
    }
  }
  stats.close();
  if (TIMING == 1) times.close();
}

void FW_frna::RunDFS() {
  stats << "Target_index"
        << "\t"
        << "Source_index"
        << "\t"
        << "Source"
        << "\t"
        << "Fitness_source"
        << "\t"
        << "Target"
        << "\t"
        << "Fitness_target"
        << "\t"
        << "Phenos_found"
        << "\t"
        << "Fittest_found?"
        << "\t"
        << "u_size"
        << "\t"
        << "v_size" << std::endl;

  target_index = 0;
  for (target_index = 0; target_index < targets; target_index++) {
    std::cout << targets << std::endl;
    // Choose which fRNA phenotype is target and which is source.
    target = rand() % phenos_start.size();

    // Timing
    std::string geno;
    source_index = 0;
    for (source_index = 0; source_index < sources; source_index++) {
      std::cout << target_index << "." << source_index << std::endl;
      if (TIMING == 1) start = clock();

      do {
        source = rand() % phenos_start.size();
      } while (source == target);

      // If load_genos_, then get source from vector
      if (load_genos_) {
        geno = genos_start[source];
      }
      // Else find a genotype with inverse_fold
      else {
        geno = FindGenotype(phenos_start[source]);
      }

      if (TIMING == 1) middle = clock();

      SingleRunDFS(geno);

      if (TIMING == 1) end = clock();

      if (TIMING == 1)
        times << target_index << "." << source_index << "\t"
              << "Geno time: " << (double)(middle - start) / CLOCKS_PER_SEC
              << "\t"
              << "Search time: " << (double)(end - middle) / CLOCKS_PER_SEC
              << "\t"
              << "Source time: " << (double)(end - start) / CLOCKS_PER_SEC
              << std::endl;
    }
  }
  stats.close();
  if (TIMING == 1) times.close();
}

bool DoubleEqual(double lhs, double rhs) {
  return std::fabs(lhs - rhs) < std::numeric_limits<double>::epsilon();
}

double ProbFixationMutant(double f_p, double f_q, int N) {
  // Kimura fixation probability
  // Probability a given mutant will fix in population.
  if (DoubleEqual(f_p, f_q)) {
    return 1. / N;
  }
  // Selection coefficient calculated of p relative to q
  double s = f_p / f_q - 1;
  return (1 - exp(-2 * s)) / (1 - exp(-2 * N * s));
}

void FW_frna::SingleRunMono(std::string geno) {
  // Clear each of the sets used for set information
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  u.clear();
  v.clear();
  p.clear();
  genos.clear();
  fitness.clear();

  std::string source_shape = shape(phenos_start[source], shape_level);
  std::string target_shape = shape(phenos_start[target], shape_level);

  fitness[trivial_shape] = 0.0;
  fitness[target_shape] = 1.0;
  SetUpFitness(source_shape);

  if (debug) {
    for (auto &it : fitness) {
      std::cout << it.first << "\n";
      std::cout << it.second << "\n";
    }
    std::cout << std::endl;
  }

  double f_g = 0.0;
  double f_n = 0.0;

  std::string g;
  g.resize(L);
  std::string g_p;
  g_p.resize(L);
  std::string g_p_prev;
  g_p_prev.resize(L);
  std::string n;
  n.resize(L);
  std::string n_p;
  n_p.resize(L);

  u.insert(geno);

  std::strcpy(char_seq, geno.c_str());
  char_seq[L] = '\0';
  Fold(char_seq, phys_seq);
  n_p.assign(phys_seq);
  n_p = shape(n_p, shape_level);
  genos[geno] = std::string(phys_seq);

  neighbour_fitnesses.clear();
  neighbour_fitnesses.reserve((int)(K - 1) * L);

  // If u is empty, then run aground.
  // If u or v are too big, then
  if (source_index == -1 and target_index == -1)
    print_component = 1;
  else
    print_component = 0;

  if (print_component == 1) fStartComponent(target_index, source_index);

  int generation = 0;

  // Make initial genotype, phenotype and fitness
  g = geno;
  std::strcpy(char_seq, g.c_str());
  char_seq[L] = '\0';
  Fold(char_seq, phys_seq);
  g_p.assign(phys_seq);
  g_p = shape(g_p, shape_level);
  f_g = SetUpFitness(g_p);

  while (((generation < threshold) and (record_diffs_only == 0)) or
         ((generation < threshold) and (record_diffs_only == 1) and
          (g_p != target_shape))) {
    // Set g_p_prev at start of loop to g_p and update g_p inside loop
    g_p_prev = g_p;

    std::strcpy(char_seq, g.c_str());
    char_seq[L] = '\0';
    Fold(char_seq, phys_seq);
    g_p.assign(phys_seq);
    g_p = shape(g_p, shape_level);
    f_g = SetUpFitness(g_p);

    // Print to file
    if (((generation == 0 or g_p != g_p_prev) and (record_diffs_only == 1)) or
        (record_diffs_only == 0)) {
      stats << target_index << "\t" << source_index << "\t" << generation
            << "\t" << g << "\t" << g_p << "\t" << f_g << "\t" << target_shape
            << "\t" << fitness[target_shape] << "\t" << p.size() << "\t"
            << (p.find(target_shape) != p.end()) << "\t" << u.size() << "\t"
            << v.size() << std::endl;
    }

    if (g_p == target_shape) {
      break;
    }

    p.insert(g_p);

    // Set-up neighbour at g
    n = g;

    // Clear neighbourbood vectors
    neighbour_fitnesses.clear();
    if (debug) {
      std::cout << "Neigh comps: " << std::endl;
    }
    // Target every 1-mutant of genotype "trial"
    for (int j = 0; j < L; j++) {
      char save = n[j];
      for (int k = 0; k < K; k++) {
        if (save == alphabet[k]) continue;

        n[j] = alphabet[k];

        // Convert string into phenotype using GP program.
        std::strcpy(char_seq, n.c_str());
        char_seq[L] = '\0';
        Fold(char_seq, phys_seq);
        n_p.assign(phys_seq);
        n_p = shape(n_p, shape_level);
        genos[n] = n_p;

        // Get fitness of neighbour or assign it
        f_n = SetUpFitness(n_p);

        // Set the fitness to zero if neutral and no neutral mutations allowed
        if (neutral_mutations_ != 1 and n_p == g_p) {
          f_n = 0.;
        }

        // Add phenotype to set of phenotypes discovered
        p.insert(n_p);

        // Info for stdout if debugging
        if (debug) {
          std::cout << "l=" << j << ", "
                    << "k=" << alphabet[k] << "\t"
                    << "Phenos: " << g_p << "\t" << n_p << "\t"
                    << "Fitness: " << f_g << "\t" << f_n << "\t"
                    << "Prob: " << ProbFixationMutant(f_n, f_g, population_size)
                    << "\t"
                    << "Same: " << (g_p == n_p) << std::endl;
        }

        neighbour_fitnesses.push_back(f_n);
      }
      // if (print_component == 1) fAddEdge(g, n);

      // Make jth postion the original base again.
      n[j] = save;
    }

    // Loop to determine which genotype of the neighbours to
    // transition to next. Markov Chain.
    // Convert fitnss vector to fixation vector:
    // 1. Assume all mutants equally likely to be produced
    // 2. Use Eq. 2.14 SS: P(s) = (1-exp(-2*s))/(1-exp(-2NS))
    //    as the probability that a mutant goes to fixation
    //
    // The above formulation is described in:
    // https://doi.org/10.1086/677571
    //
    // We use the standard definition of selection coefficient for type
    // p compared to q by measuring the proportional fitness of p to q
    if (debug) {
      std::cout << "Kimura probs:" << std::endl;
    }
    for (unsigned int i = 0; i < neighbour_fitnesses.size(); i++) {
      if (neighbour_fitnesses[i] > std::numeric_limits<double>::min()) {
        neighbour_fitnesses[i] =
            ProbFixationMutant(neighbour_fitnesses[i], f_g, population_size);

        neighbour_fitnesses[i] = std::max(std::numeric_limits<double>::min(),
                                          neighbour_fitnesses[i]);
      }

      if (debug) {
        std::cout << neighbour_fitnesses[i] << std::endl;
      }
    }
    // Convert to cumulative vector
    double total = neighbour_fitnesses[0];
    if (debug) {
      std::cout << "Cumulative" << std::endl;
    }
    for (unsigned int i = 1; i < neighbour_fitnesses.size(); i++) {
      neighbour_fitnesses[i] += neighbour_fitnesses[i - 1];
      total = neighbour_fitnesses[i];
    }

    // Assert total greater than numeric limits
    assert(total > std::numeric_limits<double>::min());

    if (debug) {
      std::cout << "Normed cumulative probs:" << std::endl;
    }
    for (unsigned int i = 0; i < neighbour_fitnesses.size(); i++) {
      neighbour_fitnesses[i] = neighbour_fitnesses[i] / total;
      if (debug) {
        std::cout << neighbour_fitnesses[i] << std::endl;
      }
    }

    // Pick a neighbour at random according to convert fixation prob from
    // fitness vector
    double f_r = (double)rand() / RAND_MAX;
    double previous = 0;
    int counter = 0;
    n = g;
    for (int j = 0; j < L; j++) {
      char save = n[j];
      for (int k = 0; k < K; k++) {
        if (save == alphabet[k]) continue;

        n[j] = alphabet[k];

        if (f_r < neighbour_fitnesses[counter] and f_r >= previous) {
          if (debug) {
            std::cout << "Old geno: " << g << std::endl;
          }
          g = n;
          if (debug) {
            std::cout << "New geno: " << g << std::endl;
          }
        }
        previous = neighbour_fitnesses[counter];
        counter++;
      }
      // Make jth postion the original base again.
      n[j] = save;
    }

    // Assert that the next genotype g must be different to n, otherwise above
    // loop
    // failed
    assert(g != n);

    // Increment generation
    generation++;
    if (debug) {
      std::cin.get();
    }
  }

  free(char_seq);
  free(phys_seq);

  if (print_component == 1) fEndComponent();
}

void FW_frna::RunMono() {
  stats << "Target_index"
        << "\t"
        << "Source_index"
        << "\t"
        << "Generation"
        << "\t"
        << "Genotype"
        << "\t"
        << "Source"
        << "\t"
        << "Fitness_source"
        << "\t"
        << "Target"
        << "\t"
        << "Fitness_target"
        << "\t"
        << "Phenos_found"
        << "\t"
        << "Fittest_found?"
        << "\t"
        << "u_size"
        << "\t"
        << "v_size" << std::endl;

  target_index = 0;
  for (target_index = 0; target_index < targets; target_index++) {
    std::cout << targets << std::endl;
    // Choose which fRNA phenotype is target and which is source.
    target = rand() % phenos_start.size();

    // Timing
    std::string geno;
    source_index = 0;
    for (source_index = 0; source_index < sources; source_index++) {
      std::cout << target_index << "." << source_index << std::endl;
      if (TIMING == 1) start = clock();

      do {
        source = rand() % phenos_start.size();
      } while (source == target);

      // If load_genos_, then get source from vector
      if (load_genos_) {
        geno = genos_start[source];
      }
      // Else find a genotype with inverse_fold
      else {
        geno = FindGenotype(phenos_start[source]);
      }

      if (TIMING == 1) middle = clock();

      SingleRunMono(geno);

      if (TIMING == 1) end = clock();

      if (TIMING == 1)
        times << target_index << "." << source_index << "\t"
              << "Geno time: " << (double)(middle - start) / CLOCKS_PER_SEC
              << "\t"
              << "Search time: " << (double)(end - middle) / CLOCKS_PER_SEC
              << "\t"
              << "Source time: " << (double)(end - start) / CLOCKS_PER_SEC
              << std::endl;
    }
  }
  stats.close();
  if (TIMING == 1) times.close();
}
