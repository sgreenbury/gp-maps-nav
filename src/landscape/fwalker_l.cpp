#include "fwalker_l.hpp"

FW_l::FW_l(int BASE, int GENOME_LENGTH, int ND, int TARGETS, int SOURCES,
           int ALL, int NEUTRAL_MUTATIONS, int THRESHOLD, int ASSEMBLY_TESTS,
           std::string FITNESS_ASSIGNMENT, std::string OUTPATH,
           std::string FILE_NAME, std::string PHENO_FNAME, std::string GP_MAP,
           int VERBOSE) {
  print_component = 0;
  sparse = 1;

  timing = 1;
  all = ALL;
  L = GENOME_LENGTH;
  K = BASE;
  nd_ref = ND;
  sources = SOURCES;
  targets = TARGETS;

  threshold = THRESHOLD;
  neutral_mutations_ = NEUTRAL_MUTATIONS;

  fitness_assignment = FITNESS_ASSIGNMENT;

  assembly_tests = ASSEMBLY_TESTS;

  outpath_ = OUTPATH;
  file_name_ = FILE_NAME;

  pheno_fname = PHENO_FNAME;
  gp_map = GP_MAP;

  verbose = VERBOSE;

  Init();
}

FW_l::~FW_l() { Clear(); }

void FW_l::SetAlphabet() {
  alphabet.resize(K);
  if (gp_map == "rna") {
    if (K == 4) {
      alphabet[0] = 'C';
      alphabet[1] = 'G';
      alphabet[2] = 'A';
      alphabet[3] = 'U';
    }
  } else if (gp_map == "polyomino") {
    for (int i = 0; i < K; i++) {
      alphabet[i] = '0' + (char)i;
    }
  }
}

void FW_l::Init() {
  SetAlphabet();

  if (gp_map == "rna") {
    std::ifstream in;
    in.open(pheno_fname);
    std::string line;
    int p_count = 0;
    while (std::getline(in, line)) {
      phenos[line] = p_count;
      p_count++;
    }
    in.close();

    N_p = p_count;

  } else if (gp_map == "polyomino") {
    // Set the number of tiles
    Nt = (int)(L / 4);

    // Load the phenotypes
    N_p = 0;
    while (1) {
      Phenotype P;
      std::string buffer = absl::StrFormat("%s/pheno%d.txt", pheno_fname, N_p);
      std::ifstream in(buffer);
      // If cannot open the file, then finished loading phenotypes
      if (!in.good()) {
        break;
      } else {
        in.close();
      }
      P.fLoad(buffer.c_str());
      phenos_vector.push_back(P);
      std::cout << N_p << std::endl;
      phenos_vector[N_p].PrintPic();
      N_p++;
    }
  }

  fitness.resize(N_p, -1);

  // File name format:
  // <outpath>/config_<file_name>.txt : this is done in nav_disk.cpp
  // <outpath>/out_<file_name>.txt
  // <outpath>/stats_<file_name>.txt
  // <outpath>/times_<file_name>.txt
  // <outpath>/components/component_<file_name>_t-<target>_s-<source>.txt
  std::string statsfile =
      absl::StrFormat("%sstats_%s.txt", outpath_.c_str(), file_name_.c_str());
  stats.open(statsfile.c_str());
  if (timing == 1) {
    std::string timesfile =
        absl::StrFormat("%stimes_%s.txt", outpath_.c_str(), file_name_.c_str());
    times.open(timesfile.c_str());
  }
}

int FW_l::Lookup(Phenotype P) {
  for (unsigned int i = 0; i < phenos_vector.size(); i++) {
    if (phenos_vector[i].MatchRots(P) == 1) return i;
  }
  return nd_ref;
}

void FW_l::Clear() {}

void FW_l::fStartComponent(int t, int s) {
  std::string componentfile =
      absl::StrFormat("%scomponents/component_%s_t-%d_s-%d.txt",
                      outpath_.c_str(), file_name_.c_str(), t, s);
  fp1.open(componentfile);
}

void FW_l::fAddEdge(std::string x, std::string y) {
  fp1 << x << "\t" << y << "\n";
}

void FW_l::fEndComponent() { fp1.close(); }

void FW_l::Transfer(std::set<int> &from, std::set<int> &to, int x) {
  from.erase(x);
  to.insert(x);
}

void FW_l::Transfer(std::unordered_set<int> &from, std::unordered_set<int> &to,
                    int x) {
  from.erase(x);
  to.insert(x);
}

void FW_l::Transfer(std::set<std::string> &from, std::set<std::string> &to,
                    std::string x) {
  from.erase(x);
  to.insert(x);
}

void FW_l::Transfer(std::unordered_set<std::string> &from,
                    std::unordered_set<std::string> &to, std::string x) {
  from.erase(x);
  to.insert(x);
}

void FW_l::AssignFitness() {
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

void FW_l::AssignFitness(int max) {
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

void FW_l::SetFitness() {
  if (all == 1)
    if (target != nd_ref)
      AssignFitness(target);
    else
      return;
  else
    AssignFitness();
}

void FW_l::IntSeq2String(int *&int_seq, std::string &s) {
  for (unsigned int i = 0; i < s.size(); i++) s[i] = '0' + (char)int_seq[i];
}

void FW_l::String2IntSeq(std::string &s, int *&int_seq) {
  for (unsigned int i = 0; i < s.size(); i++) int_seq[i] = s[i] - '0';
}

void FW_l::SingleRun(std::string geno) {
  int *int_seq;
  char *phys_seq;
  char *char_seq;

  // Assign memory for genotype and phenotypes for RNA or polyomino
  if (gp_map == "rna") {
    char_seq = (char *)malloc(sizeof(char) * (L + 1));
    char_seq[L] = '\0';
    phys_seq = (char *)malloc(sizeof(char) * (L + 1));
    phys_seq[L] = '\0';
  } else if (gp_map == "polyomino") {
    int_seq = new int[L];
  }

  // Clear and initialize the containers for performing search
  u.clear();
  v.clear();
  p.clear();
  genos.clear();

  std::string g;
  g.resize(L);

  std::string n;
  n.resize(L);

  std::string n_p;
  n_p.resize(L);

  int g_p_i;
  int n_p_i;

  u.insert(geno);

  if (gp_map == "rna") {
    std::strcpy(char_seq, geno.c_str());
    char_seq[L] = '\0';
    Fold(char_seq, phys_seq);
    n_p.assign(phys_seq);
    genos[geno] = phenos[n_p];
  } else if (gp_map == "polyomino") {
    String2IntSeq(geno, int_seq);
    Phenotype A = Assemble(int_seq, Nt, assembly_tests);
    g_p_i = Lookup(A);
    genos[geno] = g_p_i;
  }

  // Conditionally write explored component to file
  if (source == -1 and target == -1)
    print_component = 1;
  else
    print_component = 0;

  if (print_component == 1) fStartComponent(target, source);

  // If u is empty, then no more to search OR if u or v are too big OR fittest
  // found, exit loop
  while (u.size() > 0) {
    if ((u.size() + v.size()) > threshold or p.find(max_fit_p) != p.end())
      break;

    // Take first element of u and attempt all 1-mutants.
    // Copy beginning string of u to g.
    g = *u.begin();
    g_p_i = genos[g];

    n = g;

    // Test every 1-mutant of genotype "trial"
    for (int j = 0; j < L; j++) {
      char save = n[j];
      for (int k = 0; k < K; k++) {
        if (save == alphabet[k]) continue;
        n[j] = alphabet[k];

        if (genos.find(n) != genos.end()) {
          n_p_i = genos[n];
        } else {
          if (gp_map == "rna") {
            // Convert string into phenotype using GP program.
            std::strcpy(char_seq, n.c_str());
            char_seq[L] = '\0';

            Fold(char_seq, phys_seq);
            n_p.assign(phys_seq);
            n_p_i = phenos[n_p];

          } else if (gp_map == "polyomino") {
            String2IntSeq(n, int_seq);
            Phenotype X = Assemble(int_seq, Nt, assembly_tests);
            n_p_i = Lookup(X);
          }

          // Assign n_p_i to genos map
          genos[n] = n_p_i;
        }

        // If fitness n_p_i is < fitness of g_p_i or if neutral mutations not
        // allowed and phenotype is the same, continue
        if ((fitness[n_p_i] < fitness[g_p_i]) or
            ((neutral_mutations_ != 1) and n_p_i == g_p_i)) {
          continue;
        }

        p.insert(n_p_i);

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

  if (gp_map == "rna") {
    free(char_seq);
    free(phys_seq);
  } else if (gp_map == "polyomino") {
    delete[] int_seq;
  }
  if (print_component == 1) fEndComponent();

  // Print file
  stats << target << "\t" << source << "\t" << max_fit_p << "\t" << max_fit
        << "\t" << geno << "\t" << genos[geno] << "\t" << p.size() << "\t"
        << (p.find(max_fit_p) != p.end()) << "\t" << u.size() << "\t"
        << v.size() << std::endl;
}

void FW_l::Run() {
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
        << "v_size" << std::endl;

  int *int_seq;
  int_seq = new int[L];
  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  // If fitness landscape "fixed" for all source-target pairs, set fitness
  if (fitness_assignment == "fixed") {
    SetFitness();
  }

  // If all, test all phenotypes as the target, otherwise random sampling
  target = 0;
  if (all == 1) targets = N_p;
  for (target = 0; target < targets; target++) {
    if (verbose) std::cout << targets << endl;

    // If target is nd_ref and testing all phenotypes
    if (target == nd_ref and all == 1) continue;

    // If assignment mechanism for every target, set fitness
    if (fitness_assignment == "target") {
      SetFitness();
    }

    // Timing
    std::string geno(L, 'A');

    int geno_p;
    source = 0;
    for (source = 0; source < sources; source++) {
      if (verbose) std::cout << target << "." << source << endl;

      if (timing == 1) start = clock();

      // If assignment mechanism for every source, set fitness
      if (fitness_assignment == "source") {
        SetFitness();
      }

      do
        geno_p = rand() % N_p;
      while (geno_p == max_fit_p);

      if (gp_map == "rna") {
        for (std::map<std::string, int>::iterator it = phenos.begin();
             it != phenos.end(); ++it)
          if (it->second == geno_p) {
            std::strcpy(phys_seq, (it->first).c_str());
            break;
          }

        double found;
        do {
          for (int i = 0; i < L; i++) {
            int ri = rand() % K;
            char_seq[i] = alphabet[ri];
          }

          found = InverseFold(char_seq, phys_seq);
        } while (found != 0);
        geno.assign(char_seq);
      } else if (gp_map == "polyomino") {
        int found = -1;
        do {
          for (int i = 0; i < L; i++) {
            int ri = rand() % K;
            int_seq[i] = ri;
          }

          Phenotype A = Assemble(int_seq, Nt, assembly_tests);
          if (Lookup(A) == geno_p) found = 1;
        } while (found == -1);

        // Convert int_seq to string
        IntSeq2String(int_seq, geno);
      }

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

  free(char_seq);
  free(phys_seq);
  delete[] int_seq;
}
