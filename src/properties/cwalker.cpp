#include "cwalker.hpp"

nc::nc(int base, int genome_length, int ND_REF, int ND_INCLUDE,
       int ND_THRESHOLD, int SWAPS, int WRITE_COMPONENTS,
       int ESTIMATE_COMPONENT_STATS, std::string GENO_FNAME,
       std::string FILE_NAME, std::string OUTPATH, int VERBOSE, int DEBUG) {
  r_c = 0.0;
  e_wagner_c = 0.0;
  e_cowperthwaite_c = 0.0;
  f_p = 0;
  r_p = 0.0;
  e_wagner_p = 0.0;
  e_cowperthwaite_p = 0.0;

  write_components = WRITE_COMPONENTS;
  estimate_component_stats = ESTIMATE_COMPONENT_STATS;

  // ND integer in geno_list
  nd_ref = ND_REF;
  // Calculate stats for ND phenotype as well
  nd_include = ND_INCLUDE;
  // Maximum number of degenerate phenotypes for measurement
  threshold_nd = ND_THRESHOLD;

  sparse = 1;

  swaps = SWAPS;

  // File names
  geno_fname_ = GENO_FNAME;
  file_name_ = FILE_NAME;
  outpath_ = OUTPATH;

  // Verbosity for printing
  verbose_ = VERBOSE;

  // Debug mode
  debug_ = DEBUG;

  // Init call
  Init(base, genome_length);
}

nc::~nc() { Clear(); }

void nc::Init(int base, int genome_length) {
  L = genome_length;
  K = base;
  S = (int)pow(K, L);

  genotype_component_list.resize(S, -1);

  std::vector<int> temp;
  temp.reserve(5);

  std::ifstream in;
  in.open(geno_fname_);
  std::string line;
  int max = 0;
  int old_max = 0;
  int count = 0;
  int g_count = 0;
  genotype_list.resize(S * (threshold_nd + 1));
  while (std::getline(in, line)) {
    std::stringstream ls(line);
    std::string cell;
    count = 0;
    old_max = max;
    temp.clear();
    while (std::getline(ls, cell, ',')) {
      int phenotype = String2Number(cell);
      temp.push_back(phenotype);
      if (phenotype > max) max = phenotype;
    }

    if (temp.size() > threshold_nd) {
      genotype_list[o2(0, g_count, threshold_nd + 1)] = nd_ref;
      genotype_list[o2(threshold_nd, g_count, threshold_nd + 1)] = 1;
      max = old_max;
    } else {
      genotype_list[o2(threshold_nd, g_count, threshold_nd + 1)] = temp.size();
      for (unsigned int i = 0; i < temp.size(); i++)
        genotype_list[o2(i, g_count, threshold_nd + 1)] = temp[i];
    }
    g_count++;
  }
  in.close();

  N_p = max + 1;

  if (verbose_) {
    std::cout << N_p << std::endl;
  }

  // Only works for nd=1
  int scount = 0;
  int temp1 = -1;
  int temp_nd = -1;
  while (scount < swaps) {
    int a = rand() % S;
    int b = rand() % S;
    if (genotype_list[o2(0, a, threshold_nd + 1)] !=
        genotype_list[o2(0, b, threshold_nd + 1)]) {
      temp1 = genotype_list[o2(0, a, threshold_nd + 1)];
      temp_nd = genotype_list[o2(threshold_nd, a, threshold_nd + 1)];

      genotype_list[o2(0, a, threshold_nd + 1)] =
          genotype_list[o2(0, b, threshold_nd + 1)];
      genotype_list[o2(threshold_nd, a, threshold_nd + 1)] =
          genotype_list[o2(threshold_nd, b, threshold_nd + 1)];

      genotype_list[o2(0, b, threshold_nd + 1)] = temp1;
      genotype_list[o2(threshold_nd, b, threshold_nd + 1)] = temp_nd;

      scount++;
    }
  }

  component_transitions.resize(N_p, 0);
  phenotype_transitions.resize(N_p, 0);
}

void nc::Clear() {}

void nc::fStartComponent(int p, int n) {
  char *String;
  String = new char[1000];
  sprintf(String, "%scomponent_%s_p-%d_n-%d.txt", outpath_.c_str(),
          file_name_.c_str(), p, n);
  fp1.open(String);
  delete[] String;
}

void nc::fAddEdge(int x, int y) { fp1 << x << "\t" << y << "\n"; }

void nc::fEndComponent() { fp1.close(); }

void nc::Transfer(std::set<int> &from, std::set<int> &to, int x) {
  from.erase(x);
  to.insert(x);
}

void nc::CalculateComponentStatistics() {
  r_c = (double)component_transitions[phenotype] / (v.size() * (K - 1) * L);
  e_wagner_c = 0.0;
  e_wagner_c_nond = 0.0;

  int sum_nu = 0;
  int sum_nu_nond = 0;
  for (int i = 0; i < N_p; i++) {
    if (component_transitions[i] > 0 and phenotype != i) {
      e_wagner_c++;
      sum_nu += component_transitions[i];
      if (i != nd_ref) {
        e_wagner_c_nond++;
        sum_nu_nond += component_transitions[i];
      }
    }
  }

  double sum_f = 0;
  double sum_f_nond = 0;
  for (int i = 0; i < N_p; i++) {
    if (component_transitions[i] > 0 and phenotype != i) {
      sum_f += pow((double)component_transitions[i] / sum_nu, 2);
      if (i != nd_ref)
        sum_f_nond += pow((double)component_transitions[i] / sum_nu, 2);
    }
  }
  e_cowperthwaite_c = 1 - sum_f;
  e_cowperthwaite_c_nond = 1 - sum_f_nond;
}

void nc::CalculatePhenotypeStatistics() {
  f_p = VectorSum(component_volumes);
  r_p = (double)phenotype_transitions[phenotype] / (f_p * (K - 1) * L);
  e_wagner_p = 0.0;
  e_wagner_p_nond = 0.0;

  int sum_nu = 0;
  int sum_nu_nond = 0;
  for (int i = 0; i < N_p; i++) {
    if (phenotype_transitions[i] > 0 and phenotype != i) {
      e_wagner_p++;
      sum_nu += phenotype_transitions[i];
      if (i != nd_ref) {
        e_wagner_p_nond++;
        sum_nu_nond += phenotype_transitions[i];
      }
    }
  }

  double sum_f = 0.0;
  double sum_f_nond = 0.0;
  for (int i = 0; i < N_p; i++) {
    if (phenotype_transitions[i] > 0 and phenotype != i) {
      sum_f += pow((double)phenotype_transitions[i] / sum_nu, 2);
      if (i != nd_ref)
        sum_f_nond += pow((double)phenotype_transitions[i] / sum_nu_nond, 2);
    }
  }
  e_cowperthwaite_p = 1 - sum_f;
  e_cowperthwaite_p_nond = 1 - sum_f_nond;
}

void nc::WriteCompStatsHeader(std::ofstream &file) {
  /*
  Component stats header
  0: phenotype
  1: component,
  2: Frequency,
  3: surface_0,
  4: volume_1,
  5: surface_1,
  6: robustness,
  7: evolvability_wagner,
  8: evolvability_cowperthwaite,
  9: evolvability_wagner_nond,
  10:evolvability_cowperthwaite_nond,
  11:S_1, # S_1 = f * (1 - r_c)
  12:S_2  # S_1 = V_1 * 2 * (1 / r_c - 1)
  */
  file << "phenotype"
       << "\t"
       << "component"
       << "\t"
       << "Frequency"
       << "\t"
       << "surface_0"
       << "\t"
       << "volume_1"
       << "\t"
       << "surface_1"
       << "\t"
       << "robustness"
       << "\t"
       << "evolvability_wagner"
       << "\t"
       << "evolvability_cowperthwaite"
       << "\t"
       << "evolvability_wagner_nond"
       << "\t"
       << "evolvability_cowperthwaite_nond"
       << "\t"
       << "S_1"
       << "\t";                // S_1 = f * (1 - r_c)
  file << "S_2" << std::endl;  // S_1 = V_1 * 2 * (1 / r_c - 1)
}

void nc::WritePhenoStatsHeader(std::ofstream &file) {
  /*
  Phenotype stats header
  0: phenotype
  1:n_components,
  2:Frequency, # volume_0
  3:surface_0, # Sum of component surfaces
  4:volume_1, # ( f_p * (K-1) * L * r_p )/(double)2
  5:surface_1,
  6:robustness,
  7:evolvability_wagner,
  8:evolvability_cowperthwaite,
  9:evolvability_wagner_nond,
  10:evolvability_cowperthwaite_nond,
  11:S_1, # f_p * (1 - r_p)
  12:S_2, # S_1 = V_1 * 2(1/r_p - 1)
  13:surface_unique
  */
  file << "phenotype"
       << "\t"
       << "n_components"
       << "\t"
       << "Frequency"
       << "\t";  // volume_0
  file << "surface_0"
       << "\t";  // Sum of component surfaces
  file << "volume_1"
       << "\t";  // (f_p * (K-1) * L * r_p )/(double)2
  file << "surface_1"
       << "\t"
       << "robustness"
       << "\t"
       << "evolvability_wagner"
       << "\t"
       << "evolvability_cowperthwaite"
       << "\t"
       << "evolvability_wagner_nond"
       << "\t"
       << "evolvability_cowperthwaite_nond"
       << "\t"
       << "S_1"
       << "\t";  // S_1 = f * (1-r_c)
  file << "S_2"
       << "\t";  // S_1 = V_1 * 2(1/r_c - 1)
  file << "surface_unique" << std::endl;
}

void nc::fPrintComponentStatistics(std::ofstream &file) {
  CalculateComponentStatistics();
  // phenotype - component - volume_0 - surface_0 - volume_1 - surface_1 -
  // robustness_c - E_w_c - E_cowperthwaite_c - E_w_c_nond -
  // E_cowperthwaite_c_nond - S_1 = f * (1-r_c) -  S_1 = V_1 * 2(1/r_c - 1)
  double v_1 = (v.size() * (K - 1) * L * r_c) / (double)2;
  file << phenotype << "\t" << component_volumes.size() - 1 << "\t" << v.size()
       << "\t" << s.size() << "\t" << v_1 << "\t"
       << (K - 1) * L * (1 - r_c) * v.size() << "\t" << r_c << "\t"
       << e_wagner_c << "\t" << e_cowperthwaite_c << "\t" << e_wagner_c_nond
       << "\t" << e_cowperthwaite_c_nond << "\t" << v.size() * (1 - r_c) << "\t"
       << v_1 * 2 * ((double)1 / r_c - 1) << std::endl;
}

void nc::fPrintPhenotypeStatistics(std::ofstream &file) {
  CalculatePhenotypeStatistics();
  // phenotype - components - volume_0 - surface_0 - volume_1 - surface_1 -
  // robustness - E_w_p - E_cowperthwaite_p - E_w_p_nond -
  // E_cowperthwaite_p_nond - S_1 = ((a-1)*L) * f * (1-r_p) - S_1 = V_1 *
  // 2(1/r_p - 1) - unique phenotype surface (sp.size)

  double v_1 = (f_p * (K - 1) * L * r_p) / (double)2;
  file << phenotype << "\t" << component_volumes.size() << "\t" << f_p << "\t"
       << VectorSum(component_surfaces) << "\t" << v_1 << "\t"
       << (K - 1) * L * (1 - r_p) * f_p << "\t" << r_p << "\t" << e_wagner_p
       << "\t" << e_cowperthwaite_p << "\t" << e_wagner_p_nond << "\t"
       << e_cowperthwaite_p_nond << "\t" << f_p * (1 - r_p) << "\t"
       << v_1 * 2 * ((double)1 / r_p - 1) << "\t" << sp.size() << std::endl;
}

void nc::Run() {
  int *sequence;
  sequence = new int[L];

  std::stringstream ss;
  std::ofstream f_component_genos;
  std::ofstream f_comp;
  if (estimate_component_stats) {
    // GP map statistics for each component (each row is a component)
    ss.str("");
    ss << outpath_ << "component_stats_" << file_name_ << ".txt";
    f_comp.open(ss.str().c_str());

    if (verbose_) {
      std::cout << "Component stats path: " << ss.str() << std::endl;
    }

    // The component index for each genotype
    ss.str("");
    ss << outpath_ << "component_genotypes_" << file_name_ << ".txt";
    f_component_genos.open(ss.str().c_str());
    if (verbose_) {
      std::cout << "Component genotypes path: " << ss.str() << std::endl;
    }
  }

  // GP map statistics for each phenotype (each row is a phenotype)
  ss.str("");
  ss << outpath_ << "phenotype_stats_" << file_name_ << ".txt";
  std::ofstream f_pheno(ss.str().c_str());

  if (verbose_) {
    std::cout << "Phenotype stats path: " << ss.str() << std::endl;
  }

  // Edge list of phenotype connectivity
  ss.str("");
  ss << outpath_ << "phenotype_transition_stats_" << file_name_ << ".txt";
  std::ofstream f_pheno_trans(ss.str().c_str());
  if (sparse) {
    f_pheno_trans << "from"
                  << "\t"
                  << "to"
                  << "\t"
                  << "count" << std::endl;
  }
  if (verbose_) {
    std::cout << "Phenotype transitions path: " << ss.str() << std::endl;
  }

  // Write headers for files
  WriteCompStatsHeader(f_comp);
  WritePhenoStatsHeader(f_pheno);

  // Run over all the phenotypes
  phenotype = 0;
  for (phenotype = 0; phenotype < N_p; phenotype++) {
    if (phenotype == nd_ref and nd_include == 0) continue;

    if (verbose_) {
      std::cout << "phenotype: " << phenotype << std::endl;
    }

    // Clear each of the sets used for set information
    w.clear();
    u.clear();
    v.clear();
    s.clear();

    for (unsigned int j = 0; j < S; j++) {
      // Relaxed
      if (genotype_list[o2(threshold_nd, j, threshold_nd + 1)] > threshold_nd)
        continue;
      for (unsigned int i = 0;
           i < genotype_list[o2(threshold_nd, j, threshold_nd + 1)]; i++)
        if (genotype_list[o2(i, j, threshold_nd + 1)] == phenotype) w.insert(j);
    }

    int total_p = w.size();
    if (verbose_) {
      std::cout << "w: " << w.size() << std::endl;
    }

    clock_t runner = clock();
    sp.clear();

    while (w.size() > 0) {
      Transfer(w, u, *w.begin());

      if (write_components)
        fStartComponent(phenotype, component_volumes.size());

      while (u.size() > 0) {
        // Mark down new component phenotype
        if (u.size() == 1 and v.size() == 0) {
          phenotype_component_list.push_back(phenotype);
        }

        // Readout at each u genotype visited
        if (debug_) {
          clock_t now = clock();
          std::cout << "Time from last = "
                    << (double)(now - runner) / CLOCKS_PER_SEC << std::endl;
          runner = now;
          std::cout << "w: " << w.size() << "\tu: " << u.size()
                    << "\tv: " << v.size() << "\ts: " << s.size() << std::endl;
        }

        // Take first element of u and attempt all 1-mutants.
        int trial = *u.begin();
        int_to_basek(trial, sequence, L, K);

        // Test every 1-mutant of genotype "trial"
        for (int j = 0; j < L; j++) {
          int save = sequence[j];
          for (int k = 0; k < K; k++) {
            if (save == k) continue;
            sequence[j] = k;

            int robust = 0;
            int neighbour = basek_to_int(sequence, L, K);

            // Relaxed
            int neighbour_phenotype = -1;
            for (unsigned int i = 0;
                 i <
                 genotype_list[o2(threshold_nd, neighbour, threshold_nd + 1)];
                 i++) {
              // Get neighbour phenotype
              neighbour_phenotype =
                  genotype_list[o2(i, neighbour, threshold_nd + 1)];

              // Is it a neutral mutations
              if (neighbour_phenotype == phenotype) robust = 1;

              // Count transitions
              component_transitions[neighbour_phenotype]++;
              phenotype_transitions[neighbour_phenotype]++;
            }
            // Assert that phenotype is assigned
            assert(neighbour_phenotype != -1);
            // if (phenotype == -1) {
            //   std::cout << "u: " << u.size() << std::endl;
            //   std::cout << "v: " << v.size() << std::endl;
            //   std::cout << "w: " << w.size() << std::endl;
            // }

            // Updated surface sets
            if (robust == 0) {
              s.insert(neighbour);
              // TODO: fix to neighbour_phenotype?
              sp.insert(neighbour);
              continue;
            }

            // Add edge to the component graph
            if (write_components) {
              fAddEdge(trial, neighbour);
            }

            if (u.find(neighbour) != u.end()) continue;
            if (v.find(neighbour) != v.end()) continue;
            if (w.find(neighbour) == w.end()) continue;

            // Move neighbour from w to u for inclusion in component stats
            Transfer(w, u, neighbour);
          }
          // if (phenotype == -1) std::cout << "j: " << j << "\tL: " << L <<
          // std::endl;
          sequence[j] = save;
        }
        // Move trial from u to v
        Transfer(u, v, trial);
      }

      component_volumes.push_back(v.size());
      component_surfaces.push_back(s.size());
      if (write_components) fEndComponent();

      // Print stats to file...
      if (estimate_component_stats == 1) fPrintComponentStatistics(f_comp);

      // Update genotype_component_list
      for (std::set<int>::iterator it = v.begin(); it != v.end(); it++) {
        genotype_component_list[*it] = phenotype_component_list.size() - 1;
      }

      // Clear s,v
      v.clear();
      s.clear();
      for (unsigned int el = 0; el < component_transitions.size(); el++)
        component_transitions[el] = 0;
    }
    // Print stats about phenotype
    fPrintPhenotypeStatistics(f_pheno);

    // component_volumes, surfaces, transitions clear.
    component_volumes.clear();
    component_surfaces.clear();
    for (unsigned int el = 0; el < phenotype_transitions.size(); el++) {
      if (sparse == 0)
        f_pheno_trans << phenotype_transitions[el] << "\t";
      else if (phenotype_transitions[el] != 0)
        f_pheno_trans << phenotype << "\t" << el << "\t"
                      << phenotype_transitions[el] << std::endl;

      phenotype_transitions[el] = 0;
    }
    if (sparse == 0) f_pheno_trans << std::endl;
  }

  if (estimate_component_stats == 1 and nd_include) {
    for (unsigned int i = 0; i < genotype_component_list.size(); i++) {
      f_component_genos << genotype_component_list[i] << "\n";
    }
    CalculateComponentTransitions();
  }

  delete[] sequence;

  f_component_genos.close();
  f_pheno.close();
  f_comp.close();
  f_pheno_trans.close();
}

void nc::CalculateComponentTransitions() {
  int components = phenotype_component_list.size();

  int *sequence;
  sequence = new int[L];

  std::stringstream ss;
  // Component transition stats
  // Format: edge list of component transition counts between components
  //   from_component_idx \t to_coponent_idx \t count
  ss.str("");
  ss << outpath_ << "component_transition_stats_" << file_name_ << ".txt";
  std::ofstream f_comp_trans(ss.str().c_str());
  if (sparse) {
    f_comp_trans << "from"
                 << "\t"
                 << "to"
                 << "\t"
                 << "count" << std::endl;
  }
  if (verbose_) {
    std::cout << "Component transition stats path: " << ss.str() << std::endl;
  }
  // Component phenotypes
  // Format: list of length number of components with the phenotype index of
  //         each component
  ss.str("");
  ss << outpath_ << "component_phenotypes_" << file_name_ << ".txt";
  std::ofstream f_comp_pheno(ss.str().c_str());
  if (verbose_) {
    std::cout << "Component phenotypes path: " << ss.str() << std::endl;
  }

  // Resize based on number of components
  component_ind_transitions.resize(components);
  if (verbose_) {
    std::cout << "Components: " << components << std::endl;
  }

  // Run over all components
  for (int component = 0; component < components; component++) {
    if (verbose_) {
      std::cout << component << std::endl;
    }

    for (unsigned int el = 0; el < component_ind_transitions.size(); el++)
      component_ind_transitions[el] = 0;

    // Run over all genotypes
    for (unsigned int genotype = 0; genotype < genotype_component_list.size();
         genotype++)
      if (genotype_component_list[genotype] == component) {
        int_to_basek(genotype, sequence, L, K);

        // Test every 1-mutant of genotype "trial"
        for (int j = 0; j < L; j++) {
          int save = sequence[j];
          for (int k = 0; k < K; k++) {
            if (save == k) continue;
            sequence[j] = k;

            int neighbour = basek_to_int(sequence, L, K);
            int neighbour_component = genotype_component_list[neighbour];

            component_ind_transitions[neighbour_component]++;
          }
          sequence[j] = save;
        }
      }

    for (unsigned int i = 0; i < component_ind_transitions.size(); i++) {
      if (sparse == 0)
        f_comp_trans << component_ind_transitions[i] << "\t";
      else if (component_ind_transitions[i] > 0)
        f_comp_trans << component << "\t" << i << "\t"
                     << component_ind_transitions[i] << std::endl;
    }
    if (sparse == 0) f_comp_trans << std::endl;
  }

  for (unsigned int i = 0; i < phenotype_component_list.size(); i++) {
    f_comp_pheno << phenotype_component_list[i] << std::endl;
  }

  delete[] sequence;

  f_comp_trans.close();
  f_comp_pheno.close();
}
