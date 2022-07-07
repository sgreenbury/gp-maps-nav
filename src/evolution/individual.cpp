#include "individual.hpp"

Individual::Individual() {
  K = 0;
  L = 0;
  genotype = "";
  phenotype = "";
  genotype_prev = "";
  phenotype_prev = "";
  mem_assigned = 0;
  fitness = 0.0;
}

Individual::Individual(const Individual& ind) { Init(ind); }

Individual::Individual(int alphabet_size, int sequence_length) {
  Init(alphabet_size, sequence_length);
}

void Individual::MakeAlphabet(int alphabet_size) {
  alphabet.resize(alphabet_size);
  for (int i = 0; i < alphabet_size; i++) {
    switch (i) {
      case 0:
        alphabet[i] = 'G';
        break;
      case 1:
        alphabet[i] = 'C';
        break;
      case 2:
        alphabet[i] = 'A';
        break;
      case 3:
        alphabet[i] = 'U';
        break;
    }
  }
}

void Individual::Init(int alphabet_size, int sequence_length) {
  K = alphabet_size;
  L = sequence_length;
  genotype.resize(L, '-');
  phenotype.resize(L, '-');
  genotype_prev.resize(L, '-');
  phenotype_prev.resize(L, '-');
  g_seq = (char*)malloc(sizeof(char) * (L + 1));
  p_seq = (char*)malloc(sizeof(char) * (L + 1));
  g_seq[L] = '\0';
  p_seq[L] = '\0';
  mem_assigned = 1;
  MakeAlphabet(K);
}

void Individual::AssignGenotype(const char* genotype_sequence) {
  for (int i = 0; i < L; i++) genotype[i] = genotype_sequence[i];
  G2P();
}

void Individual::AssignGenotype(const std::string genotype_sequence) {
  for (int i = 0; i < L; i++) genotype[i] = genotype_sequence[i];
  G2P();
}

void Individual::Init(const Individual& ind) {
  K = ind.K;
  L = ind.L;
  genotype = ind.genotype;
  phenotype = ind.phenotype;
  genotype_prev = ind.genotype_prev;
  phenotype_prev = ind.phenotype_prev;
  g_seq = (char*)malloc(sizeof(char) * (L + 1));
  p_seq = (char*)malloc(sizeof(char) * (L + 1));
  std::strcpy(g_seq, genotype.c_str());
  std::strcpy(p_seq, phenotype.c_str());
  fitness = ind.fitness;
  mem_assigned = 1;
  MakeAlphabet(K);
}

Individual::~Individual() {
  if (mem_assigned) {
    free(g_seq);
    free(p_seq);
  }
  // std::cout << "Destroyed..." << std::endl;
}

Individual& Individual::operator=(const Individual& other) {
  if (this != &other) {
    Clear();
    Init(other);
  }
  return *this;
}

void Individual::Mutate() {
  //  std::cout << genotype << std::endl;

  int pos = rand() % L;
  int base;
  do {
    base = rand() % K;
  } while (genotype[pos] == alphabet[base]);

  genotype[pos] = alphabet[base];
  G2P();
  // std::cout << genotype << std::endl;
}

void Individual::MutateBase(int pos) {
  int base;
  do {
    base = rand() % K;
  } while (genotype[pos] == alphabet[base]);

  genotype[pos] = alphabet[base];
}

void Individual::G2P() {
  std::strcpy(g_seq, genotype.c_str());
  fold(g_seq, p_seq);
  for (int i = 0; i < genotype.size(); i++) phenotype[i] = p_seq[i];
}

void Individual::Clear() {
  if (mem_assigned) {
    free(g_seq);
    free(p_seq);
  }
  //  std::cout << "Cleared..." << std::endl;
}

void Individual::Print() {
  std::cout << "Genotype: " << genotype << std::endl;
  std::cout << "Phenotype: " << phenotype << std::endl;
}

void Individual::RandomGenotype() {
  char ch = 'b';
  for (int i = 0; i < L; i++) {
    int el = rand() % 4;
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
    g_seq[i] = ch;
    genotype[i] = ch;
  }
}
