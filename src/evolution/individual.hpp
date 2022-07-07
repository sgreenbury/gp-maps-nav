#ifndef _individual_h_included_
#define _individual_h_included_

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "base/utilities.hpp"
#include "landscape/rna_fold.hpp"

using namespace std;

class Individual {
 public:
  int L;
  int K;
  int mem_assigned;
  double fitness;

  std::vector<char> alphabet;
  std::string alph;
  std::string genotype;
  std::string genotype_prev;
  std::string phenotype;
  std::string phenotype_prev;

  char* g_seq;
  char* p_seq;

  void Init(int alphabet_size, int sequence_length);
  void Init(const Individual& ind);

  Individual();
  Individual(int alphabet_size, int sequence_length);
  Individual(const Individual& ind);
  ~Individual();

  void MakeAlphabet(int alphabet_size);

  void AssignGenotype(const char* genotype_sequence);
  void AssignGenotype(const std::string genotype);
  Individual& operator=(const Individual& other);

  void G2P();
  void Mutate();
  void MutateBase(int pos);
  void Clear();
  void Print();
  void RandomGenotype();
};

#endif
