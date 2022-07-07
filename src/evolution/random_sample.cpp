#include "random_sample.hpp"

using namespace std;

const int SEED = 1;
const int ATTEMPTS = 10;
const int SAMPLES = 100000;
const int L = 50;
const char *OUTPUT0 = "out.txt";

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

int main() {
  std::ofstream out0(OUTPUT0);

  // No need Vienna2
  initialize_fold(L);
  //
  // give_up = 1;

  clock_t start, end;
  start = clock();
  srand(SEED);
  int randomise_rand = rand();

  std::map<std::string, int> phenos;

  char *char_seq;
  char_seq = (char *)malloc(sizeof(char) * (L + 1));
  char_seq[L] = '\0';
  char *phys_seq;
  phys_seq = (char *)malloc(sizeof(char) * (L + 1));
  phys_seq[L] = '\0';

  std::string temp;
  for (int j = 0; j < SAMPLES; j++) {
    for (int i = 0; i < L; i++) {
      char_seq[i] = RandomBase();
    }

    Fold(char_seq, phys_seq);
    temp = std::string(phys_seq);

    if (phenos.find(temp) != phenos.end()) {
      phenos[temp] += 1;
    } else {
      phenos[temp] = 1;
    }
  }
  free(char_seq);
  free(phys_seq);

  std::ofstream out;
  char *buffer0;
  buffer0 = new char[50];
  sprintf(buffer0, "phenotype_sample_L%d.txt", L);
  out.open(buffer0);
  delete[] buffer0;

  for (std::map<std::string, int>::iterator it = phenos.begin();
       it != phenos.end(); it++) {
    out << it->first << "\t" << it->second << std::endl;
  }
  out.close();

  // Finish program timer
  end = clock();
  out0 << "Program time = " << (double)(end - start) / CLOCKS_PER_SEC
       << " seconds" << endl;
  out0.close();

  return 0;
}
