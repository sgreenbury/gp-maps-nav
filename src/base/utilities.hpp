#ifndef _utilities_h_included_
#define _utilities_h_included_
#include <math.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int RelativeCoordinatePosition(int x1, int y1, int x2, int y2);

int Round(double x);

std::string MakeString(int *s, int n);

// Array Utilities

inline int Step(int value, int step_point) {
  if (value >= step_point) return 1;
  return 0;
}

inline int o2(int x, int y, int k) { return x + y * k; }

inline int o3(int x, int y, int z, int k, int m) {
  return x + y * k + z * k * m;
}

template <typename T>
void ClearArray(T *&array) {
  if (array != 0) {
    delete[] array;
    array = 0;
  }
}

template <typename T>
void resize_array(T *&array, size_t old_size, size_t new_size) {
  if (new_size <= old_size) return;
  T *temp = new (std::nothrow) T[new_size];
  if (temp == NULL) throw((const char *)"Memory allocation failed!\n");
  for (size_t i = 0; i < old_size; i++) temp[i] = array[i];
  if (array != NULL) delete[] array;
  array = temp;
  return;
}

template <typename T>
void Clean(T *&array, size_t n) {
  for (size_t i = 0; i < n; i++) array[i] = 0;
  return;
}

int String2Number(const std::string &text);

std::string Number2String(int number);

int Search(int *array, int v, int l, int r);

// Maths utilities
inline int Mod(int x, int y) {
  int result = x % y;
  if (result < 0) return result + y;
  return result;
}

double my_log(double base, double y);

int Delta(int i, int j);

template <typename T>
double Mean(T String, int N) {
  double count = 0;
  for (int i = 0; i < N; i++) count += String[i];
  count = count / N;
  return count;
}

template <typename T>
double SE(T String, int N) {
  double mean = Mean(String, N);
  double count = 0;
  for (int i = 0; i < N; i++)
    count += (double)((String[i] - mean) * (String[i] - mean));
  count = sqrt(count) / N;
  return count;
}

template <typename T>
T VectorSum(std::vector<T> data) {
  T count = 0;
  for (unsigned int i = 0; i < data.size(); i++) count += data[i];
  return count;
}

template <typename T>
double VectorMean(std::vector<T> data) {
  double count = 0;
  for (unsigned int i = 0; i < data.size(); i++) count += data[i];
  count = count / data.size();
  return count;
}

template <typename T>
double VectorSE(std::vector<T> data) {
  double mean = VectorMean(data);
  double count = 0;
  for (unsigned int i = 0; i < data.size(); i++)
    count += (double)((data[i] - mean) * (data[i] - mean));
  count = sqrt(count) / sqrt(data.size() * (data.size() - 1));
  return count;
}

template <typename T>
double VectorSD(std::vector<T> data) {
  double mean = VectorMean(data);
  double count = 0;
  for (unsigned int i = 0; i < data.size(); i++)
    count += (double)((data[i] - mean) * (data[i] - mean));
  count = sqrt(count / data.size());
  return count;
}

template <typename T>
double VectorSampleSD(std::vector<T> data) {
  if (data.size() == 1) return -1;

  double mean = VectorMean(data);
  double count = 0;
  for (unsigned int i = 0; i < data.size(); i++)
    count += (double)((data[i] - mean) * (data[i] - mean));
  count = sqrt(count / (data.size() - 1));
  return count;
}

template <typename T>
T VectorMax(std::vector<T> data) {
  T max;
  for (unsigned int i = 0; i < data.size(); i++) {
    if (data[i] > max or i == 0) max = data[i];
  }
  return max;
}

template <typename T>
double VectorMedian(std::vector<T> data) {
  sort(data.begin(), data.end());
  if (data.size() % 2 == 0) {
    int mid = int((double)data.size() / 2);
    return (double)(data[mid - 1] + data[mid]) / 2;
  } else {
    int mid = int((double)data.size() / 2) - 1;
    return (double)data[mid];
  }
}

template <typename T>
double VectorMeanOmit(std::vector<T> data, unsigned int omit) {
  double count = 0;
  for (unsigned int i = 0; i < data.size(); i++) {
    if (i == omit) continue;
    count += (double)data[i];
  }
  count = count / (data.size() - 1);
  return count;
}

template <typename T>
double JackKnifeError(std::vector<T> samples) {
  std::vector<double> means;
  means.reserve(samples.size());
  for (unsigned int i = 0; i < samples.size(); i++) {
    means.push_back(VectorMeanOmit(samples, i));
  }

  //  for( unsigned int i=0; i<samples.size(); i++ ) {
  //  std::cout << samples[ i ] << "\t" << means[i] << std::endl;
  // }

  return sqrt(samples.size() - 1) * VectorSampleSD(means);
}

// Taken from:
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
template <typename T>
typename T::size_type GenLevDist(const T &source, const T &target,
                                 typename T::size_type insert_cost = 1,
                                 typename T::size_type delete_cost = 1,
                                 typename T::size_type replace_cost = 1) {
  if (source.size() > target.size()) {
    return GenLevDist(target, source, delete_cost, insert_cost, replace_cost);
  }

  using TSizeType = typename T::size_type;
  const TSizeType min_size = source.size(), max_size = target.size();
  std::vector<TSizeType> lev_dist(min_size + 1);

  lev_dist[0] = 0;
  for (TSizeType i = 1; i <= min_size; ++i) {
    lev_dist[i] = lev_dist[i - 1] + delete_cost;
  }

  for (TSizeType j = 1; j <= max_size; ++j) {
    TSizeType previous_diagonal = lev_dist[0], previous_diagonal_save;
    lev_dist[0] += insert_cost;

    for (TSizeType i = 1; i <= min_size; ++i) {
      previous_diagonal_save = lev_dist[i];
      if (source[i - 1] == target[j - 1]) {
        lev_dist[i] = previous_diagonal;
      } else {
        lev_dist[i] = std::min(
            std::min(lev_dist[i - 1] + delete_cost, lev_dist[i] + insert_cost),
            previous_diagonal + replace_cost);
      }
      previous_diagonal = previous_diagonal_save;
    }
  }

  return lev_dist[min_size];
}

// base conversion
void int_to_basek(int x, int *&String, int L, int base);
int basek_to_int(int *String, int L, int base);
void long_to_basek(long x, int *&String, int L, int base);
long basek_to_long(int *String, int L, int base);

// File utilities
void Load(const char *String, int *&array);
int Load_Nu(const char *String);

#endif
