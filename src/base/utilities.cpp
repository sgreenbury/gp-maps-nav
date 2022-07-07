#include "utilities.hpp"

using namespace std;

std::string MakeString(int *s, int n) {
  std::stringstream ss;
  for (int i = 0; i < n; i++) {
    ss << s[i];
    ss << ' ';
  }
  return ss.str();
}

std::string Number2String(int number) {
  stringstream ss;  // create a stringstream
  ss << number;     // add number to the stream
  return ss.str();  // return a string with the contents of the stream
}

int String2Number(const std::string &text) {
  stringstream ss(text);
  int result;
  ss >> result;
  return result;
}

int Search(int *array, int v, int l, int r) {
  while (r >= l) {
    int m = (int)((double)(l + r) / 2);
    if (v == array[m]) return m;
    if (v < array[m])
      r = m - 1;
    else
      l = m + 1;
  }
  return -1;
}

int RelativeCoordinatePosition(int x1, int y1, int x2, int y2) {
  if (y2 - y1 == 1 and x2 - x1 == 0) return 0;
  if (y2 - y1 == 0 and x2 - x1 == 1) return 1;
  if (y2 - y1 == -1 and x2 - x1 == 0) return 2;
  if (y2 - y1 == 0 and x2 - x1 == -1) return 3;

  return -1;
}

int Round(double x) {
  if (x < 0)
    return int(x - 0.5);
  else
    return int(x + 0.5);
}

// Maths Utilities
/*
int Mod(int x, int y)
{
  int result = x % y;
  if( result < 0 )
    return result + y;
  return result;
}
*/
double my_log(double base, double y) { return log(y) / log(base); }

// Base conversion utilities
void int_to_basek(int x, int *&String, int L, int base) {
  for (int j = 0; j < L; j++) String[j] = 0;

  int remainder = x;
  int remainder_n = x;
  double div;
  for (int j = L - 1; j >= 0; j--) {
    remainder_n = remainder % (int)(pow(base, j));
    div = (int)(remainder / pow(base, j));
    if (remainder_n == remainder) continue;

    if (remainder_n == 0) {
      if (remainder > 0) {
        String[j] = div;
        break;
      } else
        break;
    }

    remainder = remainder_n;
    String[j] = div;
  }
}

void long_to_basek(long x, int *&String, int L, int base) {
  for (int j = 0; j < L; j++) String[j] = 0;

  long remainder = x;
  long remainder_n = x;
  int div;
  for (int j = L - 1; j >= 0; j--) {
    remainder_n = remainder % (long)(pow(base, j));
    div = (int)(remainder / pow(base, j));
    if (remainder_n == remainder) continue;

    if (remainder_n == 0) {
      if (remainder > 0) {
        String[j] = div;
        break;
      } else
        break;
    }

    remainder = remainder_n;
    String[j] = div;
  }
}

long basek_to_long(int *String, int L, int base) {
  long count = 0;
  for (int i = 0; i < L; i++) count += (long)pow(base, i) * String[i];
  return count;
}

int Delta(int i, int j) {
  if (i == j)
    return 1;
  else
    return 0;
  return 0;
}

int basek_to_int(int *String, int L, int base) {
  int count = 0;
  for (int i = 0; i < L; i++) count += (int)pow(base, i) * String[i];
  return count;
}

// File Utilities
void Load(const char *String, int *&array) {
  ifstream fp;
  fp.open(String);
  int n = 0;

  char *ch;
  ch = new char[10];
  for (int i = 0; i < 10; i++) ch[i] = '\0';
  int x = 0;

  while (!fp.eof()) {
    char ch1 = fp.get();
    if (ch1 == '\n') {
      int t = atoi(ch);
      array[n] = t;
      n++;
      for (int i = 0; i < 10; i++) ch[i] = '\0';
      x = 0;
      continue;
    }

    if (fp.eof()) break;

    ch[x] = ch1;
    x++;
  }
  delete[] ch;
  fp.close();
  return;
}

int Load_Nu(const char *String) {
  ifstream fp;
  fp.open(String);
  int n = 0;
  while (!fp.eof()) {
    string String1;
    getline(fp, String1);
    if (fp.eof()) break;
    n++;
  }

  fp.close();
  return n;
}
