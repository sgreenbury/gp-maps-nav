#include "pheno_14.hpp"

using namespace std;

Phenotype::Phenotype() {
  grid = 0;
  Egrid = 0;
  Z = 0;
  E = 0;
  Label = 0;
  Symmetry = -1;
  complexity = -1;
  fitness = 0.0;
  x_c = 0.0;
  y_c = 0.0;
}

// New
Phenotype::~Phenotype() { Clear(); }

void Phenotype::Init(const Phenotype &phenotype) {
  Z = phenotype.Z;
  E = phenotype.E;

  Label = phenotype.Label;
  Symmetry = phenotype.Symmetry;
  complexity = phenotype.complexity;
  fitness = phenotype.fitness;
  x_c = phenotype.x_c;
  y_c = phenotype.y_c;

  grid = new int[2 * Z];
  for (int i = 0; i < 2 * Z; i++) {
    grid[i] = phenotype.grid[i];
  }

  Egrid = new int[2 * E];
  for (int i = 0; i < 2 * E; i++) {
    Egrid[i] = phenotype.Egrid[i];
  }
}

// Assignment
Phenotype::Phenotype(const Phenotype &phenotype) { Init(phenotype); }

Phenotype &Phenotype::operator=(const Phenotype &phenotype) {
  if (&phenotype != this) {
    Clear();
    Init(phenotype);
  }
  return *this;
}

void Phenotype::fLoad(std::ifstream &data) {
  std::string line;
  int count = 0;
  while (std::getline(data, line)) {
    count++;
  }
  Z = count;
  grid = new int[2 * Z];
  data.clear();
  data.seekg(0, ios::beg);

  count = 0;
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
      // Convert cell to an int, store in grid.
      int number = String2Number(cell);
      grid[count] = number;
      count++;
    }
  }

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
}

void Phenotype::fLoad(const char *file_name) {
  std::ifstream data(file_name);
  std::string line;
  int count = 0;
  while (std::getline(data, line)) {
    count++;
  }
  Z = count;
  grid = new int[2 * Z];
  data.close();

  data.open(file_name);
  count = 0;
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
      // Convert cell to an int, store in grid.
      int number = String2Number(cell);
      grid[count] = number;
      count++;
    }
  }
  data.close();

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
}

// Write phenotype to display
void Phenotype::PrintPic() {
  int x_low = -1000000;
  int x_high = -1000000;
  int y_low = -1000000;
  int y_high = -1000000;

  for (int j = 0; j < Z; j++) {
    if (grid[o2(0, j, 2)] < x_low or x_low == -1000000)
      x_low = grid[o2(0, j, 2)];
    if (grid[o2(0, j, 2)] > x_high or x_high == -1000000)
      x_high = grid[o2(0, j, 2)];
    if (grid[o2(1, j, 2)] < y_low or y_low == -1000000)
      y_low = grid[o2(1, j, 2)];
    if (grid[o2(1, j, 2)] > y_high or y_high == -1000000)
      y_high = grid[o2(1, j, 2)];
  }

  int buffer = 1;
  int g_height = y_high - y_low + 1 + buffer;
  int g_width = x_high - x_low + 1 + buffer;
  int x_c = -x_low + buffer;
  int y_c = -y_low + buffer;

  // Make a grid from which to print
  std::vector<std::vector<int> > xy_grid;
  xy_grid.resize(g_height);
  for (unsigned int i = 0; i < xy_grid.size(); i++) {
    xy_grid[i].resize(g_width, -1);
  }

  for (int j = 0; j < Z; j++) {
    int x = grid[o2(0, j, 2)] + x_c;
    int y = grid[o2(1, j, 2)] + y_c;
    xy_grid[y][x] = j;
  }

  cout << "/-";
  for (int i = 0; i < g_width; i++) cout << "--";
  cout << "\\" << endl;
  cout << "|";
  for (int i = 0; i < g_width; i++) cout << "  ";
  cout << " |" << endl;
  cout << "|";
  for (int i = 0; i < g_width; i++) cout << "  ";
  cout << " |" << endl;

  for (int i = int(xy_grid.size()) - 1; i >= 0; i--) {
    cout << "|";
    int DO = 0;
    if (i > 0)
      DO = 1;
    else
      DO = 0;

    for (unsigned int j = 0; j < xy_grid[i].size(); j++) {
      if (xy_grid[i][j] == -1)
        cout << "  ";
      else {
        cout << "O";
        if (j < xy_grid[i].size() - 1)
          if (xy_grid[i][j + 1] != -1)
            cout << "-";
          else
            cout << " ";
        else
          cout << " ";
      }
    }
    cout << " |" << endl;

    if (DO == 1) {
      cout << "|";
      for (unsigned int j = 0; j < xy_grid[i].size(); j++) {
        if (xy_grid[i - 1][j] == -1)
          cout << "  ";
        else {
          if (xy_grid[i][j] != -1) {
            if (xy_grid[i - 1][j] != -1)
              cout << "| ";
            else
              cout << "  ";
          } else
            cout << "  ";
        }
      }
      cout << " |" << endl;
    }
  }
  cout << "\\-";
  for (int i = 0; i < g_width; i++) cout << "--";
  cout << "/" << endl;
  cout << endl;
}

// Write phenotype to file
void Phenotype::fPrintPic(std::ofstream &fp) {
  int x_low = -1000000;
  int x_high = -1000000;
  int y_low = -1000000;
  int y_high = -1000000;

  for (int j = 0; j < Z; j++) {
    if (grid[o2(0, j, 2)] < x_low or x_low == -1000000)
      x_low = grid[o2(0, j, 2)];
    if (grid[o2(0, j, 2)] > x_high or x_high == -1000000)
      x_high = grid[o2(0, j, 2)];
    if (grid[o2(1, j, 2)] < y_low or y_low == -1000000)
      y_low = grid[o2(1, j, 2)];
    if (grid[o2(1, j, 2)] > y_high or y_high == -1000000)
      y_high = grid[o2(1, j, 2)];
  }

  int buffer = 1;
  int g_height = y_high - y_low + 1 + buffer;
  int g_width = x_high - x_low + 1 + buffer;
  int x_c = -x_low + buffer;
  int y_c = -y_low + buffer;

  // Make a grid from which to print
  std::vector<std::vector<int> > xy_grid;
  xy_grid.resize(g_height);
  for (unsigned int i = 0; i < xy_grid.size(); i++) {
    xy_grid[i].resize(g_width, -1);
  }

  for (int j = 0; j < Z; j++) {
    int x = grid[o2(0, j, 2)] + x_c;
    int y = grid[o2(1, j, 2)] + y_c;
    xy_grid[y][x] = j;
  }

  fp << "/-";
  for (int i = 0; i < g_width; i++) fp << "--";
  fp << "\\" << endl;
  fp << "|";
  for (int i = 0; i < g_width; i++) fp << "  ";
  fp << " |" << endl;
  fp << "|";
  for (int i = 0; i < g_width; i++) fp << "  ";
  fp << " |" << endl;

  for (int i = int(xy_grid.size()) - 1; i >= 0; i--) {
    fp << "|";
    int DO = 0;
    if (i > 0)
      DO = 1;
    else
      DO = 0;

    for (unsigned int j = 0; j < xy_grid[i].size(); j++) {
      if (xy_grid[i][j] == -1)
        fp << "  ";
      else {
        fp << "O";
        if (j < xy_grid[i].size() - 1)
          if (xy_grid[i][j + 1] != -1)
            fp << "-";
          else
            fp << " ";
        else
          fp << " ";
      }
    }
    fp << " |" << endl;

    if (DO == 1) {
      fp << "|";
      for (unsigned int j = 0; j < xy_grid[i].size(); j++) {
        if (xy_grid[i - 1][j] == -1)
          fp << "  ";
        else {
          if (xy_grid[i][j] != -1) {
            if (xy_grid[i - 1][j] != -1)
              fp << "| ";
            else
              fp << "  ";
          } else
            fp << "  ";
        }
      }
      fp << " |" << endl;
    }
  }
  fp << "\\-";
  for (int i = 0; i < g_width; i++) fp << "--";
  fp << "/" << endl;
  fp << endl;
}

void Phenotype::Copy(Phenotype A) {
  Clear();
  Z = A.Z;
  E = A.E;

  Label = A.Label;
  Symmetry = A.Symmetry;
  complexity = A.complexity;
  fitness = A.fitness;
  x_c = A.x_c;
  y_c = A.y_c;

  if (grid != 0) delete[] grid;

  grid = new int[2 * Z];

  for (int i = 0; i < 2 * Z; i++) {
    grid[i] = A.grid[i];
  }

  if (Egrid != 0) delete[] Egrid;

  Egrid = new int[2 * E];
  for (int i = 0; i < 2 * E; i++) {
    Egrid[i] = A.Egrid[i];
  }
}

void Phenotype::Copy(Phenotype A, Phenotype B) {
  Z = A.Z + B.Z;

  if (grid != 0) delete[] grid;

  grid = new int[2 * Z];

  for (int i = 0; i < 2 * A.Z; i++) {
    grid[i] = A.grid[i];
  }

  int x_hi = 0;
  int y_hi = 0;
  for (int i = 0; i < A.Z; i++) {
    if (grid[o2(0, i, 2)] > x_hi) x_hi = grid[o2(0, i, 2)];
    if (grid[o2(1, i, 2)] > y_hi) y_hi = grid[o2(1, i, 2)];
  }

  x_hi += 2;
  y_hi += 2;

  B.SetCorner(x_hi, y_hi);

  for (int i = 0; i < B.Z; i++) {
    grid[o2(0, A.Z + i, 2)] = B.grid[o2(0, i, 2)];
    grid[o2(1, A.Z + i, 2)] = B.grid[o2(1, i, 2)];
  }

  if (Egrid != 0) delete[] Egrid;

  MakePerimeter();
  FindCentre();
}

// Sets left (bottom) corner of phenotype A, at (X,Y)
void Phenotype::SetCorner(int X, int Y) {
  int xa = 0, ya = 0;
  for (int i = 0; i < Z; i++) {
    if (grid[o2(0, i, 2)] < xa or i == 0) {
      xa = grid[o2(0, i, 2)];
      ya = grid[o2(1, i, 2)];
      continue;
    }
    if (grid[o2(0, i, 2)] == xa and grid[o2(1, i, 2)] < ya)
      ya = grid[o2(1, i, 2)];
  }
  for (int i = 0; i < Z; i++) {
    grid[o2(0, i, 2)] = grid[o2(0, i, 2)] - xa + X;
    grid[o2(1, i, 2)] = grid[o2(1, i, 2)] - ya + Y;
  }

  for (int i = 0; i < E; i++) {
    Egrid[o2(0, i, 2)] = Egrid[o2(0, i, 2)] - xa + X;
    Egrid[o2(1, i, 2)] = Egrid[o2(1, i, 2)] - ya + Y;
  }

  return;
}

void Phenotype::Rotate90() {
  for (int i = 0; i < Z; i++) {
    int x = grid[o2(0, i, 2)];
    int y = grid[o2(1, i, 2)];
    grid[o2(0, i, 2)] = -y;
    grid[o2(1, i, 2)] = x;
  }

  for (int i = 0; i < E; i++) {
    int x = Egrid[o2(0, i, 2)];
    int y = Egrid[o2(1, i, 2)];
    Egrid[o2(0, i, 2)] = -y;
    Egrid[o2(1, i, 2)] = x;
  }

  return;
}

void Phenotype::ReflectX() {
  for (int i = 0; i < Z; i++) {
    int y = grid[o2(1, i, 2)];
    // grid[o2(0, i, 2)] = x;
    grid[o2(1, i, 2)] = -y;
  }

  for (int i = 0; i < E; i++) {
    int y = Egrid[o2(1, i, 2)];
    // Egrid[o2(0, i, 2)] = x;
    Egrid[o2(1, i, 2)] = -y;
  }

  return;
}

void Phenotype::ReflectY() {
  for (int i = 0; i < Z; i++) {
    int x = grid[o2(0, i, 2)];
    grid[o2(0, i, 2)] = -x;
    // grid[o2(1, i, 2)] = y;
  }

  for (int i = 0; i < E; i++) {
    int x = Egrid[o2(0, i, 2)];
    Egrid[o2(0, i, 2)] = -x;
    // Egrid[o2(1, i, 2)] = y;
  }

  return;
}

void Phenotype::ReflectXY() {
  for (int i = 0; i < Z; i++) {
    int x = grid[o2(0, i, 2)];
    int y = grid[o2(1, i, 2)];
    grid[o2(0, i, 2)] = y;
    grid[o2(1, i, 2)] = x;
  }

  for (int i = 0; i < E; i++) {
    int x = Egrid[o2(0, i, 2)];
    int y = Egrid[o2(1, i, 2)];
    Egrid[o2(0, i, 2)] = y;
    Egrid[o2(1, i, 2)] = x;
  }

  return;
}

void Phenotype::ReflectNXY() {
  for (int i = 0; i < Z; i++) {
    int x = grid[o2(0, i, 2)];
    int y = grid[o2(1, i, 2)];
    grid[o2(0, i, 2)] = -y;
    grid[o2(1, i, 2)] = -x;
  }

  for (int i = 0; i < E; i++) {
    int x = Egrid[o2(0, i, 2)];
    int y = Egrid[o2(1, i, 2)];
    Egrid[o2(0, i, 2)] = -y;
    Egrid[o2(1, i, 2)] = -x;
  }

  return;
}

void Phenotype::fPrint(std::ofstream &fp) {
  for (int i = 0; i < Z; i++)
    fp << grid[o2(0, i, 2)] << "," << grid[o2(1, i, 2)] << "\n";
}

void Phenotype::FindCentre() {
  int low_x = 0;
  int high_x = 0;
  int low_y = 0;
  int high_y = 0;
  for (int i = 0; i < Z; i++) {
    if (grid[o2(0, i, 2)] < low_x) low_x = grid[o2(0, i, 2)];
    if (grid[o2(1, i, 2)] < low_y) low_y = grid[o2(1, i, 2)];
    if (grid[o2(0, i, 2)] > high_x) high_x = grid[o2(0, i, 2)];
    if (grid[o2(1, i, 2)] > high_y) high_y = grid[o2(1, i, 2)];
  }
  x_c = 0.5 * (high_x + low_x);
  y_c = 0.5 * (high_y + low_y);
}

void Phenotype::CalculateSymmetry() {
  if (Z == 0) {
    Symmetry = 6;
    return;
  }

  int r_count = 0;
  int s_count = 0;

  Phenotype A;
  A.Copy(*this);

  for (int i = 0; i < 8; i++) {
    switch (i) {
      case 0:
      case 1:
      case 2:
      case 3:
        Rotate90();
        SetCorner(0, 0);
        r_count += Match(A);
        break;
      case 4:
        ReflectX();
        SetCorner(0, 0);
        s_count += Match(A);
        ReflectX();
        break;
      case 5:
        ReflectY();
        SetCorner(0, 0);
        s_count += Match(A);
        ReflectY();
        break;
      case 6:
        ReflectXY();
        SetCorner(0, 0);
        s_count += Match(A);
        ReflectXY();
        break;
      case 7:
        ReflectNXY();
        SetCorner(0, 0);
        s_count += Match(A);
        ReflectNXY();
        break;
    }
  }
  SetCorner(0, 0);
  FindCentre();
  A.Clear();

  if (r_count == 4 && s_count == 4) Symmetry = 0;
  if (r_count == 4 && s_count != 4) Symmetry = 1;
  if (r_count == 2 && s_count == 2) Symmetry = 2;
  if (r_count == 2 && s_count != 2) Symmetry = 3;
  if (r_count == 1 && s_count == 1) Symmetry = 4;
  if (r_count == 1 && s_count != 1) Symmetry = 5;
}

bool Phenotype::PositionOccupied(int *array, int n, int X, int Y) {
  for (int i = 0; i < n; i++)
    if (array[o2(0, i, 2)] == X and array[o2(1, i, 2)] == Y) return 1;
  return 0;
}

void Phenotype::MakePerimeter() {
  int *rot_coords;
  rot_coords = new int[8];
  rot_coords[0] = 0;
  rot_coords[1] = 1;
  rot_coords[2] = 1;
  rot_coords[3] = 0;
  rot_coords[4] = 0;
  rot_coords[5] = -1;
  rot_coords[6] = -1;
  rot_coords[7] = 0;

  resize_array(Egrid, E, 4 * 2 * Z);

  for (int i = 0; i < Z; i++)
    for (int j = 0; j < 4; j++) {
      int xt = grid[o2(0, i, 2)] + rot_coords[o2(0, j, 2)];
      int yt = grid[o2(1, i, 2)] + rot_coords[o2(1, j, 2)];
      if (PositionOccupied(grid, Z, xt, yt) == 0 and
          PositionOccupied(Egrid, E, xt, yt) == 0) {
        Egrid[o2(0, E, 2)] = xt;
        Egrid[o2(1, E, 2)] = yt;
        E++;
      }
    }
  delete[] rot_coords;

  return;
}

int Phenotype::Match(Phenotype T) {
  if (T.Z != Z) return 0;

  int count = 0;
  for (int i = 0; i < T.Z; i++) {
    int x = T.grid[o2(0, i, 2)];
    int y = T.grid[o2(1, i, 2)];
    count += PositionOccupied(grid, Z, x, y);
  }
  return count == T.Z;
}

int Phenotype::MatchNoRots(Phenotype T) {
  if (T.Z != Z) return 0;
  SetCorner(0, 0);
  if (Match(T) == 1) return 1;
  return 0;
}

int Phenotype::MatchRots(Phenotype T) {
  if (T.Z != Z) return 0;
  if (T.Z == 0 and Z == 0) return 1;
  for (int i = 0; i < 4; i++) {
    Rotate90();
    SetCorner(0, 0);
    if (Match(T) == 1) return 1;
  }
  return 0;
}

double Phenotype::CompareOld(Phenotype T) {
  int count = 0;
  for (int i = 0; i < T.Z; i++) {
    int x = T.grid[o2(0, i, 2)];
    int y = T.grid[o2(1, i, 2)];
    count += PositionOccupied(grid, Z, x, y);
  }

  for (int i = 0; i < T.E; i++) {
    int x = T.Egrid[o2(0, i, 2)];
    int y = T.Egrid[o2(1, i, 2)];
    count += PositionOccupied(Egrid, E, x, y);
  }

  return (double)count / (T.Z + T.E);
}

double Phenotype::Compare(Phenotype T) {
  // Count the number of differences
  int differences = 0;
  for (int i = 0; i < Z; i++) {
    int x = grid[o2(0, i, 2)];
    int y = grid[o2(1, i, 2)];
    if (PositionOccupied(T.grid, T.Z, x, y) == 0) differences += 1;
  }

  if (T.Z > Z) differences += T.Z - Z;

  return (double)1 / (1 + differences);
}

int Phenotype::DistanceInstance(Phenotype T) {
  // Count the number of differences
  int differences = 0;
  for (int i = 0; i < Z; i++) {
    int x = grid[o2(0, i, 2)];
    int y = grid[o2(1, i, 2)];
    if (PositionOccupied(T.grid, T.Z, x, y) == 0) {
      differences += 1;
    }
  }

  if (T.Z > Z) {
    differences += T.Z - Z;
  }

  return differences;
}

int Phenotype::Distance(Phenotype T, int max_distance) {
  int differences = max_distance;
  if (Z == 0 or T.Z == 0) {
    if (Z == 0 and T.Z == 0) {
      return 0;
    } else {
      return differences;
    }
  }
  for (int j = 0; j < 4; j++) {
    Rotate90();
    for (int i = 0; i < T.Z; i++) {
      int X = T.grid[o2(0, i, 2)];
      int Y = T.grid[o2(1, i, 2)];
      SetCorner(X, Y);
      int differences_new = DistanceInstance(T);
      if (differences_new < differences) {
        differences = differences_new;
      }
    }
  }

  SetCorner(0, 0);
  return differences;
}

double Phenotype::FitnessSize(const int SIZE) {
  if (Z < SIZE)
    return (double)Z / SIZE;
  else
    return (double)1;
}

double Phenotype::Fitness(Phenotype T) {
  double fit = 0.0;
  for (int j = 0; j < 4; j++) {
    Rotate90();
    for (int i = 0; i < T.Z; i++) {
      int X = T.grid[o2(0, i, 2)];
      int Y = T.grid[o2(1, i, 2)];
      SetCorner(X, Y);
      double fit_new = Compare(T);
      if (fit_new > fit) {
        fit = fit_new;
      }
    }
  }
  SetCorner(0, 0);
  return fit;
}

double Phenotype::Fitness(Phenotype X, Phenotype Y, double alpha) {
  double fit_X = Fitness(X);
  double fit_Y = Fitness(Y);
  double fit = 0.0;
  if (fit_Y > alpha * fit_X)
    fit = fit_Y;
  else
    fit = alpha * fit_X;

  return fit;
}

void Phenotype::Print() {
  cout << "Blocks" << endl;
  for (int i = 0; i < Z; i++)
    cout << grid[o2(0, i, 2)] << "," << grid[o2(1, i, 2)] << endl;
  cout << endl;
  if (E > 0) {
    cout << "Perimeter" << endl;
    for (int i = 0; i < E; i++)
      cout << Egrid[o2(0, i, 2)] << "," << Egrid[o2(1, i, 2)] << endl;
    cout << endl;
  }
  return;
}

void Phenotype::Clear() {
  Z = 0;
  E = 0;
  if (grid != 0) {
    delete[] grid;
    grid = 0;
  }
  if (Egrid != 0) {
    delete[] Egrid;
    Egrid = 0;
  }
  return;
}

void Phenotype::MakeCross() {
  Z = 5;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 1;
  grid[9] = 1;
  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeAsymmCross() {
  Z = 7;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 1;
  grid[9] = 1;
  grid[10] = 1;
  grid[11] = -2;
  grid[12] = 1;
  grid[13] = 2;
  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeBoxCross() {
  Z = 6;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 1;
  grid[9] = 1;
  grid[10] = 1;
  grid[11] = -2;
  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::Make2x2() {
  Z = 4;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 0;
  grid[5] = 1;
  grid[6] = 1;
  grid[7] = 1;
  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::Make3x1() {
  Z = 3;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeCorner() {
  Z = 3;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeOverlapSquare() {
  Z = 7;

  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 2;
  grid[9] = 0;
  grid[10] = 0;
  grid[11] = -1;
  grid[12] = 2;
  grid[13] = 1;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeCW() {
  Z = 8;

  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  grid[6] = 0;
  grid[7] = 1;
  grid[8] = 0;
  grid[9] = 2;
  grid[10] = 2;
  grid[11] = 1;
  grid[12] = 1;
  grid[13] = -1;
  grid[14] = -1;
  grid[15] = 0;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make2x3() {
  Z = 6;

  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 0;
  grid[5] = 1;
  grid[6] = 1;
  grid[7] = 1;
  grid[8] = 0;
  grid[9] = 2;
  grid[10] = 1;
  grid[11] = 2;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make2xTCW() {
  Z = 6;

  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  grid[6] = 0;
  grid[7] = 1;
  grid[8] = 0;
  grid[9] = 2;
  grid[10] = 1;
  grid[11] = -1;

  // O
  // |
  // O - O
  // |   |
  // O - o
  //     |
  //     O

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeECross() {
  Z = 6;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 2;
  grid[9] = 0;
  grid[10] = 3;
  grid[11] = 0;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeNS() {
  Z = 7;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 2;
  grid[9] = 0;
  grid[10] = 3;
  grid[11] = 0;
  grid[12] = 3;
  grid[13] = 1;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeComplexityTestShape() {
  Z = 18;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 1;
  grid[4] = 1;
  grid[5] = 0;
  grid[6] = 1;
  grid[7] = -1;
  grid[8] = 2;
  grid[9] = 0;
  grid[10] = 3;
  grid[11] = 0;
  grid[12] = 3;
  grid[13] = 1;
  grid[14] = 3;
  grid[15] = -1;
  grid[16] = 4;
  grid[17] = 1;
  grid[18] = 4;
  grid[19] = -1;
  grid[20] = 5;
  grid[21] = 0;
  grid[22] = 5;
  grid[23] = 1;
  grid[24] = 5;
  grid[25] = -1;
  grid[26] = 6;
  grid[27] = 0;
  grid[28] = 7;
  grid[29] = 0;
  grid[30] = 7;
  grid[31] = 1;
  grid[32] = 7;
  grid[33] = -1;
  grid[34] = 8;
  grid[35] = 0;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make1x1() {
  Z = 1;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make2x1() {
  Z = 2;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 0;
  grid[3] = 1;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeLCHS() {
  Z = 17;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 3;
  grid[5] = 0;
  grid[6] = 4;
  grid[7] = 0;
  grid[8] = 1;
  grid[9] = 1;
  grid[10] = 4;
  grid[11] = 1;
  grid[12] = 0;
  grid[13] = 2;
  grid[14] = 1;
  grid[15] = 2;
  grid[16] = 2;
  grid[17] = 2;
  grid[18] = 3;
  grid[19] = 2;
  grid[20] = 4;
  grid[21] = 2;
  grid[22] = 0;
  grid[23] = 3;
  grid[24] = 3;
  grid[25] = 3;
  grid[26] = 0;
  grid[27] = 4;
  grid[28] = 1;
  grid[29] = 4;
  grid[30] = 3;
  grid[31] = 4;
  grid[32] = 4;
  grid[33] = 4;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeHCHS() {
  Z = 17;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 3;
  grid[7] = 0;
  grid[8] = 4;
  grid[9] = 0;
  grid[10] = 2;
  grid[11] = 1;
  grid[12] = 4;
  grid[13] = 1;
  grid[14] = 1;
  grid[15] = 2;
  grid[16] = 2;
  grid[17] = 2;
  grid[18] = 3;
  grid[19] = 2;
  grid[20] = 0;
  grid[21] = 3;
  grid[22] = 1;
  grid[23] = 3;
  grid[24] = 4;
  grid[25] = 3;
  grid[26] = 1;
  grid[27] = 4;
  grid[28] = 2;
  grid[29] = 4;
  grid[30] = 3;
  grid[31] = 4;
  grid[32] = 4;
  grid[33] = 4;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeRing() {
  Z = 8;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 0;
  grid[7] = 1;
  grid[8] = 2;
  grid[9] = 1;
  grid[10] = 0;
  grid[11] = 2;
  grid[12] = 1;
  grid[13] = 2;
  grid[14] = 2;
  grid[15] = 2;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeU() {
  Z = 7;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 0;
  grid[7] = 1;
  grid[8] = 2;
  grid[9] = 1;
  grid[10] = 0;
  grid[11] = 2;
  grid[12] = 2;
  grid[13] = 2;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeEx2x2() {
  Z = 6;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 0;
  grid[5] = 1;
  grid[6] = 1;
  grid[7] = 1;
  grid[8] = 2;
  grid[9] = 0;
  grid[10] = 1;
  grid[11] = -1;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeLongerCross() {
  Z = 9;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 3;
  grid[7] = 0;
  grid[8] = 4;
  grid[9] = 0;
  grid[10] = 2;
  grid[11] = 1;
  grid[12] = 2;
  grid[13] = 2;
  grid[14] = 2;
  grid[15] = -1;
  grid[16] = 2;
  grid[17] = -2;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeDoubleS() {
  Z = 17;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 3;
  grid[5] = 0;
  grid[6] = 4;
  grid[7] = 0;
  grid[8] = 1;
  grid[9] = 1;
  grid[10] = 4;
  grid[11] = 1;
  grid[12] = 0;
  grid[13] = 2;
  grid[14] = 1;
  grid[15] = 2;
  grid[16] = 2;
  grid[17] = 2;
  grid[18] = 3;
  grid[19] = 2;
  grid[20] = 4;
  grid[21] = 2;
  grid[22] = 0;
  grid[23] = 3;
  grid[24] = 3;
  grid[25] = 3;
  grid[26] = 0;
  grid[27] = 4;
  grid[28] = 1;
  grid[29] = 4;
  grid[30] = 3;
  grid[31] = 4;
  grid[32] = 4;
  grid[33] = 4;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::MakeBigCross() {
  Z = 12;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 3;
  grid[7] = 0;
  grid[8] = 1;
  grid[9] = -1;
  grid[10] = 2;
  grid[11] = -1;
  grid[12] = 0;
  grid[13] = 1;
  grid[14] = 1;
  grid[15] = 1;
  grid[16] = 2;
  grid[17] = 1;
  grid[18] = 3;
  grid[19] = 1;
  grid[20] = 1;
  grid[21] = 2;
  grid[22] = 2;
  grid[23] = 2;

  Egrid = 0;
  MakePerimeter();

  return;
}

void Phenotype::Make4x4() {
  Z = 16;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 3;
  grid[7] = 0;
  grid[8] = 1;
  grid[9] = -1;
  grid[10] = 2;
  grid[11] = -1;
  grid[12] = 0;
  grid[13] = 1;
  grid[14] = 1;
  grid[15] = 1;
  grid[16] = 2;
  grid[17] = 1;
  grid[18] = 3;
  grid[19] = 1;
  grid[20] = 1;
  grid[21] = 2;
  grid[22] = 2;
  grid[23] = 2;
  grid[24] = 0;
  grid[25] = -1;
  grid[26] = 3;
  grid[27] = -1;
  grid[28] = 0;
  grid[29] = 2;
  grid[30] = 3;
  grid[31] = 2;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make3x3() {
  Z = 9;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 0;
  grid[7] = 1;
  grid[8] = 1;
  grid[9] = 1;
  grid[10] = 2;
  grid[11] = 1;
  grid[12] = 0;
  grid[13] = 2;
  grid[14] = 1;
  grid[15] = 2;
  grid[16] = 2;
  grid[17] = 2;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make3x3_1() {
  Z = 8;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 0;
  grid[7] = 1;
  grid[8] = 1;
  grid[9] = 1;
  grid[10] = 2;
  grid[11] = 1;
  grid[12] = 0;
  grid[13] = 2;
  grid[14] = 1;
  grid[15] = 2;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make6x6() {
  Z = 36;
  int count = 0;
  grid = new int[2 * Z];
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      grid[2 * count] = i;
      grid[2 * count + 1] = j;
      count++;
    }
  }

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make8x8() {
  Z = 64;
  int count = 0;
  grid = new int[2 * Z];
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      grid[2 * count] = i;
      grid[2 * count + 1] = j;
      count++;
    }
  }

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make4x4_1() {
  Z = 15;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 3;
  grid[7] = 0;
  grid[8] = 1;
  grid[9] = -1;
  grid[10] = 2;
  grid[11] = -1;
  grid[12] = 0;
  grid[13] = 1;
  grid[14] = 1;
  grid[15] = 1;
  grid[16] = 2;
  grid[17] = 1;
  grid[18] = 3;
  grid[19] = 1;
  grid[20] = 1;
  grid[21] = 2;
  grid[22] = 2;
  grid[23] = 2;
  grid[24] = 0;
  grid[25] = -1;
  grid[26] = 3;
  grid[27] = -1;
  grid[28] = 0;
  grid[29] = 2;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::Make4Rings() {
  Z = 32;
  grid = new int[2 * Z];
  int count = 0;
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      if ((i == 1 or i == 4) and (j == 1 or j == 4)) {
        continue;
      }
      grid[2 * count] = i;
      grid[2 * count + 1] = j;
      count++;
    }
  }

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeGrid(int length, int height) {
  Z = height * length;
  int count = 0;
  grid = new int[2 * Z];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < length; j++) {
      grid[2 * count] = i;
      grid[2 * count + 1] = j;
      count++;
    }
  }

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeBigRing() {
  Z = 14;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 2;
  grid[5] = 0;
  grid[6] = 3;
  grid[7] = 0;
  grid[8] = 0;
  grid[9] = 1;
  grid[10] = 3;
  grid[11] = 1;
  grid[12] = 0;
  grid[13] = 2;
  grid[14] = 3;
  grid[15] = 2;
  grid[16] = 0;
  grid[17] = 3;
  grid[18] = 1;
  grid[19] = 3;
  grid[20] = 2;
  grid[21] = 3;
  grid[22] = 3;
  grid[23] = 3;
  grid[24] = 1;
  grid[25] = 4;
  grid[26] = 2;
  grid[27] = 4;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeEx4x2() {
  Z = 10;
  int count = 0;
  grid = new int[2 * Z];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 4; j++) {
      grid[2 * count] = i;
      grid[2 * count + 1] = j;
      count++;
    }
  }

  grid[16] = 2;
  grid[17] = 1;
  grid[18] = 2;
  grid[19] = 2;
  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeSnake() {
  Z = 4;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 1;
  grid[3] = 0;
  grid[4] = 1;
  grid[5] = 1;
  grid[6] = 2;
  grid[7] = 1;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeT1() {
  Z = 9;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 0;
  grid[3] = 1;
  grid[4] = 0;
  grid[5] = 2;
  grid[6] = -1;
  grid[7] = 2;
  grid[8] = -2;
  grid[9] = 2;
  grid[10] = -2;
  grid[11] = 3;
  grid[12] = -2;
  grid[13] = 4;
  grid[14] = -1;
  grid[15] = 4;
  grid[16] = 0;
  grid[17] = 4;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}

void Phenotype::MakeT2() {
  Z = 9;
  grid = new int[2 * Z];
  grid[0] = 0;
  grid[1] = 0;
  grid[2] = 0;
  grid[3] = 1;
  grid[4] = 0;
  grid[5] = 2;
  grid[6] = 1;
  grid[7] = 2;
  grid[8] = 2;
  grid[9] = 2;
  grid[10] = 0;
  grid[11] = 3;
  grid[12] = 0;
  grid[13] = 4;
  grid[14] = 1;
  grid[15] = 4;
  grid[16] = 2;
  grid[17] = 4;

  Egrid = 0;
  MakePerimeter();
  SetCorner(0, 0);
  return;
}
