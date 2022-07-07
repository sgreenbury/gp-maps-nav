#include "fast_assemble_14.hpp"

// Assembly code from:
// https://github.com/StochasticBiology/polyomino-evolution/blob/main/assembly.c

using namespace std;

// Function for adding new positions to the available position
inline int offset2(int x, int y, int offset) { return offset * y + x; }

// Wrapper for assembly to construct a polyomino phenotype object
Phenotype Assemble(int *&seq, int nt, const int MAX_TESTS) {
  int size = 0;
  int ARR = 16;
  int allownd = 1;
  int allowub = 0;

  int *Grid = 0;
  Grid = new int[ARR * ARR];

  int response =
      Grow(seq, nt, Grid, &size, allownd, allowub, MAX_TESTS, 0, ARR);
  Phenotype A;

  if (response >= 0) {
    A.Z = size;
    A.grid = new int[2 * A.Z];
    int count = 0;
    for (int y = 0; y < ARR; y++) {
      for (int x = 0; x < ARR; x++) {
        if (Grid[o2(x, y, ARR)] != -1) {
          A.grid[o2(0, count, 2)] = x;
          A.grid[o2(1, count, 2)] = y;
          count++;
        }
      }
    }
  } else {
    A.Z = 0;
  }
  delete[] Grid;
  A.SetCorner(0, 0);

  return A;
}

// outputs grid to screen
void OutputGrid(int *G, int ARR) {
  int x, y;

  for (y = 0; y < ARR; y++) {
    for (x = 0; x < ARR; x++)
      printf("%c", (G[y * ARR + x] == -1 ? '.' : G[y * ARR + x] + 48));
    printf("\n");
  }
}

// outputs grid to file
void FileOutputGrid(FILE *fp, int *G, int ARR) {
  int x, y;

  for (y = 0; y < ARR; y++) {
    for (x = 0; x < ARR; x++)
      fprintf(fp, "%c", (G[y * ARR + x] == -1 ? '.' : G[y * ARR + x] + 48));
    fprintf(fp, "\n");
  }
}

// neighbourlist will store data in the form:
// x y type ori

// add a new possible move (data) to neighbourlist
void Push(int *data, int *neighbourlist, int *length) {
  int i;
  for (i = 0; i < 4; i++) neighbourlist[4 * (*length) + i] = data[i];
  (*length)++;
}

// remove a possible move from neighbourlist
void Pop(int *neighbourlist, int *length, int index) {
  int i;
  (*length)--;
  if (index < (*length))
    for (i = 0; i < 4; i++)
      neighbourlist[i + 4 * index] = neighbourlist[((*length) * 4) + i];
}

// this function adds a block to the assembly grid
// given the move description x, y, block type, orientation
// and updates neighbourlist appropriately
int Add(int x, int y, int type, int ori, int *Grid, int *sides, int n,
        int *sidelist, int *neighbourlist, int *length, int allownd,
        int allowub, int ARR) {
  int i;
  int bond, partner;
  int data[4];
  int nondet;

  // add this tile to the grid
  Grid[y * ARR + x] = type;
  nondet = 0;
  // this compares our proposed move to all other current possible moves
  // to check for non-determinism and remove all possible moves in the
  // space we've chosen
  for (i = 0; i < *length; i++) {
    if (neighbourlist[4 * i] == x && neighbourlist[4 * i + 1] == y) {
      if (allownd == 0 && (neighbourlist[4 * i + 2] != type ||
                           neighbourlist[4 * i + 3] != ori)) {
        return -1;
      }
      Pop(neighbourlist, length, i);
      i--;
    }
  }

  // the subsequent four sections compute the new possible moves that we've
  // gained on each side of the new block
  // and push these possible moves into neighbourlist
  if (x > 0) {
    bond = sides[4 * type + (7 - ori) % 4];
    if (bond != 0 && Grid[y * ARR + x - 1] == -1) {
      partner = (bond % 2 == 0 ? bond - 1 : bond + 1);
      for (i = 0; i < sidelist[partner * (4 * n * 2 + 1)]; i++) {
        data[0] = x - 1;
        data[1] = y;
        data[2] = sidelist[partner * (4 * n * 2 + 1) + 2 * i + 1];
        data[3] = (5 - sidelist[partner * (4 * n * 2 + 1) + 2 * i + 2]) % 4;
        Push(data, neighbourlist, length);
      }
    }
  }

  if (x < ARR - 1) {
    bond = sides[4 * type + (5 - ori) % 4];
    if (bond != 0 && Grid[y * ARR + x + 1] == -1) {
      partner = (bond % 2 == 0 ? bond - 1 : bond + 1);
      for (i = 0; i < sidelist[partner * (4 * n * 2 + 1)]; i++) {
        data[0] = x + 1;
        data[1] = y;
        data[2] = sidelist[partner * (4 * n * 2 + 1) + 2 * i + 1];
        data[3] = (3 - sidelist[partner * (4 * n * 2 + 1) + 2 * i + 2]) % 4;
        Push(data, neighbourlist, length);
      }
    }
  }
  if (y > 0) {
    bond = sides[4 * type + (4 - ori) % 4];
    if (bond != 0 && Grid[(y - 1) * ARR + x] == -1) {
      partner = (bond % 2 == 0 ? bond - 1 : bond + 1);
      for (i = 0; i < sidelist[partner * (4 * n * 2 + 1)]; i++) {
        data[0] = x;
        data[1] = y - 1;
        data[2] = sidelist[partner * (4 * n * 2 + 1) + 2 * i + 1];
        data[3] = (6 - sidelist[partner * (4 * n * 2 + 1) + 2 * i + 2]) % 4;
        Push(data, neighbourlist, length);
      }
    }
  }
  if (y < ARR - 1) {
    bond = sides[4 * type + (6 - ori) % 4];
    if (bond != 0 && Grid[(y + 1) * ARR + x] == -1) {
      partner = (bond % 2 == 0 ? bond - 1 : bond + 1);
      for (i = 0; i < sidelist[partner * (4 * n * 2 + 1)]; i++) {
        data[0] = x;
        data[1] = y + 1;
        data[2] = sidelist[partner * (4 * n * 2 + 1) + 2 * i + 1];
        data[3] = (4 - sidelist[partner * (4 * n * 2 + 1) + 2 * i + 2]) % 4;
        Push(data, neighbourlist, length);
      }
    }
  }
  // check if edge of grid is reached
  if ((x == 0 || x == ARR - 1 || y == 0 || y == ARR - 1) && allowub == 0)
    return -2;

  return 0;
}

// more efficient than calling pow
int mypow2(int n) {
  int v = 1;
  int i;

  for (i = 0; i < n; i++) v *= 2;
  return v;
}

// convert a bitstring to an integer ruleset
void Convert(int *p, int *q, int NTILE, int NBITCOL) {
  int i, j, k;

  for (i = 0; i < NTILE; i++) {
    for (j = 0; j < 4; j++) {
      q[4 * i + j] = 0;
      for (k = 0; k < NBITCOL; k++) {
        q[4 * i + j] +=
            (p[NBITCOL * 4 * i + NBITCOL * j + k] == 1 ? mypow2(NBITCOL - 1 - k)
                                                       : 0);
      }
    }
  }
}

// grow a polyomino
int Grow(int *sides, int n, int *Grid, int *size, int allownd, int allowub,
         int checksize, int display, int ARR) {
  int i, j, k;
  int max = 0;
  int *sidelist;
  int *neighbourlist;
  int length;
  int ref;
  int x, y, type, ori;
  int *done;
  int run, last;

  for (i = 0; i < n * 4; i++) {
    if (sides[i] > max) max = sides[i];
  }

  max++;

  // allocate memory for structures storing the sides of placed blocks and the
  // neighbour list of possible moves we can make at any given time
  sidelist = (int *)malloc(sizeof(int) * (4 * n * 2 + 1) * (max + 1));
  neighbourlist = (int *)calloc(4 * ARR * ARR * n, sizeof(int));

  for (i = 0; i < (4 * n * 2 + 1) * (max + 1); i++)
    sidelist[i] = (i % (4 * n * 2 + 1) == 0 ? 0 : -1);

  // horrible looking piece of code to populate list of sides i.e. the
  // (tile,side) coordinates for every instance of colour 1, colour 2, ...
  for (i = 0; i < n; i++) {
    for (j = 0; j < 4; j++) {
      if (sides[4 * i + j] != 0) {
        for (k = 0;; k += 2) {
          if (sidelist[sides[4 * i + j] * (4 * 2 * n + 1) + k + 1] == -1) {
            sidelist[sides[4 * i + j] * (4 * 2 * n + 1) + k + 1] = i;
            sidelist[sides[4 * i + j] * (4 * 2 * n + 1) + k + 2] = j;
            sidelist[sides[4 * i + j] * (4 * 2 * n + 1)]++;
            break;
          }
        }
      }
    }
  }

  // simulate as many assembly runs as we are using to check determinism
  for (run = 0; run <= checksize; run++) {
    // initialise assembly grid with empty space
    for (i = 0; i < ARR * ARR; i++) Grid[i] = -1;

    length = 0;

    // add the first tile; this will update the neighbourlist
    Add(ARR / 2, ARR / 2, 0, 0, Grid, sides, n, sidelist, neighbourlist,
        &length, allownd, allowub, ARR);
    *size = 1;

    // as long as the length of neighbourlist is nonzero, pick a random element
    // and use this as the next addition move
    while (length != 0) {
      ref = (int)(rand() % length);
      x = neighbourlist[ref * 4];
      y = neighbourlist[ref * 4 + 1];
      type = neighbourlist[ref * 4 + 2];
      ori = neighbourlist[ref * 4 + 3];

      // add the chosen tile and update neighbourlist appropriately
      ref = Add(x, y, type, ori, Grid, sides, n, sidelist, neighbourlist,
                &length, allownd, allowub, ARR);
      if (ref < 0) {
        free(sidelist);
        free(neighbourlist);
        return ref;
      }
      (*size)++;
    }

    // catch nondeterminism
    if (checksize && run > 0 && *size != last) {
      free(sidelist);
      free(neighbourlist);
      return -3;
    }

    // record size for checking in subsequent assembly tests
    last = *size;
  }

  if (display) OutputGrid(Grid, ARR);

  free(sidelist);
  free(neighbourlist);
  return 0;
}
