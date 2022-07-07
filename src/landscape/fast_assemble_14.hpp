#ifndef _fast_assemble_14_h_included_
#define _fast_assemble_14_h_included_

#include <math.h>
#include <iostream>

#include <cstdlib>
#include <fstream>

#include "base/utilities.hpp"
#include "pheno_14.hpp"

using namespace std;

// Phenotype Assemble(Genotype &G, System Test, const int MAX_TESTS);
Phenotype Assemble(int *&seq, int nt, const int MAX_TESTS);

// outputs grid to screen
void OutputGrid(int *G, int ARR);

// outputs grid to file
void FileOutputGrid(FILE *fp, int *G, int ARR);

// add data to neighbourlisty
void Push(int *data, int *neighbourlist, int *length);

// remove data from neighbourlist
void Pop(int *neighbourlist, int *length, int index);

int Add(int x, int y, int type, int ori, int *Grid, int *sides, int n,
        int *sidelist, int *neighbourlist, int *length, int allownd,
        int allowub, int ARR);

int mypow2(int n);

void Convert(int *p, int *q, int NTILE, int NBITCOL);

// grow a polyomino
int Grow(int *sides, int n, int *Grid, int *size, int allownd, int allowub,
         int checksize, int display, int ARR);

#endif
