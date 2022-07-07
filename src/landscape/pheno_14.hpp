#ifndef _pheno_14_h_included_
#define _pheno_14_h_included_

#include <fstream>
#include <iostream>

#include "base/utilities.hpp"

class Phenotype {
 public:
  int *grid;
  int Z;
  int *Egrid;
  int E;

  int Label;
  int Symmetry;
  double complexity;
  double fitness;

  double x_c;
  double y_c;

  Phenotype();
  ~Phenotype();
  Phenotype(const Phenotype &phenotype);
  void Init(const Phenotype &phenotype);
  // Assignment
  Phenotype &operator=(const Phenotype &phenotype);

  int DistanceInstance(Phenotype T);
  int Distance(Phenotype T, int max_distance);
  void Copy(Phenotype A);
  void Copy(Phenotype A, Phenotype B);
  void fLoad(const char *file_name);
  void fLoad(std::ifstream &data);
  void Rotate90();
  void ReflectX();
  void ReflectY();
  void ReflectXY();
  void ReflectNXY();
  void SetCorner(int X, int Y);
  bool PositionOccupied(int *array, int n, int x, int y);
  void MakePerimeter();
  void FindCentre();
  void PrintPic();
  void fPrintPic(std::ofstream &fp);
  void fPrint(std::ofstream &fp);
  void CalculateSymmetry();
  int Match(Phenotype T);
  int MatchNoRots(Phenotype T);
  int MatchRots(Phenotype T);
  double CompareOld(Phenotype T);
  double Compare(Phenotype T);
  double FitnessSize(const int SIZE);
  double Fitness(Phenotype T);
  double Fitness(Phenotype X, Phenotype Y, double alpha);
  void Print();
  void Clear();
  void Make1x1();
  void Make2x1();
  void MakeCross();
  void Make2x2();
  void Make3x1();
  void MakeOverlapSquare();
  void MakeCW();
  void Make2x3();
  void Make2xTCW();
  void MakeECross();
  void MakeNS();
  void MakeComplexityTestShape();
  void MakeLCHS();
  void MakeHCHS();
  void MakeRing();
  void MakeAsymmCross();
  void MakeBoxCross();
  void MakeU();
  void MakeEx2x2();
  void MakeLongerCross();
  void MakeDoubleS();
  void MakeBigCross();
  void Make4x4();
  void Make3x3();
  void Make3x3_1();
  void MakeCorner();
  void Make6x6();
  void Make8x8();
  void Make4x4_1();
  void Make4x4_2();
  void Make4Rings();
  void MakeGrid(int length, int height);
  void MakeBigRing();
  void MakeEx4x2();
  void MakeSnake();
  void MakeT1();
  void MakeT2();
};

#endif
