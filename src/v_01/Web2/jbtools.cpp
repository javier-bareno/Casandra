#include<iostream>
#include<fstream>
#include "sputt.h"

using namespace std;

// Tests of sputt.cpp to check that python version is OK

sputt run_calc(void)
{
  // Run sputtering calculation with standard hard coded 
  // params and return sputt object 
  
  sputt system;
  CData *ansData;

  // Gas atom params
  string Gas_symbol;
  double MG, ZG, RG, UG, QG, DENG, KAPA;
  MG = 39.948;
  ZG = 18;
  RG = 0.98;
  UG = 5;
  QG = 1;
  DENG = 5;
  KAPA = 0.0002;

  // Target element
  string Target_symbol;
  double MA, ZA, RA, UO, QZ, DEN;
  MA = 28;
  ZA = 14;
  RA = 2;
  UO = 4.63;
  QZ = 0.75;
  DEN = 2.33;

  // other params
  double disch_volt, Eavg,  disch_cur, pres, T0;
  disch_volt = 100;
  Eavg = 10;
  disch_cur = 1.1;
  pres = 1.2;
  T0 = 300;

  int TC = 0;
  bool Tcorrect = TC == 1;
  double in_rad, out_rad, N1, N2, L5, dist;
  in_rad =20;
  out_rad =30;
  N1 = 10;
  N2 = 20;
  L5 = 2;
  dist = 150;

  system.MG = MG;
  system.ZG = ZG;
  system.RG = RG;
  system.UG = UG;
  system.QG = QG;
  system.DENG = DENG;
  system.KAPPA = KAPA;
  system.MA = MA;
  system.ZA = ZA;
  system.RA = RA;
  system.UA = UO;
  system.DENA = DEN;
  system.EI = disch_volt;
  system.E9 = Eavg;
  system.CUR = disch_cur;
  system.P = pres;
  system.T0  = T0;
  system.R1 =  in_rad;
  system.R2 = out_rad;
  system.N1 = N1;
  system.N2 = N2;
  system.L5 = L5;
  system.Z = dist;
  system.TCorr = Tcorrect;
  system.QA = QZ;

  ansData = system.calculate();
  
  return(system);
};

int _main(void)
{
    cout << "Try to run a sputtering calculation\n" ;
    sputt calc = run_calc();
    double Y = calc.Y;
    cout << "Found sput yield: " << Y;
    return 0;
}
