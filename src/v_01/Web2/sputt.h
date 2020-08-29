//---------------------------------------------------------------------------

#ifndef sputtH
#define sputtH

//#include "CPlot.h"
#include "CData.h"

#include <cmath>
#include<fstream>
#include<iostream>

#ifndef PI
  #define PI 3.1415927
#endif


class sputt
{
  public:
  //sputter gas
  double MG, ZG, RG, UG, QG, DENG, KAPPA;
  char IO[2];
  //target
  double MA, ZA, RA, UA, QA, DENA;
  char EL[2];
  // magnetron settings
  double EI, CUR, P, T0, R1, R2, N1, N2;
  bool TCorr;
  double Z; // distance target-substarte in mm
  // calculations and results
  double TG; // gas temperature during/after temp correction
  double E9; // average energy of sputtered gas
  double L5; // radius of measured volume

  //private: made all public to test for python version
  // parameters from calculations and results
  double LAMBDA; // mean free path
  double ETA; // ~ number of collisions prior thermalization
  double LALA; // lambda/(10*TG)
  double DEP1[100], DEP21[100], DEP22[100], DEP[100];
  double Y; // sputtering yield;
  double R99; // collisionless rate?
  double J[100][100]; //diff current
    int JN1, JN2;

  // functions translated from Ivan's prog
  public:
    
  CData *calculate();
    //sputt();
    sputt(string log_fname); // constructor for loging
    ~sputt(); // destructor

  private:
  ofstream logfile; //used for loging
    
  double correct_T(double T);
  double atom_radius(double z);
  double mean_free_path(double T);
  double avg_collision_number_to_thermalize(double T);
  void collisionless_transport();
  double collisionless_dep_rate();
  double sputtering_yield();
  void difusion_source();
  void diffusion_component();


  };
//---------------------------------------------------------------------------
#endif
