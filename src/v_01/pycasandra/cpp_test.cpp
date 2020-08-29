#include<iostream>
#include<fstream>
#include<cmath>
#ifndef PI
  #define PI 3.1415927
#endif
using namespace std;

//public:
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
  double grid_R2 =2.0;
  

  //private:
  // parameters from calculations and results
  double LAMBDA; // mean free path
  double ETA; // ~ number of collisions prior thermalization
  double LALA; // lambda/(10*TG)
  double DEP1[100], DEP21[100], DEP22[100], DEP[100];
  double Y; // sputtering yield;
  double R99; // collisionless rate?
  double J[100][100]; //diff current
  int JN1, JN2, DEPN;

void difusion_source()
{ // want to calculate current J(r,z)
  int nr = ceil(L5*R2/grid_R2);
  double drv = L5*R2/nr;
  int nz = ceil(Z/2); // -> dz = 2mm
    //ensure dz < lambda
  double dz = Z/nz;
  if (dz > LAMBDA)
  {  nz = Z/LAMBDA;
     dz = Z/nz;
  };

  for (int ir=0; ir < nr; ir++)
  {
    double r = (0.5+ir)*drv;
    for(int iz=0; iz < nz; iz++)
    {
      double z = (iz+0.5) *dz;
      double s=0;
      for (int k1=0; k1 < N1; k1++)
      {
        double rc = R1 + (k1+0.5) *(R2-R1)/N1;
        for (int k2=0; k2 < N2; k2++)
        {
          double fi = (k2 +0.5) * PI/N2;
          double r0 = sqrt(z*z + r*r + rc*rc  - 2*r*rc*cos(fi));
          s += rc* z * exp(-r0/LAMBDA) /(LAMBDA * pow(r0, 3));
        };
      };
    J[ir][iz] = ((R2-R1)/N1) * (PI/N2) * s * 2 /PI;
    };
  };
  JN1 = nr;
  JN2 = nz;
  return;
};

void diffusion_component()
{
  int nr = ceil(L5*R2/grid_R2);
  double drv = L5*R2/nr;
  int nz = ceil(Z/2); // -> dz = 2mm
    //ensure dz < lambda
  double dz = Z/nz;
  if (dz > LAMBDA)
  {  nz = Z/LAMBDA;
     dz = Z/nz;
  };

 // First, DEP1[], Z9=Z in Ivan's code
 for (int ir=0; ir<=ceil(2*R2/grid_R2); ir++)
 {
   double r = grid_R2*ir;
   double sr=0;
   for (int iz=0; iz<nz; iz++)
   {
     double zc = (iz+0.5)*dz;
     double sz=0;
     for (int ir=0; ir<nr; ir++)
     {
       double rc = (ir+0.5)*drv;
       for (int k2=0; k2<N2; k2++)
       {
         double dfi = PI/N2;
         double fi = (k2+0.5)*dfi;
         double rop = sqrt(r*r +rc*rc - 2*r*rc*cos(fi));
         double A= (Z-zc)/(2*Z);
         double B= (Z+zc)/(2*Z);
         double C= rop /(2*Z);
         double kk1 = A/pow(A*A+C*C, 1.5);
         double kk2 = B/pow(B*B+C*C, 1.5);
         double s0= kk1 - kk2;
         double s1= (1+A) /pow((1+A)*(1+A)+C*C, 1.5) + (1-B) /pow((1-B)*(1-B)+C*C, 1.5);
         s1 -= (1-A) /pow((1-A)*(1-A)+C*C, 1.5) + (1+B) /pow((1+B)*(1+B)+C*C, 1.5);
         s0 += s1;
         double s2= J[ir][iz]*drv*dz*dfi*rc*s0/(16*Z*Z*PI);
         sz+=s2;
       };
     };
     sr+=sz;
   };
   DEP21[ir]=2*sr;
 };

 // now DEP22[] Z9=0 in Ivan's code
 for (int ir=0; ir<=ceil(2*R2/grid_R2); ir++)
 {
   double r = grid_R2*ir;
   double sr=0;
   for (int iz=0; iz<nz; iz++)
   {
     double zc = (iz+0.5)*dz;
     double sz=0;
     for (int ir=0; ir<nr; ir++)
     {
       double rc = (ir+0.5)*drv;
       for (int k2=0; k2<N2; k2++)
       {
         double dfi = PI/N2;
         double fi = (k2+0.5)*dfi;
         double rop = sqrt(r*r +rc*rc - 2*r*rc*cos(fi));
         double A= (0.0-zc)/(2*Z);
         double B= (0.0+zc)/(2*Z);
         double C= rop /(2*Z);
         double s0= A/pow(A*A+C*C, 1.5) - B/pow(B*B+C*C, 1.5);
         double s1= (1+A) /pow((1+A)*(1+A)+C*C, 1.5) + (1-B) /pow((1-B)*(1-B)+C*C, 1.5);
         s1 -= (1-A) /pow((1-A)*(1-A)+C*C, 1.5) + (1+B) /pow((1+B)*(1+B)+C*C, 1.5);
         s0 += s1;
         double s2= J[ir][iz]*drv*dz*dfi*rc*s0/(16*Z*Z*PI);
         sz+=s2;
       };
     };
     sr+=sz;
   };
   DEP22[ir]=2*sr;
   DEPN = ir +1;
 };
  return;
};

int main()
{
    ofstream ofile;

    L5 =2; R1 =10; R2 = 15; grid_R2 = 2; Z=100;
    LAMBDA = 5; N1 = 10; N2 = 20;
    cout <<"Getting started";
    difusion_source();
    cout << "Called dif_source\n";
    diffusion_component();
    cout << "Called dif_component\n";
    
    ofile.open("out.txt");
    ofile << "J " << JN1 << " x " << JN2 << "\n";
    for (int i=0; i<JN1; i++)
    {
      for(int j=0; j<JN2; j++)
      {
        ofile << J[i][j] << "\t";
      };
      ofile << "\n";
    };

    ofile << "\nDEP21 " << DEPN << "\n";
    for (int i=0; i< DEPN; i++)
    {
      ofile << DEP21[i] << "\t";
    };
    ofile << "\n";

    ofile.close();
}