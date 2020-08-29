//---------------------------------------------------------------------------


#include "sputt.h"


// for debug purpose
// step size in mm for DEP1, DEP21, R, totDep...

const int grid_R2 = 2;

sputt::sputt(string log_fname)
{
    //constructor opens a file for loging intermediate
    // results
    logfile.open(log_fname, ios::out);
}

sputt::~sputt()
{
    logfile.close();
}

double sputt::atom_radius(double z)
{
  double r;
  if (z < 39)
    r = sqrt(0.77*0.77 + ((0.93*0.93 - 0.77*0.77)/18) * (z - 18));
  else
    r = sqrt(0.93*0.93 + ((1.04*1.04 - 0.93*0.93)/18) * (z - 36));

  return r;
};

double sputt::mean_free_path(double T)
{
  double lambda;
  double m = MA/MG;
  lambda = 1/(2.276516*(P/T)*(RG+RA)*(RG+RA)*sqrt(1+m));
    logfile << "Called mean_free_path() \n " ;
    logfile << "Calculated lambda: " << lambda << "\n";
  return lambda;
};

double sputt::avg_collision_number_to_thermalize(double T)
{ // calculates eta, avg energy of sputt atoms is 10 eV
  double m = MG/MA;
  double vper;
  double eta;
  vper = log(sqrt(1+m)+sqrt(m))/(4*pow(m,1.5)*sqrt(1+m));
  vper += (2*pow(m,4) + 5*pow(m,3) + 3*m*m - m -1)/(4*m*pow(1+m,3));
  vper = (1-m)/(1+m) + 2 * m * vper/(1+m);
  eta = log(sqrt(8.6174E-5*T/(E9*m)));
  eta /= log(vper);
    logfile << "Called avg_collision_number_to_thermalize() \n " ;
    logfile << "Calc eta: " << eta <<"\n";
    
  return eta;
};

void sputt::collisionless_transport()
{
  double dr=(R2-R1)/N1;
  double dfi= PI/N2;
  double R=0;
  double kk;
  int maxi=ceil(L5*R2/grid_R2);
  for (int i=0; i<=maxi; i++)
  {
    R = grid_R2*i;
    double S=0;
    for(int k1=0; k1<ceil(N1); k1++)
    {
      double RC = R1 + (k1 + 0.5) * dr;
      for (int k2=0; k2<ceil(N2); k2++)
      {
        double fi =(k2+0.5)*dfi;
        double R0 = Z*Z + R*R + RC*RC - 2*RC*R*cos(fi);
        S += RC * exp(-sqrt(R0)/LAMBDA)/(R0*R0);
      };
    };
    kk=DEP1[i]=Z*Z*dr*dfi*S*2/PI;
  };
    logfile << "Called collisionless_transport() \n";
    logfile << "Calculated DEP1 logged in: /Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_log_DEP1.txt\n";
    
    string log_fname = "/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_log_DEP1.txt";
    ofstream log_DEP1;
    log_DEP1.open(log_fname, ios::out);
    log_DEP1 << maxi << "\n";
    for (int i=0; i<=maxi; i++){
        log_DEP1 << i << "\t" << DEP1[i] << "\n";
    };
    log_DEP1.close();
  return;
};

double sputt::collisionless_dep_rate()
{
    logfile << "Called collisionless_dep_rate() \n";
  double RTOT = PI*pow(grid_R2/2.0,2)*DEP1[0];
  int maxj=ceil(L5*R2/grid_R2);
  for (int j=1; j<=maxj; j++)
  { double R = j * grid_R2;
    double A = PI * 2 * j * grid_R2 * grid_R2;
    RTOT += DEP1[j] * A;
  };
    logfile << "Calculated RTOT = Integrate(DEP1 dA): ";
    logfile << RTOT << "\n";

  double R99 = RTOT / (PI * (R2*R2 - R1*R1));
    logfile << "Calculated R99: " << R99 << "\n";

  return R99; // RDTOT / RSTOT
};

double sputt::sputtering_yield()
{
    logfile << "Called sputtering_yield() \n";
  double H, TH,ALPHA, PT, CP, EPS, SN, Y;
  if (MA/MG > 1)
    H = 0.834;
  else
    H = 0.18;

  TH = (1.5 * UA * pow((1+ 1.38 * pow(MG/MA,H)),2)) / (4*MA*MG/pow((MA+MG),2));
  ALPHA = 0.1 + 0.155*pow(MA/MG,0.73);
    //logfile << "CALC ALPHA:" << ALPHA <<"\n";
  PT = (1+MG/MA)*ZA*ZG*sqrt(pow(ZA,2.0/3)+pow(ZG,2.0/3)) /0.0325;
  CP = 3.56 * MG * ZG * ZA * ALPHA / ((MA+MG) * UA * sqrt(pow(ZA,2.0/3)+pow(ZG,2.0/3)));
  EPS = EI/PT;
  SN = 3.441 * sqrt(EPS) * log(EPS+2.718);
  SN /= (1 + 6.355*sqrt(EPS)+EPS*(6.882*sqrt(EPS)-1.708));
  Y = QA * CP *SN * pow(1-sqrt(TH/EI),2);
   logfile << "Calculated Y: " << Y << "\n";
  return Y;
};

double sputt::correct_T(double T)
{
  double TC, AA;
    logfile << "Called correct_T(): \n";
    
  LAMBDA = mean_free_path(T);
  ETA = avg_collision_number_to_thermalize(T);
  LAMBDA *= ETA;
    logfile << "Updated LAMBDA:" << LAMBDA << "\n";
  LALA = LAMBDA /(10*T);
    logfile << "Calculated LALA: " << LALA << "\n";
  collisionless_transport(); // writes DEP1[]
  R99 = collisionless_dep_rate();
  Y = sputtering_yield();
    
    
    
  if (!TCorr)   // if don't want T correction
      logfile << "No T correction performed\n";
    return T0;

  AA = CUR * E9 * Y * (1-R99) / (4*PI*KAPPA*LALA);
  TC = (T0 + sqrt(T0*T0 + 4*AA))/2;
    logfile << "Calc TC: " << TC <<"\n";
  if (fabs(T-TC)>0.1) // Need to correct T?
    return correct_T(TC);
  else
    return TC;
};

void sputt::difusion_source()
{ // want to calculate current J(r,z)
    logfile << "Called difusion_source() \n";
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
    logfile << "Calculated J logged in:  \n";
    string fname = "/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_log_J.txt";
    ofstream J_log;
    J_log.open(fname, ios::out);
    J_log << JN1 << '\t' << JN2 << '\n';
    for (int i=0; i < JN1; i++)
    {
        for (int j=0; j < JN2; j++){
            J_log << i << '\t' << j;
            J_log << '\t' << J[i][j];
            J_log << '\n';
        }
    }
    J_log.close();
  return;
};

void sputt::diffusion_component()
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
    logfile << "\n ******************* \n";
    logfile << "Called diffusion_component() \n";
    logfile << "Using dz: " << dz << " \n";
    logfile << "Logged [rop, A, B, C, kk1, kk2";
    logfile << "Logged [s0, s1, s2] \n in: ";
    string fname= "/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_DEP21_calc.txt \n";
    logfile << fname << "\n";
    double log_c[9];
    ofstream calc_log;
    calc_log.open(fname);
    calc_log << ceil(2*R2/grid_R2) +1;
    calc_log << "\t" << nz << "\t" << nr ;
    calc_log << "\t" << N2 << "\t" << "9" << "\n";
    
    logfile << "Logged DEP21 in: ";
    string fname2= "/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_DEP21.txt";
    logfile << fname2 << "\n";
    ofstream DEP21_log;
    DEP21_log.open(fname2);
    DEP21_log << ceil(2*R2/grid_R2) +1 << "\n";
    
    logfile << "Logged DEP22 in: ";
    string fname3= "/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_DEP22.txt";
    logfile << fname3 << "\n";
    ofstream DEP22_log;
    DEP22_log.open(fname3);
    DEP22_log << ceil(2*R2/grid_R2) +1 << "\n";

 // First, DEP1[], Z9=Z in Ivan's code
 for (int iir=0; iir<=ceil(2*R2/grid_R2); iir++)
 {
   double r = grid_R2*iir;
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
           log_c[0] = rop;
           log_c[1] = A;
           log_c[2] = B;
           log_c[3] = C;
           log_c[4] = kk1;
           log_c[5] = kk2;
           log_c[6] = s0;
           log_c[7] = s1;
           log_c[8] = s2;
           for (int rec =0; rec < 9; rec++){
               calc_log <<iir << "\t" << iz << "\t" << ir;
               calc_log << "\t" << k2 << "\t"<< rec ;
               calc_log << "\t" << log_c[rec] << "\n";
               };
           //};
       };
     };
     sr+=sz;
   };
   DEP21[iir]=2*sr;
     DEP21_log << iir << "\t" << DEP21[iir] << "\n";
 };
    calc_log.close();
    DEP21_log.close();

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
     DEP22_log << ir << "\t" << DEP22[ir] << "\n";
 };
    DEP22_log.close();
  return;
};

CData *sputt::calculate()
{ double *pR, *pTot;
    logfile << "Sputter log file to compare to python\n";
  RG = atom_radius(ZG);
    logfile << "RG calculated: " << RG << "\n";
  RA = atom_radius(ZA);
    logfile << "RA calculated: " << RA << "\n";
  TG = correct_T(T0);
    logfile << "Called T correction function\n";
    logfile << "TG calculated: " << TG << "\n";
  difusion_source();
  diffusion_component();

  int ir=ceil(2*R2/grid_R2)+1;
  pR = new double[ir];
  pTot = new double[ir];
  for (int i=0; i<ir; i++)
  {
    double c = (CUR * Y / (PI * ((R2*R2 - R1*R1)/100) * 1.6022E-19))* MA / DENA * (1.0E7 / 6.022E23);
    pR[i] = grid_R2*i;
    // renormalize dep rates form fraction of sputter current at target
    // to nm per sec

    DEP1[i] *= c;
    DEP21[i] *= c;
    pTot[i] = DEP1[i] + DEP21[i];
    DEP22[i] = -DEP22[i];
  };


  CData *panswer_data = new CData(ir, pR, pTot, DEP1, DEP21, DEP22);
    
    logfile.close();
    
  delete[] pR;
  delete[] pTot;
  return panswer_data;
};
//---------------------------------------------------------------------------







