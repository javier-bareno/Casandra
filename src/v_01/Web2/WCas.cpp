//---------------------------------------------------------------------------


#include <fstream>
#include <sputt.h>
#include <string>


#pragma hdrstop

//---------------------------------------------------------------------------

#pragma argsused
int main(int argc, char* argv[])
{
  ifstream ifile;
  ofstream test;
  ofstream ofdat, ofbmp;
  sputt system;
  CPlot dplot;
  CData *ansData;
  Graphics::TBitmap *pbmp;
  string ofname;
  AnsiString ofnamebmp, ofnamedat;

  // Gas atom params
  string Gas_symbol;
  double MG, ZG, RG, UG, QG, DENG, KAPA;
  // Target element
  string Target_symbol;
  double MA, ZA, RA, UO, QZ, DEN;
  // other params
  double disch_volt, Eavg,  disch_cur, pres, T0;
  int TC;
  bool Tcorrect;
  double in_rad, out_rad, N1, N2, L5, dist;
  // output files name (.bmp and .dat)

  // open input file and load params
  test.open(argv[2]);
  test << "This is it";
  test.close();

  ifile.open(argv[1]);
  ifile >> Gas_symbol;
  ifile >> MG >> ZG >> RG >> UG >> QG >> DENG >> KAPA;
  ifile >> Target_symbol;
  ifile >> MA >> ZA >> RA >> UO >> QZ >> DEN;
  ifile >> disch_volt >> disch_cur >> pres >> T0;
  ifile >> Eavg >> in_rad >> out_rad >> N1 >> N2 >> L5 >> dist;
  ifile >> TC;
  Tcorrect = (TC == 1);
  ifile >> ofname;

  ifile.close();

  dplot.setGas(Gas_symbol.c_str());
  dplot.setTarget(Target_symbol.c_str());
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
  dplot.V = disch_volt;
  dplot.Eavg = Eavg;
  system.E9 = Eavg;
  system.CUR = disch_cur;
  dplot.I = disch_cur;
  system.P = pres;
  dplot.Preasure = pres;
  system.T0 = dplot.T0 = T0;
  system.R1 = dplot.in_radius = in_rad;
  system.R2 = dplot.out_radius = out_rad;
  system.N1 = N1;
  system.N2 = N2;
  system.L5 = L5;
  system.Z = dplot.distance = dist;
  system.TCorr = Tcorrect;
  system.QA = QZ;

  ansData = system.calculate();
  dplot.Tf = system.TG;
  dplot.initData(*ansData);
  pbmp = dplot.GetBitmap();

  ofnamedat = ofname.c_str();
  ofnamedat += ".dat";
  ofnamebmp = ofname.c_str();
  ofnamebmp += ".bmp";

  pbmp->SaveToFile(ofnamebmp);
  char *name;
  name = ofnamedat.c_str();
  dplot.SaveASCII(name);

  return 0;
}
//---------------------------------------------------------------------------
