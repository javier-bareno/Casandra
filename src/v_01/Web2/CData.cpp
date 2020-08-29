//---------------------------------------------------------------------------


#pragma hdrstop

#include "CData.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

CData::CData()
{
  // do nothing, really
  numElements =0;
  X = Y1 = Y2 = Y3 = Y4 = 0;
  return;
};

CData::~CData()
{
  clear();
  return;
};

CData::CData(int pnum, double* pX, double* pY1, double* pY2, double* pY3, double* pY4)
{
  numElements = pnum;
  X = new double[pnum];
  Y1 = new double[pnum];
  Y2 = new double[pnum];
  Y3 = new double[pnum];
  Y4 = new double[pnum];

  for (int i=0; i<pnum; i++)
  {
    X[i] = pX[i];
    Y1[i] = pY1[i];
    Y2[i] = pY2[i];
    Y3[i] = pY3[i];
    Y4[i] = pY4[i];
  };
};


CData & CData::operator=(const CData &operand)
{
  numElements = operand.numElements;
  X = new double[numElements];
  Y1  = new double[numElements];
  Y2  = new double[numElements];
  Y3 = new double[numElements];
  Y4 = new double[numElements];

  for (int i=0; i<numElements; i++)
  {
   X[i] = operand.X[i];
   Y1[i] = operand.Y1[i];
   Y2[i] = operand.Y2[i];
   Y3[i] = operand.Y3[i];
   Y4[i] = operand.Y4[i];
  };

  return *this;
};

int CData::AddData(double dX, double dY1, double dY2, double dY3, double dY4)
{
  double *pX, *pY1, *pY2, *pY3, *pY4;
  numElements ++;
  pX = new double[numElements];
  pY1  = new double[numElements];
  pY2  = new double[numElements];
  pY3 = new double[numElements];
  pY4 = new double[numElements];

  for (int i=0; i<numElements-1; i++)
  {
    pX[i] = X[i];
    pY1[i] = Y1[i];
    pY2[i] = Y2[i];
    pY3[i] = Y3[i];
    pY4[i] = Y4[i];
  };
  pX[numElements-1] = dX;
  pY1[numElements-1] = dY1;
  pY2[numElements-1] = dY2;
  pY3[numElements-1] = dY3;
  pY4[numElements-1] = dY4;

  if (numElements >1)
  {
    delete[] X;
    delete[] Y1;
    delete[] Y2;
    delete[] Y3;
    delete[] Y4;
  };
  
  X=pX;
  Y1=pY1;
  Y2=pY2;
  Y3=pY3;
  Y4=pY4;

  return 0;
};

void CData::SaveASCII(char *ofile)
{
 ofstream ofs;
 ofs.open(ofile);

 for(int i=0; i<numElements; i++)
 {
   ofs << X[i] << "\t" << Y1[i] << "\t";
   ofs << Y2[i] << "\t" << Y3[i] << "\t";
   ofs << Y4[i] << "\n";
 };

 ofs.close();
 return;
};

void CData::ReadASCII(char* ifile)
{
  double dX, dY1, dY2, dY3, dY4;

  ifstream ifs;
  ifs.open(ifile);
  ifs >> dX >> dY1 >> dY2 >> dY3 >> dY4;
  while(!ifs.eof())
  {
    AddData(dX, dY1, dY2, dY3, dY4);
    ifs >> dX >> dY1 >> dY2 >> dY3 >> dY4;
  };

  ifs.close();
  return;
};

double CData::minX()
{
  double min;
  min = X[0];
  for (int i=1; i<numElements; i++)
    if (X[i] < min)
      min = X[i];

  return min;
};

double CData::minY()
{
  double min;
  min = Y1[0];
  for (int i=0; i<numElements; i++)
  {  if (Y1[i] < min)
      min = Y1[i];
     if (Y2[i] < min)
      min = Y2[i];
     if (Y3[i] < min)
      min = Y3[i];
     if (Y4[i] < min)
      min = Y4[i];
  };

  return min;
};

double CData::maxX()
{
  double max;
  max = X[0];
  for (int i=1; i<numElements; i++)
    if (X[i] > max)
      max = X[i];

  return max;
};

double CData::maxY()
{
  double max;
  max = Y1[0];
  for (int i=0; i<numElements; i++)
  {  if (Y1[i] > max)
      max = Y1[i];
     if (Y2[i] > max)
      max = Y2[i];
     if (Y3[i] > max)
      max = Y3[i];
     if (Y4[i] > max)
      max = Y4[i];
  };

  return max;
};

void CData::clear()
{
  if (numElements!=0)
  {
    delete[] X;
    delete[] Y1;
    delete[] Y2;
    delete[] Y3;
    delete[] Y4;
    numElements =0;
    X = Y1 = Y2 = Y3 = Y4 = 0;
  };
  return;
};

