//---------------------------------------------------------------------------

#ifndef CDataH
#define CDataH

#include <fstream>


using namespace std;

class CData
{
  public:
  int numElements;
  double *X;
  double *Y1;
  double *Y2;
  double *Y3;
  double *Y4;

  CData();
  ~CData();
  CData(int, double*, double*, double*, double*, double*);
  CData & operator =(const CData & operand);
  int AddData(double dX, double dY1, double dY2, double dY3, double dY4);
  void SaveASCII(char *ofile);
  void ReadASCII(char *ifile);
  double minX();
  double maxX();
  double minY();
  double maxY();

  private:
  void clear();
};

//---------------------------------------------------------------------------
#endif


