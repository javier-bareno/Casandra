//---------------------------------------------------------------------------

#ifndef CPlotH
#define CPlotH

#include "CData.h"
#include <vcl.h>
#include <fstream>
#include <string.h>

struct sFrame
{
        int x1;
        int x2;
        int y1;
        int y2;
};

class CPlot
{
  private:
  CData data, data_scaled;
  double minX, minY, maxX, maxY;
  double XLabel[5], YLabel[5];
  Graphics::TBitmap *pBitmap;
  bool drawn;
  bool gasset, targetset;

  public:
  // Data for labels from calculation
  char* Gas; // EG "Ar"
  char* Target; // EG "Ti"
  double Preasure;
  double T0;
  double Tf;
  double I;
  double V;
  double in_radius;
  double out_radius;
  double distance;
  double Eavg;
  void setGas(const char *);
  void setTarget(const char *);

  public:
  CPlot();
  ~CPlot();
  CPlot(CData data);
  void initData(CData datapar);
  Graphics::TBitmap *GetBitmap();
  bool IsDrawn();
  void Paint();
  void SaveASCII(char *fname);

  private:
  //routines used by Paint
  void paintFrame(sFrame);
  void paintLabels(sFrame);
  int paintHeader();
  void dotdraw(double x, double y, sFrame in_frame);
  void tracedraw(double px1, double py1, double px2, double py2, sFrame in_frame);


};
//---------------------------------------------------------------------------
#endif


