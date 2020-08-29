//---------------------------------------------------------------------------


#pragma hdrstop

#include "CPlot.h"
#include <math>
//---------------------------------------------------------------------------

#pragma package(smart_init)

CPlot::CPlot()
{
  drawn =gasset = targetset = 0;
  return;
};

CPlot::~CPlot()
{
  if (drawn)
    delete pBitmap;

  if (gasset)
    delete Gas;

  if (targetset)
    delete Target;

  return;
};

CPlot::CPlot(CData datapar)
{
  data = datapar;
  minX = data.minX();
  minY = data.minY();
  maxX = data.maxX();
  maxY = data.maxY();

  XLabel[0] = minX;
  XLabel[1] = minX + (maxX - minX) /4;
  XLabel[2] = minX + 2 *(maxX - minX) /4;
  XLabel[3] = minX + 3 *(maxX - minX) /4;
  XLabel[4] = maxX;

  // redone to scale Y from 0, not minY
  minY = 0;
  
  YLabel[0] = minY;
  YLabel[1] = minY + (maxY - minY) /4;
  YLabel[2] = minY + 2 *(maxY - minY) /4;
  YLabel[3] = minY + 3 *(maxY - minY) /4;
  YLabel[4] = maxY;

  data_scaled=data;
  for (int i=0; i<data_scaled.numElements; i++)
  { // redone to scale Y from 0, not minY
    data_scaled.X[i] -= minX;
    data_scaled.X[i] /= (maxX - minX);
    //data_scaled.Y1[i] -= minY;
    data_scaled.Y1[i] /= (maxY /*- minY*/);
    //data_scaled.Y2[i] -= minY;
    data_scaled.Y2[i] /= (maxY /*- minY*/);
    //data_scaled.Y3[i] -= minY;
    data_scaled.Y3[i] /= (maxY /*- minY*/);
    //data_scaled.Y4[i] -= minY;
    data_scaled.Y4[i] /= (maxY /*- minY*/);
  };

  drawn = 0;

  return;
};

void CPlot::initData(CData datapar)
{
  if (drawn)
    delete pBitmap;

  data = datapar;
  // test
 // data.SaveASCII("C:\\windows\\escritorio\\initst.dat"); 

  minX = data.minX();
  minY = data.minY();
  maxX = data.maxX();
  maxY = data.maxY();

  XLabel[0] = minX;
  XLabel[1] = minX + (maxX - minX) /4;
  XLabel[2] = minX + 2 *(maxX - minX) /4;
  XLabel[3] = minX + 3 *(maxX - minX) /4;
  XLabel[4] = maxX;

  //redone to scale Y from 0 to maxY
  minY = 0;
  YLabel[0] = minY;
  YLabel[1] = minY + (maxY - minY) /4;
  YLabel[2] = minY + 2 *(maxY - minY) /4;
  YLabel[3] = minY + 3 *(maxY - minY) /4;
  YLabel[4] = maxY;

  data_scaled=data;
  for (int i=0; i<data_scaled.numElements; i++)
  { // redone to scale Y from 0, not minY
    data_scaled.X[i] -= minX;
    data_scaled.X[i] /= (maxX - minX);
    //data_scaled.Y1[i] -= minY;
    data_scaled.Y1[i] /= (maxY /*- minY*/);
    //data_scaled.Y2[i] -= minY;
    data_scaled.Y2[i] /= (maxY /*- minY*/);
    //data_scaled.Y3[i] -= minY;
    data_scaled.Y3[i] /= (maxY /*- minY*/);
    //data_scaled.Y4[i] -= minY;
    data_scaled.Y4[i] /= (maxY /*- minY*/);
  };


  drawn = 0;

  return;
};

Graphics::TBitmap *CPlot::GetBitmap()
{
  if(!drawn)
    Paint();
  return pBitmap;
};

bool CPlot::IsDrawn()
{
  return drawn;
};

void CPlot::paintFrame(sFrame in_frame)
{
        int x1, y1,  x2, y2;
        x1 = in_frame.x1; y1 = in_frame.y1;
        x2 = in_frame.x2; y2 = in_frame.y2;
        pBitmap->Canvas->Pen->Color = clBlack;
        pBitmap->Canvas->Pen->Width = 2;


        for (int i=1; i<20; i++)
        {  if (!(i==5 || i== 10 || i==15))
           {
           int pos = x1 + ceil(i*(x2-x1)/20.0);
           pBitmap->Canvas->MoveTo(pos,y1);
           pBitmap->Canvas->LineTo(pos, y1+10);
           pBitmap->Canvas->MoveTo(pos,y2);
           pBitmap->Canvas->LineTo(pos, y2-10);
           };
        };

        for (int i =1; i<20; i++)
        {  if (!(i==5 || i== 10 || i==15))
           {
           int pos = y1 + ceil(i*(y2-y1)/20.0);
           pBitmap->Canvas->MoveTo(x1,pos);
           pBitmap->Canvas->LineTo(x1+10,pos);
           pBitmap->Canvas->MoveTo(x2,pos);
           pBitmap->Canvas->LineTo(x2-10,pos);
           };
        };

        for (int i=1; i<4; i++)
        {
           int pos = x1+ ceil(i*(x2-x1)/4.0);
           pBitmap->Canvas->MoveTo(pos,y1);
           pBitmap->Canvas->LineTo(pos, y1+20);
           pBitmap->Canvas->MoveTo(pos,y2);
           pBitmap->Canvas->LineTo(pos, y2-20);
        };

        for (int i=1; i<4; i++)
        {
           int pos = y1 + ceil(i*(y2-y1)/4.0);
           pBitmap->Canvas->MoveTo(x1,pos);
           pBitmap->Canvas->LineTo(x1+20,pos);
           pBitmap->Canvas->MoveTo(x2,pos);
           pBitmap->Canvas->LineTo(x2-20,pos);
        };

        pBitmap->Canvas->MoveTo(x1, y1);
        pBitmap->Canvas->LineTo(x1, y2);
        pBitmap->Canvas->LineTo(x2, y2);
        pBitmap->Canvas->LineTo(x2, y1);
        pBitmap->Canvas->LineTo(x1, y1);
        return;
};

void CPlot::paintLabels(sFrame in_frame)
{
  AnsiString fname;
  int x1, y1,  x2, y2;
  x1 = in_frame.x1; y1 = in_frame.y1;
  x2 = in_frame.x2; y2 = in_frame.y2;

  pBitmap->Canvas->Font->Size = 12;
  fname = "Times New Roman";
  pBitmap->Canvas->Font->Name = fname;

  for (int i=0; i<5; i++)
  {
    int x,y;
    TSize sz;

    // X axis labels
    x = x1 + ceil(i*(x2-x1)/4.0);
    fname.printf("%.3G",XLabel[i]);
    sz = pBitmap->Canvas->TextExtent(fname);
    switch (i)
    {
      case 0:
        x=x1;
      break;
      case 4:
        x -= sz.cx;
      break;
      default:
        x -= ceil(sz.cx/2);
    };

    y = y2;
    y += ceil(0.25 * sz.cy);
    pBitmap->Canvas->TextOutA(x,y, fname);

    // Y axis labels
    y = y2 - ceil(i*(y2-y1)/4.0);
    fname.printf("%.2G",YLabel[i]);
    sz = pBitmap->Canvas->TextExtent(fname);
    x = x1;
    x -= ceil(0.5 * sz.cy);
    x -= sz.cx;
    y -= ceil(sz.cy/2);
    pBitmap->Canvas->TextOutA(x,y, fname);
  };

  return;
};

void CPlot::dotdraw(double x, double y, sFrame in_frame)
{
        int x1, y1,  x2, y2;
        x1 = in_frame.x1; y1 = in_frame.y1;
        x2 = in_frame.x2; y2 = in_frame.y2;
        int cx, cy;
        double dx, dy;

        dx = (x2-x1) * x + x1;
        dy = y2 - (y2-y1) * y;
        cx=ceil(dx);
        cy=ceil(dy);

        pBitmap->Canvas->Ellipse(cx-3, cy-3, cx+3, cy+3);
       // pBitmap->Canvas->MoveTo(cx,cy);
        return;
};

void CPlot::tracedraw(double px1, double py1, double px2, double py2, sFrame in_frame)
{
   int x1, y1,  x2, y2;
   x1 = in_frame.x1; y1 = in_frame.y1;
   x2 = in_frame.x2; y2 = in_frame.y2;
   int cx1, cx2, cy1, cy2;
   double dx1, dx2, dy1, dy2;


   dx1 = (x2-x1) * px1 + x1;
   dx2 = (x2-x1) * px2 + x1;
   dy1 = y2 - (y2-y1) * py1;
   dy2 = y2 - (y2-y1) * py2;
   cx1 = ceil(dx1);
   cx2 = ceil(dx2);
   cy1 = ceil(dy1);
   cy2 = ceil(dy2);
   pBitmap->Canvas->MoveTo(cx1, cy1);
   pBitmap->Canvas->LineTo(cx2, cy2);
   return;
};

void CPlot::Paint()
{
  AnsiString text, fname;
  TSize sz;
  sFrame graph_frame;
  Graphics::TBitmap *pBmp1, *pBmp2;
  int x1, x2, y1, y2;

  // code to paint graphic...

  if (drawn)
    delete pBitmap;
  else
    drawn =1;

  pBitmap = new Graphics::TBitmap;
  pBitmap->Width = 800;
  pBitmap->Height = 600;

  int header_botom = paintHeader();

  y1 = graph_frame.y1 = header_botom;
  x1 = graph_frame.x1 = 50;
  x2 = graph_frame.x2 = 790;
  y2 = graph_frame.y2 = 570;

  graph_frame.x1 += 100;
  graph_frame.y1 += 30;
  graph_frame.y2 -= 30;

  paintFrame(graph_frame);
  paintLabels(graph_frame);

  pBitmap->Canvas->Pen->Width = 1;

  pBitmap->Canvas->Pen->Color = clRed;
  pBitmap->Canvas->Brush->Color = clRed;
  pBitmap->Canvas->Brush->Style = bsSolid;
  for (int i=0; i<data.numElements; i++)
    dotdraw(data_scaled.X[i], data_scaled.Y1[i], graph_frame);
  for (int i=0; i<data.numElements-1; i++)
    tracedraw(data_scaled.X[i], data_scaled.Y1[i], data_scaled.X[i+1], data_scaled.Y1[i+1], graph_frame);

  pBitmap->Canvas->Pen->Color = clGreen;
  pBitmap->Canvas->Brush->Color = clGreen;
  for (int i=0; i<data.numElements; i++)
    dotdraw(data_scaled.X[i], data_scaled.Y2[i], graph_frame);
  for (int i=0; i<data.numElements-1; i++)
    tracedraw(data_scaled.X[i], data_scaled.Y2[i], data_scaled.X[i+1], data_scaled.Y2[i+1], graph_frame);

  pBitmap->Canvas->Pen->Color = clBlue;
  pBitmap->Canvas->Brush->Color = clBlue;
  for (int i=0; i<data.numElements; i++)
    dotdraw(data_scaled.X[i], data_scaled.Y3[i], graph_frame);
  for (int i=0; i<data.numElements-1; i++)
    tracedraw(data_scaled.X[i], data_scaled.Y3[i], data_scaled.X[i+1], data_scaled.Y3[i+1], graph_frame);

  pBitmap->Canvas->Pen->Color = clMaroon;
  pBitmap->Canvas->Brush->Color = clMaroon;
  for (int i=0; i<data.numElements; i++)
    dotdraw(data_scaled.X[i], data_scaled.Y4[i], graph_frame);
  for (int i=0; i<data.numElements-1; i++)
    tracedraw(data_scaled.X[i], data_scaled.Y4[i], data_scaled.X[i+1], data_scaled.Y4[i+1], graph_frame);


 // Print axis labels
 // Y label printed vertically

  text="Deposition rate (nm/s)";
  fname = "Times New Roman";
  pBmp1 = new Graphics::TBitmap;
  pBmp1->Canvas->Font->Color = clBlack;
  pBmp1->Canvas->Font->Size = 24;
  pBmp1->Canvas->Font->Name = fname;
  sz = pBmp1->Canvas->TextExtent(text);
  pBmp1->Width = sz.cx;
  pBmp1->Height = sz.cy;
  pBmp1->Canvas->TextOutA(0,0 , text);
  pBmp2 = new Graphics::TBitmap;
  pBmp2->Width = sz.cy;
  pBmp2->Height = sz.cx;
  for (int x=0; x<sz.cx; x++)
  {  for (int y=0; y<sz.cy; y++)
        pBmp2->Canvas->Pixels[y][sz.cx-x] = pBmp1->Canvas->Pixels[x][y];
  };
  delete pBmp1;
  pBitmap->Canvas->Draw(20,y1 + ceil(0.5*(y2-y1-sz.cx)), pBmp2);
  delete pBmp2;

  text="R (mm)";
  pBitmap->Canvas->Font->Style = TFontStyles();
  pBitmap->Canvas->Font->Size = 24;
  pBitmap->Canvas->Font->Color = clBlack;
  pBitmap->Canvas->Brush->Style = bsClear;
  sz = pBitmap->Canvas->TextExtent(text);
  pBitmap->Canvas->Font->Color = clBlack;
  pBitmap->Canvas->TextOutA((x1 + ceil(0.5*(100 +x2 - x1 - sz.cx))), 600 - sz.cy, text);


  return;
};

void CPlot::SaveASCII(char *fname)
{
 ofstream ofs;
 ofs.open(fname);

 ofs << "Magnetron sputtering of " << Target;
 ofs << " by " << Gas << "\n";
 ofs << "Gas pressure: " << Preasure << " Pa\n";
 ofs << "Initial T: " << T0 << " K\t\tCorrected T: " << Tf << " K\n";
 ofs << "Discharge current: " << I << " A\t\tDischarge voltage " << V << " V\n";
 ofs << "Erosion disk inner radius: " << in_radius << " mm\t\tOuter radius: " << out_radius << " mm\n";
 ofs << "Substrate to target distance: " << distance << " mm\n";
 ofs << "Average energy of sputtered atoms: " << Eavg << " eV\n\n";

 for(int i=0; i<data.numElements; i++)
 {
   ofs << data.X[i] << "\t" << data.Y1[i] << "\t";
   ofs << data.Y2[i] << "\t" << data.Y3[i] << "\t";
   ofs << data.Y4[i] << "\n";
 };

 ofs.close();
 return;
};

int CPlot::paintHeader()
{
  AnsiString font_name = "Times New Roman";
  AnsiString text;
  int left_mrg=20;
  int top_mrg=20;
  int top_pos = top_mrg;
  int interline;
  int lineHeight;
  int legend_pos;
  int ret_val;

  pBitmap->Canvas->Font->Name = font_name;
  pBitmap->Canvas->Font->Size=12;
  pBitmap->Canvas->Font->Color = clBlack;

  text.printf("Magnetron sputtering of %s by %s", Target, Gas);
  interline = ceil(0.5 * pBitmap->Canvas->TextHeight(text));
  lineHeight = pBitmap->Canvas->TextHeight(text);
  pBitmap->Canvas->TextOutA(left_mrg, top_mrg, text);

  top_pos += interline + lineHeight;
  text.printf("Gas pressure: %g Pa       Initial T: %g K       Corrected T: %g K", Preasure, T0, Tf);
  pBitmap->Canvas->TextOutA(left_mrg, top_pos, text);

  top_pos += interline + lineHeight;
  text.printf("Discharge current: %g A       Discharge voltage: %g V", I, V);
  pBitmap->Canvas->TextOutA(left_mrg, top_pos, text);

  top_pos += interline + lineHeight;
  text.printf("Erosion disk inner radius: %g mm       Outer radius: %g mm       Distance to target: %g mm", in_radius, out_radius, distance);
  pBitmap->Canvas->TextOutA(left_mrg, top_pos, text);

  top_pos += interline + lineHeight;
  text.printf("Sputtered atoms average energy: %g eV", Eavg);
  pBitmap->Canvas->TextOutA(left_mrg, top_pos, text);

  top_pos += interline + lineHeight;
  pBitmap->Canvas->MoveTo(10,10);
  pBitmap->Canvas->Pen->Width = 1;
  pBitmap->Canvas->Pen->Style = psSolid;
  pBitmap->Canvas->Pen->Color = clBlack;
  pBitmap->Canvas->LineTo(790, 10);
  pBitmap->Canvas->LineTo(790, top_pos);
  pBitmap->Canvas->LineTo(10, top_pos);
  pBitmap->Canvas->LineTo(10,10);

  ret_val = top_pos;

  top_pos = top_mrg + 4*(interline + lineHeight);
  legend_pos = 780;
  text.printf("Resputtered at target");
  pBitmap->Canvas->Font->Style = TFontStyles() << fsBold;
  legend_pos -= pBitmap->Canvas->TextWidth(text);
  pBitmap->Canvas->Font->Color = clMaroon;
  pBitmap->Canvas->TextOutA(legend_pos, top_pos, text);

  pBitmap->Canvas->Font->Style = TFontStyles() << fsUnderline;
  text.printf("Legend");
  pBitmap->Canvas->Font->Color = clBlack;
  pBitmap->Canvas->TextOutA(legend_pos, top_mrg, text);
  pBitmap->Canvas->Font->Style = TFontStyles() << fsBold;

  top_pos = top_mrg +  (interline + lineHeight);
  text.printf("Total deposition rate");
  pBitmap->Canvas->Font->Color = clRed;
  pBitmap->Canvas->TextOutA(legend_pos, top_pos, text);


  top_pos = top_mrg + 2*(interline + lineHeight);
  text.printf("Ballistic component");
  pBitmap->Canvas->Font->Color = clGreen;
  pBitmap->Canvas->TextOutA(legend_pos, top_pos, text);

  top_pos = top_mrg + 3*(interline + lineHeight);
  text.printf("Diffusive component");
  pBitmap->Canvas->Font->Color = clBlue;
  pBitmap->Canvas->TextOutA(legend_pos, top_pos, text);
  
  return ret_val;
};

void CPlot::setGas(const char *pgas)
{
  int l = strlen(pgas) +1;
  Gas = new char[l];
  strcpy(Gas, pgas);
  gasset=1;
  return;
};

void CPlot::setTarget(const char *ptarget)
{
  int l = strlen(ptarget) +1;
  Target = new char[l];
  strcpy(Target, ptarget);
  targetset=1;
  return;
};




