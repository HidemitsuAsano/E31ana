// UshiwakaFieldMapMan.h
#ifndef UshiwakaFieldMapMan_h
#define UshiwakaFieldMapMan_h 1

#include <map>
#include <string>
#include <stdlib.h>
#include <fstream>

#include <TObject.h>
#include <TVector3.h>
#include <TVirtualPad.h>
#include <TH2.h>
#include <TArrow.h>
#include <TGraph.h>

#include "GlobalVariables.h"
#include "USWKTrack.h"

static const std::string DefaultFieldMapName = "param/Run47/inoue/UshiwakaFieldMap.param";

class UshiwakaFieldMapMan : public TObject
{
 public:
  UshiwakaFieldMapMan();
  UshiwakaFieldMapMan(const std::string &filename);
  ~UshiwakaFieldMapMan();

 private:
  std::string FileName;
  //*** for runge-kutta ***//
  double DT;
  int SR;
  std::vector<USWKTrack> Tracks;

  //*** for field map ***//
  double GX, GY, GZ;
  double X_MAX, Y_MAX, Z_MAX;
  double DX, DY, DZ;
  int NX, NY, NZ;

  typedef std::map<unsigned long, TVector3> FieldMapContainer;
  FieldMapContainer fieldMap;

 public:
  void SetdT(const double &dt){ DT=dt; };
  double dt() const { return DT; };
  double samplingrate() const { return SR; };

  void GetGPOS(double &x, double &y, double &z){ x=GX, y=GY, z=GZ; };
  void GetRange(double &x, double &y, double &z){ x=X_MAX, y=Y_MAX, z=Z_MAX; };
  void GetDL(double &x, double &y, double &z){ x=DX, y=DY, z=DZ; };
  void GetNData(int &nx, int &ny, int &nz){ nx=NX, ny=NY, nz=NZ; };

 public:
  void SetFileName(const std::string filename){ FileName = filename; };
  std::string GetFileName(){ return FileName; };
  bool Initialize();
  TVector3 value(const double x, const double y, const double z);
  TVector3 gvalue(const double gx, const double gy, const double gz);
  //*** for Runge-Kutta ***//
  int NTrack() const { return Tracks.size(); };
  USWKTrack* GetUSWKTrack(const int &i){ return &Tracks[i]; };

  void Clear();
  bool RungeKutta(std::string &name, TVector3 &inpos, TVector3 &inmom);
  bool RungeKutta(const char* name, TVector3 &inpos, TVector3 &inmom);
  bool InField(const TVector3 &pos);

  //*** for check ***//
  void DumpStatus();

  void DrawXGraph(TVirtualPad *pad, const int direction, const int iy, const int iz);
  void DrawYGraph(TVirtualPad *pad, const int direction, const int ix, const int iz);
  void DrawZGraph(TVirtualPad *pad, const int direction, const int ix, const int iy);

  bool DrawZGraph(TVirtualPad *pad, const int direction, const double ix, const double iy, const int n=100);

  ClassDef(UshiwakaFieldMapMan,1);
};
#endif
