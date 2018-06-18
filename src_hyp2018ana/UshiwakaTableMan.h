#ifndef UshiwakaTableMan_h
#define UshiwakaTableMan_h 1

#include "GlobalVariables.h"
#include <fstream>
#include <map>

static const std::string DefaultTableName = "param/Run47/inoue/UshiwakaTable.param";

class UshiwakaTable : public TObject
{
 public:
  UshiwakaTable();
  UshiwakaTable(const double &x, const double &y, const double &dx, const double &dy, const double &mom, const double &ang, const double &fl); 
  virtual ~UshiwakaTable(){};

 private:
  double X, Y, dX, dY, Mom, Ang, FL;

 public:
  void SetX(const double &x){ X=x; };
  void SetY(const double &y){ Y=y; };
  void SetdX(const double &dx){ dX=dx; };
  void SetdY(const double &dy){ dY=dy; };
  void SetMom(const double &mom){ Mom=mom; };
  void SetAng(const double &ang){ Ang=ang; };
  void SetFL(const double &fl){ FL=fl; };

  double x() const { return X; };
  double y() const { return Y; };
  double dx() const { return dX; };
  double dy() const { return dY; };
  double mom() const { return Mom; };
  double ang() const { return Ang; };
  double fl() const { return FL; };

  ClassDef(UshiwakaTable, 1);
};

class UshiwakaTableMan : public TObject
{
  UshiwakaTableMan(const UshiwakaTableMan &right){};
  void operator = (const UshiwakaTableMan &right){};
 public:
  UshiwakaTableMan();
  UshiwakaTableMan(const std::string &filename);
  virtual ~UshiwakaTableMan(){};

 private:
  std::string FileName;

  double DT;
  double Zmax, Zmin;

  int NX;
  double Xmax, Xmin, dX;

  int NDX;
  double DXmax, DXmin, dDX;

  int NY;
  double Ymax, Ymin, dY;

  int NDY;
  double DYmax, DYmin, dDY;

  int NMom;
  double MomMin, MomMax, dMom;

 public:
  void SetFileName(const std::string &filename){ FileName = filename; };
  std::string GetFileName() const { return FileName; };
  typedef std::map<unsigned long, UshiwakaTable> TableMap;
  TableMap tableMap;

 public:
  UshiwakaTable* table(int ix, int iy, int idx, int idy, int iang);
  bool GetParam(const double &x, const double &y, const double &dx, const double &dy, const double &ang, double &mom, double &fl);

  double dt() const { return DT; };
  double zmax() const { return Zmax; };
  double zmin() const { return Zmin; };

  int nx() const { return NX; };
  double xmax() const { return Xmax; };
  double xmin() const { return Xmin; };
  double dx() const { return dX; };

  int ny() const { return NY; };
  double ymax() const { return Ymax; };
  double ymin() const { return Ymin; };
  double dy() const { return dY; };

  int ndx() const { return NDX; };
  double dxmax() const { return DXmax; };
  double dxmin() const { return DXmin; };
  double ddx() const { return dDX; };

  int ndy() const { return NDY; };
  double dymax() const { return DYmax; };
  double dymin() const { return DYmin; };
  double ddy() const { return dDY; };

  int nmom() const { return NMom; };
  double mommax() const { return MomMax; };
  double mommin() const { return MomMin; };
  double dmom() const { return dMom; };

 public:
  bool Initialize();
  bool DumpStatus();

  ClassDef(UshiwakaTableMan, 1);
};
#endif
