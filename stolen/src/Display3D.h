#ifndef DISPLAY3D_h
#define DISPLAY3D_h 1

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include <TSystem.h>
#include <TAtt3D.h>
#include <TVirtualViewer3D.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TGLViewer.h>
#include <TVector3.h>
#include <TAttLine.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>

#include "GlobalVariables.h"
#include "ConfMan.h"
#include "CDSHitMan.h"
#include "ChamberLikeHit.h"
#include "HodoscopeLikeHit.h"
#include "BeamLineHitMan.h"
#include "CDSTrackingMan.h"

class Shape : public TObject
{
 public:
  Shape( Int_t color, Double_t x, Double_t y, Double_t z );
  ~Shape() {};
  virtual TBuffer3D & GetBuffer3D( UInt_t reqSections ) = 0;

  void SetName( const std::string & name ) { Name = name; }
  std::string name() { return Name; }
  void RotationX( const Double_t &ang );
  void RotationY( const Double_t &ang );
  void RotationZ( const Double_t &ang );
  void Scale( const Double_t &sx, const Double_t &sy, const Double_t &sz );
  void Scale( const Double_t &sca );
  void SetTranslation( Double_t *mat );
  void GLTranslationCalc( Double_t *mat );
  void GLTranslationCalc( Double_t *mat1, Double_t *mat2, Double_t *mat3 );
  void SetColor( const Int_t &color ) { fColor = color; }
  void SetTransparency( const Int_t &trans ) { Transparency = trans; }
 protected:
  Double_t fX, fY, fZ;
  Int_t fColor;

  std::string Name;
  Int_t Transparency;
  Double_t TransMat[16];

  ClassDef(Shape,0);
};

class Sphere : public Shape
{
 public:
  Sphere( Int_t  color, Double_t x, Double_t y, Double_t z, Double_t r );
  ~Sphere() {};

  virtual TBuffer3D & GetBuffer3D( UInt_t reqSections );
 private:
  Double_t fRadius;

  ClassDef(Sphere,0);
};

class Box : public Shape
{
 public:
  Box( Int_t color, Double_t x, Double_t y, Double_t z,
       Double_t dx, Double_t dy, Double_t dz );
  ~Box() {};

  virtual TBuffer3D & GetBuffer3D( UInt_t reqSections );
 private:
  Double_t fDX, fDY, fDZ;

  ClassDef(Box,0);
};

class Tube : public Shape
{
 public:
  Tube( Int_t  color, Double_t x, Double_t y, Double_t z, 
	Double_t len, Double_t irad, Double_t orad, 
	Double_t phi_start, Double_t phi_end );
  ~Tube() {};

  virtual TBuffer3D & GetBuffer3D( UInt_t reqSections );
 private:
  Double_t fLength;
  Double_t fInnerRadius, fOuterRadius;
  Double_t fPhiStart, fPhiEnd;

  ClassDef(Tube,0);
};

class Line : public TObject
{
 public:
  Line();
  ~Line() {};

  void SetName( const std::string &name ) { Name = name; }
  std::string name() { return Name; }
  void SetPoints( std::string name, Int_t npoints, Double_t *x, Double_t *y, Double_t *z,
		  Int_t color=1, Int_t width=1 );
  Int_t np() { return pCon.size(); }
  TVector3 vec( Int_t i ) { return pCon[i]; }
  Int_t color() { return Color; }
  Int_t width() { return Width; }
 private:
  std::string Name;
  std::vector <TVector3> pCon;
  Int_t Color, Width;

  ClassDef(Line,0);
};

class Mark : public TObject
{
 public:
  Mark();
  ~Mark() {};

  void SetName( const std::string &name ) { Name = name; }
  std::string name() { return Name; }
  void SetPoints( std::string name, Int_t npoints, Double_t *x, Double_t *y, Double_t *z,
		  Int_t color=1, Int_t size=1);
  Int_t np() { return mCon.size(); }
  TVector3 vec( Int_t i ) { return mCon[i]; }
  Int_t color() { return Color; }
  Int_t size() { return Size; }
 private:
  std::string Name;
  std::vector <TVector3> mCon;
  Int_t Color, Size;

  ClassDef(Mark,0);
};

class Display3D : public TObject, public TAtt3D
{
 public:
  Display3D();
  Display3D( TVirtualPad *p );
  ~Display3D();

 private:
  //Display3D( Display3D & );

 public:
  void Draw(Option_t *option);
  void Paint(Option_t *option);

  // CDS
  void SetCDSSolenoidBody( ConfMan *conf, Double_t phimin=0, Double_t phimax=360, Int_t transparency = 50 );
  void SetCDSSolenoidCap( ConfMan *conf, Int_t pattern = 0xf, Int_t transparency = 50 );
  void ReSetCDH( ConfMan *conf, Int_t transparency = 50 );
  void SetCDH( ConfMan *conf, Int_t transparency = 50 );
  void SetCDHColor( ConfMan *conf, Int_t seg, Int_t color, Int_t transparency = 50 );
  void SetCDCFrame( ConfMan *conf,  Int_t transparency = 50 );
  void SetTarget( ConfMan *conf, Int_t transparency = 50 );
  void SetCDCHitWire( ConfMan *conf, CDSHitMan *cds );
  void SetCDCHit( ConfMan *conf, CDSHitMan *cds );
  void SetCDHHit( ConfMan *conf, CDSHitMan *cds );
  void SetCDCTrack( ConfMan *conf, CDSTrackingMan *trackMan );
  void SetCDCTrackHit( ConfMan *conf, CDSTrackingMan *trackMan );
  void SetCDCVertexPoint( ConfMan *conf, CDSTrackingMan *trackMan );

  // special
  void SetTargetOct2010( ConfMan *conf, Int_t transparency = 50 );

  void DeleteCDSSolenoidBody();
  void DeleteCDSSolenoidCap();
  void DeleteCDH( Int_t seg );
  void DeleteCDH();
  void DeleteTarget();

  // BeamLine
  void SetCounterFrame( ConfMan *conf, Int_t CID, Bool_t fLocal, Int_t transparency = 50 );
  void SetCounterHit( ConfMan *conf, BeamLineHitMan *bl, Int_t CID, Bool_t fLocal, Int_t transparency = 50 );
  
 public:

  void SetBox( const std::string &name, const Int_t &color, 
	       const Double_t &x, const Double_t &y, const Double_t &z,
	       const Double_t &dx, const Double_t &dy, const Double_t &dz,
	       const Int_t &transparency = 50 );

  void SetSphere( const std::string &name, const Int_t &color,
		  const Double_t &x, const Double_t &y, const Double_t &z,
		  const Double_t &r, const Int_t &transparency = 50 );

  void SetTube( const std::string &name, const Int_t &color,
		const Double_t &x, const Double_t &y, const Double_t &z,
		const Double_t &len, const Double_t &irad, const Double_t &orad,
		const Double_t &phi_start, const Double_t &phi_end, const Int_t &transparency = 50 );
  
  void SetPLine( const std::string &name, const Int_t &npoints, Double_t *x, Double_t *y, Double_t *z,
		 const Int_t &color, const Int_t &width );

  void SetPMark( const std::string &name, const Int_t &npoints, Double_t *x, Double_t *y, Double_t *z,
		 const Int_t &color, const Int_t &size );

  void RotationX( const std::string &name, const double &theta ); // unit:degree
  void RotationY( const std::string &name, const double &theta ); // unit:degree
  void RotationZ( const std::string &name, const double &theta ); // unit:degree
  void Scale(     const std::string &name, 
		  const double &sx, const double &sy, const double &sz );
  void Scale(     const std::string &name, const double &sca );

  bool Delete( const std::string &name );
  bool DeleteAll();

  void PrintObj();

 private:
  std::map < std::string, Mark * > PMarkCon;
  std::map < std::string, Line * > PLineCon;
  std::map < std::string, Shape * > ShapesCon;

  TVirtualPad *pad;

  ClassDef( Display3D, 1 );
};


#endif
