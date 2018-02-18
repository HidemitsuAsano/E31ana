// CDCWireMapMan.h

#ifndef CDCWireMapMan_h
#define CDCWireMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "GlobalVariables.h"
#include "TObject.h"

class CDCWireMap : public TObject
{
 private:
  int Layer, Wire;
  int SuperLayer, ASDNum, TransType, ASDch;
  double Radius, Phi, Tilt;
 public:
  CDCWireMap();
  CDCWireMap( int layer, int wire );
  CDCWireMap( int layer, int wire, 
	     int slayer, int asdnum, int ttype, int asdfch,
	     double rad, double phi, double tilt );
  ~CDCWireMap() {}

  CDCWireMap( const CDCWireMap &right );
  CDCWireMap & operator = ( const CDCWireMap &right );

 public:
  void SetLayer( const int &layer ) { Layer = layer; }
  void SetWire( const int &wire ) { Wire = wire; }
  void SetSuperLayer( const int &slayer ) { SuperLayer = slayer; }
  void SetASDNum( const int &asd ) { ASDNum = asd; }
  void SetTransType( const int &trans ) { TransType = trans; }
  void SetASDChannel( const int &ch ) { ASDch = ch; }
  void SetRadius( const double &rad ) { Radius = rad; }
  void SetPhi( const double &phi ) { Phi = phi; }
  void SetTilt( const double &tilt ) { Tilt = tilt; }

  int layer()       const { return Layer; }
  int wire()        const { return Wire; }
  int slayer()      const { return SuperLayer; }
  int asdnum()      const { return ASDNum; }
  int transtype()   const { return TransType; }
  int asdch()       const { return ASDch; }
  double radius()   const { return Radius; }
  double phi()      const { return Phi; }
  double tilt()     const { return Tilt; }
  
  ClassDef( CDCWireMap, 1 );
};

class RCDCWireMap : public TObject
{
 private:
  int Crate, Slot, Channel;
  int SuperLayer, ASDNum, TransType, ASDch;
  double Radius, Phi, Tilt;
 public:
  RCDCWireMap();
  explicit RCDCWireMap( int cr, int sl, int ch );
  explicit RCDCWireMap( int cr, int sl, int ch,
			int slayer, int asdnum, int ttype, int asdfch,
			double rad, double phi, double tilt );
  ~RCDCWireMap() {}
  
  RCDCWireMap( const RCDCWireMap &right );
  RCDCWireMap & operator = ( const RCDCWireMap &right );

 public:

  void SetCrate( const int &cr ) { Crate = cr; }
  void SetSlot( const int &sl ) { Slot = sl; }
  void SetChannel( const int &ch ) { Channel = ch; }
  void SetSuperLayer( const int &slayer ) { SuperLayer = slayer; }
  void SetASDNum( const int &asd ) { ASDNum = asd; }
  void SetTransType( const int &trans ) { TransType = trans; }
  void SetASDChannel( const int &ch ) { ASDch = ch; }
  void SetRadius( const double &rad ) { Radius = rad; }
  void SetPhi( const double &phi ) { Phi = phi; }
  void SetTilt( const double &tilt ) { Tilt = tilt; }

  int cr()          const { return Crate; }
  int sl()          const { return Slot; }
  int ch()          const { return Channel; }
  int slayer()      const { return SuperLayer; }
  int asdnum()      const { return ASDNum; }
  int transtype()   const { return TransType; }
  int asdch()       const { return ASDch; }
  double radius()   const { return Radius; }
  double phi()      const { return Phi; }
  double tilt()     const { return Tilt; }

  ClassDef( RCDCWireMap, 1 );
};

class CDCWireMapMan : public TObject
{
 public:
  CDCWireMapMan();
  explicit CDCWireMapMan( const std::string &filename1,
			 const std::string &filename2,
			 const std::string &filename3);
  virtual ~CDCWireMapMan();
 private:
  std::string WireMapFileName;
  std::string ASDMapFileName;
  std::string ChannelMapFileName;
  std::string GeometryFileName;

 public:
  void SetFileName( const std::string &filename1,
		    const std::string &filename2,
		    const std::string &filename3);
  bool Initialize();

  void Clear();

 private:
  typedef std::map <unsigned int,  CDCWireMap> WireMapContainer;
  typedef std::map <unsigned int, RCDCWireMap> RWireMapContainer;
  WireMapContainer wContainer;
  RWireMapContainer rwContainer;

  //double Xoffset, Yoffset, Zoffset, Thetaoffset, Phioffset;
  double GX,GY,GZ,dGX,dGY,dGZ;
  double Radius[NumOfCDCLayers+1], Phi0[NumOfCDCLayers+1];
  double dPhi[NumOfCDCLayers+1], Tilt[NumOfCDCLayers+1];

  double RotationAngle;
  double ZLengthOfWire;
  double InnerRadius, OuterRadius;

 public:
  void SetWireMap( const std::map <unsigned int, CDCWireMap> &container )
    { wContainer = container; }
  void SetWireMap( const unsigned int i, const CDCWireMap &map )
    { wContainer[i] = map; }
  void SetRWireMap( const std::map <unsigned int, RCDCWireMap> &container )
    { rwContainer = container; }
  void SetRWireMap( const unsigned int i, const RCDCWireMap &map )
    { rwContainer[i] = map; }

  std::string GetFileName1() { return ASDMapFileName; }
  std::string GetFileName2() { return ChannelMapFileName; }
  std::string GetFileName3() { return GeometryFileName; }

  bool GetWire( const int &cr, const int &sl, const int &ch,
		int &slayer, int &asdnum, int &ttype, int &asdch,
		double &rad, double &phi, double &tilt,
		int &layer, int &wire );

  bool GetWire( const int &cr, const int &sl, const int &ch,
		int &layer, int &wire );

  bool GetWire( const int &cr, const int &sl, const int &ch,
		double &rad, double &phi, double &tilt );

  bool GetWire( const int &layer, const int &wire,
		int &slayer, int &asdnum, int &ttype, int &asdch,
		double &rad, double &phi, double &tilt,
		int &cr, int &sl, int &ch );

  bool GetWire( const int &layer, const int &wire,
		double &rad, double &phi, double &tilt );

  bool GetGeom( const int &layer, double &radius, 
		double &phi0, double &dphi, double &tilt );

  bool GetGPOS( double &gx, double &gy, double &gz, 
		//double &thetaoff, double &phioff );
		double &dgx, double &dgy, double &dgz );
  bool GetFrame( double &zlen, double &rin, double &rout );
  
/*   double xoffset() const { return Xoffset; } */
/*   double yoffset() const { return Yoffset; } */
/*   double zoffset() const { return Zoffset; } */
  double gx() const { return GX; }
  double gy() const { return GY; }
  double gz() const { return GZ; }
  double dgx() const { return dGX; }
  double dgy() const { return dGY; }
  double dgz() const { return dGZ; }
/*   double thetaoffset() const { return Thetaoffset; } */
/*   double phioffset() const { return Phioffset; } */
  double radius( const int &layer ) const { return Radius[layer]; }
  double phi0( const int &layer ) const { return Phi0[layer]; }
  double dphi( const int &layer ) const { return dPhi[layer]; }
  double tilt( const int &layer ) const { return Tilt[layer]; }
  int    nw( const int &layer ) const { return (int)(360/dPhi[layer]+0.2); }

  double zlen() const { return ZLengthOfWire; }
  double rin()  const { return InnerRadius; }
  double rout() const { return OuterRadius; }
  double rot()  const { return RotationAngle; }

  bool PrintSimpleWireMap( std::ostream &p_out = std::cout );
  bool PrintWireMap();

  ClassDef( CDCWireMapMan, 1 );
};

#endif
