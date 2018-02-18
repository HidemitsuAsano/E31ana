#ifndef SDDHit_h
#define SDDHit_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include <TF1.h>

#include "ConfMan.h"
#include "vme_common.h"

class SDDHit : public TObject
{
 private:
  int CounterID;
  int HitID; 
  
 public:
  SDDHit();
  SDDHit( const SDDHit &hit );
  SDDHit( const double &rndm );
  virtual ~SDDHit() {};

 private:
  // raw data
  int Seg;
  int Out,Fout,Xout;
  int TDCsdd,TDCt0;
  int Thigh,Treset;

  int CrateOut, SlotOut, ChannelOut;
  int CrateXout, SlotXout, ChannelXout;
  int CrateFout, SlotFout, ChannelFout;
  int CrateTDCsdd, SlotTDCsdd, ChannelTDCsdd;
  int CrateTDCt0, SlotTDCt0, ChannelTDCt0;
  int CrateThigh, SlotThigh, ChannelThigh;
  int CrateTreset,SlotTreset,ChannelTreset;

  double RANDOM;

 public:
  void SetRandom(const double &v){ RANDOM=v; } 
  double rndm() const { return RANDOM; }
  int cid()  const { return CounterID; }
  int hid()  const { return HitID; }
  int seg()  const { return Seg; }
  int out()  const { return Out; }
  int xout() const { return Xout; }
  int fout() const { return Fout; }
  int thigh() const { return Thigh; }
  int treset() const { return Treset; }
  int tdcsdd() const { return TDCsdd; }
  int tdct0() const { return TDCt0; }

  int crout() const { return CrateOut; }
  int slout() const { return SlotOut; }
  int chout() const { return ChannelOut; }
  int crxout() const { return CrateXout; }
  int slxout() const { return SlotXout; }
  int chxout() const { return ChannelXout; }
  int crfout() const { return CrateFout; }
  int slfout() const { return SlotFout; }
  int chfout() const { return ChannelFout; }

  int crtdcsdd() const { return CrateTDCsdd; }
  int sltdcsdd() const { return SlotTDCsdd; }
  int chtdcsdd() const { return ChannelTDCsdd; }
  int crtdct0()  const { return CrateTDCt0; }
  int sltdct0()  const { return SlotTDCt0; }
  int chtdct0()  const { return ChannelTDCt0; }

  int crthigh()  const { return CrateThigh; }
  int slthigh()  const { return SlotThigh; }
  int chthigh()  const { return ChannelThigh; }
  int crtreset() const { return CrateTreset; }
  int sltreset() const { return SlotTreset; }
  int chtreset() const { return ChannelTreset; }

  void SetCounterID( const int & cid ) { CounterID = cid; }
  void SetHitID(     const int & hid ) { HitID = hid; }
  void SetSegment(   const int & seg ) { Seg = seg; }
  void SetOut(       const int & v   ) { Out = v; }
  void SetXout(      const int & v   ) { Xout = v; }
  void SetFout(      const int & v   ) { Fout = v; }
  void SetTDCsdd(    const int & v   ) { TDCsdd = v; }
  void SetTDCt0(     const int & v   ) { TDCt0 = v; }
  void SetThigh(     const int & v   ) { Thigh = v; }
  void SetTreset(    const int & v   ) { Treset = v; }

  void SetCrateOut(   const int & cr  ) { CrateOut = cr; }
  void SetSlotOut(    const int & sl  ) { SlotOut = sl; }
  void SetChannelOut( const int & ch  ) { ChannelOut = ch; }
  void SetCrateXout(   const int & cr  ) { CrateXout = cr; }
  void SetSlotXout(    const int & sl  ) { SlotXout = sl; }
  void SetChannelXout( const int & ch  ) { ChannelXout = ch; }
  void SetCrateFout(   const int & cr  ) { CrateFout = cr; }
  void SetSlotFout(    const int & sl  ) { SlotFout = sl; }
  void SetChannelFout( const int & ch  ) { ChannelFout = ch; }

  void SetCrateTDCsdd(   const int & cr  ) { CrateTDCsdd = cr; }
  void SetSlotTDCsdd(    const int & sl  ) { SlotTDCsdd = sl; }
  void SetChannelTDCsdd( const int & ch  ) { ChannelTDCsdd = ch; }
  void SetCrateTDCt0(   const int & cr  ) { CrateTDCt0 = cr; }
  void SetSlotTDCt0(    const int & sl  ) { SlotTDCt0 = sl; }
  void SetChannelTDCt0( const int & ch  ) { ChannelTDCt0 = ch; }

  void SetCrateThigh(   const int & cr  ) { CrateThigh = cr; }
  void SetSlotThigh(    const int & sl  ) { SlotThigh = sl; }
  void SetChannelThigh( const int & ch  ) { ChannelThigh = ch; }
  void SetCrateTreset(   const int & cr  ) { CrateTreset = cr; }
  void SetSlotTreset(    const int & sl  ) { SlotTreset= sl; }
  void SetChannelTreset( const int & ch  ) { ChannelTreset= ch; }

  
  bool SetData( const int &c, const int &n, const int &a,
		const int &cid, const int &seg, const int &at, const int &ud,
		const int &data );
  bool SetData( const int &seg,	EventStruct &vme);
  
  bool CheckRange();
  
 private:
  // conversion param
  double Gain,Offset;

  // conversion data
  double Ene,Timesdd,Timet0;
  // corrected time
  double CTimesdd,CTimet0;
  double CFout;

  double  X,  Y,  Z; // position
  double dX, dY, dZ; // direction
  double  GX,  GY,  GZ; // origin position
  double dGX, dGY, dGZ; // origin direction
  double Length, Width, Thick;

 public:
  double gain() const { return Gain; }
  double offset() const { return Offset; }
  double ene()  const { return Ene; }
  double cene() const { return Gain*(Out-Offset); }
  double cene(const double &g,const double &o) const { return g*(Out-o); }
  double timesdd() const { return Timesdd; }
  double timet0()  const { return Timet0; }
  double cfout() const { return CFout; }

  double x()   const { return X; }
  double y()   const { return Y; }
  double z()   const { return Z; }
  double dx()  const { return dX; }
  double dy()  const { return dY; }
  double dz()  const { return dZ; }
  double gx()  const { return GX; }
  double gy()  const { return GY; }
  double gz()  const { return GZ; }
  double dgx() const { return dGX; }
  double dgy() const { return dGY; }
  double dgz() const { return dGZ; }
  double len() const { return Length; }
  double wid() const { return Width; }
  double thick() const { return Thick; }
  void   pos( double &x, double &y, double &z ) { x=X;  y=Y;  z=Z; }
  void   dir( double &x, double &y, double &z ) { x=dX; y=dY; z=dZ; }
  void   gpos( double &x, double &y, double &z ) { x=GX;  y=GY;  z=GZ; }
  void   gdir( double &x, double &y, double &z ) { x=dGX; y=dGY; z=dGZ; }

  void SetEne(     const double &e ) { Ene = e; }
  void SetGain(     const double &e ) { Gain = e; }
  void SetOffset(     const double &e ) { Offset = e; }
  void SetTimesdd( const double &t ) { Timesdd = t; }
  void SetTimet0(  const double &t ) { Timet0 = t; }
  void SetCFout( const double &f) { CFout = f; }

  void SetParEne( const double &g, const double &o) { Gain = g; Offset = o; }
  void GetParEne( double &g, double &o) { g = Gain; o = Offset; }


 private:
  int PADCl, PADCh;
  int Nsample;
  double Interval;
  int FADC[MaxFADCPoint];

  double VGain,VOffset;
  double VEne;

  double FBaseHeight,FBaseSlope;
  double FPeakHeight,FPeakTime;
  double FPostSlope;

  int FBaseNpar;
  double FBasePar[5];
  int FPeakNpar;
  double FPeakPar[5];
  int FPostNpar;
  double FPostPar[5];

  double FGain,FOffset;
  double FEne;
  
 public:

  int padcl() const { return PADCl; }
  int padch() const { return PADCh; }
  int nsample() const { return Nsample; }
  double interval() const { return Interval; }
  int fadc(const int &nth) const { return FADC[nth]; }
  
  double fbase() const { return FBaseHeight; }
  double fbaseslope() const { return FBaseSlope; }
  double fpeakh() const { return FPeakHeight; }
  double fpeakt() const { return FPeakTime; }
  double fpostslope() const { return FPostSlope; }

  int fpeaknpar() const { return FPeakNpar; }
  double fpeakpar(const int &i) const { return FPeakPar[i]; }
  double *fpeakpar()  { return FPeakPar; }
  int fbasenpar() const { return FBaseNpar; }
  double fbasepar(const int &i) const { return FBasePar[i]; }
  double *fbasepar()  { return FBasePar; }
  int fpostnpar() const { return FPostNpar; }
  double fpostpar(const int &i) const { return FPostPar[i]; }
  double *fpostpar()  { return FPostPar; }

  double vgain() const { return VGain; }
  double voffset() const { return VOffset; }
  double vene()  const { return VEne; }
  double cvene() const { return VGain*(PADCh-VOffset); }
  double cvene(const double &g,const double &o) const { return g*(PADCh-o); }

  double fgain() const { return FGain; }
  double foffset() const { return FOffset; }
  double fene()  const { return FEne; }
  double cfene() const { return FGain*(FPeakHeight-FOffset); }
  double cfene(const double &g,const double &o) const { return g*(FPeakHeight-o); }
  
  void SetParVEne( const double &g, const double &o) { VGain = g; VOffset = o; }
  void GetParVEne( double &g, double &o) { g = VGain; o = VOffset; }

  void SetParFEne( const double &g, const double &o) { FGain = g; FOffset = o; }
  void GetParFEne( double &g, double &o) { g = FGain; o = FOffset; }
  
  bool Calc( ConfMan *conf );
  void Clear();

 public:
  //  void Reverse( ConfMan *conf );

  ClassDef(SDDHit, 1 );
};

#endif
