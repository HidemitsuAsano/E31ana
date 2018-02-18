#ifndef CDCHit_h
#define CDCHit_h 1

#include <string>
#include <iostream>
#include <vector>

#include "ConfMan.h"
#include "GlobalVariables.h"
#include "ChamberLikeHit.h"

class CDCHit : public ChamberLikeHit
{
 private:
  
 public:
  CDCHit();
  CDCHit( const CDCHit &hit );

 private:
  int SuperLayer, ASDNum, TransType, ASDChannel;

  double Radius, Phi; // <--> (WirePosX,WirePosY,WirePosZ)
  double WirePosXp, WirePosYp, WirePosZp; // wire position opposite side of readout

 public:
  int slayer()    const { return SuperLayer; }
  int asdnum()    const { return ASDNum; }
  int transtype() const { return TransType; }
  int asdch()     const { return ASDChannel; }
  double radius() const { return Radius; }
  double phi()    const { return Phi; }
  double wxp()    const { return WirePosXp; }
  double wyp()    const { return WirePosYp; }
  double wzp()    const { return WirePosZp; }

  void SetSuperLayer( const int &lay ) { SuperLayer = lay; }
  void SetASDNum(     const int &num ) { ASDNum = num; }
  void SetTransType(  const int &typ ) { TransType = typ; }
  void SetASDChannel( const int &ch  ) { ASDChannel = ch; }
  void SetRadius(     const double &rad ) { Radius = rad; }
  void SetPhi(        const double &phi ) { Phi = phi; }
  void SetWireOppositePosition( const double &x, const double &y, const double &z )
    { WirePosXp = x; WirePosYp = y; WirePosZp = z; }

  void Calc( ConfMan *conf, const double &toffset=0.);

  void Clear();

  ClassDef( CDCHit, 1);
};

#endif
