#include "CDCHit.h"
#include <cmath>

ClassImp(CDCHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCHit::CDCHit():
  SuperLayer(-1),ASDNum(-1),TransType(-1),ASDChannel(-1),
  Radius(-999.),Phi(-999.),
  WirePosp(DEFVECT)
{
  //  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CDCHit::Calc( ConfMan *conf, const double &toffset )
{
  double Tilt;
  int tmp1,tmp2;
  CDCWireMapMan* wireman=conf->GetCDCWireMapManager();
  if(wireman){
    //     wireman->GetWire( cr(), sl(), ch(), Radius, Phi, Tilt );
    wireman->GetWire( cr(), sl(), ch(),
		      SuperLayer, ASDNum, TransType, ASDChannel,
		      Radius, Phi, Tilt,
		      tmp1,tmp2 );
    
    SetTiltAngle(Tilt);
    SetWireLength( wireman->zlen() );
    TVector3 tmppos =wireman->GetWirePos(  cr(), sl(), ch() );
    TVector3 tmpposp=wireman->GetWirePosp( cr(), sl(), ch() );
    SetWirePosition(tmppos.X(),tmppos.Y(),tmppos.Z());
    SetWireOppositePosition(tmpposp.X(),tmpposp.Y(),tmpposp.Z());

    SetHitPosition( wx(), wy(), 0 ); // just intial value
    
    double gx,gy,gz,dgx,dgy,dgz;
    wireman->GetGPOS( gx,gy,gz,dgx,dgy,dgz );
    SetGPos( gx, gy, gz );
    SetGDir( dgx, dgy, dgz );
    SetRotationAngle( wireman->rot() );
  }

  tdchit.Calc(conf);
  SetTimeOffset(conf,toffset);
    return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDCHit::Clear()
{
  ChamberLikeHit::Clear();
  SuperLayer = ASDNum = TransType = ASDChannel = -1;
  Radius = Phi = -999;
  WirePosp.SetXYZ(-999,-999,-999);
}
