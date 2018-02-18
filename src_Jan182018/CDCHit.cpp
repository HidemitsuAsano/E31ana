#include "CDCHit.h"
#include <cmath>

ClassImp(CDCHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCHit::CDCHit()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
CDCHit::CDCHit( const CDCHit &hit )
{
  *this = hit; // I cannot understand why I need this copy constructor.
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDCHit::Calc( ConfMan *conf, const double &toffset )
{
  double Tilt;
  int tmp1,tmp2;
  if( conf->GetCDCWireMapManager() ){
//     conf->GetCDCWireMapManager()->GetWire( cr(), sl(), ch(), Radius, Phi, Tilt );
    conf->GetCDCWireMapManager()->GetWire( cr(), sl(), ch(),
					   SuperLayer, ASDNum, TransType, ASDChannel,
					   Radius, Phi, Tilt,
					   tmp1,tmp2 );

    SetTiltAngle(Tilt);
    SetWireLength( conf->GetCDCWireMapManager()->zlen() );
    SetWirePosX( Radius*cos(Phi*TMath::DegToRad()) );
    SetWirePosY( Radius*sin(Phi*TMath::DegToRad()) );
    SetWirePosZ(0.);

    if( 0<fabs(Tilt) ){
      double S = length()*tan(Tilt*TMath::DegToRad());
      double rad2 = Radius*Radius;
      double cost = 1.-S*S/(2.*rad2);
      double sint = S/(2.*rad2)*sqrt(4.*rad2-S*S);
      WirePosXp = wx()*cost - wy()*sint; WirePosYp = wx()*sint + wy()*cost;
    }
    else{
      WirePosXp = wx(); WirePosYp = wy();
    }
    SetWirePosZ( length()/2. ); WirePosZp = -length()/2.;

    SetHitPosition( wx(), wy(), 0 ); // just intial value

    double gx,gy,gz,dgx,dgy,dgz;
    conf->GetCDCWireMapManager()->GetGPOS( gx,gy,gz,dgx,dgy,dgz );
    SetGPos( gx, gy, gz );
    SetGDir( dgx, dgy, dgz );
    SetRotationAngle( conf->GetCDCWireMapManager()->rot() );
  }

  if( conf->GetGainMapManager() ){
    double dtime;
    if( 0<tdc() ){
      dtime = conf->GetGainMapManager()->CalcCValue( cr(), sl(), ch(), 1, tdc() );
      SetDriftTime(dtime);
    }
  }

  SetTimeOffset(conf,toffset);

  //  if( conf->GetXTMapManager() ){
//     double time = dt() - toffs();
//     //if( 0<TDC0 ) DriftLength = conf->GetXTMapManager()->CalcDriftLength( CounterID, Layer, Wire, DriftTime );
//     double dlen;
//     if( 0<tdc() ){
//       dlen = conf->GetXTMapManager()->CalcDriftLength( cis(), layer(), wire(), time );
//       SetDriftLength(dlen);
//     }
//  }

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CDCHit::Clear()
{
  SuperLayer = ASDNum = TransType = ASDChannel = -1;
  Radius = Phi = -999;
  WirePosXp = WirePosYp = WirePosZp = -999;
}
