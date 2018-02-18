#include <iomanip>
#include <new>

#include "ChamberLikeHit.h"
#include "TVector2.h"

ClassImp(ChamberLikeHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ChamberLikeHit::ChamberLikeHit()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ChamberLikeHit::ChamberLikeHit( const ChamberLikeHit &right )
{
  //  *this = right; // I cannot understand why I need this copy constructor.
  CounterID=right.CounterID; HitID=right.HitID; XY=right.XY;
  Layer=right.Layer; Wire=right.Wire; TDC0=right.TDC0;
  Crate=right.Crate; Slot=right.Slot; Channel=right.Channel;
  DriftTime=right.DriftTime; DriftLength=right.DriftLength; TimeOffset=right.TimeOffset; CorrDriftLength=right.CorrDriftLength;
  WirePosX=right.WirePosX; WirePosY=right.WirePosY; WirePosZ=right.WirePosZ;
  WireLength=right.WireLength; TiltAngle=right.TiltAngle; Resolution=right.Resolution;
  HitPosX=right.HitPosX; HitPosY=right.HitPosY; HitPosZ=right.HitPosZ; LR=right.LR;
  RotationAngle=right.RotationAngle;
  GX=right.GX; GY=right.GY; GZ=right.GZ;
  dGX=right.dGX; dGY=right.dGY; dGZ=right.dGZ;
  dXY=right.dXY;
  ResolutionTrue=right.ResolutionTrue; TimeOffsetTrue=right.TimeOffsetTrue;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ChamberLikeHit::Clear()
{
  CounterID = HitID = XY = -1;
  Layer = Wire = TDC0 = -1;
  //  XY = LR = -1;
  Crate = Slot = Channel = -1;
  
  LR=-1;

  DriftTime = DriftLength = CorrDriftLength = -999;

  TimeOffset = 0.;

  WirePosX = WirePosY = WirePosZ = -999.; // at readout
  WireLength = TiltAngle = Resolution = -999;
  HitPosX = HitPosY = HitPosZ = -999;
  
  GX = GY = GZ = dGX = dGY = dGZ = -999;
  dXY= -999;
  ResolutionTrue = 0.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ChamberLikeHit::SetData( const int &c, const int &n, const int &a,
			      const int &cid, const int &layer, const int &wire, const int &data )
{
#if 0
  std::cout << c << " " << n << " " << a << " " << cid << " " << layer << "  " << wire << " " << data 
	    << std::endl;
#endif
  CounterID = cid;
  Layer = layer;
  Wire  = wire;
  TDC0   = data;
  Crate = c; Slot = n; Channel = a;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ChamberLikeHit::SetSimulatedResolution( ConfMan *conf )
{
  if( conf->GetReslMapManager() ){
    //ResolutionTrue = conf->GetReslMapManager()->GetResolution( CounterID, Layer, Wire );
    double dummy;
    conf->GetReslMapManager()->GetResolution( CounterID, Layer, Wire, ResolutionTrue, dummy );
    DriftLength += ResolutionTrue;
    //    std::cout<<"CDC resl :"<< ResolutionTrue<< std::endl;
    Reverse(conf);
    return true;
  }
  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ChamberLikeHit::Reverse( ConfMan *conf )
{
  if( conf->GetXTMapManager() ){
    DriftTime = conf->GetXTMapManager()->CalcDriftTime( CounterID, Layer, Wire, DriftLength );
    DriftTime += TimeOffsetTrue;
#if 0
    std::cout << " DriftLength:" << DriftLength << " ResolutionTrue:" << ResolutionTrue
	      << " DriftTime:" << DriftTime << std::endl;
#endif
  }
  if( conf->GetGainMapManager() ){
    TDC0 = conf->GetGainMapManager()->CalcDATValue( Crate, Slot, Channel, 1, DriftTime );
#if 0
    std::cout << " TDC:" << TDC0 << std::endl;
#endif
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ChamberLikeHit::Calc( ConfMan *conf, const double &toffs )
{
#if 0
  std::cout << "cid: " << CounterID << "  ly: " << Layer << "  wire: " << Wire 
	    << std::endl;
#endif

  if( conf->GetBLDCWireMapManager() ){
    int nw,xy;
    double z, xy0, dxy;
    conf->GetBLDCWireMapManager()->GetParam( CounterID, Layer, nw, z, xy, xy0, dxy,
					     WireLength, TiltAngle, RotationAngle );
    if( nw<Wire ){
      std::cout << " too many wires !!! " << std::endl;
      return false;
    }

    WirePosZ = z;  //cm
    //    XY = xy;
    dXY = TMath::Abs(dxy);

    if( xy==0 ){
      XY = 0;
      WirePosY = 0;
      WirePosX = xy0 + (Wire-1)*dxy;
    }
    else{
      XY = 1;
      WirePosY = xy0 + (Wire-1)*dxy;
      WirePosX = 0;
    }
#if 0
    std::cout << " wx:" << WirePosX << "  xy: " << WirePosY << "  wz: " << WirePosZ << std::endl;
#endif
    SetHitPosition( WirePosX, WirePosY, WirePosZ ); // just intial value
    conf->GetBLDCWireMapManager()->GetGParam( CounterID, GX, GY, GZ, dGX, dGY, dGZ );
#if 0
    std::cout << " GX:" << GX << "  GY: " << GY << "  GZ: " << GZ << std::endl;
#endif
  }

  if( conf->GetGainMapManager() ){
    if( 0<TDC0 ) DriftTime = conf->GetGainMapManager()->CalcCValue( Crate, Slot, Channel, 1, TDC0 );
  }

  SetTimeOffset(conf,toffs);
  //  double tmptime = 0;
  //  SetTimeOffset(conf,tmptime);

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ChamberLikeHit::SetTimeOffset( ConfMan *conf, const double &toffs )
{
  if( conf->GetXTMapManager() ){
    TimeOffset = toffs;
    //    std::cout<< "DriftTime: "<<DriftTime << "\tT0Time: "<<toffs;
    DriftTime -= TimeOffset;
    //    std::cout<< "DriftTimeNew: "<<DriftTime << std::endl;
    if( 0<TDC0 ) DriftLength = conf->GetXTMapManager()->CalcDriftLength( CounterID, Layer, Wire, DriftTime );
  }
  if( DriftLength < 0. ) DriftLength=0.;

  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ChamberLikeHit::CheckRange(const double &low, const double &high){
  //  std::cout<<"ChamberLikeHit::CheckRange()"<<std::endl;
  //  std::cout<<DriftLength<<"\t"<<dXY<<std::endl;
  if(DriftLength > low*dXY/2. && DriftLength < high *dXY/2. ) return true;
  else return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ChamberLikeHit::gwpos(double &tempx,double &theta,double &tempz,const bool &TILT,const bool &GPOS)
{
  if(XY) tempx=wy(); else tempx=wx();
  tempz=wz();
  theta=rot()+90.*xy();
  if(TILT){ // local tilt
    double zcenter;
    if(cid()==CID_BLC1a) zcenter=-16.;
    if(cid()==CID_BLC1b) zcenter=16.;
    if(cid()==CID_BLC2a) zcenter=-13.75;
    if(cid()==CID_BLC2b) zcenter=13.75;
    TVector2 org(tempz,tempx);
    TVector2 tilted=org.Rotate(tilt()/180*TMath::Pi());
    tempx=tilted.Y();
    tempz=tilted.X();
  }

  if(GPOS){
    TVector2 gpos(gx(),gy());
    TVector2 rotgpos=gpos.Rotate(-theta);
    TVector2 gdir(dgx(),dgy());
    TVector2 rotgdir=gdir.Rotate(-theta);
    TVector2 org2(tempz,tempx);
    TVector2 dirmod=org2.Rotate(TMath::ATan(rotgdir.X()));
    tempx=dirmod.Y()+rotgpos.X();
    double gzpos=gz();
    if(cid()==CID_BLC2||cid()==CID_BLC2a||cid()==CID_BLC2b) gzpos += 130;
    tempz=dirmod.X()+gzpos;
  }
}

void ChamberLikeHit::gwpos2(double &tempx,double &tempz, 
			    const double &a,const double &b, 
			    const bool &TILT,const bool &YPOS)
{
  if(XY) tempx=wy(); else tempx=wx();
  tempz=wz();
  if(TILT){ // local tilt
    double zcenter;
    if(cid()==CID_BLC1a) zcenter=-16.;
    if(cid()==CID_BLC1b) zcenter=16.;
    if(cid()==CID_BLC2a) zcenter=-13.75;
    if(cid()==CID_BLC2b) zcenter=13.75;
    TVector2 org(tempz,tempx);
    TVector2 tilted=org.Rotate(tilt()/180*TMath::Pi());
    tempx=tilted.Y();
    tempz=tilted.X();
  }

  if(YPOS){
    double theta=rot();
    if(cid()==CID_BLC1a) theta-=45;
    if(cid()==CID_BLC1b) theta-=45;
    if(cid()==CID_BLC2a) theta-=135;
    if(cid()==CID_BLC2b) theta-=135;
    theta=theta/180*TMath::Pi();
    double tempy=a*tempz+b;
    TVector2 gpos(tempx,tempy);
    TVector2 rotgpos=gpos.Rotate(-theta);
    tempx=rotgpos.X();
  }
}
