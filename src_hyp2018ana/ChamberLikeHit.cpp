#include <iomanip>
#include <new>

#include "ChamberLikeHit.h"
#include "TVector2.h"

ClassImp(ChamberLikeHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ChamberLikeHit::ChamberLikeHit():
  CounterID(-1),HitID(-1),Layer(-1),Wire(-1),XY(-1),STATUS(1),
  DriftLength(-999),TimeOffset(-999),CorrDriftLength(-999),
  WirePos(DEFVECT),
  WireLength(-999.),TiltAngle(-999.),Resolution(-999.),
  HitPos(DEFVECT),
  dXY(-999),LR(-1),
  RotationAngle(-999),
  GPos(DEFVECT),GDir(DEFVECT)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ChamberLikeHit::Clear()
{
  STATUS=1;
  CounterID = HitID = XY = -1;
  Layer = Wire = -1;

  LR=-1;

  DriftLength = CorrDriftLength = -999;

  TimeOffset = 0.;

  WirePos.SetXYZ(-999,-999,-999);
  HitPos.SetXYZ(-999,-999,-999);
  GPos.SetXYZ(-999,-999,-999);
  GDir.SetXYZ(-999,-999,-999);

  WireLength = TiltAngle = Resolution = -999;  

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
  tdchit.SetHit(c,n,a,data);
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
    double drifttime = conf->GetXTMapManager()->CalcDriftTime( CounterID, Layer, Wire, DriftLength ) + TimeOffsetTrue;
    tdchit.SetTime(drifttime);
#if 0
    std::cout << " DriftLength:" << DriftLength << " ResolutionTrue:" << ResolutionTrue
	      << " DriftTime:" << DriftTime << std::endl;
#endif
  }
  tdchit.Reverse(conf);
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

    WirePos.SetZ( z );  //cm
    //    XY = xy;
    dXY = TMath::Abs(dxy);

    if( xy==0 ){
      XY = 0;
      WirePos.SetY(0);
      WirePos.SetX( xy0 + (Wire-1)*dxy );
    }
    else{
      XY = 1;
      WirePos.SetY( xy0 + (Wire-1)*dxy );
      WirePos.SetX( 0 );
    }
#if 0
    std::cout << " wx:" << WirePos.X() << "  xy: " << WirePos.Y() << "  wz: " << WirePos.Z() << std::endl;
#endif
    SetHitPosition( WirePos ); // just intial value
    conf->GetBLDCWireMapManager()->GetGParam( CounterID, GPos, GDir );
#if 0
    std::cout << " GX:" << GX << "  GY: " << GY << "  GZ: " << GZ << std::endl;
#endif
  }

  tdchit.Calc(conf);
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
    //    DriftTime -= TimeOffset;
    //    std::cout<< "DriftTimeNew: "<<DriftTime << std::endl;
    if( 0<tdc() ) DriftLength = conf->GetXTMapManager()->CalcDriftLength( CounterID, Layer, Wire, dt() - TimeOffset );
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
  tempx=wpos(XY);
  tempz=wz();
  theta=rot()+90.*xy()+dgz();
  if(TILT){ // local tilt
    TVector2 org(tempz,tempx);
    TVector2 tilted=org.Rotate(tilt()/180*TMath::Pi());
    tempx=tilted.Y();
    tempz=tilted.X();
  }
  if(GPOS){
    TVector2 gpos(gx(),gy());
    TVector2 rotgpos=gpos.Rotate(-theta);
    TVector2 gdir(dgx(),dgy());
    TVector2 rotgdir=gdir.Rotate(-theta+dgz());
    TVector2 org2(tempz,tempx);
    TVector2 dirmod=org2.Rotate(TMath::ATan(rotgdir.X()));
    tempx=dirmod.Y()*tempz+rotgpos.X();
    double gzpos=gz();
    if(cid()==CID_BLC2||cid()==CID_BLC2a||cid()==CID_BLC2b) gzpos += 130;
    tempz=dirmod.X()+gzpos;
  }
}
