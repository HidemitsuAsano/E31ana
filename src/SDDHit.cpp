// SDDHit.cpp

#include <iomanip>
#include <new>
#include <cstring>
#include "SDDHit.h"
#include <stdlib.h>
#include <time.h>
#include "TRandom.h"
#include "TGraph.h"

ClassImp(SDDHit);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SDDHit::SDDHit()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SDDHit::SDDHit( const SDDHit &hit )
{
  *this = hit; // I cannot understand why I need this copy constructor.
}
SDDHit::SDDHit( const double &rndm )
{
  RANDOM=rndm;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SDDHit::Clear()
{
  CrateOut = SlotOut = ChannelOut = -1;
  CrateXout = SlotXout = ChannelXout = -1;
  CrateFout = SlotFout = ChannelFout = -1;
  CrateTDCsdd = SlotTDCsdd = ChannelTDCsdd = -1;
  CrateTDCt0 = SlotTDCt0 = ChannelTDCt0 = -1;
  CrateTreset = SlotTreset = ChannelTreset = -1;
  CrateThigh = SlotThigh = ChannelThigh = -1;

  Seg = -1;
  Out = Fout = Xout = -1;
  TDCsdd = TDCt0 = -1;
  Treset = Thigh = -1;

  X=Y=Z=dX=dY=dZ=-999;
  GX=GY=dGX=dGY=dGZ=-999;
  Length=Width=Thick=-999;

  Gain=Offset=0.;
  Ene = Timesdd = Timet0 = -999;
  CTimesdd = CTimet0 = -999;
  CFout = -999;

  VGain=VOffset=VEne=0.;
  PADCl=PADCh=Nsample=-1;
  FBaseHeight=FBaseSlope=-1.;
  for(int i=0;i<MaxFADCPoint;i++) FADC[i]=-1;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SDDHit::CheckRange()
{
  if( 0<TDCsdd && TDCsdd<4095
      )
    return true;
  else
    return false;
}

bool SDDHit::SetData( const int &c, const int &n, const int &a,
		      const int &cid, const int &seg, const int &at, const int &ud,
		      const int &data )
{
#if 0
  std::cout << c << " " << n << " " << a << " " << cid << " " << seg << " " << at << " " << ud << " " << data 
	    << std::endl;
#endif
  CounterID = cid;
  Seg = seg;
  if( at==0 ){
    if( ud==0 ){ Out=data; CrateOut=c; SlotOut=n; ChannelOut=a; }
    else if ( ud==1 ){ Fout=data; CrateFout=c; SlotFout=n; ChannelFout=a; }
    else if ( ud==2 ){ Xout=data; CrateXout=c; SlotXout=n; ChannelXout=a; }
  }
  else{
    if( ud==0 ){ TDCsdd=data; CrateTDCsdd=c; SlotTDCsdd=n; ChannelTDCsdd=a; }
    else if(ud==1){   TDCt0=data; CrateTDCt0=c; SlotTDCt0=n; ChannelTDCt0=a; }
    else if(ud==2){   Treset=data; CrateTreset=c; SlotTreset=n; ChannelTreset=a; }
    else if(ud==3){   Thigh=data; CrateThigh=c; SlotThigh=n; ChannelThigh=a; }
  }  
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SDDHit::SetData( const int &seg, EventStruct &vme )
{
  Seg = seg;
  PADCl=vme.padcl[seg-1];
  PADCh=vme.padch[seg-1];
  Nsample=vme.nsample[seg-1];
  for(int i=0;i<vme.nsample[seg-1];i++) FADC[i]=vme.fadc[seg-1][i];
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool SDDHit::Calc( ConfMan *conf )
{
  RANDOM=gRandom->Uniform(-5,5);
  /*
v  if( conf->GetGeomMapManager() ){
    double x,y,z,dx,dy,dz,len,wid,thi,lv;
    if( conf->GetGeomMapManager()->GetParam( CounterID, Seg, x, y, z, dx, dy, dz, len, wid, thi, lv ) ){
      X = x; Y = y; Z = z; dX = dx; dY = dy; dZ = dz;
      Length = len; Width = wid; Thick = thi;
      //      LightVelocity = lv;
    }
    if( conf->GetGeomMapManager()->GetGParam( CounterID, x, y, z, dx, dy, dz ) ){
      GX = x; GY = y; GZ = z; dGX = dx; dGY = dy; dGZ = dz;
    }
  }
  */
  
  if( conf->GetFADCParamManager() && Nsample>0 && CheckRange() ){
    Interval=conf->GetFADCParamManager()->GetInterval();
    int time[MaxFADCPoint];
    for(int i=0;i<Nsample;i++){
      time[i]=i;
    }
    TGraph *tmpgr=new TGraph(Nsample,time,FADC);
    tmpgr->Fit("pol0","Q0","",
	       conf->GetFADCParamManager()->GetLowerLimitBase(),
	       conf->GetFADCParamManager()->GetUpperLimitBase());
    FBaseHeight=tmpgr->GetFunction("pol0")->GetParameter(0);
    tmpgr->Fit("pol1","Q0","",
	       conf->GetFADCParamManager()->GetLowerLimitBase(),
	       conf->GetFADCParamManager()->GetUpperLimitBase());
    FBaseSlope=tmpgr->GetFunction("pol1")->GetParameter(1);
    FBaseNpar=2;
    tmpgr->GetFunction("pol1")->GetParameters(FBasePar);
    tmpgr->Fit("expo","Q0","",
	       conf->GetFADCParamManager()->GetLowerLimitPost(),
	       conf->GetFADCParamManager()->GetUpperLimitPost());
    FPostSlope=tmpgr->GetFunction("expo")->GetParameter(1);
    FPostNpar=2;
    tmpgr->GetFunction("expo")->GetParameters(FPostPar);
    TF1 *fpeak=new TF1("fpeak",
		       conf->GetFADCParamManager()->GetPeakFunc().c_str(),
		       conf->GetFADCParamManager()->GetLowerLimitPeak(),
		       conf->GetFADCParamManager()->GetUpperLimitPeak());
    tmpgr->Fit(fpeak,"RQ0");
    FPeakNpar=fpeak->GetNpar();
    fpeak->GetParameters(FPeakPar);
    FPeakHeight=fpeak->GetMaximum();
    FPeakTime=fpeak->GetMaximumX()*Interval;
    delete fpeak;
    delete tmpgr;
  }
  
  if( conf->GetGainMapManager() ){
    if( 0<Out ){
      Ene  = conf->GetGainMapManager()->CalcCValue( CrateOut, SlotOut, ChannelOut, 0, Out+RANDOM );
      conf->GetGainMapManager()->GetParam( CrateOut,SlotOut, ChannelOut, Gain, Offset);
    }
    
    if( 0<PADCh ){
      VEne  = conf->GetGainMapManager()->CalcCValue( 9 , 1 , Seg , 0, PADCh+RANDOM );
      conf->GetGainMapManager()->GetParam( 9 , 1 , Seg , VGain, VOffset);
      FEne  = conf->GetGainMapManager()->CalcCValue( 9 , 2 , Seg , 0, FPeakHeight );
      conf->GetGainMapManager()->GetParam( 9 , 2 , Seg , FGain, FOffset);
    }
    
    if( 0<TDCsdd ) Timesdd  = conf->GetGainMapManager()->CalcCValue( CrateTDCsdd, SlotTDCsdd, ChannelTDCsdd, 1, TDCsdd );
    if( 0<TDCt0 ) Timet0  = conf->GetGainMapManager()->CalcCValue( CrateTDCt0, SlotTDCt0, ChannelTDCt0, 1, TDCt0 );
  }
  if( conf->GetSlewingMapManager() ){
    if( 0<TDCsdd ) CTimesdd  = conf->GetSlewingMapManager()->CalcCValue( CounterID, Seg, 1, Timesdd, Out );
    if( 0<TDCsdd ) CTimet0  = conf->GetSlewingMapManager()->CalcCValue( CounterID, Seg, 1, Timet0, Out );
    if( 0<Fout ) CFout  = conf->GetSlewingMapManager()->CalcCValue( CounterID, Seg, 2, Fout, Out );
  }
  return true;
}
