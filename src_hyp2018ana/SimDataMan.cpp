// SimDataMan.cpp
#include <string>
#include <cstdio>
#include <iostream>
#include <new>

using namespace std;

#include "SimDataMan.h"

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SimDataMan::SimDataMan()
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::Convert(DetectorData *detectorData, ConfMan *confMan, BeamLineHitMan *blMan, CDSHitMan *cdsMan)
{
  int nCDH=0, nIH=0, nDEF=0, nT0=0, nBPD=0, nBVC=0, nNC=0, nCVC=0, nPC=0;
  int nCDC[NumOfCDCLayers]; for( int i=0; i<NumOfCDCLayers; i++ ) nCDC[i]=0;
  int nBLC2a[8], nBLC2b[8], nBPC[8], nFDC1[8];
  for( int i=0; i<8; i++ ){ nBLC2a[i]=0; nBLC2b[i]=0; nBPC[i]=0; nFDC1[i]=0; }

  for(  int i=0; i<detectorData->detectorHitSize(); i++ ){
    DetectorHit *mchit = detectorData-> detectorHit(i);
    int cid = mchit->detectorID();
    int channel   = mchit->channelID();
    int pdg_id    = mchit->pdg();
    int parent_id = mchit->parentID();
    double time = mchit->time();
    double dE = mchit->de();
    TVector3 pos = mchit->pos();

    if( pdg_id==321 && parent_id==0 ) time *= -1;

    if( cid==CID_CDH || cid==CID_IH || cid==CID_BHD || cid==CID_DEF || cid==CID_T0 || cid==CID_BPD ||
        cid==CID_LC1 || cid==CID_LC2 || cid==CID_BVC || cid==CID_NC || cid==CID_CVC || cid==CID_PC ){

      HodoscopeLikeHit hit;
      hit.SetCounterID(cid);
      hit.SetSegment(channel+1);

      TVector3 hitposition;
      if( confMan->GetGeomMapManager()->GetGPos( cid, channel+1 , hitposition) ){
	hit.SetPos(hitposition);
      }
      else{
	cout<<"!!!!! Hit Position not found !!!!!"<<endl;
	exit(0);
      }


      int nsenser=2;
      if( cid==CID_IH ){ nsenser=1; hit.SetNSensor(1); }
      for( int at=0; at<2; at++ ){
	for( int i=0; i<nsenser; i++ ){
	  int c, n, a;
	  if( !confMan-> GetCounterMapManager()-> GetCNA(cid, channel+1, at, i, c, n, a) ){
	    std::cout<<"  !!! SimDataMan Error !!!"<<std::endl;
	    exit(0);
	  }
	  hit.SetCrate(at, i, c); hit.SetSlot(at, i, n); hit.SetChannel(at, i, a);
	}
      }
      double lv=DBL_MIN; TVector3 gpos;
      confMan-> GetGeomMapManager()-> GetGPos(cid, channel+1, lv, gpos);
      double hitpos;
      if( cid==CID_CDH || cid==CID_IH ) hitpos=pos.Z()/10.;
      else if( cid==CID_BVC ) hitpos=pos.X()/10.;
      else hitpos=pos.Y()/10.;

      hit.SetCTMean(time);
      hit.SetEMean(dE);
      hit.SetHitPosition(hitpos);
      hit.SetTSub(hitpos/lv);
      hit.SetSimulatedResolution(confMan);
      hit.SetHitPosition(hit.hitpos()+hit.dt()*lv);

      // if( confMan->GetGateMapManager()->GetParam(cid, channel+1, 0, 0, type, low, high) && type==2 ){
      // 	if( dE<low ) continue;
      // }
      // else if( dE<0.1 ) continue;
      if( dE<0.1 ) continue;

      if( cid==CID_CDH      ) nCDH++;
      else if( cid==CID_IH  ) nIH++;
      else if( cid==CID_BPD ) nBPD++;
      else if( cid==CID_DEF ) nDEF++;
      else if( cid==CID_T0  ) nT0++;
      else if( cid==CID_BVC ) nBVC++;
      else if( cid==CID_NC  ) nNC++;
      else if( cid==CID_CVC ) nCVC++;
      else if( cid==CID_PC  ) nPC++;

      if( cid==CID_CDH || cid==CID_IH ) cdsMan-> AddHit(hit);
      else blMan-> AddHit(hit);
    }
    else if( cid==CID_CDC ){
      if( pdg_id==22 || pdg_id==2112 ) continue;
      int layer = mchit->layerID()+1;
      int wire = mchit->channelID()+1;
      double dl=fabs(mchit->dx()/10.);

      int c,n,a;
      int slayer, asdnum, ttype, asdch;
      double rad, phi, tilt;
      confMan->GetCDCWireMapManager()->GetWire( layer, wire, slayer, asdnum, ttype, asdch,
						rad, phi, tilt, c, n, a );
      CDCHit hit;
      hit.SetCounterID( cid );
      hit.SetLayer( layer );
      hit.SetWire( wire );
      hit.SetDriftLength( dl );
      hit.SetTimeOffsetTrue( 0.0 );
      hit.SetCrate(c); hit.SetSlot(n); hit.SetChannel(a);
      hit.SetSuperLayer(slayer); hit.SetASDNum(asdnum); hit.SetTransType(ttype); hit.SetASDChannel(asdch);
      hit.SetRadius(rad); hit.SetPhi(phi); hit.SetTiltAngle(tilt);
      hit.SetWirePosition(rad*cos(phi/180.*3.141592),rad*sin(phi/180.*3.141592),0);

      hit.SetSimulatedResolution(confMan);
      hit.Reverse( confMan );
      double adl=hit.dl();
      hit.Calc(confMan);
      hit.SetDriftLength(adl);
      hit.SetCorrDriftLength(mchit->momentum().Pt());
      cdsMan->AddHit(hit);

      nCDC[layer-1]++;
    }
    else if( cid==CID_BPC || cid==CID_FDC1 || cid==CID_BLC2a || cid==CID_BLC2b || cid==CID_BLC1a || cid==CID_BLC1b){
      if( pdg_id==22 || pdg_id==2112 ) continue;
      int layer = mchit->layerID()+1;
      int wire = mchit->channelID()+1;
      double dl=fabs(mchit->dx()/10.);
      int c,n,a;
      confMan->GetCounterMapManager()->GetCNA( cid, layer, wire, 1, 0, c, n,a );

      ChamberLikeHit hit;
      hit.SetCounterID( cid );
      hit.SetLayer( layer );
      hit.SetWire( wire );
      hit.SetDriftLength( dl );
      hit.SetTimeOffsetTrue( 0.0 );
      hit.SetCrate(c); hit.SetSlot(n); hit.SetChannel(a);

      hit.SetSimulatedResolution(confMan);
      hit.Reverse(confMan );
      double adl=hit.dl();
      hit.Calc(confMan);
      hit.SetDriftLength(adl);
      hit.SetCorrDriftLength(mchit->momentum().Pt());
      blMan->AddHit(hit);

      if( cid==CID_BLC2a ) nBLC2a[layer-1]++;
      if( cid==CID_BLC2b ) nBLC2b[layer-1]++;    
      if( cid==CID_BPC   ) nBPC[layer-1]++;
      if( cid==CID_FDC1  ) nFDC1[layer-1]++;    
    }
  }

#if 0
  cout<<" T0  DEF BPD BVC CVC NC  PC"<<endl;
  cout<<" "<<nT0<<"   "<<nDEF<<"   "<<nBPD<<"   "<<nBVC<<"   "<<nCVC<<"   "<<nNC<<"   "<<nPC<<endl;
  cout<<"BLC2a : ";
  for( int i=0; i<8; i++ ) cout<<nBLC2a[i]<<" ";
  cout<<endl;

  cout<<"BLC2b : ";
  for( int i=0; i<8; i++ ) cout<<nBLC2b[i]<<" ";
  cout<<endl;

  cout<<"BPC   : ";
  for( int i=0; i<8; i++ ) cout<<nBPC[i]<<" ";
  cout<<endl;

  cout<<"FDC1  : ";
  for( int i=0; i<6; i++ ) cout<<nFDC1[i]<<" ";
  cout<<endl;

  cout<<"CDH : "<<nCDH<<endl;
  cout<<"CDC : ";
  for( int i=0; i<NumOfCDCLayers; i++ ) cout<<nCDC[i]<<" ";
  cout<<endl;
#endif
}

