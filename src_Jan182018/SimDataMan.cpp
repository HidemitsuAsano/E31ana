// SimDataMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <new>

#include "SimDataMan.h"

ClassImp(SimDataMan);

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SimDataMan::SimDataMan()
{
  Conf = 0;
  CDSHit = new CDSHitMan();
  BeamLineHit = new BeamLineHitMan();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
SimDataMan::SimDataMan( ConfMan *conf )
{
  Conf = conf;
  CDSHit = new CDSHitMan();
  BeamLineHit = new BeamLineHitMan();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::SetHodoscopeLikeHit( int cid, int seg, double time, double dene, double hitpos )
{
  int cds_bl=-1; //cdsMan:0 ,blMan:1

  if( cid==CID_CDH ) cds_bl=0;
  else if(cid==CID_BHD ) cds_bl=1;
  else if(cid==CID_PA ) cds_bl=1;
  else if(cid==CID_T0 ) cds_bl=1;
  else if(cid==CID_LC1 ) cds_bl=1;
  else if(cid==CID_LC2 ) cds_bl=1;
  else if(cid==CID_TOFstop ) cds_bl=1;
  else if(cid==CID_NC ) cds_bl=1;
  else if(cid==CID_TOFstop ) cds_bl=1;
  else if(cid==CID_PC ) cds_bl=1;
  else return;

//  if(cds_bl==0)
//    {
//      if( CDSHit==0 )
//	CDSHit = new CDSHitMan();
//    }
//  if(cds_bl==1)
//    {
//      if( BeamLineHit==0 )
//	BeamLineHit = new BeamLineHitMan();
//    }

  HodoscopeLikeHit *hit = new HodoscopeLikeHit();
  hit->SetCounterID( cid );
  hit->SetSegment( seg );
  hit->SetCTMean(time);
  hit->SetEMean(dene);
  hit->SetHitPosition(hitpos);

  if( Conf!=0 ){
    if( Conf->GetCounterMapManager() ){
      int c,n,a;
      int at, ud;
      at=0;ud=0; 
      Conf->GetCounterMapManager()->GetCNA( cid, seg, at, ud, c, n, a );
      hit->SetCrateAu(c); hit->SetSlotAu(n); hit->SetChannelAu(a);
      at=0;ud=1;
      Conf->GetCounterMapManager()->GetCNA( cid, seg, at, ud, c, n, a );
      hit->SetCrateAd(c); hit->SetSlotAd(n); hit->SetChannelAd(a);
      at=1;ud=0;
      Conf->GetCounterMapManager()->GetCNA( cid, seg, at, ud, c, n, a );
      hit->SetCrateTu(c); hit->SetSlotTu(n); hit->SetChannelTu(a);
      at=1;ud=1;
      Conf->GetCounterMapManager()->GetCNA( cid, seg, at, ud, c, n, a );
      hit->SetCrateTd(c); hit->SetSlotTd(n); hit->SetChannelTd(a);
    }
    else{
      std::cout << " counter map was not found " << std::endl;
    }

    if( Conf->GetGeomMapManager() ){
      double x,y,z,dx,dy,dz,len,wid,th,lv;
      Conf->GetGeomMapManager()->GetParam( cid, seg, x, y, z, dx, dy, dz, len, wid, th, lv );
      hit->SetPos(x,y,z); hit->SetDir(dx,dy,dz); hit->SetLength(len); hit->SetWidth(wid); hit->SetThick(th); hit->SetLightVelocity(lv);
      hit->SetCTSub(hitpos/lv);
    }
    else{
      std::cout << " geom map was not found " << std::endl;
    }
    hit->SetSimulatedResolution(Conf);
    hit->Reverse(Conf);
    double atime=hit->ctmean();
    double adene=hit->emean();
    double ahitpos=hit->hitpos();
    hit->Calc(Conf);
    hit->SetCTMean(atime);
    hit->SetEMean(adene);
    hit->SetHitPosition(ahitpos);

  }
  if(cds_bl==0)  CDSHit->AddHit(*hit);
  if(cds_bl==1)  BeamLineHit->AddHit(*hit);
  delete hit;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::SetChamberLikeHit( int cid, int layer, int wire, double dl, double toffset )
{
  int cds_bl=-1; //cdsMan:0 ,blMan:1

  if( cid==CID_CDC ) cds_bl=0;
  else if(cid==CID_PDC1 ) cds_bl=1;
  else if(cid==CID_PDC2 ) cds_bl=1;
  else if(cid==CID_BLC1 ) cds_bl=1;
  else if(cid==CID_BLC2 ) cds_bl=1;
  else return;

//  if(cds_bl==0)
//    {
//      if( CDSHit==0 )
//	CDSHit = new CDSHitMan();
//    }
//  if(cds_bl==1)
//    {
//      if( BeamLineHit==0 )
//	BeamLineHit = new BeamLineHitMan();
//    }

  if(cid==CID_CDC)
    {    
      CDCHit *hit = new CDCHit();
      hit->SetCounterID( cid );
      hit->SetLayer( layer );
      hit->SetWire( wire );
      hit->SetDriftLength( dl );
      hit->SetTimeOffsetTrue( toffset );

      if( Conf!=0 ){
	if(cid==CID_CDC && Conf->GetCDCWireMapManager() ){
	  int c,n,a;
	  int slayer, asdnum, ttype, asdch;
	  double rad, phi, tilt;
	  Conf->GetCDCWireMapManager()->GetWire( layer, wire, slayer, asdnum, ttype, asdch,
						 rad, phi, tilt, c, n, a );
	  hit->SetCrate(c); hit->SetSlot(n); hit->SetChannel(a);
	  hit->SetSuperLayer(slayer); hit->SetASDNum(asdnum); hit->SetTransType(ttype); hit->SetASDChannel(asdch);
	  hit->SetRadius(rad); hit->SetPhi(phi); hit->SetTiltAngle(tilt);
	  hit->SetWirePosition(rad*cos(phi/180.*3.141592),rad*sin(phi/180.*3.141592),0);
	}
	else{
	  std::cout << " wire map was not found " << std::endl;
	}
	hit->SetSimulatedResolution(Conf);
	hit->Reverse( Conf );
	double adl=hit->dl();
	hit->Calc(Conf);
	hit->SetDriftLength(adl);
      }

      CDSHit->AddHit(*hit);
      delete hit;

    }
  else  
  {
    ChamberLikeHit *hit = new ChamberLikeHit();
    hit->SetCounterID( cid );
    hit->SetLayer( layer );
    hit->SetWire( wire );
    hit->SetDriftLength( dl );
    hit->SetTimeOffsetTrue( toffset );
    int c,n,a;
    if( Conf!=0 ){
      if(  Conf->GetCounterMapManager()->GetCNA( cid, layer, wire, 1, 0, c, n,a ) )
	{     
	  //std::cout << " c n a "<<c<<" "<<n<<" "<<a<< std::endl;
	  hit->SetCrate(c); hit->SetSlot(n); hit->SetChannel(a);
	}
      else{
	std::cout << " wire map was not found " << std::endl;
      }
      hit->SetSimulatedResolution(Conf);
      hit->Reverse( Conf );
      double adl=hit->dl();
      hit->Calc(Conf);
      hit->SetDriftLength(adl);
    }

    //    if(cds_bl==0)  CDSHit->AddHit(*hit);
    if(cds_bl==1)  BeamLineHit->AddHit(*hit);
    delete hit;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::SetCDHHit( int seg, double time, double dene, double hitpos )
{
//  if( CDSHit==0 ){
//    CDSHit = new CDSHitMan();
//  }

  HodoscopeLikeHit *hit = new HodoscopeLikeHit();
  hit->SetCounterID( CID_CDH );
  hit->SetSegment( seg );
  hit->SetCTMean(time);
  hit->SetEMean(dene);
  hit->SetHitPosition(hitpos);

  if( Conf!=0 ){
    if( Conf->GetCounterMapManager() ){
      int c,n,a;
      int at, ud;
      at=0;ud=0; 
      Conf->GetCounterMapManager()->GetCNA( CID_CDH, seg, at, ud, c, n, a );
      hit->SetCrateAu(c); hit->SetSlotAu(n); hit->SetChannelAu(a);
      at=0;ud=1;
      Conf->GetCounterMapManager()->GetCNA( CID_CDH, seg, at, ud, c, n, a );
      hit->SetCrateAd(c); hit->SetSlotAd(n); hit->SetChannelAd(a);
      at=1;ud=0;
      Conf->GetCounterMapManager()->GetCNA( CID_CDH, seg, at, ud, c, n, a );
      hit->SetCrateTu(c); hit->SetSlotTu(n); hit->SetChannelTu(a);
      at=1;ud=1;
      Conf->GetCounterMapManager()->GetCNA( CID_CDH, seg, at, ud, c, n, a );
      hit->SetCrateTd(c); hit->SetSlotTd(n); hit->SetChannelTd(a);
    }
    else{
      std::cout << " counter map was not found " << std::endl;
    }

    if( Conf->GetGeomMapManager() ){
      double x,y,z,dx,dy,dz,len,wid,th,lv;
      Conf->GetGeomMapManager()->GetParam( CID_CDH, seg, x, y, z, dx, dy, dz, len, wid, th, lv );
      hit->SetPos(x,y,z); hit->SetDir(dx,dy,dz); hit->SetLength(len); hit->SetWidth(wid); hit->SetThick(th); hit->SetLightVelocity(lv);
      hit->SetCTSub(hitpos/lv);
    }
    else{
      std::cout << " geom map was not found " << std::endl;
    }
    hit->SetSimulatedResolution(Conf);
    hit->Reverse(Conf);

  }
  CDSHit->AddHit(*hit);
  delete hit;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::SetCDCHit( int layer, int wire, double dl, double toffset )
{
//  if( CDSHit==0 ){
//    CDSHit = new CDSHitMan();
//  }

  CDCHit *hit = new CDCHit();
  hit->SetCounterID( CID_CDC );
  hit->SetLayer( layer );
  hit->SetWire( wire );
  hit->SetDriftLength( dl );
  hit->SetTimeOffsetTrue( toffset );

  if( Conf!=0 ){
    if( Conf->GetCDCWireMapManager() ){
      int c,n,a;
      int slayer, asdnum, ttype, asdch;
      double rad, phi, tilt;
      Conf->GetCDCWireMapManager()->GetWire( layer, wire, slayer, asdnum, ttype, asdch,
					     rad, phi, tilt, c, n, a );
      hit->SetCrate(c); hit->SetSlot(n); hit->SetChannel(a);
      hit->SetSuperLayer(slayer); hit->SetASDNum(asdnum); hit->SetTransType(ttype); hit->SetASDChannel(asdch);
      hit->SetRadius(rad); hit->SetPhi(phi); hit->SetTiltAngle(tilt);
      hit->SetWirePosition(rad*cos(phi/180.*3.141592),rad*sin(phi/180.*3.141592),0);
    }
    else{
      std::cout << " cdc wire map was not found " << std::endl;
    }
    hit->SetSimulatedResolution(Conf);
    hit->Reverse( Conf );
  }

  CDSHit->AddHit(*hit);
  delete hit;

}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::SetBLCHit( int tag, int layer, int wire, double dl, double toffset )
{

//  if( BeamLineHit==0 ){
//    BeamLineHit = new BeamLineHitMan();
//  }

  int CID_BLC=CID_BLC1;
  if(tag==0) CID_BLC=CID_BLC1;
  else if(tag==1) CID_BLC=CID_BLC2;

  ChamberLikeHit *hit = new ChamberLikeHit();
  hit->SetCounterID( CID_BLC );
  hit->SetLayer( layer );
  hit->SetWire( wire );
  hit->SetDriftLength( dl );
  hit->SetTimeOffsetTrue( toffset );

  if( Conf!=0 ){
    int c=1,n=1,a=1; // temporal
    if(Conf->GetCounterMapManager()->GetCNA( CID_BLC, layer, wire, 1, 0, c, n,a ) )
      {     
	//std::cout << " c n a "<<c<<" "<<n<<" "<<a<< std::endl;
	hit->SetCrate(c); hit->SetSlot(n); hit->SetChannel(a);
      }
    //********* 2011 10 4 #########//

    //    if( Conf->GetBLDCWireMapManager() ){
    //int c=1,n=1,a=1; // temporal
      //hit->SetCrate(c); hit->SetSlot(n); hit->SetChannel(a); // temporal

      // *** not understand set/get c,n,a to BLC data 2011/08/12 ***
     
    //int slayer, asdnum, ttype, asdch;
    //double rad, phi, tilt;
     //     std::cout << " cid layer wire "<<CID_BLC<<" "<<layer<<" "<<wire<< std::endl;
     //     Conf->GetBLDCWireMapManager()->GetWire( layer, wire, slayer, asdnum, ttype, asdch,
     //				     rad, phi, tilt, c, n, a );
     //hit->SetCrate(c); hit->SetSlot(n); hit->SetChannel(a);
     //     hit->SetSuperLayer(slayer); hit->SetASDNum(asdnum); hit->SetTransType(ttype); hit->SetASDChannel(asdch);
     // hit->SetRadius(rad); hit->SetPhi(phi); hit->SetTiltAngle(tilt);
     //hit->SetWirePosition(rad*cos(phi/180.*3.141592),rad*sin(phi/180.*3.141592),0);
    //}
    else{
      std::cout << " blc wire map was not found " << std::endl;
    }
    hit->SetSimulatedResolution(Conf);
    hit->Reverse( Conf );
  }

  BeamLineHit->AddHit(*hit);
  delete hit;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::SetT0Hit( int seg, double time, double dene, double hitpos )
{
//  if( BeamLineHit==0 ){
//    BeamLineHit = new BeamLineHitMan();
//  }

  HodoscopeLikeHit *hit = new HodoscopeLikeHit();
  hit->SetCounterID( CID_T0 );
  hit->SetSegment( seg );
  hit->SetCTMean(time);
  hit->SetEMean(dene);
  hit->SetHitPosition(hitpos);

  if( Conf!=0 ){
    if( Conf->GetCounterMapManager() ){
      int c,n,a;
      int at, ud;
      at=0;ud=0; 
      Conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
      hit->SetCrateAu(c); hit->SetSlotAu(n); hit->SetChannelAu(a);
      at=0;ud=1;
      Conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
      hit->SetCrateAd(c); hit->SetSlotAd(n); hit->SetChannelAd(a);
      at=1;ud=0;
      Conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
      hit->SetCrateTu(c); hit->SetSlotTu(n); hit->SetChannelTu(a);
      at=1;ud=1;
      Conf->GetCounterMapManager()->GetCNA( CID_T0, seg, at, ud, c, n, a );
      hit->SetCrateTd(c); hit->SetSlotTd(n); hit->SetChannelTd(a);
    }
    else{
      std::cout << " counter map was not found " << std::endl;
    }

    if( Conf->GetGeomMapManager() ){
      double x,y,z,dx,dy,dz,len,wid,th,lv;
      Conf->GetGeomMapManager()->GetParam( CID_T0, seg, x, y, z, dx, dy, dz, len, wid, th, lv );
      hit->SetPos(x,y,z); hit->SetDir(dx,dy,dz); hit->SetLength(len); hit->SetWidth(wid); hit->SetThick(th); hit->SetLightVelocity(lv);
      hit->SetCTSub(hitpos/lv);
    }
    else{
      std::cout << " geom map was not found " << std::endl;
    }
    hit->SetSimulatedResolution(Conf);
    hit->Reverse(Conf);
  }
  BeamLineHit->AddHit(*hit);
  delete hit;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void SimDataMan::Clear()
{
  CDSHit->Clear();
  if( CDSHit!=0 ){    
    delete CDSHit;
    CDSHit=0;
  }
  BeamLineHit->Clear();
  if( BeamLineHit!=0 ){
    delete BeamLineHit;
    BeamLineHit=0;
  }
}


