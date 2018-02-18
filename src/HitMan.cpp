#include "HitMan.h"

ClassImp(HitMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HitMan::HitMan() : TObject()
{
}


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void HitMan::AddHit( const HodoscopeLikeHit &hit )
{
  int cid = hit.cid();
  HodoscopeLikeContainer *container=HodoContainer(cid);
  container->push_back(hit);
  (*container)[ container->size() - 1].SetHitID( container->size() - 1 );
  HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = container->size() -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void HitMan::AddHit( const CherenkovLikeHit &hit )
{
  int cid = hit.cid();
  CherenkovLikeContainer *container=ChereContainer(cid);
  //  if(!container) return;
  container->push_back(hit);
  (*container)[ container->size() - 1].SetHitID( container->size() - 1 );
  CherenkovSegContainer[ KEY(hit.cid(),hit.seg()) ] = container->size() -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void HitMan::AddHit( const ChamberLikeHit &hit )
{
  int cid = hit.cid();
  int layer = hit.layer();
  ChamberLikeContainer *container=ChmContainer(cid,layer);
  if(container){
    container->push_back(hit);
    (*container)[ container->size() - 1].SetHitID( container->size() - 1 );
    return;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void HitMan::AddHit( const CDCHit &hit )
{
  int cid = hit.cid();
  int layer = hit.layer();
  CDCLikeContainer *container=CDCHitContainer(cid,layer);
  container->push_back(hit);
  (*container)[ container->size() - 1].SetHitID( container->size() - 1 );
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void HitMan::AddHit( const MTDCHit &hit )
{
  int cid = hit.cid();
  MTDCLikeContainer *container=MTDCContainer(cid);
  container->push_back(hit);
  (*container)[ container->size() - 1].SetHitID( container->size() - 1 );
  HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = container->size() -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::SetData( const int &c, const int &n, const int &a,
			      const int &cid, const int &seg, const int &at, const int &ud,
			      const int &data, HodoscopeLikeContainer *container )
{
  HodoscopeLikeHit *hit = 0;
  bool pabsent=false;
  for( int i=0; i<(int)container->size(); i++ ){
    if( (*container)[i].seg() == seg ){ hit = &(*container)[i]; }
  }
  if( hit==0 ){ hit = new HodoscopeLikeHit(); pabsent=true; }
  hit->SetData( c, n, a, cid, seg, at, ud, data );
  if(cid==CID_IH) hit->SetNSensor(1);
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::SetData( const int &c, const int &n, const int &a,
			      const int &cid, const int &lay, const int &wire,
			      const int &data, ChamberLikeContainer *container )
{
  ChamberLikeHit *hit = 0;
  bool pabsent=false;
  for( int i=0; i<(int)container->size(); i++ ){
    // Now, Chambers have only one TDC information.
    // And, There will be multi-hit TDC.
    //if( (*container)[i].seg() == seg ){ hit = &(*container)[i]; }
  }
  if( hit==0 ){ hit = new ChamberLikeHit(); pabsent=true; }
  hit->SetData( c, n, a, cid, lay, wire, data );
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::SetData( const int &c, const int &n, const int &a,
			      const int &cid, const int &lay, const int &wire,
			      const int &data, CDCLikeContainer *container )
{
  CDCHit *hit = 0;
  bool pabsent=false;
  for( int i=0; i<(int)container->size(); i++ ){
    // Now, Chambers have only one TDC information.
    // And, There will be multi-hit TDC.
    //if( (*container)[i].seg() == seg ){ hit = &(*container)[i]; }
  }
  if( hit==0 ){ hit = new CDCHit(); pabsent=true; }
  hit->SetData( c, n, a, cid, lay, wire, data );
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::SetData( const int &c, const int &n, const int &a,
			      const int &cid, const int &seg, const int &at, const int &ud,
			      const int &data, MTDCLikeContainer *container )
{
  MTDCHit *hit = 0;
  bool pabsent=false;
  for( int i=0; i<(int)container->size(); i++ ){
    if( (*container)[i].seg() == seg ){ hit = &(*container)[i]; }
  }
  if( hit==0 ){ hit = new MTDCHit(); pabsent=true; }
  hit->SetData( c, n, a, cid, seg, ud, data );
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::SetData( const int &c, const int &n, const int &a,
			      const int &cid, const int &seg, const int &at, const int &ud,
			      const int &data, CherenkovLikeContainer *container )
{
  CherenkovLikeHit *hit = 0;
  bool pabsent=false;
  for( int i=0; i<(int)container->size(); i++ ){
    if( (*container)[i].seg() == seg ){ hit = &(*container)[i]; }
  }
  if( hit==0 ){ hit = new CherenkovLikeHit(); pabsent=true; }
  hit->SetData( c, n, a, cid, seg, at, ud, data );
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::Convert( TKOHitCollection *tko, ConfMan *conf, double toffs )
{
#if 0
  std::cout << " Enter HitMan::Convert(TKO,conf) " << std::endl;
#endif
  for( int ih=0; ih<tko->entries(); ih++ ){
    int c, n, a, data;
    c = tko->hit(ih)->cr(); n = tko->hit(ih)->sl(); a = tko->hit(ih)->ch(); data = tko->hit(ih)->data();
    Convert( c, n, a, data, conf );
  }
  Calc( conf, toffs );
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool HitMan::Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf )
{
#if 0
  std::cout << " Enter HitMan::Convert(c,n,a,data,conf) " << std::endl;
#endif
  int cid, lay, seg, at, ud;
  conf->GetCounterMapManager()->GetInfo( c, n, a, cid, lay, seg, at, ud );

#if 0
  std::cout << " c:" << c << " n:" << n << " a:" << a << " cid:" << cid << " lay:" << lay << " seg:" << seg
	    << " at:" << at << " ud:" << ud << " data:" << data << std::endl;
#endif  
  if(cid==CID_CDC){
    int wire;
    if( conf->GetCDCWireMapManager()->GetWire( c, n, a, lay, wire ) ){
      CDCLikeContainer *chmcon=CDCHitContainer(cid,lay);
      if(chmcon)
	SetData( c, n, a, cid, lay, wire, data, chmcon);
    }
  }else{    
    HodoscopeLikeContainer *hodocon=HodoContainer(cid);
    if(hodocon){
      SetData( c, n, a, cid, seg, at, ud, data, hodocon );
    }else{
      ChamberLikeContainer *chmcon=ChmContainer(cid,lay);
      if(chmcon){
	SetData( c, n, a, cid, lay, seg, data, chmcon );
      }else{
	CherenkovLikeContainer *cherecon=ChereContainer(cid);
	if(cherecon){
	  SetData( c, n, a, cid, seg, at, ud, data, cherecon );
	}
      }
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeHit * HitMan::Hodo( const int &cid, const int &seg )
{
  unsigned int key;
  key = KEY(cid,seg);
  unsigned int i;
  SegmentContainer::iterator is = HodoscopeSegContainer.find(key);
  if( is != HodoscopeSegContainer.end()){
    i=is->second;
  }else{
    return 0;
  }
  if(HodoContainer(cid))
    return (i<HodoContainer(cid)->size()) ? &(*HodoContainer(cid))[i] : 0;
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int HitMan::nHodo( const int &cid )
{
  if(HodoContainer(cid))
    return HodoContainer(cid)->size();
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeHit * HitMan::Hodoi( const int &cid, const unsigned int &i )
{
  HodoscopeLikeContainer *container=HodoContainer(cid);
  return (i<container->size()) ? &(*container)[i] : 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int HitMan::nChm( const int &cid, const int &layer){
  ChamberLikeContainer *container= ChmContainer(cid,layer);
  if(container)
    return container->size();
  CDCLikeContainer *container1= CDCHitContainer(cid,layer);
  if(container1)
    return container1->size();
  return -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
ChamberLikeHit* HitMan::Chm( const int &cid, const int &layer, const int &i){
  ChamberLikeContainer *container= ChmContainer(cid,layer);
  if(container)
    return &(*container)[i];
  CDCLikeContainer *container1= CDCHitContainer(cid,layer);
  if(container1)
    return &(*container1)[i];
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CherenkovLikeHit * HitMan::Chere( const int &cid, const int &seg )
{
  unsigned int i = CherenkovSegContainer[ KEY(cid,seg) ];
  CherenkovLikeContainer *container=ChereContainer(cid);
  return (i<container->size()) ? &(*container)[i] : 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int HitMan::nChere( const int &cid )
{
  CherenkovLikeContainer *container=ChereContainer(cid);
  if(container)
    return container->size();
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CherenkovLikeHit* HitMan::Cherei( const int &cid, const unsigned int &i )
{
  CherenkovLikeContainer *container=ChereContainer(cid);
  return (i<container->size()) ? &(*container)[i] : 0;
}

void HitMan::Clear()
{
  HodoscopeSegContainer.clear();
  CherenkovSegContainer.clear();
}
