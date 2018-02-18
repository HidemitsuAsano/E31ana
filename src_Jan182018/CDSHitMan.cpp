#include "CDSHitMan.h"

ClassImp(CDSHitMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CDSHitMan::CDSHitMan() : TObject()
{
  Clear();
  //  CheckContainerSize();
}

const unsigned int KEYMASK  = 0x000F;
const unsigned int SEGMASK    = 0x01FF;      /* SEG Mask 9 Bits (0-511) */
const unsigned int CIDMASK    = 0x007F;      /* CID Mask 7 Bits (0-127) */
const int          SEGSHIFT   = 8;
const int          CIDSHIFT   = 24;
const unsigned int KEYFLAG  = 0x0003;

#define KEY(cid,seg) \
((((cid)&CIDMASK)<<CIDSHIFT) | (((seg)&SEGMASK)<<SEGSHIFT) | KEYFLAG )

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void CDSHitMan::AddHit( const HodoscopeLikeHit &hit )
{
  int cid = hit.cid();
  switch( cid ){
  case CID_CDH:
    CDHContainer.push_back(hit);
    CDHContainer[ CDHContainer.size()-1 ].SetHitID( CDHContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = CDHContainer.size()-1;
    break;
  case CID_IH:
    IHContainer.push_back(hit);
    IHContainer[ IHContainer.size()-1 ].SetHitID( IHContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = IHContainer.size()-1;
    break;
  default:
    break;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void CDSHitMan::AddHit( const CDCHit &hit )
{
  int cid = hit.cid();
  int layer = hit.layer();
#if 0
  std::cout << " in Add Hit, layer:" << hit.layer() << " wire:" << hit.wire() << " CID:" << cid << std::endl;
#endif
  switch( cid ){
  case CID_CDC:
    CDCContainer[layer-1].push_back(hit);
    CDCContainer[layer-1][ CDCContainer[layer-1].size()-1 ].SetHitID( CDCContainer[layer-1].size()-1 );
#if 0
    std::cout << "?? wire:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].wire()
	      << " hid:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].hid()
	      << " slayer:" << CDCContainer[layer-1][CDCContainer[layer-1].size()-1].slayer()
	      << std::endl;
#endif
    break;
  default:
    break;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::Convert( TKOHitCollection *tko, ConfMan *conf )
{
#if 0
  std::cout << " Enter CDSHitMan::Convert(TKO,conf) " << std::endl;
#endif
  //CheckContainerSize();
  for( int ih=0; ih<tko->entries(); ih++ ){
    int c, n, a, data;
    c = tko->hit(ih)->cr(); n = tko->hit(ih)->sl(); a = tko->hit(ih)->ch(); data = tko->hit(ih)->data();
    Convert( c, n, a, data, conf );
  }
  //CheckContainerSize();

  Calc( conf );

#if 0
  std::cout << "AAAAA" << std::endl;
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    std::cout << "Layer" << layer << "  #ev=" << CDCContainer[layer-1].size() << std::endl;
    for( int i=0; i<CDCContainer[layer-1].size(); i++ ){
      int wire = CDCContainer[layer-1][i].wire();
      int tdc = CDCContainer[layer-1][i].tdc();
      int slayer = CDCContainer[layer-1][i].slayer();
      std::cout << "   wire=" << wire << " tdc=" << tdc << " slayer=" << slayer << std::endl;
    }      
  }
#endif
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf )
{
#if 0
  std::cout << " Enter CDSHitMan::Convert(c,n,a,data,conf) " << std::endl;
#endif
  int cid, lay, seg, at, ud;
  conf->GetCounterMapManager()->GetInfo( c, n, a, cid, lay, seg, at, ud );

#if 0
  std::cout << " c:" << c << " n:" << n << " a:" << a << " cid:" << cid << " lay:" << lay << " seg:" << seg
	    << " at:" << at << " ud:" << ud << " data:" << data << std::endl;
#endif

  //CheckContainerSize();
  
  switch( cid ){
  case CID_CDH:
    SetData( c, n, a, cid, seg, at, ud, data, &CDHContainer );
    break;
  case CID_IH:
    SetData( c, n, a, cid, seg, at, ud, data, &IHContainer );
    break;
  case CID_CDC:
    int wire;
    if( conf->GetCDCWireMapManager()->GetWire( c, n, a, lay, wire ) ){
#if 0
      std::cout << " layer:" << lay << " wire:" << wire << std::endl;
#endif
      SetData( c, n, a, cid, lay, wire, data,    &CDCContainer[lay-1] );
    }
    else{
#if 0
      std::cout << " cannot find wires info ... " << std::endl;
#endif
    }
    break;
  }

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::SetData( const int &c, const int &n, const int &a,
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
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::SetData( const int &c, const int &n, const int &a,
			 const int &cid, const int &lay, const int &wire,
			 const int &data, CDCHitContainer *container )
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
#if 0
  std::cout << " in Set Data, layer:" << hit->layer() << " wire:" << hit->wire() << " pabsent:" << pabsent << std::endl;
#endif
  if( pabsent ){ AddHit( *hit ); delete hit; }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::SetSimulatedResolution( ConfMan *conf )
{
  for( int i=0; i<(int)CDHContainer.size(); i++ ){
    CDHContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)IHContainer.size(); i++ ){
    IHContainer[i].SetSimulatedResolution(conf);
  }
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    for( int i=0; i<(int)CDCContainer[layer-1].size(); i++ ){
      CDCContainer[layer-1][i].SetSimulatedResolution(conf);
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::Calc( ConfMan *conf )
{

  for( int i=0; i<(int)CDHContainer.size(); i++ ){
    CDHContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)IHContainer.size(); i++ ){
    IHContainer[i].Calc(conf);
  }

  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    for( int i=0; i<(int)CDCContainer[layer-1].size(); i++ ){
      CDCContainer[layer-1][i].Calc(conf);
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool CDSHitMan::RemoveNoTDCData()
{
  unsigned int key;

  // CDH
  for( HodoscopeLikeContainer::iterator ih=CDHContainer.begin();
       ih!=CDHContainer.end(); ){
    key = KEY( (*ih).cid(), (*ih).seg() );
    SegmentContainer::iterator is = HodoscopeSegContainer.find(key);
    if( is != HodoscopeSegContainer.end() ){ HodoscopeSegContainer.erase(is); }

    if( (*ih).CheckRange()==false ){
      ih = CDHContainer.erase(ih);
    }
    else{
      ih++;
    }
  }
  for( unsigned int i=0; i<CDHContainer.size(); i++ ){
    key = KEY( CDHContainer[i].cid(), CDHContainer[i].seg() );
    HodoscopeSegContainer[key] = i;
  }

  // IH
  for( HodoscopeLikeContainer::iterator ih=IHContainer.begin();
       ih!=IHContainer.end(); ){
    key = KEY( (*ih).cid(), (*ih).seg() );
    SegmentContainer::iterator is = HodoscopeSegContainer.find(key);
    if( is != HodoscopeSegContainer.end() ){ HodoscopeSegContainer.erase(is); }

    if( (*ih).CheckRange()==false ){
      ih = IHContainer.erase(ih);
    }
    else{
      ih++;
    }
  }
  for( unsigned int i=0; i<IHContainer.size(); i++ ){
    key = KEY( IHContainer[i].cid(), IHContainer[i].seg() );
    HodoscopeSegContainer[key] = i;
  }

  //

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int CDSHitMan::nHodo( const int &cid )
{
  switch( cid ){
  case CID_CDH:
    return CDHContainer.size();
  case CID_IH:
    return IHContainer.size();
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeHit * CDSHitMan::Hodo( const int &cid, const int &seg )
{
  unsigned int i = HodoscopeSegContainer[ KEY(cid,seg) ];
  switch( cid ){
  case CID_CDH:
    return (0<=i&&i<CDHContainer.size()) ? &CDHContainer[i] : 0;
  case CID_IH:
    return (0<=i&&i<IHContainer.size()) ? &IHContainer[i] : 0;
  }
  return 0;
}

HodoscopeLikeHit * CDSHitMan::Hodoi( const int &cid, const unsigned int &i )
{
  switch( cid ){
  case CID_CDH:
    return (0<=i&&i<CDHContainer.size()) ? &CDHContainer[i] : 0;
  case CID_IH:
    return (0<=i&&i<IHContainer.size()) ? &IHContainer[i] : 0;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void CDSHitMan::Clear()
{
  CDHContainer.clear(); 
  IHContainer.clear(); 
  for( int i=0; i<NumOfCDCLayers; i++ ){
    CDCContainer[i].clear();
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void CDSHitMan::CheckContainerSize()
{
  std::cout << " Size -> "
	    << " CDH:" << CDHContainer.size()
	    << std::endl;
  std::cout << " Size -> "
	    << " IH:" << IHContainer.size()
	    << std::endl;
  std::cout << "          "; std::cout << " CDC:";
  for( int i=0; i<NumOfCDCLayers;  i++ ){ std::cout << " , " << CDCContainer[i].size(); }
  std::cout << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int CDSHitMan::nCDC()
{
  int n=0;
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    n += nCDC(layer);
  }
  return n;
}
