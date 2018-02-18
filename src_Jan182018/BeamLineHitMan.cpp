#include "BeamLineHitMan.h"

ClassImp(BeamLineHitMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
BeamLineHitMan::BeamLineHitMan() : TObject()
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
void BeamLineHitMan::AddHit( const HodoscopeLikeHit &hit )
{
  int cid = hit.cid();
  switch( cid ){
  case CID_BHD:
    BHDContainer.push_back(hit);
    BHDContainer[ BHDContainer.size()-1 ].SetHitID( BHDContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = BHDContainer.size()-1;
    break;
  case CID_BHDpost:
    BHDpostContainer.push_back(hit);
    BHDpostContainer[ BHDpostContainer.size()-1 ].SetHitID( BHDpostContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = BHDpostContainer.size()-1;
    break;
  case CID_PA:
    PAContainer.push_back(hit);
    PAContainer[ PAContainer.size()-1 ].SetHitID( PAContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = PAContainer.size()-1;
    break;
  case CID_T0:
    T0Container.push_back(hit);
    T0Container[ T0Container.size()-1 ].SetHitID( T0Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = T0Container.size()-1;
    break;
  case CID_T0pre:
    T0preContainer.push_back(hit);
    T0preContainer[ T0preContainer.size()-1 ].SetHitID( T0preContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = T0preContainer.size()-1;
    break;
  case CID_T0post:
    T0postContainer.push_back(hit);
    T0postContainer[ T0postContainer.size()-1 ].SetHitID( T0postContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = T0postContainer.size()-1;
    break;
  case CID_E0:
    E0Container.push_back(hit);
    E0Container[ E0Container.size()-1 ].SetHitID( E0Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = E0Container.size()-1;
    break;
  case CID_B1:
    B1Container.push_back(hit);
    B1Container[ B1Container.size()-1 ].SetHitID( B1Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = B1Container.size()-1;
    break;
  case CID_B2:
    B2Container.push_back(hit);
    B2Container[ B2Container.size()-1 ].SetHitID( B2Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = B2Container.size()-1;
    break;
  case CID_Range:
    RangeContainer.push_back(hit);
    RangeContainer[ RangeContainer.size()-1 ].SetHitID( RangeContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = RangeContainer.size()-1;
    break;
  case CID_CVC:
    CVCContainer.push_back(hit);
    CVCContainer[ CVCContainer.size()-1 ].SetHitID( CVCContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = CVCContainer.size()-1;
    break;
  case CID_NC:
    NCContainer.push_back(hit);
    NCContainer[ NCContainer.size()-1 ].SetHitID( NCContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = NCContainer.size()-1;
    break;
  case CID_PC:
    PCContainer.push_back(hit);
    PCContainer[ PCContainer.size()-1 ].SetHitID( PCContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = PCContainer.size()-1;
    break;
  case CID_BPD:
    BPDContainer.push_back(hit);
    BPDContainer[ BPDContainer.size()-1 ].SetHitID( BPDContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = BPDContainer.size()-1;
    break;
  case CID_BVC:
    BVContainer.push_back(hit);
    BVContainer[ BVContainer.size()-1 ].SetHitID( BVContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = BVContainer.size()-1;
    break;
  case CID_WVC:
    WVContainer.push_back(hit);
    WVContainer[ WVContainer.size()-1 ].SetHitID( WVContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = WVContainer.size()-1;
    break;
  case CID_HVC1:
    HVC1Container.push_back(hit);
    HVC1Container[ HVC1Container.size()-1 ].SetHitID( HVC1Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = HVC1Container.size()-1;
    break;
  case CID_HVC2:
    HVC2Container.push_back(hit);
    HVC2Container[ HVC2Container.size()-1 ].SetHitID( HVC2Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = HVC2Container.size()-1;
    break;
  case CID_HVC3:
    HVC3Container.push_back(hit);
    HVC3Container[ HVC3Container.size()-1 ].SetHitID( HVC3Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = HVC3Container.size()-1;
    break;
  case CID_LB:
    LBContainer.push_back(hit);
    LBContainer[ LBContainer.size()-1 ].SetHitID( LBContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = LBContainer.size()-1;
    break;
  case CID_BD:
    BDContainer.push_back(hit);
    BDContainer[ BDContainer.size()-1 ].SetHitID( BDContainer.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = BDContainer.size()-1;
    break;
  case CID_TEMP1:
    Temp1Container.push_back(hit);
    Temp1Container[ Temp1Container.size()-1 ].SetHitID( Temp1Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = Temp1Container.size()-1;
  case CID_TEMP2:
    Temp2Container.push_back(hit);
    Temp2Container[ Temp2Container.size()-1 ].SetHitID( Temp2Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = Temp2Container.size()-1;
    break;
  case CID_TEMP3:
    Temp3Container.push_back(hit);
    Temp3Container[ Temp3Container.size()-1 ].SetHitID( Temp3Container.size()-1 );
    HodoscopeSegContainer[ KEY(hit.cid(),hit.seg()) ] = Temp3Container.size()-1;
    break;
    break;
  default:
    break;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::AddHit( const CherenkovLikeHit &hit )
{
#if 0
  std::cout << " Enter Add Hit (Cherenkov) " << std::endl;
#endif
  int cid = hit.cid();
  switch( cid ){
  case CID_LC1:
    LC1Container.push_back(hit);
    LC1Container[ LC1Container.size()-1 ].SetHitID( LC1Container.size()-1 );
    CherenkovSegContainer[ KEY(hit.cid(),hit.seg()) ] = LC1Container.size()-1;
    break;
  case CID_LC2:
    LC2Container.push_back(hit);
    LC2Container[ LC2Container.size()-1 ].SetHitID( LC2Container.size()-1 );
    CherenkovSegContainer[ KEY(hit.cid(),hit.seg()) ] = LC2Container.size()-1;
    break;
  case CID_AC:
    ACContainer.push_back(hit);
    ACContainer[ ACContainer.size()-1 ].SetHitID( ACContainer.size()-1 );
    CherenkovSegContainer[ KEY(hit.cid(),hit.seg()) ] = ACContainer.size()-1;
    break;
  case CID_WC:
    WCContainer.push_back(hit);
    WCContainer[ WCContainer.size()-1 ].SetHitID( WCContainer.size()-1 );
    CherenkovSegContainer[ KEY(hit.cid(),hit.seg()) ] = WCContainer.size()-1;
    break;
  case CID_GC:
    GCContainer.push_back(hit);
    GCContainer[ GCContainer.size()-1 ].SetHitID( GCContainer.size()-1 );
    CherenkovSegContainer[ KEY(hit.cid(),hit.seg()) ] = GCContainer.size()-1;
    break;
  default:
    break;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::AddHit( const ChamberLikeHit &hit )
{
#if 0
  std::cout << " Enter Add Hit (Chamber) " << std::endl;
#endif
  int cid = hit.cid();
  int layer = hit.layer();
  switch( cid ){
  case CID_BLC1a:
    BLC1aContainer[layer-1].push_back(hit);
    BLC1aContainer[layer-1][ BLC1aContainer[layer-1].size()-1 ].SetHitID( BLC1aContainer[layer-1].size()-1 );
    break;
  case CID_BLC1b:
    BLC1bContainer[layer-1].push_back(hit);
    BLC1bContainer[layer-1][ BLC1bContainer[layer-1].size()-1 ].SetHitID( BLC1bContainer[layer-1].size()-1 );
    break;
  case CID_BLC2a:
    BLC2aContainer[layer-1].push_back(hit);
    BLC2aContainer[layer-1][ BLC2aContainer[layer-1].size()-1 ].SetHitID( BLC2aContainer[layer-1].size()-1 );
    break;
  case CID_BLC2b:
    BLC2bContainer[layer-1].push_back(hit);
    BLC2bContainer[layer-1][ BLC2bContainer[layer-1].size()-1 ].SetHitID( BLC2bContainer[layer-1].size()-1 );
    break;
  case CID_BPC:
    BPCContainer[layer-1].push_back(hit);
    BPCContainer[layer-1][ BPCContainer[layer-1].size()-1 ].SetHitID( BPCContainer[layer-1].size()-1 );
    break;
  case CID_FDC1:
    FDC1Container[layer-1].push_back(hit);
    FDC1Container[layer-1][ FDC1Container[layer-1].size()-1 ].SetHitID( FDC1Container[layer-1].size()-1 );
    break;
  default:
    break;
  }
  //CheckContainerSize();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::Convert( TKOHitCollection *tko, ConfMan *conf )
{
#if 0
  std::cout << " Enter BeamLineHitMan::Convert(TKO,conf) " << std::endl;
#endif
  //CheckContainerSize();
  for( int ih=0; ih<tko->entries(); ih++ ){
    int c, n, a, data;
    c = tko->hit(ih)->cr(); n = tko->hit(ih)->sl(); a = tko->hit(ih)->ch(); data = tko->hit(ih)->data();
    Convert( c, n, a, data, conf );
  }
  //CheckContainerSize();
  //  std::cout << " c n a" <<c<<" "<<n<<" "<<a<< std::endl;

  Calc( conf );

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf )
{
#if 0
  std::cout << " Enter BeamLineHitMan::Convert(c,n,a,data,conf) " << std::endl;
#endif
  int cid, lay, seg, at, ud;
  conf->GetCounterMapManager()->GetInfo( c, n, a, cid, lay, seg, at, ud );

#if 0
  if(c==10)
    std::cout << " c:" << c << " n:" << n << " a:" << a << " cid:" << cid << " lay:" << lay << " seg:" << seg
	      << " at:" << at << " ud:" << ud << " data:" << data << std::endl;
#endif

  //Checkcontainersize();

  switch( cid ){
  case CID_BHD:
    SetData( c, n, a, cid, seg, at, ud, data, &BHDContainer );
    break;
  case CID_BHDpost:
    SetData( c, n, a, cid, seg, at, ud, data, &BHDpostContainer );
    break;
  case CID_PA:
    SetData( c, n, a, cid, seg, at, ud, data, &PAContainer );
    break;
  case CID_T0:
    SetData( c, n, a, cid, seg, at, ud, data, &T0Container );
    break;
  case CID_T0pre:
    SetData( c, n, a, cid, seg, at, ud, data, &T0preContainer );
    break;
  case CID_T0post:
    SetData( c, n, a, cid, seg, at, ud, data, &T0postContainer );
    break;
  case CID_E0:
    SetData( c, n, a, cid, seg, at, ud, data, &E0Container );
    break;
  case CID_B1:
    SetData( c, n, a, cid, seg, at, ud, data, &B1Container );
    break;
  case CID_B2:
    SetData( c, n, a, cid, seg, at, ud, data, &B2Container );
    break;
  case CID_Range:
    SetData( c, n, a, cid, seg, at, ud, data, &RangeContainer );
    break;
  case CID_CVC:
    SetData( c, n, a, cid, seg, at, ud, data, &CVCContainer );
    break;
  case CID_NC:
    SetData( c, n, a, cid, seg, at, ud, data, &NCContainer );
    break;
  case CID_PC:
    SetData( c, n, a, cid, seg, at, ud, data, &PCContainer );
    break;
  case CID_BPD:
    SetData( c, n, a, cid, seg, at, ud, data, &BPDContainer );
    break;
  case CID_BVC:
    SetData( c, n, a, cid, seg, at, ud, data, &BVContainer );
    break;
  case CID_WVC:
    SetData( c, n, a, cid, seg, at, ud, data, &WVContainer );
    break;
  case CID_HVC1:
    SetData( c, n, a, cid, seg, at, ud, data, &HVC1Container );
    break;
  case CID_HVC2:
    SetData( c, n, a, cid, seg, at, ud, data, &HVC2Container );
    break;
  case CID_HVC3:
    SetData( c, n, a, cid, seg, at, ud, data, &HVC3Container );
    break;
  case CID_LB:
    SetData( c, n, a, cid, seg, at, ud, data, &LBContainer );
    break;
  case CID_BD:
    SetData( c, n, a, cid, seg, at, ud, data, &BDContainer );
    break;
  case CID_TEMP1:
    SetData( c, n, a, cid, seg, at, ud, data, &Temp1Container );
    break;
  case CID_TEMP2:
    SetData( c, n, a, cid, seg, at, ud, data, &Temp2Container );
    break;
  case CID_TEMP3:
    SetData( c, n, a, cid, seg, at, ud, data, &Temp3Container );
    break;
  case CID_LC1:
    SetData( c, n, a, cid, seg, at, ud, data, &LC1Container );
    break;
  case CID_LC2:
    SetData( c, n, a, cid, seg, at, ud, data, &LC2Container );
    break;
  case CID_AC:
    SetData( c, n, a, cid, seg, at, ud, data, &ACContainer );
    break;
  case CID_WC:
    SetData( c, n, a, cid, seg, at, ud, data, &WCContainer );
    break;
  case CID_GC:
    SetData( c, n, a, cid, seg, at, ud, data, &GCContainer );
    break;
  case CID_BLC1a:
    SetData( c, n, a, cid, lay, seg, data,    &(BLC1aContainer[lay-1]) );
    break;
  case CID_BLC1b:
    SetData( c, n, a, cid, lay, seg, data,    &(BLC1bContainer[lay-1]) );
    break;
  case CID_BLC2a:
    SetData( c, n, a, cid, lay, seg, data,    &(BLC2aContainer[lay-1]) );
    break;
  case CID_BLC2b:
    SetData( c, n, a, cid, lay, seg, data,    &(BLC2bContainer[lay-1]) );
    break;
  case CID_BPC:
    SetData( c, n, a, cid, lay, seg, data,    &(BPCContainer[lay-1]) );
    break;
  case CID_FDC1:
    SetData( c, n, a, cid, lay, seg, data,    &(FDC1Container[lay-1]) );
    break;
  default:
    break;
  }

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::SetData( const int &c, const int &n, const int &a,
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
bool BeamLineHitMan::SetData( const int &c, const int &n, const int &a,
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
bool BeamLineHitMan::SetData( const int &c, const int &n, const int &a,
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
bool BeamLineHitMan::SetSimulatedResolution( ConfMan *conf )
{
  for( int i=0; i<(int)BHDContainer.size(); i++ ){
    BHDContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)BHDpostContainer.size(); i++ ){
    BHDpostContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)PAContainer.size(); i++ ){
    PAContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)T0Container.size(); i++ ){
    T0Container[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)E0Container.size(); i++ ){
    E0Container[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)B1Container.size(); i++ ){
    B1Container[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)B2Container.size(); i++ ){
    B2Container[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)RangeContainer.size(); i++ ){
    RangeContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)CVCContainer.size(); i++ ){
    CVCContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)NCContainer.size(); i++ ){
    NCContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)PCContainer.size(); i++ ){
    PCContainer[i].SetSimulatedResolution(conf);
  }
  for( int i=0; i<(int)BPDContainer.size(); i++ ){
    BPDContainer[i].SetSimulatedResolution(conf);
  }

  for( int i=0; i<(int)LC1Container.size(); i++ ){
  }
  for( int i=0; i<(int)LC2Container.size(); i++ ){
  }
  for( int i=0; i<(int)ACContainer.size(); i++ ){
  }
  for( int i=0; i<(int)WCContainer.size(); i++ ){
  }
  for( int i=0; i<(int)GCContainer.size(); i++ ){
  }

 for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC1aContainer[lay-1].size(); i++ ){
      BLC1aContainer[lay-1][i].SetSimulatedResolution(conf);
    }
  }
 for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC1bContainer[lay-1].size(); i++ ){
      BLC1bContainer[lay-1][i].SetSimulatedResolution(conf);
    }
  }
  for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC2aContainer[lay-1].size(); i++ ){
      BLC2aContainer[lay-1][i].SetSimulatedResolution(conf);
    }
  }
  for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC2bContainer[lay-1].size(); i++ ){
      BLC2bContainer[lay-1][i].SetSimulatedResolution(conf);
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::Calc( ConfMan *conf )
{

  for( int i=0; i<(int)BHDContainer.size(); i++ ){
    BHDContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)BHDpostContainer.size(); i++ ){
    BHDpostContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)PAContainer.size(); i++ ){
    PAContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)T0Container.size(); i++ ){
    T0Container[i].Calc(conf);
  }
  for( int i=0; i<(int)T0preContainer.size(); i++ ){
    T0preContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)T0postContainer.size(); i++ ){
    T0postContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)E0Container.size(); i++ ){
    E0Container[i].Calc(conf);    
  }
  for( int i=0; i<(int)B1Container.size(); i++ ){
    B1Container[i].Calc(conf);
  }
  for( int i=0; i<(int)B2Container.size(); i++ ){
    B2Container[i].Calc(conf);
  }
  for( int i=0; i<(int)RangeContainer.size(); i++ ){
    RangeContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)CVCContainer.size(); i++ ){
    CVCContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)NCContainer.size(); i++ ){
    NCContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)PCContainer.size(); i++ ){
    PCContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)BPDContainer.size(); i++ ){
    BPDContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)BVContainer.size(); i++ ){
    BVContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)WVContainer.size(); i++ ){
    WVContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)HVC1Container.size(); i++ ){
    HVC1Container[i].Calc(conf);
  }
  for( int i=0; i<(int)HVC2Container.size(); i++ ){
    HVC2Container[i].Calc(conf);
  }
  for( int i=0; i<(int)HVC3Container.size(); i++ ){
    HVC3Container[i].Calc(conf);
  }
  for( int i=0; i<(int)LBContainer.size(); i++ ){
    LBContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)BDContainer.size(); i++ ){
    BDContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)Temp1Container.size(); i++ ){
    Temp1Container[i].Calc(conf);
  }
  for( int i=0; i<(int)Temp2Container.size(); i++ ){
    Temp2Container[i].Calc(conf);
  }
  for( int i=0; i<(int)Temp3Container.size(); i++ ){
    Temp3Container[i].Calc(conf);
  }
  for( int i=0; i<(int)LC1Container.size(); i++ ){
    LC1Container[i].Calc(conf);
  }
  for( int i=0; i<(int)LC2Container.size(); i++ ){
    LC2Container[i].Calc(conf);
  }
  for( int i=0; i<(int)ACContainer.size(); i++ ){
    ACContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)WCContainer.size(); i++ ){
    WCContainer[i].Calc(conf);
  }
  for( int i=0; i<(int)GCContainer.size(); i++ ){
    GCContainer[i].Calc(conf);
  }

  CalcT0Time();
  CalcBLDC(conf,T0Time);
  return true;
}

void BeamLineHitMan::CalcBLDC(ConfMan *conf, const double &toffs){
  for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC1aContainer[lay-1].size(); i++ ){
      BLC1aContainer[lay-1][i].Calc(conf, toffs);
    }
  }
  for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC1bContainer[lay-1].size(); i++ ){
      BLC1bContainer[lay-1][i].Calc(conf, toffs);
    }
  }
  for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC2aContainer[lay-1].size(); i++ ){
      BLC2aContainer[lay-1][i].Calc(conf ,toffs);
    }
  }
  for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
    for( int i=0; i<(int)BLC2bContainer[lay-1].size(); i++ ){
      BLC2bContainer[lay-1][i].Calc(conf, toffs);
    }
  }
  for( int lay=1; lay<=NumOfBPCLayers; lay++ ){
    for( int i=0; i<(int)BPCContainer[lay-1].size(); i++ ){
      BPCContainer[lay-1][i].Calc(conf, toffs);
    }
  }
  for( int lay=1; lay<=NumOfFDC1Layers; lay++ ){
    for( int i=0; i<(int)FDC1Container[lay-1].size(); i++ ){
      FDC1Container[lay-1][i].Calc(conf, toffs);
    }
  }
  //  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeHit * BeamLineHitMan::Hodo( const int &cid, const int &seg )
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
#if 0
  std::cout << " in Hodo, cid:" << cid << " seg:" << seg << " i:" << i 
	    << " KEY:" << std::hex << KEY(cid,seg) << std::dec
	    << std::endl;
#endif
  switch( cid ){
  case CID_BHD:
    return (0<=i&&i<BHDContainer.size()) ? &BHDContainer[i] : 0;
  case CID_BHDpost:
    return (0<=i&&i<BHDpostContainer.size()) ? &BHDpostContainer[i] : 0;
  case CID_PA:
    return (0<=i&&i<PAContainer.size())  ? &PAContainer[i] : 0;
  case CID_T0:
    return (0<=i&&i<T0Container.size())  ? &T0Container[i] : 0;
  case CID_T0pre:
    return (0<=i&&i<T0preContainer.size())  ? &T0preContainer[i] : 0;
  case CID_T0post:
    return (0<=i&&i<T0postContainer.size())  ? &T0postContainer[i] : 0;
  case CID_E0:
    return (0<=i&&i<E0Container.size())  ? &E0Container[i] : 0;
  case CID_B1:
    return (0<=i&&i<B1Container.size())  ? &B1Container[i] : 0;
  case CID_B2:
    return (0<=i&&i<B2Container.size())  ? &B2Container[i] : 0;
  case CID_Range:
    return (0<=i&&i<RangeContainer.size()) ? &RangeContainer[i] : 0;
  case CID_CVC:
    return (0<=i&&i<CVCContainer.size()) ? &CVCContainer[i] : 0;
  case CID_NC:
    return (0<=i&&i<NCContainer.size()) ? &NCContainer[i] : 0;
  case CID_PC:
    return (0<=i&&i<PCContainer.size()) ? &PCContainer[i] : 0;
  case CID_BPD:
    return (0<=i&&i<BPDContainer.size()) ? &BPDContainer[i] : 0;
  case CID_BVC:
    return (0<=i&&i<BVContainer.size()) ? &BVContainer[i] : 0;
  case CID_WVC:
    return (0<=i&&i<BVContainer.size()) ? &WVContainer[i] : 0;
  case CID_HVC1:
    return (0<=i&&i<HVC1Container.size()) ? &HVC1Container[i] : 0;
  case CID_HVC2:
    return (0<=i&&i<HVC2Container.size()) ? &HVC2Container[i] : 0;
  case CID_HVC3:
    return (0<=i&&i<HVC3Container.size()) ? &HVC3Container[i] : 0;
  case CID_LB:
    return (0<=i&&i<LBContainer.size()) ? &LBContainer[i] : 0;
  case CID_BD:
    return (0<=i&&i<BDContainer.size()) ? &BDContainer[i] : 0;
  case CID_TEMP1:
    return (0<=i&&i<Temp1Container.size()) ? &Temp1Container[i] : 0;
  case CID_TEMP2:
    return (0<=i&&i<Temp2Container.size()) ? &Temp2Container[i] : 0;
  case CID_TEMP3:
    return (0<=i&&i<Temp3Container.size()) ? &Temp3Container[i] : 0;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nHodo( const int &cid )
{
  switch( cid ){
  case CID_BHD:
    return BHDContainer.size();
  case CID_BHDpost:
    return BHDpostContainer.size();
  case CID_PA:
    return PAContainer.size();
  case CID_T0:
    return T0Container.size();
  case CID_T0pre:
    return T0preContainer.size();
  case CID_T0post:
    return T0postContainer.size();
  case CID_E0:
    return E0Container.size();
  case CID_B1:
    return B1Container.size();
  case CID_B2:
    return B2Container.size();
  case CID_Range:
    return RangeContainer.size();
  case CID_CVC:
    return CVCContainer.size();
  case CID_NC:
    return NCContainer.size();
  case CID_PC:
    return PCContainer.size();
  case CID_BPD:
    return BPDContainer.size();
  case CID_BVC:
    return BVContainer.size();
  case CID_WVC:
    return WVContainer.size();
  case CID_HVC1:
    return HVC1Container.size();
  case CID_HVC2:
    return HVC2Container.size();
  case CID_HVC3:
    return HVC3Container.size();
  case CID_LB:
    return LBContainer.size();
  case CID_BD:
    return BDContainer.size();
  case CID_TEMP1:
    return Temp1Container.size();
  case CID_TEMP2:
    return Temp2Container.size();
  case CID_TEMP3:
    return Temp3Container.size();
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeHit * BeamLineHitMan::Hodoi( const int &cid, const unsigned int &i )
{
  switch( cid ){
  case CID_BHD:
    return (0<=i&&i<BHDContainer.size()) ? &BHDContainer[i] : 0;
  case CID_BHDpost:
    return (0<=i&&i<BHDpostContainer.size()) ? &BHDpostContainer[i] : 0;
  case CID_PA:
    return (0<=i&&i<PAContainer.size())  ? &PAContainer[i] : 0;
  case CID_T0:
    return (0<=i&&i<T0Container.size())  ? &T0Container[i] : 0;
  case CID_T0pre:
    return (0<=i&&i<T0preContainer.size())  ? &T0preContainer[i] : 0;
  case CID_T0post:
    return (0<=i&&i<T0postContainer.size())  ? &T0postContainer[i] : 0;
  case CID_E0:
    return (0<=i&&i<E0Container.size())  ? &E0Container[i] : 0;
  case CID_B1:
    return (0<=i&&i<B1Container.size())  ? &B1Container[i] : 0;
  case CID_B2:
    return (0<=i&&i<B2Container.size())  ? &B2Container[i] : 0;
  case CID_Range:
    return (0<=i&&i<RangeContainer.size()) ? &RangeContainer[i] : 0;
  case CID_CVC:
    return (0<=i&&i<CVCContainer.size()) ? &CVCContainer[i] : 0;
  case CID_NC:
    return (0<=i&&i<NCContainer.size())  ? &NCContainer[i] : 0;
  case CID_PC:
    return (0<=i&&i<PCContainer.size())  ? &PCContainer[i] : 0;
  case CID_BPD:
    return (0<=i&&i<BPDContainer.size())  ? &BPDContainer[i] : 0;
  case CID_BVC:
    return (0<=i&&i<BVContainer.size())  ? &BVContainer[i] : 0;
  case CID_WVC:
    return (0<=i&&i<WVContainer.size())  ? &WVContainer[i] : 0;
  case CID_HVC1:
    return (0<=i&&i<HVC1Container.size())  ? &HVC1Container[i] : 0;
  case CID_HVC2:
    return (0<=i&&i<HVC2Container.size())  ? &HVC2Container[i] : 0;
  case CID_HVC3:
    return (0<=i&&i<HVC3Container.size())  ? &HVC3Container[i] : 0;
  case CID_LB:
    return (0<=i&&i<LBContainer.size())  ? &LBContainer[i] : 0;
  case CID_BD:
    return (0<=i&&i<BDContainer.size())  ? &BDContainer[i] : 0;
  case CID_TEMP1:
    return (0<=i&&i<Temp1Container.size())  ? &Temp1Container[i] : 0;
  case CID_TEMP2:
    return (0<=i&&i<Temp2Container.size())  ? &Temp2Container[i] : 0;
  case CID_TEMP3:
    return (0<=i&&i<Temp3Container.size())  ? &Temp3Container[i] : 0;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CherenkovLikeHit * BeamLineHitMan::Chere( const int &cid, const int &seg )
{
  unsigned int i = CherenkovSegContainer[ KEY(cid,seg) ];
  switch( cid ){
  case CID_LC1:
    return (0<=i&&i<LC1Container.size()) ? &LC1Container[i] : 0;
  case CID_LC2:
    return (0<=i&&i<LC2Container.size()) ? &LC2Container[i] : 0;
  case CID_AC:
    return (0<=i&&i<ACContainer.size())  ? &ACContainer[i] : 0;
  case CID_WC:
    return (0<=i&&i<WCContainer.size())  ? &WCContainer[i] : 0;
  case CID_GC:
    return (0<=i&&i<GCContainer.size())  ? &GCContainer[i] : 0;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nChere( const int &cid )
{
  switch( cid ){
  case CID_LC1:
    return LC1Container.size();    
  case CID_LC2:
    return LC2Container.size();
  case CID_AC:
    return ACContainer.size();
  case CID_WC:
    return WCContainer.size();
  case CID_GC:
    return GCContainer.size();
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CherenkovLikeHit * BeamLineHitMan::Cherei( const int &cid, const unsigned int &i )
{
  switch( cid ){
  case CID_LC1:
    return (0<=i&&i<LC1Container.size()) ? &LC1Container[i] : 0;
  case CID_LC2:
    return (0<=i&&i<LC2Container.size()) ? &LC2Container[i] : 0;
  case CID_AC:
    return (0<=i&&i<ACContainer.size())  ? &ACContainer[i] : 0;
  case CID_WC:
    return (0<=i&&i<WCContainer.size())  ? &WCContainer[i] : 0;
  case CID_GC:
    return (0<=i&&i<GCContainer.size())  ? &GCContainer[i] : 0;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::FindTrackPDC()
{
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::Clear()
{
  BHDContainer.clear(); 
  BHDpostContainer.clear(); 
  PAContainer.clear();  
  T0Container.clear();  
  T0preContainer.clear();  
  T0postContainer.clear();  
  E0Container.clear();  
  B1Container.clear();  
  B2Container.clear();  
  RangeContainer.clear();
  CVCContainer.clear(); 
  NCContainer.clear();  
  PCContainer.clear();  
  BPDContainer.clear();  
  BVContainer.clear();  
  WVContainer.clear();  
  HVC1Container.clear();  
  HVC2Container.clear();  
  HVC3Container.clear();  
  LBContainer.clear();  
  BDContainer.clear();  
  Temp1Container.clear();  
  Temp2Container.clear();  
  Temp3Container.clear();  
  LC1Container.clear(); 
  LC2Container.clear(); 
  ACContainer.clear();  
  WCContainer.clear();  
  GCContainer.clear();  
  for( int i=0; i<NumOfBLCLayers; i++ ){
    BLC1aContainer[i].clear();
    BLC1bContainer[i].clear();
    BLC2aContainer[i].clear();
    BLC2bContainer[i].clear();
  }
  for( int i=0; i<NumOfBPCLayers; i++ ){
    BPCContainer[i].clear();
  }
  for( int i=0; i<NumOfFDC1Layers; i++ ){
    FDC1Container[i].clear();
  }
  T0Time=0.;
  T0Seg=-1;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::CheckContainerSize()
{
  std::cout << " Size -> "
	    << " BHD:" << BHDContainer.size()
	    << " PA:"  << PAContainer.size()
	    << " T0:"  << T0Container.size()
	    << " E0:"  << E0Container.size()
	    << " B1:"  << B1Container.size()
	    << " B2:"  << B2Container.size()
	    << " Range:"  << RangeContainer.size()
	    << " CVC:"  << CVCContainer.size()
	    << " NC:"  << NCContainer.size()
	    << " PC:"  << PCContainer.size()
	    << " BPD:"  << BPDContainer.size()
    //	    << " CV:"  << CVContainer.size()
	    << " BVC:"  << BVContainer.size()
	    << " Longbar:"  << LBContainer.size()
	    << " Beamdump:"  << BDContainer.size()
	    << " Temp1:"  << Temp1Container.size()
	    << std::endl
	    << "          "
	    << " LC1:"  << LC1Container.size()
	    << " LC2:"  << LC2Container.size()
	    << " AC:"   << ACContainer.size()
	    << " WC:"   << WCContainer.size()
	    << " GC:"   << GCContainer.size()
	    << std::endl;
  std::cout << "          "; std::cout << " BLC1a:";
  for( int i=0; i<NumOfBLCLayers;  i++ ){ std::cout << " , " << BLC1aContainer[i].size(); }
  std::cout << std::endl;
  std::cout << "          "; std::cout << " BLC1b:";
  for( int i=0; i<NumOfBLCLayers;  i++ ){ std::cout << " , " << BLC1bContainer[i].size(); }
  std::cout << std::endl;
  std::cout << "          "; std::cout << " BLC2a:";
  for( int i=0; i<NumOfBLCLayers;  i++ ){ std::cout << " , " << BLC2aContainer[i].size(); }
  std::cout << std::endl;
  std::cout << "          "; std::cout << " BLC2b:";
  for( int i=0; i<NumOfBLCLayers;  i++ ){ std::cout << " , " << BLC2bContainer[i].size(); }
  std::cout << std::endl;
  std::cout << "          "; std::cout << " BPC:";
  for( int i=0; i<NumOfBPCLayers;  i++ ){ std::cout << " , " << BPCContainer[i].size(); }
  std::cout << std::endl;
  std::cout << "          "; std::cout << " FDC1:";
  for( int i=0; i<NumOfFDC1Layers;  i++ ){ std::cout << " , " << FDC1Container[i].size(); }
  std::cout << std::endl;
}

void BeamLineHitMan::CalcT0Time(){
  for(int i=0;i<nT0();i++)
    if(T0(i)->CheckRange2())
      if(T0(i)->ctmean()<T0Time){
	T0Time=T0(i)->ctmean();
	T0Seg=T0(i)->seg();
      }
  if(T0Time>100)
    for(int i=0;i<nT0();i++)
      if(T0(i)->CheckRange())
	if(T0(i)->ctmean()<T0Time){
	  T0Time=T0(i)->ctmean();      
	  T0Seg=T0(i)->seg();
	}
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nBLC1a()
{
  int n=0;
  for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
    n += nBLC1a(layer);
  }
  return n;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nBLC1b()
{
  int n=0;
  for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
    n += nBLC1b(layer);
  }
  return n;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nBLC2a()
{
  int n=0;
  for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
    n += nBLC2a(layer);
  }
  return n;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nBLC2b()
{
  int n=0;
  for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
    n += nBLC2b(layer);
  }
  return n;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nBPC()
{
  int n=0;
  for( int layer=1; layer<=NumOfBPCLayers; layer++ ){
    n += nBPC(layer);
  }
  return n;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nFDC1()
{
  int n=0;
  for( int layer=1; layer<=NumOfFDC1Layers; layer++ ){
    n += nFDC1(layer);
  }
  return n;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int BeamLineHitMan::nBLDC( const int &cid, const int &layer){
  if(cid==CID_BPC)
    return nBPC(layer);
  if(cid==CID_FDC1)
    return nFDC1(layer);
  if(cid==CID_BLC1a)
    return nBLC1a(layer);
  if(cid==CID_BLC1b)
    return nBLC1b(layer);
  if(cid==CID_BLC2a)
    return nBLC2a(layer);
  if(cid==CID_BLC2b)
    return nBLC2b(layer);
  return -1;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
ChamberLikeHit* BeamLineHitMan::BLDC( const int &cid, const int &layer, const int &i){
  if(cid==CID_BPC)
    return BPC(layer,i);
  if(cid==CID_FDC1)
    return FDC1(layer,i);
  if(cid==CID_BLC1a)
    return BLC1a(layer,i);
  if(cid==CID_BLC1b)
    return BLC1b(layer,i);
  if(cid==CID_BLC2a)
    return BLC2a(layer,i);
  if(cid==CID_BLC2b)
    return BLC2b(layer,i);
  return 0;
}
