#include "SDDHitMan.h"

ClassImp(SDDHitMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
SDDHitMan::SDDHitMan() : TObject()
{
  Clear();
  CheckContainerSize();
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
void SDDHitMan::AddHit( const SDDHit &hit )
{
  int cid = hit.cid();
  switch( cid ){
  case CID_SDD:
    SDDContainer.push_back(hit);
    SDDContainer[ SDDContainer.size()-1 ].SetHitID( SDDContainer.size()-1 );
    //    SDDSegContainer[ KEY(hit.cid(),hit.seg()) ] = SDDContainer.size()-1;
    break;
  default:
    break;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void SDDHitMan::AddHit( const SDDHitCom &hit )
{
  int cid = hit.cid();
  switch( cid ){
  case CID_SDD:
    SDDCommon.push_back(hit);
    break;
  default:
    break;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
/*
void SDDHitMan::AddHit( const HodoscopeLikeHit &hit )
{
  int cid = hit.cid();
  switch( cid ){
  case CID_DEF:
    Define.push_back(hit);
    Define[ Define.size()-1 ].SetHitID( Define.size()-1 );
    break;
  default:
    break;
  }
}
*/
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::Convert( TKOHitCollection *tko, ConfMan *conf )
{
#if 0
  std::cout << " Enter SDDHitMan::Convert(TKO,conf) " << std::endl;
#endif
  //CheckContainerSize();

  for( int ih=0; ih<tko->entries(); ih++ ){
    int c, n, a, data;
    c = tko->hit(ih)->cr(); n = tko->hit(ih)->sl(); a = tko->hit(ih)->ch(); data = tko->hit(ih)->data();
#if 0
  std::cout << " c:" << c << " n:" << n << " a:" << a << " data:" << data << std::endl;
#endif
    Convert( c, n, a, data, conf );
  }
  //CheckContainerSize();

  Calc( conf );

  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::Convert( const int &c, const int &n, const int &a, const int &data, ConfMan *conf )
{
#if 0
  std::cout << " Enter SDDHitMan::Convert(c,n,a,data,conf) " << std::endl;
#endif
  int cid, lay, seg, at, ud;
  conf->GetCounterMapManager()->GetInfo( c, n, a, cid, lay, seg, at, ud );

#if 1
  std::cout << " c:" << c << " n:" << n << " a:" << a << " cid:" << cid << " lay:" << lay << " seg:" << seg
	    << " at:" << at << " ud:" << ud << " data:" << data << std::endl;
#endif

  //CheckContainerSize();
  
  switch( cid ){
  case CID_SDD:
    if(seg==0){
      SetData( c, n, a, cid, seg, at, ud, data, &SDDCommon );
    }else{
      SetData( c, n, a, cid, seg, at, ud, data, &SDDContainer );
    }
    break;
  case CID_DEF:
    if(at==0) ADef=data;
    else if(at==1) TDef=data;
    break;
  default:
    break;
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
/*
bool SDDHitMan::SetData( const int &c, const int &n, const int &a,
			 const int &cid, const int &seg, const int &at, const int &ud,
			 const int &data, DefineContainer *container )
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
*/
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::SetData( const int &c, const int &n, const int &a,
			 const int &cid, const int &seg, const int &at, const int &ud,
			 const int &data, SDDHitContainer *container )
{

  SDDHit *hit= 0;
  bool pabsent=false;
  
  for( int i=0; i<(int)container->size(); i++ ){
    if( (*container)[i].seg() == seg ){ hit = &(*container)[i];} 
  }
  
  if( hit==0 ){ 
    hit = new SDDHit(); pabsent=true;
#if 0
    std::cout<<"c: "<< c <<" n: "<< n <<" a: "<< a
	     <<" seg: "<<seg<<" data: "<<data<<std::endl;
#endif
  }
  hit->SetData( c, n, a, cid, seg, at, ud, data );
  if( pabsent ){ AddHit( *hit ); delete hit; }
  
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::SetData( EventStruct &vme )
{
  for(int seg=1;seg<=NumOfSDDs;seg++) SetData(seg,vme,&SDDContainer);   
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::SetData(const int &seg, EventStruct &vme, SDDHitContainer *container )
{
  SDDHit *hit= 0;
  bool pabsent=false;
  
  for( int i=0; i<(int)container->size(); i++ ){
    if( (*container)[i].seg() == seg ){ hit = &(*container)[i];} 
  }
  
  if( hit==0 ){ 
    hit = new SDDHit(); pabsent=true;
  }
  hit->SetData(seg, vme);
  if( pabsent ){ AddHit( *hit ); delete hit; }
  
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::SetData( const int &c, const int &n, const int &a,
			 const int &cid, const int &seg, const int &at, const int &ud,
			 const int &data, SDDComContainer *container )
{

  SDDHitCom *hit= 0;
  bool pabsent=false;
  for( int i=0; i<(int)container->size(); i++ ){
    if( (*container)[i].cid() == cid ){ hit = &(*container)[i];} 
  }  
  if( hit==0 ){ 
    hit = new SDDHitCom(); pabsent=true;
#if 0
    std::cout<<"c: "<< c <<" n: "<< n <<" a: "<< a
	     <<" seg: "<<seg<<" data: "<<data<<std::endl;
#endif
  }
  hit->SetData( c, n, a, cid, seg, at, ud, data );
  if( pabsent ){ AddHit( *hit ) ; delete hit; }
  
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool SDDHitMan::Calc( ConfMan *conf )
{
  for( int i=0; i<(int)SDDContainer.size(); i++ ){
    SDDContainer[i].Calc(conf);
  }

  //  Define.Calc(conf);
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
//
//  HodoscopeLikeHit * SDDHitMan::SDD( const int &cid, const int &seg )
//  {
//  unsigned int i = HodoscopeSegContainer[ KEY(cid,seg) ];
//  switch( cid ){
//  case CID_SDD:
//    return (0<=i&&i<SDDContainer.size()) ? &SDDContainer[i] : 0;
//  }
//  return 0;
//}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int SDDHitMan::nHit(int threshold)
{
  int nhit=0;
  for(int i=0;i<(int)SDDContainer.size();i++)
    if(SDDContainer[i].tdcsdd()>threshold)
      if(SDDContainer[i].tdcsdd()<threshold+800) nhit++;
 
  return nhit;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int SDDHitMan::reset(int nth)
{
  int nn=0;
  int* temp=0;
  for(int i=0;i<(int)SDDContainer.size();i++){
    if(SDDContainer[i].treset()>0){
      temp[nn]=i+1;
      nn++;
    }
  }
  if(nn==0) return 0;
  else return temp[nth];
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
int SDDHitMan::highth(int nth)
{
  int nn=0;
  int* temp=0;
  for(int i=0;i<(int)SDDContainer.size();i++){
    if(SDDContainer[i].thigh()>0){
      temp[nn]=i+1;
      nn++;
    }
  }
  if(nn==0) return 0;
  else return temp[nth];
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void SDDHitMan::Clear()
{
  SDDContainer.clear(); 
  SDDCommon.clear();
  ADef=TDef=-1;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void SDDHitMan::CheckContainerSize()
{
  std::cout << " Size -> "
	    << " SDD:" << SDDContainer.size()
	    << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
