#include "CDSHitMan.h"

ClassImp(CDSHitMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CDSHitMan::CDSHitMan() : HitMan()
{
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeContainer *CDSHitMan::HodoContainer( const int &cid )
{
  switch( cid ){
  case CID_IH:      return &IHContainer;
  case CID_CDH:     return &CDHContainer;    
  default:          return 0;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CDCLikeContainer *CDSHitMan::CDCHitContainer( const int &cid, const int &layer )
{
  switch( cid ){
  case CID_CDC:  
    return (0<layer&&layer<=NumOfCDCLayers) ? &CDCContainer[layer-1] : 0;
  default:    
    return 0;
  }
}
// -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
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
bool CDSHitMan::Calc( ConfMan *conf, double toffs )
{
  for( int i=0; i<(int)CDHContainer.size(); i++ )    CDHContainer[i].Calc(conf);
  for( int i=0; i<(int)IHContainer.size(); i++ )     IHContainer[i].Calc(conf);

  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    for( int i=0; i<(int)CDCContainer[layer-1].size(); i++ ){
      CDCContainer[layer-1][i].Calc(conf,toffs);
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void CDSHitMan::RemoveNoTDCData()
{
  HodoscopeSegContainer.clear();
  unsigned int key;
  const int CDSHodoIDList[]={CID_CDH,CID_IH};
  const int nCDSHodoIDList=sizeof(CDSHodoIDList)/sizeof(int);
  for(int i=0;i<nCDSHodoIDList;i++){
    HodoscopeLikeContainer::iterator ittr;
    HodoscopeLikeContainer *container=HodoContainer(CDSHodoIDList[i]);
    if(container){
      for( ittr=container->end()-1; ittr!=container->begin()-1; --ittr ){
	if(!ittr->CheckRange()) container->erase(ittr);
      }
    }
    for( unsigned int i=0; i<container->size(); i++ ){
      key = KEY( (*container)[i].cid(), (*container)[i].seg() );
      HodoscopeSegContainer[key] = i;
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void CDSHitMan::Clear()
{
  CDHContainer.clear(); 
  IHContainer.clear(); 
  for( int i=0; i<NumOfCDCLayers; i++ ){
    CDCContainer[i].clear();
  }
  HitMan::Clear();
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
