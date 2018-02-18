#include "BeamLineHitMan.h"


ClassImp(BeamLineHitMan);
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
BeamLineHitMan::BeamLineHitMan() : HitMan()
{
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
HodoscopeLikeContainer *BeamLineHitMan::HodoContainer( const int &cid )
{
  switch( cid ){
  case CID_BHD:      return &BHDContainer;
  case CID_BHDpost:  return &BHDpostContainer;    
  case CID_T0:       return &T0Container;
  case CID_T0pre:    return &T0preContainer;
  case CID_T0post:   return &T0postContainer;
  case CID_DEF:      return &E0Container;
    //  case CID_E0:       return &E0Container;
  case CID_CVC:      return &CVCContainer;
  case CID_NC:       return &NCContainer;
  case CID_PC:       return &PCContainer;
  case CID_BPD:      return &BPDContainer;
  case CID_BVC:      return &BVContainer;
  case CID_WVC:      return &WVContainer;
  case CID_HVC1:     return &HVC1Container;
  case CID_HVC2:     return &HVC2Container;
  case CID_LB:       return &LBContainer;
  case CID_BD:       return &BDContainer;
  case CID_TEMP1:    return &Temp1Container;
  case CID_TEMP2:    return &Temp2Container;
  case CID_START:    return &StartContainer;
  default:           return 0;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
CherenkovLikeContainer *BeamLineHitMan::ChereContainer( const int &cid )
{
  switch( cid ){
  case CID_LC1:  return &LC1Container;
  case CID_LC2:  return &LC2Container;
  case CID_AC:   return &ACContainer;
  case CID_WC:   return &WCContainer;
  case CID_GC:   return &GCContainer;
  default:           return 0;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
MTDCLikeContainer *BeamLineHitMan::MTDCContainer( const int &cid )
{
  switch( cid ){
  case CID_BHDmul:  return &BHDMTDCContainer;
  case CID_T0mul:   return &T0MTDCContainer;
  case CID_BVCmul:  return &BVCMTDCContainer;
  case CID_HVC1mul: return &HVC1MTDCContainer;
  case CID_HVC2mul: return &HVC2MTDCContainer;
  case CID_REFmul:  return &REFMTDCContainer;
  default:          return 0;
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
ChamberLikeContainer *BeamLineHitMan::ChmContainer( const int &cid, const int &layer )
{
  switch( cid ){
  case CID_BLC1a:  
    return (0<layer&&layer<=NumOfBLCLayers) ? &BLC1aContainer[layer-1] : 0;
  case CID_BLC1b:
    return (0<layer&&layer<=NumOfBLCLayers) ? &BLC1bContainer[layer-1] : 0;
  case CID_BLC2a:  
    return (0<layer&&layer<=NumOfBLCLayers) ? &BLC2aContainer[layer-1] : 0;
  case CID_BLC2b:
    return (0<layer&&layer<=NumOfBLCLayers) ? &BLC2bContainer[layer-1] : 0;
  case CID_BPC:  
    return (0<layer&&layer<=NumOfBLCLayers) ? &BPCContainer[layer-1] : 0;
  case CID_FDC1:
    return (0<layer&&layer<=NumOfBLCLayers) ? &FDC1Container[layer-1] : 0;
  default:    
    return 0;
  }
}
// -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::SetSimulatedResolution( ConfMan *conf )
{
  for(int i=0;i<nBHodoIDList;i++){
    HodoscopeLikeContainer *container=HodoContainer(BHodoIDList[i]);
    if(container){
      for( int j=0; j<(int)container->size(); j++ ){
	(*container)[j].SetSimulatedResolution(conf);
      }
    }
  }
  for(int i=0;i<nChereIDList;i++){
    CherenkovLikeContainer *container=ChereContainer(ChereIDList[i]);
    if(container){
      for( int j=0; j<(int)container->size(); j++ ){
	//	(*container)[j].SetSimulatedResolution(conf);
      }
    }
  }
  for(int i=0;i<nBLDCIDList;i++){
    for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
      ChamberLikeContainer *container=BLDCContainer(BLDCIDList[i],lay);
      if(container){
	for( int j=0; j<(int)container->size(); j++ ){
	  //	(*container)[j].SetSimulatedResolution(conf);
	}
      }
    }
  }
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
bool BeamLineHitMan::Calc( ConfMan *conf, double toffs )
{
  for(int i=0;i<nBHodoIDList;i++){
    HodoscopeLikeContainer *container=HodoContainer(BHodoIDList[i]);
    if(container){
      for( int j=0; j<(int)container->size(); j++ ){
	(*container)[j].Calc(conf);
      }
    }
  }
  for(int i=0;i<nChereIDList;i++){
    CherenkovLikeContainer *container=ChereContainer(ChereIDList[i]);
    if(container){
      for( int j=0; j<(int)container->size(); j++ ){
	(*container)[j].Calc(conf);
      }
    }
  }
  CalcT0Time();
  CalcBLDC(conf,T0Time);
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::RemoveNoTDCData()
{
  HodoscopeSegContainer.clear();
  CherenkovSegContainer.clear();
  unsigned int key;
  for(int i=0;i<nBHodoIDList;i++){
    HodoscopeLikeContainer::iterator ittr;
    HodoscopeLikeContainer *container=HodoContainer(BHodoIDList[i]);
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
  for(int i=0;i<nChereIDList;i++){
    CherenkovLikeContainer::iterator ittr;
    CherenkovLikeContainer *container=ChereContainer(ChereIDList[i]);
    if(container){
      for( ittr=container->end()-1; ittr!=container->begin()-1; --ittr ){
	if(!ittr->CheckRange()) container->erase(ittr);
      }
    }
    for( unsigned int i=0; i<container->size(); i++ ){
      key = KEY( (*container)[i].cid(), (*container)[i].seg() );
      CherenkovSegContainer[key] = i;
    }
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::CalcBLDC(ConfMan *conf, const double &toffs){
  for(int i=0;i<nBLDCIDList;i++){
    for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
      ChamberLikeContainer *container=BLDCContainer(BLDCIDList[i],lay);
      if(container){
	for( int j=0; j<(int)container->size(); j++ ){
	  (*container)[j].Calc(conf,toffs);
	}
      }
    }
  }
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::Clear()
{
  for(int i=0;i<nBHodoIDList;i++){
    HodoscopeLikeContainer *container=HodoContainer(BHodoIDList[i]);
    if(container) container->clear();
  }
  for(int i=0;i<nChereIDList;i++){
    CherenkovLikeContainer *container=ChereContainer(ChereIDList[i]);
    if(container) container->clear();
  }
  for(int i=0;i<nMTDCIDList;i++){
    MTDCLikeContainer *container=MTDCContainer(MTDCIDList[i]);
    if(container) container->clear();
  }
  for(int i=0;i<nBLDCIDList;i++){
    for( int lay=1; lay<=NumOfBLCLayers; lay++ ){
      ChamberLikeContainer *container=BLDCContainer(BLDCIDList[i],lay);
      if(container) container->clear();
    }
  }
  T0Time=0.;
  T0Seg=-1;
  HitMan::Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- //
void BeamLineHitMan::CheckContainerSize()
{
  std::cout << " Size -> "
	    << " BHD:"  << BHDContainer.size()
	    << " T0: "  << T0Container.size()
	    << " E0: "  << E0Container.size()
	    << " CVC:"  << CVCContainer.size()
	    << " NC: "  << NCContainer.size()
	    << " PC: "  << PCContainer.size()
	    << " BPD:"  << BPDContainer.size()
	    << " BVC:"  << BVContainer.size()
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
  std::cout << "          "; std::cout << " BPC:  ";
  for( int i=0; i<NumOfBLCLayers;  i++ ){ std::cout << " , " << BPCContainer[i].size(); }
  std::cout << std::endl;
  std::cout << "          "; std::cout << " FDC1: ";
  for( int i=0; i<NumOfBLCLayers;  i++ ){ std::cout << " , " << FDC1Container[i].size(); }
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

