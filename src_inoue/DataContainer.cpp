/* DataContainer.cpp */
#include "DataContainer.h"

ClassImp(DataContainer);

DataContainer::DataContainer()
{
  Clear();
}

DataContainer::~DataContainer()
{
}

bool DataContainer::DoTrackingBL(BeamLineHitMan *blMan, ConfMan *confMan)
{
  bool result = bltrackMan.DoTracking(blMan, confMan, true, false);
  return result;
}

bool DataContainer::ExecuteCDS(CDSHitMan *cdsMan, ConfMan *confMan)
{
  bool result = trackingMan.Execute(cdsMan, confMan);
  return result;
}

void DataContainer::CalcVertex_beam(ConfMan *confMan)
{
  for( int i=0; i<trackingMan.nGoodTrack(); i++ ){
    trackingMan.CalcVertex_beam(trackingMan.GoodTrackID(i), &bltrackMan, confMan);
  }
}

void DataContainer::SetHodoscope(BeamLineHitMan *blMan, CDSHitMan *cdsMan)
{
  // for BeamLine Hodoscope
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      BHDContainer.push_back(*blMan->BHD(i));
    }
  }
  for( int i=0; i<blMan->nBHDpost(); i++ ){
    if( blMan->BHDpost(i)->CheckRange() ){
      BHDpostContainer.push_back(*blMan->BHDpost(i));
    }
  }
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      T0Container.push_back(*blMan->T0(i));
    }
  }
  for( int i=0; i<blMan->nT0pre(); i++ ){
    if( blMan->T0pre(i)->CheckRange() ){
      T0preContainer.push_back(*blMan->T0pre(i));
    }
  }
  for( int i=0; i<blMan->nT0post(); i++ ){
    if( blMan->T0post(i)->CheckRange() ){
      T0postContainer.push_back(*blMan->T0post(i));
    }
  }
  for( int i=0; i<blMan->nE0(); i++ ){
    if( blMan->E0(i)->CheckRange() ){
      E0Container.push_back(*blMan->E0(i));
    }
  }
  for( int i=0; i<blMan->nCVC(); i++ ){
    if( blMan->CVC(i)->CheckRange() ){
      CVCContainer.push_back(*blMan->CVC(i));
    }
  }
  for( int i=0; i<blMan->nBPD(); i++ ){
    if( blMan->BPD(i)->CheckRange() ){
      BPDContainer.push_back(*blMan->BPD(i));
    }
  }
  for( int i=0; i<blMan->nPC(); i++ ){
    if( blMan->PC(i)->CheckRange() ){
      PCContainer.push_back(*blMan->PC(i));
    }
  }
  for( int i=0; i<blMan->nBVC(); i++ ){
    if( blMan->BVC(i)->CheckRange() ){
      BVCContainer.push_back(*blMan->BVC(i));
    }
  }
  for( int i=0; i<blMan->nBD(); i++ ){
    if( blMan->BD(i)->CheckRange() ){
      BDContainer.push_back(*blMan->BD(i));
    }
  }
  for( int i=0; i<blMan->nLB(); i++ ){
    if( blMan->LB(i)->CheckRange() ){
      LBContainer.push_back(*blMan->LB(i));
    }
  }
  for( int i=0; i<blMan->nWVC(); i++ ){
    if( blMan->WVC(i)->CheckRange() ){
      WVCContainer.push_back(*blMan->WVC(i));
    }
  }
  for( int i=0; i<blMan->nHVC1(); i++ ){
    if( blMan->HVC1(i)->CheckRange() ){
      HVC1Container.push_back(*blMan->HVC1(i));
    }
  }
  for( int i=0; i<blMan->nHVC2(); i++ ){
    if( blMan->HVC2(i)->CheckRange() ){
      HVC2Container.push_back(*blMan->HVC2(i));
    }
  }
  for( int i=0; i<blMan->nTemp1(); i++ ){
    if( blMan->Temp1(i)->CheckRange() ){
      Temp1Container.push_back(*blMan->Temp1(i));
    }
  }
  for( int i=0; i<blMan->nTemp2(); i++ ){
    if( blMan->Temp2(i)->CheckRange() ){
      Temp2Container.push_back(*blMan->Temp2(i));
    }
  }

  // for Neutron Counter
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      int seg = blMan-> NC(i)-> seg();
      int id = (seg-1)/NumOfNCSegmentsInLayer;
      NCContainer[id].push_back(*blMan->NC(i));
    }
  }

  // for CDS Hodoscope
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    if( cdsMan->CDH(i)-> CheckRange() ){
      CDHContainer.push_back(*cdsMan->CDH(i));
    }
  }
  for( int i=0; i<cdsMan->nIH(); i++ ){
    int tdc = cdsMan->IH(i)-> tdcu();
    if( 0<tdc && tdc<4095 ){
      IHContainer.push_back(*cdsMan->IH(i));
    }
  }

  // for Cherenkov Counter
  for( int i=0; i<blMan->nLC1(); i++ ){
    LC1Container.push_back(*blMan->LC1(i));
  }
  for( int i=0; i<blMan->nLC2(); i++ ){
    LC1Container.push_back(*blMan->LC2(i));
  }
  for( int i=0; i<blMan->nAC(); i++ ){
    LC1Container.push_back(*blMan->AC(i));
  }
  for( int i=0; i<blMan->nWC(); i++ ){
    LC1Container.push_back(*blMan->WC(i));
  }
  for( int i=0; i<blMan->nGC(); i++ ){
    LC1Container.push_back(*blMan->GC(i));
  }
}

int DataContainer::nNC() const
{
  int num = 0;
  for( int i=0; i<NumOfNCLayer; i++ ){
    num += NCContainer[i].size();
  }
  return num;
}

HodoscopeLikeHit* DataContainer::NC(const int &i)
{
  int num=0;
  for( int j=0; j<NumOfNCLayer; j++ ){
    if( i<num+(int)NCContainer[j].size() ){
      return &NCContainer[j][i-num];
    }
    num += NCContainer[j].size();
  }
  return 0;
}

void DataContainer::Clear()
{
  trackingMan.Clear();
  bltrackMan.Clear();

  BHDContainer.clear();
  BHDpostContainer.clear();
  T0Container.clear();
  T0preContainer.clear();
  T0postContainer.clear();
  E0Container.clear();
  CVCContainer.clear();
  BPDContainer.clear();
  PCContainer.clear();
  BVCContainer.clear();
  BDContainer.clear();
  LBContainer.clear();
  WVCContainer.clear();
  HVC1Container.clear();
  HVC2Container.clear();
  Temp1Container.clear();
  Temp2Container.clear();

  for( int i=0; i<NumOfNCLayer; i++ ){
    NCContainer[i].clear();
  }

  CDHContainer.clear();
  IHContainer.clear();

  LC1Container.clear();
  LC2Container.clear();
  ACContainer.clear();
  WCContainer.clear();
  GCContainer.clear();

}
