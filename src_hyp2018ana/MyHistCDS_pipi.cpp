#include "MyHistCDS_pipi.h"

void initHistCDS_pipi()
{
  new TH1F("CDS_IM_pipi", "CDS #pi^{+} #pi^{-} IM", 10000, 0.0, 1.0);
}

void fillCDS_pipi(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_PiPlus); i++ ){
    CDS2Info *cds2=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, i);
    if( cds2->flag() && GeomTools::GetID(cds2->vertexBeam())==CID_Fiducial ){
      MyHistTools::fillTH("CDS_IM_pipi", cds2->im());
    }  
  }
}
