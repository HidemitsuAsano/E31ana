#include "MyHistCDS_ppim.h"

using namespace std;

void initHistCDS_ppim()
{
  new TH1F("CDS_IM_ppim", "CDS p #pi^{-} IM", 10000, 1.0, 2.0);
}

void fillCDS_ppim(AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_Proton); i++ ){
    CDS2Info *cds2=anaInfo->CDS2(CDS_PiMinus, CDS_Proton, i);
    if( cds2->flag() && GeomTools::GetID(cds2->vertexBeam())==CID_Fiducial ){
      MyHistTools::fillTH("CDS_IM_ppim", cds2->im());
    }
  }
}
