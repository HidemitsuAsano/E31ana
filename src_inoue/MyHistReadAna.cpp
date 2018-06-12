#include "MyHistReadAna.h"

void initHistReadAna()
{
  new TH2F("CDS_mass2_mom",     "CDS mass2 mom", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_pim", "CDS mass2 mom", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_km",  "CDS mass2 mom", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_pip", "CDS mass2 mom", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_p",   "CDS mass2 mom", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_d",   "CDS mass2 mom", 1000, -1.0, 9.0, 3000, -1.5, 1.5);

  new TH2F("Vtx_ZX", "Vertex Z-X Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY", "Vertex Z-Y Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_XY", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);
  new TH2F("Vtx_ZX_fp", "Vertex Z-X Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY_fp", "Vertex Z-Y Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_XY_fp", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);

  new TH2F("Vtx_ZX_f", "Vertex Z-X Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY_f", "Vertex Z-Y Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_XY_f", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);
}

void fillReadAna_CDS(AnaInfo *anaInfo)
{
  if( anaInfo->minDCA() ){
    CDSInfo *minDCA=anaInfo->minDCA();
    TVector3 vtx=minDCA->vertexBeam();

    MyHistTools::fillTH("Vtx_ZX", vtx.Z(), vtx.X());
    MyHistTools::fillTH("Vtx_ZY", vtx.Z(), vtx.Y());
    MyHistTools::fillTH("Vtx_XY", vtx.X(), vtx.Y());
    if( anaInfo->nForwardCharged()==1 && FC_P_MIN<anaInfo->forwardCharge(0)->Mass2byAng() && anaInfo->forwardCharge(0)->Mass2byAng()<FC_P_MAX ){
      MyHistTools::fillTH("Vtx_ZX_fp", vtx.Z(), vtx.X());
      MyHistTools::fillTH("Vtx_ZY_fp", vtx.Z(), vtx.Y());
      MyHistTools::fillTH("Vtx_XY_fp", vtx.X(), vtx.Y());
    }

    if( GeomTools::GetID(vtx)==CID_Fiducial ){
      MyHistTools::fillTH("Vtx_ZX_f", vtx.Z(), vtx.X());
      MyHistTools::fillTH("Vtx_ZY_f", vtx.Z(), vtx.Y());
      MyHistTools::fillTH("Vtx_XY_f", vtx.X(), vtx.Y());
    }
  } 

  for( int i=0; i<anaInfo->nCDS(); i++ ){
    if( !anaInfo->CDS(i)->flag() ) continue;
    if( anaInfo->CDS(i)->vertexBeam()!=CID_Fiducial ) continue;
    CDSInfo *cds=anaInfo->CDS(i);
    MyHistTools::fillTH("CDS_mass2_mom", cds->mass2(), cds->mom());
    if( cds->pid()==CDS_PiMinus ){
      MyHistTools::fillTH("CDS_mass2_mom_pim", cds->mass2(), cds->mom());
    }
    if( cds->pid()==CDS_PiPlus ){
      MyHistTools::fillTH("CDS_mass2_mom_pip", cds->mass2(), cds->mom());
    }
  }
}
