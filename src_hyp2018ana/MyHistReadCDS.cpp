#include "MyHistReadCDS.h"

using namespace std;

void initHistReadCDS()
{
  new TH2F("CDS_mass2_mom",     "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pim", "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_km",  "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pip", "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_p",   "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_d",   "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);

  new TH2F("CDS_mass2_mom_Z40",     "CDS mass2 mom",     1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_Zdiff10",     "CDS mass2 mom", 1000, -1.0, 9.0, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_Zmatch",     "CDS mass2 mom",  1000, -1.0, 9.0, 1000, -1.0, 1.0);

  new TH2F("Vtx_ZX", "Vertex Z-X Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY", "Vertex Z-Y Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_XY", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);

  new TH2F("Vtx_ZX_fp", "Vertex Z-X Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY_fp", "Vertex Z-Y Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_XY_fp", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);
  new TH2F("Vtx_ZX_f", "Vertex Z-X Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY_f", "Vertex Z-Y Plane", 1000,  -25,   25, 500, -12.5, 12.5);
  new TH2F("Vtx_XY_f", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);

  new TH1F("Vtx_Z_center", "Vertex Z", 2000, -25, 25);
  new TH1F("Vtx_X_center",  "Vertex X Plane",  5000, -12.5, 12.5);
  new TH1F("Vtx_Y_center",  "Vertex X Plane",  5000, -12.5, 12.5);

  new TH2F("Vtx_XY_z6",     "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);
  new TH2F("Vtx_XY_z0",     "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);

  new TH1F("CDS_IM_ppim", "CDS p #pi^{-} IM",       10000, 1.0, 2.0);
  new TH1F("CDS_IM_pipi", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);

  new TH1F("CDS_DCA",    "CDS-beam DCA", 1000, 0.0, 10);
  new TH1F("CDS_DCA_fp", "CDS-beam DCA", 1000, 0.0, 10);

  for( int seg=1; seg<=36; seg++ ){
    new TH2F(Form("CDC_Z_CDH%d_tsub",seg), Form("CDC Z vs CDH seg%d tsub", seg), 1000, -50, 50, 1000, -5, 5);
    new TH2F(Form("CDC_Z_CDH%d_diff",seg), Form("CDC Z vs Z pos diff CDH seg%d", seg), 1000, -50, 50, 1000, -25, 25);
  }

  new TNtuple("tup_CDS_mass2_mom", "", "beta:mass2:mom");
}

void fillHistReadCDS(AnaInfo *anaInfo, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan)
{
  TNtuple *tup;

  if( anaInfo->minDCA() ){
    CDSInfo *minDCA=anaInfo->minDCA();
    TVector3 vtx=minDCA->vertexBeam();

    MyHistTools::fillTH("Vtx_ZX", vtx.Z(), vtx.X());
    MyHistTools::fillTH("Vtx_ZY", vtx.Z(), vtx.Y());
    MyHistTools::fillTH("Vtx_XY", vtx.X(), vtx.Y());
    if( anaInfo->nFCharge()==1 && 
	FC_P_MIN<anaInfo->forwardCharge(0)->mass2byAng() && anaInfo->forwardCharge(0)->mass2byAng()<FC_P_MAX ){
      MyHistTools::fillTH("Vtx_ZX_fp", vtx.Z(), vtx.X());
      MyHistTools::fillTH("Vtx_ZY_fp", vtx.Z(), vtx.Y());
      MyHistTools::fillTH("Vtx_XY_fp", vtx.X(), vtx.Y());
    }

    if( GeomTools::GetID(vtx)==CID_Fiducial ){
      MyHistTools::fillTH("Vtx_ZX_f", vtx.Z(), vtx.X());
      MyHistTools::fillTH("Vtx_ZY_f", vtx.Z(), vtx.Y());
      MyHistTools::fillTH("Vtx_XY_f", vtx.X(), vtx.Y());
    }

    double r=sqrt(pow(vtx.X()-0.4, 2)+pow(vtx.Y()-0.5,2));
    if( r<3.0 ){
      MyHistTools::fillTH("Vtx_Z_center", vtx.Z());
    }
    if( -0.2<vtx.Z() && vtx.Z()<0.2 ){
      MyHistTools::fillTH("Vtx_XY_z0", vtx.X(), vtx.Y());
    }
    if( -6.2<vtx.Z() && vtx.Z()<-5.8 ){
      MyHistTools::fillTH("Vtx_XY_z6", vtx.X(), vtx.Y());
    }
    if( -8.0<vtx.Z() && vtx.Z()<-4.8 && -0.75<vtx.Y() && vtx.Y()<-0.25 ){
      MyHistTools::fillTH("Vtx_X_center", vtx.X());
    }
    if( -8.0<vtx.Z() && vtx.Z()<-4.8 && -0.65<vtx.X() && vtx.X()<-0.15 ){
      MyHistTools::fillTH("Vtx_Y_center", vtx.Y());
    }
  } 

  for( int i=0; i<anaInfo->nCDS(); i++ ){
    CDSInfo *cds=anaInfo->CDS(i);
    MyHistTools::fillTH("CDS_DCA", cds->dca());
    if( anaInfo->nFCharge()==1 ){
      ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
      if( 0.4<fcInfo->mass2byAng() && fcInfo->mass2byAng()<2.0 ){
	MyHistTools::fillTH("CDS_DCA_fp", cds->dca());
      }
    }
  }

  for( int i=0; i<anaInfo->nCDS(); i++ ){
    if( !anaInfo->CDS(i)->flag() ) continue;
    if( GeomTools::GetID(anaInfo->CDS(i)->vertexBeam())!=CID_Fiducial ) continue;
    CDSInfo *cds=anaInfo->CDS(i);
    MyHistTools::fillTH("CDS_mass2_mom", cds->mass2(), cds->mom());
    if( cds->pid()==CDS_PiMinus ){
      MyHistTools::fillTH("CDS_mass2_mom_pim", cds->mass2(), cds->mom());
    }
    if( cds->pid()==CDS_PiPlus ){
      MyHistTools::fillTH("CDS_mass2_mom_pip", cds->mass2(), cds->mom());
    }
    if( cds->pid()==CDS_Kaon ){
      MyHistTools::fillTH("CDS_mass2_mom_km", cds->mass2(), cds->mom());
    }
    if( cds->pid()==CDS_Proton ){
      MyHistTools::fillTH("CDS_mass2_mom_p", cds->mass2(), cds->mom());
    }
    if( cds->pid()==CDS_Deuteron ){
      MyHistTools::fillTH("CDS_mass2_mom_d", cds->mass2(), cds->mom());
    }

    HodoscopeLikeHit *hit=cds->CDH(cdsMan);
    CDSTrack *track=cds->track(cdstrackMan);
    TVector3 vertexCDH=track->CDHVertex();
    double ctsub=hit->ctsub();
    double hitpos=hit->hitpos();
    int seg=hit->seg();
    MyHistTools::fillTH(Form("CDC_Z_CDH%d_tsub", seg), vertexCDH.Z(), ctsub);
    MyHistTools::fillTH(Form("CDC_Z_CDH%d_diff", seg), vertexCDH.Z(), vertexCDH.Z()-hitpos);

    tup=(TNtuple*)gFile->Get("tup_CDS_mass2_mom"); tup-> Fill(cds->beta(), cds->mass2(), cds->mom());

    if( fabs(vertexCDH.Z())>40 ){
      MyHistTools::fillTH("CDS_mass2_mom_Z40", cds->mass2(), cds->mom());
    }
    if( fabs(vertexCDH.Z()-hitpos)>10 && fabs(vertexCDH.Z())<40 ){
      MyHistTools::fillTH("CDS_mass2_mom_Zdiff10", cds->mass2(), cds->mom());
    }
    if( fabs(vertexCDH.Z()-hitpos)<10 && fabs(vertexCDH.Z())<40 ){
      MyHistTools::fillTH("CDS_mass2_mom_Zmatch", cds->mass2(), cds->mom());
    }
  }

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_PiPlus); i++ ){
    CDS2Info *cds2=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, i);
    if( GeomTools::GetID(cds2-> vertexBeam())==CID_Fiducial && cds2->flag() ){
      MyHistTools::fillTH("CDS_IM_pipi", cds2->im());
    }
  }

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_Proton); i++ ){
    CDS2Info *cds2=anaInfo->CDS2(CDS_PiMinus, CDS_Proton, i);
    if( GeomTools::GetID(cds2-> vertexBeam())==CID_Fiducial && cds2->flag() ){
      MyHistTools::fillTH("CDS_IM_ppim", cds2->im());
    }
  }
}
