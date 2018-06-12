#include "MyHistTools.h"

void MyHistTools::initCDS()
{
  new TH1F("CDS_chi2", "CDS chi-square", 1000, 0.0, 500);

  new TH2F("CDS_mass2_mom", "CDS mass^{2} vs mom", 1000, -0.5, 4.5, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pim", "CDS mass^{2} vs mom", 1000, -0.5, 4.5, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pip", "CDS mass^{2} vs mom", 1000, -0.5, 4.5, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_km",  "CDS mass^{2} vs mom", 1000, -0.5, 4.5, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_p",   "CDS mass^{2} vs mom", 1000, -0.5, 4.5, 1000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_d",   "CDS mass^{2} vs mom", 1000, -0.5, 4.5, 1000, -1.0, 1.0);

  new TH1F("CDS_ppim_IM", "p #pi^{-} IM",       1000, 1.0, 2.0);
  new TH1F("CDS_pipi_IM", "#pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);
}

void MyHistTools::fillCDS(AnaInfo *anaInfo, CDSTrackingMan *cdstrackMan)
{
  TH1F *h1;
  TH2F *h2;
  for( int i=0; i<anaInfo->nCDS(); i++ ){
    CDSInfo *info=anaInfo->CDS(i);
    CDSTrack *track = info->track(cdstrackMan);
    h1 = (TH1F*)gFile-> Get("CDS_chi2"), h1-> Fill(track->Chi());
    if( info->flag() && MyTools::isFiducial(info) ){
      h2 = (TH2F*)gFile-> Get("CDS_mass2_mom"), h2-> Fill(info->mass2(), info->mom());
      if( info->pid()==CDS_PiMinus ){ h2 = (TH2F*)gFile-> Get("CDS_mass2_mom_pim"), h2-> Fill(info->mass2(), info->mom()); }
      if( info->pid()==CDS_PiPlus  ){ h2 = (TH2F*)gFile-> Get("CDS_mass2_mom_pip"), h2-> Fill(info->mass2(), info->mom()); }
      if( info->pid()==CDS_Kaon    ){ h2 = (TH2F*)gFile-> Get("CDS_mass2_mom_km"),  h2-> Fill(info->mass2(), info->mom()); }
      if( info->pid()==CDS_Proton  ){ h2 = (TH2F*)gFile-> Get("CDS_mass2_mom_p"),   h2-> Fill(info->mass2(), info->mom()); }
      if( info->pid()==CDS_Deuteron){ h2 = (TH2F*)gFile-> Get("CDS_mass2_mom_d"),   h2-> Fill(info->mass2(), info->mom()); }
    }
  }

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_PiPlus); i++ ){
    CDS2Info *info=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, i);
    if( info->flag() && GeomTools::GetID(info->vertexBeam())==CID_Fiducial ){
      h1 = (TH1F*)gFile-> Get("CDS_pipi_IM"), h1-> Fill(info->im());
    }
  }

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_Proton); i++ ){
    CDS2Info *info=anaInfo->CDS2(CDS_PiMinus, CDS_Proton, i);
    if( info->flag() && GeomTools::GetID(info->vertexBeam())==CID_Fiducial ){
      h1 = (TH1F*)gFile-> Get("CDS_ppim_IM"), h1-> Fill(info->im());
    }
  }
}
