#include "MyHistReadVtxDEF.h"

void initHistReadVtxDEF()
{
  new TH2F("Vtx_XY_fp_DEF", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);

  new TH2F("Vtx_ZX_fp_km", "Vertex X-Y Plane", 500, -25.0, 25.0, 500, -12.5, 12.5);
  new TH2F("Vtx_ZY_fp_km", "Vertex X-Y Plane", 500, -25.0, 25.0, 500, -12.5, 12.5);
  new TH2F("Vtx_XY_fp_DEF_km", "Vertex X-Y Plane", 500, -12.5, 12.5, 500, -12.5, 12.5);

  new TH1F("KN_MM_DEF",      "p(K^{-}, n)\"X\"",     2000, 0.0, 2.0);
  new TH1F("KN_MM_DEF_pipi", "p(K^{-}, n)\"X\"",     2000, 0.0, 2.0);
  new TH1F("KN_MM_DEF_pim",  "p(K^{-}, n)\"X\"",     2000, 0.0, 2.0);
  new TH1F("KN_MM_DEF_pip",  "p(K^{-}, n)\"X\"",     2000, 0.0, 2.0);
  new TH1F("KN_MM_DEF_km",   "p(K^{-}, n)\"X\"",     2000, 0.0, 2.0);
  new TH1F("KN_MM_DEF_p",    "p(K^{-}, n)\"X\"",     2000, 0.0, 2.0);
  new TH1F("CDS_IM_pipi_DEF", "", 1000, 0.0, 1.0);
  new TH1F("KCDSK0_MM_DEF", "", 1000, 0.5, 1.5);

  new TH1F("KP_MM_DEF",     "p(K^{-}, p)\"X\"",     2000, 0.0, 2.0);
  new TH1F("KP_MM_DEF_pim", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_DEF_pip", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_DEF_km",  "p(K^{-}, p)\"X\"",  2000, 0.0, 2.0);
  new TH1F("KP_MM_DEF_p",   "p(K^{-}, p)\"X\"",   2000, 0.0, 2.0);

  new TH1F("KPpim_MM_DEF", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("Ppim_IM_DEF", "p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("KP_MM_DEF_2pim_pip", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);

  new TH1F("KP_MM_DEF_km_1hit",      "p(K^{-}, p)\"X\"",  2000, 0.0, 2.0);
  new TH1F("KP_MM_DEF_km_1hit_true", "p(K^{-}, p)\"X\"",  2000, 0.0, 2.0);
  new TH1F("Kkm_MM_DEF_fp", "p(K^{-}, K^{-})\"X\"",  2000, 0.0, 2.0);
  new TH1F("Kkm_MM_DEF_fp_mmKm", "p(K^{-}, K^{-})\"X\"",  2000, 0.0, 2.0);

  new TH1F("Ene_kp_kp_DEF", "Energy Conversation", 2000, -1.0, 1.0);

  new TH1F("KP_KP_E_DEF", "K^{-} p #rightarrow K^{-} p Energy", 2000, -1.0, 1.0);

  new TH2F("mmKm_kin_DEF", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("Km_kin_DEF", "K^{-} cos#theta vs mom",   1000, -1.0, 1.0, 1000, 0.0, 1.0);

  new TH2F("KP_MM_Kkm_MM_DEF", "", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH1F("CDS_chi2_vtxDEF_fp", "", 1000, 0, 100);
  new TH1F("CDS_chi2_vtxDEF_fn", "", 1000, 0, 100);
  new TH1F("CDS_DCA_vtxDEF_fp", "",  1000, 0, 100);
  new TH1F("CDS_DCA_vtxDEF_fn", "",  1000, 0, 100);
  new TH1F("CDS_anaLevel_DEF_fn", "", 20, 0, 20);
  new TH1F("CDS_anaLevel_DEF_fp", "", 20, 0, 20);
  new TH2F("CDS_mass2_mom_vtxDEF_fp", "CDS mass^{2} vs mom", 1000, -1.0, 4.0, 300, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_vtxDEF_fn", "CDS mass^{2} vs mom", 1000, -1.0, 4.0, 300, -1.5, 1.5);
}

void fillHistReadVtxDEF(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo, 
			BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan,
			CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan)
{
  if( anaInfo->minDCA() ){
    CDSInfo *minDCA=anaInfo->minDCA();
    TVector3 vtx=minDCA->vertexBeam();

    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector target_lmom=TLorentzVector(0.0, 0.0, 0.0, pMass);

    if( anaInfo->nFNeutral()==1 && anaInfo->forwardNeutral(0)->pid()==F_Neutron ){
      ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
      TLorentzVector fn_lmom=fnInfo->lmom();
      if( -17.5<vtx.Z() && vtx.Z()<-14 ){
	MyHistTools::fillTH("CDS_DCA_vtxDEF_fn", minDCA->dca());

	for( int i=0; i<anaInfo->nCDS(); i++ ){
	  CDSInfo *cds=anaInfo->CDS(i);
	  CDSTrack *track=cds->track(cdstrackMan);
	  int ana_level=MyTools::CDSanaLevel(track, cdsMan, anaInfo->beam(0), bltrackMan, conf, false);
	  MyHistTools::fillTH("CDS_anaLevel_DEF_fn", ana_level);
	  MyHistTools::fillTH("CDS_chi2_vtxDEF_fn", track->Chi());
	  MyHistTools::fillTH("CDS_mass2_mom_vtxDEF_fn", cds->mass2(), cds->mom());
	}

	double mm=(beam_lmom+target_lmom-fn_lmom).M();
	MyHistTools::fillTH("KN_MM_DEF", mm);
	if( anaInfo->nCDS(CDS_PiMinus)>0 ) MyHistTools::fillTH("KN_MM_DEF_pim", mm);
	if( anaInfo->nCDS(CDS_PiPlus)>0  ) MyHistTools::fillTH("KN_MM_DEF_pip", mm);
	if( anaInfo->nCDS(CDS_Kaon)>0    ) MyHistTools::fillTH("KN_MM_DEF_km", mm);
	if( anaInfo->nCDS(CDS_Proton)>0  ) MyHistTools::fillTH("KN_MM_DEF_p", mm);
	
	if( anaInfo->nCDS(CDS_PiMinus)==1 && anaInfo->nCDS(CDS_PiPlus)==1 ){
	  CDS2Info *pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0);
	  MyHistTools::fillTH("KN_MM_DEF_pipi", mm);
	  MyHistTools::fillTH("CDS_IM_pipi_DEF", pipi->im());
	  
	  if( K0_MIN<pipi->im() && pipi->im()<K0_MAX ){
	    MyHistTools::fillTH("KCDSK0_MM_DEF", (beam_lmom+target_lmom-pipi->lmom()).M());
	  }
	}
      }
    }
  
    if( anaInfo->nFCharge()==1 &&
        FC_P_MIN<anaInfo->forwardCharge(0)->mass2byAng() && anaInfo->forwardCharge(0)->mass2byAng()<FC_P_MAX ){
      ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);

      TLorentzVector fc_lmom=fcInfo->lmom();

      double mm=(beam_lmom+target_lmom-fc_lmom).M();
      if( anaInfo->nCDS(CDS_Kaon)>0 ){
	MyHistTools::fillTH("Vtx_ZX_fp_km", vtx.Z(), vtx.X());
	MyHistTools::fillTH("Vtx_ZY_fp_km", vtx.Z(), vtx.Y());
      }

      if( -17.5<vtx.Z() && vtx.Z()<-14 ){
	MyHistTools::fillTH("CDS_DCA_vtxDEF_fp", minDCA->dca());
	for( int i=0; i<anaInfo->nCDS(); i++ ){
	  CDSInfo *cds=anaInfo->CDS(i);
	  CDSTrack *track=anaInfo->CDS(i)->track(cdstrackMan);
	  int ana_level=MyTools::CDSanaLevel(track, cdsMan, anaInfo->beam(0), bltrackMan, conf, false);
	  MyHistTools::fillTH("CDS_chi2_vtxDEF_fp", track->Chi());
	  MyHistTools::fillTH("CDS_mass2_mom_vtxDEF_fp", cds->mass2(), cds->mom());
	  MyHistTools::fillTH("CDS_anaLevel_DEF_fp", ana_level);
	}

	MyHistTools::fillTH("Vtx_XY_fp_DEF", vtx.X(), vtx.Y());
	MyHistTools::fillTH("KP_MM_DEF", mm);
	if( anaInfo->nCDS(CDS_PiMinus)>0 ){
	  CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, 0);
	  TLorentzVector pim_lmom=pim->lmom();
	  double kppim_mm=(beam_lmom+target_lmom-fc_lmom-pim_lmom).M();
	  MyHistTools::fillTH("KP_MM_DEF_pim", mm);
	  MyHistTools::fillTH("KPpim_MM_DEF", kppim_mm);
	  MyHistTools::fillTH("Ppim_IM_DEF", (fc_lmom+pim_lmom).M());
	}
	if( anaInfo->nCDS(CDS_PiPlus)>0 ){
	  MyHistTools::fillTH("KP_MM_DEF_pip", mm);
	}
	if( anaInfo->nCDS(CDS_PiMinus)>1 && anaInfo->nCDS(CDS_PiPlus)>0 ){
	  MyHistTools::fillTH("KP_MM_DEF_2pim_pip", mm);
	}
      
	if( anaInfo->nCDS(CDS_Kaon)>0 ){
	  MyHistTools::fillTH("KP_MM_DEF_km", mm);
	  if( anaInfo->nCDS(CDS_Kaon)==1 && anaInfo->nCDS()==1 ){
	    CDSInfo *km = anaInfo-> CDS(CDS_Kaon, 0);
	    TLorentzVector km_lmom=km->lmom();
	    TLorentzVector mm_lmom = beam_lmom+target_lmom-fc_lmom;

	    MyHistTools::fillTH("KP_MM_DEF_km_1hit", mm);
	    MyHistTools::fillTH("Kkm_MM_DEF_fp", (beam_lmom+target_lmom-km_lmom).M());

	    MyHistTools::fillTH("KP_MM_Kkm_MM_DEF", mm, (beam_lmom+target_lmom-km_lmom).M());
	    MyHistTools::fillTH("KP_KP_E_DEF", (beam_lmom+target_lmom-fc_lmom-km_lmom).E());

	    if( 0.45<mm && mm<0.54 ){
	      MyHistTools::fillTH("Kkm_MM_DEF_fp_mmKm", (beam_lmom+target_lmom-km_lmom).M());
	      MyHistTools::fillTH("KP_MM_DEF_km_1hit_true", mm);
	      MyHistTools::fillTH("mmKm_kin_DEF", mm_lmom.CosTheta(), mm_lmom.Vect().Mag());
	      MyHistTools::fillTH("Km_kin_DEF", km_lmom.CosTheta(), km_lmom.Vect().Mag());
	    }
	  }
	}
	if( anaInfo->nCDS(CDS_Proton)>0 ) MyHistTools::fillTH("KP_MM_DEF_p", mm);
      }
    }
  }
}


