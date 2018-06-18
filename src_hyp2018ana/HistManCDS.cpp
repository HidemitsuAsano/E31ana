#include "HistManBeamAna.h"

void HistManBeamAna::initCDS()
{
  std::cout<<"CDS Histogram Initialize ... "<<std::endl;
  std::cout<<"  for CDH"<<std::endl;
  new TH1F("hitpat_CDH", "CDH hit pattern", 36, 0.5, 36.5);
  new TH1F("mul_CDH", "CDH multiplicity", 37, -0.5, 36.5);
  for( int seg=1; seg<=36; seg++ ){
    new TH1F(Form("CDH_offset_%d",seg), Form("CDH offset seg%d", seg), 3000, -50, 250);

    new TH1F(Form("CDH_ADC_up_%d", seg),   Form("CDH ADC up seg%d", seg),          4000, 0.0, 4000);
    new TH1F(Form("CDH_ADC_down_%d", seg), Form("CDH ADC down seg%d", seg),        4000, 0.0, 4000);
    new TH1F(Form("CDH_AwT_up_%d", seg),   Form("CDH ADC up seg%d w/ TDC", seg),   4000, 0.0, 4000);
    new TH1F(Form("CDH_AwT_down_%d", seg), Form("CDH ADC down seg%d w/ TDC", seg), 4000, 0.0, 4000);

    new TH1F(Form("CDH_dE_up_%d", seg),     Form("CDH dE up seg%d", seg),          500, -1.0, 29);
    new TH1F(Form("CDH_dE_down_%d", seg),   Form("CDH dE down seg%d", seg),        500, -1.0, 29);
    new TH1F(Form("CDH_dEwT_up_%d", seg),   Form("CDH dE up seg%d w/ TDC", seg),   500, -1.0, 29);
    new TH1F(Form("CDH_dEwT_down_%d", seg), Form("CDH dE down seg%d w/ TDC", seg), 500, -1.0, 29);   

    new TH1F(Form("CDH_ADC_up_%d_cut", seg),   Form("CDH ADC up seg%d", seg),          4000, 0.0, 4000);
    new TH1F(Form("CDH_ADC_down_%d_cut", seg), Form("CDH ADC down seg%d", seg),        4000, 0.0, 4000);
    new TH1F(Form("CDH_AwT_up_%d_cut", seg),   Form("CDH ADC up seg%d w/ TDC", seg),   4000, 0.0, 4000);
    new TH1F(Form("CDH_AwT_down_%d_cut", seg), Form("CDH ADC down seg%d w/ TDC", seg), 4000, 0.0, 4000);

    new TH1F(Form("CDH_dE_up_%d_cut", seg),     Form("CDH dE up seg%d", seg),          500, -1.0, 29);
    new TH1F(Form("CDH_dE_down_%d_cut", seg),   Form("CDH dE down seg%d", seg),        500, -1.0, 29);
    new TH1F(Form("CDH_dEwT_up_%d_cut", seg),   Form("CDH dE up seg%d w/ TDC", seg),   500, -1.0, 29);
    new TH1F(Form("CDH_dEwT_down_%d_cut", seg), Form("CDH dE down seg%d w/ TDC", seg), 500, -1.0, 29);   
  }
#if CDC
  std::cout<<"  for CDC"<<std::endl;
  for(int lay=1; lay<=NumOfCDCLayers; lay++ ){
    new TH1F(Form("CDC_resi_%d", lay),    Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH1F(Form("CDC_resi_pi_%d", lay), Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH1F(Form("CDC_resi_p_%d", lay),  Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH1F(Form("CDC_resi_k_%d", lay),  Form("CDC residual layer%d", lay), 200, -0.1, 0.1);
    new TH2F(Form("CDC_dt_resi_%d", lay), Form("CDC dt vs residual layer%d", lay), 300, -20, 280, 200, -0.1, 0.1);
    for( int wire=1; wire<=NumOfCDCWiresInLayer[lay-1]; wire++ ){
      new TH1F(Form("CDC_dt_%d_%d", lay, wire),      Form("CDC dt l%d w%d", lay, wire), 600, -200, 400);
      new TH1F(Form("CDC_resi_%d_%d", lay, wire),    Form("CDC residual l%d w%d", lay, wire), 200, -0.1, 0.1);
      new TH1F(Form("CDC_resi_pi_%d_%d", lay, wire), Form("CDC residual l%d %d", lay, wire), 200, -0.1, 0.1);
      new TH1F(Form("CDC_resi_p_%d_%d", lay, wire),  Form("CDC residual l%d w%d", lay, wire), 200, -0.1, 0.1);
      new TH1F(Form("CDC_resi_k_%d_%d", lay, wire),  Form("CDC residual l%d w%d", lay, wire), 200, -0.1, 0.1);
      new TH2F(Form("CDC_dt_resi_%d_%d", lay, wire), Form("CDC dt vs residual l%d w%d", lay, wire), 300, -20, 280,  200, -0.1, 0.1);
      new TH2F(Form("CDC_dt_dl_%d_%d", lay, wire)  , Form("CDC dt vs dl l%d w%d", lay, wire),       300, -20, 280, 1000, -0.1, 0.9);
    }
  }
#endif
  std::cout<<"  for CDS"<<std::endl;
  new TH1F("CDS_nGoodTrack", "n Good Track", 10, -0.5, 9.5);

  new TH1F("CDS_FindMass_status", "CDS Find Mass Status",  12, -0.5, 11.5);
  new TH2F("CDS_mass2_mom", "CDS mass^{2} vs mom", 3000, -0.1, 5.9, 600, -1.5, 1.5); 
  new TH2F("CDS_mass2_mom_wtar", "CDS mass^{2} vs mom", 3000, -0.1, 5.9, 600, -1.5, 1.5); 
  new TH1F("CDS_mass2_plus",  "CDS mass2 plus", 3000, -0.5, 5.5);
  new TH1F("CDS_mass2_minus", "CDS mass2 plus", 3000, -0.5, 5.5);

  for( int seg=1; seg<=36; seg++ ){
    new TH2F(Form("CDS_mass2_mom_CDH%d", seg), "CDS mass^{2} vs mom", 3000, -0.1, 5.9, 600, -1.5, 1.5); 
    new TH2F(Form("CDS_mass2_mom_wtar_CDH%d", seg), "CDS mass^{2} vs mom", 3000, -0.1, 5.9, 600, -1.5, 1.5); 
  }

  for( int seg=1; seg<=5; seg++ ){
    new TH2F(Form("CDS_mass2_mom_T0%d", seg), "CDS mass^{2} vs mom", 3000, -0.1, 5.9, 600, -1.5, 1.5); 
    new TH2F(Form("CDS_mass2_mom_wtar_T0%d", seg), "CDS mass^{2} vs mom", 3000, -0.1, 5.9, 600, -1.5, 1.5); 
  }

  new TH2F("Vtx_XY_beam", "Vertex XY Plane",  600, -15, 15, 600, -15, 15);
  new TH2F("Vtx_ZX_beam", "Vertex ZX Plane", 1200, -30, 30, 600, -15, 15);
  new TH2F("Vtx_ZY_beam", "Vertex ZY Plane", 1200, -30, 30, 600, -15, 15);

  new TH2F("Vtx_XY_CDS", "Vertex XY Plane",  600, -15, 15, 600, -15, 15);
  new TH2F("Vtx_ZX_CDS", "Vertex ZX Plane", 1200, -30, 30, 600, -15, 15);
  new TH2F("Vtx_ZY_CDS", "Vertex ZY Plane", 1200, -30, 30, 600, -15, 15);

  new TH1F("Vtx_dis", "Vertex distance", 1000, 0.0, 100);

  new TH2F("Vtx_XY_beam_w", "Vertex XY Plane",  600, -15, 15, 600, -15, 15);
  new TH2F("Vtx_ZX_beam_w", "Vertex ZX Plane", 1200, -30, 30, 600, -15, 15);
  new TH2F("Vtx_ZY_beam_w", "Vertex ZY Plane", 1200, -30, 30, 600, -15, 15);

  new TH2F("Vtx_XY_CDS_w", "Vertex XY Plane",  600, -15, 15, 600, -15, 15);
  new TH2F("Vtx_ZX_CDS_w", "Vertex ZX Plane", 1200, -30, 30, 600, -15, 15);
  new TH2F("Vtx_ZY_CDS_w", "Vertex ZY Plane", 1200, -30, 30, 600, -15, 15);


  new TH2F("Vtx_XY_wtar", "Vertex XY Plane",  600, -15, 15, 600, -15, 15);
  new TH2F("Vtx_ZX_wtar", "Vertex ZX Plane", 1200, -30, 30, 600, -15, 15);
  new TH2F("Vtx_ZY_wtar", "Vertex ZY Plane", 1200, -30, 30, 600, -15, 15);

  std::cout<<"  for CDS IM"<<std::endl;
  new TH1F("CDS_IM_ppim", "CDS p #pi^{-} IM",       1000, 1.0, 2.0);
  new TH1F("CDS_IM_pipi", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);

  new TH1F("CDS_IM_ppim_wtar", "CDS p #pi^{-} IM",       1000, 1.0, 2.0);
  new TH1F("CDS_IM_pipi_wtar", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);

  std::cout<<"CDS Histogram Initialize ... finish"<<std::endl;
}

void HistManBeamAna::fillCDS(ConfMan *conf, EventHeader *header, CDSHitMan *cdsMan, BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan, BeamLineTrackMan *bltrackMan)
{
  TH1 *h1;
  TH2 *h2;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan-> CDH(i);
    int seg=hit->seg();
    int adcu=hit-> adcu(), adcd=hit->adcd();
    double eu=hit-> eu(), ed=hit->ed();

    h1 = getTH(Form("CDH_ADC_up_%d", seg)), h1-> Fill(adcu);
    h1 = getTH(Form("CDH_ADC_down_%d", seg)), h1-> Fill(adcd);
    h1 = getTH(Form("CDH_dE_up_%d", seg)), h1-> Fill(eu);
    h1 = getTH(Form("CDH_dE_down_%d", seg)), h1-> Fill(ed);
    if( hit-> CheckRange() ){
      fCDH_hit.push_back(hit);
      h1 = getTH("hitpat_CDH"), h1->Fill(hit->seg());

      h1 = getTH(Form("CDH_AwT_up_%d", seg)), h1-> Fill(adcu);
      h1 = getTH(Form("CDH_AwT_down_%d", seg)), h1-> Fill(adcd);
      h1 = getTH(Form("CDH_dEwT_up_%d", seg)), h1-> Fill(eu);
      h1 = getTH(Form("CDH_dEwT_down_%d", seg)), h1-> Fill(ed);
    }
  }
  h1 = getTH("mul_CDH"), h1-> Fill(fCDH_hit.size());

#if CDC
  for( int i=0;i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *track = cdstrackMan-> GoodTrack(i);
    bool single_flag = true;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      if( track-> nTrackHit(layer)!=1 ) single_flag = false;
    }
    double chi;
    double param[5];
    track-> GetParameters(param);
    double drho=param[0], phi0=param[1], rho=1./param[2], dz=param[3], tlam=param[4];
    double mom=track-> Momentum();
    if( single_flag && fabs(param[3])<12 && fabs(drho)<6 && rho>0 ){
      for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
	CDCHit *cdc = track-> TrackHit(cdsMan, lay, 0);
	double resi = cdc-> resl();
	double dt = cdc-> dt();
	double dl = cdc->dl();
	double dlr = dl-resi;
	int wire = cdc-> wire();
	h1 = getTH(Form("CDC_resi_%d", lay)), h1-> Fill(resi);
	h1 = getTH(Form("CDC_resi_%d_%d", lay, wire)), h1-> Fill(resi);
	h1 = getTH(Form("CDC_dt_%d_%d", lay, wire)), h1-> Fill(dt);
	h2 = getTH2(Form("CDC_dt_resi_%d", lay)), h2-> Fill(dt, resi);
	h2 = getTH2(Form("CDC_dt_resi_%d_%d", lay, wire)), h2-> Fill(dt, resi);
	h2 = getTH2(Form("CDC_dt_dl_%d_%d", lay, wire)), h2-> Fill(dt, dlr);
      }
    }
  }
#endif
  if( (fBeamPID==Beam_Kaon || fBeamPID==Beam_Pion) && fBeamMom>0.001 ){
    if( bltrackMan->ntrackBLC2()==1 && bltrackMan->ntrackBPC()==1 && fT0_hit.size()==1 ){
      int T0seg = fT0_hit[0]-> seg();
      LocalTrack *trackBLC2 = bltrackMan->trackBLC2(0);
      LocalTrack *trackBPC = bltrackMan->trackBPC(0);
      TVector3 gpos;
      conf-> GetGeomMapManager()-> GetGPos(CID_T0, T0seg, gpos);
      TVector3 T0pos = trackBLC2-> GetPosatZ(gpos.Z()-0.5);

      h1 = getTH("CDS_nGoodTrack"), h1-> Fill(cdstrackMan->nGoodTrack());   
      for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
	CDSTrack* cdsTrack = cdstrackMan->GoodTrack(i);
	double param[5];
	cdsTrack-> GetParameters(param);

	if( cdsTrack-> nCDHHit()==1 ){
	  HodoscopeLikeHit *CDHhit = cdsTrack-> CDHHit(cdsMan, 0);
	  int CDHseg = CDHhit-> seg();
	  TVector3 CDHvertex = cdsTrack-> CDHVertex();
	  TVector3 vtxb, vtxh;
	  double adcu = CDHhit-> adcu();
	  double adcd = CDHhit-> adcd();
	  double eu = CDHhit-> eu();
	  double ed = CDHhit-> ed();
	  double tmpdis;
	  TrackTools::CalcLineHelixVertex(trackBPC, cdsTrack, vtxb, vtxh, tmpdis);
	  double cdc_dis = MathTools::CalcHelixArc(param, vtxh, CDHvertex);
	  double mom = cdsTrack-> Momentum();
	  double v = 100.*Const*fabs(mom)/piMass/sqrt(1.0+mom*mom/(piMass*piMass));
	  double cdc_calc_tof = cdc_dis/v;
	  double beam_tof, out_mom;
	  ELossTools::CalcElossBeamTGeo(T0pos, vtxb, fBeamMom, parMass[fBeamPID], out_mom, beam_tof);
	  double offset = CDHhit->ctmean()-fT0_hit[0]->ctmean()-cdc_calc_tof-beam_tof;
	  h1 = getTH(Form("CDH_offset_%d", CDHseg)), h1-> Fill(offset);

	  if( -0.5<offset && offset<0.5 && fabs(mom)>0.3 ){
	    h1 = getTH(Form("CDH_ADC_up_%d_cut", CDHseg)), h1-> Fill(adcu);
	    h1 = getTH(Form("CDH_ADC_down_%d_cut", CDHseg)), h1-> Fill(adcd);
	    h1 = getTH(Form("CDH_dE_up_%d_cut", CDHseg)), h1-> Fill(eu);
	    h1 = getTH(Form("CDH_dE_down_%d_cut", CDHseg)), h1-> Fill(ed);
	  }
	}

	if( cdsTrack->CDHFlag() ){
	  double CDHtime = 99999;

	  int CDHseg = -1;
	  for( int i=0; i<cdsTrack->nCDHHit(); i++ ){
	    if( cdsTrack-> CDHHit(cdsMan,i)-> ctmean() <CDHtime ){
	      CDHtime = cdsTrack->CDHHit(cdsMan, i)->ctmean();
	      CDHseg = cdsTrack->CDHHit(cdsMan, i)->seg();
	    }
	  }

	  double tof = CDHtime-fT0_hit[0]->ctmean();
	  double beta_calc, mass2, tmp_tof;

       	  if( TrackTools::FindMass2C(cdsTrack, trackBPC, tof, fBeamMom, fBeamPID, beta_calc, mass2, tmp_tof) ){
	    double mom = cdsTrack->Momentum();

	    if( mom>0 ){
	      h1 = getTH("CDS_mass2_plus"), h1-> Fill(mass2);
	    }
	    else{
	      h1 = getTH("CDS_mass2_minus"), h1-> Fill(mass2);
	    }
	    h2 = getTH2("CDS_mass2_mom"), h2-> Fill(mass2, mom);
	    h2 = getTH2(Form("CDS_mass2_mom_CDH%d", CDHseg)), h2-> Fill(mass2, mom);
	    h2 = getTH2(Form("CDS_mass2_mom_T0%d", T0seg)), h2-> Fill(mass2, mom);
	    
	    int pid = TrackTools::PID(mom, mass2);
	    cdsTrack-> SetPID(pid);
	    if( pid==CDS_PiMinus       ) fCDSpim.push_back(cdsTrack); 
	    else if( pid==CDS_Kaon     ) fCDSkm.push_back(cdsTrack); 
	    else if( pid==CDS_PiPlus   ) fCDSpip.push_back(cdsTrack); 
	    else if( pid==CDS_Proton   ) fCDSp.push_back(cdsTrack); 
	    else if( pid==CDS_Deuteron ) fCDSd.push_back(cdsTrack); 

	    TVector3 vtx_beam, vtx_cds;
	    double tofvtxcdc, tmpl;
	    TrackTools::CalcLineHelixVertex(trackBPC, cdsTrack, vtx_beam, vtx_cds, tmpl);
	    bool near_flag = false;
	    if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	      h2 = getTH2("CDS_mass2_mom_wtar"), h2-> Fill(mass2, mom);
	      h2 = getTH2(Form("CDS_mass2_mom_wtar_CDH%d", CDHseg)), h2-> Fill(mass2, mom);
	      h2 = getTH2(Form("CDS_mass2_mom_wtar_T0%d", T0seg)), h2-> Fill(mass2, mom);
	    }
	    if( (vtx_beam-vtx_cds).Mag()<fVtxDis ){
	      fVtxDis = (vtx_beam-vtx_cds).Mag();
	      fVtxBeam=vtx_beam;
	      fVtxCDS=vtx_cds;
	      near_flag = true;
	    }
	    if( cdsTrack->CalcVertexTimeLength(T0pos, trackBPC->GetMomDir(), cdsMass[pid], vtx_beam, vtx_cds, tofvtxcdc, tmpl, true) ){
	      fVtxBeam_w = vtx_beam;
	      fVtxCDS_w = vtx_cds;
	    }
	  }
	}
      }
      h1 = getTH("Vtx_dis"), h1-> Fill(fVtxDis);
      if( fVtxDis<10){
	h2 = getTH2("Vtx_XY_beam"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = getTH2("Vtx_ZY_beam"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
	h2 = getTH2("Vtx_ZX_beam"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());

	h2 = getTH2("Vtx_XY_CDS"), h2-> Fill(fVtxCDS.X(), fVtxCDS.Y());
	h2 = getTH2("Vtx_ZY_CDS"), h2-> Fill(fVtxCDS.Z(), fVtxCDS.Y());
	h2 = getTH2("Vtx_ZX_CDS"), h2-> Fill(fVtxCDS.Z(), fVtxCDS.X());

	h2 = getTH2("Vtx_XY_beam_w"), h2-> Fill(fVtxBeam_w.X(), fVtxBeam_w.Y());
	h2 = getTH2("Vtx_ZY_beam_w"), h2-> Fill(fVtxBeam_w.Z(), fVtxBeam_w.Y());
	h2 = getTH2("Vtx_ZX_beam_w"), h2-> Fill(fVtxBeam_w.Z(), fVtxBeam_w.X());

	h2 = getTH2("Vtx_XY_CDS_w"), h2-> Fill(fVtxCDS_w.X(), fVtxCDS_w.Y());
	h2 = getTH2("Vtx_ZY_CDS_w"), h2-> Fill(fVtxCDS_w.Z(), fVtxCDS_w.Y());
	h2 = getTH2("Vtx_ZX_CDS_w"), h2-> Fill(fVtxCDS_w.Z(), fVtxCDS_w.X());

	if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	  h2 = getTH2("Vtx_XY_wtar"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	  h2 = getTH2("Vtx_ZY_wtar"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
	  h2 = getTH2("Vtx_ZX_wtar"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	}
      }    

      for( int i=0; i<fCDSpim.size(); i++ ){
	for( int j=0; j<fCDSp.size(); j++ ){
	  TVector3 vtx1, vtx2;
	  if( TrackTools::Calc2HelixVertex(fCDSpim[i], fCDSp[j], vtx1, vtx2) ){
	    TVector3 pim_mom, p_mom;
	    if( fCDSpim[i]-> GetMomentum(vtx1, pim_mom, true, true) ){
	      if( fCDSp[j]-> GetMomentum(vtx2, p_mom, true, true) ){
		TLorentzVector pim_lmom, p_lmom;
		pim_lmom.SetVectM(pim_mom, piMass);
		p_lmom.SetVectM(p_mom, pMass);
		TVector3 vtx_mean = 0.5*(vtx1+vtx2);
		double dist, dltmp=0;
		TVector3 vtxb, vtxp;
		MathTools::LineToLine( vtx_mean, (pim_mom+p_mom).Unit(), T0pos, trackBPC->GetMomDir(), dltmp, dist, vtxb, vtxp);

		double im = (pim_lmom+p_lmom).M();
		h1 = getTH("CDS_IM_ppim"), h1-> Fill(im);
		if( L_MIN<im && im<L_MAX ) fLflag = true;
		if( MyTools::IsTarget(vtxb) ){
		  h1 = getTH("CDS_IM_ppim_wtar"), h1-> Fill(im);
		}
	      }
	    }
	  }
	}
      }

      for( int i=0; i<fCDSpim.size(); i++ ){
	for( int j=0; j<fCDSpip.size(); j++ ){
	  TVector3 vtx1, vtx2;
	  if( TrackTools::Calc2HelixVertex(fCDSpim[i], fCDSpip[j], vtx1, vtx2) ){
	    TVector3 pim_mom, pip_mom;
	    if( fCDSpim[i]-> GetMomentum(vtx1, pim_mom, true, true) ){
	      if( fCDSpip[j]-> GetMomentum(vtx2, pip_mom, true, true) ){
		TLorentzVector pim_lmom, pip_lmom;
		pim_lmom.SetVectM(pim_mom, piMass);
		pip_lmom.SetVectM(pip_mom, piMass);
		TVector3 vtx_mean = 0.5*(vtx1+vtx2);
		double dist, dltmp=0;
		TVector3 vtxb, vtxp;
		MathTools::LineToLine( vtx_mean, (pim_mom+pip_mom).Unit(), T0pos, trackBPC->GetMomDir(), dltmp, dist, vtxb, vtxp);

		double im = (pim_lmom+pip_lmom).M();
		if( K0_MIN<im && im<K0_MAX ) fK0flag = true;
		h1 = getTH("CDS_IM_pipi"), h1-> Fill(im);
		if( MyTools::IsTarget(vtxb) ){
		  h1 = getTH("CDS_IM_pipi_wtar"), h1-> Fill(im);
		}
	      }
	    }
	  }
	}

      }
    }
  }
}


