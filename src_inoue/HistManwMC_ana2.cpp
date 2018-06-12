bool HistManwMC::ana(EventHeader *header)
{
  cdstrackMan-> Calc(cdsMan, confMan);

  TTree *tree;
  TH2F *h2;
  TH1F *h1;

  TH1F *h1_N  = (TH1F*)rtFile-> Get("N_Reduction");
  TH1F *h1_ER = (TH1F*)rtFile-> Get("EventReduction");
  if( fT0time==DBL_MIN ) return false;
  if( fD5mom==DBL_MIN ) return false;
  if( !fTrackBPC ) return false;

  for( int i=0; i<blMan->nBHD(); i++  ) if( blMan->BHD(i)->CheckRange()  ) fBHD_hit.push_back(blMan->BHD(i));
  for( int i=0; i<blMan->nT0(); i++   ) if( blMan->T0(i)->CheckRange()   ) fT0_hit.push_back(blMan->T0(i));
  for( int i=0; i<blMan->nBPD(); i++  ) if( blMan->BPD(i)->CheckRange()  ) fBPD_hit.push_back(blMan->BPD(i));
  for( int i=0; i<cdsMan->nCDH(); i++ ) if( cdsMan->CDH(i)->CheckRange() ) fCDH_hit.push_back(cdsMan->CDH(i));
  for( int i=0; i<blMan->nBVC(); i++  ) if( blMan->BVC(i)->CheckRange()  ) fBVC_hit.push_back(blMan->BVC(i));
  for( int i=0; i<blMan->nCVC(); i++  ) if( blMan->CVC(i)->CheckRange()  ) fCVC_hit.push_back(blMan->CVC(i));
  for( int i=0; i<blMan->nPC(); i++   ) if( blMan->PC(i)->CheckRange()   ) fPC_hit.push_back(blMan->PC(i));
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      int seg = blMan->NC(i)->seg();
      int lay = (seg-1)/16;
      fNC_hit[lay].push_back(blMan->NC(i));
      h1 = (TH1F*)rtFile-> Get("NC_hitpos"), h1-> Fill(blMan->NC(i)->hitpos());
    }
  }
  fT0pos = fTrackBPC->GetPosatZ(-110.5);

  double good_dis=DBL_MAX;
  bool good_vtx =false;
  //  std::cout<<"Event_Number:"<<header->ev()<<"  "<<cdstrackMan->nGoodTrack()<<std::endl;
  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *cdc = cdstrackMan->GoodTrack(i);
    double CDHtime, dis;
    double CDHdE=0;
    int CDHseg;
    if( !cdc-> GetCDHHit(cdsMan, CDHseg, CDHtime) ){
      //      std::cout<<"  !!! CDSTrack::GetCDHHit return false !!!"<<std::endl;
      continue;
    }
    for( int i=0; i<cdc->nCDHHit(); i++ ) CDHdE += cdc->CDHHit(cdsMan, i)->emean();

    double par[5];
    TVector3 vtxCDS, vtxBeam;
    cdc-> GetParameters(CID_CDC, par, vtxCDS);
    double mom = cdc-> Momentum();
    cdc-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxBeam, vtxCDS);
    //    TrackTools::CalcLineHelixVertex(fTrackBPC, cdc, vtxBeam, vtxCDS, dis);
    if( dis<good_dis ){
      good_vtx = true;
      good_dis=dis;
      fVtxBeam = vtxBeam;
      fVtxCDS = vtxCDS;
    }

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, vtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);

    TVector3 vtxCDH = cdc-> CDHVertex();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
    double beta = cdc_dis/(CDHtime-fT0time-beam_tof)/(100.*Const);
    double mass2 = mom*mom*(1./(beta*beta)-1);

    double calc_beta, tofvtxcdc;
    if( !TrackTools::FindMass2(cdc, fTrackBPC, CDHtime-fT0time, fD5mom, Beam_Kaon, calc_beta, mass2, tofvtxcdc) ){
      //      std::cout<<"  !!! return false TrackTools::FindMass2 !!!"<<std::endl;
      continue;
    }
    if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
      h1 = (TH1F*)rtFile-> Get("T0CDH_TOF"), h1-> Fill(CDHtime-fT0time);
    }
    int pid=TrackTools::PID(mom, mass2);
    //    int pid=MyTools::PID(mom, mass2);

    cdc-> SetPID(pid);
    if( MyTools::IsElectron(calc_beta, mom) ){
      if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
	h2 = (TH2F*)rtFile-> Get("CDS_overbeta_mom_e"), h2-> Fill(1./calc_beta, mom);
	h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_e"), h2-> Fill(mass2, mom);
      }
      if( mom<0 ) pid=98;
      else pid=99;
      //      pid=CDS_Other;
    }
    fCDSPID.push_back(pid);
    fCDSbeta.push_back(calc_beta);
    fCDSmass2.push_back(mass2);
    fCDSmom.push_back(mom);

    if( pid==CDS_PiMinus  ) fTrackPim.push_back(cdc);
    if( pid==CDS_Kaon     ) fTrackKm.push_back(cdc);
    if( pid==CDS_PiPlus   ) fTrackPip.push_back(cdc);
    if( pid==CDS_Proton   ) fTrackP.push_back(cdc);
    if( pid==CDS_Deuteron ) fTrackD.push_back(cdc);

    if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
      if( fBeamPID==Beam_Kaon ){
	h1 = (TH1F*)rtFile-> Get("CDS_chi2"), h1-> Fill(cdc->Chi());
	h2 = (TH2F*)rtFile-> Get("CDS_overbeta_mom"), h2-> Fill(1./calc_beta, mom);
	h2 = (TH2F*)rtFile-> Get("CDS_beta_overmom"), h2-> Fill(calc_beta, 1./mom);
	h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
	if( mom<0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_minus"), h1-> Fill(mass2);
	if( mom>0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_plus"), h1-> Fill(mass2);
      
	if( pid==CDS_PiMinus  ){
	  h1 = (TH1F*)rtFile-> Get("CDS_chi2_pim"), h1-> Fill(cdc->Chi());
	  h1 = (TH1F*)rtFile-> Get("CDS_mass2_pim"), h1-> Fill(mass2);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pim"), h2-> Fill(mass2, mom);
	}
	if( pid==CDS_Kaon     ){
	  h1 = (TH1F*)rtFile-> Get("CDS_chi2_K"), h1-> Fill(cdc->Chi());
	  h1 = (TH1F*)rtFile-> Get("CDS_mass2_K"), h1-> Fill(mass2);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_K"), h2-> Fill(mass2, mom);
	}
	if( pid==CDS_PiPlus   ){
	  h1 = (TH1F*)rtFile-> Get("CDS_chi2_pip"), h1-> Fill(cdc->Chi());
	  h1 = (TH1F*)rtFile-> Get("CDS_mass2_pip"), h1-> Fill(mass2);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pip"), h2-> Fill(mass2, mom);
	}
	if( pid==CDS_Proton   ){
	  h1 = (TH1F*)rtFile-> Get("CDS_chi2_p"), h1-> Fill(cdc->Chi());
	  h1 = (TH1F*)rtFile-> Get("CDS_mass2_p"), h1-> Fill(mass2);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_p"), h2-> Fill(mass2, mom);
	}
	if( pid==CDS_Deuteron ){
	  h1 = (TH1F*)rtFile-> Get("CDS_mass2_d"), h1-> Fill(mass2);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_d"), h2-> Fill(mass2, mom);
	}
      }
      else if( fBeamPID==Beam_Pion ){
	h1 = (TH1F*)rtFile-> Get("CDS_chi2_pi"), h1-> Fill(cdc->Chi());
	h2 = (TH2F*)rtFile-> Get("CDS_overbeta_mom_pi"), h2-> Fill(1./calc_beta, mom);
	h2 = (TH2F*)rtFile-> Get("CDS_beta_overmom_pi"), h2-> Fill(calc_beta, 1./mom);
	h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pi"), h2-> Fill(mass2, mom);
	if( mom<0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_minus_pi"), h1-> Fill(mass2);
	if( mom>0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_plus_pi"), h1-> Fill(mass2);
      }
    }

    TVector3 vtxCDS2, vtxBeam2;
    double tof, tmpl;
    if( !cdc->CalcVertexTimeLength(fT0pos, fTrackBPC->GetMomDir(), cdsMass[pid], vtxCDS2, vtxBeam2, tof, tmpl, true) ) continue;
    if( pid==CDS_PiMinus || pid==CDS_PiPlus || pid==CDS_Kaon || pid==CDS_Proton ){
      double offset = CDHtime-fT0time-beam_tof-tof;
      h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d", CDHseg)), h1-> Fill(offset);
    }

    if( simReader ){
      if( pid==CDS_PiMinus || pid==CDS_PiPlus || pid==CDS_Kaon || pid==CDS_Proton ){
	TVector3 vtx_b, vtx_cds;
	if( cdc-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_b, vtx_cds) ){
	  TVector3 p;
	  if( cdc-> GetMomentum(vtx_cds, p, true, true) ){
	    Track *track = simReader->trace(cdc, cdsMan);
	    double diff = p.Mag()-0.001*track->momentum().Mag();
	    //	    std::cout<<" ana mom : "<<p.Mag()<<" mc mom : "<<0.001*track->momentum().Mag()<<std::endl;
	    if( pid==CDS_PiPlus && track->pdgID()==211 ) h2 = (TH2F*)rtFile-> Get("CDS_mom_diff_pip_MC"), h2-> Fill(p.Mag(), diff);
	    if( pid==CDS_Proton && track->pdgID()==2212 )    h2 = (TH2F*)rtFile-> Get("CDS_mom_diff_p_MC"), h2-> Fill(p.Mag(), diff);
	    if( pid==CDS_PiMinus && track->pdgID()==-211 ) h2 = (TH2F*)rtFile-> Get("CDS_mom_diff_pim_MC"), h2-> Fill(p.Mag(), diff);
	    if( pid==CDS_Kaon && track->pdgID()==-321 )    h2 = (TH2F*)rtFile-> Get("CDS_mom_diff_km_MC"), h2-> Fill(p.Mag(), diff);
	  }
	}
      }
    }
  }

  if( good_vtx ){
    if( simReader ){
      TVector3 vtxMC = 0.1*simReader->getMCTrack(1)->vertex();
      double diff_x = vtxMC.X()-fVtxBeam.X();
      double diff_y = vtxMC.Y()-fVtxBeam.Y();
      double diff_z = vtxMC.Z()-fVtxBeam.Z();

      h2 = (TH2F*)rtFile-> Get("Vtx_XY_diff_MC"), h2-> Fill(diff_x, diff_y);
      h1 = (TH1F*)rtFile-> Get("Vtx_Z_diff_MC"), h1-> Fill(diff_z);
    }

    if( header ){
      if( fBeamPID==Beam_Kaon ){
	//	std::cout<<"  CDS 1track"<<std::endl;
	if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(7);
	h1_ER-> Fill(7);
      }
    }
    else{
      h1_ER-> Fill(4);
    }

    TVector3 beam_mom = fTrackBPC->GetMomDir();
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
    beam_mom.SetMag(beam_out);
    if( simReader ){
      ReactionData *reacData = simReader-> getReactionData();
      double beam_mom_MC = 0.001*reacData->GetInitParticle(0).Vect().Mag();
      h1 = (TH1F*)rtFile-> Get("BeamMom_diff_MC"), h1-> Fill(beam_mom_MC-beam_mom.Mag());
    }
    fBeamLmom.SetVectM(beam_mom, parMass[fBeamPID]);
    if( fBeamPID==Beam_Kaon ){
      h2 = (TH2F*)rtFile->Get("Vtx_XY"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
      h2 = (TH2F*)rtFile->Get("Vtx_ZX"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
      h2 = (TH2F*)rtFile->Get("Vtx_ZY"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	h2 = (TH2F*)rtFile->Get("Vtx_XY_wtar"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = (TH2F*)rtFile->Get("Vtx_ZX_wtar"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	h2 = (TH2F*)rtFile->Get("Vtx_ZY_wtar"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      }
      if( MyTools::IsTarget1(fVtxBeam, confMan) ){
	h2 = (TH2F*)rtFile->Get("Vtx_XY_wtar1"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = (TH2F*)rtFile->Get("Vtx_ZX_wtar1"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	h2 = (TH2F*)rtFile->Get("Vtx_ZY_wtar1"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      }
      if( MyTools::IsTarget2(fVtxBeam, confMan) ){
	h2 = (TH2F*)rtFile->Get("Vtx_XY_wtar2"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = (TH2F*)rtFile->Get("Vtx_ZX_wtar2"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	h2 = (TH2F*)rtFile->Get("Vtx_ZY_wtar2"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      }
      if( MyTools::IsTarget3(fVtxBeam, confMan) ){
	h2 = (TH2F*)rtFile->Get("Vtx_XY_wtar3"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = (TH2F*)rtFile->Get("Vtx_ZX_wtar3"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	h2 = (TH2F*)rtFile->Get("Vtx_ZY_wtar3"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      }
      if( MyTools::IsTarget4(fVtxBeam, confMan) ){
	h2 = (TH2F*)rtFile->Get("Vtx_XY_wtar4"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = (TH2F*)rtFile->Get("Vtx_ZX_wtar4"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	h2 = (TH2F*)rtFile->Get("Vtx_ZY_wtar4"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      }
      if( MyTools::IsTarget5(fVtxBeam, confMan) ){
	h2 = (TH2F*)rtFile->Get("Vtx_XY_wtar5"), h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
	h2 = (TH2F*)rtFile->Get("Vtx_ZX_wtar5"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.X());
	h2 = (TH2F*)rtFile->Get("Vtx_ZY_wtar5"), h2-> Fill(fVtxBeam.Z(), fVtxBeam.Y());
      }
    }
  }

  if( good_vtx && fCVC_hit.size()==0 && fBVC_hit.size()==0 && nNC()>0 ){
    //    if( !simReader ) anaNC(header);
    //    else anaNC_MC();
    anaNC(header);
  }
}

void HistManwMC::fill(EventHeader *header)
{
  rtFile-> cd();
  if( simReader ) simReader->fill(rtFile);

  TH1F *h1;
  TH2F *h2;
  TTree *tree;
  TNtuple *tup;
  if( header ){
    for( int i=0; i<20; i++ ){
