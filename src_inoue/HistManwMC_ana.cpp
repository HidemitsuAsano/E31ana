#include "HistManwMC.h"

using namespace std;

bool HistManwMC::ana(EventHeader *header)
{
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

  cdstrackMan-> Calc(cdsMan, confMan);
  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *cdc = cdstrackMan->GoodTrack(i);
    double CDHtime, dis;
    double CDHdE=0;
    int CDHseg; int CDH_clus_size;
    if( !MyTools::getCDHHit(cdsMan, cdstrackMan, cdc, CDHseg, CDH_clus_size, CDHtime) ){
      //    if( !cdc-> GetCDHHit(cdsMan, CDHseg, CDHtime) ){
      //      std::cout<<"  !!! CDSTrack::GetCDHHit return false !!!"<<std::endl;
      continue;
    }


    h1 = (TH1F*)rtFile-> Get("CDH_clus_size"), h1-> Fill(cdc->nCDHHit());

    for( int i=0; i<cdc->nCDHHit(); i++ ){
      CDHdE += cdc->CDHHit(cdsMan, i)->emean();
      for( int j=i+1; j<cdc->nCDHHit(); j++ ){
	HodoscopeLikeHit *hit1 = cdc->CDHHit(cdsMan, i);
	HodoscopeLikeHit *hit2 = cdc->CDHHit(cdsMan, j);

	double time1 = hit1->ctmean(), time2=hit2->ctmean();
	h1 = (TH1F*)rtFile-> Get("CDH_time_diff"), h1-> Fill(time2-time1);
	if( simReader ){
	  if( simReader->isSameHit(hit1, hit2) ) h1 = (TH1F*)rtFile-> Get("CDH_time_diff_true"), h1-> Fill(time2-time1);
	  else	h1 = (TH1F*)rtFile-> Get("CDH_time_diff_false"), h1-> Fill(time2-time1);
	}
      }
    }

    double par[5];
    TVector3 vtxCDS, vtxBeam;
    cdc-> GetParameters(CID_CDC, par, vtxCDS);
    double mom = cdc-> Momentum();
    if( !cdc-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxBeam, vtxCDS) ) continue;

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, vtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);

    TVector3 vtxCDH = cdc-> CDHVertex();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxCDS);
    double beta = cdc_dis/(CDHtime-fT0time-beam_tof)/(100.*Const);
    double mass2 = mom*mom*(1./(beta*beta)-1);

    double calc_beta, tofvtxcdc;
    if( !TrackTools::FindMass2(cdc, fTrackBPC, CDHtime-fT0time, fD5mom, Beam_Kaon, calc_beta, mass2, tofvtxcdc) ){
      continue;
    }
    //*** Add 2015/11/25 ******************************//
    if( !simReader ){
      cdc-> Retiming(cdsMan, confMan, calc_beta, true);
      //      for(int j=0;j<3;j++) cdc->HelixFitting(cdsMan);
    }
    //*************************************************//
  }
  cdstrackMan-> Calc(cdsMan, confMan);
  h2 = (TH2F*)rtFile->Get("nCDC");
  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    h2-> Fill(cdsMan->nCDC(lay), lay);
  }

  double good_dis=DBL_MAX;
  bool good_vtx =false;

  // if( simReader ){
  //   cout<<"===== CDS Check ====="<<endl;
  //   for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
  //     cout<<"> lay : "<<lay<<" nCDC : "<<cdsMan->nCDC(lay)<<endl;
  //   }
  //   cout<<"> CDSTrackMan   nGoodTrack : "<<cdstrackMan->nGoodTrack()<<endl;
  // }

  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *cdc = cdstrackMan->GoodTrack(i);
    double CDHtime, dis;
    double CDHdE=0;
    int CDHseg;
    int CDH_clus_size;
    if( !MyTools::getCDHHit(cdsMan, cdstrackMan, cdc, CDHseg, CDH_clus_size, CDHtime) ){
      //    if( !cdc-> GetCDHHit(cdsMan, CDHseg, CDHtime) ){
      //      std::cout<<"  !!! CDSTrack::GetCDHHit return false !!!"<<std::endl;
      continue;
    }
    for( int i=0; i<cdc->nCDHHit(); i++ ) CDHdE += cdc->CDHHit(cdsMan, i)->emean();

    double par[5];
    TVector3 vtxCDS, vtxBeam;
    cdc-> GetParameters(CID_CDC, par, vtxCDS);
    double mom = cdc-> Momentum();
    if( !cdc-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxBeam, vtxCDS) ) continue;
    dis = (vtxBeam-vtxCDS).Mag();
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
    bool is_share_CDH = MyTools::isShareCDH(cdsMan, cdstrackMan, i);
    // if( is_share_CDH ){
    //   int tmp_seg, tmp_size; double tmp_time;
    //   MyTools::getCDHHit(confMan, cdsMan, cdstrackMan, cdc, tmp_seg, tmp_size, tmp_time);
    // }

    if( !TrackTools::FindMass2(cdc, fTrackBPC, CDHtime-fT0time, fD5mom, Beam_Kaon, calc_beta, mass2, tofvtxcdc) ){
      //      std::cout<<"  !!! return false TrackTools::FindMass2 !!!"<<std::endl;
      continue;
    }

    if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
      h1 = (TH1F*)rtFile-> Get("T0CDH_TOF"), h1-> Fill(CDHtime-fT0time);
    }
    int pid=MyTools::PIDcorr_wide(mom, mass2); // 2D mass2 vs mom
    //    int pid=MyTools::PIDcorr_wide(mom, mass2); // 2D mass2 vs mom
    //    int pid=MyTools::PID(mom, mass2); // mass2 w/o electron-like cut
    int pid2=TrackTools::PID(mom, mass2);
    if( is_share_CDH ) pid=CDS_Other;
    //    if( cdc->nCDHHit()>1 ) pid=CDS_Other;

    // if( pid==CDS_Other ) cdc-> SetPID(pid2);
    // else cdc-> SetPID(pid);
    cdc->SetPID(pid);

    if( simReader ){
      //      if( (pid2==CDS_PiMinus && pid!=CDS_PiMinus) || (pid2==CDS_PiPlus && pid!=CDS_PiPlus ) ){

      //      }
       
      h2 = (TH2F*)rtFile-> Get("nCDC_tail");
      for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
	h2-> Fill(cdsMan->nCDC(lay), lay);
      }
    }

    if( MyTools::IsElectron(calc_beta, mom) ){
      if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
	h2 = (TH2F*)rtFile-> Get("CDS_overbeta_mom_e"), h2-> Fill(1./calc_beta, mom);
	h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_e"), h2-> Fill(mass2, mom);
      }
      if( mom<0 ) pid=98;
      else pid=99;
      //      pid=CDS_Other;
    }

    h1 = (TH1F*)rtFile-> Get("CDS_chi2"), h1-> Fill(cdc->Chi());

    if( cdc->Chi()<30. ){ // off tempolary
      fCDSPID.push_back(pid);
      fCDSbeta.push_back(calc_beta);
      fCDSmass2.push_back(mass2);
      fCDSmom.push_back(mom);

      fillCDC(cdc, calc_beta, cdsMan, confMan);

      for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
	for( int i=0; i<cdc->nTrackHit(lay); i++ ){
	  CDCHit *hit = cdc->hit(cdsMan, lay, i);
	  h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d", lay)), h1-> Fill(hit->resl());
	  if( cdc->PID()==CDS_PiMinus ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_pim", lay)), h1-> Fill(hit->resl());
	  if( cdc->PID()==CDS_PiPlus  ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_pip", lay)), h1-> Fill(hit->resl());
	  if( cdc->PID()==CDS_Kaon    ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_km", lay)), h1-> Fill(hit->resl());
	  if( cdc->PID()==CDS_Proton  ) h1 = (TH1F*)rtFile-> Get(Form("CDC_res_%d_p", lay)), h1-> Fill(hit->resl());
	}
      }
    
      if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
	if( fBeamPID==Beam_Kaon ){
	  if( is_share_CDH ) h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_CDHshare"), h2-> Fill(mass2, mom);
	  if( cdc->nCDHHit()>1 ) h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_CDH2"), h2-> Fill(mass2, mom);
	
	  h2 = (TH2F*)rtFile-> Get("CDS_overbeta_mom"), h2-> Fill(1./calc_beta, mom);
	  h2 = (TH2F*)rtFile-> Get("CDS_beta_overmom"), h2-> Fill(calc_beta, 1./mom);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
	  h2 = (TH2F*)rtFile-> Get(Form("CDS_mass2_mom_CDH_seg%d", CDHseg)), h2->Fill(mass2, mom);

	  if( simReader ){
	    //	    cout<<" CDHseg="<<CDHseg<<endl;
	    DetectorHit *det_hit = simReader-> getHit(CID_CDH, CDHseg);
	    int CDScheck = simReader-> checkCDS(cdc, cdstrackMan, cdsMan, mass2, vtxCDS, vtxCDH);
	    //	    std::cout<<"CDH Hit mass2:"<<mass2<<"  mom:"<<mom<<std::endl;
	    //	    std::cout<<"    "<<simReader->parName(det_hit->pdg())<<std::endl;

	    if( det_hit->pdg()==11 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_e"), h2-> Fill(mass2, mom);
	      if( CDScheck==2 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_e_2"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_e_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_e_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else if( det_hit->pdg()==-11 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_ep"), h2-> Fill(mass2, mom);
	      if( CDScheck==3 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_ep_3"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_ep_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_ep_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else if( det_hit->pdg()==13 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mum"), h2-> Fill(mass2, mom);
	      if( CDScheck==10 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mum_10"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==11 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mum_11"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mum_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mum_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else if( det_hit->pdg()==-13 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mup"), h2-> Fill(mass2, mom);
	      if( CDScheck==15 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mup_15"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==16 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mup_16"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mup_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_mup_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else if( det_hit->pdg()==211 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip"), h2-> Fill(mass2, mom);
	      if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip_0"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==30 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip_30"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==31 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip_31"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==32 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip_32"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pip_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else if( det_hit->pdg()==-211 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim"), h2-> Fill(mass2, mom);
	      if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim_0"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==20 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim_20"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==21 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim_21"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==22 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim_22"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_pim_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else if( det_hit->pdg()==-321 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_km"), h2-> Fill(mass2, mom);
	      if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_km_0"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_km_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_km_2tra"), h2-> Fill(mass2, mom); }
	    }
	    else if( det_hit->pdg()==321 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_kp"), h2-> Fill(mass2, mom);
	    }
	    else if( det_hit->pdg()==2212 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p"), h2-> Fill(mass2, mom);
	      if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_0"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==50 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_50"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==51 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_51"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==52 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_52"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==53 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_53"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_2hit"), h2-> Fill(mass2, mom); }
	      else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_p_2tra"), h2-> Fill(mass2, mom); }
	      //	      else simReader->checkCDS(cdc, cdsMan, true);
	    }
	    else{
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_MC_other"), h2-> Fill(mass2, mom);
	    }
	  }

	  if( mom<0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_minus"), h1-> Fill(mass2);
	  if( mom>0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_plus"), h1-> Fill(mass2);
	  
	  if( pid==CDS_PiMinus  ){
	    h1 = (TH1F*)rtFile-> Get("CDS_mass2_pim"), h1-> Fill(mass2);
	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pim"), h2-> Fill(mass2, mom);
	  }
	  if( pid==CDS_Kaon     ){
	    h1 = (TH1F*)rtFile-> Get("CDS_mass2_K"), h1-> Fill(mass2);
	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_K"), h2-> Fill(mass2, mom);
	  }
	  if( pid==CDS_PiPlus   ){
	    h1 = (TH1F*)rtFile-> Get("CDS_mass2_pip"), h1-> Fill(mass2);
	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pip"), h2-> Fill(mass2, mom);
	  }
	  if( pid==CDS_Proton   ){
	    h1 = (TH1F*)rtFile-> Get("CDS_mass2_p"), h1-> Fill(mass2);
	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_p"), h2-> Fill(mass2, mom);
	  }
	  if( pid==CDS_Deuteron ){
	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_d"), h2-> Fill(mass2, mom);
	  }
	}
	else if( fBeamPID==Beam_Pion ){
	  h2 = (TH2F*)rtFile-> Get("CDS_overbeta_mom_pi"), h2-> Fill(1./calc_beta, mom);
	  h2 = (TH2F*)rtFile-> Get("CDS_beta_overmom_pi"), h2-> Fill(calc_beta, 1./mom);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pi"), h2-> Fill(mass2, mom);
	  if( mom<0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_minus_pi"), h1-> Fill(mass2);
	  if( mom>0 ) h1 = (TH1F*)rtFile-> Get("CDS_mass2_plus_pi"), h1-> Fill(mass2);
	}
      }
    
      TVector3 vtxCDS2, vtxBeam2;
      double tof, tmpl;
      if( !cdc->CalcVertexTimeLength(fT0pos, fTrackBPC->GetMomDir(), cdsMass[pid], vtxCDS2, vtxBeam2, tof, tmpl, true) ){
	//	std::cout<<" !!! Calc Vertex Length return false !!!"<<std::endl; 
	continue;
      }
      //      else{ std::cout<<" !!! Calc Vertex Length return true !!!"<<std::endl; }

      if( pid==CDS_PiMinus  ) fTrackPim.push_back(cdc);
      if( pid==CDS_Kaon     ) fTrackKm.push_back(cdc);
      if( pid==CDS_PiPlus   ) fTrackPip.push_back(cdc);
      if( pid==CDS_Proton   ) fTrackP.push_back(cdc);
      if( pid==CDS_Deuteron ) fTrackD.push_back(cdc);

      // if( !simReader ){
      // 	TVector3 vtxCDH = cdc->CDHVertex();
      // 	double param[5];
      // 	cdc-> GetParameters(param);
      // 	double fl = MathTools::CalcHelixArc(param, vtxCDS2, vtxCDH);
      // 	double mom = cdc->Momentum();
      // 	if( mom>0.3 ){
      // 	  double beta = fabs(mom/(piMass*sqrt(1+mom*mom/(piMass*piMass))));
      // 	  double calc_tof = fl/(beta*Const*100.);
      // 	  double offset = CDHtime-fT0time-calc_tof-beam_tof;
      // 	  h1 = (TH1F*)rtFile-> Get(Form("CDH_offset0_%d", CDHseg)), h1-> Fill(offset);
      // 	}
      // }
    
      if( pid==CDS_PiMinus || pid==CDS_PiPlus || pid==CDS_Kaon || pid==CDS_Proton ){
	TVector3 vtxCFRP;
	double param[5];
	cdc-> GetParameters(CID_CDCCFRP, param, vtxCFRP);
	//	std::cout<<" CDCCFRP r="<<vtxCFRP.Perp()<<std::endl;
	TVector3 vtxCDH = cdc->CDHVertex();
	cdc-> GetParameters(param);
	double fl = MathTools::CalcHelixArc(param, vtxCFRP, vtxCDH);
	double mom = fabs(cdc-> Momentum(CID_CDCCFRP));
	double mass = cdsMass[pid];
	double momout;
	double tmptof;
	ELossTools::CalcdE(mom, mass, fl, "CDCGas", momout, -1, tmptof);
	double offset = CDHtime-fT0time-beam_tof-tof-tmptof;
	//	MyTools::printCDSParams(cdc);
	  //	std::cout<<"Mat : "<<mat<<" mom CFRP : "<<mom<< " mom CDH : "<<momout<<std::endl;
	  //	std::cout<<"CDH offset : "<<offset<<" beam_tof="<<beam_tof<<"  Vtx_CDCin tof="<<tof<<"  CDCin CDH tof="<<tmptof<<std::endl;
	HodoscopeLikeHit *hit=0;
	for( int i=0; i<cdc->nCDHHit(); i++ ){
	  HodoscopeLikeHit *hit2 = cdc->CDHHit(cdsMan, i);
	  if( hit2 && hit2->seg()==CDHseg ) hit=hit2;
	}
	double eu = hit->eu();
	double ed = hit->ed();
	double tu = hit->tu()-fT0time-beam_tof-tof-tmptof;
	double td = hit->td()-fT0time-beam_tof-tof-tmptof;
	double ctu = hit->ctu()-fT0time-beam_tof-tof-tmptof;
	double ctd = hit->ctd()-fT0time-beam_tof-tof-tmptof;
	double tsub = hit->ctsub();

	if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
	  if( !simReader ){	
	    // TNtuple *tup = (TNtuple*)rtFile-> Get(Form("CDH%d_slewing", CDHseg));
	    // tup-> Fill(eu, ed, tu, td, ctu, ctd, offset);
	  }
	  
	  h2 = (TH2F*)rtFile-> Get(Form("CDH_Z_tsub_%d", CDHseg)), h2-> Fill(vtxCDH.Z(), tsub);
	  
	  h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_T0%d", fT0_hit[0]->seg())), h1-> Fill(offset);
	  h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d_T0%d", CDHseg, fT0_hit[0]->seg())), h1-> Fill(offset);

	  h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d", CDHseg)), h1-> Fill(offset);
	  h2 = (TH2F*)rtFile-> Get("CDH_offset_mom"), h2-> Fill(offset, mom);
	  if( pid==CDS_PiMinus ){
	    h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d_pim", CDHseg)), h1-> Fill(offset);
	    h2 = (TH2F*)rtFile-> Get("CDH_offset_mom_pim"), h2-> Fill(offset, mom);
	  }
	  if( pid==CDS_PiPlus  ){
	    h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d_pip", CDHseg)), h1-> Fill(offset);
	    h2 = (TH2F*)rtFile-> Get("CDH_offset_mom_pip"), h2-> Fill(offset, mom);
	  }
	  if( pid==CDS_Kaon    ){
	    h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d_km", CDHseg)), h1-> Fill(offset);
	    h2 = (TH2F*)rtFile-> Get("CDH_offset_mom_km"), h2-> Fill(offset, mom);
	  }
	  if( pid==CDS_Proton  ){
	    h1 = (TH1F*)rtFile-> Get(Form("CDH_offset_%d_p", CDHseg)), h1-> Fill(offset);
	    h2 = (TH2F*)rtFile-> Get("CDH_offset_mom_p"), h2-> Fill(offset, mom);
	  }
	
	  double IH_time;
	  int IH_seg;
	  if( cdc->GetIHHit(cdsMan, IH_seg, IH_time) ){
	    TVector3 vtxIH = cdc-> IHVertex();
	    cdc-> GetParameters(param);
	    double fl = MathTools::CalcHelixArc(param, vtxCFRP, vtxIH);
	    mom = fabs(cdc-> Momentum(CID_CDCCFRP));
	    ELossTools::CalcdE(mom, mass, fl, "CDCGas", momout, -1, tmptof);
	    HodoscopeLikeHit *IHhit=0;
	    for( int i=0; i<cdsMan->nIH(); i++ ){
	      HodoscopeLikeHit *hit2 = cdsMan->IH(i);
	      if( hit2->seg()==IH_seg ) IHhit=hit2;
	    }
	    int adc = IHhit-> adcu();
	    double dE = IHhit-> eu();
	    double IHoffset = IH_time-fT0time-beam_tof-tof-tmptof;
	    
	    if( mom>0.3 && (pid==CDS_PiMinus || pid==CDS_PiPlus) ){
	      h1 = (TH1F*)rtFile-> Get(Form("IH_ADC_%d_cut", IH_seg)), h1-> Fill(adc);
	      h1 = (TH1F*)rtFile-> Get(Form("IH_dE_%d_cut", IH_seg)), h1-> Fill(dE);
	      h2 = (TH2F*)rtFile-> Get(Form("IH_Z_dE_%d_cut", IH_seg)), h2->Fill(vtxIH.Z(), dE);
	    }
	    //	    std::cout<<" IH offset seg"<<IH_seg<<" time : "<<IH_time<<"  offset : "<<IHoffset<<std::endl;
	    h2 = (TH2F*)rtFile-> Get(Form("IH_dE_offset_%d", IH_seg)), h2-> Fill(dE, offset);
	    h1 = (TH1F*)rtFile-> Get(Form("IH_offset_%d", IH_seg)), h1-> Fill(IHoffset);
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
      }
    }// if{ cdc->chi2()<30) }
    else{
      if( GeomTools::GetID(vtxBeam)==CID_Fiducial ){
    	h2 = (TH2F*)rtFile->Get("CDS_mass2_mom_chi30"), h2-> Fill(mass2, mom);
      }
    }
  }

  h1 = (TH1F*)rtFile-> Get("CDC_nGoodTrack_wo"), h1-> Fill(fCDSPID.size());
  if( good_vtx ){
    h1 = (TH1F*)rtFile-> Get("CDC_nGoodTrack"), h1-> Fill(fCDSPID.size());
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
      h1_ER-> Fill(7);
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
      if( -0.1<fVtxBeam.Z()  && fVtxBeam.z()<0.1 ){
	h2 = (TH2F*)rtFile-> Get("Vtx_XY_z0"),h2-> Fill(fVtxBeam.X(), fVtxBeam.Y());
      }
      if( fVtxBeam.Perp()<2.0 ){
	h1 = (TH1F*)rtFile-> Get("Vtx_Z_r2"), h1-> Fill(fVtxBeam.Z());
      }
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
  //  anaFC(header);

  if( cdstrackMan->nGoodTrack()==3 && fCDSPID.size()==3 ){
    CDSTrack *track0 = cdstrackMan-> GoodTrack(0);
    CDSTrack *track1 = cdstrackMan-> GoodTrack(1);
    CDSTrack *track2 = cdstrackMan-> GoodTrack(2);
    int nm=0;
    int p_id=-1;
    if( track0->Momentum()>0 ) nm++, p_id=0;
    if( track1->Momentum()>0 ) nm++, p_id=1;
    if( track2->Momentum()>0 ) nm++, p_id=2;

    if( nm==2 ){
      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	h2 = (TH2F*)rtFile->Get("CDS_mass2_mom_1p2m");
	h2-> Fill(fCDSmass2[0], fCDSmom[0]); h2-> Fill(fCDSmass2[1], fCDSmom[1]); h2-> Fill(fCDSmass2[2], fCDSmom[2]);
      }
    }
  }

  if( cdstrackMan->nGoodTrack()==2 && fCDSPID.size()==2 ){
    CDSTrack *track0 = cdstrackMan-> GoodTrack(0);
    CDSTrack *track1 = cdstrackMan-> GoodTrack(1);
    if( ( track0->Momentum()<0 && track1->Momentum()>0 ) || ( track0->Momentum()>0 && track1->Momentum()<0 ) ){
      if( fFPID==F_Neutron ){
	// if( simReader ){
	//   std::cout<<"====== CDS +- charge & neutron ====="<<std::endl;
	//   simReader-> trace(CID_CDH, seg0, true);
	//   simReader-> trace(CID_CDH, seg1, true);
	// }

	if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	  h2 = (TH2F*)rtFile-> Get("hitpatNC");
	  for( int lay=0; lay<7; lay++ ){
	    for( int i=0; i<fNC_hit[lay].size(); i++ ){
	      int seg2 = fNC_hit[lay][i]-> seg()%16;
	      if( seg2==0 ) seg2=16;
	      //	      cout<<" lay : "<<lay<<"  "<<seg2<<endl;
	      h2-> Fill(seg2, lay+1);
	    }
	  }

	  double mass2_plus=fCDSmass2[0], mass2_minus=fCDSmass2[1], mom_plus=fCDSmom[0], mom_minus=fCDSmom[1];
	  if( mom_plus<0 ){
	    std::swap(mom_minus, mom_plus);
	    std::swap(mass2_minus, mass2_plus);
	  }
	  int plus_id=0, minus_id=1;
	  if( track0-> Momentum()<0 ){
	    std::swap(track0, track1);
	    plus_id=1, minus_id=0;
	  }

	  if( MyTools::isShareCDH(cdsMan, cdstrackMan, plus_id) ) h2=(TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_CDHshare"),h2-> Fill(mass2_plus,  mom_plus);
	  if( MyTools::isShareCDH(cdsMan, cdstrackMan, minus_id) ) h2=(TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_CDHshare"),h2-> Fill(mass2_minus, mom_minus);
	  if( track0->nCDHHit()>1 ) h2=(TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_CDH2"),h2-> Fill(mass2_plus,  mom_plus);
	  if( track1->nCDHHit()>1 ) h2=(TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_CDH2"),h2-> Fill(mass2_minus, mom_minus);

	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm"), h2-> Fill(fCDSmass2[0], fCDSmom[0]);
	  h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm"), h2-> Fill(fCDSmass2[1], fCDSmom[1]);
	  h2 = (TH2F*)rtFile-> Get("nCDC_wn_pm");
	  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){ h2-> Fill(cdsMan->nCDC(lay), lay); };
	
	  //	  cout<<" mom+ : "<<track0->Momentum()<<" "<<mom_plus<<endl;
	  //	  cout<<" mom- : "<<track1->Momentum()<<" "<<mom_minus<<endl;
	
	  int seg0, seg1;
	  double time0, time1;
	  bool flag0=track0-> GetCDHHit(cdsMan, seg0, time0);
	  bool flag1=track1-> GetCDHHit(cdsMan, seg1, time1);
	  if( flag0 && flag1 ){
	  // if( seg0==0 ) cout<<" CDH seg0 : "<<seg0<<endl;
	  // if( seg1==0 ) cout<<" CDH seg1 : "<<seg1<<endl;
	
	    int pid_plus=MyTools::PIDcorr2(mom_plus, mass2_plus);
	    int pid2_plus=TrackTools::PID(mom_plus, mass2_plus);
	    int pid_minus =MyTools::PIDcorr2(mom_minus, mass2_minus);
	    int pid2_minus=TrackTools::PID(mom_minus, mass2_minus);
	    if( (pid2_plus==CDS_PiMinus && pid_plus!=CDS_PiMinus) || (pid2_minus==CDS_PiMinus && pid_minus!=CDS_PiMinus) ){
	      h2 = (TH2F*)rtFile-> Get("nCDC_wn_pm_tail");
	      for( int lay=1; lay<=NumOfCDCLayers; lay++ ){ h2-> Fill(cdsMan->nCDC(lay), lay); };
	    }
	    
	    double kn_mm = (fBeamLmom+D_LMOM-fFLmom).M();
	    if( fTrackPip.size()>0 ) h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_pip"), h2-> Fill(mass2_minus, mom_minus);
	    if( fTrackP.size()>0 ) h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_p"), h2-> Fill(mass2_minus, mom_minus);
	    if( fTrackKm.size()>0  ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_km"), h2-> Fill(mass2_plus, mom_plus);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_km_2track"), h1-> Fill(kn_mm);
	    }
	 
	    TVector3 beam_plus_vtx, beam_minus_vtx, plus_beam_vtx, minus_beam_vtx, plus_vtx, minus_vtx;
	    TVector3 beam_plsu_mom, beam_minus_mom, plus_mom, minus_mom;
	    bool plus_vtx_flag  = track0-> GetVertex(fT0pos, fTrackBPC-> GetMomDir(), beam_plus_vtx, plus_beam_vtx);
	    bool minus_vtx_flag = track1-> GetVertex(fT0pos, fTrackBPC-> GetMomDir(), minus_beam_vtx, beam_minus_vtx);
	    bool CDS_vtx_flag = TrackTools::Calc2HelixVertex(track0, track1, plus_vtx, minus_vtx);
	    bool plus_mom_flag = track0-> GetMomentum(plus_beam_vtx, plus_mom, true, true);
	    bool minus_mom_flag = track0-> GetMomentum(minus_beam_vtx, minus_mom, true, true);
	    TVector3 plus_mm_mom = fBeamLmom.Vect()-fFLmom.Vect()-plus_mom;
	    TVector3 minus_mm_mom = fBeamLmom.Vect()-fFLmom.Vect()-minus_mom;
	    TVector3 plus_mm_vtx, minus_mm_vtx;
	    TVector3 mm_vtx_plus_l, mm_vtx_minus_h, mm_vtx_minus_l, mm_vtx_plus_h;
	    bool plus_mm_vtx_flag = track1-> GetVertex(beam_plus_vtx, plus_mm_mom.Unit(), mm_vtx_plus_l,  mm_vtx_minus_h);
	    bool minus_mm_vtx_flag = track0-> GetVertex(beam_minus_vtx, minus_mm_mom.Unit(), mm_vtx_minus_l,  mm_vtx_plus_h);
	    
	    double dca_bp = (beam_plus_vtx-plus_beam_vtx).Mag();
	    double dca_bm = (beam_minus_vtx-minus_beam_vtx).Mag();

	    double dca_mm_lp = (mm_vtx_plus_l-mm_vtx_plus_h).Mag();
	    double dca_mm_lm = (mm_vtx_plus_l-mm_vtx_plus_h).Mag();

	    h2 = (TH2F*)rtFile-> Get("CDS_DCA_wn_pm"), h2-> Fill(dca_bp, dca_bm);
	    h2 = (TH2F*)rtFile-> Get("CDS_mmDCA_wn_pm"), h2-> Fill(dca_mm_lp, dca_mm_lm);

	    int chk_plus  = MyTools::checkArea(mom_plus, mass2_plus);
	    int chk_minus = MyTools::checkArea(mom_minus, mass2_minus);
	    if( chk_minus==-1 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_DCA_tailm"), h2-> Fill(dca_bp, dca_bm);
	      h2 = (TH2F*)rtFile-> Get("CDS_mmDCA_tailm"), h2-> Fill(dca_mm_lp, dca_mm_lm);
	    }
	    if( chk_plus==1   ){
	      h2 = (TH2F*)rtFile-> Get("CDS_DCA_tailp"), h2-> Fill(dca_bp, dca_bm);
	      h2 = (TH2F*)rtFile-> Get("CDS_mmDCA_tailp"), h2-> Fill(dca_mm_lp, dca_mm_lm);
	    }
	    if( fTrackPim.size()==1 && fTrackPip.size()==1 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_DCA_npipi"), h2-> Fill(dca_bp, dca_bm);
	      h2 = (TH2F*)rtFile-> Get("CDS_mmDCA_npipi"), h2-> Fill(dca_mm_lp, dca_mm_lm);
	    }
	    
	    
	    if( simReader ){
	      DetectorHit *hit0 = simReader->getHit(CID_CDH, seg0);
	      DetectorHit *hit1 = simReader->getHit(CID_CDH, seg1);
	      if( track0->Momentum()<0 ){ std::swap(hit0, hit1); }
	    // std::cout<<"===== CDS +- charge & forward n ====="<<std::endl;
	    // std::cout<<"> + par ; "<<simReader->parName(hit0->pdg())<<std::endl;
	    // std::cout<<"> - par ; "<<simReader->parName(hit1->pdg())<<std::endl;
	    // if( hit0->pdg()!=211 ){ simReader->printTrackTrace(hit0->trackID()); };
	    // if( hit1->pdg()!=-211  ){ simReader->printTrackTrace(hit1->trackID()); };
	    //	    simReader->printTrackTrace(hit0->trackID());
	    //	    simReader->printTrackTrace(hit1->trackID());
	      TVector3 vtxBeam, vtxCDS;
	      
	      track0-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxBeam, vtxCDS);
	      int CDScheck = simReader-> checkCDS(track0, cdstrackMan, cdsMan, mass2_plus, vtxBeam, vtxCDS);
	      if( hit0->pdg()==11 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==2 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e_2"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==-11 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==3 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep_3"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==13 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==10 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_10"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==11 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_11"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==-13 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==15 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_15"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==16 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_16"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==-211 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_0"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==20 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_20"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==21 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_21"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==22 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_22"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==211 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_0"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==30 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_30"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==31 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_31"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==32 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_32"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==-321 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km_0"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else if( hit0->pdg()==321 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_kp"), h2-> Fill(mass2_plus, mom_plus);
	      }
	      else if( hit0->pdg()==2212 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p"), h2-> Fill(mass2_plus, mom_plus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_0"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==50 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_50"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==51 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_51"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==52 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_52"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==53 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_53"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_2hit"), h2-> Fill(mass2_plus, mom_plus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_2tra"), h2-> Fill(mass2_plus, mom_plus); }
	      }
	      else{
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_other"), h2-> Fill(mass2_plus, mom_plus);
	      }

	      track1-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxBeam, vtxCDS);						
	      CDScheck = simReader-> checkCDS(track1, cdstrackMan, cdsMan, mass2_minus, vtxBeam, vtxCDS);
	      if( hit1->pdg()==11 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==2 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e_2"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_e_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==-11 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==3 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep_3"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_ep_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==13 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==10 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_10"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==11 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_11"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mum_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==-13 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==15 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_15"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==16 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_16"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_mup_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==-211 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_0"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==20 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_20"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==21 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_21"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==22 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_22"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pim_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==211 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_0"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==30 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_30"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==31 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_31"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==32 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_32"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_pip_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==-321 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km_0"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_km_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else if( hit1->pdg()==321 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_kp"), h2-> Fill(mass2_minus, mom_minus);
	      }
	      else if( hit1->pdg()==2212 ){
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p"), h2-> Fill(mass2_minus, mom_minus);
		if( CDScheck==0 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_0"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==50 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_50"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==51 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_51"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==52 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_52"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==53 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_53"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==100 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_2hit"), h2-> Fill(mass2_minus, mom_minus); }
		else if( CDScheck==200 ){ h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_p_2tra"), h2-> Fill(mass2_minus, mom_minus); }
	      }
	      else{
		h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_MC_other"), h2-> Fill(mass2_minus, mom_minus);
	      }
	    }
	  }
	}
      }
    }
  }
  if( fFPID==F_Neutron ){
    if( cdstrackMan->nGoodTrack()>2 ){
      bool pim_flag=false, pip_flag=false;
      for( int i=0; i<fCDSmass2.size(); i++ ){
	if( TrackTools::PID(fCDSmom[i], fCDSmass2[i])==CDS_PiMinus ) pim_flag=true;
	if( TrackTools::PID(fCDSmom[i], fCDSmass2[i])==CDS_PiPlus ) pip_flag=true;
      }

      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_wn_pm_3");
      for( int i=0; i<fCDSmom.size(); i++ ){
	h2-> Fill(fCDSmass2[i], fCDSmom[i]);
      }
    }
  }
}
