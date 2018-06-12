#include "HistManwMC.h"
#include "MyParam.h"

static const TLorentzVector TGT_LMOM = TLorentzVector(0., 0., 0., pMass);
static const TLorentzVector P_LMOM = TLorentzVector(0., 0., 0., pMass);
static const double NC_THRE = 8.0;

HistManwMC::HistManwMC(TFile *f, ConfMan *conf)
  : rtFile(f), confMan(conf), runHeaderMC(0), evHeaderMC(0), mcData(0), detData(0), reacData(0), blMan(0), cdsMan(0), bltrackMan(0), cdstrackMan(0)
{
  Clear();
}

bool HistManwMC::ana(EventHeader *header)
{
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
      int lay = (seg-1)/15;
      fNC_hit[lay].push_back(blMan->NC(i));
      h1 = (TH1F*)rtFile-> Get("NC_hitpos"), h1-> Fill(blMan->NC(i)->hitpos());
    }
  }
  fT0pos = fTrackBPC->GetPosatZ(-110.5);

  double good_dis=DBL_MAX;
  bool good_vtx =false;

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
    }

    TVector3 vtxCDS2, vtxBeam2;
    double tof, tmpl;
    if( !cdc->CalcVertexTimeLength(fT0pos, fTrackBPC->GetMomDir(), cdsMass[pid], vtxCDS2, vtxBeam2, tof, tmpl, true) ) continue;
  }

  if( good_vtx ){
    if( header ){
      if( fBeamPID==Beam_Kaon ){
	if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(7);
	h1_ER-> Fill(7);
      }
    } 
    TVector3 beam_mom = fTrackBPC->GetMomDir();
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
    beam_mom.SetMag(beam_out);
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
    bool NC_flag = false;
    int nlay=0;
    for( int lay=0; lay<8; lay++ ){
      for( int i=0; i<fNC_hit[lay].size(); i++ ){
	if( fNC_hit[lay][i]->emean()>NC_THRE && fNC_hit[lay][i]->emean()>fNCdE ){ 
	  NC_flag =true;
	  fNCseg  = fNC_hit[lay][i]->seg();
	  fNCtime = fNC_hit[lay][i]->ctmean();
	  fNCdE   = fNC_hit[lay][i]->emean();
	  confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNC_hit[lay][i]->seg() , fNCpos);

	  if( !mcData ) fNCpos.SetY(fNC_hit[lay][i]->hitpos());
	  //	  fNCpos.SetY(fNC_hit[lay][i]->hitpos());
	}
      }
      if( NC_flag ) break;
    }

    if( fNCdE>0 ){
      double fl = (fNCpos-fVtxBeam).Mag();
      double beam_out, beam_tof;
      ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      fNCbeta = fl/((fNCtime-fT0time-beam_tof)*100.*Const);

      if( fBeamPID==Beam_Kaon ){
	h1 = (TH1F*)rtFile-> Get("NC_overbeta"), h1-> Fill(1./fNCbeta);
	h2 = (TH2F*)rtFile-> Get("NC_overbeta_dE"), h2-> Fill(1./fNCbeta, fNCdE);
	if( fNCdE>NC_THRE  ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee"), h1-> Fill(1./fNCbeta);

	if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	  h1 = (TH1F*)rtFile-> Get("NC_overbeta_wtar"), h1-> Fill(1./fNCbeta);
	  if( fNCdE>NC_THRE) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_wtar"), h1-> Fill(1./fNCbeta);
	  h2 = (TH2F*)rtFile-> Get("NC_overbeta_dE_wtar"), h2-> Fill(1./fNCbeta, fNCdE);
	}
      }
    }

    if( NC_flag ){
      if( header && fBeamPID==Beam_Kaon ){
	h1_ER-> Fill(8);
	if( header->IsTrig(Trig_Neutral) ){
	  h1_N-> Fill(8);
	  h1_ER-> Fill(9);	
	}
      }
    
      if( fNCbeta>0.9 ){
	fFPID=F_Gamma;
	if( header ){
	  if( header->IsTrig(Trig_Neutral) && fBeamPID==Beam_Kaon ) h1_N-> Fill(9);
	}
      }
      else{
	fFPID=F_Neutron;
	double NCmom = nMass*fNCbeta/sqrt(1-fNCbeta*fNCbeta);
	TVector3 n_mom = fNCpos-fVtxBeam;
	n_mom.SetMag(NCmom);
	fFLmom.SetVectM(n_mom, nMass);
	if( header ){
	  if( header->IsTrig(Trig_Neutral) && fBeamPID==Beam_Kaon){
	    h1_N-> Fill(10);
	    h1_ER->Fill(10);
	  }
	}
      }
    }
  }

  // get MC Data
  if( mcData ) anaMC();
}

void HistManwMC::fill(EventHeader *header)
{
  rtFile-> cd();
  TH1F *h1;
  TH2F *h2;
  TNtuple *tup;
  TH1F *h1_N  = (TH1F*)rtFile-> Get("N_Reduction");
  TH1F *h1_ER  = (TH1F*)rtFile->Get("EventReduction");
  bool trigNC_flag = false;
  if( header ){
    if( header->IsTrig(Trig_Neutral) ) trigNC_flag=true;
  }
  else{
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit = detData-> detectorHit(i);
      if( hit->detectorID()==CID_NC ) trigNC_flag=true;
    }
  }
  
  //  std::cout<<"===== HistManwMC::fill START  ====="<<std::endl;
  //******************//
  //*** for CDS IM ***//
  //******************//
  for( int i=0; i<fTrackPim.size(); i++ ){
    //*** CDS pi+ pi- ***//
    for( int j=0; j<fTrackPip.size(); j++ ){
      TVector3 vtx_pim, vtx_pip;
      if( TrackTools::Calc2HelixVertex(fTrackPim[i], fTrackPip[j], vtx_pim, vtx_pip) ){
	TVector3 pim_mom, pip_mom;
	if( fTrackPim[i]-> GetMomentum(vtx_pim, pim_mom, true, true) &&
	    fTrackPip[j]-> GetMomentum(vtx_pip, pip_mom, true, true) ){
	  TLorentzVector pim_lmom, pip_lmom;
	  pim_lmom.SetVectM(pim_mom, piMass);
	  pip_lmom.SetVectM(pip_mom, piMass);
	  TVector3 vtx_mean = 0.5*(vtx_pim+vtx_pip);
	  TVector3 mom_sum = pim_mom+pip_mom;

	  TVector3 vtxb, vtxcds;
	  double dltmp=0;
	  double dca;
	  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtxcds, vtxb);

	  double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(fT0pos, vtxb, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	  TVector3 beam_mom = fTrackBPC->GetMomDir();
	  beam_mom.SetMag(beam_out);
	  TLorentzVector beam_lmom;
	  beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);

	  if( GeomTools::GetID(vtxb)==CID_Fiducial ){
	    if( fBeamPID==Beam_Kaon ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi"), h1-> Fill((pim_lmom+pip_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KCDSpipi_MM"), h1-> Fill((beam_lmom+TGT_LMOM-pim_lmom-pip_lmom).M());
	    }
	  }
	}
      }
    }
  
    //*** CDS p pi- ***//
    for( int j=0; j<fTrackP.size(); j++ ){
      TVector3 vtx_pim, vtx_p;
      if( TrackTools::Calc2HelixVertex(fTrackPim[i], fTrackP[j], vtx_pim, vtx_p) ){
	TVector3 pim_mom, p_mom;
	if( fTrackPim[i]-> GetMomentum(vtx_pim, pim_mom, true, true) &&
	    fTrackP[j]-> GetMomentum(vtx_p, p_mom, true, true) ){
	  TLorentzVector pim_lmom, p_lmom;
	  pim_lmom.SetVectM(pim_mom, piMass);
	  p_lmom.SetVectM(p_mom, pMass);
	  TVector3 vtx_mean = 0.5*(vtx_pim+vtx_p);
	  TVector3 mom_sum = pim_mom+p_mom;

	  TVector3 vtxb, vtxcds;
	  double dltmp=0;
	  double dca;
	  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtxcds, vtxb);

	  double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(fT0pos, vtxb, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	  TVector3 beam_mom = fTrackBPC->GetMomDir();
	  beam_mom.SetMag(beam_out);
	  TLorentzVector beam_lmom;
	  beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);

	  if( GeomTools::GetID(vtxb)==CID_Fiducial ){
	    if( fBeamPID==Beam_Kaon ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_ppim"), h1-> Fill((pim_lmom+p_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KCDSppim_MM"), h1-> Fill((beam_lmom+TGT_LMOM-pim_lmom-p_lmom).M());
	    }
	  }
	}
      }
    }
  }

  //*** CDS p K- ***//
  for( int i=0; i<fTrackP.size(); i++ ){
    for( int j=0; j<fTrackKm.size(); j++ ){
      TVector3 vtx_p, vtx_km;
      if( TrackTools::Calc2HelixVertex(fTrackP[i], fTrackKm[j], vtx_p, vtx_km) ){
	TVector3 p_mom, km_mom;
	if( fTrackP[i]-> GetMomentum(vtx_p, p_mom, true, true) &&
	    fTrackKm[j]-> GetMomentum(vtx_km, km_mom, true, true) ){
	  TLorentzVector p_lmom, km_lmom;
	  p_lmom.SetVectM(p_mom, pMass);
	  km_lmom.SetVectM(km_mom, kpMass);
	  TVector3 vtx_mean = 0.5*(vtx_p+vtx_km);
	  TVector3 mom_sum = p_mom+km_mom;

	  TVector3 vtxb, vtxcds;
	  double dltmp=0;
	  double dca;
	  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtxcds, vtxb);

	  double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(fT0pos, vtxb, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	  TVector3 beam_mom = fTrackBPC->GetMomDir();
	  beam_mom.SetMag(beam_out);
	  TLorentzVector beam_lmom;
	  beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);

          if( fTrackP.size()==1 && fTrackKm.size()==1 ){
            h1 = (TH1F*)rtFile-> Get("DCA_pkm"), h1-> Fill((vtx_p-vtx_km).Mag());
            if( (vtx_p-vtx_km).Mag()<1.0 ){
              TVector3 vtx_mean = 0.5*(vtx_p+vtx_km);
              TVector3 BPC_calc_pos = fTrackBPC->GetPosatZ(vtx_mean.Z());
              TVector3 diff = vtx_mean-BPC_calc_pos;
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_dx_dy"), h2-> Fill(diff.X(), diff.Y());
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_z_dx"),  h2-> Fill(vtx_mean.Z(), diff.X());
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_z_dy"),  h2-> Fill(vtx_mean.Z(), diff.Y());
            }
          }

	  if( GeomTools::GetID(vtxb)==CID_Fiducial ){
	    if( fBeamPID==Beam_Kaon ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_pkm"), h1-> Fill((p_lmom+km_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KCDSpkm_MM"), h1-> Fill((beam_lmom+TGT_LMOM-p_lmom-km_lmom).M());
	    }
	  }
	}
      }
    }
  }

  //**********************************//
  //*** MC K- d -> L(1405) n study ***//
  //**********************************//
  if( fFNFlag && fD5mom!=DBL_MIN && fT0time!=DBL_MIN && fTrackBPC!=0 ){
    if( mcData ){
      if( reacData->ReactionID()==1098 || reacData-> ReactionID()==3000 ){
	if( fBeamPID==Beam_Kaon ){
	  h1 = (TH1F*)rtFile->Get("L1405_mass_MC_wNChit"), h1->Fill(fYstarLmom.M());
	  if( fDecayMode==-1 ) h1 = (TH1F*)rtFile-> Get("L1405_mass_MC_wNChit_Sm"), h1-> Fill(fYstarLmom.M());
	  if( fDecayMode==0  ) h1 = (TH1F*)rtFile-> Get("L1405_mass_MC_wNChit_S0"), h1-> Fill(fYstarLmom.M());
	  if( fDecayMode==1  ) h1 = (TH1F*)rtFile-> Get("L1405_mass_MC_wNChit_Sp"), h1-> Fill(fYstarLmom.M());
	}
      }
    }
  }

  //************************//
  //*** Forward Neutron  ***//
  //************************//
  if( fFPID==F_Neutron && trigNC_flag ){
    TLorentzVector kn_lmom = fBeamLmom+TGT_LMOM-fFLmom;
    if( mcData ){
      //*** for resolution estimation by MC ***//
      if( fFNFlag && (reacData->ReactionID()==1098 || reacData->ReactionID()==3000 ) ){
	if( fBeamPID==Beam_Kaon ){
	  h1 = (TH1F*)rtFile->Get("n_mom_res"), h1->Fill(1000.*(fFLmom.Vect().Mag()-fNLmom1.Vect().Mag()));
	  h1 = (TH1F*)rtFile->Get("KN_MM_res"), h1->Fill(1000.*(kn_lmom.M()-fYstarLmom.M()));
	  h2 = (TH2F*)rtFile->Get("L1405_KN_MM_res"), h2->Fill(fYstarLmom.M(), 1000.*(kn_lmom.M()-fYstarLmom.M()));
	}
      }
    }

    if( fBeamPID==Beam_Kaon ){
      h1 = (TH1F*)rtFile-> Get("KN_MM_wo"), h1->Fill(kn_lmom.M());
      h1 = (TH1F*)rtFile-> Get("KN_MM_wo_B"), h1->Fill(kn_lmom.M());
      if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim_wo"), h1->Fill(kn_lmom.M());
      if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip_wo"), h1->Fill(kn_lmom.M());
      if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km_wo"), h1->Fill(kn_lmom.M());
      if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p_wo"), h1->Fill(kn_lmom.M());
      if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim_wo_B"), h1->Fill(kn_lmom.M());
      if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip_wo_B"), h1->Fill(kn_lmom.M());
      if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km_wo_B"), h1->Fill(kn_lmom.M());
      if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p_wo_B"), h1->Fill(kn_lmom.M());

      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	if( header ){
	  if( header->IsTrig(Trig_Neutral) ){
	    h1_N-> Fill(11);
	  }
	  h1_ER->Fill(11);
	}
	
	h1 = (TH1F*)rtFile-> Get("KN_MM"), h1->Fill(kn_lmom.M());
	h1 = (TH1F*)rtFile-> Get("KN_MM_B"), h1->Fill(kn_lmom.M());
	if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim"), h1->Fill(kn_lmom.M());
	if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip"), h1->Fill(kn_lmom.M());
	if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km"), h1->Fill(kn_lmom.M());
	if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p"), h1->Fill(kn_lmom.M());
	if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim_B"), h1->Fill(kn_lmom.M());
	if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip_B"), h1->Fill(kn_lmom.M());
	if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km_B"), h1->Fill(kn_lmom.M());
	if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p_B"), h1->Fill(kn_lmom.M());
      }
    }

    //***************************//
    //*** d(K-, n) w/ pi+ pi- ***//
    //***************************//
    if( fTrackPim.size()==1 && fTrackPip.size()==1 ){  
      if( header && fBeamPID==Beam_Kaon ){
	h1_ER-> Fill(12);
	if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(12);
      }
      double pim_mass2;
      double pip_mass2;
      double pim_mom2;
      double pip_mom2;
      for( int i=0; i<fCDSPID.size(); i++ ){
	if( fCDSPID[i]==CDS_PiMinus ){
	  pim_mass2 = fCDSmass2[i];
	  pim_mom2 = fCDSmom[i];
	}
	else if( fCDSPID[i]==CDS_PiPlus ){
	  pip_mass2 = fCDSmass2[i];
	  pip_mom2 = fCDSmom[i];
	}
      }


      TVector3 vtx_pim, vtx_pip;
      TVector3 vtx_beam_pim, vtx_beam_pip, vtx_pim_beam, vtx_pip_beam;
      double dis_beam_pim, dis_beam_pip;
      //      bool pim_vtx_flag = TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPim[0], vtx_beam_pim, vtx_pim_beam, dis_beam_pim);
      //      bool pip_vtx_flag = TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPip[0], vtx_beam_pip, vtx_pip_beam, dis_beam_pip);
      bool pim_vtx_flag = fTrackPim[0]-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam_pim, vtx_pim_beam);
      bool pip_vtx_flag = fTrackPip[0]-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam_pip, vtx_pip_beam);
      dis_beam_pim = (vtx_beam_pim-vtx_pim_beam).Mag();
      dis_beam_pip = (vtx_beam_pip-vtx_pip_beam).Mag();
      bool pipi_vtx_flag = TrackTools::Calc2HelixVertex(fTrackPim[0], fTrackPip[0], vtx_pim, vtx_pip);
      TVector3 vtx_pim_beam_mean = 0.5*(vtx_pim_beam+vtx_beam_pim);
      TVector3 vtx_pip_beam_mean = 0.5*(vtx_pip_beam+vtx_beam_pip);

      TVector3 pim_mom, pip_mom, pim_mom_beam, pip_mom_beam;
      bool pim_mom_flag0 = fTrackPim[0]-> GetMomentum(vtx_pim, pim_mom, true, true);
      bool pip_mom_flag0 = fTrackPip[0]-> GetMomentum(vtx_pip, pip_mom, true, true);
      bool pim_mom_flag1 = fTrackPim[0]-> GetMomentum(vtx_pim_beam, pim_mom_beam, true, true);
      bool pip_mom_flag1 = fTrackPip[0]-> GetMomentum(vtx_pip_beam, pip_mom_beam, true, true);
      TVector3 kn_pim_mom = fBeamLmom.Vect()-fFLmom.Vect()-pim_mom_beam;
      TVector3 kn_pip_mom = fBeamLmom.Vect()-fFLmom.Vect()-pip_mom_beam;
      TLorentzVector pim_lmom, pip_lmom, pim_lmom_beam , pip_lmom_beam;
      pim_lmom.SetVectM(pim_mom, piMass);
      pip_lmom.SetVectM(pip_mom, piMass);
      pim_lmom_beam.SetVectM(pim_mom_beam, piMass);
      pip_lmom_beam.SetVectM(pip_mom_beam, piMass);

      TVector3 vtx_cds, vtx_beam;
      double dltmp=0;
      double dca;
      MathTools::LineToLine(0.5*(vtx_pim+vtx_pip), (pim_mom+pip_mom).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtx_cds, vtx_beam);

      TVector3 vtx_pim_mm, vtx_mm_pim;
      TVector3 vtx_pip_mm, vtx_mm_pip;
      bool pim_vtx_mm_flag = fTrackPim[0]-> GetVertex(vtx_pim_beam_mean,  kn_pim_mom.Unit(), vtx_mm_pim, vtx_pim_mm);
      bool pip_vtx_mm_flag = fTrackPip[0]-> GetVertex(vtx_pip_beam_mean,  kn_pip_mom.Unit(), vtx_mm_pip, vtx_pip_mm);

      TVector3 pim_mom_mm, pip_mom_mm;
      bool pim_mom_flag2 = fTrackPim[0]-> GetMomentum(vtx_pim_mm, pim_mom_mm, true, true);
      bool pip_mom_flag2 = fTrackPip[0]-> GetMomentum(vtx_pip_mm, pip_mom_mm, true, true);
      TLorentzVector pim_lmom_mm, pip_lmom_mm;
      pim_lmom_mm.SetVectM(pim_mom_mm, piMass);
      pip_lmom_mm.SetVectM(pip_mom_mm, piMass);

      bool fiducial_flag = false;
      double beam_out, beam_tof;
      TVector3 vtx_c;
      if( dis_beam_pim<dis_beam_pip ){
	if( GeomTools::GetID(vtx_beam_pim)==CID_Fiducial ) fiducial_flag=true; 
	ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam_pim, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	vtx_c = vtx_pim_beam;
      }
      else{ 
	if( GeomTools::GetID(vtx_beam_pip)==CID_Fiducial ) fiducial_flag=true; 
	ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam_pip, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	vtx_c = vtx_pip_beam;
      }
      TVector3 beam_mom = fTrackBPC->GetMomDir();
      beam_mom.SetMag(beam_out);
      TLorentzVector beam_lmom;
      beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);

      if( fBeamPID==Beam_Kaon ){
	h1 = (TH1F*)rtFile-> Get("KNpipi_status"), h1-> Fill(0);
	if( !pipi_vtx_flag ) h1-> Fill(1);
	if( !pim_vtx_flag ) h1-> Fill(2);
	if( !pip_vtx_flag ) h1-> Fill(3);
	if( !pim_vtx_flag ) h1-> Fill(4);
	if( !pim_mom_flag0 ) h1-> Fill(5);
	if( !pip_mom_flag0 ) h1-> Fill(6);
	if( !pim_mom_flag1 ) h1-> Fill(7);
	if( !pip_mom_flag1 ) h1-> Fill(8);
	if( !pim_vtx_mm_flag ) h1-> Fill(9);
	if( !pip_vtx_mm_flag ) h1-> Fill(10);
      }

      if( pim_mom_flag0 && pip_mom_flag0 && pim_mom_flag1 && pip_mom_flag1 ){
	TVector3 n_mom = (fNCpos-vtx_c);
	double NC_beta_c=(fNCpos-vtx_c).Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
	double nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	n_mom.SetMag(nc_mom);
	TLorentzVector n_lmom;
	n_lmom.SetVectM(n_mom, nMass);
	TLorentzVector knpipi_lmom = beam_lmom+TGT_LMOM-n_lmom-pim_lmom_beam-pip_lmom_beam;
	const double im_npim = (n_lmom+pim_lmom_beam).M();
	const double im_npip = (n_lmom+pip_lmom_beam).M();
	const double mm_knpim = (beam_lmom+TGT_LMOM-n_lmom-pim_lmom_beam).M();
	const double mm_knpip = (beam_lmom+TGT_LMOM-n_lmom-pip_lmom_beam).M();
	const double mm_knpipi = (beam_lmom+TGT_LMOM-n_lmom-pim_lmom_beam-pip_lmom_beam).M();
	const double im_pipi = (pim_lmom+pip_lmom).M();
	const double mm_kn = (fBeamLmom+TGT_LMOM-n_lmom).Mag();
	const double im_npipi = (n_lmom+pim_lmom_beam+pip_lmom_beam).M();
	const double mm_kp_npipi = (beam_lmom+P_LMOM-n_lmom-pim_lmom_beam-pip_lmom_beam).M();

	//*** +-2sigma (default)
	if( Npipi_K0_MIN<im_pipi && im_pipi<Npipi_K0_MAX ) fKNpipi_K0_flag=true;
	if( Npipi_Sm_MIN<im_npim && im_npim<Npipi_Sm_MAX ) fKNpipi_Sm_flag=true;
	if( Npipi_Sp_MIN<im_npip && im_npip<Npipi_Sp_MAX ) fKNpipi_Sp_flag=true;
	//*** +- 2.5sigma
	if( Npipi_K0_MIN25<im_pipi && im_pipi<Npipi_K0_MAX25 ) fKNpipi_K0_flag25=true;
	if( Npipi_Sm_MIN25<im_npim && im_npim<Npipi_Sm_MAX25 ) fKNpipi_Sm_flag25=true;
	if( Npipi_Sp_MIN25<im_npip && im_npip<Npipi_Sp_MAX25 ) fKNpipi_Sp_flag25=true;
	//*** +- 3sigma
	if( Npipi_K0_MIN3<im_pipi && im_pipi<Npipi_K0_MAX3 ) fKNpipi_K0_flag3=true;
	if( Npipi_Sm_MIN3<im_npim && im_npim<Npipi_Sm_MAX3 ) fKNpipi_Sm_flag3=true;
	if( Npipi_Sp_MIN3<im_npip && im_npip<Npipi_Sp_MAX3 ) fKNpipi_Sp_flag3=true;
	//*** +- 4sigma
	if( Npipi_K0_MIN4<im_pipi && im_pipi<Npipi_K0_MAX4 ) fKNpipi_K0_flag4=true;
	if( Npipi_Sm_MIN4<im_npim && im_npim<Npipi_Sm_MAX4 ) fKNpipi_Sm_flag4=true;
	if( Npipi_Sp_MIN4<im_npip && im_npip<Npipi_Sp_MAX4 ) fKNpipi_Sp_flag4=true;

	if( KNpim_MM_Sp_MIN<mm_knpim && mm_knpim<KNpim_MM_Sp_MAX ) fKNpim_Sp_flag=true;
	if( KNpip_MM_Sm_MIN<mm_knpip && mm_knpip<KNpip_MM_Sm_MAX ) fKNpip_Sm_flag=true;
	if( Npipi_Sm_MIN<im_npim && im_npim<Npipi_Sm_MAX ) fKNpipi_Sm_flag=true;
	if( Npipi_Sp_MIN<im_npip && im_npip<Npipi_Sp_MAX ) fKNpipi_Sp_flag=true;
	if( Npipi_N_MIN<mm_knpipi && mm_knpipi<Npipi_N_MAX ) fKNpipi_N_flag=true;
	if( Npipi_N_MIN0<mm_knpipi && mm_knpipi<Npipi_N_MAX0 ) fKNpipi_N0_flag=true;
	if( Npipi_N_MIN1<mm_knpipi && mm_knpipi<Npipi_N_MAX1 ) fKNpipi_N1_flag=true;
	if( Npipi_N_MIN2<mm_knpipi && mm_knpipi<Npipi_N_MAX2 ) fKNpipi_N2_flag=true;

	//*** set data to NpipiData
	fNpipiData-> setBeamLmom(beam_lmom);
	fNpipiData-> setFNLmom(n_lmom);
	fNpipiData-> setPimLmom(pim_lmom);
	fNpipiData-> setPipLmom(pip_lmom);
	fNpipiData-> setVtxCDS(vtx_cds);
	fNpipiData-> setVtxBeam(vtx_beam);
	fNpipiData-> setVtxPim(vtx_pim);
	fNpipiData-> setVtxPim(vtx_pip);
	fNpipiData-> setVtxPimBeam(vtx_pim_beam);
	fNpipiData-> setVtxPipBeam(vtx_pip_beam);
	fNpipiData-> setVtxBeamPim(vtx_beam_pim);
	fNpipiData-> setVtxBeamPip(vtx_beam_pip);
      
	if( fBeamPID==Beam_Kaon ){
	  if( fiducial_flag ){
	    h1= (TH1F*)rtFile-> Get("CDS_GoodTrack_npipi"), h1-> Fill(cdstrackMan->nGoodTrack());
	    h1= (TH1F*)rtFile-> Get("KNpipi_event"), h1-> Fill(0);
	    if( fKNpipi_K0_flag ) h1->Fill(1);
	    if( fKNpipi_Sm_flag ) h1->Fill(2);
	    if( fKNpipi_Sp_flag ) h1->Fill(3);
	    if( fKNpipi_K0_flag && fKNpipi_Sm_flag ) h1-> Fill(4);
	    if( fKNpipi_K0_flag && fKNpipi_Sp_flag ) h1-> Fill(5);
	    if( fKNpipi_Sp_flag && fKNpipi_Sm_flag ) h1-> Fill(6);
	    if( fKNpipi_K0_flag && fKNpipi_Sm_flag && fKNpipi_Sp_flag ) h1-> Fill(7);

	    h1= (TH1F*)rtFile-> Get("KNpipi_event25"), h1-> Fill(0);
	    if( fKNpipi_K0_flag25 ) h1->Fill(1);
	    if( fKNpipi_Sm_flag25 ) h1->Fill(2);
	    if( fKNpipi_Sp_flag25 ) h1->Fill(3);
	    if( fKNpipi_K0_flag25 && fKNpipi_Sm_flag25 ) h1-> Fill(4);
	    if( fKNpipi_K0_flag25 && fKNpipi_Sp_flag25 ) h1-> Fill(5);
	    if( fKNpipi_Sp_flag25 && fKNpipi_Sm_flag25 ) h1-> Fill(6);
	    if( fKNpipi_K0_flag25 && fKNpipi_Sm_flag25 && fKNpipi_Sp_flag25 ) h1-> Fill(7);

	    h1= (TH1F*)rtFile-> Get("KNpipi_event3"), h1-> Fill(0);
	    if( fKNpipi_K0_flag3 ) h1->Fill(1);
	    if( fKNpipi_Sm_flag3 ) h1->Fill(2);
	    if( fKNpipi_Sp_flag3 ) h1->Fill(3);
	    if( fKNpipi_K0_flag3 && fKNpipi_Sm_flag3 ) h1-> Fill(4);
	    if( fKNpipi_K0_flag3 && fKNpipi_Sp_flag3 ) h1-> Fill(5);
	    if( fKNpipi_Sp_flag3 && fKNpipi_Sm_flag3 ) h1-> Fill(6);
	    if( fKNpipi_K0_flag3 && fKNpipi_Sm_flag3 && fKNpipi_Sp_flag3 ) h1-> Fill(7);

	    h1= (TH1F*)rtFile-> Get("KNpipi_event4"), h1-> Fill(0);
	    if( fKNpipi_K0_flag4 ) h1->Fill(1);
	    if( fKNpipi_Sm_flag4 ) h1->Fill(2);
	    if( fKNpipi_Sp_flag4 ) h1->Fill(3);
	    if( fKNpipi_K0_flag4 && fKNpipi_Sm_flag4 ) h1-> Fill(4);
	    if( fKNpipi_K0_flag4 && fKNpipi_Sp_flag4 ) h1-> Fill(5);
	    if( fKNpipi_Sp_flag4 && fKNpipi_Sm_flag4 ) h1-> Fill(6);
	    if( fKNpipi_K0_flag4 && fKNpipi_Sm_flag4 && fKNpipi_Sp_flag4 ) h1-> Fill(7);

	    if( 0.98<mm_knpipi && mm_knpipi<1.06 ){
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pim_test"), h2-> Fill(pim_mass2, pim_mom2);
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pip_test"), h2-> Fill(pip_mass2, pip_mom2);
	      h1 = (TH1F*)rtFile-> Get("pim_dmom_test"), h1-> Fill(pim_lmom_beam.Mag()-pim_mom2);
	      h1 = (TH1F*)rtFile-> Get("pip_dmom_test"), h1-> Fill(pip_lmom_beam.Mag()-pip_mom2);
	    }
	    
	    if( header ){
	      //	    std::cout<<">>> Event Number : "<<header->ev()<<" n pi+ pi- hit fiducial"<<std::endl;
	      tup = (TNtuple*)rtFile->Get("tupNpipi"), tup-> Fill(header->ev());
	    }
	    if( header ){
	      h1_ER-> Fill(13);
	      if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(13);
	    }

	    h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN"), h1-> Fill(im_pipi);
	    if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wSm"), h1-> Fill(im_pipi);
	    if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wSp"), h1-> Fill(im_pipi);
	    if( fKNpipi_Sm_flag || fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wSf"), h1-> Fill(im_pipi);
	    
	    h1 = (TH1F*)rtFile-> Get("Npim_IM_pip"), h1-> Fill(im_npim);
	    if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wK0"), h1-> Fill(im_npim);
	    if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wSp"), h1-> Fill(im_npim);
	    if( fKNpipi_K0_flag || fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wSp_K0"), h1-> Fill(im_npim);
	    
	    h1 = (TH1F*)rtFile-> Get("Npip_IM_pim"), h1-> Fill(im_npip);
	    if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wK0"), h1-> Fill(im_npip);
	    if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wSm"), h1-> Fill(im_npip);
	    if( fKNpipi_K0_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wSm_K0"), h1-> Fill(im_npip);
	    
	    h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM"), h2-> Fill(im_npim, im_npip);
	    
	    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi"), h1->Fill(mm_kn);
	    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_B"), h1->Fill(mm_kn);
	    h1 = (TH1F*)rtFile-> Get("KNpipi_MM"), h1->Fill(mm_knpipi);
	    
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM"), h2-> Fill(mm_knpim, mm_knpip);
	    
	    h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
	    h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip"), h2-> Fill(dis_beam_pim, dis_beam_pip);
	    
	    if( fKNpipi_K0_flag ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_B"), h1->Fill(mm_kn);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_wK0"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wK0"), h1->Fill(mm_knpipi);
	      if( fKNpipi_N_flag ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_B"), h1->Fill(mm_kn);
	      }
	    }
	    if( fKNpipi_K0_flag25 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_25"), h1->Fill(mm_kn);
	    if( fKNpipi_K0_flag3  ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_3"), h1->Fill(mm_kn);
	    if( fKNpipi_K0_flag4  ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_4"), h1->Fill(mm_kn);

	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_B"), h1->Fill(mm_kn);
	      if( fKNpipi_N_flag ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_wN_B"), h1->Fill(mm_kn);
	      }
	    }
	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_25"), h1->Fill(mm_kn);
	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_3"), h1->Fill(mm_kn);
	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_4"), h1->Fill(mm_kn);
	    
	    if( !fKNpipi_K0_flag && !fKNpipi_Sp_flag && !fKNpipi_Sm_flag ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_B"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll"), h1->Fill(mm_knpipi);
	      if( fBeamPion ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_wPi"), h1->Fill(mm_knpipi);
	      if( fBHD_hit.size()==1 ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_BHD1hit"), h1->Fill(mm_knpipi);
	      if( MyTools::IsTarget1(vtx_c, confMan) ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx1"), h1->Fill(mm_knpipi);
	      if( MyTools::IsTarget2(vtx_c, confMan) ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx2"), h1->Fill(mm_knpipi);
	      if( MyTools::IsTarget3(vtx_c, confMan) ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx3"), h1->Fill(mm_knpipi);
	      if( MyTools::IsTarget4(vtx_c, confMan) ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx4"), h1->Fill(mm_knpipi);
	      if( MyTools::IsTarget5(vtx_c, confMan) ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx5"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_pim_mom_woAll"), h2->Fill(mm_knpipi, pim_lmom_beam.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_pip_mom_woAll"), h2->Fill(mm_knpipi, pip_lmom_beam.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_n_mom_woAll"), h2->Fill(mm_knpipi, n_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_NC_seg_woAll"), h2->Fill(mm_knpipi, fNCseg);
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll"), h2-> Fill(mm_knpim, mm_knpip);

	      h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip_woAll"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
	      h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip_woAll"), h2-> Fill(dis_beam_pim, dis_beam_pip);
	      
	      double chi2_pim = fTrackPim[0]-> Chi();
	      double chi2_pip = fTrackPip[0]-> Chi();
	      h1 = (TH1F*)rtFile-> Get("CDC_pim_chi2_Npipi"), h1-> Fill(chi2_pim);
	      h1 = (TH1F*)rtFile-> Get("CDC_pip_chi2_Npipi"), h1-> Fill(chi2_pip);
	      if( chi2_pim<50 && chi2_pip<50 ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_chi2_50"), h1->Fill(mm_knpipi);
	      if( chi2_pim<30 && chi2_pip<30 ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_chi2_30"), h1->Fill(mm_knpipi);
	      
	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_B"), h1->Fill(mm_kn);
		
		h2 = (TH2F*)rtFile-> Get("KNpim_mom"), h2-> Fill(kn_pim_mom.CosTheta(), kn_pim_mom.Mag());
		h2 = (TH2F*)rtFile-> Get("KNpip_mom"), h2-> Fill(kn_pip_mom.CosTheta(), kn_pip_mom.Mag());
		
		h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip_woAll_wN"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
		h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip_woAll_wN"), h2-> Fill(dis_beam_pim, dis_beam_pip);
		
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wS_B"), h1->Fill(mm_kn);

		  h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip_woAll_wN_wS"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
		  h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip_woAll_wN_wS"), h2-> Fill(dis_beam_pim, dis_beam_pip);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wS_B"), h1->Fill(mm_kn);
		}
	      }

	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wS_B"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( !fKNpipi_K0_flag25 && !fKNpipi_Sp_flag25 && !fKNpipi_Sm_flag25 ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_B"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll25"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll25"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25"), h2-> Fill(mm_knpim, mm_knpip);

	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN_B"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN0"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN0_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN0_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN0_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN1"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN1_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN1_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN1_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN2_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN2_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN2_wS_B"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( !fKNpipi_K0_flag3 && !fKNpipi_Sp_flag3 && !fKNpipi_Sm_flag3 ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_B"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll3"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll3"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3"), h2-> Fill(mm_knpim, mm_knpip);

	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN_B"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN0"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN0_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN0_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN0_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN1"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN1_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN1_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN1_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN2_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN2_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN2_wS_B"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( !fKNpipi_K0_flag4 && !fKNpipi_Sp_flag4 && !fKNpipi_Sm_flag4 ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_B"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll4"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll4"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4"), h2-> Fill(mm_knpim, mm_knpip);

	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN_B"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN0"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN0_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN0_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN0_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN1"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN1_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN1_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN1_wS_B"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN2_B"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN2_wS"), h1->Fill(mm_kn);
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN2_wS_B"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_wN"), h2-> Fill(mm_knpim, mm_knpip);
	    }
	    if( !fKNpipi_K0_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woK0"), h1->Fill(mm_knpipi);
	    }

	    //*** n pi- pi+ IM Analysis ***//
	    h1 = (TH1F*)rtFile-> Get("Npipi_IM"), h1-> Fill(im_npipi);
	    h1 = (TH1F*)rtFile-> Get("KP_Npipi_MM"), h1-> Fill(mm_kp_npipi);
	    h2 = (TH2F*)rtFile-> Get("Npipi_IM_KP_Npipi_MM"), h2-> Fill(im_npipi, mm_kp_npipi);
	    h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM"), h2-> Fill(im_npipi, mm_knpipi);
	  }
	
	  if( mcData ){
	    if( reacData->ReactionID()==1098 || reacData->ReactionID()==3000 ){
	      h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit"), h1-> Fill(fYstarLmom.M());
	      if( fDecayMode==-1     )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_Sm"), h1-> Fill(fYstarLmom.M());
	      else if( fDecayMode==1 )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_Sp"), h1-> Fill(fYstarLmom.M());
	      else std::cout<<"  !!! n pi+ pi- event decay mode : "<<fDecayMode<<std::endl;
	      if( fKNpipi_K0_flag ){
		h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_K0"), h1-> Fill(fYstarLmom.M());
		if( fDecayMode==-1     )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_K0_Sm"), h1-> Fill(fYstarLmom.M());
		else if( fDecayMode==1 )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_K0_Sp"), h1-> Fill(fYstarLmom.M());
	      }
	      if( fKNpipi_Sm_flag || fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_Sf"), h1-> Fill(fYstarLmom.M());
		if( fDecayMode==-1     )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_Sf_Sm"), h1-> Fill(fYstarLmom.M());
		else if( fDecayMode==1 )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_Sf_Sp"), h1-> Fill(fYstarLmom.M());
	      }
	      if( !fKNpipi_K0_flag && !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll"), h1-> Fill(fYstarLmom.M());
		if( fDecayMode==-1     )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_Sm"), h1-> Fill(fYstarLmom.M());
		else if( fDecayMode==1 )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_Sp"), h1-> Fill(fYstarLmom.M());
		if( fKNpipi_N_flag ){
		  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_wN"), h1-> Fill(fYstarLmom.M());
		  if( fDecayMode==-1     )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_wN_Sm"), h1-> Fill(fYstarLmom.M());
		  else if( fDecayMode==1 )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_wN_Sp"), h1-> Fill(fYstarLmom.M());
		  if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		    h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_wN_wS"), h1-> Fill(fYstarLmom.M());
		    if( fDecayMode==-1     )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_wN_wS_Sm"), h1-> Fill(fYstarLmom.M());
		    else if( fDecayMode==1 )  h1 = (TH1F*)rtFile-> Get("L1405_pipi_hit_woAll_wN_wS_Sp"), h1-> Fill(fYstarLmom.M());
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    bool IM_Sp_flag2 = false;
    bool IM_Sm_flag2 = false;
    for( int i=0; i<fTrackPim.size(); i++ ){
      TVector3 vtxb, vtxcds;
      double dis;
      //      if( TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPim[i], vtxb, vtxcds, dis) ){
      if( fTrackPim[i]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxb, vtxcds) ){
	TVector3 p;
	if( fTrackPim[0]->GetMomentum(vtxcds, p, true, true) ){
	  TLorentzVector lmom;
	  lmom.SetVectM(p, piMass);
	  double im = (fFLmom+lmom).M();
	  if( GeomTools::GetID(fVtxBeam)==CID_Fiducial && fBeamPID==Beam_Kaon ){
	    h1 = (TH1F*)rtFile-> Get("Npim_IM"), h1-> Fill(im);
	  }
	  if( Sm_MIN<im && im<Sm_MAX ) IM_Sm_flag2=true;
	}
      }
    }

    for( int i=0; i<fTrackPip.size(); i++ ){
      TVector3 vtxb, vtxcds;
      double dis;
      if( TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPip[i], vtxb, vtxcds, dis) ){
      //      if( fTrackPip[i]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxb, vtxcds) ){
	TVector3 p;
	if( fTrackPip[0]->GetMomentum(vtxcds, p, true, true) ){
	  TLorentzVector lmom;
	  lmom.SetVectM(p, piMass);
	  double im = (fFLmom+lmom).M();
	  if( GeomTools::GetID(fVtxBeam)==CID_Fiducial && fBeamPID==Beam_Kaon ){
	    h1 = (TH1F*)rtFile-> Get("Npip_IM"), h1-> Fill(im);
	  }
	  if( Sp_MIN<im && im<Sp_MAX ) IM_Sp_flag2=true;
	}
      }
    }

    for( int i=0; i<fTrackKm.size(); i++ ){
      TVector3 vtxb, vtxcds;
      double dis;
      if( TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackKm[i], vtxb, vtxcds, dis) ){
      //      if( fTrackKm[i]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxb, vtxcds) ){
	TVector3 p;
	if( fTrackKm[0]->GetMomentum(vtxcds, p, true, true) ){
	  TLorentzVector lmom;
	  lmom.SetVectM(p, kpMass);
	  if( GeomTools::GetID(fVtxBeam)==CID_Fiducial && fBeamPID==Beam_Kaon ){
	    h1 = (TH1F*)rtFile-> Get("Nkm_IM"), h1-> Fill((fFLmom+lmom).M());
	  }
	}
      }
    }
    if( GeomTools::GetID(fVtxBeam)==CID_Fiducial && fBeamPID==Beam_Kaon ){
      if( IM_Sp_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wSp"), h1->Fill(kn_lmom.M());
      if( IM_Sm_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wSm"), h1->Fill(kn_lmom.M());
      if( IM_Sp_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wSp_B"), h1->Fill(kn_lmom.M());
      if( IM_Sm_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wSm_B"), h1->Fill(kn_lmom.M());
      if( IM_Sp_flag2 || IM_Sm_flag2 ){
	if( IM_Sm_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wS"), h1->Fill(kn_lmom.M());
	if( IM_Sm_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wS_B"), h1->Fill(kn_lmom.M());
      }
    }
  }

  // for MC
  if( mcData ){
    if( reacData->ReactionID()==1098 || reacData->ReactionID()==3000 ){
      if( fBeamPID==Beam_Kaon ){
	h1 = (TH1F*)rtFile-> Get("L1405_mass_MC"), h1-> Fill(fYstarLmom.M());
	if( fDecayMode==-1 ){
	  h1 = (TH1F*)rtFile-> Get("FL_Sm_MC"), h1-> Fill((fVertex-fDecayPos).Mag());
	}
	if( fDecayMode==0  ){
	  h1 = (TH1F*)rtFile-> Get("FL_S0_MC"), h1-> Fill((fVertex-fDecayPos).Mag());
	}
	if( fDecayMode==1  ){
	  h1 = (TH1F*)rtFile-> Get("FL_Sp_MC"), h1-> Fill((fVertex-fDecayPos).Mag());
	}
      }
    }
  }
  //  std::cout<<"===== HistManwMC::fill FINISH ====="<<std::endl;
}

void HistManwMC::anaMC()
{
  TH1F *h1;
  TH2F *h2;
  fFNFlag = false;

  int nT0_MC=0;
  DetectorHit *T0hit = 0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==CID_T0 && hit->parentID()==0 && hit->pdg()==321){
      nT0_MC++;
      T0hit = hit;
    }
  }

  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==CID_NC ){
      if( hit->pdg()==2112 && hit->parentID()==0 ){
	h1 = (TH1F*)rtFile-> Get("NC_dE_MC"), h1-> Fill(hit->adc());

	bool hit_flag = false;
	for( int j=0; j<blMan->nNC(); j++ ){
	  if( blMan->NC(j)->seg()==hit->channelID()+1 ){
	    HodoscopeLikeHit *hit2 = blMan->NC(j);
	    hit_flag = true;
	    h1 = (TH1F*)rtFile-> Get("NC_dE_diff_MC"), h1-> Fill(hit->adc()-hit2->emean());
	    if( fT0time!=DBL_MIN && T0hit ){
	      double tof_data = hit2->ctmean()-fT0time;
	      double tof_MC = hit-> tdc()+T0hit->tdc();
	      h1 = (TH1F*)rtFile-> Get("T0NC_diff_MC"), h1-> Fill(tof_MC-tof_data);
	    }
	  }
	
	  if( blMan->NC(j)->emean()>NC_THRE ){
	    fFNFlag = true;
	  }
	}
	if( !hit_flag ){
	  h1 = (TH1F*)rtFile-> Get("NC_wo_hit_dE_MC"), h1-> Fill(hit->adc());
	}
      }
    }
  }

  // 1098 K- d -> L(1405) n
  for( int i=0; i<reacData->InitParticleSize(); i++ ){
    if( reacData->InitPDG(i)==-321 ) fBeamLmomMC=reacData->GetInitParticle(i)*MeV;
    else fTargetLmom=reacData->GetInitParticle(i)*MeV;
  }

  //  std::cout<<"  ReactionID : "<<reacData->ReactionID()<<std::endl;
  if( reacData->ReactionID()==1098 || reacData->ReactionID()==3000 ){
    for( int i=0; i<reacData->ParticleSize(); i++ ){
      if( reacData-> PDG(i)==2112  ) fNLmom1 = reacData->GetParticle(i)*MeV;
      if( reacData-> PDG(i)==13122 ) fYstarLmom = reacData->GetParticle(i)*MeV;
    }

    int YstarID;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->pdgID()==13122 ){
	YstarID=track->trackID();
	fVertex = track->vertex()*mm;
      }
    }

    int SigmaID;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track=mcData->track(i);
      if( track->parentTrackID()==YstarID ){
	if( track->pdgID()==3222 ){
	  fDecayMode=1;
	  fSLmom.SetVectM(track->momentum()*MeV, spMass);
	  SigmaID=track->trackID();
	}
	if( track->pdgID()==3112 ){
	  fDecayMode=-1;
	  fSLmom.SetVectM(track->momentum()*MeV, smMass);
	  SigmaID=track->trackID();
	}
	if( track->pdgID()==3212 ){
	  fDecayMode=0;
	  fSLmom.SetVectM(track->momentum()*MeV, s0Mass);
	  SigmaID=track->trackID();
	}
	if( track->pdgID()==211  ) fPiLmom1.SetVectM(track->momentum(), piMass);
	if( track->pdgID()==-211 ) fPiLmom1.SetVectM(track->momentum(), piMass);
	if( track->pdgID()==111  ) fPiLmom1.SetVectM(track->momentum(), 0.1349766);
      }      
    }

    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track=mcData->track(i);
      if( track->parentTrackID()==SigmaID ){
	fDecayPos = track->vertex()*mm;
	if( track->pdgID()==211  ) fPiLmom2.SetVectM(track->momentum()*MeV, piMass);
	if( track->pdgID()==-211 ) fPiLmom2.SetVectM(track->momentum()*MeV, piMass);
	if( track->pdgID()==111  ) fPiLmom2.SetVectM(track->momentum()*MeV, 0.1349766);

	if( track->pdgID()==2112 ) fNLmom2.SetVectM(track->momentum()*MeV, nMass);
      }
    }
  }
}

//*** after this only test ***//
void HistManwMC::dump()
{
  std::cout<<" detData="<<detData<<std::endl;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==CID_T0 ){
      std::cout<<">>> Detector Hit  T0 seg"<<hit->channelID()+1<<std::endl;
      std::cout<<">>>    time : "<<hit->time()<<std::endl; 
      std::cout<<">>>    dE   : "<<hit->de()<<"[MeV/c]"<<std::endl; 
      std::cout<<">>>    hitZ : "<<hit->pos().Z()<<std::endl;
    }
  }

  for( int i=0; i<reacData->InitParticleSize(); i++ ){
    int pdgId = reacData-> InitPDG(i);
    TLorentzVector lmom = reacData-> GetInitParticle(i);
    if( pdgId==-321 ){
      std::cout<<">>> Init K- ====="<<std::endl;
      std::cout<<">>>   lmom ("<<lmom.X()<<", "<<lmom.Y()<<", "<<lmom.Z()<<", "<<lmom.T()<<")"<<std::endl;
    }
  }
}

void HistManwMC::finit()
{
  rtFile-> Write();
  // This method occure abort if param file don't exist
  //  confMan-> SaveParams();
  //  confMan-> SaveCode();
  rtFile-> Close();

}

void HistManwMC::print(const int &evnum)
{
  if( evnum<0 ) std::cout<<"===== HistManwMC::print unknown evnum ====="<<std::endl;
  else  std::cout<<"===== HistManwMC::print Event Number : "<<evnum<<" ====="<<std::endl;
  if( fD5mom>0 ) std::cout<<"> Beam Mom : "<<fD5mom<<"[GeV/c]"<<std::endl;
  if( fFPID==F_Gamma ) std::cout<<"> Forward Gamma"<<std::endl;
  else if( fFPID==F_Neutron ) std::cout<<"> Forward Neutron"<<std::endl;

  std::cout<<"> CDS GoodTrack : "<<cdstrackMan->nGoodTrack()<<std::endl;
  for( int i=0; i<fCDSPID.size(); i++ ){
    std::string str;
    char c_str[512];
    if( fCDSPID[i]==CDS_PiMinus  ) str += ">> pi- : ";
    if( fCDSPID[i]==CDS_PiPlus   ) str += ">> pi+ : ";
    if( fCDSPID[i]==CDS_Kaon     ) str += ">> K-  : ";
    if( fCDSPID[i]==CDS_Proton   ) str += ">> p   : ";
    if( fCDSPID[i]==CDS_Deuteron ) str += ">> d   : ";
    if( fCDSPID[i]==98           ) str += ">> e-  : ";
    if( fCDSPID[i]==99           ) str += ">> e+  : ";
    sprintf(c_str, "mass2=%5.4lf :", fCDSmass2[i]);
    str += c_str;
    sprintf(c_str, "mom=%5.4lf :", fCDSmom[i]);
    str += c_str;

    std::cout<<str<<std::endl;
  }



}

bool HistManwMC::Clear(const bool flag)
{
  if( blMan ) blMan->Clear();
  if( cdsMan ) cdsMan->Clear();
  if( bltrackMan ) bltrackMan->Clear();
  if( cdstrackMan ) cdstrackMan->Clear();
  fBeamPID = Beam_Other;
  fBeamKaon=false;
  fBeamPion=false;
  fT0time = DBL_MIN;
  fD5mom = DBL_MIN;
  fTrackBPC = 0;

  fNCtime = DBL_MAX;
  fNCseg = -1;
  fNCdE = 0.0;
  fBeamLmom.SetXYZT(DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX);

  fT0pos.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
  fNCpos.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
  fVtxCDS.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
  fVtxBeam.SetXYZ(DBL_MAX, DBL_MAX, DBL_MAX);

  fBHD_hit.clear();
  fT0_hit.clear();
  fBPD_hit.clear();
  fCDH_hit.clear();
  fBVC_hit.clear();
  fCVC_hit.clear();
  fPC_hit.clear();
  for( int i=0; i<8; i++ ) fNC_hit[i].clear();
  fBD_hit.clear();

  fCDSPID.clear();
  fCDSbeta.clear();
  fCDSmass2.clear();
  fCDSmom.clear();

  fTrackPim.clear();
  fTrackKm.clear();
  fTrackPip.clear();
  fTrackP.clear();
  fTrackD.clear();

  fFPID = F_Other;
  fNCbeta = DBL_MIN;
  fFLmom.SetXYZT(DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX);
  fKNpipi_K0_flag = false;
  fKNpipi_Sp_flag = false;
  fKNpipi_Sm_flag = false;
  fKNpipi_K0_flag25 = false;
  fKNpipi_Sp_flag25 = false;
  fKNpipi_Sm_flag25 = false;
  fKNpipi_K0_flag3 = false;
  fKNpipi_Sp_flag3 = false;
  fKNpipi_Sm_flag3 = false;
  fKNpipi_K0_flag4 = false;
  fKNpipi_Sp_flag4 = false;
  fKNpipi_Sm_flag4 = false;
  fKNpipi_N_flag = false;
  fKNpipi_N0_flag = false;
  fKNpipi_N1_flag = false;
  fKNpipi_N2_flag = false;
  fKNpim_Sp_flag = false;
  fKNpip_Sm_flag = false;

  return flag;
}

