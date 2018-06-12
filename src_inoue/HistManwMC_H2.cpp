#include "HistManwMC.h"
#include "MyParam.h"

static const int TGT_A=1;
static const TLorentzVector TGT_LMOM=P_LMOM;
//static const TLorentzVector TGT_LMOM=TLorentzVector(0., 0., 0., ThreeHeMass);

HistManwMC::HistManwMC(TFile *f, ConfMan *conf)
  : rtFile(f), confMan(conf), blMan(0), cdsMan(0), bltrackMan(0), cdstrackMan(0)
{
  Clear();
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
      h1 = (TH1F*)rtFile-> Get(Form("trig_%d", i)), h1-> Fill(header->pattern(i));
    }
  }

  TH1F *h1_N  = (TH1F*)rtFile-> Get("N_Reduction");
  TH1F *h1_ER  = (TH1F*)rtFile->Get("EventReduction");
  bool trigNC_flag = false;

  if( header ){
    if( header->IsTrig(Trig_Neutral) ) trigNC_flag=true;
  }
  else{
    trigNC_flag=simReader->isNC();
  }
  bool trigCDH2 = false;
  if( header ){
    if( header-> trigmode2(Mode_KCDH2, confMan) ) trigCDH2=true;
  }
  else{
    trigCDH2=simReader->isCDH2();
  }

  //******************//
  //*** for CDS IM ***//
  //******************//
  //*** Common Part *********************************************************************************************************************************//
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
	  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtxcds\
                                , vtxb);

          double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(fT0pos, vtxb, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
          TVector3 beam_mom = fTrackBPC->GetMomDir();
          beam_mom.SetMag(beam_out);
          TLorentzVector beam_lmom;
          beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);

          if( GeomTools::GetID(vtxb)==CID_Fiducial ){
            if( fBeamPID==Beam_Kaon ){
              TLorentzVector pipi_missing_lmom = beam_lmom+TGT_LMOM-pim_lmom-pip_lmom;
              TVector3 pipi_missing_mom = pipi_missing_lmom.Vect();
              TVector3 mm_dir = pipi_missing_mom.Unit();
              mm_dir.SetMag(1/mm_dir.Z());
              TVector3 nc_calc_pos = fVtxBeam+mm_dir*(1490.-fVtxBeam.Z());

              double im_pipi = (pim_lmom+pip_lmom).M();
              double mm_kn_pipi = (beam_lmom+TGT_LMOM-pim_lmom-pip_lmom).M();
              double mm_kn_pim = (beam_lmom+TGT_LMOM-pim_lmom).M();
              double mm_kn_pip = (beam_lmom+TGT_LMOM-pip_lmom).M();
              h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi"), h1-> Fill(im_pipi);
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

	  bool L_flag=false;
	  if( L_MIN<(pim_lmom+p_lmom).M() && (pim_lmom+p_lmom).M()<L_MAX ) L_flag=true;

	  if( GeomTools::GetID(vtxb)==CID_Fiducial ){
	    if( fBeamPID==Beam_Kaon ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_ppim"), h1-> Fill((pim_lmom+p_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KCDSppim_MM"), h1-> Fill((beam_lmom+TGT_LMOM-pim_lmom-p_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KCDSppim_MM2"), h1-> Fill((beam_lmom+TGT_LMOM-pim_lmom-p_lmom).M2());

	      if( L_flag && cdstrackMan->nGoodTrack()==2 ){
		TVector3 ppim_mom=pim_mom+p_mom;
		TLorentzVector l_lmom; l_lmom.SetVectM(ppim_mom, lMass);
		TLorentzVector kppim_lmom = beam_lmom+TGT_LMOM-pim_lmom-p_lmom;
		TLorentzVector kl_lmom = beam_lmom+TGT_LMOM-l_lmom;
		h1 = (TH1F*)rtFile-> Get("KCDSppim_MM_wL"), h1-> Fill(kppim_lmom.M());
		h1 = (TH1F*)rtFile-> Get("KCDSppim_MM2_wL"), h1-> Fill(kppim_lmom.M2());
		h1 = (TH1F*)rtFile-> Get("KL_MM"), h1-> Fill(kl_lmom.M());
		h1 = (TH1F*)rtFile-> Get("KL_MM2"), h1-> Fill(kl_lmom.M2());
		if( KP_L_MM2_MIN<kl_lmom.M2() && kl_lmom.M2()<KP_L_MM2_MAX ){
		  TLorentzVector l_lmom_CM = l_lmom;
		  TVector3 LabToCM = -(beam_lmom+TGT_LMOM).BoostVector();
		  l_lmom_CM.Boost(LabToCM);
		  TVector3 pi_mom = -l_lmom_CM.Vect();
		  h2 = (TH2F*)rtFile-> Get("KL_pi0_kin"), h2-> Fill(pi_mom.CosTheta(), pi_mom.Mag());
		  h2 = (TH2F*)rtFile-> Get("Kpi0_L_kin"), h2-> Fill(l_lmom.Vect().CosTheta(), l_lmom.Vect().Mag());
		}
	      }	      
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
              TVector3 mean_CDC_beam = 0.5*(vtx_mean+BPC_calc_pos);
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_dx_dy"), h2-> Fill(diff.X(), diff.Y());
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_z_dx"),  h2-> Fill(vtx_mean.Z(), diff.X());
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_z_dy"),  h2-> Fill(vtx_mean.Z(), diff.Y());
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_x_dx"),  h2-> Fill(mean_CDC_beam.X(), diff.X());
              h2 = (TH2F*)rtFile-> Get("CDS2_BPC_y_dy"),  h2-> Fill(mean_CDC_beam.Y(), diff.Y());
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

  //************************//
  //*** Forward Neutron  ***//
  //************************//
  if( fFPID==F_Neutron && trigNC_flag ){
    TLorentzVector kn_lmom = fBeamLmom+TGT_LMOM-fFLmom;

    if( fBeamPID==Beam_Kaon ){
      if( simReader ){
	Track *n_track = simReader->trace(CID_NC, fNCseg);
	ReactionData *reacData = simReader->getReactionData();
	TLorentzVector beam_lmom_MC = 0.001*reacData->GetInitParticle(0);
	TLorentzVector n_lmom_MC;
	n_lmom_MC.SetVectM(0.001*n_track->momentum(), nMass);
	double mm_kn_MC = (beam_lmom_MC+TGT_LMOM-n_lmom_MC).Mag();
	h1 = (TH1F*)rtFile-> Get("KN_MM_diff_MC"), h1-> Fill(mm_kn_MC-kn_lmom.M());
	h1 = (TH1F*)rtFile-> Get("KN_MM_MC"), h1-> Fill(mm_kn_MC);
      }
      h1 = (TH1F*)rtFile-> Get("KN_MM_wo"), h1->Fill(kn_lmom.M());
      if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim_wo"), h1->Fill(kn_lmom.M());
      if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip_wo"), h1->Fill(kn_lmom.M());
      if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km_wo"), h1->Fill(kn_lmom.M());
      if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p_wo"), h1->Fill(kn_lmom.M());

      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	if( header ){
	  if( header->IsTrig(Trig_Neutral) ){
	    h1_N-> Fill(11);
	  }
	  //	  std::cout<<"  on Target"<<std::endl;
	  h1_ER->Fill(11);
	}
	else {
	  h1_ER->Fill(11);
	}
	
	h1 = (TH1F*)rtFile-> Get("KN_MM"), h1->Fill(kn_lmom.M());
	if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim"), h1->Fill(kn_lmom.M());
	if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip"), h1->Fill(kn_lmom.M());
	if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km"), h1->Fill(kn_lmom.M());
	if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p"), h1->Fill(kn_lmom.M());

	int NClay = (fNCseg-1)/16+1;
	int NCseg2 = fNCseg%16;
	if( NCseg2==0 ) NCseg2=16;
	h2 = (TH2F*)rtFile-> Get("KN_MM_NCseg"), h2-> Fill(kn_lmom.M(), NCseg2);
	h2 = (TH2F*)rtFile-> Get("KN_MM_NClay"), h2-> Fill(kn_lmom.M(), NClay);

	h2 = (TH2F*)rtFile-> Get("overbeta_NCseg"), h2-> Fill(1./fNCbeta, NCseg2);
	h2 = (TH2F*)rtFile-> Get("overbeta_NClay"), h2-> Fill(1./fNCbeta, NClay);
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
	  double pim_mm = (fBeamLmom+P_LMOM-lmom).M();
	  double kn_mm = (fBeamLmom+P_LMOM-fFLmom).M();
	  double im = (fFLmom+lmom).M();
	  double mm = (fBeamLmom+P_LMOM-fFLmom-lmom).M();
	  double mm2 = (fBeamLmom+P_LMOM-fFLmom-lmom).M2();
	  if( GeomTools::GetID(fVtxBeam)==CID_Fiducial && fBeamPID==Beam_Kaon ){
	    h1 = (TH1F*)rtFile-> Get("Npim_IM"), h1-> Fill(im);
	    h1 = (TH1F*)rtFile-> Get("KNpim_MM"), h1-> Fill(mm);
	    h1 = (TH1F*)rtFile-> Get("KNpim_MM2"), h1-> Fill(mm2);
	    h2 = (TH2F*)rtFile-> Get("KNpim_MM_Npim_IM"), h2-> Fill(mm, im);
	    h2 = (TH2F*)rtFile-> Get("Kpim_MM_Npim_IM"), h2-> Fill(pim_mm, im);
	    if( Sm_MIN<im && im<Sm_MAX ){
	      h1 = (TH1F*)rtFile-> Get("KNpim_MM_wSm"), h1-> Fill(mm);
	    }
	    if( KP_Npim_MM2_MIN<mm2 && mm2<KP_Npim_MM2_MAX ){
	      h2 = (TH2F*)rtFile-> Get("Kpim_MM_Npim_IM_mm_pip"), h2-> Fill(pim_mm, im);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_Npim_IM_mm_pip"), h2-> Fill(kn_mm, im);
	      if( kn_mm<KP_K0_MIN && KP_K0_MAX<kn_mm && Sm_MIN<im && im<Sm_MAX ){
		TLorentzVector Sm_lmom_CM=fFLmom+lmom;
		TVector3 LabToCM = -(fBeamLmom+P_LMOM).BoostVector();
		Sm_lmom_CM.Boost(LabToCM);
		h2 = (TH2F*)rtFile-> Get("KN_MM_Npim_IM_mm_pip_wSm_woK0"), h2-> Fill(kn_mm, im);
		h1 = (TH1F*)rtFile-> Get("Sm_ang_Npim"), h1-> Fill(cos(fBeamLmom.Angle(Sm_lmom_CM.Vect())));
	      }
	    }
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
	  double kn_mm = (fBeamLmom+P_LMOM-fFLmom).M();
	  double pip_mm = (fBeamLmom+P_LMOM-lmom).M();
	  double im = (fFLmom+lmom).M();
	  double mm = (fBeamLmom+P_LMOM-fFLmom-lmom).M();
	  double mm2 = (fBeamLmom+P_LMOM-fFLmom-lmom).M2();
	  if( GeomTools::GetID(fVtxBeam)==CID_Fiducial && fBeamPID==Beam_Kaon ){
	    h1 = (TH1F*)rtFile-> Get("Npip_IM"), h1-> Fill(im);
	    h1 = (TH1F*)rtFile-> Get("KNpip_MM"), h1-> Fill(mm);
	    h1 = (TH1F*)rtFile-> Get("KNpip_MM2"), h1-> Fill(mm2);
	    h2 = (TH2F*)rtFile-> Get("KNpip_MM_Npip_IM"), h2-> Fill(mm, im);
	    h2 = (TH2F*)rtFile-> Get("Kpip_MM_Npip_IM"), h2-> Fill(pip_mm, im);
	    if( Sm_MIN<im && im<Sm_MAX ){
	      h1 = (TH1F*)rtFile-> Get("KNpip_MM_wSp"), h1-> Fill(mm);
	    }
	    if( KP_Npip_MM2_MIN<mm2 && mm2<KP_Npip_MM2_MAX ){
	      h2 = (TH2F*)rtFile-> Get("Kpip_MM_Npip_IM_mm_pim"), h2-> Fill(pip_mm, im);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_Npip_IM_mm_pim"), h2-> Fill(kn_mm, im);
	      if( kn_mm<KP_K0_MIN && KP_K0_MAX<kn_mm && Sm_MIN<im && im<Sm_MAX ){
		TLorentzVector Sp_lmom_CM=fFLmom+lmom;
		TVector3 LabToCM = -(fBeamLmom+P_LMOM).BoostVector();
		Sp_lmom_CM.Boost(LabToCM);
		h2 = (TH2F*)rtFile-> Get("KN_MM_Npip_IM_mm_pim_wSp_woK0"), h2-> Fill(kn_mm, im);
		h1 = (TH1F*)rtFile-> Get("Sp_ang_Npip"), h1-> Fill(cos(fBeamLmom.Angle(Sp_lmom_CM.Vect())));
	      }
	    }
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
      if( IM_Sp_flag2 || IM_Sm_flag2 ){
	if( IM_Sm_flag2 ) h1 = (TH1F*)rtFile->Get("KN_MM_wS"), h1->Fill(kn_lmom.M());
      }
    }
  }

  if( fTrackKm.size()==1 ){
    TVector3 vtx_beam, vtx_cds;
    TVector3 km_mom;
    bool vtx_flag = fTrackKm[0]-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_cds, vtx_beam);
    bool mom_flag = fTrackKm[0]-> GetMomentum(vtx_cds, km_mom, true, true);
    TLorentzVector km_lmom;
    km_lmom.SetVectM(km_mom, kpMass);

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
    TVector3 beam_mom = fTrackBPC-> GetMomDir();
    beam_mom.SetMag(beam_out);
    TLorentzVector beam_lmom;
    beam_lmom.SetVectM(beam_mom, kpMass);

    if( fFPID==F_Neutron && trigNC_flag ){
      if( vtx_flag && mom_flag ){
	TVector3 n_mom = (fNCpos-vtx_beam);
	double NC_beta_c;
	if( simReader ) NC_beta_c=((fNCpos-vtx_beam).Mag()-2.5)/((fNCtime-fT0time-beam_tof)*100.*Const);
	else NC_beta_c=(fNCpos-vtx_beam).Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
	double nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	n_mom.SetMag(nc_mom);
	TLorentzVector n_lmom;
	n_lmom.SetVectM(n_mom, nMass);

	double kn_kn_mm = (beam_lmom+TGT_LMOM-n_lmom-km_lmom).M();
	h1 = (TH1F*)rtFile-> Get("KNkm_MM"),h1-> Fill(kn_kn_mm);
      }
    }
  }
  //*************************************************************************************************************************************************//
  
  if( fTrackPim.size()==1 && fTrackPip.size()==1 ){  
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
    //    bool pim_vtx_flag = TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPim[0], vtx_beam_pim, vtx_pim_beam, dis_beam_pim);
    //    bool pip_vtx_flag = TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPip[0], vtx_beam_pip, vtx_pip_beam, dis_beam_pip);
    bool pim_vtx_flag = fTrackPim[0]-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam_pim, vtx_pim_beam);
    bool pip_vtx_flag = fTrackPip[0]-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam_pip, vtx_pip_beam);
    dis_beam_pim = (vtx_beam_pim-vtx_pim_beam).Mag();
    dis_beam_pip = (vtx_beam_pip-vtx_pip_beam).Mag();
    bool pipi_vtx_flag = TrackTools::Calc2HelixVertex(fTrackPim[0], fTrackPip[0], vtx_pim, vtx_pip);
    TVector3 vtx_pim_beam_mean = 0.5*(vtx_pim_beam+vtx_beam_pim);
    TVector3 vtx_pip_beam_mean = 0.5*(vtx_pip_beam+vtx_beam_pip);
    TVector3 vtx_pipi = 0.5*(vtx_pim+vtx_pip);

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
    double dis;
    MathTools::LineToLine(0.5*(vtx_pim+vtx_pip), (pim_mom+pip_mom).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtx_cds, vtx_beam);
    const double dca_pipi_beam = dca;
    
    TVector3 vtx_pim_mm, vtx_mm_pim;
    TVector3 vtx_pip_mm, vtx_mm_pip;
    //    bool pim_vtx_mm_flag = TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPim[0], vtx_mm_pim, vtx_pim_mm, dis);
    //    bool pip_vtx_mm_flag = TrackTools::CalcLineHelixVertex(fTrackBPC, fTrackPip[0], vtx_mm_pip, vtx_pip_mm, dis);
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
    double min_dca;
    if( dis_beam_pim<dis_beam_pip ){
      if( GeomTools::GetID(vtx_beam_pim)==CID_Fiducial ) fiducial_flag=true; 
      ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam_pim, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      vtx_c = vtx_pim_beam;
      min_dca=dis_beam_pim;
    }
    else{ 
      if( GeomTools::GetID(vtx_beam_pip)==CID_Fiducial ) fiducial_flag=true; 
      ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam_pip, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      vtx_c = vtx_pip_beam;
      min_dca=dis_beam_pip;
    }
    TVector3 beam_mom = fTrackBPC->GetMomDir();
    beam_mom.SetMag(beam_out);
    TLorentzVector beam_lmom;
    beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);

    if( pim_mom_flag0 && pip_mom_flag0 && pim_mom_flag1 && pip_mom_flag1 ){
      if( -0.1<vtx_c.Z() && vtx_c.Z()<0.1 ) h2 = (TH2F*)rtFile-> Get("Vtx_XY_min_z0"), h2-> Fill(vtx_c.X(), vtx_c.Y());
      if( vtx_c.Perp()<2.0 ) h1 = (TH1F*)rtFile-> Get("Vtx_Z_min"), h1-> Fill(vtx_c.Z());
      h1 = (TH1F*)rtFile-> Get("DCA_pipi_beam_woF"), h1-> Fill(dca_pipi_beam);
      if( fiducial_flag ) h1 = (TH1F*)rtFile-> Get("DCA_pipi_beam"), h1-> Fill(dca_pipi_beam);

      h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi"), h2-> Fill(vtx_cds.X(), vtx_cds.Y());
      h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi"), h1-> Fill(vtx_cds.Z());

      //****************************//
      //*** p(K-, pi+ pi-) study ***//
      //****************************//
      TLorentzVector pipi_missing_lmom = beam_lmom+TGT_LMOM-pim_lmom-pip_lmom;
      TVector3 LabToCM = -(beam_lmom+TGT_LMOM).BoostVector();
      TVector3 pipi_missing_mom = pipi_missing_lmom.Vect();
      TVector3 mm_dir = pipi_missing_mom.Unit();
      mm_dir.SetMag(1/mm_dir.Z());
      TVector3 nc_calc_pos = fVtxBeam+mm_dir*(1490.-fVtxBeam.Z());
      
      double im_pipi = (pim_lmom+pip_lmom).M();
      double mm_kp_pipi = pipi_missing_lmom.M();
      double mm_kp_pim = (beam_lmom+TGT_LMOM-pim_lmom).M();
      double mm_kp_pip = (beam_lmom+TGT_LMOM-pip_lmom).M();

      TVector3 n_mom = (fNCpos-vtx_c);
      double NC_beta_c;
      if( simReader ) NC_beta_c=((fNCpos-vtx_c).Mag()-2.5)/((fNCtime-fT0time-beam_tof)*100.*Const);
      else NC_beta_c=(fNCpos-vtx_c).Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
      double nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
      n_mom.SetMag(nc_mom);
      TLorentzVector n_lmom;
      n_lmom.SetVectM(n_mom, nMass);

      int NClay = (fNCseg-1)/16+1;
      int NCseg2 = fNCseg%16;
      if( NCseg2==0 ) NCseg2=16;

      bool missing_n_flag=false;
      bool missing_L_flag=false;
      bool K0_flag2=false;
      bool K0_flag3=false;
      bool K0_flag4=false;
      bool K0_flag5=false;
      bool mm_Sm_flag=false;
      bool mm_Sp_flag=false;

      if( Kpipi_N_MIN<mm_kp_pipi && mm_kp_pipi<Kpipi_N_MAX ) missing_n_flag=true;
      if( Kpipi_L_MIN<mm_kp_pipi && mm_kp_pipi<Kpipi_L_MAX ) missing_L_flag=true;
      if( K0_MIN <im_pipi && im_pipi<K0_MAX  ) K0_flag2=true;
      if( K0_MIN3<im_pipi && im_pipi<K0_MAX3 ) K0_flag3=true;
      if( K0_MIN4<im_pipi && im_pipi<K0_MAX4 ) K0_flag4=true;
      if( K0_MIN5<im_pipi && im_pipi<K0_MAX5 ) K0_flag5=true;
 
      if( fBeamPID==Beam_Kaon && trigCDH2 && fiducial_flag){
	h1 = (TH1F*)rtFile-> Get("KPpipi_MM_noCut"), h1-> Fill(mm_kp_pipi);
	if( fTrackP.size()>0 ) h1 = (TH1F*)rtFile-> Get("KPpipi_MM_noCut_wp"), h1-> Fill(mm_kp_pipi);

	if( fTrackP.size()==1 ){
	  CDSTrack *p_track = fTrackP[0];
	  TVector3 vtx_pipi_p, vtx_p_pipi;
	  if( p_track-> GetVertex(0.5*(vtx_pim+vtx_pip), pipi_missing_mom.Unit(), vtx_pipi_p, vtx_p_pipi) ){
	    TVector3 p_mom;
	    if( p_track-> GetMomentum(vtx_p_pipi, p_mom, true, true) ){
	      TLorentzVector p_lmom; p_lmom.SetVectM(p_mom, pMass);
	      TVector3 vtx_LP = 0.5*(vtx_pipi_p+vtx_p_pipi);
	      double fl_L = (vtx_pipi-vtx_LP).Mag();
	      if( vtx_pipi.Z()>vtx_LP.Z() ) fl_L *= -1.0;
	      double dca_pipi_p = (vtx_pipi_p-vtx_p_pipi).Mag();
	      double kpipip_mm = (pipi_missing_lmom-p_lmom).M();
	      double kpipip_mm2 = (pipi_missing_lmom-p_lmom).M2();
	      h1 = (TH1F*)rtFile-> Get("KPpipip_MM2"), h1-> Fill(kpipip_mm2);
	      if( missing_L_flag ){
		TLorentzVector l_lmom; l_lmom.SetVectM(pipi_missing_mom, lMass);
		double Lp_mm = (l_lmom+TGT_LMOM-p_lmom).M();
		h1 = (TH1F*)rtFile-> Get("DCA_LP"), h1-> Fill(dca_pipi_p);
		h1 = (TH1F*)rtFile-> Get("FL_L"), h1-> Fill(fl_L);
		h2 = (TH2F*)rtFile-> Get("KPpipip_MM2_KPpipi_PL_MM_mmL"), h2-> Fill(kpipip_mm2, Lp_mm);

		h2 = (TH2F*)rtFile-> Get("KPpipi_mom_wp_mmL"), h2-> Fill(pipi_missing_lmom.CosTheta(), pipi_missing_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("LP_vtx_XY"), h2-> Fill(vtx_LP.X(), vtx_LP.Y());
		h2 = (TH2F*)rtFile-> Get("LP_vtx_ZX"), h2-> Fill(vtx_LP.Z(), vtx_LP.X());
		h2 = (TH2F*)rtFile-> Get("LP_vtx_ZY"), h2-> Fill(vtx_LP.Z(), vtx_LP.Y());

		if( GeomTools::GetID(vtx_LP)==CID_Fiducial ){
		  h1 = (TH1F*)rtFile-> Get("KPpipi_PL_MM_vtx"), h1-> Fill(Lp_mm);
		}
	      }
	    }
	  }
	}
	if( !K0_flag3 ){
	  h1 = (TH1F*)rtFile-> Get("KPpipi_MM_noCut_woK0"), h1-> Fill(mm_kp_pipi);
	}
	if( missing_L_flag ){
	  h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wL"), h1-> Fill(im_pipi);
	  h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM_mmL"), h2-> Fill(mm_kp_pip, mm_kp_pim);
	  h2 = (TH2F*)rtFile-> Get("KPpipi_mom_mmL"), h2-> Fill(pipi_missing_lmom.CosTheta(), pipi_missing_lmom.Vect().Mag());
	  if( !K0_flag3 ){
	    h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM_mmL_woK0"), h2-> Fill(mm_kp_pip, mm_kp_pim);
	  }
	}
      }
 
      if( fBeamPID==Beam_Kaon && trigCDH2 && cdstrackMan->nGoodTrack()==2 ){
	//	if( fiducial_flag && GeomTools::GetID(vtx_pipi)==CID_Fiducial ){
	if( fiducial_flag ){
	  h1 = (TH1F*)rtFile-> Get("KPpipi_MM"), h1-> Fill(mm_kp_pipi);
	  if( K0_flag2 ){
	    h1 = (TH1F*)rtFile-> Get("KPpipi_MM_wK0"), h1-> Fill(mm_kp_pipi);
	  }
	  h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM"), h2-> Fill(mm_kp_pip, mm_kp_pim);

	  h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_NCseg"), h2-> Fill((beam_lmom+TGT_LMOM-n_lmom).M(), NCseg2);
	  h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_NClay"), h2-> Fill((beam_lmom+TGT_LMOM-n_lmom).M(), NClay);

	  if( missing_n_flag ){
	    h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_NCseg"), h2-> Fill((beam_lmom+TGT_LMOM-n_lmom).M(), NCseg2);
	    h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_NClay"), h2-> Fill((beam_lmom+TGT_LMOM-n_lmom).M(), NClay);

	    TLorentzVector K0_lmom = pim_lmom+pip_lmom;
	    TLorentzVector K0_lmom_CM = K0_lmom;
	    TLorentzVector beam_lmom_CM = beam_lmom;
	    beam_lmom_CM.Boost(LabToCM);
	    K0_lmom_CM.Boost(LabToCM);
	    double K0_ang = beam_lmom.Angle(K0_lmom.Vect());
	    double K0_ang_CM = beam_lmom_CM.Angle(K0_lmom_CM.Vect());


	    TVector3 nc_hitpos_sim;
	    TVector3 n_mom_sim;
	    if( simReader ){
	      MCData *mcData = simReader-> getMCData();
	      for( int i=0; i<mcData->trackSize(); i++ ){
		Track *track = mcData->track(i);
		if( track->parentTrackID()==0 && track->pdgID()==2112 ){
		  TVector3 n_mom_sim = 0.001*track-> momentum();
		  TVector3 vtx = 0.1*track-> vertex();
		  TVector3 n_mom_dir = (1./n_mom_sim.Z())*n_mom_sim;
		  nc_hitpos_sim = vtx+n_mom_dir*(1490.-vtx.Z());
		}
	      }
	    }
	  	  
	    if( simReader ){
	      if( fabs(nc_hitpos_sim.X())<160 && fabs(nc_hitpos_sim.Y())<75 ){
		TVector3 diff = nc_calc_pos- nc_hitpos_sim;
		h2 = (TH2F*)rtFile-> Get("diff_n_pos"), h2-> Fill(diff.X(), diff.Y());
	      }
	    }

	    h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN"), h1-> Fill(im_pipi);

	    if( fFPID==F_Neutron ){
	      TVector3 n_mom = (fNCpos-vtx_beam);
	      double NC_beta_c;
	      if( simReader ) NC_beta_c = (n_mom.Mag()-2.5)/((fNCtime-fT0time-beam_tof)*100.*Const);
	      else NC_beta_c = n_mom.Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
	      double nc_mom = NC_beta_c*nMass/sqrt(1-NC_beta_c*NC_beta_c);
	      n_mom.SetMag(nc_mom);
	      TLorentzVector n_lmom;
	      n_lmom.SetVectM(n_mom, nMass);
	      double mom_diff = (pipi_missing_mom - n_mom).Mag();
	      double kn_mm = (beam_lmom+TGT_LMOM-n_lmom).M();
	    }

	    if( mm_kp_pim>MM_Sp_MAX3 && mm_kp_pip>MM_Sm_MAX3 ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_woS"), h1-> Fill(im_pipi);
	      if( K0_flag2 ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_CM = K0_lmom;
		K0_lmom_CM.Boost(LabToCM);
		h1 = (TH1F*)rtFile-> Get("K0_ang"), h1-> Fill(cos(K0_ang));
		h1 = (TH1F*)rtFile-> Get("K0_ang_CM"), h1-> Fill(cos(K0_ang_CM));
		h2 = (TH2F*)rtFile-> Get("n_hitpos"), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		h1 = (TH1F*)rtFile-> Get("NC_trig");
		if( fabs(nc_calc_pos.X())<160 && fabs(nc_calc_pos.Y())<75 ) h1-> Fill(0);
		if( fabs(nc_calc_pos.X())<150 && fabs(nc_calc_pos.Y())<65 ) h1-> Fill(1);
		if( fabs(nc_calc_pos.X())<140 && fabs(nc_calc_pos.Y())<55 ) h1-> Fill(2);
		if( fabs(nc_calc_pos.X())<130 && fabs(nc_calc_pos.Y())<45 ) h1-> Fill(3);
		if( fabs(nc_calc_pos.X())<120 && fabs(nc_calc_pos.Y())<35 ) h1-> Fill(4);
		if( fabs(nc_calc_pos.X())<110 && fabs(nc_calc_pos.Y())<25 ) h1-> Fill(5);
		if( fabs(nc_calc_pos.X())<100 && fabs(nc_calc_pos.Y())<15 ) h1-> Fill(6);
		if( fabs(nc_calc_pos.X())<90  && fabs(nc_calc_pos.Y())<5  ) h1-> Fill(7);
	      
		if( simReader ){
		  if( fBVC_hit.size()==0 && fCVC_hit.size()==0 ){
		    DetectorHit *nc_hit=0;
		    h2 = (TH2F*)rtFile-> Get("n_hitpos_MC"), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		    DetectorData *dData=simReader->getDetectorData();
		    for( int i=0; i<dData->detectorHitSize(); i++ ){
		      DetectorHit *hit=dData->detectorHit(i);
		      if( hit->detectorID()==CID_NC && hit->pdg()==2112 ){
			nc_hit=hit;
		      }
		    }
		    if( nc_hit ){
		      h2 = (TH2F*)rtFile-> Get("n_hitpos_wMChit"), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		    }
		  }
		}

		HodoscopeLikeHit* nc_hit[50];
		for( int i=0; i<50; i++ ){
		  nc_hit[i] = 0;
		  for( int j=0; j<8; j++ ){
		    for( int k=0; k<fNC_hit[j].size(); k++ ){
		      double dE = fNC_hit[j][k]-> emean();
		      double time = fNC_hit[j][k]-> ctmean();
		      if( dE>2*i ){
			if( !nc_hit[i] ) nc_hit[i]=fNC_hit[j][k];
			else if( nc_hit[i]->ctmean()<time )  nc_hit[i]=fNC_hit[j][k];
		      }
		    }
		  }
		}

		h2 = (TH2F*)rtFile-> Get("NC_eff2");
		for( int i=0; i<50; i++ ){
		  if( nc_hit[i] ){
		    TVector3 nc_pos;
		    confMan-> GetGeomMapManager()-> GetGPos(CID_NC, nc_hit[i]->seg() , nc_pos);
		    nc_pos.SetY(nc_hit[i]->hitpos());

		    double fl = (nc_pos-fVtxBeam).Mag();
		    double beam_out, beam_tof;
		    ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
		    double NC_beta = fl/((nc_hit[i]->ctmean()-fT0time-beam_tof)*100.*Const);

		    if( fabs(nc_calc_pos.X())<160 && fabs(nc_calc_pos.Y())<75 ) h2-> Fill(0.0, (double)i);
		    if( fabs(nc_calc_pos.X())<150 && fabs(nc_calc_pos.Y())<65 ) h2-> Fill(1.0, (double)i);
		    if( fabs(nc_calc_pos.X())<140 && fabs(nc_calc_pos.Y())<55 ) h2-> Fill(2.0, (double)i);
		    if( fabs(nc_calc_pos.X())<130 && fabs(nc_calc_pos.Y())<45 ){
		      h2-> Fill(3.0, (double)i);
		      h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_dE%d", i)), h1-> Fill(1./NC_beta);
		    }
		    if( fabs(nc_calc_pos.X())<120 && fabs(nc_calc_pos.Y())<35 ) h2-> Fill(4.0, (double)i);
		    if( fabs(nc_calc_pos.X())<110 && fabs(nc_calc_pos.Y())<25 ) h2-> Fill(5.0, (double)i);
		    if( fabs(nc_calc_pos.X())<100 && fabs(nc_calc_pos.Y())<15 ) h2-> Fill(6.0, (double)i);
		    if( fabs(nc_calc_pos.X())<90  && fabs(nc_calc_pos.Y())<5  ) h2-> Fill(7.0, (double)i);
		    //		    h2 = (TH2F*)rtFile-> Get(Form("n_hitpos_whit_%d", i)), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		  }
		}

		if( fabs(nc_calc_pos.X())<150 && fabs(nc_calc_pos.Y())<75 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev1"), h1->Fill(6);
		}
		if( fabs(nc_calc_pos.X())<150 && fabs(nc_calc_pos.Y())<65 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev2"), h1->Fill(6);
		}
		if( fabs(nc_calc_pos.X())<140 && fabs(nc_calc_pos.Y())<55 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(4);
		  if( nc_hit[4] && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(5);
		  if( nc_hit[4] && fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(6);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(7);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 && nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(8);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev3"), h1->Fill(9);
		}
		if( fabs(nc_calc_pos.X())<130 && fabs(nc_calc_pos.Y())<45 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev4"), h1->Fill(6);
		}
		if( fabs(nc_calc_pos.X())<120 && fabs(nc_calc_pos.Y())<35 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev5"), h1->Fill(6);
		}
		if( fabs(nc_calc_pos.X())<110 && fabs(nc_calc_pos.Y())<25 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev6"), h1->Fill(6);
		}
		if( fabs(nc_calc_pos.X())<100 && fabs(nc_calc_pos.Y())<15 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev7"), h1->Fill(6);
		}
		if( fabs(nc_calc_pos.X())<90 && fabs(nc_calc_pos.Y())<5 ){
		  h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(0);
		  if( nNC()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(1);
		  if( nc_hit[4] ) h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(2);
		  if( fCVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(3);
		  if( fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(4);
		  if( fCVC_hit.size()>0 && fBVC_hit.size()>0 ) h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(5);
		  if( fFPID==F_Neutron ) h1=(TH1F*)rtFile->Get("NC_eff_ev8"), h1->Fill(6);
		}
	      
		if( fCVC_hit.size()>0 ){
		  h2 = (TH2F*)rtFile-> Get("n_hitpos_wCVC"), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		}	      
		if( fBVC_hit.size()>0 ){
		  h2 = (TH2F*)rtFile-> Get("n_hitpos_wBVC"), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		}

		if( fFPID==F_Neutron ){
		  TVector3 n_mom = (fNCpos-vtx_beam);
		  double NC_beta_c;
		  if( simReader ) NC_beta_c = (n_mom.Mag()-2.5)/((fNCtime-fT0time-beam_tof)*100.*Const);
		  else NC_beta_c = n_mom.Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
		  double nc_mom = NC_beta_c*nMass/sqrt(1-NC_beta_c*NC_beta_c);
		  n_mom.SetMag(nc_mom);
		  TLorentzVector n_lmom;
		  n_lmom.SetVectM(n_mom, nMass);
		  double mom_diff = (pipi_missing_mom - n_mom).Mag();
		  double kn_mm = (beam_lmom+TGT_LMOM-n_lmom).M();

		  h1 = (TH1F*)rtFile-> Get("K0_ang_nhit"), h1-> Fill(cos(K0_ang));
		  h1 = (TH1F*)rtFile-> Get("K0_ang_CM_nhit"), h1-> Fill(cos(K0_ang_CM));
		
		  h2 = (TH2F*)rtFile-> Get("n_mom_diff"), h2-> Fill(pipi_missing_mom.Mag(), mom_diff);

		  h2 = (TH2F*)rtFile-> Get("n_hitpos_whit"), h2-> Fill(nc_calc_pos.X(), nc_calc_pos.Y());
		  h1 = (TH1F*)rtFile-> Get("NC_eff");
		  if( fabs(nc_calc_pos.X())<160 && fabs(nc_calc_pos.Y())<80 ) h1-> Fill(0);
		  if( fabs(nc_calc_pos.X())<140 && fabs(nc_calc_pos.Y())<70 ) h1-> Fill(1);
		  if( fabs(nc_calc_pos.X())<120 && fabs(nc_calc_pos.Y())<60 ) h1-> Fill(2);
		  if( fabs(nc_calc_pos.X())<100 && fabs(nc_calc_pos.Y())<50 ) h1-> Fill(3);
		  if( fabs(nc_calc_pos.X())<80  && fabs(nc_calc_pos.Y())<40 ) h1-> Fill(4);
		  if( fabs(nc_calc_pos.X())<60  && fabs(nc_calc_pos.Y())<30 ) h1-> Fill(5);
		  if( fabs(nc_calc_pos.X())<40  && fabs(nc_calc_pos.Y())<20 ) h1-> Fill(6);
		  if( fabs(nc_calc_pos.X())<20  && fabs(nc_calc_pos.Y())<10 ) h1-> Fill(7);
		  
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN"), h1-> Fill(kn_mm);
		  if( 160<fabs(nc_calc_pos.X()) && 75<fabs(nc_calc_pos.Y()) ){
		    h1 = (TH1F*)rtFile->Get("KN_MM_pipi_wK0_wN_NC0"), h1-> Fill(kn_mm);
		  }
		  else if( 130<fabs(nc_calc_pos.X()) ){
		    h1 = (TH1F*)rtFile->Get("KN_MM_pipi_wK0_wN_NC1"), h1-> Fill(kn_mm);
		  }
		  else if( 100<fabs(nc_calc_pos.X()) ){
		    h1 = (TH1F*)rtFile->Get("KN_MM_pipi_wK0_wN_NC2"), h1-> Fill(kn_mm);
		  }
		  else if( 70<fabs(nc_calc_pos.X()) ){
		    h1 = (TH1F*)rtFile->Get("KN_MM_pipi_wK0_wN_NC3"), h1-> Fill(kn_mm);
		  }
		  else if( 40<fabs(nc_calc_pos.X()) ){
		    h1 = (TH1F*)rtFile->Get("KN_MM_pipi_wK0_wN_NC4"), h1-> Fill(kn_mm);
		  }
		  else{
		    h1 = (TH1F*)rtFile->Get("KN_MM_pipi_wK0_wN_NC5"), h1-> Fill(kn_mm);
		  }
		}
	      }
	    }

	    h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM_wN"), h2-> Fill(mm_kp_pip, mm_kp_pim);
	    if( !K0_flag2 ){
	      h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM_wN_woK0"), h2-> Fill(mm_kp_pip, mm_kp_pim);
	    }
	    if( !K0_flag3 ){
	      h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM_wN_woK0_3"), h2-> Fill(mm_kp_pip, mm_kp_pim);

	    }
	    if( !K0_flag5 ){
	      h2 = (TH2F*)rtFile-> Get("KPpip_KPpim_MM_wN_woK0_5"), h2-> Fill(mm_kp_pip, mm_kp_pim);

	      if( MM_Sp_MIN2<mm_kp_pim && mm_kp_pim<MM_Sp_MAX2 ){
		TLorentzVector Sp_lmom = beam_lmom+TGT_LMOM-pim_lmom;
		TLorentzVector Sp_lmom_CM = Sp_lmom;
		Sp_lmom_CM.Boost(LabToCM);
		h1 = (TH1F*)rtFile-> Get("Sp_ang"), h1-> Fill(Sp_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("Sp_ang_CM"), h1-> Fill(Sp_lmom_CM.CosTheta());
		if( fFPID==F_Neutron ){
		  h1 = (TH1F*)rtFile-> Get("Npip_IM_Sp"), h1-> Fill((pip_lmom+n_lmom).M());
		}
	      }

	      if( MM_Sm_MIN2<mm_kp_pip && mm_kp_pip<MM_Sm_MAX2 ){
		TLorentzVector Sm_lmom = beam_lmom+TGT_LMOM-pip_lmom;
		TLorentzVector Sm_lmom_CM = Sm_lmom;
		Sm_lmom_CM.Boost(LabToCM);
		h1 = (TH1F*)rtFile-> Get("Sm_ang"), h1-> Fill(Sm_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("Sm_ang_CM"), h1-> Fill(Sm_lmom_CM.CosTheta());
		if( fFPID==F_Neutron ){
		  h1 = (TH1F*)rtFile-> Get("Npim_IM_Sm"), h1-> Fill((pim_lmom+n_lmom).M());
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //  fill_pppim(header);
  //  std::cout<<"===== HistManwMC::fill FINISH ====="<<std::endl;
#if CALIB
  if( !simReader ) fillCalib(header);
#endif
}

void HistManwMC::finit()
{
  std::cout<<"HistManwMC::finit()"<<std::endl;
  rtFile-> cd();
  rtFile-> Write();
  std::cout<<"rtFile Write"<<std::endl;
  // This method occure abort if param file don't exist
  //  confMan-> SaveParams();
  //  confMan-> SaveCode();
  rtFile-> Close();
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

  fNC_eff_hit = 0;
  fFPID = F_Other;
  fNCbeta = DBL_MIN;
  fFLmom.SetXYZT(DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX);
  fKNpipi_K0_flag = false;
  fKNpipi_K0_SB0_flag = false;
  fKNpipi_K0_SB1_flag = false;
  fKNpipi_K0_SB2_flag = false;
  fKNpipi_K0_SB3_flag = false;
  fKNpipi_K0_SB4_flag = false;
  fKNpipi_K0_SB5_flag = false;
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

  fSm_mass_rc_flag = false;
  fSp_mass_rc_flag = false;
  fSm_mass_rc_flag2 = false;
  fSp_mass_rc_flag2 = false;
  fSm_mass_rc_flag3 = false;
  fSp_mass_rc_flag3 = false;
  fSm_mass_rc_flag4 = false;
  fSp_mass_rc_flag4 = false;
  fSm_mass_rc_flag5 = false;
  fSp_mass_rc_flag5 = false;

  if( simReader ) simReader-> clear();

  fFCflag=false;
  fFChit       = 0;
  fFDC1track   = 0;
  fFC_start    = DEFVECT;
  fFC_FDC1pos  = DEFVECT;
  fFC_hitpos   = DEFVECT;
  fFC_Angle    = DEFAULTD;
  fFC_Mom_USWK = DEFAULTD;
  fFC_Mom_TOF  = DEFAULTD;

  return flag;
}

