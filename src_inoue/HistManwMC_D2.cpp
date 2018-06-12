#include "HistManwMC.h"

#include "MyParam.h"

using namespace std;

static const int TGT_A=2;
static const TLorentzVector TGT_LMOM=D_LMOM;
//static const TLorentzVector TGT_LMOM=TLorentzVector(0., 0., 0., ThreeHeMass);

HistManwMC::HistManwMC(TFile *f, ConfMan *conf)
  : rtFile(f), confMan(conf), blMan(0), cdsMan(0), bltrackMan(0), cdstrackMan(0)
{
  Clear();
}

void HistManwMC::fill(EventHeader *header)
{
  //  rtFile-> cd();
  if( simReader ) simReader->fill(rtFile);

  TH1F *h1;
  TH2F *h2;
  TTree *tree;
  TNtuple *tup;
  if( header ){
    for( int i=0; i<20; i++ ){
      if( header->pattern(i)>0 ){
	h1 = (TH1F*)rtFile-> Get(Form("trig_%d", i)), h1-> Fill(header->pattern(i));
      }
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
  bool trigC = false;
  if( header ){
    if( header-> trigmode2(Mode_KCDH2, confMan) ) trigCDH2=true;
    if( header-> IsTrig(Trig_Charged) )trigC=true;
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
	  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtxcds, vtxb);

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
	      h2 = (TH2F*)rtFile-> Get("pim_mom_CDS_IM_pipi"), h2-> Fill(pim_lmom.Vect().Mag(), im_pipi);
	      h2 = (TH2F*)rtFile-> Get("pip_mom_CDS_IM_pipi"), h2-> Fill(pip_lmom.Vect().Mag(), im_pipi);
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
	      //	      std::cout<<" CDS p pi- evnum: "<<header->ev()<<std::endl; 
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_ppim"), h1-> Fill((pim_lmom+p_lmom).M());
	      h2 = (TH2F*)rtFile-> Get("p_mom_CDS_IM_ppim"), h2-> Fill(p_lmom.Vect().Mag(), (pim_lmom+p_lmom).M());
	      h2 = (TH2F*)rtFile-> Get("pim_mom_CDS_IM_ppim"), h2-> Fill(pim_lmom.Vect().Mag(), (pim_lmom+p_lmom).M());
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

  if( cdstrackMan->nGoodTrack()==3 ){
    if( fTrackPim.size()==2 && fTrackP.size()==1 ){
      //      cout<<" p pi- pi- event"<<endl;
      bool L_flag=false;
      CDSTrack *pim_track0=fTrackPim[0];
      CDSTrack *pim_track1=fTrackPim[1];
      CDSTrack *p_track=fTrackP[0];

      TVector3 vtx_pim0, vtx_pim1, vtx_p0, vtx_p1;
      TVector3 pim_mom0, pim_mom1, p_mom0, p_mom1;
      TLorentzVector pim_lmom0, pim_lmom1, p_lmom0, p_lmom1;
      bool flag0 = TrackTools::Calc2HelixVertex(pim_track0, p_track, vtx_pim0, vtx_p0);
      double dis0 = (vtx_pim0-vtx_p0).Mag();
      if( !pim_track0-> GetMomentum(vtx_pim0, pim_mom0, true, true) ) flag0=false;
      if( !p_track-> GetMomentum(vtx_p0, p_mom0, true, true) ) flag0=false;
      pim_lmom0.SetVectM(pim_mom0, piMass);
      p_lmom0.SetVectM(p_mom0, pMass);
      double ppim_im0 = (p_lmom0+pim_lmom0).M();

      bool flag1 = TrackTools::Calc2HelixVertex(pim_track1, p_track, vtx_pim1, vtx_p1);
      double dis1 = (vtx_pim1-vtx_p1).Mag();
      if( !pim_track1-> GetMomentum(vtx_pim1, pim_mom1, true, true) ) flag1=false;
      if( !p_track-> GetMomentum(vtx_p1, p_mom1, true, true) ) flag1=false;
      pim_lmom1.SetVectM(pim_mom1, piMass);
      p_lmom1.SetVectM(p_mom1, pMass);
      double ppim_im1 = (p_lmom1+pim_lmom1).M();

      TVector3 vtx_b0, vtx_b1, vtx_b_pim0, vtx_b_pim1;
      pim_track1-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_b1, vtx_b_pim1);
      pim_track0-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_b0, vtx_b_pim0);
    
      if( flag1 && flag0 ){
	//	cout<<" IM0 : "<<ppim_im0<<"  IM1 : "<<ppim_im1<<endl;
	if( GeomTools::GetID(vtx_b1)==CID_Fiducial || GeomTools::GetID(vtx_b0)==CID_Fiducial ){
	  //*** dumping ***//
	  if( header ){
	    std::cout<<header->ev()<<"  ppimpim_event"<<std::endl;
	  }
	  //*** dumping ***//

	  if( dis0<dis1 ) h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim"), h2->Fill(ppim_im0, ppim_im1);
	  else h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim"), h2->Fill(ppim_im1, ppim_im0);
	  if( trigCDH2 ){
	    if( dis0<dis1 ) h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim_CDH2"), h2->Fill(ppim_im0, ppim_im1);
	    else h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim_CDH2"), h2->Fill(ppim_im1, ppim_im0);
	  }
	  if( trigC ){
	    if( dis0<dis1 ) h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim_C"), h2->Fill(ppim_im0, ppim_im1);
	    else h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim_C"), h2->Fill(ppim_im1, ppim_im0);
	  }
	}
	
	bool L_flag0=false, L_flag1=false;
	if( L_MIN<ppim_im0 && ppim_im0<L_MAX ) L_flag0=true;
	if( L_MIN<ppim_im1 && ppim_im1<L_MAX ) L_flag1=true;
	
	TVector3 vtxb, vtxcds, vtxL, pim_mom;
	TLorentzVector L_lmom, pim_lmom;
	bool flag=true;
	if( !L_flag1 && !L_flag0 ) flag=false;
	
	if( L_flag0 || L_flag1 ){
	  if( L_flag0 && L_flag1 ){
	    double diff0 = fabs(ppim_im1-L_peak);
	    double diff1 = fabs(ppim_im1-L_peak);
	    if( diff0<diff1 ) L_flag1=false;
	    else L_flag0=false;
	  }

	  if( L_flag0 ){
	    if( !pim_track1-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxb, vtxcds) ) flag=false;
	    if( !pim_track1-> GetMomentum(vtxcds, pim_mom) ) flag=false;
	    pim_lmom.SetVectM(pim_mom, piMass);
	    L_lmom = p_lmom0+pim_lmom0;
	    vtxL = 0.5*(vtx_p0+vtx_pim0);
	  }
	  else if( L_flag1 ){
	    if( !pim_track0-> GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtxb, vtxcds) ) flag=false;
	    if( !pim_track0-> GetMomentum(vtxcds, pim_mom) ) flag=false;
	    pim_lmom.SetVectM(pim_mom, piMass);
	    L_lmom = p_lmom1+pim_lmom1;
	    vtxL = 0.5*(vtx_p1+vtx_pim1);
	  }
       
	  if( flag && GeomTools::GetID(vtxcds)==CID_Fiducial ){
	    h2 =  (TH2F*)rtFile-> Get("CDS_ppimpim_vtxL_XY"), h2-> Fill(vtxL.X(), vtxL.Y());
	    h2 =  (TH2F*)rtFile-> Get("CDS_ppimpim_vtxL_ZX"), h2-> Fill(vtxL.Z(), vtxL.X());
	    h2 =  (TH2F*)rtFile-> Get("CDS_ppimpim_vtxL_ZY"), h2-> Fill(vtxL.Z(), vtxL.Y());

	    TLorentzVector Lpim_missing_lmom = fBeamLmom+D_LMOM-L_lmom-pim_lmom;
	    TVector3 LabToCM = -(fBeamLmom+D_LMOM).BoostVector();
	    h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim"), h1-> Fill((L_lmom+pim_lmom).M());
	    h1 = (TH1F*)rtFile-> Get("KLpim_MM"), h1-> Fill(Lpim_missing_lmom.M());   
	    h2 = (TH2F*)rtFile-> Get("CDS_IM_Lpim_KLpim_MM"), h2-> Fill((L_lmom+pim_lmom).M(), Lpim_missing_lmom.M());
	    if( trigCDH2 ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_CDH2"), h1-> Fill((L_lmom+pim_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KLpim_MM_CDH2"), h1-> Fill(Lpim_missing_lmom.M());   
	      h2 = (TH2F*)rtFile-> Get("CDS_IM_Lpim_KLpim_MM_CDH2"), h2-> Fill((L_lmom+pim_lmom).M(), Lpim_missing_lmom.M());
	    }
	    if( trigC ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_C"), h1-> Fill((L_lmom+pim_lmom).M());
	      h1 = (TH1F*)rtFile-> Get("KLpim_MM_C"), h1-> Fill(Lpim_missing_lmom.M());   
	      h2 = (TH2F*)rtFile-> Get("CDS_IM_Lpim_KLpim_MM_C"), h2-> Fill((L_lmom+pim_lmom).M(), Lpim_missing_lmom.M());
	    }

	    if( simReader ){
	      simReader-> fillLpim(rtFile);
	    }

	    if( 0.89<Lpim_missing_lmom.M() && Lpim_missing_lmom.M()<0.99 ){
	      TLorentzVector p_lmom_CM;
	      p_lmom_CM.SetVectM(Lpim_missing_lmom.Vect(), pMass);
	      TLorentzVector beam_lmom_CM=fBeamLmom;
	      p_lmom_CM.Boost(LabToCM);
	      beam_lmom_CM.Boost(LabToCM);

	      double p_ang = cos(beam_lmom_CM.Angle(p_lmom_CM.Vect()));
	      double p_ang_lab = cos(fBeamLmom.Angle(Lpim_missing_lmom.Vect()));
	      h1 = (TH1F*)rtFile-> Get("CDS_Lpim_p_ang"), h1-> Fill(p_ang);
	      h1 = (TH1F*)rtFile-> Get("CDS_Lpim_p_ang_lab"), h1-> Fill(p_ang_lab);
	      h2 = (TH2F*)rtFile-> Get("CDS_Lpim_p_ang_CDS_IM_Lpim"), h2-> Fill(p_ang, (L_lmom+pim_lmom).M());
	      if( trigCDH2 ){
		h1 = (TH1F*)rtFile-> Get("CDS_Lpim_p_ang_CDH2"), h1-> Fill(p_ang);
		h1 = (TH1F*)rtFile-> Get("CDS_Lpim_p_ang_lab_CDH2"), h1-> Fill(p_ang_lab);
		h2 = (TH2F*)rtFile-> Get("CDS_Lpim_p_ang_CDS_IM_Lpim_CDH2"), h2-> Fill(p_ang, (L_lmom+pim_lmom).M());
	      }
	      if( trigC ){
		h1 = (TH1F*)rtFile-> Get("CDS_Lpim_p_ang_C"), h1-> Fill(p_ang);
		h1 = (TH1F*)rtFile-> Get("CDS_Lpim_p_ang_lab_C"), h1-> Fill(p_ang_lab);
		h2 = (TH2F*)rtFile-> Get("CDS_Lpim_p_ang_CDS_IM_Lpim_C"), h2-> Fill(p_ang, (L_lmom+pim_lmom).M());
	      
		h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_C_mm_p"), h1-> Fill((L_lmom+pim_lmom).M());
		if( header ){
		  cout<<"===== L pi- \"p\" event : "<<header->ev()<<" ====="<<endl;
		  cout<<"> nBVC : "<<fBVC_hit.size()<<endl;
		  cout<<"> nCVC : "<<fCVC_hit.size()<<endl;
		  cout<<"> nPC  : "<<fPC_hit.size()<<endl;
		  cout<<"> nFDC1 : "<<bltrackMan->ntrackFDC1()<<endl;
		}
		if( fCVC_hit.size()>0 ||  fPC_hit.size()>0 ){
		  h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_C_mm_p_wCVCPC"), h1-> Fill((L_lmom+pim_lmom).M());
		}
		if( bltrackMan->ntrackFDC1()==1 ){
		  h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_C_mm_p_wFDC1"), h1-> Fill((L_lmom+pim_lmom).M());
		  if( fFCflag ){
		    h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_C_mm_p_FC"), h1-> Fill((L_lmom+pim_lmom).M());
		    h1 = (TH1F*)rtFile-> Get("Lpim_C_mm_p_mass2"), h1-> Fill(fFC_mass2);
		    if( fFPID==F_Proton ){
		      h1 = (TH1F*)rtFile-> Get("CDS_IM_Lpim_C_mm_p_fp"), h1-> Fill((L_lmom+pim_lmom).M());
		    }
		  }
		}
	      }
	      if( simReader ){
		simReader-> fillLpim_p(rtFile);
	      }
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
	
	int NClay = (fNCseg-1)/16+1;
	int NCseg2 = fNCseg%16;
	if( NCseg2==0 ) NCseg2=16;
	h2 = (TH2F*)rtFile-> Get("KN_MM_NCseg"), h2-> Fill(kn_lmom.M(), NCseg2);
	h2 = (TH2F*)rtFile-> Get("KN_MM_NClay"), h2-> Fill(kn_lmom.M(), NClay);

	h1 = (TH1F*)rtFile-> Get("KN_MM"), h1->Fill(kn_lmom.M());
	if( fTrackPim.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pim"), h1->Fill(kn_lmom.M());
	if( fTrackPip.size()>0 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pip"), h1->Fill(kn_lmom.M());
	if( fTrackKm.size()>0 )  h1 = (TH1F*)rtFile-> Get("KN_MM_km"), h1->Fill(kn_lmom.M());
	if( fTrackP.size()>0 )   h1 = (TH1F*)rtFile-> Get("KN_MM_p"), h1->Fill(kn_lmom.M());
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
  if( fTrackP.size()==1 && fTrackPim.size()==1 ){
    TVector3 p_vtx, pim_vtx;
    bool ppim_flag = TrackTools::Calc2HelixVertex(fTrackP[0], fTrackPim[0], p_vtx, pim_vtx);
    TVector3 vtx_ppim = 0.5*(p_vtx+pim_vtx);
    
    TVector3 p_mom, pim_mom;
    bool p_mom_flag = fTrackP[0]->GetMomentum(p_vtx, p_mom, true, true);
    bool pim_mom_flag = fTrackPim[0]->GetMomentum(pim_vtx, pim_mom, true, true);
    
    TVector3 vtxb, vtxcds;
    double dltmp=0;
    double dca;
    if( ppim_flag && p_mom_flag && pim_mom_flag ){
      MathTools::LineToLine(vtx_ppim, (p_mom+pim_mom).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtxcds, vtxb);
    
      double beam_out, beam_tof;
      ELossTools::CalcElossBeamTGeo(fT0pos, vtxb, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      TVector3 beam_mom = fTrackBPC->GetMomDir();
      beam_mom.SetMag(beam_out);
      TLorentzVector beam_lmom;
      beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);
      
      if( fFPID==F_Neutron && trigNC_flag && fBeamPID==Beam_Kaon ){
	TVector3 n_mom = (fNCpos-vtxb);
	double NC_beta_c;
	if( simReader ) NC_beta_c=((fNCpos-vtxb).Mag()-2.5)/((fNCtime-fT0time-beam_tof)*100.*Const);
	else NC_beta_c=(fNCpos-vtxb).Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
	double nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	n_mom.SetMag(nc_mom);
	TLorentzVector n_lmom;
	n_lmom.SetVectM(n_mom, nMass);
	
	TLorentzVector p_lmom, pim_lmom;
	p_lmom.SetVectM(p_mom, pMass), pim_lmom.SetVectM(pim_mom, piMass);
	TLorentzVector ppim_lmom=p_lmom+pim_lmom;

	TLorentzVector kn_lmom = beam_lmom+TGT_LMOM-n_lmom;
	TLorentzVector knppim_lmom = beam_lmom+TGT_LMOM-n_lmom-p_lmom-pim_lmom;

	bool L_flag=false;
	if( L_MIN<ppim_lmom.M() && ppim_lmom.M()<L_MAX ) L_flag=true;
	
	bool fiducial_flag=false;
	if( GeomTools::GetID(vtxb)==CID_Fiducial ) fiducial_flag=true;

	if( fiducial_flag && cdstrackMan->nGoodTrack()==2 ){
	  h1 = (TH1F*)rtFile-> Get("CDS_IM_ppim_wN"), h1-> Fill(ppim_lmom.M());
	  h2 = (TH2F*)rtFile-> Get("KN_MM_KNppim_MM"), h2-> Fill(kn_lmom.M(), knppim_lmom.M());
	  h2 = (TH2F*)rtFile-> Get("KN_MM_KNppim_MM2"), h2-> Fill(kn_lmom.M(), knppim_lmom.M2());
	  if( L_flag ){
	    TVector3 ppim_mom = p_mom+pim_mom;
	    TLorentzVector l_lmom;
	    l_lmom.SetVectM(ppim_mom, lMass);
	    TLorentzVector knl_lmom = beam_lmom+TGT_LMOM-n_lmom-l_lmom;
	    h1 = (TH1F*)rtFile-> Get("KN_MM_wL"), h1-> Fill(kn_lmom.M());
	    h2 = (TH2F*)rtFile-> Get("KN_MM_KNppim_MM2_wL"), h2-> Fill(kn_lmom.M(), knppim_lmom.M2());
	    h2 = (TH2F*)rtFile-> Get("KN_MM_KNL_MM"), h2-> Fill(kn_lmom.M(), knl_lmom.M());
	    h2 = (TH2F*)rtFile-> Get("KN_MM_KNL_MM2"), h2-> Fill(kn_lmom.M(), knl_lmom.M2());
	    if( knppim_lmom.M()>0 ){
	      h2 = (TH2F*)rtFile-> Get("KN_MM_KNppim_MM_wL"), h2-> Fill(kn_lmom.M(), knppim_lmom.M());
	    }
	  }
	}
      }
    }
  }

  if( fTrackPim.size()==1 && fTrackPip.size()==1 ){  
    double pim_mass2, pip_mass2, pim_mom2, pip_mom2, pim_time, pip_time;
    int pim_seg, pip_seg;
    fTrackPim[0]-> GetCDHHit(cdsMan, pim_seg, pim_time);
    fTrackPip[0]-> GetCDHHit(cdsMan, pip_seg, pip_time);
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
    bool vtx_pim_flag=false;
    if( dis_beam_pim<dis_beam_pip ){
      if( GeomTools::GetID(vtx_beam_pim)==CID_Fiducial ) fiducial_flag=true; 
      ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam_pim, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      vtx_c = vtx_pim_beam;
      min_dca=dis_beam_pim;
      vtx_pim_flag=true;
    }
    else{ 
      if( GeomTools::GetID(vtx_beam_pip)==CID_Fiducial ) fiducial_flag=true; 
      ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam_pip, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      vtx_c = vtx_pip_beam;
      min_dca=dis_beam_pip;
      vtx_pim_flag=false;
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
    }

    if( fFPID==F_Neutron && trigNC_flag ){
      if( header && fBeamPID==Beam_Kaon ){
	//	std::cout<<"  pi+ pi- tag"<<std::endl;
	h1_ER-> Fill(12);
	if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(12);
      }
      else{
	h1_ER-> Fill(12);
      }
   
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
      	//***************************//
	//*** d(K-, n) w/ pi+ pi- ***//
	//***************************//
	if( simReader ){
	  simReader-> fillpipi(rtFile);
	  //	  if( simReader->getReactionData()->ReactionID()==821 ) simReader->printReaction();
	}

	TVector3 n_mom = (fNCpos-vtx_c);
	double NC_beta_c;
	if( simReader ) NC_beta_c=((fNCpos-vtx_c).Mag()-2.5)/((fNCtime-fT0time-beam_tof)*100.*Const);
	else NC_beta_c=(fNCpos-vtx_c).Mag()/((fNCtime-fT0time-beam_tof)*100.*Const);
	double nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	n_mom.SetMag(nc_mom);
	TLorentzVector n_lmom;
	n_lmom.SetVectM(n_mom, nMass);

	TLorentzVector knpipi_lmom = beam_lmom+TGT_LMOM-n_lmom-pim_lmom_beam-pip_lmom_beam;
	TLorentzVector missing_n_lmom;
	missing_n_lmom.SetVectM(missing_n_lmom.Vect(), nMass);
	const double missing_ene = beam_lmom.E()+TGT_LMOM.E()-n_lmom.E()-pim_lmom.E()-pip_lmom.E()-missing_n_lmom.E();
	const TLorentzVector kn_lmom = beam_lmom+TGT_LMOM-n_lmom;
	const TVector3 kn_mom = kn_lmom.Vect();

	const double im_npim = (n_lmom+pim_lmom_beam).M();
	const double im_npip = (n_lmom+pip_lmom_beam).M();
	const double mm_knpim = (beam_lmom+TGT_LMOM-n_lmom-pim_lmom_beam).M();
	const double mm_knpip = (beam_lmom+TGT_LMOM-n_lmom-pip_lmom_beam).M();
	const double mm_knpipi = (beam_lmom+TGT_LMOM-n_lmom-pim_lmom_beam-pip_lmom_beam).M();
	const double im_pipi = (pim_lmom+pip_lmom).M();
	const double mm_kn = (fBeamLmom+TGT_LMOM-n_lmom).Mag();
	const double im_npipi = (n_lmom+pim_lmom_beam+pip_lmom_beam).M();
	const double mm_kp_npipi = (beam_lmom+P_LMOM-n_lmom-pim_lmom_beam-pip_lmom_beam).M();
	const double mm_kp_pipi  = (beam_lmom+P_LMOM-pim_lmom_beam-pip_lmom_beam).M();
	TLorentzVector MM_pipi_lmom = beam_lmom+P_LMOM-pim_lmom-pip_lmom;
	TLorentzVector MM_N_lmom;
	MM_N_lmom.SetVectM(MM_pipi_lmom.Vect(), nMass);
	const double Sm_mass_rc = (MM_N_lmom+pim_lmom).M();
	const double Sp_mass_rc = (MM_N_lmom+pip_lmom).M();

	TLorentzVector knpim_lmom = beam_lmom+TGT_LMOM-n_lmom-pim_lmom;
	TLorentzVector knpip_lmom = beam_lmom+TGT_LMOM-n_lmom-pip_lmom;
	const TLorentzVector tot_lmom = beam_lmom+TGT_LMOM;
	const TVector3 LabToCM = -tot_lmom.BoostVector();
	const TVector3 LabToY = -kn_lmom.BoostVector();
	const TLorentzVector Kp_lmom = TGT_LMOM + beam_lmom;
	const TVector3 LabToKp = -Kp_lmom.BoostVector();

	TLorentzVector kn_lmom_CM = beam_lmom+TGT_LMOM-n_lmom;
	TLorentzVector knpim_lmom_CM = beam_lmom+TGT_LMOM-n_lmom-pim_lmom;
	TLorentzVector knpip_lmom_CM = beam_lmom+TGT_LMOM-n_lmom-pip_lmom;
	kn_lmom_CM.Boost(LabToCM);
	knpim_lmom_CM.Boost(LabToCM);
	knpip_lmom_CM.Boost(LabToCM);
	const TVector3 kn_mom_CM = kn_lmom_CM.Vect();
	const TVector3 knpim_mom_CM = knpim_lmom_CM.Vect();
	const TVector3 knpip_mom_CM = knpip_lmom_CM.Vect();
	const double KNpim_cos_CM = kn_mom_CM.Dot(knpim_mom_CM)/(kn_mom_CM.Mag()*knpim_mom_CM.Mag());
	const double KNpip_cos_CM = kn_mom_CM.Dot(knpip_mom_CM)/(kn_mom_CM.Mag()*knpip_mom_CM.Mag());

	TLorentzVector knpim_lmom_Y = beam_lmom+TGT_LMOM-n_lmom-pim_lmom;
	TLorentzVector knpip_lmom_Y = beam_lmom+TGT_LMOM-n_lmom-pip_lmom;
	TLorentzVector pim_lmom_Y = pim_lmom;
	TLorentzVector pip_lmom_Y = pip_lmom;
	knpim_lmom_Y.Boost(LabToY);
	knpip_lmom_Y.Boost(LabToY);
	pim_lmom_Y.Boost(LabToY);
	pip_lmom_Y.Boost(LabToY);
	const TVector3 knpim_mom_Y = knpim_lmom_Y.Vect();
	const TVector3 knpip_mom_Y = knpip_lmom_Y.Vect();
	const TVector3 pim_mom_Y = pim_lmom_Y.Vect();
	const TVector3 pip_mom_Y = pip_lmom_Y.Vect();
	const double KNpim_cos_Y = knpim_mom_Y.Dot(kn_mom)/(kn_mom.Mag()*knpim_mom_Y.Mag());
	const double KNpip_cos_Y = knpip_mom_Y.Dot(kn_mom)/(kn_mom.Mag()*knpip_mom_Y.Mag());
	const double pim_cos_Y = pim_mom_Y.Dot(kn_mom)/(kn_mom.Mag()*pim_mom_Y.Mag());
	const double pip_cos_Y = pip_mom_Y.Dot(kn_mom)/(kn_mom.Mag()*pip_mom_Y.Mag());

	KinFitMan kinFitMan;
	if( simReader ) kinFitMan.pipiMMSA_n_spec(beam_lmom, TGT_LMOM, pip_lmom, pim_lmom, n_lmom, rtFile, simReader);
	else kinFitMan.pipiMMSA_n_spec(beam_lmom, TGT_LMOM, pip_lmom, pim_lmom, n_lmom, rtFile);
	TLorentzVector mmsa_n_lmom = kinFitMan.getReacLmom();
	const double npim_im_mmsa = (mmsa_n_lmom+pim_lmom).M();
	const double npip_im_mmsa = (mmsa_n_lmom+pip_lmom).M();
	const double mm_knn = (N_LMOM+beam_lmom-n_lmom).M();

        int NClay = (fNCseg-1)/16+1;
        int NCseg2 = fNCseg%16;
        if( NCseg2==0 ) NCseg2=16;

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

	if( Npipi_K0_SB0_MIN<im_pipi && im_pipi<Npipi_K0_SB0_MAX ) fKNpipi_K0_SB0_flag=true;
	if( Npipi_K0_SB1_MIN<im_pipi && im_pipi<Npipi_K0_SB1_MAX ) fKNpipi_K0_SB1_flag=true;
	if( Npipi_K0_SB0_MIN<im_pipi && im_pipi<Npipi_K0_SB0_MAX ) fKNpipi_K0_SB2_flag=true;
	if( Npipi_K0_SB1_MIN<im_pipi && im_pipi<Npipi_K0_SB1_MAX ) fKNpipi_K0_SB3_flag=true;
	if( Npipi_K0_SB0_MIN<im_pipi && im_pipi<Npipi_K0_SB0_MAX ) fKNpipi_K0_SB4_flag=true;
	if( Npipi_K0_SB1_MIN<im_pipi && im_pipi<Npipi_K0_SB1_MAX ) fKNpipi_K0_SB5_flag=true;

	if( KNpim_MM_Sp_MIN<mm_knpim && mm_knpim<KNpim_MM_Sp_MAX ) fKNpim_Sp_flag=true;
	if( KNpip_MM_Sm_MIN<mm_knpip && mm_knpip<KNpip_MM_Sm_MAX ) fKNpip_Sm_flag=true;
	if( Npipi_Sm_MIN<im_npim && im_npim<Npipi_Sm_MAX ) fKNpipi_Sm_flag=true;
	if( Npipi_Sp_MIN<im_npip && im_npip<Npipi_Sp_MAX ) fKNpipi_Sp_flag=true;

	if( MM_Npim_MIN <Sm_mass_rc && Sm_mass_rc<MM_Npim_MAX  ) fSm_mass_rc_flag = true;
	if( MM_Npip_MIN <Sp_mass_rc && Sp_mass_rc<MM_Npip_MAX  ) fSp_mass_rc_flag = true;
	if( MM_Npim_MIN2<Sm_mass_rc && Sm_mass_rc<MM_Npim_MAX2 ) fSm_mass_rc_flag2 = true;
	if( MM_Npip_MIN2<Sp_mass_rc && Sp_mass_rc<MM_Npip_MAX2 ) fSp_mass_rc_flag2 = true;
	if( MM_Npim_MIN3<Sm_mass_rc && Sm_mass_rc<MM_Npim_MAX3 ) fSm_mass_rc_flag3 = true;
	if( MM_Npip_MIN3<Sp_mass_rc && Sp_mass_rc<MM_Npip_MAX3 ) fSp_mass_rc_flag3 = true;
	if( MM_Npim_MIN4<Sm_mass_rc && Sm_mass_rc<MM_Npim_MAX4 ) fSm_mass_rc_flag4 = true;
	if( MM_Npip_MIN4<Sp_mass_rc && Sp_mass_rc<MM_Npip_MAX4 ) fSp_mass_rc_flag4 = true;
	if( MM_Npim_MIN5<Sm_mass_rc && Sm_mass_rc<MM_Npim_MAX5 ) fSm_mass_rc_flag5 = true;
	if( MM_Npip_MIN5<Sp_mass_rc && Sp_mass_rc<MM_Npip_MAX5 ) fSp_mass_rc_flag5 = true;

	if( Npipi_N_MIN<mm_knpipi && mm_knpipi<Npipi_N_MAX ) fKNpipi_N_flag=true;
	if( Npipi_N_MIN0<mm_knpipi && mm_knpipi<Npipi_N_MAX0 ) fKNpipi_N0_flag=true;
	if( Npipi_N_MIN1<mm_knpipi && mm_knpipi<Npipi_N_MAX1 ) fKNpipi_N1_flag=true;
	if( Npipi_N_MIN2<mm_knpipi && mm_knpipi<Npipi_N_MAX2 ) fKNpipi_N2_flag=true;
      	
	if( fBeamPID==Beam_Kaon ){
	  h2 = (TH2F*)rtFile-> Get("Vtx_XY_npipi"), h2-> Fill(vtx_c.X(), vtx_c.Y());
	  h2 = (TH2F*)rtFile-> Get("Vtx_ZX_npipi"), h2-> Fill(vtx_c.Z(), vtx_c.X());
	  h2 = (TH2F*)rtFile-> Get("Vtx_ZY_npipi"), h2-> Fill(vtx_c.Z(), vtx_c.Y());

	  TVector3 vtx_pipi = 0.5*(vtx_pim+vtx_pip);
	  h1 = (TH1F*)rtFile-> Get("DCA_pipi_beam_npipi_woF"), h1-> Fill(dca_pipi_beam);

	  h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi_n_woF"), h2-> Fill(vtx_pipi.X(), vtx_pipi.Y());
	  h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi_n_woF"), h1-> Fill(vtx_pipi.Z());
	  if( fKNpipi_N_flag ){
	    h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi_n_wN_woF"), h2-> Fill(vtx_pipi.X(), vtx_pipi.Y());
	    h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi_n_wN_woF"), h1-> Fill(vtx_pipi.Z());
	  }
	  if( 1.0<mm_knpipi && mm_knpipi<1.05 ){
	    h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi_n_tail_woF"), h2-> Fill(vtx_pipi.X(), vtx_pipi.Y());
	    h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi_n_tail_woF"), h1-> Fill(vtx_pipi.Z());
	  }

	  if( cdstrackMan->nGoodTrack()==2 ){
	    h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_woF"), h1->Fill(mm_knpipi);
	    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_woF"), h1-> Fill(mm_kn);
	  
	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile-> Get("CDS_chi2_npipin");
	      if( vtx_pim_flag ){ h2-> Fill(fTrackPim[0]->Chi(), fTrackPip[0]->Chi()); }
	      else{ h2-> Fill(fTrackPip[0]->Chi(), fTrackPim[0]->Chi()); }
	    }

	    if( MyTools::IsTarget0(vtx_c, confMan) && MyTools::IsTarget0(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx0"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx0"), h1-> Fill(mm_kn);
	    }
	  }
	
	  if( cdstrackMan->nGoodTrack()==2 ){
	    h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_woF"), h1->Fill(mm_knpipi);
	    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_woF"), h1-> Fill(mm_kn);
	    if( MyTools::IsTarget0(vtx_c, confMan) && MyTools::IsTarget0(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx0"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx0"), h1-> Fill(mm_kn);
	    }
	  }
	
	  //	  if( fiducial_flag && GeomTools::GetID(vtx_pipi)==CID_Fiducial && cdstrackMan->nGoodTrack()==2 ){
	  if( fiducial_flag && cdstrackMan->nGoodTrack()==2 ){
	    MyTools::checkCDC2(cdsMan, fTrackPim[0], fTrackPip[0]);
	    //	  if( fiducial_flag ){
	    //*** set data to NpipiData
	    // fNpipiData-> setBeamLmom(beam_lmom);
	    // fNpipiData-> setFNLmom(n_lmom);
	    // fNpipiData-> setPimLmom(pim_lmom);
	    // fNpipiData-> setPipLmom(pip_lmom);
	    // fNpipiData-> setVtxCDS(vtx_cds);
	    // fNpipiData-> setVtxBeam(vtx_beam);
	    // fNpipiData-> setVtxPim(vtx_pim);
	    // fNpipiData-> setVtxPim(vtx_pip);
	    // fNpipiData-> setVtxPimBeam(vtx_pim_beam);
	    // fNpipiData-> setVtxPipBeam(vtx_pip_beam);
	    // fNpipiData-> setVtxBeamPim(vtx_beam_pim);
	    // fNpipiData-> setVtxBeamPip(vtx_beam_pip);
	    // tree = (TTree*)rtFile-> Get("NpipiTree"), tree-> Fill();      

	    // 2015/12/15 Add **************************************************************//
	    //	    tup = (TNtuple*)rtFile-> Get("tupNpipi");
	    //	    tup-> Fill(mm_kn, mm_knpipi, mm_knpim, mm_knpip, im_pipi, im_npim, im_npip, 
	    //		       n_lmom.Vect().Mag(), pim_lmom.Vect().Mag(), pip_lmom.Vect().Mag(), (pim_lmom+pip_lmom).Vect().Mag(), knpipi_lmom.CosTheta());
	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile-> Get("CDS_chi2_npipin");
	      if( vtx_pim_flag ){ h2-> Fill(fTrackPim[0]->Chi(), fTrackPip[0]->Chi()); }
	      else{ h2-> Fill(fTrackPip[0]->Chi(), fTrackPim[0]->Chi()); }
	    }

	    h2 = (TH2F*)rtFile->Get("KN_MM_pipi_NCseg"), h2-> Fill(mm_kn, NCseg2);
	    h2 = (TH2F*)rtFile->Get("KN_MM_pipi_NClay"), h2-> Fill(mm_kn, NClay);

	    h2 = (TH2F*)rtFile->Get("KNpipi_MM_pipi_mom"), h2-> Fill(mm_knpipi, (pim_lmom+pip_lmom).Vect().Mag());

	    // 2015/12/12 Add **************************************************************//
	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile->Get("KN_MM_pipi_wN_NCseg"), h2-> Fill(mm_kn, NCseg2);
	      h2 = (TH2F*)rtFile->Get("KN_MM_pipi_wN_NClay"), h2-> Fill(mm_kn, NClay);

	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_Npim_IM"), h2->Fill(mm_kn, im_npim);
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_Npip_IM"), h2->Fill(mm_kn, im_npip);
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_CDS_IM_pipi"), h2->Fill(mm_kn, im_pipi);
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_KNpip_MM"), h2-> Fill(mm_kn, mm_knpip);
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_KNpim_MM"), h2-> Fill(mm_kn, mm_knpim);
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_n_theta"), h2-> Fill(mm_kn, knpipi_lmom.CosTheta());
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_mmN_mom"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_pim_mom"), h2-> Fill(mm_kn, pim_lmom.Vect().Mag());
	      // h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_pip_mom"), h2-> Fill(mm_kn, pip_lmom.Vect().Mag());

	      TVector3 LabToKP = -(beam_lmom+P_LMOM).BoostVector();
	      if( fKNpipi_K0_flag ){
		TLorentzVector K0_lmom_KP=pim_lmom+pip_lmom;
		K0_lmom_KP.Boost(LabToKP);
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wK0_K0_theta"), h2-> Fill(mm_kn, K0_lmom_KP.CosTheta());
	      }
	      if( fKNpipi_Sm_flag ){
		TLorentzVector pip_lmom_KP = pip_lmom;
		pip_lmom_KP.Boost(LabToKP);
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wSm_pip_theta"), h2-> Fill(mm_kn, pip_lmom_KP.CosTheta());
	      }
	      if( fKNpipi_Sp_flag ){
		TLorentzVector pim_lmom_KP = pim_lmom;
		pim_lmom_KP.Boost(LabToKP);
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wSp_pim_theta"), h2-> Fill(mm_kn, pim_lmom_KP.CosTheta());
	      }

	      if( !fKNpipi_K0_flag && !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ){
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_wN_KNpip_MM"), h2-> Fill(mm_kn, mm_knpip);
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_wN_KNpim_MM"), h2-> Fill(mm_kn, mm_knpim);
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_wN_n_theta"), h2-> Fill(mm_kn, knpipi_lmom.CosTheta());
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_wN_mmN_mom"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_wN_pim_mom"), h2-> Fill(mm_kn, pim_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_wN_pip_mom"), h2-> Fill(mm_kn, pip_lmom.Vect().Mag());
	      }

	      if( !fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woK0"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSm"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSp"), h1-> Fill(mm_kn);

	      if( !fKNpipi_K0_flag25 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woK0_25"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sm_flag25 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSm_25"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sp_flag25 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSp_25"), h1-> Fill(mm_kn);

	      if( !fKNpipi_K0_flag3 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woK0_3"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sm_flag3 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSm_3"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sp_flag3 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSp_3"), h1-> Fill(mm_kn);

	      if( !fKNpipi_K0_flag4 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woK0_4"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sm_flag4 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSm_4"), h1-> Fill(mm_kn);
	      if( !fKNpipi_Sp_flag4 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woSp_4"), h1-> Fill(mm_kn);

	      if( !fKNpipi_K0_flag25 && !fKNpipi_Sm_flag25 && !fKNpipi_Sp_flag25 ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woAll25"), h1-> Fill(mm_kn);
		//		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_woAll25"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      }
	      if( !fKNpipi_K0_flag3 && !fKNpipi_Sm_flag3 && !fKNpipi_Sp_flag3 ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woAll3"), h1-> Fill(mm_kn);
		//		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_woAll3"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      }
	      if( !fKNpipi_K0_flag4 && !fKNpipi_Sm_flag4 && !fKNpipi_Sp_flag4 ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN_woAll4"), h1-> Fill(mm_kn);
		//		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_woAll4"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      }
	    }
	    // 2015/12/12 Add **************************************************************//
	  
	    if( MyTools::IsTarget1(vtx_c, confMan) && MyTools::IsTarget1(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx1"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx1"), h1-> Fill(mm_kn);
	    }
	    if( MyTools::IsTarget2(vtx_c, confMan) && MyTools::IsTarget2(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx2"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx2"), h1-> Fill(mm_kn);
	    }
	    if( MyTools::IsTarget3(vtx_c, confMan) && MyTools::IsTarget3(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx3"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx3"), h1-> Fill(mm_kn);
	    }
	    if( MyTools::IsTarget4(vtx_c, confMan) && MyTools::IsTarget4(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx4"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx4"), h1-> Fill(mm_kn);
	    }
	    if( MyTools::IsTarget5(vtx_c, confMan) && MyTools::IsTarget5(vtx_pipi, confMan) ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_vtx5"), h1->Fill(mm_knpipi);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_vtx5"), h1-> Fill(mm_kn);
	    }

	    h1 = (TH1F*)rtFile-> Get("DCA_pipi_beam_npipi"), h1-> Fill(dca_pipi_beam);

	    h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi_n"), h2-> Fill(vtx_pipi.X(), vtx_pipi.Y());
	    h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi_n"), h1-> Fill(vtx_pipi.Z());
	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi_n_wN"), h2-> Fill(vtx_pipi.X(), vtx_pipi.Y());
	      h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi_n_wN"), h1-> Fill(vtx_pipi.Z());
	    }
	    if( 1.0<mm_knpipi && mm_knpipi<1.05 ){
	      h2 = (TH2F*)rtFile-> Get("Vtx_XY_pipi_n_tail"), h2-> Fill(vtx_pipi.X(), vtx_pipi.Y());
	      h1 = (TH1F*)rtFile-> Get("Vtx_Z_pipi_n_tail"), h1-> Fill(vtx_pipi.Z());
	    }
#if 0
	    if( simReader ){
	      ReactionData *reacData = simReader->getReactionData();
	      std::cout<<"===== K- d -> n pi+ pi- Event in Fiducial ====="<<std::endl;
	      std::cout<<"> Reaction ID : "<<reacData->ReactionID()<<std::endl;
	      std::cout<<">>> Analyzied Data"<<std::endl;
	      std::cout<<"> d(K-, n) : "<<mm_kn<<" [GeV/c2]"<<std::endl;
	      std::cout<<"> d(K-, n pi+ pi-) : "<<mm_knpipi<<" [GeV/c2]"<<std::endl;
	      std::cout<<">>> NC dumpping "<<std::endl;
	      simReader->trace(CID_NC, fNCseg, true);
	      std::cout<<"> pi+   ana mom : "<<pip_lmom.Vect().Mag()<<" [GeV/c]"<<std::endl;
	      simReader->trace(fTrackPip[0], cdsMan, true);
	      std::cout<<"> pi- dummpping : "<<pim_lmom.Vect().Mag()<<" [GeV/c]"<<std::endl;
	      simReader->trace(fTrackPim[0], cdsMan, true);
	    }
#endif
	    if( simReader ){
	      rtFile-> cd();
	      simReader-> fillNpipi(rtFile, cdsMan, beam_lmom, fTrackPim[0], pim_lmom, fTrackPip[0], pip_lmom, fNCseg, n_lmom, mm_knpipi);
	      simReader-> fillHist_knEl_data(fBeamLmom, n_lmom, pim_lmom, pip_lmom, vtx_pim, vtx_pip, 
					     vtx_pim_beam, vtx_pip_beam, vtx_beam_pim, vtx_beam_pip,
					     fKNpipi_N_flag, fKNpipi_K0_flag, fKNpipi_Sm_flag, fKNpipi_Sp_flag);
	      if( fKNpipi_N_flag ){
		if( !fKNpipi_K0_flag ) simReader-> fill_woK0(rtFile); 
		if( !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ) simReader-> fill_woSf(rtFile);
		simReader-> fillMM_N(rtFile);
		if( !fKNpipi_K0_flag ) simReader-> fillMM_N_woK0(rtFile);
		if( !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ) simReader-> fillMM_N_woSf(rtFile);
		if( !fKNpipi_K0_flag && !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ){
		  simReader-> fillMM_N_wo(rtFile);
		  if( fKNpim_Sp_flag || fKNpip_Sm_flag ) simReader-> fillCharged(rtFile);
		  if( !fSm_mass_rc_flag && !fSp_mass_rc_flag ) simReader-> fillMM_N_woRcS(rtFile);
		  if( fKNpim_Sp_flag && !fKNpip_Sm_flag ) simReader-> fillSp(rtFile);
		  if( fKNpip_Sm_flag && !fKNpim_Sp_flag ) simReader-> fillSm(rtFile);
		}
	      }

	      simReader-> fillKN_MM(rtFile, mm_kn);
	      simReader-> fillKNpipiMM(rtFile, mm_knpipi);

	      if( fKNpipi_K0_flag ) simReader-> fillKNpipiMM_wK0(rtFile, mm_knpipi);
	      if( fKNpipi_Sm_flag ) simReader-> fillKNpipiMM_wSm(rtFile, mm_knpipi);
	      if( fKNpipi_Sp_flag ) simReader-> fillKNpipiMM_wSp(rtFile, mm_knpipi);
	      if( !fKNpipi_K0_flag && !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ) simReader-> fillKNpipiMM_woAll(rtFile, mm_knpipi);
	    }

	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pim_npipi"), h2-> Fill(pim_mass2, pim_mom2);
	    h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pip_npipi"), h2-> Fill(pip_mass2, pip_mom2);

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

	    if( fKNpipi_N_flag ){
	      h1= (TH1F*)rtFile-> Get("KNpipi_n_event"), h1-> Fill(0);
	      if( fKNpipi_K0_flag ) h1->Fill(1);
	      if( fKNpipi_Sm_flag ) h1->Fill(2);
	      if( fKNpipi_Sp_flag ) h1->Fill(3);
	      if( fKNpipi_K0_flag && fKNpipi_Sm_flag ) h1-> Fill(4);
	      if( fKNpipi_K0_flag && fKNpipi_Sp_flag ) h1-> Fill(5);
	      if( fKNpipi_Sp_flag && fKNpipi_Sm_flag ) h1-> Fill(6);
	      if( fKNpipi_K0_flag && fKNpipi_Sm_flag && fKNpipi_Sp_flag ) h1-> Fill(7);

	      h1= (TH1F*)rtFile-> Get("KNpipi_n_event25"), h1-> Fill(0);
	      if( fKNpipi_K0_flag25 ) h1->Fill(1);
	      if( fKNpipi_Sm_flag25 ) h1->Fill(2);
	      if( fKNpipi_Sp_flag25 ) h1->Fill(3);
	      if( fKNpipi_K0_flag25 && fKNpipi_Sm_flag25 ) h1-> Fill(4);
	      if( fKNpipi_K0_flag25 && fKNpipi_Sp_flag25 ) h1-> Fill(5);
	      if( fKNpipi_Sp_flag25 && fKNpipi_Sm_flag25 ) h1-> Fill(6);
	      if( fKNpipi_K0_flag25 && fKNpipi_Sm_flag25 && fKNpipi_Sp_flag25 ) h1-> Fill(7);

	      h1= (TH1F*)rtFile-> Get("KNpipi_n_event3"), h1-> Fill(0);
	      if( fKNpipi_K0_flag3 ) h1->Fill(1);
	      if( fKNpipi_Sm_flag3 ) h1->Fill(2);
	      if( fKNpipi_Sp_flag3 ) h1->Fill(3);
	      if( fKNpipi_K0_flag3 && fKNpipi_Sm_flag3 ) h1-> Fill(4);
	      if( fKNpipi_K0_flag3 && fKNpipi_Sp_flag3 ) h1-> Fill(5);
	      if( fKNpipi_Sp_flag3 && fKNpipi_Sm_flag3 ) h1-> Fill(6);
	      if( fKNpipi_K0_flag3 && fKNpipi_Sm_flag3 && fKNpipi_Sp_flag3 ) h1-> Fill(7);

	      h1= (TH1F*)rtFile-> Get("KNpipi_n_event4"), h1-> Fill(0);
	      if( fKNpipi_K0_flag4 ) h1->Fill(1);
	      if( fKNpipi_Sm_flag4 ) h1->Fill(2);
	      if( fKNpipi_Sp_flag4 ) h1->Fill(3);
	      if( fKNpipi_K0_flag4 && fKNpipi_Sm_flag4 ) h1-> Fill(4);
	      if( fKNpipi_K0_flag4 && fKNpipi_Sp_flag4 ) h1-> Fill(5);
	      if( fKNpipi_Sp_flag4 && fKNpipi_Sm_flag4 ) h1-> Fill(6);
	      if( fKNpipi_K0_flag4 && fKNpipi_Sm_flag4 && fKNpipi_Sp_flag4 ) h1-> Fill(7);

	    }

	    if( 1.00<mm_knpipi && mm_knpipi<1.05 ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_test"), h1-> Fill(im_pipi);
	      h1 = (TH1F*)rtFile-> Get("Npim_IM_test"), h1-> Fill(im_npim);
	      h1 = (TH1F*)rtFile-> Get("Npip_IM_test"), h1-> Fill(im_npip);
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_test"), h2-> Fill(mm_knpim, mm_knpip);
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pim_test"), h2-> Fill(pim_mass2, pim_mom2);
	      h2 = (TH2F*)rtFile-> Get("CDS_mass2_mom_pip_test"), h2-> Fill(pip_mass2, pip_mom2);
	      h1 = (TH1F*)rtFile-> Get("pim_dmom_test"), h1-> Fill(pim_lmom_beam.Mag()-pim_mom2);
	      h1 = (TH1F*)rtFile-> Get("pip_dmom_test"), h1-> Fill(pip_lmom_beam.Mag()-pip_mom2);
	    }
	  
	    if( header ){
	      //	    std::cout<<">>> Event Number : "<<header->ev()<<" n pi+ pi- hit fiducial"<<std::endl;
	      //	      TNtuple *tup;

	      //	      tup = (TNtuple*)rtFile->Get("tupNpipi"), tup-> Fill(header->ev());
	    }
	    if( header ){
	      //	      std::cout<<"  pi+ pi- tag on Target"<<std::endl;
	      h1_ER-> Fill(13);
	      if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(13);
	    }
	    else{
	      h1_ER-> Fill(13);
	    }

	    h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN"), h1-> Fill(im_pipi);
	    if( fKNpipi_N_flag ){
	      h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wN"), h1-> Fill(im_pipi);
	      if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wN_wSm"), h1-> Fill(im_pipi);
	      if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wN_wSp"), h1-> Fill(im_pipi);
	      if( fKNpipi_Sm_flag || fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wN_wS"), h1-> Fill(im_pipi);
	    }
	    if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wSm"), h1-> Fill(im_pipi);
	    if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wSp"), h1-> Fill(im_pipi);
	    if( fKNpipi_Sm_flag || fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi_wN_wSf"), h1-> Fill(im_pipi);
	  
	    h1 = (TH1F*)rtFile-> Get("Npim_IM_pip"), h1-> Fill(im_npim);
	    h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_MMSA"), h1-> Fill(npim_im_mmsa);
	    if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_rcSm"), h1-> Fill(im_npim);

	    if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wK0"), h1-> Fill(im_npim);
	    if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wSp"), h1-> Fill(im_npim);
	    if( fKNpipi_K0_flag || fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wSp_K0"), h1-> Fill(im_npim);
	    if( fKNpipi_N_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wN"), h1-> Fill(im_npim);
	      if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_wN_rcSm"), h1-> Fill(im_npim);

	      if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wN_wK0"), h1-> Fill(im_npim);
	      if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wN_wSp"), h1-> Fill(im_npim);
	      if( fKNpipi_K0_flag || fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npim_IM_pip_wN_wK0_Sp"), h1-> Fill(im_npim);
	    }
	 
	    h1 = (TH1F*)rtFile-> Get("Npip_IM_pim"), h1-> Fill(im_npip);
	    h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_MMSA"), h1-> Fill(npip_im_mmsa);
	    if( fSp_mass_rc_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npip_IM_rcSp"), h1-> Fill(im_npip);
	    }
	    if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wK0"), h1-> Fill(im_npip);
	    if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wSm"), h1-> Fill(im_npip);
	    if( fKNpipi_K0_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wSm_K0"), h1-> Fill(im_npip);
	    if( fKNpipi_N_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wN"), h1-> Fill(im_npip);
	      if( fSp_mass_rc_flag ){
		h1 = (TH1F*)rtFile-> Get("Npip_IM_wN_rcSp"), h1-> Fill(im_npip);
	      }
	      if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wN_wK0"), h1-> Fill(im_npip);
	      if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wN_wSm"), h1-> Fill(im_npip);
	      if( fKNpipi_K0_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npip_IM_pim_wN_wK0_Sm"), h1-> Fill(im_npip);
	    }

	    
	    h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM"), h2-> Fill(im_npim, im_npip);
	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM_wN"), h2-> Fill(im_npim, im_npip);
	      if( fKNpipi_K0_flag ) h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM_wN_wK0"), h2-> Fill(im_npim, im_npip);
	      else h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM_wN_woK0"), h2-> Fill(im_npim, im_npip);
	    }
	    if( fKNpipi_K0_flag ){
	      h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wK0_KNpipi_MM"), h2->Fill(mm_kn, mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM_wK0"), h2-> Fill(im_npim, im_npip);
	    }
	    else h2 = (TH2F*)rtFile-> Get("Npim_Npip_IM_woK0"), h2-> Fill(im_npim, im_npip);

	    if( simReader ){
	      Track *n_track = simReader->trace(CID_NC, fNCseg);
	      ReactionData *reacData = simReader->getReactionData();
	      TLorentzVector beam_lmom_MC = 0.001*reacData->GetInitParticle(0);
	      TLorentzVector n_lmom_MC;
	      n_lmom_MC.SetVectM(0.001*n_track->momentum(), nMass);
	      double mm_kn_MC = (beam_lmom_MC+TGT_LMOM-n_lmom_MC).Mag();

	      Track *nc_hit_track = simReader->trace(CID_NC, fNCseg);
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_diff_MC"), h1-> Fill(mm_kn_MC-mm_kn);
	      if( fKNpipi_K0_flag )h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_diff_MC"), h1-> Fill(mm_kn_MC-mm_kn);
	      if( simReader-> getReactionData()->ReactionID()==735 ){
		//		if( nc_hit_track->trackID()==3 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_n_spec_MC"), h1-> Fill(mm_kn);
	      }

	      int NCstatus = simReader->NChitNeutronStatus(fNCseg);
	    }

	    //*** dumping ***//
	    if( header ){
	      std::cout<<header->ev()<<"  npipi_event"<<std::endl;
	    }
	  
	    //*** dumping ***//
	
	    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi"), h1->Fill(mm_kn);
	    h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_KNpipi_MM"), h2->Fill(mm_kn, mm_knpipi);
	    h1 = (TH1F*)rtFile-> Get("KNpipi_MM"), h1->Fill(mm_knpipi);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_KP_pipi_MM"), h2-> Fill(mm_knpipi, mm_kp_pipi);

	    h2 = (TH2F*)rtFile-> Get("Npim_IM_MM_Npim_IM_no_cut"), h2-> Fill(im_npim, Sm_mass_rc);
	    h2 = (TH2F*)rtFile-> Get("Npip_IM_MM_Npip_IM_no_cut"), h2-> Fill(im_npip, Sp_mass_rc);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_MM_Npim_IM_no_cut"), h2-> Fill(mm_knpipi, Sm_mass_rc);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_MM_Npip_IM_no_cut"), h2-> Fill(mm_knpipi, Sp_mass_rc);

	    h1 = (TH1F*)rtFile-> Get("MM_Npim_IM_no_cut"), h1-> Fill(Sm_mass_rc);
	    h1 = (TH1F*)rtFile-> Get("MM_Npip_IM_no_cut"), h1-> Fill(Sp_mass_rc);

	    if( mm_knpipi<1.05 ){
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_MM_Npim_IM"), h2-> Fill(mm_knpipi, Sm_mass_rc);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_MM_Npip_IM"), h2-> Fill(mm_knpipi, Sp_mass_rc);

	      h1 = (TH1F*)rtFile-> Get("MM_Npim_IM"), h1-> Fill(Sm_mass_rc);
	      if( !fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("MM_Npim_IM_woSm"), h1-> Fill(Sm_mass_rc);
	      else h1 = (TH1F*)rtFile-> Get("MM_Npim_IM_wSm"), h1-> Fill(Sm_mass_rc);

	      h1 = (TH1F*)rtFile-> Get("MM_Npip_IM"), h1-> Fill(Sp_mass_rc);
	      if( !fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("MM_Npip_IM_woSp"), h1-> Fill(Sp_mass_rc);
	      else h1 = (TH1F*)rtFile-> Get("MM_Npip_IM_wSp"), h1-> Fill(Sp_mass_rc);
	    }
	    if( fSm_mass_rc_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_rcSm"), h1-> Fill(mm_knpipi);
	      if( !fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_rcSm_woSm"), h1-> Fill(mm_knpipi);
	    }
	    if( fSp_mass_rc_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_rcSp"), h1-> Fill(mm_knpipi);
	      if( !fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_rcSp_woSp"), h1-> Fill(mm_knpipi);
	    
	    }	   
	    if( fSm_mass_rc_flag || fSp_mass_rc_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_rcS"), h1-> Fill(mm_knpipi);
	      if( !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_rcS_woS"), h1-> Fill(mm_knpipi);
	    }

	    if( !fSm_mass_rc_flag || !fSp_mass_rc_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wo_rcS"), h1-> Fill(mm_knpipi);
	    }
	    if( !fSm_mass_rc_flag2 || !fSp_mass_rc_flag2 ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wo_rcS2"), h1-> Fill(mm_knpipi);
	    }
	    if( !fSm_mass_rc_flag3 || !fSp_mass_rc_flag3 ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wo_rcS3"), h1-> Fill(mm_knpipi);
	    }
	    if( !fSm_mass_rc_flag4 || !fSp_mass_rc_flag4 ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wo_rcS4"), h1-> Fill(mm_knpipi);
	    }
	    if( !fSm_mass_rc_flag5 || !fSp_mass_rc_flag5 ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wo_rcS5"), h1-> Fill(mm_knpipi);
	    }

	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_missing_ene"), h2->Fill(mm_knpipi, missing_ene);

	    if( fKNpipi_Sp_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wSp"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wSp_KNpipi_MM"), h2->Fill(mm_kn, mm_knpipi);
	      if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wSp_rcSp"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_missing_ene_wSp"), h2->Fill(mm_knpipi, missing_ene);
	    }
	    if( fKNpipi_Sm_flag ){
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wSm"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wSm_KNpipi_MM"), h2->Fill(mm_kn, mm_knpipi);
	      if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wSm_rcSm"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_missing_ene_wSm"), h2->Fill(mm_knpipi, missing_ene);
	    }
	  
	    h2 = (TH2F*)rtFile-> Get("KNpipi_KN_MM"), h2->Fill(mm_knpipi, mm_kn);
	    if( fKNpipi_K0_flag ) h2 = (TH2F*)rtFile-> Get("KNpipi_KN_MM_wK0"), h2->Fill(mm_knpipi, mm_kn);
	    if( fKNpipi_Sm_flag ) h2 = (TH2F*)rtFile-> Get("KNpipi_KN_MM_wSm"), h2->Fill(mm_knpipi, mm_kn);
	    if( fKNpipi_Sp_flag ) h2 = (TH2F*)rtFile-> Get("KNpipi_KN_MM_wSp"), h2->Fill(mm_knpipi, mm_kn);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_KNpim_MM"), h2->Fill(mm_knpipi, mm_knpim);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_KNpip_MM"), h2->Fill(mm_knpipi, mm_knpip);	    
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_beam_mom"), h2-> Fill(mm_knpipi, beam_lmom.Vect().Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_NC_seg"), h2->Fill(mm_knpipi, fNCseg);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_NC_dE"), h2->Fill(mm_knpipi, fNCdE);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_n_mom"), h2->Fill(mm_knpipi, n_lmom.Vect().Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_nc_beta"), h2-> Fill(mm_knpipi, fNCbeta);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_nc_overbeta"), h2-> Fill(mm_knpipi, 1/fNCbeta);
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_pim_mom"), h2-> Fill(mm_knpipi, pim_lmom_beam.Vect().Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_pip_mom"), h2-> Fill(mm_knpipi, pip_lmom_beam.Vect().Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_KNpim_mom"), h2-> Fill(mm_knpipi, kn_pim_mom.Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_KNpip_mom"), h2-> Fill(mm_knpipi, kn_pip_mom.Mag());

	    h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	    h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM"), h2-> Fill(mm_knpim, mm_knpip);
	    
	    h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
	    h1 = (TH1F*)rtFile-> Get("MinDCA_npipi"), h1-> Fill(min_dca);
	    h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip"), h2-> Fill(dis_beam_pim, dis_beam_pip);
 	    
	    h1 = (TH1F*)rtFile-> Get("KNn_MM_npipi"), h1-> Fill(mm_knn);
	    if( mm_knpipi<1.05 ){
	      h1 = (TH1F*)rtFile-> Get("KNn_MM_npipi_cut"), h1-> Fill(mm_knn);
	      if( fKNpipi_K0_flag ) h1 = (TH1F*)rtFile-> Get("KNn_MM_npipi_cut_wK0"), h1-> Fill(mm_knn);
	      if( fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KNn_MM_npipi_cut_wSm"), h1-> Fill(mm_knn);
	      if( fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("KNn_MM_npipi_cut_wSp"), h1-> Fill(mm_knn);
	      if( !fKNpipi_K0_flag && !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("KNn_MM_npipi_cut_woAll"), h1-> Fill(mm_knn);
	      }
	    }
       
	    if( fKNpipi_N_flag ){
	      h2 = (TH2F*)rtFile-> Get("KN_MM_KNpim_MM_wN"), h2-> Fill(mm_kn, mm_knpim);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_KNpip_MM_wN"), h2-> Fill(mm_kn, mm_knpip);

	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wN"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get(Form("KN_MM_pipi_wN_T0%d", fT0_hit[0]->seg())), h1->Fill(mm_kn);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wN_dE"), h2->Fill(mm_kn, fNC_eff_hit->emean());
	      h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());

	      h1 = (TH1F*)rtFile-> Get("KNpipi_mom"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      if( fKNpipi_K0_flag ){
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0"), h1-> Fill(knpipi_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_wK0"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      }
	      if( fKNpipi_Sm_flag ){
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wSm"), h1-> Fill(knpipi_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_wSm"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      }
	      if( fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wSp"), h1-> Fill(knpipi_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_wSp"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
	      }
	      if( !fKNpipi_K0_flag && !fKNpipi_Sm_flag && !fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_woAll"), h1-> Fill(knpipi_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpipi_mom_woAll"), h2-> Fill(mm_kn, knpipi_lmom.Vect().Mag());
		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpim_MM_woAll_wN"), h2-> Fill(mm_kn, mm_knpim);
		h2 = (TH2F*)rtFile-> Get("KN_MM_KNpip_MM_woAll_wN"), h2-> Fill(mm_kn, mm_knpip);
		if( knpipi_lmom.Vect().Mag()>0.2 ){
		  h2 = (TH2F*)rtFile-> Get("KN_MM_KNpim_MM_woAll_wN_200"), h2-> Fill(mm_kn, mm_knpim);
		  h2 = (TH2F*)rtFile-> Get("KN_MM_KNpip_MM_woAll_wN_200"), h2-> Fill(mm_kn, mm_knpip);
		}
	      }

	      if( simReader) simReader-> fillKN_MM_wN(rtFile, mm_kn);
	      if( fKNpipi_Sp_flag || fKNpipi_Sm_flag || fKNpipi_K0_flag){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wAll_wN"), h1->Fill(mm_kn);
	      }
	      if( fKNpipi_Sp_flag ){
		TLorentzVector pim_lmom_Kp = pim_lmom;
		pim_lmom_Kp.Boost(LabToKp);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wSp_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("pim_ang_npimSp"), h1-> Fill(pim_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("pim_ang_npimSp_Kp"), h1-> Fill(pim_lmom_Kp.CosTheta());
		if( !fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wSp_woSm_wN"), h1->Fill(mm_kn);
	      }
	      if( fSp_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_rcSp_wN"), h1->Fill(mm_kn);

	      if( fKNpipi_Sm_flag ){
		TLorentzVector pip_lmom_Kp = pip_lmom;
		pip_lmom_Kp.Boost(LabToKp);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wSm_wN"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("pip_ang_npipSm"), h1-> Fill(pip_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("pip_ang_npipSm_Kp"), h1-> Fill(pip_lmom_Kp.CosTheta());
	       	if( !fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wSm_woSp_wN"), h1->Fill(mm_kn);
	      }
	      if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_rcSm_wN"), h1->Fill(mm_kn);
	    }

	    if( fKNpipi_K0_flag ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0"), h1->Fill(mm_kn);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_wK0"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_wK0"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_missing_ene_wK0"), h2->Fill(mm_knpipi, missing_ene);
	      if( fKNpipi_N_flag ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN"),    h1->Fill(mm_kn);
		h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_wK0_wN_dE"), h2->Fill(mm_kn, fNC_eff_hit->emean());

		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp"), h1-> Fill(K0_lmom_Kp.CosTheta());
	      }
	    }
	    if( fKNpipi_K0_SB0_flag ){
	      if( fKNpipi_N_flag ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_SB0"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_SB0"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp_SB0"), h1-> Fill(K0_lmom_Kp.CosTheta());
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0_SB0"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      }
	    }
	    if( fKNpipi_K0_SB1_flag ){
	      if( fKNpipi_N_flag ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_SB1"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_SB1"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp_SB1"), h1-> Fill(K0_lmom_Kp.CosTheta());
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0_SB1"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      }
	    }
	    if( fKNpipi_K0_SB2_flag ){
	      if( fKNpipi_N_flag ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_SB2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_SB2"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp_SB2"), h1-> Fill(K0_lmom_Kp.CosTheta());
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0_SB2"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      }
	    }
	    if( fKNpipi_K0_SB3_flag ){
	      if( fKNpipi_N_flag ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_SB3"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_SB3"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp_SB3"), h1-> Fill(K0_lmom_Kp.CosTheta());
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0_SB3"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      }
	    }
	    if( fKNpipi_K0_SB4_flag ){
	      if( fKNpipi_N_flag ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_SB4"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_SB4"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp_SB4"), h1-> Fill(K0_lmom_Kp.CosTheta());
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0_SB4"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      }
	    }
	    if( fKNpipi_K0_SB5_flag ){
	      if( fKNpipi_N_flag ){
		TLorentzVector K0_lmom = pim_lmom+pip_lmom;
		TLorentzVector K0_lmom_Kp = K0_lmom;
		K0_lmom.Boost(LabToKp);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_wN_SB5"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_SB5"), h1-> Fill(K0_lmom.CosTheta());
		h1 = (TH1F*)rtFile-> Get("K0_ang_npipin_Kp_SB5"), h1-> Fill(K0_lmom_Kp.CosTheta());
		h1 = (TH1F*)rtFile-> Get("KNpipi_mom_wK0_SB5"), h1-> Fill(knpipi_lmom.Vect().Mag());
	      }
	    }
	  
	    if( fKNpipi_K0_flag25 ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_25"), h1->Fill(mm_kn);
	    if( fKNpipi_K0_flag3  ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_3"), h1->Fill(mm_kn);
	    if( fKNpipi_K0_flag4  ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wK0_4"), h1->Fill(mm_kn);

	    if( fKNpipi_Sm_flag  ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wSm"), h1-> Fill(mm_kn);
	    if( fKNpipi_Sp_flag  ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wSp"), h1-> Fill(mm_kn);
	    if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_rcSm"), h1-> Fill(mm_kn);
	    if( fSp_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_rcSp"), h1-> Fill(mm_kn);

	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS"), h1->Fill(mm_kn);
	      if( fKNpipi_N_flag ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_wN"), h1->Fill(mm_kn);
	      }
	    }
	  
	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_25"), h1->Fill(mm_kn);
	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_3"), h1->Fill(mm_kn);
	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_wS_4"), h1->Fill(mm_kn);

	  
	    if( !fKNpipi_K0_flag && !fKNpipi_Sp_flag && !fKNpipi_Sm_flag ){
	      if( fKNpipi_N_flag ){
		h1 = (TH1F*)rtFile-> Get("KNpim_ang_CM_npipi"), h1-> Fill(KNpim_cos_CM);
		h1 = (TH1F*)rtFile-> Get("KNpip_ang_CM_npipi"), h1-> Fill(KNpip_cos_CM);

		h1 = (TH1F*)rtFile-> Get("KNpim_ang_piS_npipi"), h1-> Fill(KNpim_cos_Y);
		h1 = (TH1F*)rtFile-> Get("KNpip_ang_piS_npipi"), h1-> Fill(KNpip_cos_Y);
		h1 = (TH1F*)rtFile-> Get("pim_ang_piS_npipi"), h1-> Fill(pim_cos_Y);
		h1 = (TH1F*)rtFile-> Get("pip_ang_piS_npipi"), h1-> Fill(pip_cos_Y);

		if( fKNpim_Sp_flag ){
		  h1 = (TH1F*)rtFile-> Get("KNpim_ang_CM_npipi_mmSp"),  h1-> Fill(KNpim_cos_CM);
		  h1 = (TH1F*)rtFile-> Get("KNpim_ang_piS_npipi_mmSp"), h1-> Fill(KNpim_cos_Y);
		  h1 = (TH1F*)rtFile-> Get("pim_ang_piS_npipi_mmSp"),   h1-> Fill(pim_cos_Y);
		  if( mm_kn<1.44 ){
		    h1 = (TH1F*)rtFile-> Get("KNpim_ang_piS_npipi_mmSp_144"), h1-> Fill(KNpim_cos_Y);
		    h1 = (TH1F*)rtFile-> Get("pim_ang_piS_npipi_mmSp_144"),   h1-> Fill(pim_cos_Y);
		  }
		  if( !fKNpip_Sm_flag ){
		    h1 = (TH1F*)rtFile-> Get("KNpim_ang_CM_npipi_mmSp_woSm"),  h1-> Fill(KNpim_cos_CM);
		    h1 = (TH1F*)rtFile-> Get("KNpim_ang_piS_npipi_mmSp_woSm"), h1-> Fill(KNpim_cos_Y);
		    h1 = (TH1F*)rtFile-> Get("pim_ang_piS_npipi_mmSp_woSm"),   h1-> Fill(pim_cos_Y);
		    if( mm_kn<1.44 ){
		      h1 = (TH1F*)rtFile-> Get("KNpim_ang_piS_npipi_mmSp_woSm_144"), h1-> Fill(KNpim_cos_Y);
		      h1 = (TH1F*)rtFile-> Get("pim_ang_piS_npipi_mmSp_woSm_144"),   h1-> Fill(pim_cos_Y);
		    }
		  }
		}
		if( fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KNpip_ang_CM_npipi_mmSm"),  h1-> Fill(KNpip_cos_CM);
		  h1 = (TH1F*)rtFile-> Get("KNpip_ang_piS_npipi_mmSm"), h1-> Fill(KNpip_cos_Y);
		  h1 = (TH1F*)rtFile-> Get("pip_ang_piS_npipi_mmSm"),   h1-> Fill(pip_cos_Y);
		  if( mm_kn<1.44 ){
		    h1 = (TH1F*)rtFile-> Get("KNpip_ang_piS_npipi_mmSm_144"), h1-> Fill(KNpip_cos_Y);
		    h1 = (TH1F*)rtFile-> Get("pip_ang_piS_npipi_mmSm_144"),   h1-> Fill(pip_cos_Y);
		  }
		  if( !fKNpim_Sp_flag ){
		    h1 = (TH1F*)rtFile-> Get("KNpip_ang_CM_npipi_mmSm_woSp"), h1-> Fill(KNpip_cos_CM);
		    h1 = (TH1F*)rtFile-> Get("KNpip_ang_piS_npipi_mmSm_woSp"), h1-> Fill(KNpip_cos_Y);
		    h1 = (TH1F*)rtFile-> Get("pip_ang_piS_npipi_mmSm_woSp"),   h1-> Fill(pip_cos_Y);
		    if( mm_kn<1.44 ){
		      h1 = (TH1F*)rtFile-> Get("KNpip_ang_piS_npipi_mmSm_woSp_144"), h1-> Fill(KNpip_cos_Y);
		      h1 = (TH1F*)rtFile-> Get("pip_ang_piS_npipi_mmSm_woSp_144"),   h1-> Fill(pip_cos_Y);
		    }
		  }
		}
	      }
	    
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll"), h1->Fill(mm_kn);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_woAll_MM_Npim_IM"), h2->Fill(mm_kn, Sm_mass_rc);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_woAll_MM_Npip_IM"), h2->Fill(mm_kn, Sp_mass_rc);
	      if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_rcSm"), h1->Fill(mm_kn);
	      if( fSp_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_rcSp"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll"), h1->Fill(mm_knpipi);
	      if( fSm_mass_rc_flag || fSp_mass_rc_flag ){
		h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_rcS"), h1->Fill(mm_knpipi);
	      }
	      if( fSm_mass_rc_flag2 || fSp_mass_rc_flag2 ){
		h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_rcS2"), h1->Fill(mm_knpipi);
	      }
	      if( !fSm_mass_rc_flag2 && !fSp_mass_rc_flag2 ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_woS2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_woS2"), h1->Fill(mm_knpipi);
	      }
	      if( !fSm_mass_rc_flag && !fSp_mass_rc_flag  ){
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll2"), h1->Fill(mm_kn);
		h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll2"), h1->Fill(mm_knpipi);
	      }

	      for( int mm=130; mm<160; mm++ ){
		if( 0.01*mm<mm_kn && mm_kn<0.01*(mm+1) ) h1 = (TH1F*)rtFile-> Get(Form("KNpipi_MM_woAll_%d", mm)), h1->Fill(mm_knpipi);
	      }
	      if( 1.6<mm_kn ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_160"), h1->Fill(mm_knpipi);
	      if( mm_kn<kpMass+pMass ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_down"), h1->Fill(mm_knpipi);
	      else h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_up"), h1->Fill(mm_knpipi);
	    
	      h2 = (TH2F*)rtFile-> Get("KN_MM_pipi_woAll_KNpipi_MM"), h2-> Fill(mm_kn, mm_knpipi);

	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_woAll_DCA"), h2->Fill(mm_knpipi, min_dca);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_missing_ene_woAll"), h2->Fill(mm_knpipi, missing_ene);

	      if( min_dca<3.) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_DCA_3"), h1->Fill(mm_knpipi);
	      if( min_dca<2.) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_DCA_2"), h1->Fill(mm_knpipi);
	      if( min_dca<1.) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_DCA_1"), h1->Fill(mm_knpipi);

	      if( fBeamPion ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_wPi"), h1->Fill(mm_knpipi);
	      if( fBHD_hit.size()==1 ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_BHD1hit"), h1->Fill(mm_knpipi);
	   
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_n_mom_woAll"), h2->Fill(mm_knpipi, n_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_NC_seg_woAll"), h2->Fill(mm_knpipi, fNCseg);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_NC_dE_woAll"), h2->Fill(mm_knpipi, fNCdE);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_beam_mom_woAll"), h2-> Fill(mm_knpipi, beam_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_nc_beta_woAll"), h2-> Fill(mm_knpipi, fNCbeta);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_nc_overbeta_woAll"), h2-> Fill(mm_knpipi, 1/fNCbeta);
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll"), h2-> Fill(mm_knpim, mm_knpip);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_KN_MM_woAll"), h2->Fill(mm_knpipi, mm_kn);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_KNpim_MM_woAll"), h2->Fill(mm_knpipi, mm_knpim);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_KNpip_MM_woAll"), h2->Fill(mm_knpipi, mm_knpip);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_pim_mom_woAll"), h2->Fill(mm_knpipi, pim_lmom_beam.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_pip_mom_woAll"), h2->Fill(mm_knpipi, pip_lmom_beam.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_KNpim_mom_woAll"), h2-> Fill(mm_knpipi, kn_pim_mom.Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_KNpip_mom_woAll"), h2-> Fill(mm_knpipi, kn_pip_mom.Mag());
	      if( !fSm_mass_rc_flag && !fSp_mass_rc_flag  ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll2"), h2-> Fill(mm_knpim, mm_knpip);
	      }	    

	      h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip_woAll"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
	      h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip_woAll"), h2-> Fill(dis_beam_pim, dis_beam_pip);
	      
	      double chi2_pim = fTrackPim[0]-> Chi();
	      double chi2_pip = fTrackPip[0]-> Chi();
	      h1 = (TH1F*)rtFile-> Get("CDC_pim_chi2_Npipi"), h1-> Fill(chi2_pim);
	      h1 = (TH1F*)rtFile-> Get("CDC_pip_chi2_Npipi"), h1-> Fill(chi2_pip);
	      if( chi2_pim<50 && chi2_pip<50 ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_chi2_50"), h1->Fill(mm_knpipi);
	      if( chi2_pim<30 && chi2_pip<30 ) h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll_chi2_30"), h1->Fill(mm_knpipi);
	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KN_MM_woAll_wN_MM_Npim_IM"), h2->Fill(mm_kn, Sp_mass_rc);
		h2 = (TH2F*)rtFile-> Get("KN_MM_woAll_wN_MM_Npim_IM"), h2->Fill(mm_kn, Sp_mass_rc);
		if( fSm_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_rcSm"), h1->Fill(mm_kn);
		if( fSp_mass_rc_flag ) h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_rcSp"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSp_down"), h1-> Fill(mm_kn);
		}
		if( fKNpim_Sp_flag && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSp_up"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim<KNpim_MM_Sp_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSm_down"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim>KNpim_MM_Sp_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSm_up"), h1-> Fill(mm_kn);
		}

		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wA"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wB"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wC"), h1-> Fill(mm_kn);
		}
		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wD"), h1-> Fill(mm_kn);
		}
	      }

	      if( fKNpipi_N0_flag ){
		if( fKNpim_Sp_flag && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wSp_down"), h1-> Fill(mm_kn);
		}
		if( fKNpim_Sp_flag && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wSp_up"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim<KNpim_MM_Sp_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wSm_down"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim>KNpim_MM_Sp_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wSm_up"), h1-> Fill(mm_kn);
		}

		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wA"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wB"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wC"), h1-> Fill(mm_kn);
		}
		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wD"), h1-> Fill(mm_kn);
		}
	      }

	      if( fKNpipi_N1_flag ){
		if( fKNpim_Sp_flag && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wSp_down"), h1-> Fill(mm_kn);
		}
		if( fKNpim_Sp_flag && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wSp_up"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim<KNpim_MM_Sp_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wSm_down"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim>KNpim_MM_Sp_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wSm_up"), h1-> Fill(mm_kn);
		}

		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wA"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wB"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wC"), h1-> Fill(mm_kn);
		}
		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wD"), h1-> Fill(mm_kn);
		}
	      }

	      if( fKNpipi_N2_flag ){
		if( fKNpim_Sp_flag && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wSp_down"), h1-> Fill(mm_kn);
		}
		if( fKNpim_Sp_flag && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wSp_up"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim<KNpim_MM_Sp_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wSm_down"), h1-> Fill(mm_kn);
		}
		if( fKNpip_Sm_flag && mm_knpim>KNpim_MM_Sp_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wSm_up"), h1-> Fill(mm_kn);
		}

		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wA"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip>KNpip_MM_Sm_MAX ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wB"), h1-> Fill(mm_kn);
		}
		if( mm_knpim>KNpim_MM_Sp_MAX && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wC"), h1-> Fill(mm_kn);
		}
		if( mm_knpim<KNpim_MM_Sp_MIN && mm_knpip<KNpip_MM_Sm_MIN ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wD"), h1-> Fill(mm_kn);
		}
	      }
	      
	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN"), h2-> Fill(mm_knpim, mm_knpip);
		if( !fSm_mass_rc_flag && !fSp_mass_rc_flag ){
		  h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll2_wN"), h2-> Fill(mm_knpim, mm_knpip);
		}

		for( int mm=1300; mm<1600; mm+=5 ){
		  if( 0.001*mm<mm_kn && mm_kn<0.001*mm+0.005 ) h2 = (TH2F*)rtFile-> Get(Form("KNpim_KNpip_MM_woAll_wN_%d", mm)), h2-> Fill(mm_knpim, mm_knpip);
		}
		if( mm_kn<kpMass+pMass ) h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN_down"), h2-> Fill(mm_knpim, mm_knpip);
		else h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN_up"), h2-> Fill(mm_knpim, mm_knpip);

		//		if( 1.6<mm_kn ) h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN_160"), h2-> Fill(mm_knpim, mm_knpip);

		if( mm_kn<1.43 ) h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN_lt143"), h2-> Fill(mm_knpim, mm_knpip);
		if( mm_kn<1.44 ) h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN_lt144"), h2-> Fill(mm_knpim, mm_knpip);

		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN"), h1->Fill(mm_kn);
		if( !fSm_mass_rc_flag && !fSp_mass_rc_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll2_wN"), h1->Fill(mm_kn);
		}
		if( !fSm_mass_rc_flag2 && !fSp_mass_rc_flag2 ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_woS2_wN"), h1->Fill(mm_kn);
		}
		
		h2 = (TH2F*)rtFile-> Get("KNpim_mom"), h2-> Fill(kn_pim_mom.CosTheta(), kn_pim_mom.Mag());
		h2 = (TH2F*)rtFile-> Get("KNpip_mom"), h2-> Fill(kn_pip_mom.CosTheta(), kn_pip_mom.Mag());

		h2 = (TH2F*)rtFile-> Get("KNpim_mom_CM"), h2-> Fill(knpim_mom_CM.CosTheta(), knpim_mom_CM.Mag());
		h2 = (TH2F*)rtFile-> Get("KNpip_mom_CM"), h2-> Fill(knpip_mom_CM.CosTheta(), knpip_mom_CM.Mag());

		h2 = (TH2F*)rtFile-> Get("KN_mom"), h2-> Fill(kn_mom.CosTheta(), kn_mom.Mag());
		h2 = (TH2F*)rtFile-> Get("KN_mom"), h2-> Fill(kn_mom.CosTheta(), kn_mom.Mag());

		h2 = (TH2F*)rtFile-> Get("KN_mom_CM"), h2-> Fill(kn_mom_CM.CosTheta(), kn_mom_CM.Mag());
		h2 = (TH2F*)rtFile-> Get("KN_mom_CM"), h2-> Fill(kn_mom_CM.CosTheta(), kn_mom_CM.Mag());
		
		h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip_woAll_wN"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
		h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip_woAll_wN"), h2-> Fill(dis_beam_pim, dis_beam_pip);

		if( fKNpim_Sp_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSp"), h1->Fill(mm_kn);
		  if( !fKNpip_Sm_flag ){
		    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSp_woSm"), h1->Fill(mm_kn);
		  }
		}

		if( fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSm"), h1->Fill(mm_kn);
		  if( !fKNpim_Sp_flag ){
		    h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wSm_woSp"), h1->Fill(mm_kn);
		  }
		}
		
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN_wS"), h1->Fill(mm_kn);

		  h1 = (TH1F*)rtFile-> Get("Vtxz_pim_pip_woAll_wN_wS"), h1-> Fill(vtx_pim_beam_mean.Z()-vtx_pip_beam_mean.Z());
		  h2 = (TH2F*)rtFile-> Get("Vtx_bpim_bpip_woAll_wN_wS"), h2-> Fill(dis_beam_pim, dis_beam_pip);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN0_wS"), h1->Fill(mm_kn);
		}
	      }

	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN1_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll_wN2_wS"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( !fKNpipi_K0_flag25 && !fKNpipi_Sp_flag25 && !fKNpipi_Sm_flag25 ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll25"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll25"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25"), h2-> Fill(mm_knpim, mm_knpip);

	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN0"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN0_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN1"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN1_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll25_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN2"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll25_wN2_wS"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( !fKNpipi_K0_flag && !fKNpipi_Sp_flag3 && !fKNpipi_Sm_flag3 ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll3"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll3"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3"), h2-> Fill(mm_knpim, mm_knpip);

	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN0"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN0_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN1"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN1_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll3_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN2"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll3_wN2_wS"), h1->Fill(mm_kn);
		}
	      }
	    }

	    if( !fKNpipi_K0_flag4 && !fKNpipi_Sp_flag4 && !fKNpipi_Sm_flag4 ){
	      h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4"), h1->Fill(mm_kn);
	      h1 = (TH1F*)rtFile-> Get("KNpipi_MM_woAll4"), h1->Fill(mm_knpipi);
	      h2 = (TH2F*)rtFile-> Get("KNpipi_MM_mom_woAll4"), h2->Fill(mm_knpipi, knpipi_lmom.Vect().Mag());
	      h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4"), h2-> Fill(mm_knpim, mm_knpip);

	      if( fKNpipi_N_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN"), h1->Fill(mm_kn);

		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N0_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN0"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN0"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN0_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N1_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN1"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN1"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN1_wS"), h1->Fill(mm_kn);
		}
	      }
	      if( fKNpipi_N2_flag ){
		h2 = (TH2F*)rtFile-> Get("KNpim_KNpip_MM_woAll4_wN2"), h2-> Fill(mm_knpim, mm_knpip);
		h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN2"), h1->Fill(mm_kn);
		if( fKNpim_Sp_flag || fKNpip_Sm_flag ){
		  h1 = (TH1F*)rtFile-> Get("KN_MM_pipi_woAll4_wN2_wS"), h1->Fill(mm_kn);
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
	    h2 = (TH2F*)rtFile-> Get("KN_MM_Npipi_IM"), h2-> Fill(mm_kn, im_npipi);
	    h1 = (TH1F*)rtFile-> Get("KP_Npipi_MM"), h1-> Fill(mm_kp_npipi);
	    h2 = (TH2F*)rtFile-> Get("Npipi_IM_KP_Npipi_MM"), h2-> Fill(im_npipi, mm_kp_npipi);
	    h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM"), h2-> Fill(im_npipi, mm_knpipi);

	    if( fKNpipi_K0_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npipi_IM_wK0"), h1-> Fill(im_npipi);
	      h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM_wK0"), h2-> Fill(im_npipi, mm_knpipi);
	    }

	    if( fKNpipi_Sp_flag || fKNpipi_Sm_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npipi_IM_wSf"), h1-> Fill(im_npipi);
	      h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM_wSf"), h2-> Fill(im_npipi, mm_knpipi);
	    }

	    if( fKNpipi_Sm_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npipi_IM_wSm"), h1-> Fill(im_npipi);
	      h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM_wSm"), h2-> Fill(im_npipi, mm_knpipi);
	      if( !fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wSm_woSp"), h1-> Fill(im_npipi);
		h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM_wSm_woSp"), h2-> Fill(im_npipi, mm_knpipi);
	      }
	    }

	    if( fKNpipi_Sp_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npipi_IM_wSp"), h1-> Fill(im_npipi);
	      h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM_wSp"), h2-> Fill(im_npipi, mm_knpipi);
	      if( !fKNpipi_Sm_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wSp_woSm"), h1-> Fill(im_npipi);
		h2 = (TH2F*)rtFile-> Get("Npipi_IM_KNpipi_MM_wSp_woSm"), h2-> Fill(im_npipi, mm_knpipi);
	      }
	    }

	    if( fKNpipi_N_flag ){
	      h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN"), h1-> Fill(im_npipi);
	      h2 = (TH2F*)rtFile-> Get("KN_MM_Npipi_IM_wN"), h2-> Fill(mm_kn, im_npipi);
	      if( fKNpipi_K0_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0"), h1-> Fill(im_npipi);
		h2 = (TH2F*)rtFile-> Get("KN_MM_Npipi_IM_wN_K0"), h2-> Fill(mm_kn, im_npipi);
	      }
	      if( fKNpipi_K0_SB0_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0_SB0"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_K0_SB1_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0_SB1"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_K0_SB2_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0_SB2"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_K0_SB3_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0_SB3"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_K0_SB4_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0_SB4"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_K0_SB5_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wK0_SB5"), h1-> Fill(im_npipi);
	      }

	      if( fKNpipi_Sm_flag || fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wSf"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_Sm_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wSm"), h1-> Fill(im_npipi);
		h2 = (TH2F*)rtFile-> Get("KN_MM_Npipi_IM_wN_Sm"), h2-> Fill(mm_kn, im_npipi);
		if( !fKNpipi_Sp_flag ) h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wSm_woSp"), h1-> Fill(im_npipi);
	      }
	      if( fKNpipi_Sp_flag ){
		h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wSp"), h1-> Fill(im_npipi);
		h2 = (TH2F*)rtFile-> Get("KN_MM_Npipi_IM_wN_Sp"), h2-> Fill(mm_kn, im_npipi);
		if( !fKNpipi_Sm_flag ) h1 = (TH1F*)rtFile-> Get("Npipi_IM_wN_wSp_woSm"), h1-> Fill(im_npipi);
	      }
	    }
	  }
	}
      }
    }
  }

  //****************************************//
  //** K- d -> n p(BPD) pi- tagged event ***//
  //****************************************//
  if( fFPID=F_Neutron && trigNC_flag && fBeamPID==Beam_Kaon && fTrackPim.size()==1 ){
    if( simReader){
      simReader-> fillppim_MC(rtFile);
    }
  }

  fill_pppim(header);
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

  fFCflag=false;
  fFChit       = 0;
  fFDC1track   = 0;
  fFC_start    = DEFVECT;
  fFC_FDC1pos  = DEFVECT;
  fFC_hitpos   = DEFVECT;
  fFC_Angle    = DEFAULTD;
  fFC_Mom_USWK = DEFAULTD;
  fFC_Mom_TOF  = DEFAULTD;



  if( simReader ) simReader-> clear();

  return flag;
}

