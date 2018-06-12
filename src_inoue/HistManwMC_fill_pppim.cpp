#include "HistManwMC.h"
#include "MyParam.h"

void HistManwMC::fill_pppim(EventHeader *header)
{
  TH1F *h1;
  TH2F *h2;

  if( cdstrackMan->nGoodTrack()==3 ){
    if( fTrackP.size()==2 && fTrackPim.size()==1 ){
      TVector3 pim_mom0, pim_mom1, p_mom0, p_mom1;
      TLorentzVector pim_lmom0, pim_lmom1, p_lmom0, p_lmom1;
      TVector3 vtx_pim0, vtx_pim1, vtx_p0, vtx_p1;
      if( !TrackTools::Calc2HelixVertex(fTrackPim[0], fTrackP[0], vtx_pim0, vtx_p0) ) return;
      if( !fTrackPim[0]->GetMomentum(vtx_pim0, pim_mom0, true, true) ) return;
      if( !fTrackP[0]->GetMomentum(vtx_p0, p_mom0, true, true) ) return;
      pim_lmom0.SetVectM(pim_mom0, piMass);
      p_lmom0.SetVectM(p_mom0, pMass);
      double ppim_im0 = (pim_lmom0+p_lmom0).M();

      if( !TrackTools::Calc2HelixVertex(fTrackPim[0], fTrackP[1], vtx_pim1, vtx_p1) ) return;
      if( !fTrackPim[0]->GetMomentum(vtx_pim1, pim_mom1, true, true) ) return;
      if( !fTrackP[1]->GetMomentum(vtx_p1, p_mom1, true, true) ) return;
      pim_lmom1.SetVectM(pim_mom1, piMass);
      p_lmom1.SetVectM(p_mom1, pMass);
      double ppim_im1 = (pim_lmom1+p_lmom1).M();
      TVector3 vtx0 = 0.5*(vtx_pim0+vtx_p0);
      TVector3 vtx1 = 0.5*(vtx_pim1+vtx_p1);
      double dltmp=0;
      double dca;
      TVector3 vtx_b0, vtx_cds0, vtx_b1, vtx_cds1;
      MathTools::LineToLine(vtx0, (p_mom0+pim_mom0).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtx_cds0, vtx_b0);
      double dis0 = (vtx_b0-vtx0).Mag();
      MathTools::LineToLine(vtx1, (p_mom1+pim_mom1).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtx_cds1, vtx_b1);
      double dis1 = (vtx_b1-vtx1).Mag();
      if( GeomTools::GetID(vtx_b0)==CID_Fiducial && GeomTools::GetID(vtx_b1)==CID_Fiducial ){
	if( fBeamPID==Beam_Kaon ){
	  //	  h2 = (TH2F*)rtFile-> Get("CDS_IM_ppim_ppim"), h2-> Fill(ppim_im0, ppim_im1);
	  h1 = (TH1F*)rtFile-> Get("CDS_IM_ppim_2p"), h1-> Fill(ppim_im0), h1-> Fill(ppim_im1);
	}
      }

      bool L_flag0 = false, L_flag1 = false;
      if( L_MIN<ppim_im0 && ppim_im0<L_MAX ) L_flag0 = true;
      if( L_MIN<ppim_im1 && ppim_im1<L_MAX ) L_flag1 = true;

      if( !L_flag0 && !L_flag1 ) return;

      TLorentzVector L_lmom, p_lmom;
      TVector3 p_mom, vtx_cds;
      TVector3 vtx_p, vtx_beam0, vtx_beam1, vtx_L;
      if( L_flag0 && !L_flag1 ){
	L_lmom = pim_lmom0+p_lmom0;
	vtx_L = vtx0;
	fTrackP[1]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam0, vtx_p);
	if( !fTrackP[1]-> GetMomentum(vtx_p, p_mom, true, true) ) return;
	MathTools::LineToLine(vtx0, (p_mom0+pim_mom0).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtx_cds, vtx_beam1);
      }
      else if( L_flag1 && !L_flag0 ){
	L_lmom = pim_lmom1+p_lmom1;
	vtx_L = vtx1;
	fTrackP[0]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam0, vtx_p);
	if( !fTrackP[0]-> GetMomentum(vtx_p, p_mom, true, true) ) return;
	MathTools::LineToLine(vtx1, (p_mom1+pim_mom1).Unit(), fT0pos, fTrackBPC->GetMomDir(), dltmp, dca, vtx_cds, vtx_beam1);
      }
      else if( L_flag1 && L_flag0 ){
	if( dis0>dis1 ){
	  L_lmom = p_lmom0+pim_lmom0;
	  vtx_L = vtx0;
	  vtx_beam1 = vtx_b0;
	  fTrackP[1]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam0, vtx_p);
	  if( !fTrackP[1]-> GetMomentum(vtx_p, p_mom, true, true) ) return;
	}
	else{
	  L_lmom = p_lmom1+pim_lmom1;
	  vtx_L = vtx1;
	  vtx_beam1 = vtx_b1;
	  fTrackP[0]->GetVertex(fT0pos, fTrackBPC->GetMomDir(), vtx_beam0, vtx_p);
	  if( !fTrackP[1]-> GetMomentum(vtx_p, p_mom, true, true) ) return;
	}
      }
      p_lmom.SetVectM(p_mom, pMass);

      TLorentzVector Lp_lmom = L_lmom+p_lmom;

      double beam_out, beam_tof;
      ELossTools::CalcElossBeamTGeo(fT0pos, vtx_beam0, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
      TVector3 beam_mom = fTrackBPC-> GetMomDir();
      beam_mom.SetMag(beam_out);
      TLorentzVector beam_lmom;
      beam_lmom.SetVectM(beam_mom, parMass[fBeamPID]);
      TLorentzVector kn_Lp_lmom = beam_lmom+D_LMOM-L_lmom-p_lmom;

      if( GeomTools::GetID(vtx_beam0)==CID_Fiducial && GeomTools::GetID(vtx_beam1)==CID_Fiducial ){
	if( fBeamPID==Beam_Kaon ){
	  h1 = (TH1F*)rtFile-> Get("CDS_IM_Lp"), h1-> Fill(Lp_lmom.M());
	  h2 = (TH2F*)rtFile-> Get("CDS_IM_Lp_KCDSLp_MM"), h2-> Fill(Lp_lmom.M(), kn_Lp_lmom.M());
	  h1 = (TH1F*)rtFile-> Get("Missing_mom_Lp"), h1-> Fill(kn_Lp_lmom.Vect().Mag());
	  h1 = (TH1F*)rtFile-> Get("Missing_ang_Lp"), h1-> Fill(kn_Lp_lmom.Vect().CosTheta());
	  h2 = (TH2F*)rtFile-> Get("Missing_ang_mom_Lp"), h2-> Fill(kn_Lp_lmom.Vect().CosTheta(), kn_Lp_lmom.Vect().Mag());
	}
      }
    }
  }
}
