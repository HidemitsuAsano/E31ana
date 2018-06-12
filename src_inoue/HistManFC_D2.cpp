#include "HistManBeamAna.h"

//static const double pi_v = 0.992046*Const*m; // 1.1GeV/c   
static const double pi_v = 0.9904*Const*m;   // 1.0GeV/c
//static const double pi_v = 0.988188*Const*m; // 0.9GeV/c
static const double USWK_Z = 250; //cm
static const double NC_beta_thre = 0.9;
//static const TLorentzVector TargetLmom(0., 0., 0., dMass);
static const TLorentzVector TargetLmom(0., 0., 0., dMass);

void HistManBeamAna::initFC()
{
  fFile-> cd();

  new TH1F("KCDH1N_Reduction", "K#timesCDH1#timesN Event Reduction", 20, -0.5, 19.5);
  new TH1F("Neutral_Reduction", "Neuterl Event Reduction", 20, -0.5, 19.5);

  new TH1F("NC_overbeta",  "NC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("PC_overbeta",  "CVC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("CVC_overbeta", "PC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("NC_overbeta_wtar",  "NC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("PC_overbeta_wtar",  "CVC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("CVC_overbeta_wtar", "PC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("hitpat_NC_n", "hitpat NC neutral", 112, 0.5, 112.5);
  new TH1F("hitpat_NC_n_wtar", "hitpat NC neutral", 112, 0.5, 112.5);

  new TH1F("T0NC_offset", "NC offset", 500, -10, 10);
  new TH1F("T0NC_offset_8MeVee", "NC offset", 500, -10, 10);
  new TH1F("T0NC_offset_8MeVee_wtar", "NC offset", 500, -10, 10);
  new TH1F("NC_overbeta_8MeVee",  "NC 1/#beta", 5000, -0.5, 4.5);
  new TH1F("NC_overbeta_8MeVee_wtar",  "NC 1/#beta", 5000, -0.5, 4.5);

  new TH2F("NC_overbeta_dE",       "NC 1/#beta vs dE", 500, -0.5, 4.5, 1000, 0.0, 100.);
  new TH2F("NC_overbeta_dE_wtar",  "NC 1/#beta vs dE", 500, -0.5, 4.5, 1000, 0.0, 100.);

  new TH1F("KN_MM",     "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_km",  "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pim", "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pip", "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_p",   "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_d",   "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi","d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_L",   "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_K0",  "d(K^{-}, n) MM", 3000, 0.0, 3.0);

  new TH1F("KN_MM_wtar",     "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_km_wtar",  "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pim_wtar", "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pip_wtar", "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_p_wtar",   "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_d_wtar",   "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wtar","d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_L_wtar",   "d(K^{-}, n) MM", 3000, 0.0, 3.0);
  new TH1F("KN_MM_K0_wtar",  "d(K^{-}, n) MM", 3000, 0.0, 3.0);

  new TH1F("Npim_IM",  "n #pi^{-} IM",           2000, 1.0, 3.0);
  new TH1F("Npip_IM",  "n #pi^{+} IM",           2000, 1.0, 3.0);
  new TH1F("Nkm_IM",   "n K^{-} IM",             2000, 1.0, 3.0);
  new TH1F("Npipi_IM", "n #pi^{+} #pi^{-} IM",   2000, 1.0, 3.0);
  new TH1F("NK0_IM",   "n K^{0} IM",             2000, 1.0, 3.0);

  new TH1F("Npim_IM_wtar",  "n #pi^{-} IM",           2000, 1.0, 3.0);
  new TH1F("Npip_IM_wtar",  "n #pi^{+} IM",           2000, 1.0, 3.0);
  new TH1F("Nkm_IM_wtar",   "n K^{-} IM",             2000, 1.0, 3.0);
  new TH1F("Npipi_IM_wtar", "n #pi^{+} #pi^{-} IM",   2000, 1.0, 3.0);
  new TH1F("NK0_IM_wtar",   "n K^{0} IM",             2000, 1.0, 3.0);

  new TH1F("KNpim_MM",  "d(K^{-}, n #pi^{-}) MM",         2000, 0.5, 2.5);
  new TH1F("KNpip_MM",  "d(K^{-}, n #pi^{+}) MM",         2000, 0.5, 2.5);
  new TH1F("KNpipi_MM", "d(K^{-}, n #pi^{+} #pi^{-}) MM", 2000, 0.0, 2.0);
  new TH1F("KNkm_MM",   "d(K^{-}, n #pi^{+} #pi^{-}) MM", 2000, 0.0, 2.0);
  new TH1F("KNK0_MM",   "d(K^{-}, n K^{0}) MM",           2000, 0.0, 2.0);

  new TH1F("KNpim_MM_wtar",  "d(K^{-}, n #pi^{-}) MM",         2000, 0.5, 2.5);
  new TH1F("KNpip_MM_wtar",  "d(K^{-}, n #pi^{+}) MM",         2000, 0.5, 2.5);
  new TH1F("KNkm_MM_wtar",   "d(K^{-}, n #pi^{+} #pi^{-}) MM", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wtar", "d(K^{-}, n #pi^{+} #pi^{-}) MM", 2000, 0.0, 2.0);
  new TH1F("KNK0_MM_wtar",   "d(K^{-}, n K^{0}) MM",           2000, 0.0, 2.0);

  new TH2F("KNpip_KNpim_MM",      "d(K^{-}, n #pi^{+}) vs d(K^{-}, n #pi^{-})", 200, 1.0, 2.0,200, 1.0, 2.0);
  new TH2F("KNpip_KNpim_MM_woK0", "d(K^{-}, n #pi^{+}) vs d(K^{-}, n #pi^{-})", 200, 1.0, 2.0,200, 1.0, 2.0);
  new TH2F("KNpip_KNpim_MM_wtar",      "d(K^{-}, n #pi^{+}) vs d(K^{-}, n #pi^{-})", 200, 1.0, 2.0,200, 1.0, 2.0);
  new TH2F("KNpip_KNpim_MM_woK0_wtar", "d(K^{-}, n #pi^{+}) vs d(K^{-}, n #pi^{-})", 200, 1.0, 2.0,200, 1.0, 2.0);

  new TH1F("CDS_IM_K0_n", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0); 
  new TH1F("CDS_IM_K0_n_wtar", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0); 
  new TH1F("CDS_IM_pipi_n", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0); 
  new TH1F("CDS_IM_pipi_n_wtar", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0); 

  fFile->cd();
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("T0CVC%d_offset", seg), Form("CVC seg%d offset by T0-CVC", seg), 20000, -100, 100);
    new TH1F(Form("T0CVC%d_offset_wtar", seg), Form("CVC seg%d offset by T0-CVC", seg), 20000, -100, 100);
  }
  for( int seg=1; seg<=27; seg++ ){
    new TH1F(Form("T0PC%d_offset", seg), Form("PC seg%d offset by T0-PC", seg), 20000, -100, 100);
    new TH1F(Form("T0PC%d_offset_wtar", seg), Form("PC seg%d offset by T0-PC", seg), 20000, -100, 100);
  }
  for( int seg=1; seg<=112; seg++ ){
    new TH1F(Form("T0NC%d_offset", seg), Form("NC seg%d offset by T0-NC", seg), 20000, -100, 100);
    new TH1F(Form("T0NC%d_offset_wtar", seg), Form("NC seg%d offset by T0-NC", seg), 20000, -100, 100);
    new TH1F(Form("T0NC%d_offset_8MeVee", seg), Form("NC seg%d offset by T0-NC", seg), 20000, -100, 100);
    new TH1F(Form("T0NC%d_offset_8MeVee_wtar", seg), Form("NC seg%d offset by T0-NC", seg), 20000, -100, 100);

    new TH1F(Form("NC_overbeta_%d", seg), Form("NC seg%d 1/#beta", seg), 5000, -0.5, 4.5);
    new TH1F(Form("NC_overbeta_%d_wtar", seg), Form("NC seg%d 1/#beta", seg), 5000, -0.5, 4.5);

    new TH1F(Form("NC_overbeta_%d_8MeVee", seg), Form("NC seg%d 1/#beta", seg), 5000, -0.5, 4.5);
    new TH1F(Form("NC_overbeta_%d_8MeVee_wtar", seg), Form("NC seg%d 1/#beta", seg), 5000, -0.5, 4.5);

    new TH2F(Form("euctm_NC_%d", seg), Form("NC seg%d up dE vs time", seg),   1000, 0, 50, 500, -10, 10);
    new TH2F(Form("edctm_NC_%d", seg), Form("NC seg%d down dE vs time", seg), 1000, 0, 50, 500, -10, 10);

    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("T0%dNC%d_offset", seg2, seg), Form("NC seg%d - T0 seg%d offset", seg, seg2), 20000, -100, 100);
      new TH1F(Form("T0%dNC%d_offset_wtar", seg2, seg), Form("NC seg%d - T0 seg%d offset", seg, seg2), 20000, -100, 100);
      new TH1F(Form("T0%dNC%d_offset_8MeVee", seg2, seg), Form("NC seg%d - T0 seg%d offset", seg, seg2), 20000, -100, 100);
      new TH1F(Form("T0%dNC%d_offset_8MeVee_wtar", seg2, seg), Form("NC seg%d - T0 seg%d offset", seg, seg2), 20000, -100, 100);

      new TH1F(Form("NC_overbeta_%d_T0%d", seg, seg2),      Form("T0 seg%d - NC seg%d 1/#beta", seg2, seg), 5000, -0.5, 4.5);
      new TH1F(Form("NC_overbeta_%d_T0%d_wtar", seg, seg2), Form("T0 seg%d - NC seg%d 1/#beta", seg2, seg), 5000, -0.5, 4.5);

      new TH1F(Form("NC_overbeta_%d_8MeVee_T0%d", seg, seg2),      Form("T0 seg%d - NC seg%d 1/#beta", seg2, seg), 5000, -0.5, 4.5);
      new TH1F(Form("NC_overbeta_%d_8MeVee_T0%d_wtar", seg, seg2), Form("T0 seg%d - NC seg%d 1/#beta", seg2, seg), 5000, -0.5, 4.5);
    }
  }
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("T0%dNC_offset", seg), Form("NC - T0 seg%d offset", seg), 20000, -100, 100);
    new TH1F(Form("T0%dNC_offset_wtar", seg), Form("NC - T0 seg%d offset", seg), 20000, -100, 100);
    new TH1F(Form("T0%dNC_offset_8MeVee", seg), Form("NC - T0 seg%d offset", seg), 20000, -100, 100);
    new TH1F(Form("T0%dNC_offset_8MeVee_wtar", seg), Form("NC - T0 seg%d offset", seg), 20000, -100, 100);

    new TH2F(Form("euctm_T0_%d_NC", seg), Form("T0 seg%d up dE vs offset byNC", seg),   500, 0., 10, 500, -10, 10);
    new TH2F(Form("edctm_T0_%d_NC", seg), Form("T0 seg%d down dE vs offset byNC", seg), 500, 0., 10, 500, -10, 10);
  }
}

void HistManBeamAna::fillFC(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan)
{
  TH1 *h1;
  TH2 *h2;
  TTree *tree;

  if( fT0_hit.size()==1 && bltrackMan->ntrackFDC1()==1 && bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ){
    LocalTrack *trackBLC2 = bltrackMan->trackBLC2(0);
    LocalTrack *trackBPC  = bltrackMan->trackBPC(0);
    LocalTrack *trackFDC  = bltrackMan->trackFDC1(0);
    TVector3 gBPCpos, gFDCpos, gBPCdir, gFDCdir;
    conf-> GetBLDCWireMapManager()-> GetGParam(CID_FDC1, gFDCpos, gFDCdir);
    conf-> GetBLDCWireMapManager()-> GetGParam(CID_BPC, gBPCpos, gBPCdir);
    TVector3 pos;
    conf->GetGeomMapManager()-> GetGPos(CID_T0, fT0_hit[0]->seg(), pos);
    TVector3 T0pos = trackBLC2-> GetPosatZ(pos.Z()-0.5);
    TVector3 FDCpos = trackFDC-> GetPosatZ(gFDCpos.Z());
    TVector3 BPCpos = trackBPC-> GetPosatZ(gBPCpos.Z());
    TVector3 dir = (FDCpos-BPCpos).Unit();
    TVector3 USWKpos = FDCpos+(250-FDCpos.Z())*dir;

    if( fCVC_hit.size()+fPC_hit.size()==1 ){
      HodoscopeLikeHit *charge_hit =0;
      TVector3 charge_pos;
      double tof;
      if( fCVC_hit.size()==1 ){
        tof = fCVC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
        conf-> GetGeomMapManager()-> GetGPos(CID_CVC, fCVC_hit[0]->seg(), charge_pos);
      }
      else{
        tof = fPC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
        conf-> GetGeomMapManager()-> GetGPos(CID_PC, fPC_hit[0]->seg(), charge_pos);
      }
      double bend_angle = (USWKpos-BPCpos).Angle(charge_pos-USWKpos);

      double leff = 100;
      double r = leff/sqrt(2*(1-cos(bend_angle)));
      double larc = r*bend_angle;
      double line = 2*leff/sqrt(2*(1+cos(bend_angle)));
      double diff = line-larc;
      double fl = (T0pos-BPCpos).Mag()+(BPCpos-USWKpos).Mag()+(USWKpos-charge_pos).Mag()-diff;
      fl -= 1.5;
      double offset = tof - fl/pi_v;

      if( fCVC_hit.size()==1 ){
#if SLEWING_DATA
	h1 = getTH(Form("T0CVC%d_offset", fCVC_hit[0]->seg())), h1->Fill(offset);
	fSlewingT0CVC-> Set(fT0_hit[0], fCVC_hit[0]);
	fSlewingT0CVC-> SetOffset(2, fl/pi_v);
	tree = (TTree*)fFile-> Get("TreeT0CVC");
	tree-> Fill();
#endif
      }
      else{
#if SLEWING_DATA
	h1 = getTH(Form("T0PC%d_offset", fPC_hit[0]->seg())), h1->Fill(offset);
	fSlewingT0PC-> Set(fT0_hit[0], fPC_hit[0]);
	fSlewingT0PC-> SetOffset(2, fl/pi_v);
	tree = (TTree*)fFile-> Get("TreeT0PC");
	tree-> Fill();
#endif
      }
      if( fCVC_hit.size()==1 ){
	for( int lay=0; lay<8; lay++ ){
	  if( fNC_hit[lay].size()==1 ){
	    HodoscopeLikeHit *nc_hit = fNC_hit[lay][0];
	    tof = nc_hit->ctmean()-fT0_hit[0]->ctmean();
	    TVector3 nc_pos;
	    conf-> GetGeomMapManager()-> GetGPos(CID_NC, nc_hit->seg(), nc_pos);
	    bend_angle = (USWKpos-BPCpos).Angle(nc_pos-USWKpos);
	    r = leff/sqrt(2*(1-cos(bend_angle)));
	    larc = r*bend_angle;
	    line = 2*leff/sqrt(2*(1+cos(bend_angle)));
	    diff = line-larc;
	    fl = (T0pos-BPCpos).Mag()+(BPCpos-USWKpos).Mag()+(USWKpos-nc_pos).Mag()-diff;
	    fl -= 2.5;
	    double offse = tof - fl/pi_v;
	    h1 = getTH(Form("T0NC%d_offset", nc_hit->seg())), h1-> Fill(offset);
	    h1 = getTH(Form("T0%dNC%d_offset", fT0_hit[0]->seg(), nc_hit->seg())), h1-> Fill(offset);
	    h1 = getTH(Form("T0%dNC_offset", fT0_hit[0]->seg())), h1-> Fill(offset);
	    h2 = getTH2(Form("euctm_NC_%d", nc_hit->seg())), h2-> Fill(nc_hit->eu(), offset);
	    h2 = getTH2(Form("edctm_NC_%d", nc_hit->seg())), h2-> Fill(nc_hit->ed(), offset);
	    
	    h2 = getTH2(Form("euctm_T0_%d_NC", fT0_hit[0]->seg())), h2->Fill(fT0_hit[0]->eu(), -offset);
	    h2 = getTH2(Form("edctm_T0_%d_NC", fT0_hit[0]->seg())), h2->Fill(fT0_hit[0]->ed(), -offset);
	    
	    if( nc_hit->emean()>8 ){
	      h1 = getTH(Form("T0NC%d_offset_8MeVee", nc_hit->seg())), h1-> Fill(offset);
	      h1 = getTH(Form("T0%dNC%d_offset_8MeVee", fT0_hit[0]->seg(), nc_hit->seg())), h1-> Fill(offset);
	      h1 = getTH(Form("T0%dNC_offset_8MeVee", fT0_hit[0]->seg())), h1-> Fill(offset);
	    }
#if SLEWING_DATA
	    fSlewingT0NC-> Set(fT0_hit[0], nc_hit);
	    fSlewingT0NC-> SetOffset(2, fl/pi_v);
	    tree = (TTree*)fFile-> Get("TreeT0NC");
	    tree-> Fill();
#endif
	  }
	}
      }
    }
  }
}

void HistManBeamAna::fillFC(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan)
{
  TH1 *h1;
  TH2 *h2;
  TTree *tree;

  bool vtx_flag = false;
  bool forward_n = false;
  bool forward_g = false;
  if( fBeamPID==Beam_Pion || fBeamPID==Beam_Kaon ){
    if( fT0_hit.size()==1 && bltrackMan->ntrackBPC()==1 && bltrackMan->ntrackBLC2()==1 ){
      if( fBeamMom>0.0001 ){
	LocalTrack *trackBPC = bltrackMan->trackBPC(0);
	LocalTrack *trackBLC2 = bltrackMan->trackBLC2(0);
	int T0seg = fT0_hit[0]-> seg();
	TVector3 gpos, gdir;
	conf-> GetGeomMapManager()-> GetGPos(CID_T0, T0seg, gpos);
	TVector3 T0pos = trackBPC->GetPosatZ(gpos.Z()-0.5);
	double out_mom;
	double beam_tof;

	double tmpdis = 99999.;
	TVector3 vtx;
	for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
	  CDSTrack *cdsTrack = cdstrackMan->GoodTrack(id);
	  TVector3 vtxb, vtxcds;
	  double tmpl;
	  if( TrackTools::CalcLineHelixVertex(trackBPC, cdsTrack, vtxb, vtxcds, tmpl) ){
	    if( tmpl<tmpdis ){
	      vtx_flag = true;
	      tmpdis=tmpl;
	      vtx = vtxb;
	    }
	  }
	}

	double beam_fl = (T0pos-vtx).Mag();
	ELossTools::CalcElossBeamTGeo(T0pos, vtx, fBeamMom, parMass[fBeamPID], out_mom, beam_tof);

	if( fCVC_hit.size()==0 ){

	  bool first_flag = true;
	  for( int lay=0; lay<8; lay++ ){
	    for( int i=0; i<fNC_hit[lay].size(); i++ ){
	      int NCseg = fNC_hit[lay][i]-> seg();
	      double tof = fNC_hit[lay][i]->ctmean()-fT0_hit[0]->ctmean();
	      tof -= beam_tof;
		
	      conf-> GetGeomMapManager()-> GetGPos(CID_NC, NCseg, gpos);
	      TVector3 NCpos = gpos;
	      double fl = (vtx-NCpos).Mag();
	      double beta = fl/(tof*100.*Const);	      
	      double offset = tof - fl/(100*Const);
#if SLWEING_DATA
	      fSlewingT0NC-> Set(fT0_hit[0], fNC_hit[lay][i]);
	      fSlewingT0NC-> SetOffset(2, fl/(100*Const));
	      tree = (TTree*)fFile-> Get("TreeT0NC");
	      tree-> Fill();
#endif
	      h1 = getTH("hitpat_NC_n"), h1->Fill(fNC_hit[lay][i]->seg());
	      h1 = getTH("T0NC_offset"), h1-> Fill(offset);
	      h1 = getTH(Form("T0NC%d_offset", fNC_hit[lay][i]->seg())), h1-> Fill(offset);
	      h1 = getTH(Form("T0%dNC%d_offset", T0seg, fNC_hit[lay][i]->seg())), h1-> Fill(offset);
	      h1 = getTH(Form("T0%dNC_offset", T0seg)), h1-> Fill(offset);

	      h2 = getTH2(Form("euctm_NC_%d", fNC_hit[lay][i]->seg())), h2-> Fill(fNC_hit[lay][i]->eu(), offset);
	      h2 = getTH2(Form("edctm_NC_%d", fNC_hit[lay][i]->seg())), h2-> Fill(fNC_hit[lay][i]->ed(), offset);
	      if( fNC_hit[lay][i]->emean()>8 ){
		h1 = getTH("T0NC_offset_8MeVee"), h1-> Fill(offset);
		h1 = getTH(Form("T0NC%d_offset_8MeVee", fNC_hit[lay][i]->seg())), h1-> Fill(offset);
		h1 = getTH(Form("T0%dNC%d_offset_8MeVee", T0seg, fNC_hit[lay][i]->seg())), h1-> Fill(offset);
		h1 = getTH(Form("T0%dNC_offset_8MeVee", T0seg)), h1-> Fill(offset);
		if( MyTools::IsTarget(vtx) ){
		  h1 = getTH("T0NC_offset_8MeVee_wtar"), h1-> Fill(offset);
		  h1 = getTH(Form("T0NC%d_offset_8MeVee_wtar", fNC_hit[lay][i]->seg())), h1-> Fill(offset);
		  h1 = getTH(Form("T0%dNC%d_offset_8MeVee_wtar", T0seg, fNC_hit[lay][i]->seg())), h1-> Fill(offset);
		  h1 = getTH(Form("T0%dNC_offset_8MeVee_wtar", T0seg)), h1-> Fill(offset);
		}
	      }

	      if( GeomTools::GetID(vtx)==CID_Fiducial ){
		h1 = getTH("hitpat_NC_n_wtar"), h1->Fill(fNC_hit[lay][i]->seg());
		h1 = getTH(Form("T0NC%d_offset_wtar", fNC_hit[lay][i]->seg())), h1-> Fill(offset);
		h1 = getTH(Form("T0%dNC%d_offset_wtar", T0seg, fNC_hit[lay][i]->seg())), h1-> Fill(offset);
		h1 = getTH(Form("T0%dNC_offset_wtar", T0seg)), h1-> Fill(offset);
	      }

	      if( first_flag ){
		//		std::cout<<"Fill NC"<<std::endl;
		h1 = getTH("NC_overbeta"), h1-> Fill(1./beta);
		h1 = getTH(Form("NC_overbeta_%d", fNC_hit[lay][i]->seg())), h1-> Fill(1./beta);
		h1 = getTH(Form("NC_overbeta_%d_T0%d", fNC_hit[lay][i]->seg(), T0seg)), h1-> Fill(1./beta);
		if( fNC_hit[lay][i]->emean()>8 ){
		  h1 = getTH("NC_overbeta_8MeVee"), h1-> Fill(1./beta);
		  h1 = getTH(Form("NC_overbeta_%d_8MeVee", fNC_hit[lay][i]->seg())), h1-> Fill(1./beta);
		  h1 = getTH(Form("NC_overbeta_%d_8MeVee_T0%d", fNC_hit[lay][i]->seg(), T0seg)), h1-> Fill(1./beta);
		}
		h2 = getTH2("NC_overbeta_dE"), h2-> Fill(1./beta, fNC_hit[lay][i]->emean());
		if( GeomTools::GetID(vtx)==CID_Fiducial ){
		  h1 = getTH("NC_overbeta_wtar"), h1-> Fill(1./beta);
		  h1 = getTH(Form("NC_overbeta_%d_wtar", fNC_hit[lay][i]->seg())), h1-> Fill(1./beta);
		  h1 = getTH(Form("NC_overbeta_%d_T0%d_wtar", fNC_hit[lay][i]->seg(), T0seg)), h1-> Fill(1./beta);
		  if( fNC_hit[lay][i]->emean()>8 ){
		    h1 = getTH("NC_overbeta_8MeVee_wtar"), h1-> Fill(1./beta);
		    h1 = getTH(Form("NC_overbeta_%d_8MeVee_wtar", fNC_hit[lay][i]->seg())), h1-> Fill(1./beta);
		    h1 = getTH(Form("NC_overbeta_%d_8MeVee_T0%d_wtar", fNC_hit[lay][i]->seg(), T0seg)), h1-> Fill(1./beta);
		  }
		  h2 = getTH2("NC_overbeta_dE_wtar"), h2-> Fill(1./beta, fNC_hit[lay][i]->emean());
		}
		
		TLorentzVector fn_lmom;
		//		if( beta<NC_beta_thre ){
		if( beta<NC_beta_thre && fNC_hit[lay][i]->emean()>8 ){
		  forward_n = true;
		  
		  double mom = nMass*beta/sqrt(1-beta*beta);
		  TVector3 n_mom = (NCpos-vtx);
		  n_mom.SetMag(mom);
		  fn_lmom.SetVectM(n_mom, nMass);
		  
		  //		  if( fBeamPID==Beam_Kaon && header-> trigmode2(Mode_KCDH1N) ){
		  if( fBeamPID==Beam_Kaon && header-> IsTrig(Trig_Neutral) ){
		    TLorentzVector beam_lmom;
		    TVector3 beam_mom = trackBPC-> GetMomDir();
		    beam_mom.SetMag(out_mom);
		    beam_lmom.SetVectM(beam_mom, kpMass);
		    
		    double mm = (beam_lmom+TargetLmom-fn_lmom).M();
		    h1 = getTH("KN_MM"), h1-> Fill(mm);
		    if( fCDSpim.size()>0 ) h1 = getTH("KN_MM_pim"), h1-> Fill(mm);
		    if( fCDSkm.size()>0  ) h1 = getTH("KN_MM_km"),  h1-> Fill(mm);
		    if( fCDSpip.size()>0 ) h1 = getTH("KN_MM_pip"), h1-> Fill(mm);
		    if( fCDSp.size()>0   ) h1 = getTH("KN_MM_p"),   h1-> Fill(mm);
		    if( fCDSd.size()>0   ) h1 = getTH("KN_MM_d"),  h1-> Fill(mm);
		    if( fLflag  ) h1 = getTH("KN_MM_L"),  h1-> Fill(mm);
		    if( fK0flag ) h1 = getTH("KN_MM_K0"), h1-> Fill(mm);
		    if( fCDSpim.size()>0 && fCDSpip.size()>0 ){
		      h1 = getTH("KN_MM_pipi"), h1-> Fill(mm);
		    
		      if( GeomTools::GetID(vtx)==CID_Fiducial ){
			h1 = getTH("KN_MM_pipi_wtar"), h1-> Fill(mm);
		      }
		    }
		    
		    for( int i=0; i<fCDSpim.size(); i++ ){
		      CDSTrack *cdsTrack = fCDSpim[i];
		      TVector3 vtxb, vtxh;
		      double tmpdis;
		      if( TrackTools::CalcLineHelixVertex(trackBPC, cdsTrack, vtxb, vtxh, tmpdis) ){
			TVector3 p;
			if( cdsTrack-> GetMomentum(vtxh, p, true, true) ){
			  TLorentzVector lmom;
			  lmom.SetVectM(p, piMass);
			  
			  double mm2 = (beam_lmom+TargetLmom-fn_lmom-lmom).M();
			  double im = (fn_lmom+lmom).M();
			  h1 = getTH("KNpim_MM"), h1-> Fill(mm2);
			  h1 = getTH("Npim_IM"), h1-> Fill(im);
			  if( GeomTools::GetID(vtx)==CID_Fiducial ){
			    h1 = getTH("Npim_IM_wtar"), h1-> Fill(im);
			    h1 = getTH("KNpim_MM_wtar"), h1-> Fill(mm2);
			  }
			}
		      }
		    }

		    for( int i=0; i<fCDSpip.size(); i++ ){
		      CDSTrack *cdsTrack = fCDSpip[i];
		      TVector3 vtxb, vtxh;
		      double tmpdis;
		      if( TrackTools::CalcLineHelixVertex(trackBPC, cdsTrack, vtxb, vtxh, tmpdis) ){
			TVector3 p;
			if( cdsTrack-> GetMomentum(vtxh, p, true, true) ){
			  TLorentzVector lmom;
			  lmom.SetVectM(p, piMass);
			  
			  double mm2 = (beam_lmom+TargetLmom-fn_lmom-lmom).M();
			  double im = (fn_lmom+lmom).M();
			  h1 = getTH("KNpip_MM"), h1-> Fill(mm2);
			  h1 = getTH("Npip_IM"), h1-> Fill(im);
			  if( GeomTools::GetID(vtx)==CID_Fiducial ){
			    h1 = getTH("KNpip_MM_wtar"), h1-> Fill(mm2);
			    h1 = getTH("Npip_IM_wtar"), h1-> Fill(im);
			  }
			}
		      }
		    }
		    
		    for( int i=0; i<fCDSkm.size(); i++ ){
		      CDSTrack *cdsTrack = fCDSkm[i];
		      TVector3 vtxb, vtxh;
		      double tmpdis;
		      if( TrackTools::CalcLineHelixVertex(trackBPC, cdsTrack, vtxb, vtxh, tmpdis) ){
			TVector3 p;
			if( cdsTrack-> GetMomentum(vtxh, p, true, true) ){
			  TLorentzVector lmom;
			  lmom.SetVectM(p, kpMass);

			  double mm2 = (beam_lmom+TargetLmom-fn_lmom-lmom).M();
			  double im = (fn_lmom+lmom).M();
			  h1 = getTH("KNkm_MM"), h1-> Fill(mm2);
			  h1 = getTH("Nkm_IM"), h1-> Fill(im);
			  if( GeomTools::GetID(vtx)==CID_Fiducial ){
			    h1 = getTH("KNkm_MM_wtar"), h1-> Fill(mm2);
			    h1 = getTH("Nkm_IM_wtar"), h1-> Fill(im);
			  }
			}
		      }
		    }

		    for( int i=0; i<fCDSpim.size(); i++ ){
		      for( int j=0; j<fCDSpip.size(); j++ ){
			CDSTrack *cdsTrack1 = fCDSpim[i];
			CDSTrack *cdsTrack2 = fCDSpip[j];
			TVector3 vtxh1, vtxh2;
			if( TrackTools::Calc2HelixVertex(cdsTrack1, cdsTrack2, vtxh1, vtxh2) ){
			  TVector3 pim_mom, pip_mom;
			  if( cdsTrack1->GetMomentum(vtxh1, pim_mom, true, true) ){
			    if( cdsTrack2->GetMomentum(vtxh2, pip_mom, true, true) ){
			      TLorentzVector pim_lmom, pip_lmom;
			      pim_lmom.SetVectM(pim_mom, piMass);
			      pip_lmom.SetVectM(pip_mom, piMass);
			      TVector3 vtx_mean = 0.5*(vtxh1+vtxh2);
			      double dist, dltmp=0;
			      TVector3 vtxb, vtxp;
			      MathTools::LineToLine(vtx_mean, (pim_mom+pip_mom).Unit(), T0pos, trackBPC->GetMomDir(), dltmp, dist, vtxb, vtxp);
			      
			      double im = (pim_lmom+pip_lmom).M();
			      double im2 = (pim_lmom+pip_lmom+fn_lmom).M();
			      double mm2 = (beam_lmom+TargetLmom-fn_lmom-pim_lmom-pip_lmom).M();
			      h1 = getTH("CDS_IM_pipi_n"), h1-> Fill(im);
			      h1 = getTH("KNpipi_MM"), h1-> Fill(mm2);
			      h1 = getTH("Npipi_IM"), h1-> Fill(im2);
			      
			      h2 = getTH2("KNpip_KNpim_MM"), h2-> Fill((beam_lmom+TargetLmom-fn_lmom-pip_lmom).M(), 
								       (beam_lmom+TargetLmom-fn_lmom-pim_lmom).M());
			      if( K0_MIN<im && im<K0_MAX ){
				h1 = getTH("CDS_IM_K0_n"), h1-> Fill(im);
				h1 = getTH("KNK0_MM"), h1-> Fill(mm2);
				h1 = getTH("NK0_IM"), h1-> Fill(im2);
				h2 = getTH2("KNpip_KNpim_MM_woK0"), h2-> Fill((beam_lmom+TargetLmom-fn_lmom-pip_lmom).M(), 
									      (beam_lmom+TargetLmom-fn_lmom-pim_lmom).M());
			      }
			      if( GeomTools::GetID(vtx)==CID_Fiducial ){
				h1 = getTH("CDS_IM_pipi_n_wtar"), h1-> Fill(im);
				h1 = getTH("KNpipi_MM_wtar"), h1-> Fill(mm2);
				h1 = getTH("Npipi_IM_wtar"), h1-> Fill(im2);
				h2 = getTH2("KNpip_KNpim_MM_wtar"), h2-> Fill((beam_lmom+TargetLmom-fn_lmom-pip_lmom).M(), 
									      (beam_lmom+TargetLmom-fn_lmom-pim_lmom).M());
				if( K0_MIN<im && im<K0_MAX ){
				  h1 = getTH("CDS_IM_K0_n_wtar"), h1-> Fill(im);
				  h1 = getTH("KNK0_MM_wtar"), h1-> Fill(mm2);
				  h1 = getTH("NK0_IM_wtar"), h1-> Fill(im2);
				  h2 = getTH2("KNpip_KNpim_MM_woK0_wtar"), h2-> Fill((beam_lmom+TargetLmom-fn_lmom-pip_lmom).M(), 
										     (beam_lmom+TargetLmom-fn_lmom-pim_lmom).M());
				}
			      }
			    }
			  }
			}
		      }
		    }

		  
		    if( GeomTools::GetID(vtx)==CID_Fiducial ){
		      h1 = getTH("KN_MM_wtar"), h1-> Fill(mm);
		      if( fCDSpim.size()>0 ) h1 = getTH("KN_MM_pim_wtar"), h1-> Fill(mm);
		      if( fCDSkm.size()>0  ) h1 = getTH("KN_MM_km_wtar"),  h1-> Fill(mm);
		      if( fCDSpip.size()>0 ) h1 = getTH("KN_MM_pip_wtar"), h1-> Fill(mm);
		      if( fCDSp.size()>0   ) h1 = getTH("KN_MM_p_wtar"),   h1-> Fill(mm);
		      if( fCDSd.size()>0   ) h1 = getTH("KN_MM_d_wtar"),  h1-> Fill(mm);
		      if( fLflag  ) h1 = getTH("KN_MM_L_wtar"),  h1-> Fill(mm);
		      if( fK0flag ) h1 = getTH("KN_MM_K0_wtar"), h1-> Fill(mm);
		    }
		  }
		}
		else forward_g=true;

		first_flag = false;
	      }
	    }
	    if( fNC_hit[lay].size()>0 ) break;
	  }
	}
	
	if( bltrackMan-> ntrackFDC1()==1 && fCVC_hit.size()==1 && fPC_hit.size()==0 ){
	  LocalTrack *trackFDC = bltrackMan->trackFDC1(0);
	  conf-> GetBLDCWireMapManager()-> GetGParam(CID_FDC1, gpos, gdir);
	  TVector3 FDCpos = trackFDC->GetPosatZ(gpos.Z());
	  TVector3 dir = (FDCpos-vtx).Unit(); 
	  TVector3 USWKpos = FDCpos+(250-FDCpos.Z())*dir;
	  int CVCseg = fCVC_hit[0]-> seg();
	  
	  double tof = fCVC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
	  tof -= beam_tof;
	  
	  TVector3 CVCpos;
	  conf-> GetGeomMapManager()-> GetGPos(CID_CVC, CVCseg, CVCpos);
	  double angle = (USWKpos-vtx).Angle(CVCpos-USWKpos);
	  double leff = 100.;
	  double r = leff/sqrt(2*(1-cos(angle)));
	  double larc = r*angle;
	  double line = 2*leff/sqrt(2*(1+cos(angle)));
	  double diff = line-larc;
	  double fl = (USWKpos-vtx).Mag()+(CVCpos-USWKpos).Mag()-diff;
	  fl -= 1.5;
	  double beta = fl/(tof*100.*Const);
	  //	std::cout<<"CVC seg"<<CVCseg<<" beta="<<beta<<std::endl;
	  h1 = getTH("CVC_overbeta"), h1-> Fill(1./beta);
	  if( GeomTools::GetID(vtx)==CID_Fiducial ){
	    h1 = getTH("CVC_overbeta_wtar"), h1-> Fill(1./beta);
	  }
	}
	
	if( bltrackMan-> ntrackFDC1()==1 && fCVC_hit.size()==0 && fPC_hit.size()==1 ){
	  LocalTrack *trackFDC = bltrackMan->trackFDC1(0);
	  conf-> GetBLDCWireMapManager()-> GetGParam(CID_FDC1, gpos, gdir);
	  TVector3 FDCpos = trackFDC->GetPosatZ(gpos.Z());
	  TVector3 dir = (FDCpos-vtx).Unit(); 
	  TVector3 USWKpos = FDCpos+(250-FDCpos.Z())*dir;
	  int PCseg = fPC_hit[0]-> seg();
	  
	  double tof = fPC_hit[0]->ctmean()-fT0_hit[0]->ctmean();
	  tof -= beam_tof;
	  
	  TVector3 PCpos;
	  conf-> GetGeomMapManager()-> GetGPos(CID_PC, PCseg, PCpos);
	  double angle = (USWKpos-vtx).Angle(PCpos-USWKpos);
	  double leff = 100.;
	  double r = leff/sqrt(2*(1-cos(angle)));
	  double larc = r*angle;
	  double line = 2*leff/sqrt(2*(1+cos(angle)));
	  double diff = line-larc;
	  double fl = (USWKpos-vtx).Mag()+(PCpos-USWKpos).Mag()-diff;
	  fl -= 1.5;
	  double beta = fl/(tof*100.*Const);
	  //	std::cout<<"PC seg"<<PCseg<<" beta="<<beta<<std::endl;
	  h1 = getTH("PC_overbeta"), h1-> Fill(1./beta);
	  if( GeomTools::GetID(vtx)==CID_Fiducial ){
	    h1 = getTH("PC_overbeta_wtar"), h1-> Fill(1./beta);
	  }
	}
      }
    }
  }

  if( header->trigmode2(Mode_KCDH1N) && fT0_hit.size()==1 && fBeamMom>0.1 && fBeamPID==Beam_Kaon && bltrackMan->ntrackBPC()==1 ){
    h1 = getTH("Neutral_Reduction"), h1-> Fill(0);
    if( vtx_flag ){
      h1-> Fill(5);
      if( forward_n){
	h1-> Fill(6);
      }
      if( forward_g){
	h1-> Fill(7);
      }
    }
  }

  if( header->IsTrig(Trig_Neutral) && fT0_hit.size()==1 && fBeamMom>0.1 && fBeamPID==Beam_Kaon && bltrackMan->ntrackBPC()==1 ){
    h1 = getTH("Neutral_Reduction"), h1-> Fill(0);
    if( vtx_flag ){
      h1-> Fill(5);
      if( forward_n){
	h1-> Fill(6);
      }
      if( forward_g){
	h1-> Fill(7);
      }
    }
  }
}
