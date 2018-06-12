#include "kn_El.h"

void initHist()
{
  gFile->cd();
  new TH1F("KNpipi_MM_InEl_hit_wD", "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_Cap_hit_wD", "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("km_FL_InEl_npipi_hit_wD",     "K^{-} Fl capture", 10000, 0.0, 100);
  new TH1F("km_dE_InEl_npipi_hit_wD",  "K^{-} dE Inelastic", 1000,  0.0, 1.0);
  new TH1F("km_fmom_InEl_wD", "K^{-} mom", 1000, 0.0, 1.0);

  new TH1F("pi_mom1_Cap_wD", "#pi mom", 1000, 0.0, 1.0);
  new TH1F("pi_mom2_Cap_wD", "#pi mom", 1000, 0.0, 1.0);

  new TH1F("km_FL_Cap_npipi_hit", "K^{-} FL",  1000, 0.0, 100.);
  new TH1F("km_FL_InEl_npipi_hit", "K^{-} FL", 1000, 0.0, 100.);
  new TH1F("km_FL_InEl_npipi_hit_2tra", "K^{-} FL", 1000, 0.0, 100.);
  new TH1F("km_FL_InEl_npipi_hit_Spi", "K^{-} FL", 1000, 0.0, 100.);
  new TH1F("km_FL_Decay_npipi_hit", "K^{-} FL", 1000, 0.0, 100.);

  new TH2F("km_fpos_XY_Cap_npipi_hit", "K^{-} fnail pos", 2000, -100, 100, 2000, -100, 100);
  new TH1F("km_fpos_Z_Cap_npipi_hit", "K^{-} fnail pos",  2000, -100, 100);

  new TH2F("km_fpos_XY_InEl_npipi_hit", "K^{-} fnail pos", 2000, -100, 100, 2000, -100, 100);
  new TH1F("km_fpos_Z_InEl_npipi_hit", "K^{-} fnail pos",  2000, -100, 100);

  new TH1F("KNpipi_MM_Cap_wD", "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_InEl_wD", "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);

  new TH1F("beam_mom", "beam mom", 1000, 0.5, 1.5);

  new TH2F("n_theta_mom",  "n cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("p_theta_mom",  "p cos#theta vs mom",     1000, -1.0, 1.0, 1000, 0.0, 1.0);

  new TH2F("n_theta_mom_n_hit",  "n cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom_n_hit", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);

  new TH2F("n_theta_mom_npipi_hit",  "n cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom_npipi_hit", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);

  new TH1F("npipi_hit_ev",  "n #pi^{+} #pi^{-} hit event evaluation", 1000, 0, 1000);
  new TH1F("ntraNpipiCap",  "n track n pi pi hit K^{-} Capture", 100, 0, 100);
  new TH1F("ntraNpipiInEl", "n track n pi pi hit K^{-} Inelastic", 100, 0, 100);
  new TH1F("KNpipi_MM_InEl_hit", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0, 2.0);
  new TH1F("KNpipi_MM_InEl_hit_2tra", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0, 2.0);
  new TH1F("KNpipi_MM_InEl_hit_Spi",  "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0, 2.0);
  new TH1F("KNpipi_MM_Cap_hit",  "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0, 2.0);
  new TH1F("KNpipi_MM_Decay_hit",  "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0, 2.0);

  new TH2F("km_theta_mom_InEl", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("km_theta_mom_kCap", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);

  new TH2F("km_theta_mom_InEl2", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("km_theta_mom_kCap_S", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("km_theta_mom_kCap_Spi", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH1F("km_mom_Cap", "K^{-} mom", 1000, 0.0, 0.5);
  new TH1F("km_mom_Cap_S", "K^{-} mom", 1000, 0.0, 0.5);
  new TH1F("km_mom_Cap_Spi", "K^{-} mom", 1000, 0.0, 0.5);

  new TH1F("km_dE_Cap",     "K^{-} dE capture", 1000,  0.0, 1.0);
  new TH1F("km_dE_Cap_S",   "K^{-} dE capture", 1000,  0.0, 1.0);
  new TH1F("km_dE_Cap_Spi", "K^{-} dE capture", 1000,  0.0, 1.0);
  new TH1F("km_dE_Cap_wD", "K^{-} dE capture", 1000,  0.0, 1.0);

  new TH1F("km_dE_InEl",     "K^{-} dE Inleastic", 1000,  0.0, 1.0);
  new TH1F("km_dE_InEl_S",   "K^{-} dE Inelastic", 1000,  0.0, 1.0);
  new TH1F("km_dE_InEl_Spi", "K^{-} dE Inelastic", 1000,  0.0, 1.0);
  new TH1F("km_dE_InEl_wD",  "K^{-} dE Inelastic", 1000,  0.0, 1.0);

  new TH1F("km_FL_Cap",     "K^{-} Fl capture", 10000, 0.0, 100);
  new TH1F("km_FL_Cap_S",   "K^{-} Fl capture", 10000, 0.0, 100);
  new TH1F("km_FL_Cap_Spi", "K^{-} Fl capture", 10000, 0.0, 100);

  new TH1F("km_FL_InEl",     "K^{-} Fl capture", 10000, 0.0, 100);
  new TH1F("km_FL_InEl_S",   "K^{-} Fl capture", 10000, 0.0, 100);
  new TH1F("km_FL_InEl_Spi", "K^{-} Fl capture", 10000, 0.0, 100);

  new TH2F("Sm_theta_mom", "#Sigma^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("S0_theta_mom", "#Sigma^{0} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("Sp_theta_mom", "#Sigma^{p} cos#theta vs mom", 1000, -1.0, 1.0, 1000, 0.0, 1.0);

  new TH1F("Km_Reac", "K^{-} reaction", 100, 0., 100);
  new TH1F("ntraKInEl", "ntrack K^{-} Inelastic", 100, 0, 100);
  new TH1F("ntraKInEl_S", "ntrack K^{-} Inelastic", 100, 0, 100);
  new TH1F("ntraKInEl_Spi", "ntrack K^{-} Inelastic", 100, 0, 100);
  new TH1F("ntraKInEl_wD",  "ntrack K^{-} Inelastic", 100, 0, 100);

  new TH1F("ntraKCap",  "ntrack K^{-} Capture", 100, 0, 100);
  new TH1F("ntraKCap_S",  "ntrack K^{-} Capture", 100, 0, 100);
  new TH1F("ntraKCap_Spi",  "ntrack K^{-} Capture", 100, 0, 100);
  new TH1F("ntraKCap_wD",  "ntrack K^{-} Capture", 100, 0, 100);
}

int fillHist(DetectorData* detData, MCData *mcData, ReactionData *reacData)
{
  TH1F *h1;
  TH2F *h2;
  if( reacData->ReactionID()!=821 ){
    std::cout<<" !!! not k-n elastic event !!!"<<std::endl;
    return 0;
  }
  int pim_id = 0;
  int pip_id = 0;
  bool n_hit = false;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData-> detectorHit(i);
    if( hit->parentID()==0 && hit->pdg()==2112 && hit->detectorID()==CID_NC  ) n_hit=true;
    if( hit->pdg()==-211 && hit->detectorID()==CID_CDH ) pim_id=hit->trackID();
    if( hit->pdg()== 211 && hit->detectorID()==CID_CDH ) pip_id=hit->trackID();
  }

  Track *sigma_track = 0;
  Track *km_track    = 0;
  Track *n_track     = 0;
  Track *p_track     = 0;
  Track *beam_track  = 0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track = mcData-> track(i);

    if( track->pdgID()==3222 ){
      TVector3 mom = track->momentum();
      sigma_track=track;
      h2 = (TH2F*)gFile-> Get("Sp_theta_mom"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
    }
    if( track->pdgID()==3112 ){
      TVector3 mom = track->momentum();
      sigma_track=track;
      h2 = (TH2F*)gFile-> Get("Sm_theta_mom"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
    }
    if( track->pdgID()==3212 ){
      TVector3 mom = track->momentum();
      sigma_track=track;
      h2 = (TH2F*)gFile-> Get("S0_theta_mom"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
    }
    if( track->pdgID()==321  && track->parentTrackID()==0 ){
      beam_track=track;
      TVector3 mom = 0.001*beam_track->momentum();
      h1 = (TH1F*)gFile-> Get("beam_mom"), h1-> Fill(mom.Mag());
    }
    if( track->pdgID()==2112 && track->parentTrackID()==0){
      n_track=track;
      TVector3 mom = track->momentum();
      h2 = (TH2F*)gFile-> Get("n_theta_mom"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
      if( n_hit ) h2 = (TH2F*)gFile-> Get("n_theta_mom_n_hit"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
      if( n_hit && pim_id>0 && pip_id>0 ) h2 = (TH2F*)gFile-> Get("n_theta_mom_npipi_hit"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
    }
    if( track->pdgID()==-321 && track->parentTrackID()==0){
      km_track=track;
      TVector3 mom = track->momentum();
      h2 = (TH2F*)gFile-> Get("km_theta_mom"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
      if( n_hit ) h2 = (TH2F*)gFile-> Get("km_theta_mom_n_hit"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
      if( n_hit && pim_id>0 && pip_id>0) h2 = (TH2F*)gFile-> Get("km_theta_mom_npipi_hit"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
    }
  
    if( track->pdgID()==2212 && track->parentTrackID()==0){
      p_track=track;
      TVector3 mom = track->momentum();
      h2 = (TH2F*)gFile-> Get("p_theta_mom"), h2-> Fill(mom.CosTheta(), 0.001*mom.Mag());
    }
  }

  if( n_track ){
    std::cout<<"> K- reaction : "<<SimTools::GetProcessName(km_track->finalProcessID())<<std::endl;
  }

  TString mat_name = "";
  int mat_id=-999;
  std::vector<Track*> tracks;
  std::vector<Track*> tracks2;
  bool reacD_flag = false;
  if( km_track ){
    TVector3 vtx = 0.1*km_track->finalPos();
    mat_id=GeomTools::GetIDMat(vtx, mat_name);
    //    std::cout<<"> K react  ID="<<mat_id<<"  "<<mat_name<<std::endl;
    if( mat_id==151 || mat_id==160 ) reacD_flag=true;
    //    std::cout<<"> K- final process : "<<SimTools::GetProcessName(km_track->finalProcessID())<<std::endl;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->parentTrackID()==km_track->trackID() ){
	tracks2.push_back(track);
	if( track->vertex()==km_track->finalPos() ) tracks.push_back(track);
	//	SimTools::Dump(track);
      }
    }
    // if( tracks.size()!=tracks2.size() ){
    //   std::cout<<" !!! tracks2.size()="<<tracks2.size()<<" tracks.size()="<<tracks.size()<<" !!!"<<std::endl;
    //   for( int i=0; i<tracks2.size(); i++ ) SimTools::Dump(tracks2[i], 1);
    // }
    
    if( km_track->finalProcessID()==102 ){
      //      std::cout<<"> K Inelastic react  ID="<<mat_id<<"  "<<mat_name<<std::endl;
      h1 = (TH1F*)gFile-> Get("ntraKInEl"), h1-> Fill(tracks.size());
      TVector3 pos = km_track->vertex();
      TVector3 fpos = km_track->finalPos();
      TVector3 mom = 0.001*km_track-> momentum();
      //      TVector3 fmom = 0.001*km_track-> finalMom();
      TVector3 fmom;
      for( int i=0; i<tracks.size(); i++ ) fmom += tracks[i]->momentum();
      fmom *= 0.001;
      TLorentzVector km_lmom, km_lmom2;
      km_lmom.SetVectM(mom, kpMass);
      km_lmom2.SetVectM(fmom, kpMass);
      double dE = km_lmom.E()-km_lmom2.E();
      double fl = 0.1*(pos-fpos).Mag();

      h2 = (TH2F*)gFile-> Get("km_theta_mom_InEl"), h2-> Fill(mom.CosTheta(), mom.Mag());
      h1 = (TH1F*)gFile-> Get("km_FL_InEl"), h1-> Fill(fl);
      if( reacD_flag ){
	//	std::cout<<"===== K- Inelastic w/ D final mom : "<<fmom.Mag()<<std::endl;
	// std::cout<<"===== K- Inelastic w/ D ====="<<std::endl;
	// for( int i=0; i<tracks.size(); i++ ) SimTools::Dump(tracks[i]);
	// std::cout<<"============================="<<std::endl;
      	h1 = (TH1F*)gFile-> Get("ntraKInEl_wD"), h1-> Fill(tracks.size());
	h1 = (TH1F*)gFile-> Get("km_dE_InEl_wD"), h1-> Fill(dE);
	TVector3 sum;
	for( int i=0; i<tracks.size(); i++ ) sum += tracks[i]->momentum();
	h1 = (TH1F*)gFile-> Get("km_fmom_InEl_wD"), h1-> Fill(0.001*sum.Mag());
      }
    }
  
    if( km_track->finalProcessID()==204 ){
      h1 = (TH1F*)gFile-> Get("ntraKCap"), h1-> Fill(tracks.size());
      TVector3 pos = km_track->vertex();
      TVector3 fpos = km_track->finalPos();
      TVector3 mom = 0.001*km_track-> momentum();
      //      TVector3 fmom = 0.001*km_track-> finalMom();
      TVector3 fmom;
      for( int i=0; i<tracks.size(); i++ ) fmom += tracks[i]->momentum();
      fmom *= 0.001;
      TLorentzVector km_lmom, km_lmom2;
      km_lmom.SetVectM(mom, kpMass);
      km_lmom2.SetVectM(fmom, kpMass);
      double dE = km_lmom.E()-km_lmom2.E();
      double fl = 0.1*(pos-fpos).Mag();

      h2 = (TH2F*)gFile-> Get("km_theta_mom_kCap"), h2-> Fill(mom.CosTheta(), mom.Mag());
      h1 = (TH1F*)gFile-> Get("km_FL_Cap"), h1-> Fill(fl);
      h1 = (TH1F*)gFile-> Get("km_mom_Cap"), h1-> Fill(mom.Mag());
      if( reacD_flag ){
	//	std::cout<<"===== K- Capture fmom : "<<fmom.Mag()<<std::endl;
	if( tracks.size()>2 ){
	  // std::cout<<"===== K- Capture w/ Deuteron  track:"<<tracks.size()<<" ====="<<std::endl;
	  // for( int i=0; i<tracks.size(); i++ ) SimTools::Dump(tracks[i], 1);
	  // std::cout<<"============================================"<<std::endl;
	}
	h1 = (TH1F*)gFile-> Get("ntraKCap_wD"), h1-> Fill(tracks.size());
	h1 = (TH1F*)gFile-> Get("km_dE_Cap_wD"), h1-> Fill(dE);

	if( sigma_track ){
	  Track *pi_track1 = 0;
	  Track *pi_track2 = 0;
	  for( int i=0; i<mcData->trackSize(); i++ ){
	    Track *track = mcData->track(i);
	    if( track->parentTrackID()==sigma_track->trackID() ){
	      if( track->pdgID()==-211 || track->pdgID()==111 || track->pdgID()==211 ) pi_track2=track;
	    }
	    if( track->parentTrackID()==sigma_track->parentTrackID() ){
	      if( track->pdgID()==-211 || track->pdgID()==111 || track->pdgID()==211 ) pi_track1=track;
	    }
	  }
	  if( pi_track1 && pi_track2 ){
	    TVector3 mom1 = 0.001*pi_track1->momentum();
	    TVector3 mom2 = 0.001*pi_track2->momentum();

	    TLorentzVector beam_lmom, pi_lmom1, pi_lmom2, n_lmom;
	    TLorentzVector tgt_lmom(0, 0, 0, dMass);
	    beam_lmom.SetVectM(-0.001*beam_track->momentum(), kpMass);
	    n_lmom.SetVectM(0.001*n_track->momentum(), nMass);
	    pi_lmom1.SetVectM(mom1, piMass);
	    pi_lmom2.SetVectM(mom2, piMass);
	    double mm_knpipi = (beam_lmom+tgt_lmom-pi_lmom1-pi_lmom2-n_lmom).M();

	    h1 = (TH1F*)gFile-> Get("pi_mom1_Cap_wD"), h1-> Fill(mom1.Mag());
	    h1 = (TH1F*)gFile-> Get("pi_mom2_Cap_wD"), h1-> Fill(mom2.Mag());
	    h1 = (TH1F*)gFile-> Get("KNpipi_MM_Cap_wD"), h1-> Fill(mm_knpipi);
	  }
	  else{
	    if( sigma_track->pdgID()!=3212 && sigma_track->finalProcessID()!=204 ){
	      SimTools::Dump(sigma_track);
	      std::cout<<" !!! pi_track1="<<pi_track1<<" pi_track2="<<pi_track2<<" !!!"<<std::endl;
	    }
	  }
	}
      }
      //      std::cout<<"> K- captrue "<<std::endl;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( km_track->trackID()==track->parentTrackID() ){
	  //	  SimTools::Dump(track);
	}
      }
    }
  }
  else{
    std::cout<<"  !!! K- track not found !!!"<<std::endl;
  }

  int n_pi=0;
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track = mcData->track(i);
    if( track->pdgID()==211 || track->pdgID()==111 || track->pdgID()==-211 ) n_pi++;
  }

  Track *pi_track=0;
  if( sigma_track ){
    if( sigma_track->pdgID()==3222 ){
      //      std::cout<<"===== S+ pi- Event ===="<<std::endl;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData-> track(i);
	if( track->pdgID()==-211 ) pi_track=track;
      }
    }
    if( sigma_track->pdgID()==3212 ){
      //      std::cout<<"===== S0 pi0 Event ===="<<std::endl;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData-> track(i);
	if( track->pdgID()==111 ) pi_track=track;
      }
    }
    if( sigma_track->pdgID()==3112 ){
      //      std::cout<<"===== S- pi+ Event ===="<<std::endl;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData-> track(i);
	if( track->pdgID()==211 ) pi_track=track;
      }
    }

    if( sigma_track->processID()==102 ){
      //      std::cout<<" K- Inelastic  /w ID="<<mat_id<<" "<<mat_name<<std::endl;
      TVector3 mom = 0.001*km_track-> momentum();
      //      TVector3 fmom = 0.001*km_track-> finalMom();
      TVector3 fmom;
      for( int i=0; i<tracks.size(); i++ ) fmom += tracks[i]->momentum();
      fmom *= 0.001;
      TLorentzVector km_lmom1, km_lmom2;
      km_lmom1.SetVectM(mom, kpMass);
      km_lmom2.SetVectM(fmom, kpMass);
      double dE = km_lmom1.E()-km_lmom2.E();
      TVector3 pos = km_track->vertex();
      TVector3 fpos = km_track->finalPos();
      double fl = 0.1*(pos-fpos).Mag();

      h2 = (TH2F*)gFile-> Get("km_theta_mom_InEl"), h2-> Fill(mom.CosTheta(), mom.Mag());
      h1 = (TH1F*)gFile-> Get("km_FL_InEl_S"), h1-> Fill(fl);
      h1 = (TH1F*)gFile-> Get("ntraKInEl_S"), h1-> Fill(tracks.size());
      if( pi_track ){
	h1 = (TH1F*)gFile-> Get("ntraKInEl_Spi"), h1-> Fill(tracks.size());
	h1 = (TH1F*)gFile-> Get("km_FL_InEl_Spi"), h1-> Fill(fl);
      }
    }
    else if( sigma_track->processID()==204 ){
      //      std::cout<<" K- Capture    /w ID="<<mat_id<<" "<<mat_name<<"  Secondary"<<tracks.size()<<std::endl;

      TVector3 mom = 0.001*km_track-> momentum();
      //      TVector3 fmom = 0.001*km_track-> finalMom();
      TVector3 fmom;
      for( int i=0; i<tracks.size(); i++ ) fmom += tracks[i]->momentum();
      fmom *= 0.001;
      TLorentzVector km_lmom1, km_lmom2;
      km_lmom1.SetVectM(mom, kpMass);
      km_lmom2.SetVectM(fmom, kpMass);
      double dE = km_lmom1.E()-km_lmom2.E();
      TVector3 pos = km_track->vertex();
      TVector3 fpos = km_track->finalPos();
      double fl = 0.1*(pos-fpos).Mag();
      h2 = (TH2F*)gFile-> Get("km_theta_mom_kCap_S"), h2-> Fill(mom.CosTheta(), mom.Mag());
      h1 = (TH1F*)gFile-> Get("km_mom_Cap_S"), h1-> Fill(fmom.Mag());
      h1 = (TH1F*)gFile-> Get("km_dE_Cap_S"), h1-> Fill(dE);
      h1 = (TH1F*)gFile-> Get("km_FL_Cap_S"), h1-> Fill(fl);
      h1 = (TH1F*)gFile-> Get("ntraKCap_S"), h1-> Fill(tracks.size());

      if( pi_track ){
	h2 = (TH2F*)gFile-> Get("km_theta_mom_kCap_Spi"), h2-> Fill(mom.CosTheta(), mom.Mag());
	h1 = (TH1F*)gFile-> Get("km_mom_Cap_Spi"), h1-> Fill(fmom.Mag());
	h1 = (TH1F*)gFile-> Get("km_dE_Cap_Spi"), h1-> Fill(dE);
	h1 = (TH1F*)gFile-> Get("km_FL_Cap_Spi"), h1-> Fill(fl);
	h1 = (TH1F*)gFile-> Get("ntraKCap_Spi"), h1-> Fill(tracks.size());
      }
    }

    //    std::cout<<">>> "<<SimTools::GetProcessName(sigma_track->processID())<<std::endl;
  }

  if( pim_id>0 && pip_id>0 && n_hit ){
    Track *pim_track=0;
    Track *pip_track=0;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->trackID()==pim_id ) pim_track=track;
      if( track->trackID()==pip_id ) pip_track=track;
    }
    h1 = (TH1F*)gFile-> Get("npipi_hit_ev");
    if( km_track->finalProcessID()==102 ) h1-> Fill(490);
    else if( km_track->finalProcessID()==204 ) h1-> Fill(491);
    else h1-> Fill(km_track->finalProcessID());

    if( km_track->finalProcessID()==1 ){
      std::cout<<"===== K- decay event pi+ pi- n hit ====="<<std::endl;
      SimTools::Dump(km_track);
      for( int i=0; i<tracks2.size(); i++ ) SimTools::Dump(tracks2[i]);
      std::cout<<"===== Dump hit pi+ pi- track ====="<<std::endl;
      SimTools::Dump(pim_track);
      SimTools::Dump(pip_track);
    }
    TVector3 pos  = 0.1*km_track->vertex();
    TVector3 fpos = 0.1*km_track->finalPos();
    double fl = (pos-fpos).Mag();

    TVector3 n_mom = 0.001*n_track-> momentum();
    TVector3 pim_mom = 0.001*pim_track-> momentum();
    TVector3 pip_mom = 0.001*pip_track-> momentum();
    TVector3 beam_mom = -0.001*beam_track-> momentum();
    TLorentzVector n_lmom, pim_lmom, pip_lmom, beam_lmom;
    n_lmom.SetVectM(n_mom, nMass);
    pim_lmom.SetVectM(pim_mom, piMass);
    pip_lmom.SetVectM(pip_mom, piMass);
    beam_lmom.SetVectM(beam_mom, kpMass);
    TLorentzVector tgt_lmom(0, 0, 0, dMass);
    double mm_knpipi=(beam_lmom+tgt_lmom-n_lmom-pim_lmom-pip_lmom).M();

    TVector3 km_mom0 = 0.001*km_track->momentum();
    TVector3 km_mom1 = 0.001*km_track->finalMom();
    TLorentzVector km_lmom0, km_lmom1;
    km_lmom0.SetVectM(km_mom0, kpMass);
    km_lmom1.SetVectM(km_mom1, kpMass);
    double dE = km_lmom0.E()- km_lmom1.E();

    if( km_track->finalProcessID()==102 ){
      h1 = (TH1F*)gFile-> Get("km_FL_InEl_npipi_hit"), h1-> Fill(fl);
      h2 = (TH2F*)gFile-> Get("km_fpos_XY_InEl_npipi_hit"), h2-> Fill(fpos.X(), fpos.Y());
      h1 = (TH1F*)gFile-> Get("km_fpos_Z_InEl_npipi_hit"), h1-> Fill(fpos.Z());
      h1 = (TH1F*)gFile-> Get("ntraNpipiInEl"), h1-> Fill(tracks.size());
      h1 = (TH1F*)gFile-> Get("KNpipi_MM_InEl_hit"), h1-> Fill(mm_knpipi);
      if( tracks.size()==2 ){
	h1 = (TH1F*)gFile-> Get("KNpipi_MM_InEl_hit_2tra"), h1-> Fill(mm_knpipi);
	h1 = (TH1F*)gFile-> Get("km_FL_InEl_npipi_hit_2tra"), h1-> Fill(fl);
      }
      if( sigma_track && pi_track ){
	h1 = (TH1F*)gFile-> Get("KNpipi_MM_InEl_hit_Spi"), h1-> Fill(mm_knpipi);
	h1 = (TH1F*)gFile-> Get("km_FL_InEl_npipi_hit_Spi"), h1-> Fill(fl);
      }
      if( reacD_flag ){
	h1 = (TH1F*)gFile-> Get("KNpipi_MM_InEl_hit_wD"), h1-> Fill(mm_knpipi);
	h1 = (TH1F*)gFile-> Get("km_FL_InEl_npipi_hit_wD"), h1-> Fill(fl);
	h1 = (TH1F*)gFile-> Get("km_dE_InEl_npipi_hit_wD"), h1-> Fill(dE);
      }
    }
    if( km_track->finalProcessID()==204 ){
      h1 = (TH1F*)gFile-> Get("km_FL_Cap_npipi_hit"), h1-> Fill(fl);
      h2 = (TH2F*)gFile-> Get("km_fpos_XY_Cap_npipi_hit"), h2-> Fill(fpos.X(), fpos.Y());
      h1 = (TH1F*)gFile-> Get("km_fpos_Z_Cap_npipi_hit"), h1-> Fill(fpos.Z());
      h1 = (TH1F*)gFile-> Get("ntraNpipiCap"), h1-> Fill(tracks.size());
      h1 = (TH1F*)gFile-> Get("KNpipi_MM_Cap_hit"), h1-> Fill(mm_knpipi);
      if( reacD_flag ){
	h1 = (TH1F*)gFile-> Get("KNpipi_MM_Cap_hit_wD"), h1-> Fill(mm_knpipi);
      }
    }
    if( km_track->finalProcessID()==1 ){
      h1 = (TH1F*)gFile-> Get("km_FL_Decay_npipi_hit"), h1-> Fill(fl);
      h1 = (TH1F*)gFile-> Get("KNpipi_MM_Decay_hit"), h1-> Fill(mm_knpipi);
    }


#if 0
    std::cout<<"===== #pi+ pi- hit event ====="<<std::endl;
    SimTools::Dump(km_track, 1);
    SimTools::Dump(p_track, 1);
    std::cout<<">>>>> K- Reaction Dumpping"<<std::endl;
    for( int i=0; i<tracks.size(); i++ ) SimTools::Dump(tracks[i], 1);
    std::cout<<">>>>> K- create track : "<<tracks.size()<<std::endl;
    std::cout<<">>> d(K-, n pi+ pi-) : "<<mm_knpipi<<" (PDG: "<<nMass<<")"<<std::endl;
    SimTools::Dump(pim_track);
    SimTools::Dump(pip_track);
    SimTools::Dump(n_track);
    std::cout<<"=============================="<<std::endl;
#endif
  }
}

void writeHist()
{
  std::cout<<"===== Write Histgram ====="<<std::endl;
  gFile-> Write();
  gFile-> Close();
}
