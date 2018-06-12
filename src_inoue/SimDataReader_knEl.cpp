#include "SimDataReader.h"
#include "KnuclRootData.h"

using namespace std;

void SimDataReader::initHist_knEl()
{
  new TH2F("knEl_event", "K^{-} n elastic scattering event", 20, 0, 20, 5, 0, 5);

  new TH2F("knEl_pipi_pos_XY_wD", "#pi^{+} #pi^{-} Vtx XY", 1000, -25, 25, 1000, -25, 25);
  new TH2F("knEl_pipi_pos_ZX_wD", "#pi^{+} #pi^{-} Vtx XZ", 1000,  50, 50, 1000, -25, 25);
  new TH2F("knEl_pipi_pos_ZY_wD", "#pi^{+} #pi^{-} Vtx YZ", 1000,  50, 50, 1000, -25, 25);

  new TH2F("knEl_pipi_pos_XY_woD", "#pi^{+} #pi^{-} Vtx XY", 1000, -25, 25, 1000, -25, 25);
  new TH2F("knEl_pipi_pos_ZX_woD", "#pi^{+} #pi^{-} Vtx XZ", 1000,  50, 50, 1000, -25, 25);
  new TH2F("knEl_pipi_pos_ZY_woD", "#pi^{+} #pi^{-} Vtx YZ", 1000,  50, 50, 1000, -25, 25);

  new TH1F("pim_mom", "#pi^{-} mom", 1000, 0, 1.0);
  new TH1F("pip_mom", "#pi^{+} mom", 1000, 0, 1.0);
  new TH1F("pim_mom_InEl_wD",  "#pi^{-} mom", 1000, 0, 1.0);
  new TH1F("pip_mom_InEl_wD",  "#pi^{+} mom", 1000, 0, 1.0);
  new TH1F("pim_mom_KCap_woD", "#pi^{-} mom", 1000, 0, 1.0);
  new TH1F("pip_mom_KCap_woD", "#pi^{+} mom", 1000, 0, 1.0);

  new TH1F("ntraKInEl_n_hit_wD",  "ntrack from K^{-} Inelastic", 10, 0, 10);
  new TH1F("ntraKCap_n_hit_woD", "ntrack from K^{-} Capture",    50, 0, 50);

  new TH2F("km_theta_mom", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("n_theta_mom",  "n cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("p_theta_mom",  "p cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);

  new TH2F("km_theta_mom_n_hit",     "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom_npipi_hit", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("n_theta_mom_n_hit",  "n cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("p_theta_mom_n_hit",  "p cos#theta vs mom",     1000, -1.0, 1.0, 1500, 0.0, 1.5);

  new TH2F("km_theta_mom_n_hit_Inel",    "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom_n_hit_Inel_wD", "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom_n_hit_Cap",     "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);
  new TH2F("km_theta_mom_n_hit_Cap_wD",  "K^{-} cos#theta vs mom", 1000, -1.0, 1.0, 1500, 0.0, 1.5);

  new TH1F("KNpipi_MM_knEl",         "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_knEl_InEl",    "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_knEl_InEl_wD", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_knEl_Cap",     "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_knEl_Cap_wD",  "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_knEl_Decay",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);

  new TH1F("KN_MM_woAll_wN_knEl",         "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_woAll_wN_knEl_InEl",    "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_woAll_wN_knEl_InEl_wD", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_woAll_wN_knEl_Cap",     "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_woAll_wN_knEl_Cap_wD",  "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_woAll_wN_knEl_Decay",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
}

void SimDataReader::fillHist_knEl()
{
  TH1F *h1;
  TH2F *h2;
  TH2F *h2_ev = (TH2F*)gFile-> Get("knEl_event");

  Track *n_track=0;
  Track *km_track=0;
  Track *p_track=0;
  Track *beam_track=0;
  std::vector<Track*> tracks;

  //****************************************************************************************************//
  //*** Serch Track & Reaction                                                                       ***//
  //****************************************************************************************************//

  Track *sigma_track=0;
  Track *pi_track1=0;
  Track *pi_track2=0;

  int km_fpos_id=-999;
  TString km_fpos_mat;
  bool reacD_flag = false;

  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track = mcData->track(i);
    if( track->parentTrackID()==0 ){
      if( track->pdgID()==-321 ) km_track=track;
      if( track->pdgID()==2112 ) n_track=track;
      if( track->pdgID()==2212 ) p_track=track;
      if( track->pdgID()==321  ) beam_track=track;
    }

    if( track->pdgID()==3222 || track->pdgID()==3112 || track->pdgID()==3212 || track->pdgID()==3122) sigma_track=track;
  }

  if( km_track ){
    km_fpos_id=GeomTools::GetIDMat(0.1*km_track->finalPos(), km_fpos_mat);
    if( km_fpos_id==151 || km_fpos_id==160 ) reacD_flag = true;

    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->parentTrackID()==km_track->trackID() ){
	if( track->vertex()==km_track->finalPos() ) tracks.push_back(track);
      }
    }
  }
  else std::cout<<"  !!! K- track not found !!!"<<std::endl;

  if( sigma_track ){
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->parentTrackID()==sigma_track->parentTrackID() ){
	if( track->pdgID()==-211 || track->pdgID()==111 || track-> pdgID()==211 ) pi_track1=track;
      }
      if( track->parentTrackID()==sigma_track->trackID() ){
	if( track->pdgID()==-211 || track->pdgID()==111 || track-> pdgID()==211 ) pi_track2=track;
      }
    }
  }

  //****************************************************************************************************//
  //*** Serch DetectorHit                                                                            ***//
  //****************************************************************************************************//

  DetectorHit *n_hit=0;
  DetectorHit *pim_hit=0;
  DetectorHit *pip_hit=0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit-> parentID()==0 && hit-> pdg()==2112 && hit-> detectorID()==CID_NC ) n_hit=hit;
    if( hit->pdg()==-211 && hit->detectorID()==CID_CDH ) pim_hit=hit;
    if( hit->pdg()== 211 && hit->detectorID()==CID_CDH ) pip_hit=hit;
  }

  bool npipi_hit = false;
  if( n_hit && pim_hit && pip_hit ) npipi_hit=true;

  TVector3 km_final_pos(DBL_MAX, DBL_MAX, DBL_MAX);
  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track = mcData->track(i);
    if( track->parentTrackID()==0 && track->pdgID()==-321 ){
      km_final_pos = 0.1*track->finalPos();
      if( npipi_hit ){
	cout<<"n pi+ pi- Reaction "<<getProcessName(track->finalProcessID())<<endl;
      }
    }
  }

  // TVector3 pip_vtx(DBL_MAX, DBL_MAX, DBL_MAX);
  // for( int i=0; i<mcData->trackSize(); i++ ){
  //   Track *track  = mcData->track(i);
  //   if( track->pdgID()==211 ){
  //     Track *parent_track = getMCTrack(track->parentTrackID());
  //     while( true ){
  // 	if( parent_track->pdgID()!=211 ){
  // 	  pip_vtx=track->vertex();
  // 	  break;
  // 	}
  // 	track = parent_track;
  // 	parent_track = getMCTrack(track->parentTrackID());
  //     }
  //     break;
  //   }
  // }

  // TVector3 pim_vtx(DBL_MAX, DBL_MAX, DBL_MAX);
  // for( int i=0; i<mcData->trackSize(); i++ ){
  //   Track *track  = mcData->track(i);
  //   if( track->pdgID()==-211 ){
  //     Track *parent_track = getMCTrack(track->parentTrackID());
  //     while( true ){
  // 	if( parent_track->pdgID()!=-211 ){
  // 	  pim_vtx=track->vertex();
  // 	  break;
  // 	}
  // 	track = parent_track;
  // 	parent_track = getMCTrack(track->parentTrackID());
  //     }
  //     break;
  //   }
  // }


  h2_ev-> Fill(1., 0.);
  if( km_track->finalProcessID()==102 ){
    h2_ev-> Fill(2., 0.);
    if( reacD_flag ) h2_ev-> Fill(3., 0.);
    else h2_ev-> Fill(4., 0.);
  }
  else if( km_track->finalProcessID()==204 ){
    h2_ev-> Fill(5., 0.);
    if( reacD_flag ) h2_ev-> Fill(6., 0.);
    else h2_ev-> Fill(7., 0.);
  }
  else if( km_track->finalProcessID()==1 ){
    h2_ev-> Fill(8., 0.);
  }
  else{
    h2_ev-> Fill(9., 0.);
  }

  h2 = (TH2F*)gFile-> Get("km_theta_mom"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
  h2 = (TH2F*)gFile-> Get("n_theta_mom"), h2-> Fill(n_track->momentum().CosTheta(), 0.001*n_track->momentum().Mag());
  h2 = (TH2F*)gFile-> Get("p_theta_mom"), h2-> Fill(p_track->momentum().CosTheta(), 0.001*p_track->momentum().Mag());

  if( n_hit ){
    h2_ev-> Fill(1., 1.);
    if( npipi_hit ) h2_ev-> Fill(1., 2.);
    h2 = (TH2F*)gFile-> Get("km_theta_mom_n_hit"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
    if( npipi_hit ){
      h2 = (TH2F*)gFile-> Get("km_theta_mom_npipi_hit"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
      TString mat_name;
      int geom_id=GeomTools::GetIDMat(km_final_pos, mat_name);
      cout<<"   mat : "<<mat_name<<endl;
      cout<<" K- final pos ("<<km_final_pos.X()<<", "<<km_final_pos.Y()<<", "<<km_final_pos.Z()<<")"<<endl; 
      if( mat_name=="LDeuterium" ){
	h2 = (TH2F*)gFile-> Get("knEl_pipi_pos_XY_wD"), h2-> Fill(km_final_pos.X(), km_final_pos.Y());
	h2 = (TH2F*)gFile-> Get("knEl_pipi_pos_ZX_wD"), h2-> Fill(km_final_pos.Z(), km_final_pos.X());
	h2 = (TH2F*)gFile-> Get("knEl_pipi_pos_ZY_wD"), h2-> Fill(km_final_pos.Z(), km_final_pos.Y());
      }
      else{
	h2 = (TH2F*)gFile-> Get("knEl_pipi_pos_XY_woD"), h2-> Fill(km_final_pos.X(), km_final_pos.Y());
	h2 = (TH2F*)gFile-> Get("knEl_pipi_pos_ZX_woD"), h2-> Fill(km_final_pos.Z(), km_final_pos.X());
	h2 = (TH2F*)gFile-> Get("knEl_pipi_pos_ZY_woD"), h2-> Fill(km_final_pos.Z(), km_final_pos.Y());
      }
    }
    h2 = (TH2F*)gFile-> Get("n_theta_mom_n_hit"), h2-> Fill(n_track->momentum().CosTheta(), 0.001*n_track->momentum().Mag());
    h2 = (TH2F*)gFile-> Get("p_theta_mom_n_hit"), h2-> Fill(p_track->momentum().CosTheta(), 0.001*p_track->momentum().Mag());

    if( km_track->finalProcessID()==102 ){
      h2_ev-> Fill(2., 1.);
      h2 = (TH2F*)gFile-> Get("km_theta_mom_n_hit_Inel"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
      if( npipi_hit ){
	h2_ev-> Fill(2., 2.);
      }
      if( reacD_flag ){
	// std::cout<<"===== K- n elastic  w/ D2 ====="<<std::endl;
	// for( int i=0; i<tracks.size(); i++ ){
	//   SimTools::Dump(tracks[i]);
	// }
	// std::cout<<"==============================="<<std::endl;
	h1 = (TH1F*)gFile-> Get("ntraKInEl_n_hit_wD"), h1-> Fill(tracks.size());

	h2 = (TH2F*)gFile-> Get("km_theta_mom_n_hit_Inel_wD"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
	h2_ev-> Fill(3., 1.);
	if( npipi_hit ){
	  h2_ev-> Fill(3., 2.);
	}
      }
      else{
	h2_ev-> Fill(4., 1.);
	if( npipi_hit ){
	  h2_ev-> Fill(4., 2.);
	}
      }
    }
    else if( km_track->finalProcessID()==204 ){
      h2 = (TH2F*)gFile-> Get("km_theta_mom_n_hit_Cap"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
      h2_ev-> Fill(5., 1.);
      if( npipi_hit ){
	h2_ev-> Fill(5., 2.);
      }
      if( reacD_flag ){
	h2 = (TH2F*)gFile-> Get("km_theta_mom_n_hit_Cap_wD"), h2-> Fill(km_track->momentum().CosTheta(), 0.001*km_track->momentum().Mag());
	h2_ev-> Fill(6., 1.);
	if( npipi_hit ){
	  h2_ev-> Fill(6., 2.);
	}
      }
      else{
	h1 = (TH1F*)gFile-> Get("ntraKCap_n_hit_woD"), h1-> Fill(tracks.size());
	// std::cout<<"===== K- Capture w/o D mat : "<<km_fpos_mat<<" ====="<<std::endl;
	// for( int i=0; i<tracks.size(); i++ ) SimTools::Dump(tracks[i]);
	// std::cout<<"===================================================="<<std::endl;
	h2_ev-> Fill(7., 1.);
	if( npipi_hit ){
	  // std::cout<<"===== K- Capture w/o D mat : "<<km_fpos_mat<<" ====="<<std::endl;
	  // for( int i=0; i<tracks.size(); i++ ) SimTools::Dump(tracks[i]);
	  // std::cout<<"===================================================="<<std::endl;

	  h2_ev-> Fill(7., 2.);
	}
      }
    }
    else if( km_track->finalProcessID()==1 ){
      h2_ev-> Fill(8., 1.);
      if( npipi_hit ){
	h2_ev-> Fill(8., 2.);
      }
    }
    else{
      h2_ev-> Fill(9., 1);
      if( npipi_hit ){
	h2_ev-> Fill(9., 2.);
      }
    }
  }
}


void SimDataReader::fillHist_knEl_data(const TLorentzVector &beam_lmom, const TLorentzVector &n_lmom, const TLorentzVector &pim_lmom, const TLorentzVector &pip_lmom,
				       const TVector3 &pim_vtx, const TVector3 &pip_vtx, 
				       const TVector3 &pim_b_vtx, const TVector3 &pip_b_vtx, const TVector3 &b_pim_vtx, const TVector3  &b_pip_vtx,
				       bool mm_n_flag, bool im_K0_flag, bool im_Sm_flag, bool im_Sp_flag)
{
  CrossSectionTable table = runHeader-> CStable();
  if( table.CSSize()>=10 ) return;

  if( reacData->ReactionID()!=821 ) return;
  const double mm_kn = (beam_lmom+D_LMOM-n_lmom).M();
  const double mm_knpipi = (beam_lmom+D_LMOM-n_lmom-pim_lmom-pip_lmom).M();
  const double mm_knpim = (beam_lmom+D_LMOM-n_lmom-pim_lmom).M();
  const double mm_knpip = (beam_lmom+D_LMOM-n_lmom-pip_lmom).M();

  TH1F *h1;
  TH2F *h2;
  TH2F *h2_ev = (TH2F*)gFile-> Get("knEl_event");

  Track *n_track=0;
  Track *km_track=0;
  Track *p_track=0;
  Track *beam_track=0;
  std::vector<Track*> tracks;

  //****************************************************************************************************//
  //*** Serch Track & Reaction                                                                       ***//
  //****************************************************************************************************//

  Track *sigma_track=0;
  Track *pi_track1=0;
  Track *pi_track2=0;

  int km_fpos_id=-999;
  TString km_fpos_mat;
  bool reacD_flag = false;

  for( int i=0; i<mcData->trackSize(); i++ ){
    Track *track = mcData->track(i);
    if( track->parentTrackID()==0 ){
      if( track->pdgID()==-321 ) km_track=track;
      if( track->pdgID()==2112 ) n_track=track;
      if( track->pdgID()==2212 ) p_track=track;
      if( track->pdgID()==321  ) beam_track=track;
    }

    if( track->pdgID()==3222 || track->pdgID()==3112 || track-> pdgID()==3212 ) sigma_track=track;
  }

  if( km_track ){
    km_fpos_id=GeomTools::GetIDMat(0.1*km_track->finalPos(), km_fpos_mat);
    if( km_fpos_id==151 || km_fpos_id==160 ) reacD_flag = true;

    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->parentTrackID()==km_track->trackID() ){
	if( track->vertex()==km_track->finalPos() ) tracks.push_back(track);
      }
    }
  }
  else std::cout<<"  !!! K- track not found !!!"<<std::endl;

  if( sigma_track ){
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( track->parentTrackID()==sigma_track->parentTrackID() ){
	if( track->pdgID()==-211 || track->pdgID()==111 || track-> pdgID()==211 ) pi_track1=track;
      }
      if( track->parentTrackID()==sigma_track->trackID() ){
	if( track->pdgID()==-211 || track->pdgID()==111 || track-> pdgID()==211 ) pi_track2=track;
      }
    }
  }

  //****************************************************************************************************//
  //*** Serch DetectorHit                                                                            ***//
  //****************************************************************************************************//

  DetectorHit *n_hit=0;
  DetectorHit *pim_hit=0;
  DetectorHit *pip_hit=0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit-> parentID()==0 && hit-> pdg()==2112 && hit-> detectorID()==CID_NC ) n_hit=hit;
    if( hit->pdg()==-211 && hit->detectorID()==CID_CDH ) pim_hit=hit;
    if( hit->pdg()== 211 && hit->detectorID()==CID_CDH ) pip_hit=hit;
  }

  h2_ev-> Fill(1., 3.);
  if( mm_n_flag && !im_K0_flag && !im_Sm_flag && !im_Sp_flag ){
    h1 = (TH1F*)gFile-> Get("KN_MM_woAll_wN_knEl"), h1-> Fill(mm_kn);
  }
  if( km_track->finalProcessID()==102 ){
    h2_ev-> Fill(2., 3.);
    h1 = (TH1F*)gFile-> Get("KNpipi_MM_knEl_InEl"), h1-> Fill(mm_knpipi);
    if( mm_n_flag && !im_K0_flag && !im_Sm_flag && !im_Sp_flag ){
      h1 = (TH1F*)gFile-> Get("KN_MM_woAll_wN_knEl_InEl"), h1-> Fill(mm_kn);
    }
    if( reacD_flag ){
      h1 = (TH1F*)gFile-> Get("KNpipi_MM_knEl_InEl_wD"), h1-> Fill(mm_knpipi);
      if( mm_n_flag && !im_K0_flag && !im_Sm_flag && !im_Sp_flag ){
	h1 = (TH1F*)gFile-> Get("KN_MM_woAll_wN_knEl_InEl_wD"), h1-> Fill(mm_kn);
      }
      h2_ev-> Fill(3., 3.);
    }
    else{
      h2_ev-> Fill(4., 3.);
    }
  }
  else if( km_track->finalProcessID()==204 ){
    h2_ev-> Fill(5., 3.);
    h1 = (TH1F*)gFile-> Get("KNpipi_MM_knEl_Cap"), h1-> Fill(mm_knpipi);
    if( mm_n_flag && !im_K0_flag && !im_Sm_flag && !im_Sp_flag ){
      h1 = (TH1F*)gFile-> Get("KN_MM_woAll_wN_knEl_Cap"), h1-> Fill(mm_kn);
    }
    if( reacD_flag ){
      h1 = (TH1F*)gFile-> Get("KNpipi_MM_knEl_Cap_wD"), h1-> Fill(mm_knpipi);
      if( mm_n_flag && !im_K0_flag && !im_Sm_flag && !im_Sp_flag ){
	h1 = (TH1F*)gFile-> Get("KN_MM_woAll_wN_knEl_Cap_wD"), h1-> Fill(mm_kn);
      }
      h2_ev-> Fill(6., 3.);
    }
    else{
      h2_ev-> Fill(7., 3.);
    }
  }
  else if( km_track->finalProcessID()==1 ){
    h2_ev-> Fill(8., 3.);
    h1 = (TH1F*)gFile-> Get("KNpipi_MM_knEl_Decay"), h1-> Fill(mm_knpipi);
    if( mm_n_flag && !im_K0_flag && !im_Sm_flag && !im_Sp_flag ){
      h1 = (TH1F*)gFile-> Get("KN_MM_woAll_wN_knEl_Decay"), h1-> Fill(mm_kn);
    }
  }
  else{
    h2_ev-> Fill(9., 3.);
  }
}
