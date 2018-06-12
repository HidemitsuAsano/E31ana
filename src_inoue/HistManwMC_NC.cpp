#include "HistManwMC.h"
#include "MyParam.h"

using namespace std;

void HistManwMC::anaNC(EventHeader *header)
{
  TH1F *h1;
  TH2F *h2;
  TH1F *h1_N  = (TH1F*)rtFile-> Get("N_Reduction");
  TH1F *h1_ER = (TH1F*)rtFile-> Get("EventReduction");

  bool NC_flag_1st = false;
  double NC_time_1st = DBL_MAX;
  double NC_dE_1st = -999;
  HodoscopeLikeHit *hit_1st = 0;

  int nlay=0;
  int nNC=0;
  for( int lay=0; lay<8; lay++ ){
    if( fNC_hit[lay].size()!=0 ){
      if( !fNC_hit[lay].empty() ) nlay++;
      nNC += fNC_hit[lay].size();

      for( int i=0; i<fNC_hit[lay].size(); i++ ){
	//      if( fNC_hit[lay][i]->emean()>NC_dE_1st ){                                                                       
	if( fNC_hit[lay][i]->ctmean()<NC_time_1st ){
	  NC_flag_1st = true;
	  NC_time_1st = fNC_hit[lay][i]-> ctmean();
	  NC_dE_1st = fNC_hit[lay][i]-> emean();
	  hit_1st = fNC_hit[lay][i];
	}

	if( fNC_hit[lay][i]->emean()>NC_THRE ){
	  TVector3 nc_pos;
	  confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNC_hit[lay][i]->seg() , nc_pos);
	  nc_pos.SetY(hit_1st->hitpos());

	  double fl = (nc_pos-fVtxBeam).Mag();
	  double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	  double NC_beta = fl/((fNC_hit[lay][i]->ctmean()-fT0time-beam_tof)*100.*Const);
	
	  if( NC_beta>0.9 ){
	    double calc_tof = fl/(Const*100.);
	    double offset = fNC_hit[lay][i]->ctmean()-fT0time-beam_tof-calc_tof;

	    // cout<<"===== Event Number : "<<header->ev()<<endl;
	    // cout<<Form("T0pos (%lf, %lf, %lf)", fT0pos.X(), fT0pos.Y(), fT0pos.Z())<<endl;
	    // cout<<Form("Vtx (%lf, %lf, %lf)", fVtxBeam.X(), fVtxBeam.Y(), fVtxBeam.Z())<<endl;
	    // cout<<" NC seg"<<fNC_hit[lay][i]->seg()<<"  offset="<<offset<<endl;
	    if( !simReader ){	
	      h1 = (TH1F*)rtFile-> Get(Form("T0NC_offset_NC%d_All", fNC_hit[lay][i]->seg())), h1-> Fill(offset);
	    }
	  }
	}
      }
    }
    if( NC_flag_1st ) break;
  }
  h1 = (TH1F*)rtFile-> Get("nlayNC"), h1-> Fill(nlay);
  h1 = (TH1F*)rtFile-> Get("nNC"), h1-> Fill(nNC);

  if( NC_flag_1st ){
    TVector3 nc_pos;
    confMan-> GetGeomMapManager()-> GetGPos(CID_NC, hit_1st->seg() , nc_pos);
    nc_pos.SetY(hit_1st->hitpos());
    double nc_time = hit_1st->ctmean();
    double fl;

    if( !simReader ) fl = (nc_pos-fVtxBeam).Mag();
    else  fl = (nc_pos-fVtxBeam).Mag()-2.5;

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
    double NC_beta_1st = fl/((nc_time-fT0time-beam_tof)*100.*Const);

    if( fBeamPID==Beam_Kaon ){
      h1 = (TH1F*)rtFile-> Get("NC_overbeta"), h1-> Fill(1./NC_beta_1st);
      h2 = (TH2F*)rtFile-> Get("NC_overbeta_dE"), h2-> Fill(1./NC_beta_1st, hit_1st->emean());
      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){


	h1 = (TH1F*)rtFile-> Get("NC_overbeta_wtar"), h1-> Fill(1./NC_beta_1st);
	h2 = (TH2F*)rtFile-> Get("NC_overbeta_dE_wtar"), h2-> Fill(1./NC_beta_1st, hit_1st->emean());
      }
    }
  }

  //************************//
  //*** for NC S/N study ***//
  //************************//
  double NC_dE_test[30];
  bool NC_test_flag[30];
  double NC_beta_test[30];
  for( int i=0; i<30; i++ ){
    NC_dE_test[i]   = -999.;
    NC_test_flag[i] = false;
    NC_beta_test[i] = 0.;
  }
  for( int lay=0; lay<8; lay++ ){
    for( int i=0; i<fNC_hit[lay].size(); i++ ){
      for( int j=0; j<30; j++ ){
	if( !NC_test_flag[j] && 2*j<fNC_hit[lay][i]->emean() ){
	  NC_dE_test[j] = fNC_hit[lay][i]->emean();
	  TVector3 nc_pos;
	  confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNC_hit[lay][i]->seg() , nc_pos);
	  nc_pos.SetY(fNC_hit[lay][i]->hitpos());
	  double nc_time = fNC_hit[lay][i]->ctmean();
	  double fl;

	  if( !simReader ) fl = (nc_pos-fVtxBeam).Mag();
	  else  fl = (nc_pos-fVtxBeam).Mag()-2.5;

	  double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
	  NC_beta_test[j] = fl/((nc_time-fT0time-beam_tof)*100.*Const);
	  NC_test_flag[j] = true;
	}
      }
    }
  }
  for( int i=0; i<30; i++ ){
    if( fBeamPID==Beam_Kaon ){
      if( NC_test_flag[i] && NC_dE_test[i]<2*i+2. ){
	h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_%d", 2*i+1)), h1-> Fill(1./NC_beta_test[i]);
	if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
	  h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_%d_wtar", 2*i+1)), h1-> Fill(1./NC_beta_test[i]);
	}
      }
    }
  }

  //************************//
  //*** NC Main Analysis ***//
  //************************//
  bool NC_flag = false;
#if 0
  //*** select yangest layer & Max dE
  for( int lay=0; lay<8; lay++ ){
    for( int i=0; i<fNC_hit[lay].size(); i++ ){
      if( fNC_hit[lay][i]->emean()>NC_THRE && fNC_hit[lay][i]->emean()>fNCdE ){
      //      if( fNC_hit[lay][i]->emean()>NC_THRE && fNC_hit[lay][i]->ctmean()<fNCtime ){                                      
	NC_flag =true;
	fNC_eff_hit = fNC_hit[lay][i];
	fNCseg  = fNC_hit[lay][i]->seg();
	fNCtime = fNC_hit[lay][i]->ctmean();
	fNCdE   = fNC_hit[lay][i]->emean();
	confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNC_hit[lay][i]->seg() , fNCpos);

	fNCpos.SetY(fNC_hit[lay][i]->hitpos());
	//    fNCpos.SetY(fNC_hit[lay][i]->hitpos());                                                                           
      }
    }
    if( NC_flag ) break;
  }
#endif
#if 1
  //*** select fastest hit
  for( int lay=0; lay<8; lay++ ){
    for( int i=0; i<fNC_hit[lay].size(); i++ ){
      //      if( fNC_hit[lay][i]->emean()>NC_THRE && fNC_hit[lay][i]->emean()>fNCdE ){
      if( fNC_hit[lay][i]->emean()>NC_THRE && fNC_hit[lay][i]->ctmean()<fNCtime ){                                      
	NC_flag =true;
	fNC_eff_hit = fNC_hit[lay][i];
	fNCseg  = fNC_hit[lay][i]->seg();
	fNCtime = fNC_hit[lay][i]->ctmean();
	fNCdE   = fNC_hit[lay][i]->emean();
	confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNC_hit[lay][i]->seg() , fNCpos);

	fNCpos.SetY(fNC_hit[lay][i]->hitpos());
	//    fNCpos.SetY(fNC_hit[lay][i]->hitpos());                                                                           
      }
    }
  }
#endif
#if 0
  //  std::cout<<">>>>> NC Clustering Analysis START "<<std::endl;
  std::vector<HodoscopeLikeHit*> tmpNC_hit[7];
  for( int lay=0; lay<7; lay++ ){
    //    std::cout<<">>> NC lay : "<<lay+1<<"  nHit : "<<fNC_hit[lay].size()<<std::endl;
    for( int i=0; i<fNC_hit[lay].size(); i++ ){
      tmpNC_hit[lay].push_back(fNC_hit[lay][i]);
      //      std::cout<<"> seg:"<<fNC_hit[lay][i]->seg()<<" tiem:"<<fNC_hit[lay][i]->ctmean()<<" dE:"<<fNC_hit[lay][i]->emean()<<std::endl;
    }
  }

  std::vector<std::vector<HodoscopeLikeHit*> > NC_clusters;
  while(true){
    bool all_empty=true;
    for( int lay=0; lay<7; lay++ ){
      if( !tmpNC_hit[lay].empty() ) all_empty=false;
    }
    if( all_empty ) break;

    std::vector<HodoscopeLikeHit*> cluster;

    //*** push_back 1st hit in cluster ***//
    for( int lay=0; lay<7; lay++ ){
      if( !tmpNC_hit[lay].empty() ){
	cluster.push_back(tmpNC_hit[lay][0]);
	tmpNC_hit[lay].erase(tmpNC_hit[lay].begin());
	break;
      }
    }
    while( true ){
      bool add=false;
      for( int lay=0; lay<7; lay++ ){
	for( int i=0; i<tmpNC_hit[lay].size(); i++ ){
	  int seg=tmpNC_hit[lay][i]->seg()%16;
	  for( int j=0; j<cluster.size(); j++ ){
	    int lay2 = (cluster[j]-> seg()-1)/16;
	    int seg2 = cluster[j]->seg()%16;
	    if( lay==lay2 ){
	      if( seg-1==seg2 || seg+1==seg2 ){
		add=true;
		cluster.push_back(tmpNC_hit[lay][i]);
		std::vector<HodoscopeLikeHit*>::iterator ite=tmpNC_hit[lay].begin();
		for( int k=0; k<i; k++ ) ++ite;
		tmpNC_hit[lay].erase(ite);
		break;
	      }
	    }
	    if( lay+1==lay2 || lay-1==lay2 ){
	      if( seg==seg2 || seg-1==seg2 || seg+1==seg2 ){
		add=true;
		cluster.push_back(tmpNC_hit[lay][i]);
		std::vector<HodoscopeLikeHit*>::iterator ite=tmpNC_hit[lay].begin();
		for( int k=0; k<i; k++ ) ++ite;
		tmpNC_hit[lay].erase(ite);
		break;
	      }
	    }
	  }
	  if( add ) break;
	}
	if( add ) break;
      }
      if( !add ) break;
    }

    NC_clusters.push_back(cluster);
  }
  //  std::cout<<"Num of NC cluster : "<<NC_clusters.size()<<std::endl;
  h1 = (TH1F*)rtFile-> Get("nNCCluster"), h1-> Fill(NC_clusters.size());
  for( int i=0; i<NC_clusters.size(); i++ ){
    h1 = (TH1F*)rtFile-> Get("NCClusterSize"), h1-> Fill(NC_clusters[i].size());
    h2 = (TH2F*)rtFile-> Get("nNCCluster_Size"), h2-> Fill(NC_clusters.size(),  NC_clusters[i].size());
  }

  HodoscopeLikeHit *fastest_hit[NC_clusters.size()];
  double fastest_time[NC_clusters.size()];
  bool NC_eff_flag[NC_clusters.size()];
  for( int i=0; i<NC_clusters.size(); i++ ){
    fastest_hit[i]=0;
    fastest_time[i]=DBL_MAX;
    NC_eff_flag[i]=false;
  }

  for( int i=0; i<NC_clusters.size(); i++ ){
    for( int j=0; j<NC_clusters[i].size(); j++ ){
      if( NC_clusters[i][j]->emean()>NC_THRE ){
	NC_eff_flag[i]=true;
	NC_flag=true;
      }
      if( NC_clusters[i][j]->ctmean()<fastest_time[i] ){
	fastest_time[i] = NC_clusters[i][j]->ctmean();
	fastest_hit[i]  = NC_clusters[i][j];
      }
    }
  }

  int fastest_index=-1;
  double fastest_timming = DBL_MAX;
  for( int i=0; i<NC_clusters.size(); i++ ){
    if( NC_eff_flag[i] ){
      if( fastest_timming>fastest_time[i] ){
	fastest_timming = fastest_time[i];
	fNC_eff_hit = fastest_hit[i];
      }
    }
  }
  if( NC_flag ){
    fNCseg  = fNC_eff_hit->seg();
    fNCtime = fNC_eff_hit->ctmean();
    fNCdE   = fNC_eff_hit->emean();
    confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNCseg , fNCpos);
    
    fNCpos.SetY(fNC_eff_hit->hitpos());
  }

  // std::string hit_str[7];
  // for( int i=0; i<7; i++) hit_str[i]="| | | | | | | | | | | | | | | | |";
  // for( int i=0; i<NC_clusters.size(); i++ ){
  //   for( int j=0; j<NC_clusters[i].size(); j++ ){
  //     int lay=(NC_clusters[i][j]->seg()-1)/16;
  //     int seg=NC_clusters[i][j]->seg()%16;
  //     if( i==0 )      hit_str[lay][2*seg+1]='1';
  //     else if( i==1 ) hit_str[lay][2*seg+1]='2';
  //     else if( i==2 ) hit_str[lay][2*seg+1]='3';
  //     else if( i==3 ) hit_str[lay][2*seg+1]='4';
  //     else if( i==4 ) hit_str[lay][2*seg+1]='5';
  //     else if( i==5 ) hit_str[lay][2*seg+1]='6';
  //     else            hit_str[lay][2*seg+1]='0';
  //   }
  // }
  // for( int i=0; i<7; i++) std::cout<<hit_str[i]<<std::endl;
  // std::string str;
  // std::cin>>str;
#endif

  if( NC_flag ){
    double fl;
    if( !simReader ) fl = (fNCpos-fVtxBeam).Mag();
    else fl = (fNCpos-fVtxBeam).Mag()-2.5;

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
    fNCbeta = fl/((fNCtime-fT0time-beam_tof)*100.*Const);

    int NClay = (fNCseg-1)/16+1;
    int NCseg2 = fNCseg%16;
    if( NCseg2==0 ) NCseg2=16;

    if( fBeamPID==Beam_Kaon ){
      h2 = (TH2F*)rtFile-> Get("overbeta_NCseg"), h2-> Fill(1./fNCbeta, NCseg2);
      h2 = (TH2F*)rtFile-> Get("overbeta_NClay"), h2-> Fill(1./fNCbeta, NClay);

      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee"), h1-> Fill(1./fNCbeta);
      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma"), h1-> Fill(1./fNCbeta);
      h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_8MeVee_gamma_seg%d", fNCseg)), h1-> Fill(1./fNCbeta);
      if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
#if CALIB
	if( !simReader ){
	  h1 = (TH1F*)rtFile-> Get(Form("NC_tsub_%d", fNCseg)), h1-> Fill(fNC_eff_hit-> tsub());
	  h1 = (TH1F*)rtFile-> Get(Form("NC_hitpos_%d", fNCseg)), h1-> Fill(fNC_eff_hit-> hitpos());
	}
#endif
	h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_wtar"), h1-> Fill(1./fNCbeta);
	h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma_wtar"), h1-> Fill(1./fNCbeta);
	h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_8MeVee_gamma_wtar_seg%d", fNCseg)), h1-> Fill(1./fNCbeta);
      }
    } 
    if( header && fBeamPID==Beam_Kaon ){
      //    std::cout<<"  NC hit"<<std::endl;                                                                                   
      h1_ER-> Fill(8);
      if( header->IsTrig(Trig_Neutral) ){
	//    std::cout<<"  Neutral trigger"<<std::endl;                                                                        
	h1_N-> Fill(8);
	h1_ER-> Fill(9);
      }
      else{
	h1_N-> Fill(8);
      }
    }

    if( fNCbeta>0.9 ){
      fFPID=F_Gamma;
      if( header ){
	if( header->IsTrig(Trig_Neutral) && fBeamPID==Beam_Kaon ) h1_N-> Fill(9);
      }

      double calc_tof = fl/(Const*100.);
      double offset = fNC_eff_hit->ctmean()-fT0time-beam_tof-calc_tof;
      int T0seg=fT0_hit[0]-> seg();
      double T0eu=fT0_hit[0]->eu(), T0ed=fT0_hit[0]->ed();
      double T0tu=fT0_hit[0]->tu()-fNC_eff_hit->ctmean()+beam_tof+calc_tof;
      double T0td=fT0_hit[0]->td()-fNC_eff_hit->ctmean()+beam_tof+calc_tof;
      double T0ctu=fT0_hit[0]->ctu()-fNC_eff_hit->ctmean()+beam_tof+calc_tof;
      double T0ctd=fT0_hit[0]->ctd()-fNC_eff_hit->ctmean()+beam_tof+calc_tof;
      int NCseg=fNC_eff_hit-> seg();
      double NCeu=fNC_eff_hit->eu(), NCed=fNC_eff_hit->ed();
      double NCtu=fNC_eff_hit->tu()-fT0time-beam_tof-calc_tof;
      double NCtd=fNC_eff_hit->td()-fT0time-beam_tof-calc_tof;
      double NCctu=fNC_eff_hit->ctu()-fT0time-beam_tof-calc_tof;
      double NCctd=fNC_eff_hit->ctd()-fT0time-beam_tof-calc_tof;
    
      if( !simReader ){
	//	h2 = (TH2F*)rtFile-> Get("NC_eu_ctm"), h2-> Fill(fNC_eff_hit->eu(), offset);
	//	h2 = (TH2F*)rtFile-> Get("NC_ed_ctm"), h2-> Fill(fNC_eff_hit->ed(), offset);
	//	std::cout<<"offset"<<offset<<" NC seg"<<NCseg<<std::endl;

	// TNtuple *tup = (TNtuple*)rtFile-> Get(Form("NC%d_slewing_info", NCseg));
	// tup-> Fill(NCeu, NCed, NCtu, NCtd, NCctu, NCctd, offset);

	// tup = (TNtuple*)rtFile-> Get(Form("T0%d_slewing_info", T0seg));
	// tup-> Fill(T0eu, T0ed, T0tu, T0td, T0ctu, T0ctd, -offset);

	// tup = (TNtuple*)rtFile-> Get(Form("NC%d_slewing_info_T0%d", NCseg, T0seg));
	// tup-> Fill(NCeu, NCed, NCtu, NCtd, NCctu, NCctd, offset);

	// tup = (TNtuple*)rtFile-> Get(Form("T0%d_slewing_info_NC%d", T0seg, NCseg));
	// tup-> Fill(T0eu, T0ed, T0tu, T0td, T0ctu, T0ctd, -offset);

	h1 = (TH1F*)rtFile-> Get(Form("T0NC_offset_T0%d", fT0_hit[0]->seg())), h1-> Fill(offset);
	h1 = (TH1F*)rtFile-> Get(Form("T0NC_offset_NC%d", fNC_eff_hit->seg())), h1-> Fill(offset);	
	h1 = (TH1F*)rtFile-> Get(Form("T0NC_offset_T0%d_NC%d", fT0_hit[0]->seg(), fNC_eff_hit->seg())), h1-> Fill(offset);
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
	  //            std::cout<<"  neutron"<<std::endl;                                                                      
	  h1_N-> Fill(10);
	  h1_ER->Fill(10);
	}
      }
      else{
	h1_ER->Fill(10);
      }
    }
  }
}

void HistManwMC::anaNC_MC()
{
  TH1F *h1;
  DetectorData *detData = simReader-> getDetectorData();

  DetectorHit *hit_1st=0;
  double first_time=DBL_MAX;
  int nNC_hit = 0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit-> detectorID()==CID_NC ){
      if( hit-> pdg()==2112 || hit-> pdg()==22 ){
        if( hit-> tdc()<first_time ){
          nNC_hit++;
          hit_1st = hit;
          first_time = hit-> tdc();
        }
      }
    }
  }

  if( !hit_1st ){
    fFPID=F_Other;
    return;
  }

  fNCseg = hit_1st-> channelID()+1;
  fNC_eff_hit = 0;
  HodoscopeLikeHit *hodo_hit = 0;
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->seg()==fNCseg ) fNC_eff_hit=blMan->NC(i);
  }
  if( !fNC_eff_hit ){
    std::cout<<"  !!! HistManwwMC::anaMC_NC not find HodoscopeLikeHit !!!"<<std::endl;
    return;
  }

  fNCdE = fNC_eff_hit-> emean();
  fNCtime = fNC_eff_hit->ctmean();
  confMan-> GetGeomMapManager()-> GetGPos(CID_NC, fNCseg , fNCpos);
  fNCpos.SetY(fNC_eff_hit->hitpos());

  //  simReader->trace(CID_NC, fNCseg, true);

  double fl = (fNCpos-fVtxBeam).Mag()-2.5;
  TVector3 vtxMC = simReader->getMCTrack(1)->vertex();
  TVector3 hitposMC = hit_1st-> pos();
  double flMC = 0.1*(hitposMC-vtxMC).Mag();
  h1 = (TH1F*)rtFile-> Get("T0NC_FL_diff_MC"), h1-> Fill(flMC-fl);

  //  double fl = (fNCpos-fVtxBeam).Mag();
  double beam_out, beam_tof;
  ELossTools::CalcElossBeamTGeo(fT0pos, fVtxBeam, fD5mom, parMass[fBeamPID], beam_out, beam_tof);
  fNCbeta = fl/((fNCtime-fT0time-beam_tof)*100.*Const);

  if( fBeamPID==Beam_Kaon ){
    h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee"), h1-> Fill(1./fNCbeta);
    if( hit_1st->pdg()==22 ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma"), h1-> Fill(1./fNCbeta);
    h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_8MeVee_gamma_seg%d", fNCseg)), h1-> Fill(1./fNCbeta);
    if( GeomTools::GetID(fVtxBeam)==CID_Fiducial ){
      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_wtar"), h1-> Fill(1./fNCbeta);
      if( hit_1st->pdg()==22 ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma_wtar"), h1-> Fill(1./fNCbeta);
      h1 = (TH1F*)rtFile-> Get(Form("NC_overbeta_8MeVee_gamma_wtar_seg%d", fNCseg)), h1-> Fill(1./fNCbeta);
      h1 = (TH1F*)rtFile-> Get("NC_overbeta_8MeVee_gamma_wtar"), h1-> Fill(1./fNCbeta);
    }
  }

  if( hit_1st->pdg()==22 ){
    h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_MC"), h1-> Fill(1./fNCbeta);
    if( simReader->isPiDecay(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_from_pi_MC"), h1-> Fill(1./fNCbeta);
    else if( simReader->isSigmaDecay(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_from_S_MC"), h1-> Fill(1./fNCbeta);
    else if( simReader->isFromE(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_from_e_MC"), h1-> Fill(1./fNCbeta);
    else if( simReader->isInitial(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_init_MC"), h1-> Fill(1./fNCbeta);
    else h1 = (TH1F*)rtFile-> Get("NC_overbeta_gamma_acci_MC"), h1-> Fill(1./fNCbeta);
    fFPID=F_Gamma;
  }
  if( hit_1st->pdg()==2112 ){
    h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_MC"), h1-> Fill(1./fNCbeta);
    if( simReader->isInitial(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_init_MC"), h1-> Fill(1./fNCbeta);
    else if( simReader->isSigmaDecay(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_from_S_MC"), h1-> Fill(1./fNCbeta);
    else if( simReader->isLambdaDecay(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_from_L_MC"), h1-> Fill(1./fNCbeta);
    else if( simReader->isL1520Decay(CID_NC, fNCseg) ) h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_from_L1520_MC"), h1-> Fill(1./fNCbeta);
    else  h1 = (TH1F*)rtFile-> Get("NC_overbeta_n_acci_MC"), h1-> Fill(1./fNCbeta);

    fFPID=F_Neutron;
    double NCmom = nMass*fNCbeta/sqrt(1-fNCbeta*fNCbeta);
    TVector3 n_mom = fNCpos-fVtxBeam;
    n_mom.SetMag(NCmom);
    Track *n_track = simReader->trace(CID_NC, fNCseg);
    double n_mom_MC = 0.001*n_track->momentum().Mag();
    h1 = (TH1F*)rtFile-> Get("n_mom_diff_MC"), h1-> Fill(n_mom_MC-n_mom.Mag());

    fFLmom.SetVectM(n_mom, nMass);
  }
}
