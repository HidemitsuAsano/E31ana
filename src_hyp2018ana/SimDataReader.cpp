#include "SimDataReader.h"

using namespace std;

SimDataReader::SimDataReader() : runHeader(0), evHeader(0), reacData(0), detData(0), mcData(0)
{
}

bool SimDataReader::isSameHit(HodoscopeLikeHit *hit1, HodoscopeLikeHit *hit2)
{
  int cid1=hit1->cid(), cid2=hit2->cid(),  seg1=hit1->seg(), seg2=hit2->seg();
  DetectorHit *d_hit1 = getHit(cid1, seg1);
  DetectorHit *d_hit2 = getHit(cid2, seg2);

  int track_id1=d_hit1->trackID();
  int track_id2=d_hit2->trackID();

  int parent_id1=d_hit1->parentID();
  int parent_id2=d_hit2->parentID();

  bool result=false;
  if( track_id1==track_id2 ) result=true;
  if( track_id1==parent_id2 ) result=true;
  if( track_id2==parent_id1 ) result=true;

  // if( result ){
  //   cout<<"SimDataReade::isSameHit"<<endl;
  //   cout<<" 1> tID : "<<track_id1<<"  pID : "<<parent_id1<<endl;
  //   cout<<" 2> tID : "<<track_id2<<"  pID : "<<parent_id2<<endl;
  //   cout<<" result : "<<boolalpha<<result<<endl;
  // }

  return result;
}

void SimDataReader::printProcess(const int &trackID)
{
  Track *track = getMCTrack(trackID);
  int parentID = track->parentTrackID();
  std::string process_name=getProcessName(track->processID());
  vector<Track*> tracks;
  for( int i=0; i<mcData->trackSize(); i++ ){ 
    Track *track2 = mcData-> track(i);
    if( track2->parentTrackID()==parentID ){ 
      tracks.push_back(track2);
    }
  }

  Track *p_track = getMCTrack(parentID);
  cout<<"===== SimDataReader::printProcess ====="<<endl;
  cout<<"> process name : "<<process_name<<"  parent par : "<<parName(p_track->pdgID())<<endl;
  for( int i=0; i<tracks.size(); i++ ){
    cout<<">                         par : "<<parName(tracks[i]->pdgID())<<endl;
  }
}

bool SimDataReader::printCDCTrack(CDSHitMan *cdsMan, CDSTrack* track)
{
  vector<DetectorHit*> CDHhits;
  for( int i=0; i<track->nCDHHit(); i++ ){
    CDHhits.push_back(getHit(CID_CDH, track->CDHHit(cdsMan, i)-> seg()));
  }

  bool all_same=true;
  int pdg1=-999999;
  if( CDHhits.size()==0 ) all_same=false;
  else{
    pdg1=CDHhits[0]->pdg();
    for( int i=1; i<CDHhits.size(); i++ ){
      if( pdg1!=CDHhits[i]->pdg() ) all_same=false;
    }
  }
  if( track->PID()==CDS_PiMinus && pdg1!=-211 ) all_same=false;
  if( track->PID()==CDS_PiPlus  && pdg1!=211  ) all_same=false;
  if( track->PID()==CDS_Kaon    && pdg1!=-321 ) all_same=false;
  if( track->PID()==CDS_Proton  && pdg1!=2122 ) all_same=false;

  cout<<"> CDH hit seg : ";
  for( int i=0; i<CDHhits.size(); i++ ) cout<<CDHhits[i]->channelID()+1<<" ";
  cout<<endl;
  cout<<"> CDH hit par : ";
  for( int i=0; i<CDHhits.size(); i++ ) cout<<parName(CDHhits[i]->pdg())<<" ";
  cout<<endl;
  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    string str=Form("> lay%2d ", lay);
    string str2=" ";
    for( int i=0; i<track->nTrackHit(lay); i++ ){
      CDCHit *hit = track->TrackHit(cdsMan, lay, i);
      DetectorHit *d_hit = getHit(CID_CDC, lay, hit->wire());
      if( pdg1!= d_hit->pdg() ) all_same=false;

      str += Form("%d ", hit->wire());
      str2 += parName(d_hit->pdg())+" ";
    }
    cout<<str<<" "<<str2<<endl;
  }

  if( !all_same ){

  }

  if( CDHhits.size()!=1 ) return false;
  if( all_same!=1 ) return false;

  return true;
}

int SimDataReader::checkCDS(CDSTrack *track, CDSTrackingMan *cdstrackingMan, CDSHitMan *cdsMan, const double &mass2, const TVector3 &vtxCDS, const TVector3 &vtxBeam, bool print)
{
  vector<HodoscopeLikeHit*> CDHHits;
  for( int i=0; i<track->nCDHHit(); i++ ){ CDHHits.push_back(track->CDHHit(cdsMan, i)); }
  vector<DetectorHit*> DetHits;
  for( int i=0; i<track->nCDHHit(); i++ ){ DetHits.push_back(getHit(CID_CDH, CDHHits[i]->seg())); }

  if( CDHHits.size()==0 ){
    cout<<"  !!! CDSTrack doesn't have CDHHit !!!"<<endl;
    return -1;
  }
  int hit_pdg=DetHits[0]->pdg();
  int pid=track->PID();
  std::string cds_par;
  if( pid==CDS_PiPlus   )      cds_par="pi+     ";
  else if( pid==CDS_PiMinus  ) cds_par="pi-     ";
  else if( pid==CDS_Kaon     ) cds_par="kaon-   ";
  else if( pid==CDS_Proton   ) cds_par="proton  ";
  else if( pid==CDS_Deuteron ) cds_par="deuteron";
  else                         cds_par="unknown ";

  vector<Track*> tracks;
  tracks.push_back(getMCTrack(DetHits[0]->trackID()));
  while( tracks[tracks.size()-1]->parentTrackID()>0 ){ tracks.push_back(getMCTrack(tracks[tracks.size()-1]->parentTrackID())); }

  bool pid_match=true;
  if( pid==CDS_PiPlus ){ if( hit_pdg!= 211 ) pid_match=false; }
  else if( pid==CDS_PiMinus ){ if( hit_pdg!=-211 ) pid_match=false; }
  else if( pid==CDS_Kaon    ){ if( hit_pdg!=-321 ) pid_match=false; }
  else if( pid==CDS_Proton  ){ if( hit_pdg!=2212 ) pid_match=false; }
  else if( pid==CDS_Deuteron ){ if( hit_pdg!=1000010020 ) pid_match=false; }
  else{ pid_match=false; }

  bool CDH2hit=false;
  bool CDC2tra=false;
  int result=9999;
  if( hit_pdg==11 ){
    for( int i=0; i<tracks.size()-1; i++ ){ 
      if( tracks[i+1]->pdgID()==22 ){
	result=2;
	break;
      }
    }
  }
  else if( hit_pdg==-11 ){
    for( int i=0; i<tracks.size()-1; i++ ){ 
      if( tracks[i+1]->pdgID()==22 ){
	result=3;
	break;
      }
    }
  }
  else if( hit_pdg==13 ){
    for( int i=0; i<tracks.size()-1; i++ ){ 
      if( tracks[i]->processID()==1 && tracks[i+1]->pdgID()==-211 ){
	result=10;
	break;
      }
      else if( tracks[i]->processID()==1 && tracks[i+1]->pdgID()==-321 ){
	result=11;
	break;
      }
    }
  }
  else if( hit_pdg==-13 ){
    for( int i=0; i<tracks.size()-1; i++ ){ 
      if( tracks[i]->processID()==1 && tracks[i+1]->pdgID()==211 ){
	result=15;
	break;
      }
      else if( tracks[i]->processID()==1 && tracks[i+1]->pdgID()==321 ){
	result=16;
	break;
      }
    }
  }
  else if( hit_pdg==-211 ){
    for( int i=0; i<tracks.size(); i++ ){ 
      if( i<tracks.size()-1 ){
	if( tracks[i]->processID()==1 && tracks[i+1]->pdgID()==-321 ){
	  result=20;
	  break;
	}
	else if( tracks[i]->processID()==100 && tracks[i+1]->pdgID()==-211 ){
	  result=21;
	  break;
	}
	else if( tracks[i+1]->pdgID()==-321 ){
	  result=22;
	  break;
	}
      }
      if( tracks[i]->pdgID()==-211 && (tracks[i]->processID()==1 || tracks[i]->processID()==0) ){
	result=0;
	break;
      }
    }
  }
  else if( hit_pdg==211 ){
    for( int i=0; i<tracks.size(); i++ ){ 
      if( i<tracks.size()-1 ){
	if( tracks[i]->processID()==1 && tracks[i+1]->pdgID()==-321 ){
	  result=30;
	  break;
	}
	else if( tracks[i]->processID()==100 && tracks[i+1]->pdgID()==211 ){
	  result=31;
	  break;
	}
	else if( tracks[i+1]->pdgID()==-321 ){
	  result=32;
	  break;
	}
      }
      if( tracks[i]->pdgID()==211 && (tracks[i]->processID()==1 || tracks[i]->processID()==0) ){
	result=0;
	break;
      }
    }
  }
  else if( hit_pdg==-321 ){ 
    for( int i=0; i<tracks.size(); i++ ){
      if( tracks[i]->processID()==0 || tracks[i]->processID()==1 ){
	result=0;
	break;
      }
    }
  }
  else if( hit_pdg==2212 ){
    for( int i=0; i<tracks.size(); i++ ){
      if( i<tracks.size()-1 ){
	if( tracks[i+1]->pdgID()==2112 ){
	  result=50;
	  break;
	}
	else if( tracks[i+1]->pdgID()==-321 ){
	  result=51;
	  break;
	}
	else if( tracks[i+1]->pdgID()==-211 ){
	  result=52;
	  break;
	}
	else if( tracks[i+1]->pdgID()==211 ){
	  result=53;
	  break;
	}
      }
      if( tracks[i]->pdgID()==2212 && (tracks[i]->processID()==0 || tracks[i]->processID()==1) ){
	result=0;
	break;
      }
    }
  }
  if( DetHits.size()>1 ){ result=100, CDH2hit=true; }
  for( int i=0; i<cdstrackingMan->nGoodTrack(); i++ ){
    CDSTrack *track2 = cdstrackingMan->GoodTrack(i);
    if( track==track2 ){ continue; }

    for( int j=0; j<track2->nCDHHit(); j++ ){
      for( int k=0; k<track->nCDHHit(); k++ ){
	int seg1= track->CDHHit(cdsMan, k)->seg();
	int seg2= track2->CDHHit(cdsMan, j)->seg();
	if( seg1==seg2 ) CDC2tra=true;
      }
    }
  }
  if( CDC2tra ) result=200;

  // if( CDC2tra ){
  //   cout<<"      CDC 2track"<<endl;
  //   cout<<"===== print All CDC Track ====="<<endl;
  //   for( int i=0; i<cdstrackingMan->nGoodTrack(); i++ ){
  //     printCDCTrack(cdsMan, cdstrackingMan->GoodTrack(i));
  //   }
  //     string input;
  //     cin>>input;
  //     if( input=="q" ) exit(0);
  // }


  if( print ){
  //  if( (hit_pdg==-211 && track->Momentum()<0 && mass2>0.3) || CDC2tra ){
    cout<<"===== SimDataReader::checkCDS ====="<<endl;
    cout<<"> nCDH : "<<track->nCDHHit()<<"  seg : ";
    for( int i=0; i<CDHHits.size(); i++ ){ cout<<CDHHits[i]->seg()<<" "; }
    cout<<endl;

    cout<<"> par cds : "<<cds_par<<"  MC : "<<parName(hit_pdg)<<"  mass2 : "<<mass2<<endl;
    for( int i=1; i<DetHits.size(); i++ ){
      cout<<"                          : "<<parName(DetHits[i]->pdg())<<endl;
    }
    cout<<"> mom : "<<track->Momentum()<<"  MC : "<<0.001*tracks[0]->momentum().Mag()<<endl;

    string header=">>";
    for( int i=0; i<tracks.size(); i++ ){
      cout<<header<<" "<<parName(tracks[i]->pdgID())<<"  reaction : "<<getProcessName(tracks[i]->processID())<<"  trackID : "<<tracks[i]->trackID()<<" parentID : "<<tracks[i]->parentTrackID()<<endl;
      if( getProcessName(tracks[i]->processID()).find("elastic")!=string::npos ){
	printProcess(tracks[i]->trackID());
      }
      if( tracks[i]->processID()==1 && (tracks[i]->pdgID()==-321 || tracks[i]->pdgID()==321) ){
	printProcess(tracks[i]->trackID());
      }
      header += ">";
    }

    printCDCTrack(cdsMan, track);

    if( !pid_match ){ cout<<">  hit par not match"<<endl; }
    if( result==2 ) cout<<"> gamma -> e- conv"<<endl;
    else if( result==3 ) cout<<"> gamma -> e+ conv"<<endl;
    else if( result==10 ) cout<<"> pi- -> mu- Decay"<<endl;
    else if( result==11 ) cout<<"> K- -> mu- Decay"<<endl;
    else if( result==15 ) cout<<"> pi+ -> mu+ Decay"<<endl;
    else if( result==16 ) cout<<"> K+ -> mu+ Decay"<<endl;
    else if( result==20 ) cout<<"> K- -> pi- Decay"<<endl;
    else if( result==21 ) cout<<"> pi-Inelastic   pi- hit"<<endl;
    else if( result==22 ) cout<<"> kaon-Inelastic   pi- hit"<<endl;
    else if( result==30 ) cout<<"> K+ -> pi+ Decay"<<endl;
    else if( result==31 ) cout<<"> pi+Inelastic   pi+ hit"<<endl;
    else if( result==32 ) cout<<"> kaon-Inelastic   pi+ hit"<<endl;
    else if( result==50 ) cout<<"> neutron -> proton hit"<<endl;
    else if( result==51 ) cout<<"> kaon- -> proton hit"<<endl;
    else if( result==52 ) cout<<"> pi- -> proton hit"<<endl;
    else if( result==53 ) cout<<"> pi+ -> proton hit"<<endl;
    else if( result==100 ) cout<<"> pi+ -> proton hit"<<endl;

    cout<<"> result : "<<result<<endl;

    if( print ){
      string input;
      cin>>input;
      if( input=="q" ) exit(0);
    }
  }

  return result;
}

void SimDataReader::setTree2(TTree *tree)
{
  if( tree-> FindBranch("RunHeaderMC") ) tree-> SetBranchAddress("RunHeaderMC", &runHeader);
}

void SimDataReader::setTree(TTree *tree)
{
  if( tree-> FindBranch("EventHeaderMC") ) tree-> SetBranchAddress("EventHeaderMC", &evHeader);
  if( tree-> FindBranch("ReactionData") ) tree-> SetBranchAddress("ReactionData", &reacData);
  if( tree-> FindBranch("DetectorData") ) tree-> SetBranchAddress("DetectorData", &detData);
  if( tree-> FindBranch("MCData") ) tree-> SetBranchAddress("MCData", &mcData);
}

void SimDataReader::checkStatus()
{
  std::cout<<"===== SimDataReader::checkStatus ====="<<std::endl;
  std::cout<<"> RunHeaderMC       : "<<runHeader<<std::endl;
  std::cout<<"> EventHeaderMC     : "<<evHeader<<std::endl;
  std::cout<<"> ReactionData      : "<<reacData<<std::endl;
  std::cout<<"> DetectorData      : "<<detData<<std::endl;
  std::cout<<"> MCData            : "<<mcData<<std::endl;
  std::cout<<"======================================"<<std::endl;
}

void SimDataReader::printSize()
{
  std::cout<<"===== SimDataReader::printSize ====="<<std::endl;
  if( evHeader ) std::cout<<"> EventID : "<<evHeader->eventID()<<std::endl;
  else std::cout<<"> No EventHeaderMC"<<std::endl;
  if( reacData ){
    std::cout<<"> ReactionID : "<<reacData-> ReactionID()<<std::endl;
    std::cout<<"> Init : "<<reacData->InitPDGSize()<<"  Particle : "<<reacData->PDGSize()<<std::endl;
  }
  else std::cout<<"> No EventHeaderMC"<<std::endl;
  std::cout<<"===================================="<<std::endl;
}

Track* SimDataReader::getMCTrack(const int &trackID)
{
  for( int i=0; i<mcData->trackSize(); i++ ){
    if( mcData->track(i)->trackID()==trackID ) return mcData->track(i);
  }
  std::cout<<"  !!! SimDataReader::getMCTrack("<<trackID<<") not found !!!"<<std::endl;
  return 0;
}

DetectorHit* SimDataReader::getHit(const int &cid, const int &seg)
{
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==cid && hit->channelID()+1==seg ) return detData-> detectorHit(i);
  }
  std::cout<<"  !!! SimDataReader::getHit("<<cid<<", "<<seg<<") not found !!!"<<std::endl;
  return 0;
}

vector<DetectorHit*> SimDataReader::getHits(const int &cid, const int &seg)
{
  vector<DetectorHit*> hits;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==cid && hit->channelID()+1==seg ) hits.push_back(hit);
  }
  return hits;
}


DetectorHit* SimDataReader::getHit(const int &cid, const int &lay, const int &wire)
{
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==cid && hit->channelID()+1==wire && hit->layerID()+1==lay ) return detData-> detectorHit(i);
  }
  std::cout<<"  !!! SimDataReader::getHit("<<cid<<", "<<lay<<", "<<wire<<") not found !!!"<<std::endl;
  return 0;
}


Track* SimDataReader::trace(const int &cid, const int &seg, const bool &print)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return trace(hit, print);
}

Track* SimDataReader::trace(CDSTrack *cds_track, CDSHitMan *cdsMan, const bool &print)
{
  int CDHseg;
  double CDHtime;
  if( !cds_track-> GetCDHHit(cdsMan, CDHseg, CDHtime) ){
    std::cout<<"  !!! SimDataReader::trace not find CDHHit !!!"<<std::endl;
    return 0;
  }
  DetectorHit *hit = getHit(CID_CDH, CDHseg);
  if( !hit ) return 0;

  return trace(hit, print);
}

Track* SimDataReader::trace(DetectorHit *hit, const bool &print)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();
  if( print ){
    TVector3 mom = hit->momentum();
    std::cout<<"===== SimDataReader::trace CID:"<<hit->detectorID()<<" seg:"<<hit->channelID()+1<<" ====="<<std::endl;
    std::cout<<"> Hit Particle : "<<parName(hit->pdg())<<"   trackID : "<<trackID<<" parentID : "<<parentID<<std::endl;
    std::cout<<">   mom : "<<0.001*mom.Mag()<<"[GeV/c]  cos(ang) : "<<mom.CosTheta()<<std::endl;
  }

  while( parentID>0 ){
    Track *track = getMCTrack(trackID);

    if( print ){
      TVector3 mom = track->momentum();
      std::cout<<">>> depth : "<<depth<<"  Create Reaction : "<<getProcessName(track->processID())<<std::endl;
      std::cout<<"> Particle : "<<parName(track->pdgID())<<"  trackID : "<<track->trackID()<<" trackID : "<<track->parentTrackID()<<std::endl;
      std::cout<<">   mom : "<<0.001*mom.Mag()<<"[GeV/c]  cos(ang) : "<<mom.CosTheta()<<std::endl;
    }
    if( track->pdgID()!=originPDG ) break;

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  Track *track = getMCTrack(trackID);
  if( print ){
    TVector3 mom=track->momentum();

    std::cout<<">>> depth : "<<depth<<std::endl;
    std::cout<<"> mom diff : "<<mom.Mag()-originMom.Mag()<<" [MeV/c]  hit mom : "<<originMom.Mag()<<" [GeV/c]   init mom : "<<mom.Mag()<<" [GeV/c]"<<std::endl;
    std::cout<<"==================================================="<<std::endl;

    std::string str;
    std::cin>>str;
  }

  return track;
}

Track* SimDataReader::trace(Track *in_track, const bool &print)
{
  const int originPDG = in_track->pdgID();
  int depth=1;
  int parentID = in_track->parentTrackID();
  int trackID = in_track->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);

    if( print ){
      TVector3 mom = track->momentum();
      std::cout<<">>> depth : "<<depth<<std::endl;
      std::cout<<"> Particle : "<<parName(track->pdgID())<<"  trackID : "<<track->trackID()<<" trackID : "<<track->parentTrackID()<<std::endl;
      std::cout<<">   mom : "<<0.001*mom.Mag()<<"[GeV/c]  cos(ang) : "<<mom.CosTheta()<<std::endl;
    }
    if( track->pdgID()!=originPDG ) break;

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  Track *track = getMCTrack(trackID);
  if( print ){
    TVector3 mom=track->momentum();

    std::cout<<">>> depth : "<<depth<<std::endl;
    std::cout<<"==================================================="<<std::endl;
    std::string str;
    std::cin>>str;
  }

  return track;
}

Track* SimDataReader::traceAll(Track *in_track, const bool &print)
{
  const int originPDG = in_track->pdgID();
  int depth=1;
  int parentID = in_track->parentTrackID();
  int trackID = in_track->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);

    if( print ){
      TVector3 mom = track->momentum();
      std::cout<<">>> depth : "<<depth<<std::endl;
      std::cout<<"> Particle : "<<parName(track->pdgID())<<"  trackID : "<<track->trackID()<<" trackID : "<<track->parentTrackID()<<std::endl;
      std::cout<<">   mom : "<<0.001*mom.Mag()<<"[GeV/c]  cos(ang) : "<<mom.CosTheta()<<std::endl;
    }

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  Track *track = getMCTrack(trackID);
  if( print ){
    TVector3 mom=track->momentum();

    std::cout<<">>> depth : "<<depth<<std::endl;
    std::cout<<"==================================================="<<std::endl;
    std::string str;
    std::cin>>str;
  }

  return track;
}

Track* SimDataReader::traceNC(const int &seg, const bool &print)
{
  DetectorHit *hit = getHit(CID_NC, seg);
  if( !hit ) return 0;
  if( print ){
    std::cout<<"===== SimDataReader::traceNC ====="<<std::endl;
    std::cout<<"> NC seg:"<<seg<<std::endl;
    std::cout<<"> particle : "<<parName(hit->pdg())<<std::endl;
    std::cout<<"> time : "<<hit->time()<<" dE : "<<hit->adc()<<std::endl;
  }
  Track *track = getMCTrack(hit->trackID());

  while( true ){
    bool NChit=false;
    for( int i=0; i<track->detectorHitLinkSize(); i++ ){
      if( detData->detectorHit(track->detectorHitLink(i))->detectorID()==CID_NC ){
	NChit = true;
      }
    }
    if( !NChit ) break;
    if( track->parentTrackID()==0 ) return track;
    track = getMCTrack(track->parentTrackID());
  }
  return trace(track, print);
}

bool SimDataReader::fillNpipi(TFile *f, CDSHitMan *cdsMan, const TLorentzVector &beam_lmom, CDSTrack *pim, const TLorentzVector &pim_lmom, 
			      CDSTrack *pip, const TLorentzVector &pip_lmom, const int &NCseg, const TLorentzVector &n_lmom, const double &mm)
{
  CrossSectionTable table = runHeader-> CStable();
  if( table.CSSize()<10 ){
    TH2F *h2;
    TH1F *h1;
    //  std::cout<<"  !!! SimDataReader::fillNpipi START !!!"<<std::endl;
    int reactionID = reacData-> ReactionID();

    if( reactionID==821 ){
      Track *sm_track = 0;
      Track *sp_track = 0;
      Track *k0_track = 0;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData-> track(i);
	if( track->pdgID()==3222 ) sp_track = track;
	if( track->pdgID()==3112 ) sm_track = track;
	if( track->pdgID()==311  ) k0_track = track;
      }

      if( sm_track ){
	h1 = (TH1F*)f-> Get("KN_MM_wN_Reac821_Sm"), h1-> Fill(mm);
	h2 = (TH2F*)f-> Get("KNpim_KNpip_MM_wN_Reac821_Sm"), h2-> Fill((beam_lmom+D_LMOM-pim_lmom).M(), (beam_lmom+D_LMOM-pip_lmom).M());
      }

      if( sp_track ){
	h1 = (TH1F*)f-> Get("KN_MM_wN_Reac821_Sp"), h1-> Fill(mm);
	h2 = (TH2F*)f-> Get("KNpim_KNpip_MM_wN_Reac821_Sp"), h2-> Fill((beam_lmom+D_LMOM-pim_lmom).M(), (beam_lmom+D_LMOM-pip_lmom).M());
      }

      if( k0_track ){
	h1 = (TH1F*)f-> Get("KN_MM_wN_Reac821_K0"), h1-> Fill(mm);
	h2 = (TH2F*)f-> Get("KNpim_KNpip_MM_wN_Reac821_K0"), h2-> Fill((beam_lmom+D_LMOM-pim_lmom).M(), (beam_lmom+D_LMOM-pip_lmom).M());
      }
    }

    DetectorHit *pim_hit=0;
    DetectorHit *pip_hit=0;
    int pim_seg, pip_seg;
    double pim_time, pip_time;
    pim-> GetCDHHit(cdsMan, pim_seg, pim_time);
    pip-> GetCDHHit(cdsMan, pip_seg, pip_time);

    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit = detData-> detectorHit(i);
      if( hit->detectorID()==CID_CDH ){
	if( hit->channelID()+1==pim_seg ) pim_hit=hit;
	if( hit->channelID()+1==pip_seg ) pip_hit=hit;
      }
    }

    if( pim_hit->pdg()!=-211 ){
      std::cout<<"!!! pi- hit not match !!!"<<std::endl;
      return false;
    }
    if( pip_hit->pdg()!=211 ){
      std::cout<<"!!! pi+ hit not match !!!"<<std::endl;
      return false;
    }

    Track *pim_track = getMCTrack(pim_hit->trackID());
    Track *pip_track = getMCTrack(pip_hit->trackID());
    while( pim_track->parentTrackID()!=0 ){
      Track *track  = getMCTrack(pim_track->parentTrackID());
      if( track-> pdgID()!=-211 ) break;
      pim_track = track;
    }
    while( pip_track->parentTrackID()!=0 ){
      Track *track = getMCTrack(pip_track->parentTrackID());
      if( track-> pdgID()!=211  ) break;
      pip_track = track;
    }

    int parentID = -999;
    for( int i=0; i<detData-> detectorHitSize(); i++ ){
      DetectorHit *hit = detData->detectorHit(i);
      if( hit->detectorID()==CID_NC && hit->channelID()+1==NCseg ){
	//      std::cout<<"> Detector Hit  dE : "<<hit->adc()<<"  time : "<<hit->time()<<std::endl;
	//      std::cout<<"  TrackID : "<<hit->trackID()<<"  particle : "<<parName(hit->pdg())<<"  reaction : "<<reacContainer->getReaction(hit->parentID())<<std::endl;
	parentID = hit->trackID();
      }
    }
    if( parentID<0 ){
      std::cout<<"!!! Detector Hit not found !!!"<<std::endl;
      return 0;
    }

    Track *n_track = getMCTrack(parentID);
    bool n_flag = false;

    bool n_scat = false;
    bool n_scat_Fiducial = false;
    bool n_scat_CDS  = false;
    bool n_scat_USWK = false;
    bool n_scat_BD   = false;
    bool n_scat_EG   = false;
    bool n_scat_CVC  = false;
    std::vector<TVector3> n_scat_pos;
    if( n_track-> pdgID()==2112 ) n_flag = true;
    while( n_track->parentTrackID()!=0 ){
      Track *track = getMCTrack(n_track->parentTrackID());
      if( track->pdgID()==2112 && n_track->pdgID()==2112 ){
	//      if( reacContainer-> getReaction(n_track->parentTrackID())=="neutronInelastic" ){
	if( track->processID() ){
	  TVector3 pos = 0.1*n_track->vertex();
	  //	std::cout<<"> neutron scatter pos("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<")"<<std::endl;
	  n_scat_pos.push_back(pos);

	  if( GeomTools::GetID(pos)!=CID_NC       ) n_scat=true;
	  if( GeomTools::GetID(pos)==CID_Fiducial ) n_scat_Fiducial=true;
	  if( GeomTools::GetID(pos)==CID_USWK     ) n_scat_USWK=true;
	  //	if( GeomTools::GetID(pos)==CID_BD       ) n_scat_BD=true;
	  if( -58.<pos.Z() && pos.Z()<58.         ) n_scat_CDS=true;
	  if( pos.Z()>1580.                       ) n_scat_BD=true;
	  if( GeomTools::GetID(pos)==CID_CVC      ) n_scat_CVC=true;
	  if( GeomTools::GetID(pos)==CID_Doraemon ){
	    if( 58.<pos.Z() && pos.Z()<74. ) n_scat_EG=true;
	  }
	}
      }
  
      if( !n_flag ){
	if( track->pdgID()==2112 ) n_flag=true;
      }
      else if( n_flag ){
	if( track->pdgID()!=2112 ) break;
      }
      n_track = track;
    }
    if( n_scat          ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_MC"), h1-> Fill(mm);
    if( n_scat_Fiducial ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_Fiducial_MC"), h1-> Fill(mm);
    if( n_scat_CDS      ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_CDS_MC"), h1-> Fill(mm);
    if( n_scat_EG       ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_EG_MC"), h1-> Fill(mm);
    if( n_scat_USWK     ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_USWK_MC"), h1-> Fill(mm);
    if( n_scat_BD       ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_BD_MC"), h1-> Fill(mm);
    if( n_scat_CVC      ) h1 = (TH1F*)f-> Get("KNpipi_MM_n_Scat_CVC_MC"), h1-> Fill(mm);

    for( int i=0; i<n_scat_pos.size(); i++ ){
      TVector3 pos = n_scat_pos[i];
      h1 = (TH1F*)f-> Get("N_Scat_Z_npipi_MC"), h1-> Fill(pos.Z());
      if( -58.<pos.Z() && pos.Z()<58. )  h2 = (TH2F*)f-> Get("N_Scat_CDS_npipi_MC"),  h2-> Fill(pos.X(), pos.Y());
      if( 58.<pos.Z() && pos.Z()<74. )   h2 = (TH2F*)f-> Get("N_Scat_EG_npipi_MC"),   h2-> Fill(pos.X(), pos.Y());
      if( 180.<pos.Z() && pos.Z()<320. ) h2 = (TH2F*)f-> Get("N_Scat_USWK_npipi_MC"), h2-> Fill(pos.X(), pos.Y());
    }

    if( !n_flag ){
      std::cout<<"!!! not found neutron !!!"<<std::endl;
      return false;
    }

    TVector3 pim_mom_MC = 0.001*pim_track->momentum();
    TVector3 pip_mom_MC = 0.001*pip_track->momentum();
    TVector3 n_mom_MC = 0.001*n_track->momentum();

    double pim_diff = pim_mom_MC.Mag()-pim_lmom.Vect().Mag();
    double pip_diff = pip_mom_MC.Mag()-pip_lmom.Vect().Mag();
    double n_diff   = n_mom_MC.Mag()-n_lmom.Vect().Mag();
#if 0
    std::cout<<"> d(K-, n pi+ pi-) : "<<mm<<" [GeV/c^{2}]"<<std::endl;
    std::cout<<">>> pi- : "<<pim_diff<<std::endl;
    std::cout<<">>> pi+ : "<<pip_diff<<std::endl;
    std::cout<<">>> n   : "<<n_diff<<std::endl;
    std::string str;
    std::cin>>str;
    if( str=="q" ) exit(-1);
    std::cout<<"  !!! SimDataReader::fillNpipi FINISH !!!"<<std::endl;
#endif
    h2 = (TH2F*)f-> Get("KNpipi_MM_pim_diff"), h2-> Fill(mm, pim_diff);
    h2 = (TH2F*)f-> Get("KNpipi_MM_pip_diff"), h2-> Fill(mm, pip_diff);
    h2 = (TH2F*)f-> Get("KNpipi_MM_n_diff"),   h2-> Fill(mm, n_diff);
  }
  return true;
}


int SimDataReader::NChitNeutronStatus(const int &NCseg)
{
  //  std::cout<<"===== SimDataReader::NChitNeutronStatus START NC seg : "<<NCseg<<" ====="<<std::endl;
  int parentID = -999;
  for( int i=0; i<detData-> detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==CID_NC && hit->channelID()+1==NCseg ){
      //      std::cout<<"> Detector Hit  dE : "<<hit->adc()<<"  time : "<<hit->time()<<std::endl;
      //      std::cout<<"  TrackID : "<<hit->trackID()<<"  particle : "<<parName(hit->pdg())<<"  reaction : "<<reacContainer->getReaction(hit->parentID())<<std::endl;
      parentID = hit->trackID();
    }
  }
  if( parentID<0 ){
    std::cout<<"!!! Detector Hit not found !!!"<<std::endl;
    return 0;
  }

  bool N_flag = false;
  bool N_Inel = false;
  bool Dora_Scat = false;
  bool USWK_Scat = false;
  bool BD_Scat = false;
  bool from_Decay = false;
  bool N_init  = false;

  while( true ){
    bool next_flag = false;
    for( int i=0; i<mcData->trackSize(); i++ ){
      Track *track = mcData->track(i);
      if( parentID==track->trackID() ){
	TVector3 vtx = 0.1*track-> vertex();
	//	std::cout<<" trackID : "<<track->trackID()<<" particle : "<<parName(track->pdgID())<<" reaction : "<<reacContainer->getReaction(track->parentTrackID())<<std::endl;
	//	std::cout<<"   Vertex ("<<vtx.X()<<", "<<vtx.Y()<<", "<<vtx.Z()<<")  CID : "<<GeomTools::GetID(vtx)<<std::endl;

	if( track->pdgID()==2112 ) N_flag = true;
	if( track->parentTrackID()==0 && track->pdgID()==2112 ) N_init=true;
	if( track->pdgID()==2112 && track->processID()==111 ){
	  N_Inel = true;
	  if( GeomTools::GetID(vtx)==CID_Doraemon ) Dora_Scat=true;
	  if( GeomTools::GetID(vtx)==CID_USWK     ) USWK_Scat=true;
	  if( 110<=GeomTools::GetID(vtx) && GeomTools::GetID(vtx)<120     ) BD_Scat=true;
	}
	if( track->pdgID()==2112 && track->processID()==1 ) from_Decay=true;

	parentID = track->parentTrackID();
	next_flag = true;
	break;
      }
    }
    if( !next_flag ){
      std::cout<<"!!! ParentTrack not found !!!"<<std::endl;
      return 0;
    }
    if( parentID==0 ) break;
  }
#if 0
  std::cout<<"> Neutron flag : "<<std::boolalpha<<N_flag<<std::endl;
  if( N_init     ) std::cout<<">  Initial Neutron"<<std::endl;
  if( from_Decay ) std::cout<<">  Neutron from Decay "<<std::endl;
  if( N_Inel    ) std::cout<<">     Neutron Inelastic"<<std::endl;
  if( Dora_Scat ) std::cout<<">     Neutron Scattered by Dora"<<std::endl;
  if( USWK_Scat ) std::cout<<">     Neutron Scattered by USWK"<<std::endl;
  if( BD_Scat   ) std::cout<<">     Neutron Scattered by Beam Dump"<<std::endl;
  std::cout<<"===== SimDataReader::NChitNeutronStatus FINISH ====="<<std::endl;
#endif
  if( !N_flag ) return -1; // not found neutron track;
  if( Dora_Scat ){
    if( N_init     ) return 1101;
    if( from_Decay ) return 1102;
    return 1103;
  }
  if( USWK_Scat ){
    if( N_init     ) return 1201;
    if( from_Decay ) return 1202;
    return 1203;
  }
  if( BD_Scat   ){
    if( N_init     ) return 1301;
    if( from_Decay ) return 1302;
    return 1303;
  }
  if( N_Inel ){
    if( N_init     ) return 1001;
    if( from_Decay ) return 1002;
    return 1003;
  }

  if( N_init ) return 1;
  if( from_Decay ) return 2;
  return 3;
}

bool SimDataReader::isHyperon(const int &pdg)
{
  bool flag = false;
  if( pdg==3122 || pdg==13122 || pdg==3124 || pdg==23122 || pdg==33122 || pdg==13124 || pdg==43122 || pdg==53122 || pdg==3126 || pdg==13126 || pdg==23124 ||
      pdg==3224 || pdg==3114  || pdg==3214 || pdg==3222  || pdg==3112  || pdg==3212 ) flag=true;
  if( flag ) std::cout<<parName(pdg)<<std::endl;

  return flag;
}

bool SimDataReader::isNeutronInElasticNC(const int &seg)
{
  DetectorHit *NChit = 0;
  double time = DBL_MAX;
  //  std::cout<<"========================"<<std::endl;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData->detectorHit(i);
    if( hit->detectorID()==CID_NC && hit->time()<time ){
      time = hit->time();
      NChit=hit;
    }
    // if( hit->detectorID()==OACID_NC ){
    //   std::cout<<" NC seg : "<<hit->channelID()+1<<std::endl;
    //   std::cout<<"   Particle : "<<parName(hit->pdg())<<std::endl;
    //   std::cout<<"   time : "<<hit->time()<<" dE : "<<hit->adc()<<std::endl;
    // }
  }
  //  std::cout<<"========================"<<std::endl;
  if( !NChit ) return false;

  return isNeutronInElastic(NChit);
}


bool SimDataReader::isNeutronInElastic(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isNeutronInElastic(hit);
}

bool SimDataReader::isNeutronInElastic(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();
  if( originPDG!=2112 ){
    std::cout<<"  !!! SimDataReader::isNeutronInElastic not neutron : "<<parName(originPDG)<<" !!!"<<std::endl;
    return false;
  }

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);
    if( track->processID()==111 ) return true;

    if( track->pdgID()!=2112 ) break;

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

bool SimDataReader::isFromGamma(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isFromGamma(hit);
}

bool SimDataReader::isFromGamma(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);

    if( track->pdgID()!=originPDG ){
      if( track->pdgID()==22 ) return true;
      else break;
    }
    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

bool SimDataReader::isFromE(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isFromE(hit);
}

bool SimDataReader::isFromE(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);

    if( track->pdgID()!=originPDG ){
      if( track->pdgID()==11 || track->pdgID()==-11 ) return true;
      else break;
    }
    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

bool SimDataReader::isInitial(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isInitial(hit);
}

bool SimDataReader::isInitial(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);
    if( track->pdgID()!=originPDG ) return false;

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return true;
}

bool SimDataReader::isSigmaDecay(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isSigmaDecay(hit);
}

bool SimDataReader::isSigmaDecay(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);
    if( track->pdgID()!=originPDG ){
      if( track->pdgID()==3112 || track->pdgID()==3222 || track->pdgID()==3212 ) return true;
      else break;
    }

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

bool SimDataReader::isPiDecay(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isPiDecay(hit);
}

bool SimDataReader::isPiDecay(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);
    if( track->pdgID()!=originPDG ){
      if( track->pdgID()==211 || track->pdgID()==111 || track->pdgID()==-211 ) return true;
      else break;
    }

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

bool SimDataReader::isLambdaDecay(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isLambdaDecay(hit);
}

bool SimDataReader::isLambdaDecay(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);
    if( track->pdgID()!=originPDG ){
      if( track->pdgID()==3122 ) return true;
      else break;
    }

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

bool SimDataReader::isL1520Decay(const int &cid, const int &seg)
{
  DetectorHit *hit = getHit(cid, seg);
  if( !hit ) return 0;

  return isL1520Decay(hit);
}

bool SimDataReader::isL1520Decay(DetectorHit *hit)
{
  const int originPDG = hit->pdg();
  const TVector3 originMom = hit->momentum();
  int depth=1;
  int parentID = hit->parentID();
  int trackID = hit->trackID();

  while( parentID>0 ){
    Track *track = getMCTrack(parentID);
    if( track->pdgID()!=originPDG ){
      if( track->pdgID()==3124 ) return true;
      else break;
    }

    trackID=parentID;
    parentID=track->parentTrackID();
    depth++;
  }

  return false;
}

void SimDataReader::get()
{
  fBeamPDG = reacData-> InitPDG(0);
  fTarPDG = reacData-> InitPDG(1);

  fBeamLmom = reacData-> GetInitParticle(0);
  fTarLmom = reacData-> GetInitParticle(1);

}

void SimDataReader::clear()
{
  //  std::cout<<"SimDataReader::clear"<<std::endl;
  fTrigger=false;
  fL1405_lmom.SetXYZT(0., 0., 0., 0.);
  fBeamPDG = INT_MIN;
  fTarPDG = INT_MIN;
  fBeamLmom.SetXYZT(0, 0, 0, 0);
  fTarLmom.SetXYZT(0, 0, 0, 0);

}

void SimDataReader::printReaction(const int &max_gen)
{
  std::cout<<"===== SimDataReader::printReaction START ====="<<std::endl;
  std::cout<<"> Event ID : "<<evHeader-> eventID()<<std::endl;
  std::cout<<"> Reaction ID : "<<reacData-> ReactionID()<<std::endl;
  std::cout<<">    "<<reacName()<<std::endl;
  std::cout<<"===== SimDataReader::printReaction START ====="<<std::endl;
}

std::string SimDataReader::reacName()
{
  std::string name;
  name += parName(fBeamPDG)+" ";
  name += parName(fTarPDG)+" -> ";
  int n1 = reacData->NParticle(0);
  int n2 = n1+reacData->NParticle(1);
  for( int i=0; i<n1; i++ ){
    name += parName(reacData->PDG(i))+" ";
  }
  for( int i=n1; i<n2; i++ ){
    name += parName(reacData->PDG(i))+"_s ";
  }
  return name;
}

void SimDataReader::initHist(TFile* f)
{
  f-> cd();

  std::cout<<"===== SimDataReader::initHist START  ====="<<std::endl;
  new TH1F("CDH_time_diff_true",  "CDH cluster time diff", 1000, -10.0, 10.0);
  new TH1F("CDH_time_diff_false", "CDH cluster time diff", 1000, -10.0, 10.0);

  new TH1F("ReactionID", "Reaction ID", 5000, -0.5, 4999.5);

  new TH1F("NC_overbeta_n_MC", "NC 1/#beta by neutron", 10000, 0.0, 10.0);

  new TH1F("NC_overbeta_n_init_MC", "NC 1/#beta by neutron", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_n_from_S_MC", "NC 1/#beta by neutron", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_n_from_L_MC", "NC 1/#beta by neutron", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_n_from_L1520_MC", "NC 1/#beta by neutron", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_n_acci_MC", "NC 1/#beta by neutron", 10000, 0.0, 10.0);

  new TH1F("NC_overbeta_gamma_MC", "NC 1/#beta by #gamma", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_gamma_init_MC", "NC 1/#beta by #gamma", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_gamma_from_S_MC", "NC 1/#beta by #gamma", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_gamma_from_pi_MC", "NC 1/#beta by #gamma", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_gamma_from_e_MC", "NC 1/#beta by #gamma", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_gamma_acci_MC", "NC 1/#beta by #gamma", 10000, 0.0, 10.0);

  new TH2F("KNpipi_MM_pim_diff", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs #pi^{-} mom diff", 400, 0.0, 2.0, 500, -0.5, 0.5);
  new TH2F("KNpipi_MM_pip_diff", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs #pi^{+} mom diff", 400, 0.0, 2.0, 500, -0.5, 0.5);
  new TH2F("KNpipi_MM_n_diff",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs n mom diff",       400, 0.0, 2.0, 750, -0.1, 1.4);

  new TH1F("N_Scat_Z_npipi_MC", "neutron scatter pos Z", 25000, -500, 2000);

  new TH2F("N_Scat_CDS_npipi_MC",  "neutron scatter pos XY at CDS",        200,  -100, 100,  200, -100, 100);
  new TH2F("N_Scat_EG_npipi_MC",   "neutron scatter pos XY at End Guard",  200 , -100, 100,  200, -100, 100);
  new TH2F("N_Scat_USWK_npipi_MC", "neutron scatter pos XY at USWK",       1000, -500, 500, 1000, -500, 500);

  new TH1F("KNpipi_MM_n_Scat_MC",          "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_n_Scat_EG_MC",       "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_n_Scat_CDS_MC",      "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_n_Scat_Fiducial_MC", "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_n_Scat_USWK_MC",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_n_Scat_BD_MC",       "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_n_Scat_CVC_MC",      "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);

  new TH2F("n_mom_2D_MMSA_MC", "n mom true vs MMSA", 500, 0., 2.0, 500, 0., 2.0);
  new TH1F("n_mom_diff_MMSA_MC", "n mom true - MMSA", 1000, -0.5, 0.5);
#if 1
  CrossSectionTable table = runHeader-> CStable();
  std::cout<<"> CSTable size : "<<table.CSSize()<<std::endl;
  for( int i=0; i<table.CSSize(); i++ ){
    const CrossSection cs=table.CS(i);
    int reacID=cs.Id();
    new TH1F(Form("KNpipi_MM_Reac%d", reacID),       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F(Form("KNpipi_MM_woAll_Reac%d", reacID), "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F(Form("KNpipi_MM_wK0_Reac%d", reacID),   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F(Form("KNpipi_MM_wSm_Reac%d", reacID),   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F(Form("KNpipi_MM_wSp_Reac%d", reacID),   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F(Form("KN_MM_Reac%d", reacID), "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
    new TH1F(Form("KN_MM_wN_Reac%d", reacID), "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
  }
  if( table.CSSize()<10 ){
    for( int i=0; i<table.CSSize(); i++ ){
      const CrossSection cs=table.CS(i);
      int reacID=cs.Id();
      std::cout<<">>> reaction ID : "<<reacID<<std::endl;
      if( reacID==2 ){
	new TH1F("Reac2_km_ang_MC",    "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac2_p_ang_MC",     "p cos#theta Lab Flame",     1000, -1.0, 1.0);
	new TH1F("Reac2_km_ang_CM_MC", "K^{-} cos#theta CM Flame",  1000, -1.0, 1.0);
	new TH1F("Reac2_p_ang_CM_MC",  "p cos#theta CM Flame",      1000, -1.0, 1.0);
      }
      else if( reacID==272 ){
	new TH1F("Reac272_L_ang_MC", "#Lambda cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_pip_ang_MC", "#pi^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_L_ang_CM_MC", "#Lambda cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_pip_ang_CM_MC", "#pi^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac272_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==280 ){
	new TH1F("Reac280_L_ang_MC", "#Lambda cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac280_pip_ang_MC", "#pi^{+} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac280_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac280_L_ang_CM_MC", "#Lambda cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac280_pip_ang_CM_MC", "#pi^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac280_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==330 ){
	new TH1F("Reac330_L_ang_MC", "#Lambda cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac330_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac330_L_ang_CM_MC", "#Lambda cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac330_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==349 ){
	new TH1F("Reac349_L_ang_MC", "#Lambda cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac349_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac349_L_ang_CM_MC", "#Lambda cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac349_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==361 ){
	new TH1F("Reac361_L_ang_MC", "#Lambda cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac361_eta_ang_MC", "#eta cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac361_L_ang_CM_MC", "#Lambda cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac361_eta_ang_CM_MC", "#eta cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==364 ){
	new TH1F("Reac364_L_ang_MC", "#Lambda cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac364_rho0_ang_MC", "#rho^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac364_L_ang_CM_MC", "#Lambda cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac364_rho0_ang_CM_MC", "#rho^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==432 ){
	new TH1F("Reac432_L1405_ang_MC", "#Lambda(1405) cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac432_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac432_L1405_ang_CM_MC", "#Lambda(1405) cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac432_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==442 ){
	new TH1F("Reac442_L1520_ang_MC", "#Lambda(1520) cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac442_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac442_L1520_ang_CM_MC", "#Lambda(1520) cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac442_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==455 ){
	new TH1F("Reac442_L1690_ang_MC", "#Lambda(1690) cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac442_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac442_L1690_ang_CM_MC", "#Lambda(1690) cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac442_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==464 ){
	new TH1F("Reac464_Sp_ang_MC", "#sigma^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac464_pip_ang_MC", "#pi^{+} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac464_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac464_Sp_ang_CM_MC", "#sigma^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac464_pip_ang_CM_MC", "#pi^{+} cos#thera CM Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac464_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==467 ){
	new TH1F("Reac467_Sp_ang_MC", "#sigma^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac467_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac467_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac467_Sp_ang_CM_MC", "#sigma^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac467_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac467_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==476 ){
	new TH1F("Reac476_n_spec_mom_MC", "n_{Spec} mom", 1000, 0.0, 1.0);
	new TH1F("Reac476_Sp_ang_MC", "#sigma^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac476_Sp_ang_CM_MC", "#sigma^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 

	new TH1F("Reac476_decay_n_ang_MC", "n decay from #Sigma cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_decay_n_mom_MC", "n decay from #Sigma mom Fab Flame", 2000, 0.0, 2.0);

	new TH1F("Reac476_Sp_ang_MC_n_hit", "#sigma^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_MC_n_hit", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac476_Sp_ang_CM_MC_n_hit", "#sigma^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_CM_MC_n_hit", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 

	new TH1F("Reac476_decay_n_ang_MC_n_hit", "n decay from #Sigma cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_decay_n_mom_MC_n_hit", "n decay from #Sigma mom Fab Flame", 2000, 0.0, 2.0);

	new TH1F("Reac476_Sp_ang_MC_n_scat", "#sigma^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_MC_n_scat", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac476_Sp_ang_CM_MC_n_scat", "#sigma^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_CM_MC_n_scat", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 

	new TH1F("Reac476_decay_n_ang_MC_n_scat", "n decay from #Sigma cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_decay_n_mom_MC_n_scat", "n decay from #Sigma mom Fab Flame", 2000, 0.0, 2.0);

	new TH1F("Reac476_Sp_ang_MC_npipi_hit", "#sigma^{+} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_MC_npipi_hit", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac476_Sp_ang_CM_MC_npipi_hit", "#sigma^{+} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_pim_ang_CM_MC_npipi_hit", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 

	new TH1F("Reac476_decay_n_ang_MC_npipi_hit", "n decay from #Sigma cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac476_decay_n_mom_MC_npipi_hit", "n decay from #Sigma mom Fab Flame", 2000, 0.0, 2.0);
      }
      else if( reacID==506 ){
	new TH1F("Reac506_S0_ang_MC", "#sigma^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac506_pip_ang_MC", "#pi^{+} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac506_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac506_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac506_S0_ang_CM_MC", "#sigma^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac506_pip_ang_CM_MC", "#pi^{+} cos#thera CM Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac506_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac506_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==508 ){
	new TH1F("Reac508_S0_ang_MC", "#sigma^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac508_pip_ang_MC", "#pi^{+} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac508_pim_ang_MC", "#pi^{-} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac508_S0_ang_CM_MC", "#sigma^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac508_pip_ang_CM_MC", "#pi^{+} cos#thera CM Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac508_pim_ang_CM_MC", "#pi^{-} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==514 ){
	new TH1F("Reac514_S0_ang_MC", "#sigma^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac514_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac514_S0_ang_CM_MC", "#sigma^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac514_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==516 ){
	new TH1F("Reac516_S0_ang_MC", "#sigma^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac516_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac516_S0_ang_CM_MC", "#sigma^{0} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac516_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==536 ){
	new TH1F("Reac536_n_spec_mom_MC",  "n_{Spec} mom", 1000, 0.0, 1.0);
	new TH1F("Reac536_Sm_ang_MC", "#sigma^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac536_Sm_ang_CM_MC", "#sigma^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
  
	new TH1F("Reac536_decay_n_ang_MC", "n decay from #Sigma cos#theta Fab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_decay_n_mom_MC", "n decay from #Sigma mom Fab Flame", 2000, 0.0, 2.0);

	new TH1F("Reac536_Sm_ang_MC_n_hit", "#sigma^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_MC_n_hit", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac536_Sm_ang_CM_MC_n_hit", "#sigma^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_CM_MC_n_hit", "#pi^{0} cos#thera CM Flame", 2000, 0.0, 2.0); 

	new TH1F("Reac536_decay_n_ang_MC_n_hit", "n decay from #Sigma cos#theta Fab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_decay_n_mom_MC_n_hit", "n decay from #Sigma mom Fab Flame", 1000, -1.0, 1.0);

	new TH1F("Reac536_Sm_ang_MC_n_scat", "#sigma^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_MC_n_scat", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac536_Sm_ang_CM_MC_n_scat", "#sigma^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_CM_MC_n_scat", "#pi^{0} cos#thera CM Flame", 2000, 0.0, 2.0); 

	new TH1F("Reac536_decay_n_ang_MC_n_scat", "n decay from #Sigma cos#theta Fab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_decay_n_mom_MC_n_scat", "n decay from #Sigma mom Fab Flame", 1000, -1.0, 1.0);

	new TH1F("Reac536_Sm_ang_MC_npipi_hit", "#sigma^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_MC_npipi_hit", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac536_Sm_ang_CM_MC_npipi_hit", "#sigma^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_pip_ang_CM_MC_npipi_hit", "#pi^{0} cos#thera CM Flame", 2000, 0.0, 2.0); 

	new TH1F("Reac536_decay_n_ang_MC_npipi_hit", "n decay from #Sigma cos#theta Fab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac536_decay_n_mom_MC_npipi_hit", "n decay from #Sigma mom Fab Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==537 ){
	new TH1F("Reac537_Sm_ang_MC", "#sigma^{-} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac537_pip_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac537_pi0_ang_MC", "#pi^{0} cos#thera Lab Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac537_Sm_ang_CM_MC", "#sigma^{-} cos#thera CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac537_pip_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
	new TH1F("Reac537_pi0_ang_CM_MC", "#pi^{0} cos#thera CM Flame", 1000, -1.0, 1.0); 
      }
      else if( reacID==541 ){
	new TH1F("Reac541_Sm_ang_MC", "#Sigma^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac541_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac541_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac541_Sm_ang_CM_MC", "#Sigma^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac541_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac541_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==556 ){
	new TH1F("Reac556_Sm_ang_MC", "#Sigma^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac556_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac556_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac556_Sm_ang_CM_MC", "#Sigma^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac556_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac556_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==579 ){
	new TH1F("Reac579_S1385p_ang_MC", "#Sigma(1385)^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac579_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac579_S1385p_ang_CM_MC", "#Sigma(1385)^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac579_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==592 ){
	new TH1F("Reac579_S13850_ang_MC", "#Sigma(1385)^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac579_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac579_S13850_ang_CM_MC", "#Sigma(1385)^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac579_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==597 ){
	new TH1F("Reac597_S1385m_ang_MC", "#Sigma(1385)^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac597_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac597_S1385m_ang_CM_MC", "#Sigma(1385)^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac597_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==635 ){
	new TH1F("Reac635_p_ang_MC", "p cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_p_ang_CM_MC", "p cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac635_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==670 ){
	new TH1F("Reac670_p_ang_MC", "p cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac670_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac670_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac670_p_ang_CM_MC", "p cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac670_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac670_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==679 ){
	new TH1F("Reac679_p_ang_MC", "p cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac679_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac679_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac679_p_ang_CM_MC", "p cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac679_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac679_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==707 ){
	new TH1F("Reac707_p_ang_MC", "p cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac707_ksm_ang_MC", "K*^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac707_p_ang_CM_MC", "p cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac707_ksm_ang_CM_MC", "K*^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==714 ){
	new TH1F("Reac714_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac714_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==716 ){
	new TH1F("Reac716_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac716_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac716_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac716_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac716_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac716_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==729 ){
	new TH1F("Reac729_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac729_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac729_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac729_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac729_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac729_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==735 ){
	new TH1F("Reac735_n_spec_mom_MC", "n_{Spec} mom", 1000, 0.0, 1.0);

	new TH1F("Reac735_k0_ang_MC_trig", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_k0_ang_CM_MC_trig", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);

	new TH1F("Reac735_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_mom_MC", "n cos#theta Lab Flame", 1000, 0.0, 2.0);
	new TH1F("Reac735_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);

	new TH1F("Reac735_n_ang_MC_hit", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_mom_MC_hit", "n cos#theta Lab Flame", 2000, 0.0, 2.0);
	new TH1F("Reac735_k0_ang_MC_hit", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_ang_CM_MC_hit", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_k0_ang_CM_MC_hit", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);

	new TH1F("Reac735_n_ang_MC_n_hit", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_mom_MC_n_hit", "n cos#theta Lab Flame", 2000, 0.0, 2.0);
	new TH1F("Reac735_k0_ang_MC_n_hit", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_ang_CM_MC_n_hit", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_k0_ang_CM_MC_n_hit", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);

	new TH1F("Reac735_n_ang_MC_n_scat", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_mom_MC_n_scat", "n cos#theta Lab Flame", 2000, 0.0, 2.0);
	new TH1F("Reac735_k0_ang_MC_n_scat", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_ang_CM_MC_n_scat", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_k0_ang_CM_MC_n_scat", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);

	new TH1F("Reac735_n_ang_MC_npipi_hit", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_mom_MC_npipi_hit", "n cos#theta Lab Flame", 2000, 0.0, 2.0);
	new TH1F("Reac735_k0_ang_MC_npipi_hit", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_n_ang_CM_MC_npipi_hit", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac735_k0_ang_CM_MC_npipi_hit", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==737 ){
	new TH1F("Reac737_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac737_ks0_ang_MC", "K*^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac737_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac737_ks0_ang_CM_MC", "K*^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==767 ){
	new TH1F("Reac767_Deltap_ang_MC", "#Delta^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac767_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac767_Deltap_ang_CM_MC", "#Delta^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac767_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==771 ){
	new TH1F("Reac771_Delta0_ang_MC", "#Delta^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac771_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac771_Delta0_ang_CM_MC", "#Delta^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac771_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==821 ){
	std::cout<<" Init Hist Reaction821"<<std::endl;
	new TH1F("Reac821_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_n_mom_MC", "n mom Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);

	new TH1F("KN_MM_wN_Reac821_Sm", "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
	new TH1F("KN_MM_wN_Reac821_Sp", "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
	new TH1F("KN_MM_wN_Reac821_K0", "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);

	new TH1F("Reac821_km_mom_MC",    "K^{-} mom Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_km_mom_MC_Sm", "K^{-} mom Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_km_mom_MC_Sp", "K^{-} mom Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac821_km_mom_MC_K0", "K^{-} mom Lab Flame", 1000, -1.0, 1.0);

	new TH2F("KNpim_KNpip_MM_wN_Reac821_Sm", "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
	new TH2F("KNpim_KNpip_MM_wN_Reac821_Sp", "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
	new TH2F("KNpim_KNpip_MM_wN_Reac821_K0", "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);

	initHist_knEl();
      }
      else if( reacID==851 ){
	new TH1F("Reac851_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac851_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac851_L_ang_MC", "#Lambda cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac851_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac851_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac851_L_ang_CM_MC", "#Lambda cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==856 ){
	new TH1F("Reac856_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac856_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac856_L_ang_MC", "#Lambda cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac856_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac856_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac856_L_ang_CM_MC", "#Lambda cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==862 ){
	new TH1F("Reac862_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac862_L_ang_MC", "#Lambda cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac862_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac862_L_ang_CM_MC", "#Lambda cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==893 ){
	new TH1F("Reac893_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac893_L1405_ang_MC", "#Lambda(1405) cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac893_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac893_L1405_ang_CM_MC", "#Lambda(1405) cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==897 ){
	new TH1F("Reac897_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac897_L1520_ang_MC", "#Lambda(1520) cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac897_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac897_L1520_ang_CM_MC", "#Lambda(1520) cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==912 ){
	new TH1F("Reac912_Sp_ang_MC", "#Sigma^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac912_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac912_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac912_Sp_ang_CM_MC", "#Sigma^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac912_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac912_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==915 ){
	new TH1F("Reac915_Sp_ang_MC", "#Sigma^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac915_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac915_Sp_ang_CM_MC", "#Sigma^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac915_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==922 ){
	new TH1F("Reac922_S0_ang_MC", "#Sigma^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac922_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac922_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac922_S0_ang_CM_MC", "#Sigma^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac922_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac922_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==923 ){
	new TH1F("Reac923_S0_ang_MC", "#Sigma^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac923_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac923_S0_ang_CM_MC", "#Sigma^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac923_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==929 ){
	new TH1F("Reac929_Sm_ang_MC", "#Sigma^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_Sm_ang_CM_MC", "#Sigma^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac929_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==931 ){
	new TH1F("Reac931_Sm_ang_MC", "#Sigma^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac931_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac931_pip_ang_MC", "#pi^{+} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac931_Sm_ang_CM_MC", "#Sigma^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac931_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac931_pip_ang_CM_MC", "#pi^{+} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==938 ){
	new TH1F("Reac938_Sm_ang_MC", "#Sigma^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac938_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac938_Sm_ang_CM_MC", "#Sigma^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac938_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==946 ){
	new TH1F("Reac946_Sm_ang_MC", "#Sigma^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac946_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac946_Sm_ang_CM_MC", "#Sigma^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac946_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==984 ){
	new TH1F("Reac984_S1385_ang_MC", "#Sigma(1385)^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac984_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac984_S1385_ang_CM_MC", "#Sigma(1385)^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac984_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==990 ){
	new TH1F("Reac990_S1385m_ang_MC", "#Sigma(1385)^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac990_pi0_ang_MC", "#pi^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac990_S1385m_ang_CM_MC", "#Sigma(1385)^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac990_pi0_ang_CM_MC", "#pi^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==1025 ){
	new TH1F("Reac1025_p_ang_MC", "p cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1025_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1025_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1025_p_ang_CM_MC", "p cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1025_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1025_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==1041 ){
	new TH1F("Reac1041_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1041_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1041_pim_ang_MC", "#pi^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1041_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1041_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1041_pim_ang_CM_MC", "#pi^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==1045 ){
	new TH1F("Reac1045_n_ang_MC", "n cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1045_ksm_ang_MC", "K*^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1045_n_ang_CM_MC", "n cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1045_ksm_ang_CM_MC", "K*^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==1055 ){
	new TH1F("Reac1055_D0_ang_MC", "D^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1055_km_ang_MC", "K^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1055_D0_ang_CM_MC", "D^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1055_km_ang_CM_MC", "K^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==1061 ){
	new TH1F("Reac1061_Dm_ang_MC", "D^{-} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1061_k0_ang_MC", "K^{0} cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1061_Dm_ang_CM_MC", "D^{-} cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1061_k0_ang_CM_MC", "K^{0} cos#theta CM Flame", 1000, -1.0, 1.0);
      }
      else if( reacID==1098 ){
	std::cout<<"> K- d -> n L(1405)     Final : "<<cs.FinlPdgSize()<<" Spec : "<<cs.SpecPdgSize()<<std::endl;
	new TH1F("Reac1098_n_ang_MC",     "neutron cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1098_l1405_ang_MC", "neutron cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1098_n_ang_CM_MC",     "neutron cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1098_l1405_ang_CM_MC", "neutron cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac1098_l1405_mass_MC",    "#Lambda(1405) mass", 1000, 1.0, 2.0);
      }
      else if( reacID==3098 ){
	new TH1F("Reac3098_n_ang_MC",     "neutron cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac3098_l1405_ang_MC", "neutron cos#theta Lab Flame", 1000, -1.0, 1.0);
	new TH1F("Reac3098_n_ang_CM_MC",     "neutron cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac3098_l1405_ang_CM_MC", "neutron cos#theta CM Flame", 1000, -1.0, 1.0);
	new TH1F("Reac3098_l1405_mass_MC",    "#Lambda(1405) mass", 1000, 1.0, 2.0);
      }
    
      if( reacID==3098 || reacID==1098 ){
	new TH1F("Trig_mass",         "Trig_mass", 1000, 1.0, 2.0);
	new TH1F("Trig_mass_hit",     "Trig_mass", 1000, 1.0, 2.0);
	new TH1F("Eff_pipi",          "Eff #pi^{+} #pi^{-} hit", 1000, 1.0, 2.0);
	new TH1F("Eff_pipi_woK0",     "Eff #pi^{+} #pi^{-} hit", 1000, 1.0, 2.0);
	new TH1F("Eff_pipi_woSf",     "Eff #pi^{+} #pi^{-} hit", 1000, 1.0, 2.0);
	new TH1F("Eff_mmN",           "Eff \"n\"",   1000, 1.0, 2.0);
	new TH1F("Eff_mmN_woK0",      "Eff \"n\"",   1000, 1.0, 2.0);
	new TH1F("Eff_mmN_woSf",      "Eff \"n\"",   1000, 1.0, 2.0);
	new TH1F("Eff_mmN_woAll",     "Eff \"n\" w/o All",   1000, 1.0, 2.0);
	new TH1F("Eff_Charged",       "Eff Charged Mode", 1000, 1.0, 2.0);
	new TH1F("Eff_mmN_woRcS",     "Eff \"n\" w/o #Sigma", 1000, 1.0, 2.0);
	new TH1F("Eff_Sp",            "Eff_Sp_mass", 1000, 1.0, 2.0);
	new TH1F("Eff_S0",            "Eff_S0_mass", 1000, 1.0, 2.0);
	new TH1F("Eff_Sm",            "Eff_Sm_mass", 1000, 1.0, 2.0);

	new TH1F("Eff_ppim_hit_MC", "Eff_Sm_mass", 1000, 1.0, 2.0);

	new TH1F("Eff_Lpim",   "Eff #Lambda #pi^{-} hit",            1000, 1.0, 2.0);
	new TH1F("Eff_Lpim_p", "Eff d(K^{-}, #Lambda #pi^{-})\"X\"", 1000, 1.0, 2.0);
	new TH1F("Eff_Lpim_p_hit", "Eff d(K^{-}, #Lambda #pi^{-})\"X\"", 1000, 1.0, 2.0);
      }
    }
  }
#endif
  std::cout<<"===== SimDataReader::initHist FINISH ====="<<std::endl;
}

void SimDataReader::fill(TFile *f)
{
  f-> cd();
  TH1F *h1;
  TH2F *h2;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get("ReactionID"), h1-> Fill(reactionID);

  bool eff_NC=false;
  bool CVC_flag=false;
  bool BVC_flag=false;
  bool PC_flag=false;
  
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData-> detectorHit(i);
    if( hit->detectorID()==CID_NC && hit->adc()>NC_THRE ) eff_NC=true;
    if( hit->detectorID()==CID_BVC && hit->adc()>0.05 ) BVC_flag=true;
    if( hit->detectorID()==CID_CVC && hit->adc()>0.1 ) CVC_flag=true;
    if( hit->detectorID()==CID_PC  && hit->adc()>0.1 ) PC_flag=true;
  }
  if( BVC_flag || CVC_flag ) eff_NC=false;

  bool S1385_flag = false;
  if( reactionID==1098 || reactionID==3098 ){
    for( int i=0; i<reacData->ParticleSize(); i++ ){
      if( reacData-> PDG(i)==13122 ) fL1405_lmom = 0.001*reacData->GetParticle(i);
      if( reacData-> PDG(i)==3114 ){
	S1385_flag=true;
	fL1405_lmom = 0.001*reacData->GetParticle(i);
      }
    }

    if( S1385_flag ){
      bool PCCVC_flag=false;
      for( int i=0; i<detData->detectorHitSize(); i++ ){
	DetectorHit *hit=detData->detectorHit(i);
	if( hit->detectorID()==CID_PC || hit->detectorID()==CID_CVC ){
	  if( hit->pdg()==2212 && hit->parentID()==0 ){
	    PCCVC_flag=true;
	  }
	}
      }
      if( PCCVC_flag ){
	TH1F *h1 = (TH1F*)f-> Get("Trig_mass_hit");
	h1-> Fill(fL1405_lmom.M());
      }

      if( PC_flag || CVC_flag ){
	h1 = (TH1F*)f-> Get("Trig_mass"), h1-> Fill(fL1405_lmom.M());
      }
    }

    if( !S1385_flag && eff_NC ){
      h1 = (TH1F*)f-> Get("Trig_mass"), h1-> Fill(fL1405_lmom.M());
    }
  }

#if 1
  CrossSectionTable table = runHeader-> CStable();
  if( table.CSSize()<10 ){
    if( reactionID==2 ){
      TLorentzVector km_lmom = reacData->GetParticle(0);
      TLorentzVector p_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(km_lmom+p_lmom).BoostVector();
      TLorentzVector km_lmom_CM = km_lmom;
      TLorentzVector p_lmom_CM = p_lmom;
      km_lmom_CM.Boost(LabToCM), p_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac2_km_ang_MC"), h1-> Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac2_p_ang_MC"), h1-> Fill(p_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac2_km_ang_CM_MC"), h1-> Fill(km_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac2_p_ang_CM_MC"), h1-> Fill(p_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==272 ){
      TLorentzVector l_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(2);
      TLorentzVector pim_lmom = reacData->GetParticle(3);
      TVector3 LabToCM = -(l_lmom+pip_lmom+pi0_lmom+pim_lmom).BoostVector();
      TLorentzVector l_lmom_CM = l_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      l_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac272_L_ang_MC"), h1-> Fill(l_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_L_ang_CM_MC"), h1-> Fill(l_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac272_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==280 ){
      TLorentzVector l_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(l_lmom+pip_lmom+pim_lmom).BoostVector();
      TLorentzVector l_lmom_CM = l_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      l_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac280_L_ang_MC"), h1-> Fill(l_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac280_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac280_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac280_L_ang_CM_MC"), h1-> Fill(l_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac280_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac280_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==330 ){
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TLorentzVector l_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(pi0_lmom+l_lmom).BoostVector();
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector l_lmom_CM = l_lmom;
      pi0_lmom_CM.Boost(LabToCM), l_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac330_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac330_L_ang_MC"), h1-> Fill(l_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac330_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac330_L_ang_CM_MC"), h1-> Fill(l_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==349 ){
      TLorentzVector l_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom2 = reacData->GetParticle(2);
      TVector3 LabToCM = -(pi0_lmom+pi0_lmom2+l_lmom).BoostVector();
      TLorentzVector pi0_lmom2_CM = pi0_lmom2;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector l_lmom_CM = l_lmom;
      pi0_lmom_CM.Boost(LabToCM), pi0_lmom2_CM.Boost(LabToCM), l_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac349_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta()), h1->Fill(pi0_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac349_L_ang_MC"), h1-> Fill(l_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac349_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta()), h1-> Fill(pi0_lmom2_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac349_L_ang_CM_MC"), h1-> Fill(l_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==361 ){
      TLorentzVector eta_lmom = reacData->GetParticle(0);
      TLorentzVector l_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(eta_lmom+l_lmom).BoostVector();
      TLorentzVector eta_lmom_CM = eta_lmom;
      TLorentzVector l_lmom_CM = l_lmom;
      eta_lmom_CM.Boost(LabToCM), l_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac361_eta_ang_MC"), h1-> Fill(eta_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac361_L_ang_MC"), h1-> Fill(l_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac361_eta_ang_CM_MC"), h1-> Fill(eta_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac361_L_ang_CM_MC"), h1-> Fill(l_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==364 ){
      TLorentzVector rho0_lmom = reacData->GetParticle(0);
      TLorentzVector l_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(rho0_lmom+l_lmom).BoostVector();
      TLorentzVector rho0_lmom_CM = rho0_lmom;
      TLorentzVector l_lmom_CM = l_lmom;
      rho0_lmom_CM.Boost(LabToCM), l_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac364_rho0_ang_MC"), h1-> Fill(rho0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac364_L_ang_MC"), h1-> Fill(l_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac364_rho0_ang_CM_MC"), h1-> Fill(rho0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac364_L_ang_CM_MC"), h1-> Fill(l_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==432 ){
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TLorentzVector l1405_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(pi0_lmom+l1405_lmom).BoostVector();
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector l1405_lmom_CM = l1405_lmom;
      pi0_lmom_CM.Boost(LabToCM), l1405_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac432_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac432_L1405_ang_MC"), h1-> Fill(l1405_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac432_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac432_L1405_ang_CM_MC"), h1-> Fill(l1405_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==442 ){
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TLorentzVector l1520_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(pi0_lmom+l1520_lmom).BoostVector();
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector l1520_lmom_CM = l1520_lmom;
      pi0_lmom_CM.Boost(LabToCM), l1520_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac442_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac442_L1520_ang_MC"), h1-> Fill(l1520_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac442_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac442_L1520_ang_CM_MC"), h1-> Fill(l1520_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==455 ){
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TLorentzVector l1690_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(pi0_lmom+l1690_lmom).BoostVector();
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector l1690_lmom_CM = l1690_lmom;
      pi0_lmom_CM.Boost(LabToCM), l1690_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac455_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac455_L1690_ang_MC"), h1-> Fill(l1690_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac455_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac455_L1690_ang_CM_MC"), h1-> Fill(l1690_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==464 ){
      TLorentzVector Sp_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TLorentzVector pim_lmom2 = reacData->GetParticle(3);
      TVector3 LabToCM = -(Sp_lmom+pip_lmom+pim_lmom).BoostVector();
      TLorentzVector Sp_lmom_CM = Sp_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector pim_lmom2_CM = pim_lmom2;
      Sp_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), pim_lmom2_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac464_Sp_ang_MC"), h1-> Fill(Sp_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac464_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac464_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta()), h1-> Fill(pim_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac464_Sp_ang_CM_MC"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac464_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac464_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta()), h1-> Fill(pim_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==467 ){
      TLorentzVector Sp_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(Sp_lmom+pi0_lmom+pim_lmom).BoostVector();
      TLorentzVector Sp_lmom_CM = Sp_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sp_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac467_Sp_ang_MC"), h1-> Fill(Sp_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac467_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac467_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac467_Sp_ang_CM_MC"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac467_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac467_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==476 ){
      if( reacData->ParticleSize()>2 ){
	TLorentzVector n_spec_lmom = reacData->GetParticle(2);
	h1 = (TH1F*)f-> Get("Reac476_n_spec_mom_MC"), h1-> Fill(n_spec_lmom.Vect().Mag());
      }
      TLorentzVector Sp_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(Sp_lmom+pim_lmom).BoostVector();
      TLorentzVector Sp_lmom_CM = Sp_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sp_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac476_Sp_ang_MC"), h1-> Fill(Sp_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac476_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac476_Sp_ang_CM_MC"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac476_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());

      bool pim_hit = false;
      bool pip_hit = false;
      bool n_hit = false;
      for( int i=0; i<detData->detectorHitSize(); i++ ){
	DetectorHit *hit = detData-> detectorHit(i);
	if( hit-> detectorID()==CID_CDH ){
	  if( hit-> pdg()== 211 ) pip_hit=true;
	  if( hit-> pdg()==-211 ) pim_hit=true;
	}
	if( hit-> detectorID()==CID_NC && hit->pdg()==2112 ) n_hit = true;
      }
      if( n_hit ){
	h1 = (TH1F*)f-> Get("Reac476_Sp_ang_MC_n_hit"), h1-> Fill(Sp_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac476_pim_ang_MC_n_hit"), h1-> Fill(pim_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac476_Sp_ang_CM_MC_n_hit"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac476_pim_ang_CM_MC_n_hit"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      }
      if( pim_hit && pip_hit && n_hit ){
	h1 = (TH1F*)f-> Get("Reac476_Sp_ang_MC_npipi_hit"), h1-> Fill(Sp_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac476_pim_ang_MC_npipi_hit"), h1-> Fill(pim_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac476_Sp_ang_CM_MC_npipi_hit"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac476_pim_ang_CM_MC_npipi_hit"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      }

      int sigmaID=-1;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( track->pdgID()==3222 && track->parentTrackID()==0 ) sigmaID=track->trackID();
      }
      if( sigmaID>0 ){
	Track *n_track=0;
	for( int i=0; i<mcData->trackSize(); i++ ){
	  Track *track = mcData->track(i);
	  if( track->pdgID()==2112 && track->parentTrackID()==sigmaID ) n_track=track;
	}
	if( n_track ){
	  bool n_scat_flag = false;
	  if( n_track->processID()==111 ){
	    for( int i=0; i<mcData->trackSize(); i++ ){
	      Track *track = mcData-> track(i);
	      if( track->parentTrackID()==n_track->trackID() && track->pdgID()==2112 ){
		TVector3 vtx = 0.1*track->vertex();
		if( -58.<vtx.Z() && vtx.Z()<58. && vtx.Perp()<15. ) n_scat_flag=true;
	      }
	    }
	  }
	  h1 = (TH1F*)f-> Get("Reac476_decay_n_ang_MC"), h1-> Fill(n_track->momentum().CosTheta());
	  h1 = (TH1F*)f-> Get("Reac476_decay_n_mom_MC"), h1-> Fill(0.001*n_track->momentum().Mag());
	  if( n_scat_flag ){
	    h1 = (TH1F*)f-> Get("Reac476_Sp_ang_MC_n_scat"), h1-> Fill(Sp_lmom.Vect().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac476_pim_ang_MC_n_scat"), h1-> Fill(pim_lmom.Vect().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac476_Sp_ang_CM_MC_n_scat"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac476_pim_ang_CM_MC_n_scat"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());

	    h1 = (TH1F*)f-> Get("Reac476_decay_n_ang_MC_n_scat"), h1-> Fill(n_track->momentum().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac476_decay_n_mom_MC_n_scat"), h1-> Fill(0.001*n_track->momentum().Mag());
	  }
	  if( n_hit ){
	    h1 = (TH1F*)f-> Get("Reac476_decay_n_ang_MC_n_hit"), h1-> Fill(n_track->momentum().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac476_decay_n_mom_MC_n_hit"), h1-> Fill(0.001*n_track->momentum().Mag());
	  }
	  if( pim_hit && pip_hit && n_hit ){
	    h1 = (TH1F*)f-> Get("Reac476_decay_n_ang_MC_npipi_hit"), h1-> Fill(n_track->momentum().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac476_decay_n_mom_MC_npipi_hit"), h1-> Fill(0.001*n_track->momentum().Mag());
	  }
	}
      }
    }
    else if( reactionID==506 ){
      TLorentzVector S0_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(2);
      TLorentzVector pim_lmom = reacData->GetParticle(3);
      TVector3 LabToCM = -(S0_lmom+pip_lmom+pi0_lmom+pim_lmom).BoostVector();
      TLorentzVector S0_lmom_CM = S0_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      S0_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac506_S0_ang_MC"), h1-> Fill(S0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_S0_ang_CM_MC"), h1-> Fill(S0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac506_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==508 ){
      TLorentzVector S0_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(S0_lmom+pip_lmom+pim_lmom).BoostVector();
      TLorentzVector S0_lmom_CM = S0_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      S0_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac508_S0_ang_MC"), h1-> Fill(S0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac508_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac508_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac508_S0_ang_CM_MC"), h1-> Fill(S0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac508_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac508_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==514 ){
      TLorentzVector S0_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(S0_lmom+pi0_lmom).BoostVector();
      TLorentzVector S0_lmom_CM = S0_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      S0_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM);

      h1= (TH1F*)f-> Get("Reac514_S0_ang_MC"), h1-> Fill(S0_lmom.Vect().CosTheta());
      h1= (TH1F*)f-> Get("Reac514_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1= (TH1F*)f-> Get("Reac514_S0_ang_CM_MC"), h1-> Fill(S0_lmom_CM.Vect().CosTheta());
      h1= (TH1F*)f-> Get("Reac514_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==516 ){
      TLorentzVector S0_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom2 = reacData->GetParticle(2);
      TVector3 LabToCM = -(S0_lmom+pi0_lmom+pi0_lmom2).BoostVector();
      TLorentzVector S0_lmom_CM = S0_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pi0_lmom2_CM = pi0_lmom2;
      S0_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pi0_lmom2_CM.Boost(LabToCM);

      h1= (TH1F*)f-> Get("Reac516_S0_ang_MC"), h1-> Fill(S0_lmom.Vect().CosTheta());
      h1= (TH1F*)f-> Get("Reac516_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta()), h1-> Fill(pi0_lmom2.Vect().CosTheta());
      h1= (TH1F*)f-> Get("Reac516_S0_ang_CM_MC"), h1-> Fill(S0_lmom_CM.Vect().CosTheta());
      h1= (TH1F*)f-> Get("Reac516_pi0_ang_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta()), h1-> Fill(pi0_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==536 ){
      if( reacData->ParticleSize()>2 ){
	TLorentzVector n_spec_lmom = reacData->GetParticle(2);
	h1 = (TH1F*)f-> Get("Reac536_n_spec_mom_MC"), h1-> Fill(n_spec_lmom.Vect().Mag());
      }

      TLorentzVector Sm_lmom = reacData->GetParticle(1);
      TLorentzVector pip_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(Sm_lmom+pip_lmom).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      Sm_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac536_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac536_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());

      bool pim_hit = false;
      bool pip_hit = false;
      bool n_hit = false;
      for( int i=0; i<detData->detectorHitSize(); i++ ){
	DetectorHit *hit = detData-> detectorHit(i);
	if( hit-> detectorID()==CID_CDH ){
	  if( hit-> pdg()== 211 ) pip_hit=true;
	  if( hit-> pdg()==-211 ) pim_hit=true;
	}
	if( hit-> detectorID()==CID_NC && hit->pdg()==2112 ) n_hit = true;
      }
      if( n_hit ){
	h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC_n_hit"), h1-> Fill(Sm_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC_n_hit"), h1-> Fill(pip_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC_n_hit"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC_n_hit"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      }

      if( pim_hit && pip_hit && n_hit ){
	h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC_npipi_hit"), h1-> Fill(Sm_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC_npipi_hit"), h1-> Fill(pip_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC_npipi_hit"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC_npipi_hit"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      }

      int sigmaID=-1;
      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData->track(i);
	if( track->pdgID()==3112 && track->parentTrackID()==0 ) sigmaID=track->trackID();
      }
      if( sigmaID>0 ){
	Track *n_track=0;
	for( int i=0; i<mcData->trackSize(); i++ ){
	  Track *track = mcData->track(i);
	  if( track->pdgID()==2112 && track->parentTrackID()==sigmaID ) n_track=track;
	}
	if( n_track ){
	  bool n_scat_flag = false;
	  if( n_track->processID()==111 ){
	    for( int i=0; i<mcData->trackSize(); i++ ){
	      Track *track = mcData-> track(i);
	      if( track->parentTrackID()==n_track->trackID() && track->pdgID()==2112 ){
		TVector3 vtx = 0.1*track->vertex();
		if( -58.<vtx.Z() && vtx.Z()<58. && vtx.Perp()<15. ) n_scat_flag=true;
	      }
	    }
	  }

	  h1 = (TH1F*)f-> Get("Reac536_decay_n_ang_MC"), h1-> Fill(n_track->momentum().CosTheta());
	  h1 = (TH1F*)f-> Get("Reac536_decay_n_mom_MC"), h1-> Fill(0.001*n_track->momentum().Mag());
	  if( n_scat_flag ){
	    h1 = (TH1F*)f-> Get("Reac536_decay_n_ang_MC_n_scat"), h1-> Fill(n_track->momentum().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac536_decay_n_mom_MC_n_scat"), h1-> Fill(0.001*n_track->momentum().Mag());

	    h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC_n_scat"), h1-> Fill(Sm_lmom.Vect().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC_n_scat"), h1-> Fill(pip_lmom.Vect().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac536_Sm_ang_MC_n_scat"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac536_pip_ang_CM_MC_n_scat"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
	  }
	  if( n_hit ){
	    h1 = (TH1F*)f-> Get("Reac536_decay_n_ang_MC_n_hit"), h1-> Fill(n_track->momentum().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac536_decay_n_mom_MC_n_hit"), h1-> Fill(0.001*n_track->momentum().Mag());
	  }
	  if( pim_hit && pip_hit && n_hit ){
	    h1 = (TH1F*)f-> Get("Reac536_decay_n_ang_MC_npipi_hit"), h1-> Fill(n_track->momentum().CosTheta());
	    h1 = (TH1F*)f-> Get("Reac536_decay_n_mom_MC_npipi_hit"), h1-> Fill(0.001*n_track->momentum().Mag());
	  }
	}
      }
    }
    else if( reactionID==537 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(Sm_lmom+pip_lmom).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      Sm_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac537_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac537_pip_ang_CM_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac537_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac537_Sm_ang_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac537_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac537_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==541 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(2);
      TLorentzVector pi0_lmom2 = reacData->GetParticle(3);
      TVector3 LabToCM = -(Sm_lmom+pip_lmom+pi0_lmom+pi0_lmom2).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pi0_lmom2_CM = pi0_lmom2;
      Sm_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pi0_lmom2_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac541_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac541_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac541_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta()), h1-> Fill(pi0_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac541_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac541_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac541_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta()), h1-> Fill(pi0_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==556 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pip_lmom2 = reacData->GetParticle(2);
      TLorentzVector pim_lmom = reacData->GetParticle(3);
      TVector3 LabToCM = -(Sm_lmom+pip_lmom+pip_lmom2+pim_lmom).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pip_lmom2_CM = pip_lmom2;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sm_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pip_lmom2_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac556_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac556_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac556_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta()), h1-> Fill(pip_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac556_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac556_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac556_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta()), h1-> Fill(pip_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==579 ){
      TLorentzVector Sp1385_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(Sp1385_lmom+pim_lmom).BoostVector();
      TLorentzVector Sp1385_lmom_CM = Sp1385_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sp1385_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac579_S1385p_ang_MC"), h1-> Fill(Sp1385_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac579_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac579_S1385p_ang_CM_MC"), h1-> Fill(Sp1385_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac579_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==592 ){
      TLorentzVector S1385_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(S1385_lmom+pi0_lmom).BoostVector();
      TLorentzVector S1385_lmom_CM = S1385_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      S1385_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac592_S13850_ang_MC"), h1-> Fill(S1385_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac592_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac592_S13850_ang_CM_MC"), h1-> Fill(S1385_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac592_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==597 ){
      TLorentzVector Sm1385_lmom = reacData->GetParticle(1);
      TLorentzVector pip_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(Sm1385_lmom+pip_lmom).BoostVector();
      TLorentzVector Sm1385_lmom_CM = Sm1385_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      Sm1385_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac597_S1385m_ang_MC"), h1-> Fill(Sm1385_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac597_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac597_S1385m_ang_CM_MC"), h1-> Fill(Sm1385_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac597_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==635 ){
      TLorentzVector p_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TLorentzVector km_lmom = reacData->GetParticle(3);
      TVector3 LabToCM = -(p_lmom+pip_lmom+pim_lmom+km_lmom).BoostVector();
      TLorentzVector p_lmom_CM = p_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector km_lmom_CM = km_lmom;
      p_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), km_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac635_p_ang_MC"), h1->Fill(p_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_pip_ang_MC"), h1->Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_pim_ang_MC"), h1->Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_km_ang_MC"), h1->Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_p_ang_CM_MC"), h1->Fill(p_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_pip_ang_CM_MC"), h1->Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_pim_ang_CM_MC"), h1->Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac635_km_ang_CM_MC"), h1->Fill(km_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==670 ){
      TLorentzVector p_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector km_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(p_lmom+pi0_lmom+km_lmom).BoostVector();
      TLorentzVector p_lmom_CM = p_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector km_lmom_CM = km_lmom;
      p_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), km_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac670_p_ang_MC"), h1->Fill(p_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac670_pi0_ang_MC"), h1->Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac670_km_ang_MC"), h1->Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac670_p_ang_CM_MC"), h1->Fill(p_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac670_pi0_ang_CM_MC"), h1->Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac670_km_ang_CM_MC"), h1->Fill(km_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==679 ){
      TLorentzVector p_lmom = reacData->GetParticle(0);
      TLorentzVector pim_lmom = reacData->GetParticle(1);
      TLorentzVector k0_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(p_lmom+pim_lmom+k0_lmom).BoostVector();
      TLorentzVector p_lmom_CM = p_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector k0_lmom_CM = k0_lmom;
      p_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), k0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac679_p_ang_MC"), h1->Fill(p_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac679_pim_ang_MC"), h1->Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac679_k0_ang_MC"), h1->Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac679_p_ang_CM_MC"), h1->Fill(p_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac679_pim_ang_CM_MC"), h1->Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac679_k0_ang_CM_MC"), h1->Fill(k0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==707 ){
      TLorentzVector kms_lmom = reacData-> GetParticle(0);
      TLorentzVector p_lmom = reacData-> GetParticle(1);
      TVector3 LabToCM = -(kms_lmom+p_lmom).BoostVector();
      TLorentzVector kms_lmom_CM = kms_lmom;
      TLorentzVector p_lmom_CM = p_lmom;
      kms_lmom_CM.Boost(LabToCM), p_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac707_ksm_ang_MC"), h1-> Fill(kms_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac707_p_ang_MC"), h1-> Fill(p_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac707_ksm_ang_CM_MC"), h1-> Fill(kms_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac707_p_ang_CM_MC"), h1-> Fill(p_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==714 ){
      TLorentzVector n_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TLorentzVector k0_lmom = reacData->GetParticle(3);
      TVector3 LabToCM = -(n_lmom+pip_lmom+pim_lmom+k0_lmom).BoostVector();
      TLorentzVector n_lmom_CM = n_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector k0_lmom_CM = k0_lmom;
      n_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), k0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac714_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_k0_ang_MC"), h1-> Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac714_k0_ang_CM_MC"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==716 ){
      TLorentzVector n_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector km_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(n_lmom+pip_lmom+km_lmom).BoostVector();
      TLorentzVector n_lmom_CM = n_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector km_lmom_CM = km_lmom;
      n_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), km_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac716_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac716_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac716_km_ang_MC"), h1-> Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac716_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac716_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac716_km_ang_CM_MC"), h1-> Fill(km_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==729 ){
      TLorentzVector n_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector k0_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(n_lmom+pi0_lmom+k0_lmom).BoostVector();
      TLorentzVector n_lmom_CM = n_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector k0_lmom_CM = k0_lmom;
      n_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), k0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac729_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac729_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac729_k0_ang_MC"), h1-> Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac729_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac729_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac729_k0_ang_CM_MC"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==735 ){
      if( reacData->ParticleSize()>2 ){
	TLorentzVector n_spec_lmom = reacData->GetParticle(2);
	h1 = (TH1F*)f-> Get("Reac735_n_spec_mom_MC"), h1-> Fill(0.001*n_spec_lmom.Vect().Mag());
      }
      TLorentzVector n_lmom = reacData->GetParticle(1);
      TLorentzVector k0_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(n_lmom+k0_lmom).BoostVector();
      TLorentzVector n_lmom_CM = n_lmom;
      TLorentzVector k0_lmom_CM = k0_lmom;
      n_lmom_CM.Boost(LabToCM), k0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac735_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac735_n_mom_MC"), h1-> Fill(0.001*n_lmom.Vect().Mag());
      h1 = (TH1F*)f-> Get("Reac735_k0_ang_MC"), h1-> Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac735_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac735_k0_ang_CM_MC"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());

      if( eff_NC ){
	h1 = (TH1F*)f-> Get("Reac735_k0_ang_MC_trig"), h1-> Fill(k0_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_k0_ang_CM_MC_trig"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
      }
      bool pim_hit = false;
      bool pip_hit = false;
      bool n_hit = false;
      for( int i=0; i<detData->detectorHitSize(); i++ ){
	DetectorHit *hit = detData-> detectorHit(i);
	if( hit-> detectorID()==CID_CDH ){
	  if( hit-> pdg()== 211 ) pip_hit=true;
	  if( hit-> pdg()==-211 ) pim_hit=true;
	}
	if( hit-> detectorID()==CID_NC && hit->pdg()==2112 ) n_hit = true;
      }


      if( n_hit ){
	h1 = (TH1F*)f-> Get("Reac735_n_ang_MC_n_hit"), h1-> Fill(n_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_n_mom_MC_n_hit"), h1-> Fill(0.001*n_lmom.Vect().Mag());
	h1 = (TH1F*)f-> Get("Reac735_k0_ang_MC_n_hit"), h1-> Fill(k0_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_n_ang_CM_MC_n_hit"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_k0_ang_CM_MC_n_hit"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
      }
      if( pim_hit && pip_hit && n_hit ){
	h1 = (TH1F*)f-> Get("Reac735_n_ang_MC_npipi_hit"), h1-> Fill(n_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_n_mom_MC_npipi_hit"), h1-> Fill(0.001*n_lmom.Vect().Mag());
	h1 = (TH1F*)f-> Get("Reac735_k0_ang_MC_npipi_hit"), h1-> Fill(k0_lmom.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_n_ang_CM_MC_npipi_hit"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
	h1 = (TH1F*)f-> Get("Reac735_k0_ang_CM_MC_npipi_hit"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
      }
    }
    else if( reactionID==737 ){
      TLorentzVector n_lmom = reacData->GetParticle(1);
      TLorentzVector ks0_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(n_lmom+ks0_lmom).BoostVector();
      TLorentzVector n_lmom_CM = n_lmom;
      TLorentzVector ks0_lmom_CM = ks0_lmom;
      n_lmom_CM.Boost(LabToCM), ks0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac737_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac737_ks0_ang_MC"), h1-> Fill(ks0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac737_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac737_ks0_ang_CM_MC"), h1-> Fill(ks0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==767 ){
      TLorentzVector km_lmom = reacData->GetParticle(0);
      TLorentzVector delta_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(km_lmom+delta_lmom).BoostVector();
      TLorentzVector km_lmom_CM = km_lmom;
      TLorentzVector delta_lmom_CM = delta_lmom;
      km_lmom.Boost(LabToCM), delta_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac767_km_ang_MC"), h1-> Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac767_Deltap_ang_MC"), h1-> Fill(delta_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac767_km_ang_CM_MC"), h1-> Fill(km_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac767_Deltap_ang_CM_MC"), h1-> Fill(delta_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==771 ){
      TLorentzVector k0_lmom = reacData->GetParticle(0);
      TLorentzVector delta_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(k0_lmom+delta_lmom).BoostVector();
      TLorentzVector k0_lmom_CM = k0_lmom;
      TLorentzVector delta_lmom_CM = delta_lmom;
      k0_lmom.Boost(LabToCM), delta_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac771_k0_ang_MC"), h1-> Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac771_Delta0_ang_MC"), h1-> Fill(delta_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac771_k0_ang_CM_MC"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac771_Delta0_ang_CM_MC"), h1-> Fill(delta_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==821 ){
      TLorentzVector km_lmom = reacData->GetParticle(0);
      TLorentzVector n_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(km_lmom+n_lmom).BoostVector();
      TLorentzVector km_lmom_CM = km_lmom;
      TLorentzVector n_lmom_CM = n_lmom;
      km_lmom_CM.Boost(LabToCM), n_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac821_km_ang_MC"), h1-> Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac821_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac821_km_ang_CM_MC"), h1-> Fill(km_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac821_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());

      h1 = (TH1F*)f-> Get("Reac821_km_mom_MC"), h1-> Fill(0.001*km_lmom.Vect().Mag());
      Track *sm_track = 0;
      Track *sp_track = 0;
      Track *k0_track = 0;

      for( int i=0; i<mcData->trackSize(); i++ ){
	Track *track = mcData-> track(i);
	if( track->pdgID()==3112 ) sm_track=track;
	if( track->pdgID()==3222 ) sp_track=track;
	if( track->pdgID()==311  ) k0_track=track;
      }
      if( sm_track ){
	h1 = (TH1F*)f-> Get("Reac821_km_mom_MC_Sm"), h1-> Fill(0.001*km_lmom.Vect().Mag());
	printReac(sm_track->parentTrackID());
      }
      if( sp_track ){
	h1 = (TH1F*)f-> Get("Reac821_km_mom_MC_Sp"), h1-> Fill(0.001*km_lmom.Vect().Mag());
	printReac(sp_track->parentTrackID());
      }
      if( k0_track ){
	h1 = (TH1F*)f-> Get("Reac821_km_mom_MC_K0"), h1-> Fill(0.001*km_lmom.Vect().Mag());
	printReac(k0_track->parentTrackID());
      }
      fillHist_knEl();
    }
    else if( reactionID==851 ){
      TLorentzVector L_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TLorentzVector pim_lmom2 = reacData->GetParticle(3);
      TVector3 LabToCM = -(L_lmom+pip_lmom+pim_lmom+pim_lmom2).BoostVector();
      TLorentzVector L_lmom_CM = L_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector pim_lmom2_CM = pim_lmom2;
      L_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), pim_lmom2_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac851_L_ang_MC"), h1-> Fill(L_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac851_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac851_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta()),  h1-> Fill(pim_lmom2.Vect().CosTheta());

      h1 = (TH1F*)f-> Get("Reac851_L_ang_CM_MC"), h1-> Fill(L_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac851_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac851_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta()),  h1-> Fill(pim_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==856 ){
      TLorentzVector L_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(L_lmom+pi0_lmom+pim_lmom).BoostVector();
      TLorentzVector L_lmom_CM = L_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      L_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac856_L_ang_MC"), h1-> Fill(L_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac856_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac856_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac856_L_ang_CM_MC"), h1-> Fill(L_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac856_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac856_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==862 ){
      TLorentzVector L_lmom = reacData->GetParticle(0);
      TLorentzVector pim_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(L_lmom+pim_lmom).BoostVector();
      TLorentzVector L_lmom_CM = L_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      L_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac862_L_ang_MC"), h1-> Fill(L_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac862_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac862_L_ang_CM_MC"), h1-> Fill(L_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac862_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==893 ){
      TLorentzVector pim_lmom = reacData->GetParticle(0);
      TLorentzVector L1405_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(pim_lmom+L1405_lmom).BoostVector();
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector L1405_lmom_CM = L1405_lmom;
      pim_lmom_CM.Boost(LabToCM), L1405_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac893_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac893_L1405_ang_MC"), h1-> Fill(L1405_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac893_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac893_L1405_ang_CM_MC"), h1-> Fill(L1405_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==897 ){
      TLorentzVector pim_lmom = reacData->GetParticle(0);
      TLorentzVector L1520_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(pim_lmom+L1520_lmom).BoostVector();
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector L1520_lmom_CM = L1520_lmom;
      pim_lmom_CM.Boost(LabToCM), L1520_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac897_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac897_L1520_ang_MC"), h1-> Fill(L1520_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac897_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac897_L1520_ang_CM_MC"), h1-> Fill(L1520_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==912 ){
      TLorentzVector Sp_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TLorentzVector pim_lmom2 = reacData->GetParticle(3);
      TVector3 LabToCM = -(Sp_lmom+pi0_lmom+pim_lmom+pim_lmom2).BoostVector();
      TLorentzVector Sp_lmom_CM = Sp_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector pim_lmom2_CM = pim_lmom2;
      Sp_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), pim_lmom2_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac912_Sp_ang_MC"), h1-> Fill(Sp_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta()), h1-> Fill(pim_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_Sp_ang_CM_MC"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta()), h1-> Fill(pim_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==915 ){
      TLorentzVector Sp_lmom = reacData->GetParticle(0);
      TLorentzVector pim_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom2 = reacData->GetParticle(2);
      TVector3 LabToCM = -(Sp_lmom+pim_lmom+pim_lmom2).BoostVector();
      TLorentzVector Sp_lmom_CM = Sp_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector pim_lmom2_CM = pim_lmom2;
      Sp_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), pim_lmom2_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac915_Sp_ang_MC"), h1-> Fill(Sp_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac915_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta()), h1-> Fill(pim_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac915_Sp_ang_CM_MC"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac915_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta()), h1-> Fill(pim_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==922 ){
      TLorentzVector Sp_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(Sp_lmom+pi0_lmom+pim_lmom).BoostVector();
      TLorentzVector Sp_lmom_CM = Sp_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sp_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac912_Sp_ang_MC"), h1-> Fill(Sp_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_Sp_ang_CM_MC"), h1-> Fill(Sp_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac912_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==923 ){
      TLorentzVector S0_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(S0_lmom+pim_lmom).BoostVector();
      TLorentzVector S0_lmom_CM = S0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      S0_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac923_S0_ang_MC"), h1-> Fill(S0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac923_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac923_S0_ang_CM_MC"), h1-> Fill(S0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac923_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==929 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(2);
      TLorentzVector pim_lmom = reacData->GetParticle(3);
      TVector3 LabToCM = -(Sm_lmom+pip_lmom+pi0_lmom+pim_lmom).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sm_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac929_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac929_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==931 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(0);
      TLorentzVector pip_lmom = reacData->GetParticle(1);
      TLorentzVector pim_lmom = reacData->GetParticle(2);
      TVector3 LabToCM = -(Sm_lmom+pip_lmom+pim_lmom).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pip_lmom_CM = pip_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      Sm_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), pip_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac931_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac931_pip_ang_MC"), h1-> Fill(pip_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac931_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac931_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac931_pip_ang_CM_MC"), h1-> Fill(pip_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac931_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==938 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom = reacData->GetParticle(0);
      TVector3 LabToCM = -(Sm_lmom+pi0_lmom).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      Sm_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac938_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac938_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac938_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac938_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==946 ){
      TLorentzVector Sm_lmom = reacData->GetParticle(0);
      TLorentzVector pi0_lmom = reacData->GetParticle(1);
      TLorentzVector pi0_lmom2 = reacData->GetParticle(2);
      TVector3 LabToCM = -(Sm_lmom+pi0_lmom+pi0_lmom2).BoostVector();
      TLorentzVector Sm_lmom_CM = Sm_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      TLorentzVector pi0_lmom2_CM = pi0_lmom2;
      Sm_lmom_CM.Boost(LabToCM), pi0_lmom_CM.Boost(LabToCM), pi0_lmom2_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac946_Sm_ang_MC"), h1-> Fill(Sm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac946_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta()), h1-> Fill(pi0_lmom2.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac946_Sm_ang_CM_MC"), h1-> Fill(Sm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac946_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta()), h1-> Fill(pi0_lmom2_CM.Vect().CosTheta());
    }
    else if( reactionID==984 ){
      TLorentzVector S1385_lmom = reacData-> GetParticle(1);
      TLorentzVector pim_lmom = reacData-> GetParticle(0);
      TVector3 LabToCM = -(S1385_lmom+pim_lmom).BoostVector();
      TLorentzVector S1385_lmom_CM = S1385_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      S1385_lmom_CM.Boost(LabToCM), pim_lmom.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac984_S1385_ang_MC"), h1-> Fill(S1385_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac984_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac984_S1385_ang_CM_MC"), h1-> Fill(S1385_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac984_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==990 ){
      TLorentzVector S1385m_lmom = reacData-> GetParticle(1);
      TLorentzVector pi0_lmom = reacData-> GetParticle(0);
      TVector3 LabToCM = -(S1385m_lmom+pi0_lmom).BoostVector();
      TLorentzVector S1385m_lmom_CM = S1385m_lmom;
      TLorentzVector pi0_lmom_CM = pi0_lmom;
      S1385m_lmom_CM.Boost(LabToCM), pi0_lmom.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac990_S1385m_ang_MC"), h1-> Fill(S1385m_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac990_pi0_ang_MC"), h1-> Fill(pi0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac990_S1385m_ang_CM_MC"), h1-> Fill(S1385m_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac990_pi0_ang_CM_MC"), h1-> Fill(pi0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==1025 ){
      TLorentzVector p_lmom = reacData-> GetParticle(0);
      TLorentzVector pim_lmom = reacData-> GetParticle(1);
      TLorentzVector km_lmom = reacData-> GetParticle(2);
      TVector3 LabToCM = -(p_lmom+pim_lmom+km_lmom).BoostVector();
      TLorentzVector p_lmom_CM = p_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector km_lmom_CM = km_lmom;
      p_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), km_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac1025_p_ang_MC"), h1-> Fill(p_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1025_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1025_km_ang_MC"), h1-> Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1025_p_ang_CM_MC"), h1-> Fill(p_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1025_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1025_km_ang_CM_MC"), h1-> Fill(km_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==1041 ){
      TLorentzVector n_lmom = reacData-> GetParticle(0);
      TLorentzVector pim_lmom = reacData-> GetParticle(1);
      TLorentzVector k0_lmom = reacData-> GetParticle(2);
      TVector3 LabToCM = -(n_lmom+pim_lmom+k0_lmom).BoostVector();
      TLorentzVector n_lmom_CM = n_lmom;
      TLorentzVector pim_lmom_CM = pim_lmom;
      TLorentzVector k0_lmom_CM = k0_lmom;
      n_lmom_CM.Boost(LabToCM), pim_lmom_CM.Boost(LabToCM), k0_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac1041_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1041_pim_ang_MC"), h1-> Fill(pim_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1041_k0_ang_MC"), h1-> Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1041_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1041_pim_ang_CM_MC"), h1-> Fill(pim_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1041_k0_ang_CM_MC"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==1045 ){
      TLorentzVector ksm_lmom = reacData->GetParticle(0);
      TLorentzVector n_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(ksm_lmom+n_lmom).BoostVector();
      TLorentzVector ksm_lmom_CM = ksm_lmom;
      TLorentzVector n_lmom_CM = n_lmom;
      ksm_lmom_CM.Boost(LabToCM), n_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac1045_ksm_ang_MC"), h1-> Fill(ksm_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1045_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1045_ksm_ang_CM_MC"), h1-> Fill(ksm_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1045_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==1055 ){
      TLorentzVector km_lmom = reacData->GetParticle(0);
      TLorentzVector D_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(km_lmom+D_lmom).BoostVector();
      TLorentzVector km_lmom_CM = km_lmom;
      TLorentzVector D_lmom_CM = D_lmom;
      km_lmom_CM.Boost(LabToCM), D_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac1055_km_ang_MC"), h1-> Fill(km_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1055_D0_ang_MC"), h1-> Fill(D_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1055_km_ang_CM_MC"), h1-> Fill(km_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1055_D0_ang_CM_MC"), h1-> Fill(D_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==1061 ){
      TLorentzVector k0_lmom = reacData->GetParticle(0);
      TLorentzVector D_lmom = reacData->GetParticle(1);
      TVector3 LabToCM = -(k0_lmom+D_lmom).BoostVector();
      TLorentzVector k0_lmom_CM = k0_lmom;
      TLorentzVector D_lmom_CM = D_lmom;
      k0_lmom_CM.Boost(LabToCM), D_lmom_CM.Boost(LabToCM);

      h1 = (TH1F*)f-> Get("Reac1061_k0_ang_MC"), h1-> Fill(k0_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1061_Dm_ang_MC"), h1-> Fill(D_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1061_k0_ang_CM_MC"), h1-> Fill(k0_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1061_Dm_ang_CM_MC"), h1-> Fill(D_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==1098 ){
      TLorentzVector n_lmom = reacData->GetParticle(0);
      TLorentzVector n_lmom_CM = reacData->GetCMParticle(0);
      TLorentzVector l1405_lmom = reacData->GetParticle(1);
      TLorentzVector l1405_lmom_CM = reacData->GetCMParticle(1);

      h1 = (TH1F*)f-> Get("Reac1098_l1405_mass_MC"), h1-> Fill(0.001*l1405_lmom.M());
      h1 = (TH1F*)f-> Get("Reac1098_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1098_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1098_l1405_ang_MC"), h1-> Fill(l1405_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac1098_l1405_ang_CM_MC"), h1-> Fill(l1405_lmom_CM.Vect().CosTheta());
    }
    else if( reactionID==3098 ){
      TLorentzVector n_lmom = reacData->GetParticle(0);
      TLorentzVector n_lmom_CM = reacData->GetCMParticle(0);
      TLorentzVector l1405_lmom = reacData->GetParticle(1);
      TLorentzVector l1405_lmom_CM = reacData->GetCMParticle(1);

      h1 = (TH1F*)f-> Get("Reac3098_l1405_mass_MC"), h1-> Fill(0.001*l1405_lmom.M());
      h1 = (TH1F*)f-> Get("Reac3098_n_ang_MC"), h1-> Fill(n_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac3098_n_ang_CM_MC"), h1-> Fill(n_lmom_CM.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac3098_l1405_ang_MC"), h1-> Fill(l1405_lmom.Vect().CosTheta());
      h1 = (TH1F*)f-> Get("Reac3098_l1405_ang_CM_MC"), h1-> Fill(l1405_lmom_CM.Vect().CosTheta());
    }
  }
#endif
}

std::string SimDataReader::parName(const int &pdg)
{
  if( pdg==11 ) return "e-";
  else if( pdg==22    ) return "gamma";
  else if( pdg==12    ) return "nu_e";
  else if( pdg==13    ) return "mu-";
  else if( pdg==14    ) return "nu_mu";
  else if( pdg==15    ) return "tau-";
  else if( pdg==16    ) return "nu_tau";
  else if( pdg==-11   ) return "e+";
  else if( pdg==-12   ) return "anti_nu_e";
  else if( pdg==-13   ) return "mu+";
  else if( pdg==-14   ) return "anti_nu_mu";
  else if( pdg==-15   ) return "tau+";
  else if( pdg==-16   ) return "anti_nu_tau";
  else if( pdg==323   ) return "k_star+";
  else if( pdg==313   ) return "k_star0";
  else if( pdg==321   ) return "kaon+";
  else if( pdg==311   ) return "kaon0";
  else if( pdg==-323  ) return "k_star-";
  else if( pdg==-313  ) return "anti_k_star0";
  else if( pdg==-321  ) return "kaon-";
  else if( pdg==-311  ) return "anti_kaon0";
  else if( pdg==130   ) return "kaon0L";
  else if( pdg==310   ) return "kaon0S";
  else if( pdg==211   ) return "pi+";
  else if( pdg==-211  ) return "pi-";
  else if( pdg==111   ) return "pi0";
  else if( pdg==113   ) return "rho0";
  else if( pdg==213   ) return "rho+";
  else if( pdg==-213  ) return "rho-";
  else if( pdg==2214  ) return "delta+";
  else if( pdg==1114  ) return "delta-";
  else if( pdg==2114  ) return "delta0";
  else if( pdg==3122  ) return "lambda";
  else if( pdg==13122 ) return "lambda(1405)";
  else if( pdg==3124  ) return "lambda(1520)";
  else if( pdg==23122 ) return "lambda(1600)";
  else if( pdg==33122 ) return "lambda(1670)";
  else if( pdg==13124 ) return "lambda(1690)";
  else if( pdg==43122 ) return "lambda(1800)";
  else if( pdg==53122 ) return "lambda(1810)";
  else if( pdg==3126  ) return "lambda(1820)";
  else if( pdg==13126 ) return "lambda(1830)";
  else if( pdg==23124 ) return "lambda(1890)";
  else if( pdg==2112  ) return "neutron";
  else if( pdg==2212  ) return "proton";
  else if( pdg==3224  ) return "sigma(1385)+";
  else if( pdg==3114  ) return "sigma(1385)-";
  else if( pdg==3214  ) return "sigma(1385)0";
  else if( pdg==3222  ) return "sigma+";
  else if( pdg==3112  ) return "sigma-";
  else if( pdg==3212  ) return "sigma0";
  else if( pdg==1000010020 ) return "deuteron";
  else if( pdg==1000010030 ) return "triton";
  else if( pdg==1000020030 ) return "He3";
  else if( pdg==1000020040 ) return "alpha";
  else if( pdg==1000060120 ) return "C12";
  else if( pdg==1000070140 ) return "N14";
  else if( pdg==1000080160 ) return "O16";
  else if( pdg==1000140280 ) return "Si28";
  else if( pdg==1000160330 ) return "S33";
  else if( pdg==1000200410 ) return "Ca41";
  else if( pdg==1000210480 ) return "Sc48";
  else if( pdg==1000220460 ) return "Ti46";
  else if( pdg==1000220470 ) return "Ti47";
  else if( pdg==1000240480 ) return "Cr48";
  else if( pdg==1000240490 ) return "Cr49";
  else if( pdg==1000240500 ) return "Cr50";
  else if( pdg==1000240510 ) return "Cr51";
  else if( pdg==1000240520 ) return "Cr52";
  else if( pdg==1000240530 ) return "Cr53";
  else if( pdg==1000240540 ) return "Cr54";
  else if( pdg==1000250530 ) return "Mn53";
  else if( pdg==1000250540 ) return "Mn54";
  else if( pdg==1000260520 ) return "Fe52";
  else if( pdg==1000260530 ) return "Fe53";
  else if( pdg==1000260540 ) return "Fe54";
  else if( pdg==1000260550 ) return "Fe55";
  else if( pdg==1000260560 ) return "Fe56";
  else if( pdg==1000260570 ) return "Fe57";
  else if( pdg==1000260580 ) return "Fe58";
  else return Form("%d", pdg);
}

bool SimDataReader::isNC()
{
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData-> detectorHit(i);
    if( hit->detectorID()==CID_NC ){
      if( hit->pdg()==22 || hit->pdg()==2112 ) return true;
    }
  }
  return false;
}

bool SimDataReader::isCDH2()
{
  int nCDH=0;
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit = detData-> detectorHit(i);
    if( hit->detectorID()==CID_CDH  ) nCDH++;
  }
  if( nCDH>1 ) return true;
  else return false;
}

void SimDataReader::fillMM_N(TFile *f)
{
  //  std::cout<<"fillMM_N"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_mmN");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillMM_N_woK0(TFile *f)
{
  //  std::cout<<"fillMM_N w/o"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_mmN_woK0");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillMM_N_woSf(TFile *f)
{
  //  std::cout<<"fillMM_N w/o"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_mmN_woSf");
    h1-> Fill(fL1405_lmom.M());
  }
}


void SimDataReader::fillMM_N_wo(TFile *f)
{
  //  std::cout<<"fillMM_N w/o"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_mmN_woAll");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillLpim(TFile *f)
{
  //  std::cout<<"fillMM_N"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_Lpim");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillLpim_p(TFile *f)
{
  //  std::cout<<"fillMM_N"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_Lpim_p");
    h1-> Fill(fL1405_lmom.M());
    bool PCCVC_flag=false;
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit=detData->detectorHit(i);
      if( hit->detectorID()==CID_PC || hit->detectorID()==CID_CVC ){
	if( hit->pdg()==2212 && hit->parentID()==0 ){
	  PCCVC_flag=true;
	}
      }
    }
    if( PCCVC_flag ){
      TH1F *h1 = (TH1F*)f-> Get("Eff_Lpim_p_hit");
      h1-> Fill(fL1405_lmom.M());
    }
  }
}

void SimDataReader::fillCharged(TFile *f)
{
  //  std::cout<<"fillCharged"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_Charged");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillpipi(TFile *f)
{
  //  std::cout<<"fillCharged"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_pipi");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fill_woK0(TFile *f)
{
  //  std::cout<<"fillCharged"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_pipi_woK0");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fill_woSf(TFile *f)
{
  //  std::cout<<"fillCharged"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_pipi_woSf");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillSp(TFile *f)
{
  //  std::cout<<"fillSp"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_Sp");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillSm(TFile *f)
{
  //  std::cout<<"fillSm"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_Sm");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillMM_N_woRcS(TFile *f)
{
  //  std::cout<<"fillSm"<<std::endl;
  if( fL1405_lmom.M()>1.0 ){
    TH1F *h1 = (TH1F*)f-> Get("Eff_mmN_woRcS");
    h1-> Fill(fL1405_lmom.M());
  }
}

void SimDataReader::fillppim_MC(TFile *f)
{
  if( fL1405_lmom.M()>1.0 ){
    int BPD_p_hit=0;
    int BPD_pim_hit=0;
    int BPD_pip_hit=0;
    int BPD_beam_hit=0;
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit = detData->detectorHit(i);
      if( hit-> detectorID()==CID_BPD ){
	if( hit-> pdg()==2212 ) BPD_p_hit++;
	if( hit-> pdg()==-211 ) BPD_pim_hit++;
	if( hit-> pdg()== 211 ) BPD_pip_hit++;
	if( hit-> pdg()== 321 && hit-> parentID()==0 ) BPD_beam_hit++;
      }
      if( BPD_p_hit>0 ){
	TH1F *h1 = (TH1F*)f-> Get("Eff_ppim_hit_MC");
	h1-> Fill(fL1405_lmom.M());
      }
    }
  }
}

void SimDataReader::fillKNpipiMM(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KNpipi_MM_Reac%d", reactionID)), h1-> Fill(mm);
}

void SimDataReader::fillKNpipiMM_woAll(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KNpipi_MM_woAll_Reac%d", reactionID)), h1-> Fill(mm);
}

void SimDataReader::fillKNpipiMM_wK0(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KNpipi_MM_wK0_Reac%d", reactionID)), h1-> Fill(mm);
}

void SimDataReader::fillKNpipiMM_wSp(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KNpipi_MM_wSp_Reac%d", reactionID)), h1-> Fill(mm);
}

void SimDataReader::fillKNpipiMM_wSm(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KNpipi_MM_wSm_Reac%d", reactionID)), h1-> Fill(mm);
}

void SimDataReader::fillKN_MM(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KN_MM_Reac%d", reactionID)), h1-> Fill(mm);
}

void SimDataReader::fillKN_MM_wN(TFile *f, const double &mm)
{
  TH1F *h1=0;
  int reactionID = reacData-> ReactionID();
  h1 = (TH1F*)f-> Get(Form("KN_MM_wN_Reac%d", reactionID)), h1-> Fill(mm);
}

std::string SimDataReader::getProcessName(const int &processID)
{
  if( processID==0        ) return "Initial";
  else if( processID==1   ) return "Decay";
  else if( processID==2   ) return "conv";
  else if( processID==3   ) return "Transportation";
  else if( processID==4   ) return "phot";
  else if( processID==5   ) return "annihil";
  else if( processID==6   ) return "compt";
  else if( processID==7   ) return "eBrem";
  else if( processID==8   ) return "hadElastic";
  else if( processID==9   ) return "CoulombScat";
  else if( processID==10  ) return "nKiller";
  else if( processID==11  ) return "photoNuclear";
  else if( processID==100 ) return "pi-Inelastic";
  else if( processID==101 ) return "pi+Inelastic";
  else if( processID==102 ) return "kaon-Inelastic";
  else if( processID==103 ) return "kaon+Inelastic";
  else if( processID==104 ) return "kaon0LInelastic";
  else if( processID==105 ) return "kaon0SInelastic";
  else if( processID==106 ) return "lambdaInelastic";
  else if( processID==107 ) return "sigma+Inelastic";
  else if( processID==108 ) return "sigma-Inelastic";
  else if( processID==109 ) return "sigma0Inelastic";
  else if( processID==110 ) return "protonInelastic";
  else if( processID==111 ) return "neutronInelastic";
  else if( processID==112 ) return "dInelastic";
  else if( processID==113 ) return "tInelastic";
  else if( processID==200 ) return "eIoni";
  else if( processID==201 ) return "hIoni";
  else if( processID==202 ) return "ionIoni";
  else if( processID==203 ) return "muIoni";
  else if( processID==204 ) return "hBertiniCAptureAtRest";
  else if( processID==205 ) return "nCaptrue";
  else if( processID==206 ) return "muMinusCaptrueAtRest";
  else return "unknown";
}

void SimDataReader::printReac(const int &parentID)
{

}

void SimDataReader::printTrackTrace(const int &trackID)
{
  Track *track=getMCTrack(trackID);
  std::string str=">";
  while( track->parentTrackID()!=0 ){
    std::cout<<str<<parName(track->pdgID())<<"  Process : "<<getProcessName(track->processID())<<std::endl;
    str += ">";
    track = getMCTrack(track->parentTrackID());
  }
}
