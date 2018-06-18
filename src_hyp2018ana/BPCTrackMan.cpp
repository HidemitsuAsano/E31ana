#include "BPCTrackMan.h"

BPCTrackMan::BPCTrackMan(ConfMan *conf, BeamLineHitMan *blMan) : fConfMan(conf), fBLMan(blMan)
{
  set();
}

void BPCTrackMan::set()
{
  int nw, xy;
  double z, xy0, dxy, wl, tilt, ra;
  BLDCWireMapMan *wireMap = fConfMan-> GetBLDCWireMapManager();
  int index_x=0;
  int index_y=0;
  for( int lay=1; lay<=8; lay++ ){
    wireMap-> GetParam(CID_BPC, lay, nw, z, xy, xy0, dxy, wl, tilt, ra);
    std::vector<ChamberLikeHit*> fHits;
    for( int i=0; i<fBLMan->nBPC(lay); i++ ){
      if( xy==0 ) fXZHits[index_x].push_back(fBLMan->BPC(lay, i));
      else if( xy==1 ) fYZHits[index_y].push_back(fBLMan->BPC(lay, i));
    }
    if( xy==0 ) index_x++;
    else if( xy==1 ) index_y++;
  }
}

void BPCTrackMan::DoTracking()
{
#if 0
  std::cout<<" X Layer 1 & 2 n Cluseter : "<<nXCluster(0)<<std::endl;
  std::cout<<" X Layer 3 & 4 n Cluseter : "<<nXCluster(1)<<std::endl;

  std::cout<<" Y Layer 1 & 2 n Cluseter : "<<nYCluster(0)<<std::endl;
  std::cout<<" Y Layer 3 & 4 n Cluseter : "<<nYCluster(1)<<std::endl;
#endif
  while( true ){
    XYTrack best_track(0);
    for( int i1=0; i1<nXCluster(0); i1++ ){
      for( int i2=0; i2<nXCluster(1); i2++ ){
	XYTrack track(0);
	setXCluster(0, i1, track);
	setXCluster(1, i2, track);
	if( track.nhit()>2 ){
	  track.calcChi2();
	  //	  track.dump();
	  if( track.chi2()<best_track.chi2() ){
	    best_track = track;
	  }
	}
      }
    }
    if( best_track.chi2()<10000 ){
      fXTracks.push_back(best_track);
      deleteHit(best_track);
    }
    else break;
  }

  while( true ){
    XYTrack best_track(1);
    for( int i1=0; i1<nYCluster(0); i1++ ){
      for( int i2=0; i2<nYCluster(1); i2++ ){
	XYTrack track(1);
	setYCluster(0, i1, track);
	setYCluster(1, i2, track);
	if( track.nhit()>2 ){
	  track.calcChi2();
	  //	  track.dump();
	  if( track.chi2()<best_track.chi2() ){
	    best_track = track;
	  }
	}
      }
    }
    if( best_track.chi2()<10000 ){
      fYTracks.push_back(best_track);
      deleteHit(best_track);
    }
    else break;
  }
}

bool BPCTrackMan::setXCluster(const int &i1, const int &i2, XYTrack &track)
{
  int index=0;
  if( !fXZHits[2*i1].empty() && !fXZHits[2*i1+1].empty() ){
    for( int j1=0; j1<fXZHits[2*i1].size(); j1++ ){
      for( int j2=0; j2<fXZHits[2*i1+1].size(); j2++ ){
	if( fabs(fXZHits[2*i1][j1]->wx()-fXZHits[2*i1+1][j2]->wx())<0.54 ){
	  if( index==i2 ){
	    track.addHit(fXZHits[2*i1][j1]);
	    track.addHit(fXZHits[2*i1+1][j2]);
	    return true;
	  }
	  else index++;
	}
      }
    }
    
    for( int j1=0; j1<fXZHits[2*i1].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fXZHits[2*i1+1].size(); j2++ ){
	if( fabs(fXZHits[2*i1][j1]->wx()-fXZHits[2*i1+1][j2]->wx())<0.54 ) flag=false;
      }
      if( flag ){
	if( index==i2 ){
	  track.addHit(fXZHits[2*i1][j1]);
	  return true;
	}
	else index++;
      }
    }
    
    for( int j1=0; j1<fXZHits[2*i1+1].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fXZHits[2*i1].size(); j2++ ){
	if( fabs(fXZHits[2*i1][j2]->wx()-fXZHits[2*i1+1][j1]->wx())<0.54 ) flag=false;
      }
      if( flag ){
	if( index==i2 ){
	  track.addHit(fXZHits[2*i1+1][j1]);
	  return true;
	}
	index++;
      }
    }
  }
  else if( fXZHits[2*i1].empty() ){
    track.addHit(fXZHits[2*i1+1][index]);
    return true;
  }
  else if( fXZHits[2*i1+1].empty() ){
    track.addHit(fXZHits[2*i1][index]);
    return true;
  }
  return false;
}

bool BPCTrackMan::setYCluster(const int &i1, const int &i2, XYTrack &track)
{
  int index=0;
  if( !fYZHits[2*i1].empty() && !fYZHits[2*i1+1].empty() ){
    for( int j1=0; j1<fYZHits[2*i1].size(); j1++ ){
      for( int j2=0; j2<fYZHits[2*i1+1].size(); j2++ ){
	if( fabs(fYZHits[2*i1][j1]->wx()-fYZHits[2*i1+1][j2]->wx())<0.54 ){
	  if( index==i2 ){
	    track.addHit(fYZHits[2*i1][j1]);
	    track.addHit(fYZHits[2*i1+1][j2]);
	    return true;
	  }
	  else index++;
	}
      }
    }
    
    for( int j1=0; j1<fYZHits[2*i1].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fYZHits[2*i1+1].size(); j2++ ){
	if( fabs(fYZHits[2*i1][j1]->wx()-fYZHits[2*i1+1][j2]->wx())<0.54 ) flag=false;
      }
      if( flag ){
	if( index==i2 ){
	  track.addHit(fYZHits[2*i1][j1]);
	  return true;
	}
	else index++;
      }
    }
    
    for( int j1=0; j1<fYZHits[2*i1+1].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fYZHits[2*i1].size(); j2++ ){
	if( fabs(fYZHits[2*i1][j2]->wx()-fYZHits[2*i1+1][j1]->wx())<0.54 ) flag=false;
      }
      if( flag ){
	if( index==i2 ){
	  track.addHit(fYZHits[2*i1+1][j1]);
	  return true;
	}
	index++;
      }
    }
  }
  else if( fYZHits[2*i1].empty() ){
    track.addHit(fYZHits[2*i1+1][index]);
    return true;
  }
  else if( fYZHits[2*i1+1].empty() ){
    track.addHit(fYZHits[2*i1][index]);
    return true;
  }
  return false;
}

int BPCTrackMan::nXCluster(const int &i)
{
  int index=0;
  if( !fXZHits[2*i].empty() && !fXZHits[2*i+1].empty() ){
    for( int j1=0; j1<fXZHits[2*i].size(); j1++ ){
      for( int j2=0; j2<fXZHits[2*i+1].size(); j2++ ){
	if( fabs(fXZHits[2*i][j1]->wx()-fXZHits[2*i+1][j2]->wx())<0.54 ) index++;
#if 0
	else{
	  std::cout<<" layer"<<2*i+1<<" & layer"<<2*i+2<<std::endl;
	  std::cout<<"    wx1="<<fXZHits[2*i][j1]->wx()<<std::endl;
	  std::cout<<"    wx2="<<fXZHits[2*i+1][j2]->wx()<<std::endl;
	  std::cout<<" can not create cluster"<<std::endl;
	}
#endif
      }
    }
    
    for( int j1=0; j1<fXZHits[2*i].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fXZHits[2*i+1].size(); j2++ ){
	if( fabs(fXZHits[2*i][j1]->wx()-fXZHits[2*i+1][j2]->wx())<0.54 ) flag=false;
      }
      if( flag ){
	std::cout<<"aaa"<<std::endl;
	index++;
      }
    }
    
    for( int j1=0; j1<fXZHits[2*i+1].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fXZHits[2*i].size(); j2++ ){
	if( fabs(fXZHits[2*i][j2]->wx()-fXZHits[2*i+1][j1]->wx())<0.54 ) flag=false;
      }
      if( flag ){
	std::cout<<"bbb"<<std::endl;
	index++;
      }
    }
  }
  else if( fXZHits[2*i].empty() ) index = fXZHits[2*i+1].size();
  else if( fXZHits[2*i+1].empty() ) index = fXZHits[2*i].size();

  return index;
}

int BPCTrackMan::nYCluster(const int &i)
{
  int index=0;
  if( !fYZHits[2*i].empty() && !fYZHits[2*i+1].empty() ){
    for( int j1=0; j1<fYZHits[2*i].size(); j1++ ){
      for( int j2=0; j2<fYZHits[2*i+1].size(); j2++ ){
	if( fabs(fYZHits[2*i][j1]->wx()-fYZHits[2*i+1][j2]->wx())<0.54 ) index++;
#if 0
	else{
	  std::cout<<" layer"<<2*i+1<<" & layer"<<2*i+2<<std::endl;
	  std::cout<<"    wx1="<<fYZHits[2*i][j1]->wx()<<std::endl;
	  std::cout<<"    wx2="<<fYZHits[2*i+1][j2]->wx()<<std::endl;
	  std::cout<<" can not create cluster"<<std::endl;
	}
#endif
      }
    }
    
    for( int j1=0; j1<fYZHits[2*i].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fYZHits[2*i+1].size(); j2++ ){
	if( fabs(fYZHits[2*i][j1]->wx()-fYZHits[2*i+1][j2]->wx())<0.54 ) flag=false;
      }
      if( flag ) index++;
    }
    
    for( int j1=0; j1<fYZHits[2*i+1].size(); j1++ ){
      bool flag = true;
      for( int j2=0; j2<fYZHits[2*i].size(); j2++ ){
	if( fabs(fYZHits[2*i][j2]->wx()-fYZHits[2*i+1][j1]->wx())<0.54 ) flag=false;
      }
      if( flag ) index++;
    }
  }
  else if( fYZHits[2*i].empty() ) index = fYZHits[2*i+1].size();
  else if( fYZHits[2*i+1].empty() ) index = fYZHits[2*i].size();

  return index;
}

void BPCTrackMan::deleteHit(XYTrack &track)
{
  if( track.xy()==0 ){
    for( int i=0; i<track.nhit(); i++ ){
      ChamberLikeHit *hit = track.getHit(i);
      int index = 0;
      if( hit-> layer()==1 ) index=0;
      else if( hit-> layer()==2 ) index=1;
      else if( hit-> layer()==5 ) index=2;
      else if( hit-> layer()==6 ) index=3;
      else{
	std::cout<<" !!! Delete X !!!"<<std::endl;
	exit(-1);
      }
      for( std::vector<ChamberLikeHit*>::iterator ite = fXZHits[index].begin(); ite!=fXZHits[index].end(); ++ite ){
	if( hit->layer()==(*ite)->layer() && hit->wire()==(*ite)->wire() ){
	  fXZHits[index].erase(ite);
	  break;
	}
      }
    }
  }
  else if( track.xy()==1 ){
    for( int i=0; i<track.nhit(); i++ ){
      ChamberLikeHit *hit = track.getHit(i);
      int index = 0;
      if( hit-> layer()==3 ) index=0;
      else if( hit-> layer()==4 ) index=1;
      else if( hit-> layer()==7 ) index=2;
      else if( hit-> layer()==8 ) index=3;
      else{
	std::cout<<" !!! Delete Y !!!"<<std::endl;
	exit(-1);
      }
      for( std::vector<ChamberLikeHit*>::iterator ite = fYZHits[index].begin(); ite!=fYZHits[index].end(); ++ite ){
	if( hit->layer()==(*ite)->layer() && hit->wire()==(*ite)->wire() ){
	  fYZHits[index].erase(ite);
	  break;
	}
      }
    }
  }

}

void BPCTrackMan::dumpHits()
{
  for( int i=0; i<4; i++ ){
    std::cout<<" X Layer "<<i+1<<" : "<<fXZHits[i].size()<<std::endl;
  }

  for( int i=0; i<4; i++ ){
    std::cout<<" Y Layer "<<i+1<<" : "<<fYZHits[i].size()<<std::endl;
  }
}

void BPCTrackMan::dumpTrack()
{
  for( int i=0; i<fXTracks.size(); i++ ){
    fXTracks[i].dump();
  }

  for( int i=0; i<fYTracks.size(); i++ ){
    fYTracks[i].dump();
  }
}
