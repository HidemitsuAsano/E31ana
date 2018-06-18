#include "MyHistTrigChk.h"

using namespace std;

void initHistTrigChk()
{
  new TH1F("trig_chk_K",    "trig_chk_K",    10, 0, 10);
  new TH1F("trig_chk_CDH1", "trig_chk_CDH1", 10, 0, 10);
  new TH1F("trig_chk_CDH2", "trig_chk_CDH2", 10, 0, 10);
  new TH1F("trig_chk_N",    "trig_chk_N",    10, 0, 10);
  new TH1F("trig_chk_N2",   "trig_chk_N2",   10, 0, 10);
  new TH1F("trig_chk_C",    "trig_chk_C",    10, 0, 10);

  new TH1F("T0CDH_tof_out_trigger", "", 1000, -100, 100);
  new TH1F("T0NC_tof_out_trigger", "",  1000, -100, 100);

  new TH1F("CDH_time", "", 1000, -100, 100);
  new TH1F("CDH_time_wTrig", "", 1000, -100, 100);

  for( int i=0; i<20; i++ ){
    new TH1F(Form("trig_%d",i), Form("trig_%d", i), 4000, 0, 4000);
  }
}

void fillHistTrigChk(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, AnaInfo *anaInfo)
{
  for( int i=0; i<20; i++ ){
    if( header->pattern(i)<0 ) continue;
    MyHistTools::fillTH(Form("trig_%d", i), header->pattern(i));
  }

  if( header->trigmode(Mode_Beam) ){
    bool AC_flag=false;
    vector<HodoscopeLikeHit*> DEFhits=MyTools::getHodo(blMan, CID_DEF);
    for( int i=0; i<blMan->nAC(); i++ ){
      CherenkovLikeHit *hit=blMan->AC(i);
      if( hit->CheckRange(1) ) AC_flag=true;
    }
    if( !AC_flag && DEFhits.size()>0 ){
      MyHistTools::fillTH("trig_chk_K", 0);
      if( header->IsTrig(Trig_Kaon) ) MyHistTools::fillTH("trig_chk_K", 1);
    }
  }

  if( header->trigmode(Mode_Kf) ){
    vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
    vector<HodoscopeLikeHit*> CDHhits=MyTools::getCDH(cdsMan);
    bool time_flag=false;
    for( int i=0; i<CDHhits.size(); i++ ){
      MyHistTools::fillTH("CDH_time", CDHhits[i]->ctmean());
      if( header-> IsTrig(Trig_KCDH1) ) MyHistTools::fillTH("CDH_time_wTrig", CDHhits[i]->ctmean());
      if(  0<CDHhits[i]->ctmean() && CDHhits[i]->ctmean()<40 ) time_flag=true;
    }
    if( time_flag ){
      MyHistTools::fillTH("trig_chk_CDH1", 8);
      if( header->IsTrig(Trig_KCDH1) ) MyHistTools::fillTH("trig_chk_CDH1", 9);
    }

    if( CDHhits.size()>0 ){
      MyHistTools::fillTH("trig_chk_CDH1", 0);
      if( header->IsTrig(Trig_KCDH1) ) MyHistTools::fillTH("trig_chk_CDH1", 1);

      if( T0hits.size()==1 ){
	MyHistTools::fillTH("trig_chk_CDH1", 2);
	if( header->IsTrig(Trig_KCDH1) ) MyHistTools::fillTH("trig_chk_CDH1", 3);
	bool tof_flag=false;
	for( int i=0; i<CDHhits.size(); i++ ){
	  double tof=CDHhits[i]->ctmean()-T0hits[0]->ctmean();
	  if( -10<tof && tof<50 ) tof_flag=true;
	  if( !header->IsTrig(Trig_KCDH1) ){
	    MyHistTools::fillTH("T0CDH_tof_out_trigger", tof);
	  }
	}
	if( tof_flag ){
	  MyHistTools::fillTH("trig_chk_CDH1", 4);
	  if( header->IsTrig(Trig_KCDH1) ) MyHistTools::fillTH("trig_chk_CDH1", 5);
	}
	if( MyAnaTools::goodBeam(anaInfo) ){
	  MyHistTools::fillTH("trig_chk_CDH1", 6);
	  if( header->IsTrig(Trig_KCDH1) ) MyHistTools::fillTH("trig_chk_CDH1", 7);
	}
      }
    }
  }


  if( header->trigmode(Mode_Kf) ){
    vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
    vector<HodoscopeLikeHit*> CDHhits=MyTools::getCDH(cdsMan);
    if( CDHhits.size()>1 ){
      MyHistTools::fillTH("trig_chk_CDH2", 0);
      if( header->IsTrig(Trig_KCDH2) ) MyHistTools::fillTH("trig_chk_CDH2", 1);
      if( T0hits.size()==1 ){
	MyHistTools::fillTH("trig_chk_CDH2", 2);
	if( header->IsTrig(Trig_KCDH2) ) MyHistTools::fillTH("trig_chk_CDH2", 3);
	int n_tof_flag=0;
	for( int i=0; i<CDHhits.size(); i++ ){
	  double tof=CDHhits[i]->ctmean()-T0hits[0]->ctmean();
	  if( -10<tof && tof<50 ) n_tof_flag++;
	  if( !header->IsTrig(Trig_KCDH2) ){
	    MyHistTools::fillTH("T0CDH_tof_out_trigger", tof);
	  }
	}
	if( n_tof_flag>1 ){
	  MyHistTools::fillTH("trig_chk_CDH2", 4);
	  if( header->IsTrig(Trig_KCDH2) ) MyHistTools::fillTH("trig_chk_CDH2", 5);
	}

	if( MyAnaTools::goodBeam(anaInfo) ){
	  MyHistTools::fillTH("trig_chk_CDH2", 6);
	  if( header->IsTrig(Trig_KCDH2) ) MyHistTools::fillTH("trig_chk_CDH2", 7);
	}
      }
    }
  }

  if( header->trigmode(Mode_Kf) ){
    vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
    vector<HodoscopeLikeHit*> CDHhits=MyTools::getCDH(cdsMan);
    vector<HodoscopeLikeHit*> NChits=MyTools::getHodo(blMan, CID_NC);
    vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
    vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
    if( NChits.size()>0 && CVChits.size()==0 && BVChits.size()==0 && CDHhits.size()>0 ){
      MyHistTools::fillTH("trig_chk_N2", 0);
      if( header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("trig_chk_N2", 1);

      if( T0hits.size()==1 ){
	MyHistTools::fillTH("trig_chk_N2", 2);
	if( header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("trig_chk_N2", 3);
	if( MyAnaTools::goodBeam(anaInfo) ){
	  MyHistTools::fillTH("trig_chk_N2", 4);
	  if( header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("trig_chk_N2", 5);
	}
      }
    }
  }

  if( header->trigmode()==Mode_KCDH1f ){
    vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
    vector<HodoscopeLikeHit*> NChits=MyTools::getHodo(blMan, CID_NC);
    vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
    vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
    if( NChits.size()>0 && CVChits.size()==0 && BVChits.size()==0 ){
      MyHistTools::fillTH("trig_chk_N", 0);
      if( header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("trig_chk_N", 1);
      if( T0hits.size()==1 ){
	for( int i=0; i<NChits.size(); i++ ){
	  double tof=NChits[i]->ctmean()-T0hits[0]->ctmean();
	  if( !header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("T0NC_tof_out_trigger", tof);
	}
	MyHistTools::fillTH("trig_chk_N", 2);
	if( header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("trig_chk_N", 3);

	if( MyAnaTools::goodBeam(anaInfo) ){
	  MyHistTools::fillTH("trig_chk_N", 4);
	  if( header->IsTrig(Trig_Neutral) ) MyHistTools::fillTH("trig_chk_N", 5);
	}
      }
    }
  }

  if( header->trigmode()==Mode_KCDH1f ){
    vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
    vector<HodoscopeLikeHit*> NChits=MyTools::getHodo(blMan, CID_NC);
    vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
    vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
    vector<HodoscopeLikeHit*> PChits=MyTools::getHodo(blMan, CID_PC);
    bool CVC_flag=false;
    for( int i=0; i<CVChits.size(); i++ ){
      if( CVChits[i]->seg()>17 ) CVC_flag=true;
    }

    if( CVC_flag || PChits.size()>0 ){
      MyHistTools::fillTH("trig_chk_C", 0);
      if( header->IsTrig(Trig_Charged) ) MyHistTools::fillTH("trig_chk_C", 1);

      if( T0hits.size()==1 ){
	MyHistTools::fillTH("trig_chk_C", 2);
	if( header->IsTrig(Trig_Charged) ) MyHistTools::fillTH("trig_chk_C", 3);
	if( MyAnaTools::goodBeam(anaInfo) ){
	  MyHistTools::fillTH("trig_chk_C", 6);
	  if( header->IsTrig(Trig_Charged) ) MyHistTools::fillTH("trig_chk_C", 7);
	}
    }
    }
  }
}
