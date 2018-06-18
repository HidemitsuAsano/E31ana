#include "MyHistT0PCCVC.h"

using namespace std;

void initHistT0PCCVC()
{
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("CVC%d_offset", seg), Form("CVC seg%d offset", seg), 2000, -50, 50);
    new TNtuple(Form("CVC%d_slewing_info", seg), Form("CVC seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");
  }

  for( int seg=1; seg<=27; seg++ ){
    new TH1F(Form("PC%d_offset", seg), Form("PC seg%d offset", seg), 2000, -50, 50);
    new TNtuple(Form("PC%d_slewing_info", seg), Form("PC seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");
  }
}

void fillT0PCCVC_BT(ConfMan *conf, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamInfo *beam)
{
  if( beam->flag() ){
    std::vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
    std::vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
    std::vector<HodoscopeLikeHit*> PChits=MyTools::getHodo(blMan, CID_PC);

    if( beam->nFDC1()==1 ){
      HodoscopeLikeHit *fc_hit=0;
      if( CVChits.size()==1 && PChits.size()==0 ) fc_hit=CVChits[0];
      if( PChits.size()==1 && CVChits.size()==0 ) fc_hit=PChits[0];

      if( fc_hit ){	
	ForwardChargeInfo fcInfo;
	fcInfo.SetFDC1(beam->FDC1(0));
	fcInfo.SetHodo(fc_hit);

	if( fcInfo.calc_simple_beam_through(*beam, conf, blMan) ){
	  HodoscopeLikeHit *T0hit=beam->T0(blMan);

	  double eu=fc_hit->eu(), ed=fc_hit->ed();
	  double tu=fc_hit->tu()-beam->T0time()-fcInfo.calc_tof();
	  double td=fc_hit->td()-beam->T0time()-fcInfo.calc_tof();
	  double ctu=fc_hit->ctu()-beam->T0time()-fcInfo.calc_tof();
	  double ctd=fc_hit->ctd()-beam->T0time()-fcInfo.calc_tof();

	  double tof=fc_hit->ctmean()-beam->T0time();
	  double offset=tof-fcInfo.calc_tof();
	  //	  cout<<"PC/CVC offset"<<offset<<endl;
	  if( fc_hit->cid()==CID_CVC ){
	    MyHistTools::fillTH(Form("CVC%d_offset",fc_hit->seg()), offset);
	    TNtuple *tup=(TNtuple*)gFile->Get(Form("CVC%d_slewing_info", fc_hit->seg()));
	    tup->Fill(eu, ed, tu, td, ctu, ctd, offset);
	  }
	  if( fc_hit->cid()==CID_PC ){
	    MyHistTools::fillTH(Form("PC%d_offset",fc_hit->seg()), offset);

	    TNtuple *tup=(TNtuple*)gFile->Get(Form("PC%d_slewing_info", fc_hit->seg()));
	    tup->Fill(eu, ed, tu, td, ctu, ctd, offset);
	  }
	}
      }
    }
  }
}
