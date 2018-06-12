#include "MyHistDEF.h"

using namespace std;

void initHistT0DEF()
{
  for( int seg=1; seg<=8; seg++ ){
    new TH1F(Form("T0DEF%d_tof", seg), Form("T0-DEF seg%d TOF", seg), 600, -100, 200);
    new TH1F(Form("DEF%d_offset", seg), Form("T0-DEF seg%d offset", seg), 600, -100, 200);
  }
}

void fillHistT0DEF(BeamLineHitMan *blMan, AnaInfo *info)
{
  if( info->nBeam()==1 && info->beam(0)->flag() && info->beam(0)->pid()==Beam_Kaon ){
    BeamInfo *beam=info->beam();
    BLDCTrackInfo *trackBPC=&beam->BPC(0);
    vector<HodoscopeLikeHit*> DEFhits=MyTools::getHodo(blMan, CID_DEF);
    if( DEFhits.size()==1 ){
      HodoscopeLikeHit *DEFhit=DEFhits[0];
      TVector3 T0pos = beam->T0pos();
      double z=DEFhit->z()-0.15;
      TVector3 DEFpos=trackBPC->GetPosatZ(z);
      double fl=(DEFpos-T0pos).Mag();
      double D5mom=beam->D5mom();
      
      double tof=DEFhit->ctmean()-beam->T0time();
      double mom_out, calc_tof;
      ELossTools::CalcElossBeamTGeo(T0pos, DEFpos, D5mom, beam->mass(), mom_out, calc_tof);

      MyHistTools::fillTH(Form("T0DEF%d_tof", DEFhit->seg()), tof);
      MyHistTools::fillTH(Form("DEF%d_offset", DEFhit->seg()), tof-calc_tof);
    }
  }
}
