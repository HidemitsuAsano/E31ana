#include "MyHistBVC.h"

using namespace std;

void initHistT0BVC()
{
  for( int seg=1; seg<=8; seg++ ){
    new TNtuple(Form("BVC%d_slewing_info", seg), Form("BVC seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");

    new TH1F(Form("T0BVC%d_tof", seg), Form("T0-BVC seg%d TOF", seg), 600, -100, 200);
    new TH1F(Form("T0BVC%d_tof_CDS_nohit", seg), Form("T0-BVC seg%d TOF", seg), 600, -100, 200);
    new TH1F(Form("BVC%d_offset", seg),           Form("T0-BVC seg%d offset", seg), 600, -100, 200);
    new TH1F(Form("BVC%d_offset_CDS_nohit", seg), Form("T0-BVC seg%d offset", seg), 600, -100, 200);
  }

}

void fillHistT0BVC(BeamLineHitMan *blMan, AnaInfo *info)
{
  if( info->nBeam()==1 && info->beam(0)->flag() && info->beam(0)->pid()==Beam_Kaon ){
    BeamInfo *beam=info->beam();
    if( beam->nFDC1()==1 ){
      BLDCTrackInfo *trackFDC1=&beam->FDC1(0);
      vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
      if( BVChits.size()==1 ){
	HodoscopeLikeHit *BVChit=BVChits[0];
	TVector3 T0pos = beam->T0pos();
	double z=BVChit->z()-0.5;;

	TVector3 BVCpos=trackFDC1->GetPosatZ(z);
	double fl=(BVCpos-T0pos).Mag();
	double D5mom=beam->D5mom();
	double tof=BVChit->ctmean()-beam->T0time();
	double mom_out, calc_tof;
	ELossTools::CalcElossBeamTGeo(T0pos, BVCpos, D5mom, beam->mass(), mom_out, calc_tof);
      
	double eu=BVChit->eu(), ed=BVChit->ed();
	double tu=BVChit->tu()-beam->T0time()-calc_tof;
	double td=BVChit->td()-beam->T0time()-calc_tof;
	double ctu=BVChit->ctu()-beam->T0time()-calc_tof;
	double ctd=BVChit->ctd()-beam->T0time()-calc_tof;

	MyHistTools::fillTH(Form("T0BVC%d_tof", BVChit->seg()), tof);
	MyHistTools::fillTH(Form("BVC%d_offset", BVChit->seg()), tof-calc_tof);
	if( info->nCDS()==0 ){
	  MyHistTools::fillTH(Form("T0BVC%d_tof_CDS_nohit", BVChit->seg()), tof);
	  MyHistTools::fillTH(Form("BVC%d_offset_CDS_nohit", BVChit->seg()), tof-calc_tof);

	  TNtuple *tup=(TNtuple*)gFile->Get(Form("BVC%d_slewing_info", BVChit->seg()));
	  tup->Fill(eu, ed, tu, td, ctu, ctd, tof-calc_tof);
	}
      }
    }
  }
}
