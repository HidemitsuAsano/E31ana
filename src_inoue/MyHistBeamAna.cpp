#include "MyHistBeamAna.h"

using namespace std;

static const double BLC1_TIME_WINDOW_MIN=-30;
static const double BLC1_TIME_WINDOW_MAX=100;

static const double BLC1_TIME_MIN=-10;
static const double BLC1_TIME_MAX=10;
static const double BLC1_CHI2_MAX=30;

void fillReadBeam_BLC1(EventHeader *header, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);
  int nBLC1;
  BLDCTrackInfo infoBLC1;
  for( int i=0; i<beam->nBLC1(); i++ ){
    BLDCTrackInfo tmpBLC1=beam->BLC1(i);
    MyHistTools::fillTH("BLC1_time", tmpBLC1.time());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC1_time_Kf", tmpBLC1.time());

    if( BLC1_TIME_WINDOW_MIN<infoBLC1.time() && infoBLC1.time()<BLC1_TIME_WINDOW_MAX ){
      nBLC1++;
      infoBLC1=tmpBLC1;
    }
  }
  if( nBLC1!=1 ) return;

  if( BLC1_TIME_MIN<infoBLC1.time() || infoBLC1.time()<BLC1_TIME_MAX ){
    if( infoBLC1.chi2()<BLC1_CHI2_MAX ) return;
  }

  return;
}

void initHistBeamAna()
{
  new TH2F("BLC2BPC_diff", "BLC2 BPC diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff", "BLC2 BPC diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);

  new TH2F("BeamProf_FF",      "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_FF_Kf",   "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_FF_CDH2", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);

  new TH1F("BeamMom", "Beam Mom", 1000, 0.5, 1.5);
  new TH1F("D5_chi2", "D5 chi-square", 1000, 0, 100);

  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("BHDT0_T0%d_offset", seg), Form("T0 seg%d offset by BHD-T0", seg), 200, -5.0, 5.0);
    new TH1F(Form("BHDT0_T0%d_offset_pi", seg), Form("T0 seg%d offset by BHD-T0", seg), 200, -5.0, 5.0);
  }

  for( int seg=1; seg<=20; seg++ ){
    new TH1F(Form("BHDT0_BHD%d_offset", seg), Form("BHD seg%d offset by BHD-T0", seg), 200, -5.0, 5.0);
    new TH1F(Form("BHDT0_BHD%d_offset_pi", seg), Form("BHD seg%d offset by BHD-T0", seg), 200, -5.0, 5.0);
    new TNtuple(Form("BHD%d_slewing_info", seg), Form("BHD seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");

    for( int seg2=1; seg2<=5; seg2++ ){
      new TH1F(Form("BHDT0_BHD%d_T0%d_offset", seg, seg2), Form("BHD seg%d T0 seg%d offset by BHD-T0", seg, seg2), 200, -5.0, 5.0);
      new TH1F(Form("BHDT0_BHD%d_T0%d_offset_pi", seg, seg2), Form("BHD seg%d T0 seg%d offset by BHD-T0", seg, seg2), 200, -5.0, 5.0);
    }
  }
}

void fillBeamProf(EventHeader *header, BeamLineTrackMan *bltrackMan)
{
  LocalTrack *BPC=MyTools::trackBPC(bltrackMan);
  if( BPC ){
    TVector3 profFF=BPC->GetPosatZ(0);
    MyHistTools::fillTH("BeamProf_FF", profFF.X(), profFF.Y());
    if( header->trigmode2(Mode_Kf)    ) MyHistTools::fillTH("BeamProf_FF_Kf", profFF.X(), profFF.Y());
    if( header->trigmode2(Mode_KCDH2) ) MyHistTools::fillTH("BeamProf_FF_CDH2", profFF.X(), profFF.Y());
  }
}
void fillBLC2BPC(BeamLineTrackMan *bltrackMan)
{
  LocalTrack *BLC2=MyTools::trackBLC2(bltrackMan);
  LocalTrack *BPC =MyTools::trackBPC(bltrackMan);
  if( BLC2 && BPC ){
    TVector3 BPCmom_dir=BPC->GetMomDir();
    TVector3 BLC2mom_dir=BLC2->GetMomDir();
    TVector3 mom_diff=BLC2mom_dir-BPCmom_dir;

    TVector3 BPC_pos=BPC->GetPosatZ(0.5*(-130-20.3));
    TVector3 BLC2_pos=BLC2->GetPosatZ(0.5*(-130-20.3));
    TVector3 BLC2BPC_diff=BLC2_pos-BPC_pos;

    MyHistTools::fillTH("BLC2BPC_diff", BLC2BPC_diff.X(), BLC2BPC_diff.Y());
    MyHistTools::fillTH("BLC2BPC_dir_diff", mom_diff.X(), mom_diff.Y());
  }
}

void fillBeamSpectrometer(BeamSpectrometer *beam)
{
  if( beam->stat()>0 ){
    MyHistTools::fillTH("BeamMom", beam->mom());
    MyHistTools::fillTH("D5_chi2", beam->chisquare());
  }
}

void fillBHDT0offset(BeamInfo *info, BeamLineHitMan *blMan)
{
  if( info->flag() ){

    double calc_tof=info->calc_tof();
    HodoscopeLikeHit *T0hit=info->T0(blMan);
    HodoscopeLikeHit *BHDhit=info->BHD(blMan);
    int T0seg=T0hit->seg();
    int BHDseg=BHDhit->seg();
    double tof=T0hit->ctmean()-BHDhit->ctmean();
    double eu=BHDhit->eu(), ed=BHDhit->ed();
    double tu=BHDhit->tu()-T0hit->ctmean()+calc_tof;
    double td=BHDhit->td()-T0hit->ctmean()+calc_tof;
    double ctu=BHDhit->ctu()-T0hit->ctmean()+calc_tof;
    double ctd=BHDhit->ctd()-T0hit->ctmean()+calc_tof;

    if( info->pid()==Beam_Kaon ){
      //    cout<<"BHD-T0 : "<<tof<<"  calc : "<<calc_tof<<endl;
      if( T0seg==3 ){
	TNtuple *tup=(TNtuple*)gFile->Get(Form("BHD%d_slewing_info", BHDseg));
	tup->Fill(eu, ed, tu, td, ctu, ctd, calc_tof-tof);
      }
      
      MyHistTools::fillTH(Form("BHDT0_T0%d_offset", T0seg),   tof-calc_tof);
      MyHistTools::fillTH(Form("BHDT0_BHD%d_offset", BHDseg), calc_tof-tof);
      MyHistTools::fillTH(Form("BHDT0_BHD%d_T0%d_offset", BHDseg, T0seg), calc_tof-tof);
    }

    if( info->pid()==Beam_Pion ){      
      MyHistTools::fillTH(Form("BHDT0_T0%d_offset_pi", T0seg),   tof-calc_tof);
      MyHistTools::fillTH(Form("BHDT0_BHD%d_offset_pi", BHDseg), calc_tof-tof);
      MyHistTools::fillTH(Form("BHDT0_BHD%d_T0%d_offset_pi", BHDseg, T0seg), calc_tof-tof);
    }
  }
}
