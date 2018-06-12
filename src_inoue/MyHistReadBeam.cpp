#include "MyHistReadBeam.h"

using namespace std;

void fillReadBeam_Profile(ConfMan *conf, EventHeader *header, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);
  int nBPC=0;
  BLDCTrackInfo infoBPC;
  for( int i=0; i<beam->nBPC(); i++ ){
    BLDCTrackInfo tmpBPC=beam->BPC(i);
    if( BPC_TIME_WINDOW_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_WINDOW_MAX ){
      nBPC++;
      infoBPC=tmpBPC;
    }
  }
  if( nBPC!=1 ) return;

  TVector3 pos;
  conf-> GetGeomMapManager()->GetPos(CID_TarCell, 0, pos);
  TVector3 prof_at_FF=infoBPC.GetPosatZ(pos.Z());
  MyHistTools::fillTH("BeamProf", prof_at_FF.X(), prof_at_FF.Y());
  if( GeomTools::GetID(prof_at_FF)==CID_Fiducial )  MyHistTools::fillTH("BeamProf_true", prof_at_FF.X(), prof_at_FF.Y());
  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("BeamProf_Kf", prof_at_FF.X(), prof_at_FF.Y());
    if( GeomTools::GetID(prof_at_FF)==CID_Fiducial )  MyHistTools::fillTH("BeamProf_Kf_true", prof_at_FF.X(), prof_at_FF.Y());
  }
  if( header->IsTrig(Trig_KCDH2f) ){
    MyHistTools::fillTH("BeamProf_KCDH2", prof_at_FF.X(), prof_at_FF.Y());
  }
}

void fillReadBeam_D5BHD(EventHeader *header, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  if( anaInfo->beam(0)->D5chi2()<-1 || D5_CHI2_MAX<anaInfo->beam(0)->D5chi2() ) return;
  BeamInfo *beam=anaInfo->beam(0);
  double D5mom=beam->D5mom();
  if( beam->nBHD()==1 ){
    int BHDseg;
    double  BHDtime;
    beam->getBHD(0, BHDseg, BHDtime);
    MyHistTools::fillTH(Form("D5_mom_BHD%d", BHDseg), D5mom);
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH(Form("D5_mom_BHD%d_Kf", BHDseg), D5mom);
  }
}

void fillReadBeam_BLC2BPC(EventHeader *header, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;

  BeamInfo *beam=anaInfo->beam(0);
  int nBLC2=0;
  BLDCTrackInfo infoBLC2;
  for( int i=0; i<beam->nBLC2(); i++ ){
    BLDCTrackInfo tmpBLC2=beam->BLC2(i);
    if( BLC2_TIME_WINDOW_MIN<tmpBLC2.time() && tmpBLC2.time()<BLC2_TIME_WINDOW_MAX ){
      nBLC2++;
      infoBLC2=tmpBLC2;
    }
  }
  if( nBLC2!=1 ) return;

  int nBPC=0;
  BLDCTrackInfo infoBPC;
  for( int i=0; i<beam->nBPC(); i++ ){
    BLDCTrackInfo tmpBPC=beam->BPC(i);
    if( BPC_TIME_WINDOW_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_WINDOW_MAX ){
      nBPC++;
      infoBPC=tmpBPC;
    }
  }
  if( nBPC!=1 ) return;

  double z = 0.5*(-130-20.3);
  TVector3 dir_diff=infoBLC2.dir()-infoBPC.dir();
  TVector3 pos_diff=infoBLC2.GetPosatZ(z)-infoBPC.GetPosatZ(z);

  MyHistTools::fillTH("BLC2BPC_diff", pos_diff.X(), pos_diff.Y());
  MyHistTools::fillTH("BLC2BPC_diff_X", pos_diff.X());
  if( BLC2BPC_X_MIN<pos_diff.X() && pos_diff.X()<BLC2BPC_X_MAX ) MyHistTools::fillTH("BLC2BPC_diff_X_true", pos_diff.X());
  MyHistTools::fillTH("BLC2BPC_diff_Y", pos_diff.Y());
  if( BLC2BPC_Y_MIN<pos_diff.Y() && pos_diff.Y()<BLC2BPC_Y_MAX ) MyHistTools::fillTH("BLC2BPC_diff_Y_true", pos_diff.Y());
  MyHistTools::fillTH("BLC2BPC_dir_diff", dir_diff.X(), dir_diff.Y());
  MyHistTools::fillTH("BLC2BPC_dir_diff_X", dir_diff.X());
  if( BLC2BPC_dX_MIN<dir_diff.X() && dir_diff.X()<BLC2BPC_dX_MAX ) MyHistTools::fillTH("BLC2BPC_dir_diff_X_true", dir_diff.X());
  MyHistTools::fillTH("BLC2BPC_dir_diff_Y", dir_diff.Y());
  if( BLC2BPC_dY_MIN<dir_diff.Y() && dir_diff.Y()<BLC2BPC_dY_MAX ) MyHistTools::fillTH("BLC2BPC_dir_diff_Y_true", dir_diff.Y());

  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("BLC2BPC_diff_Kf", pos_diff.X(), pos_diff.Y());
    MyHistTools::fillTH("BLC2BPC_diff_X_Kf", pos_diff.X());
    if( BLC2BPC_X_MIN<pos_diff.X() && pos_diff.X()<BLC2BPC_X_MAX ) MyHistTools::fillTH("BLC2BPC_diff_X_Kf_true", pos_diff.X());
    MyHistTools::fillTH("BLC2BPC_diff_Y_Kf", pos_diff.Y());
    if( BLC2BPC_Y_MIN<pos_diff.Y() && pos_diff.Y()<BLC2BPC_Y_MAX ) MyHistTools::fillTH("BLC2BPC_diff_Y_Kf_true", pos_diff.Y());
    MyHistTools::fillTH("BLC2BPC_dir_diff_Kf", dir_diff.X(), dir_diff.Y());
    MyHistTools::fillTH("BLC2BPC_dir_diff_X_Kf", dir_diff.X());
    if( BLC2BPC_dX_MIN<dir_diff.X() && dir_diff.X()<BLC2BPC_dX_MAX ) MyHistTools::fillTH("BLC2BPC_dir_diff_X_Kf_true", dir_diff.X());
    MyHistTools::fillTH("BLC2BPC_dir_diff_Y_Kf", dir_diff.Y());
    if( BLC2BPC_dY_MIN<dir_diff.Y() && dir_diff.Y()<BLC2BPC_dY_MAX ) MyHistTools::fillTH("BLC2BPC_dir_diff_Y_Kf_true", dir_diff.Y());
  }
}

void fillReadBeam_D5(EventHeader *header, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;

  MyHistTools::fillTH("D5_chi2", anaInfo->beam(0)->D5chi2());
  MyHistTools::fillTH("D5_mom", anaInfo->beam(0)->D5mom());
  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("D5_chi2_Kf", anaInfo->beam(0)->D5chi2());
    MyHistTools::fillTH("D5_mom_Kf", anaInfo->beam(0)->D5mom());
  }
  if( -1<anaInfo->beam(0)->D5chi2() && anaInfo->beam(0)->D5chi2()<D5_CHI2_MAX ){
    MyHistTools::fillTH("D5_chi2_true", anaInfo->beam(0)->D5chi2());
    MyHistTools::fillTH("D5_mom_true", anaInfo->beam(0)->D5mom());
    if( header->IsTrig(Trig_Kf) ){
      MyHistTools::fillTH("D5_chi2_Kf_true", anaInfo->beam(0)->D5chi2());
      MyHistTools::fillTH("D5_mom_Kf_true", anaInfo->beam(0)->D5mom());
    }
  }
}

void fillReadBeam_BLC1(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);
  if( MyAnaTools::trigBLC1(anaInfo, conf) ){
    if( header->IsTrig(Trig_Kf) ){
      MyHistTools::fillTH("BLC1_eff_Kf", 0);
      if( beam->nBLC1()>0 ) MyHistTools::fillTH("BLC1_eff_Kf", 1);
    }
  }

  int nBLC1=0;
  int nBLC1_2=0;
  BLDCTrackInfo infoBLC1;
  for( int i=0; i<beam->nBLC1(); i++ ){
    BLDCTrackInfo tmpBLC1=beam->BLC1(i);
    MyHistTools::fillTH("BLC1_time", tmpBLC1.time());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC1_time_Kf", tmpBLC1.time());
    if( BLC1_TIME_MIN<tmpBLC1.time() && tmpBLC1.time()<BLC1_TIME_MAX ){
      MyHistTools::fillTH("BLC1_time_true", tmpBLC1.time());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC1_time_Kf_true", tmpBLC1.time());
    }
    if( -30<tmpBLC1.time() && tmpBLC1.time()<200 ) nBLC1_2++;

    if( BLC1_TIME_WINDOW_MIN<tmpBLC1.time() && tmpBLC1.time()<BLC1_TIME_WINDOW_MAX ){
      nBLC1++;
      infoBLC1=tmpBLC1;
    }
  }
  MyHistTools::fillTH("nBLC1_all", beam->nBLC1());
  MyHistTools::fillTH("nBLC1_2", nBLC1_2);
  MyHistTools::fillTH("nBLC1", nBLC1);
  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("nBLC1_all_Kf", beam->nBLC1());
    MyHistTools::fillTH("nBLC1_2_Kf", nBLC1_2);
    MyHistTools::fillTH("nBLC1_Kf", nBLC1);
  }
  if( nBLC1!=1 ) return;

  if( BLC1_TIME_MIN<infoBLC1.time() || infoBLC1.time()<BLC1_TIME_MAX ){
    MyHistTools::fillTH("BLC1_chi2", infoBLC1.chi2());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC1_chi2_Kf", infoBLC1.chi2());
    if( infoBLC1.chi2()<BLC1_CHI2_MAX ){
      MyHistTools::fillTH("BLC1_chi2_true", infoBLC1.chi2());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC1_chi2_Kf_true", infoBLC1.chi2());
    }
  }
}

void fillReadBeam_BLC2(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);
  if( MyAnaTools::trigBLC2(anaInfo, conf) ){
    if( header->IsTrig(Trig_Kf) ){
      MyHistTools::fillTH("BLC2_eff_Kf", 0);
      if( beam->nBLC2()>0 ) MyHistTools::fillTH("BLC2_eff_Kf", 1);
    }
  }

  int nBLC2=0;
  int nBLC2_2=0;
  BLDCTrackInfo infoBLC2;
  for( int i=0; i<beam->nBLC2(); i++ ){
    BLDCTrackInfo tmpBLC2=beam->BLC2(i);
    MyHistTools::fillTH("BLC2_time", tmpBLC2.time());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC2_time_Kf", tmpBLC2.time());
    if( BLC2_TIME_MIN<tmpBLC2.time() && tmpBLC2.time()<BLC2_TIME_MAX ){
      MyHistTools::fillTH("BLC2_time_true", tmpBLC2.time());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC2_time_Kf_true", tmpBLC2.time());
    }
    if( -30<tmpBLC2.time() && tmpBLC2.time()<200 ) nBLC2_2++;

    if( BLC2_TIME_WINDOW_MIN<tmpBLC2.time() && tmpBLC2.time()<BLC2_TIME_WINDOW_MAX ){
      nBLC2++;
      infoBLC2=tmpBLC2;
    }
  }
  MyHistTools::fillTH("nBLC2_all", beam->nBLC2());
  MyHistTools::fillTH("nBLC2_2", nBLC2_2);
  MyHistTools::fillTH("nBLC2", nBLC2);
  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("nBLC2_all_Kf", beam->nBLC2());
    MyHistTools::fillTH("nBLC2_2_Kf", nBLC2_2);
    MyHistTools::fillTH("nBLC2_Kf", nBLC2);
  }
  if( nBLC2!=1 ) return;

  if( BLC2_TIME_MIN<infoBLC2.time() || infoBLC2.time()<BLC2_TIME_MAX ){
    MyHistTools::fillTH("BLC2_chi2", infoBLC2.chi2());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC2_chi2_Kf", infoBLC2.chi2());
    if( infoBLC2.chi2()<BLC2_CHI2_MAX ){
      MyHistTools::fillTH("BLC2_chi2_true", infoBLC2.chi2());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BLC2_chi2_Kf_true", infoBLC2.chi2());
    }
  }
}

void fillReadBeam_BPC(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo, BeamLineHitMan *blMan)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);
  if( header->IsTrig(Trig_Kf) ){
    if( MyAnaTools::trigBPC(anaInfo, conf, blMan) ){
      MyHistTools::fillTH("BPC_eff_Kf", 0);
      if( beam->nBPC()>0 ) MyHistTools::fillTH("BPC_eff_Kf", 1);
    }
  }

  int nBPC=0;
  int nBPC_2=0;
  BLDCTrackInfo infoBPC;
  for( int i=0; i<beam->nBPC(); i++ ){
    BLDCTrackInfo tmpBPC=beam->BPC(i);
    MyHistTools::fillTH("BPC_time", tmpBPC.time());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BPC_time_Kf", tmpBPC.time());
    if( BPC_TIME_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_MAX ){
      MyHistTools::fillTH("BPC_time_true", tmpBPC.time());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BPC_time_Kf_true", tmpBPC.time());
    }
    if( -30<tmpBPC.time() && tmpBPC.time()<200 ) nBPC_2++;

    if( BPC_TIME_WINDOW_MIN<tmpBPC.time() && tmpBPC.time()<BPC_TIME_WINDOW_MAX ){
      nBPC++;
      infoBPC=tmpBPC;
    }
  }
  MyHistTools::fillTH("nBPC_all", beam->nBPC());
  MyHistTools::fillTH("nBPC_2", nBPC_2);
  MyHistTools::fillTH("nBPC", nBPC);
  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("nBPC_all_Kf", beam->nBPC());
    MyHistTools::fillTH("nBPC_2_Kf", nBPC_2);
    MyHistTools::fillTH("nBPC_Kf", nBPC);
  }
  if( nBPC!=1 ) return;

  if( BPC_TIME_MIN<infoBPC.time() || infoBPC.time()<BPC_TIME_MAX ){
    MyHistTools::fillTH("BPC_chi2", infoBPC.chi2());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BPC_chi2_Kf", infoBPC.chi2());
    if( infoBPC.chi2()<BPC_CHI2_MAX ){
      MyHistTools::fillTH("BPC_chi2_true", infoBPC.chi2());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("BPC_chi2_Kf_true", infoBPC.chi2());
    }
  }
}

void fillReadBeam_FDC1(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo, BeamLineHitMan *blMan)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);
  int nFDC1=0;
  if( header->IsTrig(Trig_Kf) ){
    if( MyAnaTools::trigFDC1(anaInfo, conf, blMan) ){
      MyHistTools::fillTH("FDC1_eff_Kf", 0);
      if( beam->nFDC1()>0 ) MyHistTools::fillTH("FDC1_eff_Kf", 1);
    }
  }

  BLDCTrackInfo infoFDC1;
  for( int i=0; i<beam->nFDC1(); i++ ){
    BLDCTrackInfo tmpFDC1=beam->FDC1(i);
    MyHistTools::fillTH("FDC1_time", tmpFDC1.time());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("FDC1_time_Kf", tmpFDC1.time());
    if( FDC1_TIME_MIN<tmpFDC1.time() && tmpFDC1.time()<FDC1_TIME_MAX ){
      MyHistTools::fillTH("FDC1_time_true", tmpFDC1.time());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("FDC1_time_Kf_true", tmpFDC1.time());
    }
    if( FDC1_TIME_WINDOW_MIN<tmpFDC1.time() && tmpFDC1.time()<FDC1_TIME_WINDOW_MAX ){
      nFDC1++;
      infoFDC1=tmpFDC1;
    }
  }
  MyHistTools::fillTH("nFDC1", nFDC1);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("nFDC1_Kf", nFDC1);
  if( nFDC1!=1 ) return;

  if( FDC1_TIME_MIN<infoFDC1.time() || infoFDC1.time()<FDC1_TIME_MAX ){
    MyHistTools::fillTH("FDC1_chi2", infoFDC1.chi2());
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("FDC1_chi2_Kf", infoFDC1.chi2());
    if( infoFDC1.chi2()<FDC1_CHI2_MAX ){
      MyHistTools::fillTH("FDC1_chi2_true", infoFDC1.chi2());
      if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("FDC1_chi2_Kf_true", infoFDC1.chi2());
    }
  }
}

void initHistReadBeam()
{
  new TH1F("EventReduction", "Event Reduction", 20, -0.5, 19.5);
  new TH1F("Kf_Reduction",   "K/f Reduction", 20, -0.5, 19.5);

  new TH1F("nT0_Kf", "T0 multiplicity", 6, -0.5, 5.5);
  new TH1F("nT0", "T0 multiplicity", 6, -0.5, 5.5);

  new TH1F("nBHD_Kf", "BHD multiplicity", 20, -0.5, 20.5);
  new TH1F("nBHD", "BHD multiplicity", 20, -0.5, 20.5);

  new TH1F("BHDT0_Kf", "BHD-T0 TOF", 50000, 0, 50);
  new TH1F("BHDT0", "BHD-T0 TOF", 50000, 0, 50);
  new TH1F("BHDT0_Kf_true", "BHD-T0 TOF", 50000, 0, 50);
  new TH1F("BHDT0_true", "BHD-T0 TOF", 50000, 0, 50);

  new TH1F("BLC1_eff_Kf", "BLC1 eff kaon beam", 2, 0, 2);
  new TH1F("BLC1_eff",    "BLC1 eff",           2, 0, 2);
  new TH1F("BLC1_time", "BLC1 time",        600,  -200,  400);
  new TH1F("BLC1_time_Kf", "BLC1 time",        600,  -200,  400);
  new TH1F("BLC1_time_true", "BLC1 time",        600,  -200,  400);
  new TH1F("BLC1_time_Kf_true", "BLC1 time",        600,  -200,  400);
  new TH1F("nBLC1_all", "BLC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC1_2",   "BLC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC1",     "BLC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC1_all_Kf", "BLC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC1_2_Kf",   "BLC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC1_Kf",     "BLC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("BLC1_chi2", "BLC1 chi-square",  1000,    0,  100);
  new TH1F("BLC1_chi2_Kf", "BLC1 chi-square",  1000,    0,  100);
  new TH1F("BLC1_chi2_true", "BLC1 chi-square",  1000,    0,  100);
  new TH1F("BLC1_chi2_Kf_true", "BLC1 chi-square",  1000,    0,  100);

  new TH1F("BLC2_eff_Kf", "BLC2 eff kaon beam", 2, 0, 2);
  new TH1F("BLC2_eff",    "BLC2 eff",           2, 0, 2);
  new TH1F("BLC2_time_all", "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time_2",   "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time",     "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time_all_Kf", "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time_2_Kf",   "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time_Kf",     "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time_true", "BLC2 time",        600,  -200,  400);
  new TH1F("BLC2_time_Kf_true", "BLC2 time",        600,  -200,  400);
  new TH1F("nBLC2_all", "BLC2 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC2_2", "BLC2 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC2", "BLC2 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC2_all_Kf", "BLC2 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC2_2_Kf", "BLC2 ntrack",     10,   -0.5,  9.5);
  new TH1F("nBLC2_Kf", "BLC2 ntrack",     10,   -0.5,  9.5);
  new TH1F("BLC2_chi2", "BLC2 chi-square",  1000,    0,  100);
  new TH1F("BLC2_chi2_Kf", "BLC2 chi-square",  1000,    0,  100);
  new TH1F("BLC2_chi2_true", "BLC2 chi-square",  1000,    0,  100);
  new TH1F("BLC2_chi2_Kf_true", "BLC2 chi-square",  1000,    0,  100);

  new TH1F("BPC_eff_Kf", "BPC eff kaon beam", 2, 0, 2);
  new TH1F("BPC_eff",    "BPC eff",           2, 0, 2);
  new TH1F("BPC_time", "BPC time",        600,  -200,  400);
  new TH1F("BPC_time_Kf", "BPC time",        600,  -200,  400);
  new TH1F("BPC_time_true", "BPC time",        600,  -200,  400);
  new TH1F("BPC_time_Kf_true", "BPC time",        600,  -200,  400);
  new TH1F("nBPC_all", "BPC ntrack",     10,   -0.5,  9.5);
  new TH1F("nBPC_2",   "BPC ntrack",     10,   -0.5,  9.5);
  new TH1F("nBPC",     "BPC ntrack",     10,   -0.5,  9.5);
  new TH1F("nBPC_all_Kf", "BPC ntrack",     10,   -0.5,  9.5);
  new TH1F("nBPC_2_Kf",   "BPC ntrack",     10,   -0.5,  9.5);
  new TH1F("nBPC_Kf",     "BPC ntrack",     10,   -0.5,  9.5);
  new TH1F("BPC_chi2", "BPC chi-square",  1000,    0,  100);
  new TH1F("BPC_chi2_Kf", "BPC chi-square",  1000,    0,  100);
  new TH1F("BPC_chi2_true", "BPC chi-square",  1000,    0,  100);
  new TH1F("BPC_chi2_Kf_true", "BPC chi-square",  1000,    0,  100);

  new TH1F("D5_chi2", "D5 chi2",  1000, 0.0, 100);
  new TH1F("D5_mom",  "Beam mom",  4000, 0.8, 1.2);
  new TH1F("D5_chi2_Kf", "D5 chi2",  1000, 0.0, 100);
  new TH1F("D5_mom_Kf",  "Beam mom",  4000, 0.8, 1.2);
  new TH1F("D5_chi2_true", "D5 chi2",  1000, 0.0, 100);
  new TH1F("D5_mom_true",  "Beam mom",  4000, 0.8, 1.2);
  new TH1F("D5_chi2_Kf_true", "D5 chi2",  1000, 0.0, 100);
  new TH1F("D5_mom_Kf_true",  "Beam mom",  4000, 0.8, 1.2);
  for( int seg=1; seg<=20; seg++ ){
    new TH1F(Form("D5_mom_BHD%d", seg), Form("Beam mom w/ BHD seg%d", seg), 4000, 0.8, 1.2);
    new TH1F(Form("D5_mom_BHD%d_Kf", seg), Form("Beam mom w/ BHD seg%d_Kf", seg), 4000, 0.8, 1.2);
  }

  new TH2F("BLC2BPC_diff",     "BLC2-BPC pos difference", 1000,  -10,  10, 1000,  -10,  10);
  new TH2F("BLC2BPC_dir_diff", "BLC2-BPC dir difference", 1000, -0.1, 0.1, 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_diff_X", "BLC2-BPC pos difference X", 1000, -10, 10);
  new TH1F("BLC2BPC_diff_Y", "BLC2-BPC pos difference Y", 1000, -10, 10);
  new TH1F("BLC2BPC_dir_diff_X", "BLC2-BPC dir difference X", 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_dir_diff_Y", "BLC2-BPC dir difference Y", 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_diff_X_true", "BLC2-BPC pos difference X", 1000, -10, 10);
  new TH1F("BLC2BPC_diff_Y_true", "BLC2-BPC pos difference Y", 1000, -10, 10);
  new TH1F("BLC2BPC_dir_diff_X_true", "BLC2-BPC dir difference X", 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_dir_diff_Y_true", "BLC2-BPC dir difference Y", 1000, -0.1, 0.1);

  new TH2F("BLC2BPC_diff_Kf",     "BLC2-BPC pos difference", 1000,  -10,  10, 1000,  -10,  10);
  new TH2F("BLC2BPC_dir_diff_Kf", "BLC2-BPC dir difference", 1000, -0.1, 0.1, 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_diff_X_Kf", "BLC2-BPC pos difference X", 1000, -10, 10);
  new TH1F("BLC2BPC_diff_Y_Kf", "BLC2-BPC pos difference Y", 1000, -10, 10);
  new TH1F("BLC2BPC_dir_diff_X_Kf", "BLC2-BPC dir difference X", 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_dir_diff_Y_Kf", "BLC2-BPC dir difference Y", 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_diff_X_Kf_true", "BLC2-BPC pos difference X", 1000, -10, 10);
  new TH1F("BLC2BPC_diff_Y_Kf_true", "BLC2-BPC pos difference Y", 1000, -10, 10);
  new TH1F("BLC2BPC_dir_diff_X_Kf_true", "BLC2-BPC dir difference X", 1000, -0.1, 0.1);
  new TH1F("BLC2BPC_dir_diff_Y_Kf_true", "BLC2-BPC dir difference Y", 1000, -0.1, 0.1);

  new TH2F("BeamProf", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_Kf", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_KCDH2", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_true", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_Kf_true", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);
  new TH2F("BeamProf_KCDH2_true", "Beam Profile at FF", 1000, -25, 25, 1000, -25, 25);

  new TH1F("FDC1_eff_Kf", "FDC1 eff kaon beam", 2, 0, 2);
  new TH1F("FDC1_eff",    "FDC1 eff",           2, 0, 2);
  new TH1F("FDC1_time", "FDC1 time",        600,  -200,  400);
  new TH1F("FDC1_time_Kf", "FDC1 time",        600,  -200,  400);
  new TH1F("FDC1_time_true", "FDC1 time",        600,  -200,  400);
  new TH1F("FDC1_time_Kf_true", "FDC1 time",        600,  -200,  400);
  new TH1F("nFDC1", "FDC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("nFDC1_Kf", "FDC1 ntrack",     10,   -0.5,  9.5);
  new TH1F("FDC1_chi2", "FDC1 chi-square",  1000,    0,  100);
  new TH1F("FDC1_chi2_Kf", "FDC1 chi-square",  1000,    0,  100);
  new TH1F("FDC1_chi2_true", "FDC1 chi-square",  1000,    0,  100);
  new TH1F("FDC1_chi2_Kf_true", "FDC1 chi-square",  1000,    0,  100);
}

void fillReadBeam_T0(EventHeader *header, BeamLineHitMan *blMan)
{
  MyHistTools::fillTH("nT0", MyTools::getHodo(blMan, CID_T0).size()+0.5);

  if( header->IsTrig(Trig_Kf) ){
    MyHistTools::fillTH("nT0_Kf", MyTools::getHodo(blMan, CID_T0).size()+0.5);
  }
}

void fillReadBeam_BHDT0(EventHeader *header, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);

  MyHistTools::fillTH("nBHD", beam->nBHD()+0.5);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("nBHD_Kf", beam->nBHD()+0.5);

  for( int i=0; i<beam->nBHD(); i++ ){
    int seg;
    double time;
    beam->getBHD(i, seg, time);
    double tof=beam->T0time()-time;

    MyHistTools::fillTH("BHDT0", tof);
    if( MyAnaTools::isTOFKaon(tof) ) MyHistTools::fillTH("BHDT0_true", tof);
    if( header->IsTrig(Trig_Kf) ){
      MyHistTools::fillTH("BHDT0_Kf", tof);
      if( MyAnaTools::isTOFKaon(tof) ) MyHistTools::fillTH("BHDT0_Kf_true", tof);
    }
  }
}
