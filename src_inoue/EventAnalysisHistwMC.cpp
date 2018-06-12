#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
//#include "HistManBeamAna.h"
#include "CDSTrackingMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "BeamSpectrometer.h"

#include "HistManwMC.h"
#include "HistTools.h"
#include "ProtonArm.h"

#define BLDC_CHI_CUT 1
#define BLDC_WIDE_TIMING 1

static const double BLC1_CHI2_MAX=10;
static const double BLC2_CHI2_MAX=10;
static const double BPC_CHI2_MAX=10;
//static const double BLC1_CHI2_MAX=30;
//static const double BLC2_CHI2_MAX=30;
//static const double BPC_CHI2_MAX=30;
//static const double BEAM_SPEC_CHI2_MAX=20;
static const double BEAM_SPEC_CHI2_MAX=30;

static const double TOF_K_MIN=27.8588;
static const double TOF_K_MAX=29.5663;
static const double TOF_PI_MIN=24.5606;
static const double TOF_PI_MAX=26.9866;
//static const double TOF_K_MIN=27;
//static const double TOF_K_MAX=31;
//static const double TOF_PI_MIN=24;
//static const double TOF_PI_MAX=27;

// BLC BPC consistency check off
static const double BLC2BPC_X_MIN=-0.72645;
static const double BLC2BPC_X_MAX=0.770232;
static const double BLC2BPC_Y_MIN=-0.778481;
static const double BLC2BPC_Y_MAX=0.755978;
static const double BLC2BPC_dX_MIN=-0.0202092;
static const double BLC2BPC_dX_MAX=0.0201483;
static const double BLC2BPC_dY_MIN=-0.0200048;
static const double BLC2BPC_dY_MAX=0.0204583;

static const double BHD_MATCH_MIN[10]={ 0.970692, 0.9754, 0.981274, 0.986733, 0.992045, 
					0.99713, 1.00219, 1.00708, 1.01191, 1.01588 };

static const double BHD_MATCH_MAX[10]={ 1.01324, 1.01846, 1.02366, 1.0282, 1.0321,
					1.03563, 1.03903, 1.0422, 1.04523, 1.04518 };

static const double CDH_Thre=4.0;
static const double IH_Thre=0.8; 

using namespace std;

class EventAnalysisHistwMC: public EventTemp
{
public:
  EventAnalysisHistwMC();
  ~EventAnalysisHistwMC();

private:
  int Event_Number_onSpill;
  int CDC_Event_Num;
  int t0, t1;
  TFile *cdcFile;
  TTree *cdcTree;
  CDSTrackingMan *cdstrackMan;
  EventHeader *cdcHeader;

  TTree *scaTree;

  TFile *rtFile;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;

  HistManwMC *histMan;

public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Clear();
  void Finalize();
};

EventAnalysisHistwMC::EventAnalysisHistwMC()
  : EventTemp()
{
  std::cout<<"EventAnalysisHistwMC::Constractor Call"<<std::endl;
}

EventAnalysisHistwMC::~EventAnalysisHistwMC()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisHistwMC::Initialize( ConfMan *conf )
{
  std::cout<<"EventAnalysisHistwMC::Initialization START"<<std::endl;
  Event_Number_onSpill = 0;
  CDC_Event_Num = 0;
  confMan = conf;
  ProtonArm::Initialize(confMan);

  cdcFile = new TFile((TString)confMan->GetCDSTrackFileName());
  cdstrackMan = new CDSTrackingMan();
  cdcHeader = new EventHeader();
  cdcTree = (TTree*)cdcFile-> Get("EventTree");
  cdcTree-> SetBranchAddress("EventHeader", &cdcHeader);
  cdcTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  scaTree = new TTree("ScalerTree", "ScalerTree");
  scaMan = new ScalerMan();
  scaTree-> Branch("ScalerMan", &scaMan);

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(confMan);

  histMan = new HistManwMC(rtFile, confMan);
  histMan-> setK18ana(blMan, cdsMan, bltrackMan, cdstrackMan);
  histMan-> initHist();

  InitializeHistogram();
  t0 = time(0);
  std::cout<<"EventAnalysisHistwMC::Initialization FINISH"<<std::endl;
}

void EventAnalysisHistwMC::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaTree-> Fill();
  scaMan->Clear();
}

bool EventAnalysisHistwMC::UAna( TKOHitCollection *tko )
{
  TH1F *h1;
  TH2F *h2;

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%5000==0 ){
    //  if( Event_Number%10==0 ){
    t1 = time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): "<<(t1-t0)<<std::endl;
  }
  TH1F *h1_ER = (TH1F*)gFile-> Get("EventReduction");
  h1_ER-> Fill(-1);
    
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );

  for( int i=0; i<20; i++ ){
    int tdc = header->pattern(i);
    if( -1<tdc && tdc<4095 ){
      h1 = (TH1F*)gFile-> Get(Form("trig_tdc_%d", i)), h1-> Fill(tdc);
      if( header-> IsTrig(i, confMan) ) h1 = (TH1F*)gFile-> Get(Form("trig_tdc_%d_eff", i)), h1-> Fill(tdc);
    }
  }
  h1 = (TH1F*)gFile-> Get("Trigger");
  h1-> Fill(0);
  for( int i=1; i<18; i++ ){ if( header->IsTrig(i) ) h1-> Fill(i); }

  h1 = (TH1F*)gFile-> Get("TrigMode");
  h1-> Fill(0);
  h1-> Fill(header->trigmode(confMan));

  h1 = (TH1F*)gFile-> Get("TrigMode1");
  h1-> Fill(0);
  for( int i=1; i<18; i++ ){ if( header->trigmode(i, confMan) ) h1-> Fill(i); }

  h1 = (TH1F*)gFile-> Get("TrigMode2");
  h1-> Fill(0);
  for( int i=1; i<18; i++ ){ if( header->trigmode2(i, confMan) ) h1-> Fill(i); }


  if( header->IsTrig(Trig_Cosmic) ){
    Clear();
    return true;
  }
  Event_Number_onSpill++;

  TH1F *h1_Kf = (TH1F*)gFile-> Get("Kf_Reduction");
  TH1F *h1_N  = (TH1F*)gFile-> Get("N_Reduction");

  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  fillHistIH(cdsMan);

  h1_ER-> Fill(0);
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(0);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(0);

  //  if( header->IsTrig(Trig_Beam ) ) std::cout<<" Trig Beam "<<std::endl;

  //*** Check T0 1hit ***//
  int nT0=0;
  double T0time = DBL_MIN;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      T0time = blMan->T0(i)->ctmean();
    }
  }
  if( header->IsTrig(Trig_Kf) ){
    h1 = (TH1F*)gFile-> Get("nT0"), h1->Fill(nT0);
  }

  if( nT0==1 ){
    h1_ER-> Fill(1);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(1);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(1);
  }

  //*** for Trigger Efficiency ***//
  if( header->IsTrig(Trig_Kf) ){
    int tmp_nCDH=0;
    for( int i=0; i<cdsMan->nCDH(); i++ ){
      if( cdsMan->CDH(i)->CheckRange() ){
	tmp_nCDH++;
      }
    }
    int tmp_nNC=0;
    for( int i=0; i<blMan->nNC(); i++ ){
      if( blMan->NC(i)->CheckRange() ){
	tmp_nNC++;
      }
    }
    int tmp_nCVC=0;
    for( int i=0; i<blMan->nCVC(); i++ ){
      if( blMan->CVC(i)->CheckRange() ){
	tmp_nCVC++;
      }
    }
    int tmp_nBVC=0;
    for( int i=0; i<blMan->nBVC(); i++ ){
      if( blMan->BVC(i)->CheckRange() ){
	tmp_nBVC++;
      }
    }

    h1 = (TH1F*)gFile-> Get("TrigEff_CDH2");
    h1-> Fill(0);
    if( tmp_nCDH>1 ){
      h1-> Fill(1);
      if( header->IsTrig(Trig_KCDH2, confMan) ){ h1-> Fill(2); }
      if( header->IsTrig(Trig_1stMix, confMan) ){ h1-> Fill(3); }
      if( header->IsTrig(Trig_KCDH2, confMan ) && header->IsTrig(Trig_1stMix, confMan)){ h1-> Fill(4); }
    }

    if( tmp_nCDH>0 ){
      h1 = (TH1F*)gFile-> Get("TrigEff_N");
      h1-> Fill(0);
      if( header-> IsTrig(Trig_KCDH1) ){ h1-> Fill(1); }
      if( tmp_nNC>0 && tmp_nCVC==0 && tmp_nBVC==0 ){
	if( header-> IsTrig(Trig_KCDH1, confMan) ){ h1-> Fill(2); };
	if( header-> IsTrig(Trig_Neutral, confMan) ){ h1-> Fill(3); };
	if( header-> IsTrig(Mode_KCDH1N, confMan ) && header->IsTrig(Trig_Neutral, confMan) ){ h1->Fill(4); }
      }
    }
  }

  //******************************************************************************************************************************************//
  //*** return T0!=1 event                                                                                                                 ***//
  //******************************************************************************************************************************************//
  if( nT0!=1 ){
    Clear();
    return true;
  }
  //*** Check BHD-T0 TOF ***//
  bool kaon_flag = false;
  bool pion_flag = false;
  int nBHD=0;
  int nBHD_woOut=0;
  int BHDseg=-1;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      nBHD++;
      if( header->IsTrig(Trig_Kf) ){
	h1 = (TH1F*)gFile-> Get("hitpatBHD"), h1->Fill(blMan->BHD(i)->seg());
      }
      if( 5<blMan->BHD(i)->seg() && blMan-> BHD(i)->seg()<16 ){
	BHDseg = blMan->BHD(i)->seg();
	nBHD_woOut++;
	if( nT0==1 ){
	  double tof = T0time-blMan->BHD(i)->ctmean();
	  h1 = (TH1F*)gFile-> Get("BHDT0"), h1-> Fill(tof);
	  if( header->IsTrig(Trig_Kf) ) h1 = (TH1F*)gFile-> Get("BHDT0_Kf"), h1-> Fill(tof);
	  if( header->IsTrig(Trig_Kaon) ) h1 = (TH1F*)gFile-> Get("BHDT0_K"), h1-> Fill(tof);
	  if( header->IsTrig(Trig_Pion) ) h1 = (TH1F*)gFile-> Get("BHDT0_pi"), h1-> Fill(tof);
	  if( TOF_K_MIN<tof && tof<TOF_K_MAX ) kaon_flag = true;
	  if( TOF_PI_MIN<tof && tof<TOF_PI_MAX ) pion_flag = true;
	}
      }
    }
  }
  if( header->IsTrig(Trig_Kf) ){
    h1 = (TH1F*)gFile-> Get("nBHD"), h1->Fill(nBHD);
    h1 = (TH1F*)gFile-> Get("nBHD_woOut"), h1->Fill(nBHD_woOut);
  }

  histMan->setBeamKaon(kaon_flag);
  histMan->setBeamPion(pion_flag);
  int beam_pid=Beam_Other;
  if( header->IsTrig(Trig_Kaon) ){
    if( kaon_flag ) beam_pid=Beam_Kaon;
  } 
  else if( header->IsTrig(Trig_Pion) ){
    if( pion_flag ) beam_pid=Beam_Pion;
  }

  histMan->setBeamPID(beam_pid);
  if( beam_pid==Beam_Kaon ){
    h1_ER-> Fill(2);
    if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(2);
    if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(2);
  }
  histMan->setBeamPID(beam_pid);

  //  bltrackMan-> DoTracking2(blMan, confMan);
  bltrackMan-> DoTracking(blMan, confMan, true, true);
  // cout<<" ntrack : "<<bltrackMan->ntrackBPC()<<"  status : "<<bltrackMan->status(CID_BPC)<<endl;
  // for( int l=1; l<=8; l++ ){
  //   cout<<" n BLC lay"<<l<<" nhit : "<<blMan->nBPC(l)<<endl;
  // }

  bool BLC1_flag=false;
  bool BLC2_flag=false;
  bool BPC_flag=false;
  bool beam_mom_flag = false;
  bool BLC2BPC_flag = false;
  bool beam_fiducial_flag = false;

  int ntrackBLC1=0;
  int ntrackBLC1_1=0;
  int ntrackBLC1_2=0;
  LocalTrack *trackBLC1=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double tracktime = bltrackMan->trackBLC1(i)->GetTrackTime();
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC1_time"), h1-> Fill(tracktime);
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC1_time_pi"), h1-> Fill(tracktime);
    if( -50<tracktime && tracktime<200 ) ntrackBLC1_1++;
    if( -30<tracktime && tracktime<100 ) ntrackBLC1_2++;
    if( -10<tracktime && tracktime<10 ){
      ntrackBLC1++;
      trackBLC1 = bltrackMan->trackBLC1(i);
    }
  }
  if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC1"), h1-> Fill(bltrackMan->ntrackBLC1());
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_1"), h1-> Fill(ntrackBLC1_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_2"), h1-> Fill(ntrackBLC1_2);
  }
  else if( beam_pid==Beam_Pion ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_pi"), h1-> Fill(bltrackMan->ntrackBLC1());
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_1_pi"), h1-> Fill(ntrackBLC1_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC1_2_pi"), h1-> Fill(ntrackBLC1_2);
  }
  if( trackBLC1 ){
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC1_chi2"), h1-> Fill(trackBLC1->chi2all());
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC1_chi2_pi"), h1-> Fill(trackBLC1->chi2all());
  }
  if( ntrackBLC1_2==1 && ntrackBLC1==1 ){
    if( trackBLC1->chi2all()<BLC1_CHI2_MAX ){
      BLC1_flag=true;
    }
  }
  if( header->IsTrig(Trig_Kf) ){
    for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
      LocalTrack *track = bltrackMan-> trackBLC1(i);
      TVector3 local_pos = track-> GetPosatZ(track->hit(i)->gz());
      h2 = (TH2F*)gFile-> Get("hitposBLC1"), h2-> Fill(local_pos.X(), local_pos.Y());
    }
  }

  int ntrackBLC2=0;
  int ntrackBLC2_1=0;
  int ntrackBLC2_2=0;
  LocalTrack *trackBLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double tracktime = bltrackMan->trackBLC2(i)->GetTrackTime();
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC2_time"), h1-> Fill(tracktime);
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC2_time_pi"), h1-> Fill(tracktime);
    if( -50<tracktime && tracktime<200 ) ntrackBLC2_1++;
    if( -30<tracktime && tracktime<100 ) ntrackBLC2_2++;
    if( -10<tracktime && tracktime<10 ){
      ntrackBLC2++;
      trackBLC2 = bltrackMan->trackBLC2(i);
    }
  }
  if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC2"), h1-> Fill(bltrackMan->ntrackBLC2());
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_1"), h1-> Fill(ntrackBLC2_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_2"), h1-> Fill(ntrackBLC2_2);
  }
  else if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_pi"), h1-> Fill(bltrackMan->ntrackBLC2());
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_1_pi"), h1-> Fill(ntrackBLC2_1);
    h1 = (TH1F*)gFile-> Get("ntrackBLC2_2_pi"), h1-> Fill(ntrackBLC2_2);
  }
  if( trackBLC2 ){
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BLC2_chi2"), h1-> Fill(trackBLC2->chi2all());
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BLC2_chi2_pi"), h1-> Fill(trackBLC2->chi2all());
  }
  if( ntrackBLC2_2==1 && ntrackBLC2==1 ){
    if( trackBLC2->chi2all()<BLC2_CHI2_MAX ){
      BLC2_flag=true;
    }
  }
  if( header->IsTrig(Trig_Kf) ){
    for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
      LocalTrack *track = bltrackMan-> trackBLC2(i);
      TVector3 local_pos = track-> GetPosatZ(track->hit(i)->gz());
      h2 = (TH2F*)gFile-> Get("hitposBLC2"), h2-> Fill(local_pos.X(), local_pos.Y());
    }
  }

  double D5mom=-999;
  if( trackBLC1 && trackBLC2 ){
    beamSpec->TMinuitFit(trackBLC1, trackBLC2, confMan);
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("D5_chi2"), h1-> Fill(beamSpec->chisquare());
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("D5_chi2"), h1-> Fill(beamSpec->chisquare());
    D5mom = beamSpec->mom();

    if( nBHD==1 && BHDseg>0 ){
      h1 = (TH1F*)gFile-> Get(Form("D5_mom_BHD%d", BHDseg)), h1-> Fill(D5mom);
    }
    if( beam_pid==Beam_Kaon )  h1 = (TH1F*)gFile-> Get("D5_mom"), h1->Fill(beamSpec->mom());
    else if( beam_pid==Beam_Pion ) h1= (TH1F*)gFile-> Get("D5_mom_pi"), h1-> Fill(beamSpec->mom());
    if( beamSpec->chisquare()<BEAM_SPEC_CHI2_MAX ){
      //***** off BLC1 & BHD matching *****//
      //      beam_mom_flag = true;
      for( int i=0; i<blMan->nBHD(); i++ ){
	if( blMan->BHD(i)->CheckRange() ){
	  if( 5<blMan->BHD(i)->seg() && blMan-> BHD(i)->seg()<16 ){
	    int seg = blMan->BHD(i)->seg();
	    if( BHD_MATCH_MIN[seg-6]<D5mom && D5mom<BHD_MATCH_MAX[seg-6] ) beam_mom_flag=true;
	  }
	}
      }
    }
  }

  int ntrackBPC=0;
  int ntrackBPC_1=0;
  int ntrackBPC_2=0;
  LocalTrack *trackBPC=0;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    double tracktime = bltrackMan->trackBPC(i)->GetTrackTime();
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BPC_time"), h1-> Fill(tracktime);

    if( -50<tracktime && tracktime<200 ) ntrackBPC_1++;
    if( -30<tracktime && tracktime<100 ) ntrackBPC_2++;
    if( -10<tracktime && tracktime<10 ){
      ntrackBPC++;
      trackBPC = bltrackMan->trackBPC(i);
    }
  }
  if( beam_pid==Beam_Kaon ){
    h1 = (TH1F*)gFile-> Get("ntrackBPC"), h1-> Fill(bltrackMan->ntrackBPC());
    h1 = (TH1F*)gFile-> Get("ntrackBPC_1"), h1-> Fill(ntrackBPC_1);
    h1 = (TH1F*)gFile-> Get("ntrackBPC_2"), h1-> Fill(ntrackBPC_2);
  }
  else if( beam_pid==Beam_Pion ){
    h1 = (TH1F*)gFile-> Get("ntrackBPC_pi"), h1-> Fill(bltrackMan->ntrackBPC());
    h1 = (TH1F*)gFile-> Get("ntrackBPC_1_pi"), h1-> Fill(ntrackBPC_1);
    h1 = (TH1F*)gFile-> Get("ntrackBPC_2_pi"), h1-> Fill(ntrackBPC_2);
  }

  if( trackBPC ){
    if( beam_pid==Beam_Kaon ) h1 = (TH1F*)gFile-> Get("BPC_chi2"), h1-> Fill(trackBPC->chi2all());
    else if( beam_pid==Beam_Pion ) h1 = (TH1F*)gFile-> Get("BPC_chi2_pi"), h1-> Fill(trackBPC->chi2all());
  }
  if( ntrackBPC_2==1 && ntrackBPC==1 ){
    if( trackBPC->chi2all()<BPC_CHI2_MAX ){
      BPC_flag=true;
    }
  }

  TVector3 BLC2BPC_diff;
  TVector3 BLC2BPC_dir_diff;

  if( trackBLC2 && trackBPC ){
    double center_z = 0.5*(-130-20.3);
    TVector3 BLC2_track_pos = trackBLC2-> GetPosatZ(center_z);
    TVector3 BPC_track_pos  = trackBPC-> GetPosatZ(center_z);
    TVector3 BLC2_mom_dir = trackBLC2-> GetMomDir();
    TVector3 BPC_mom_dir = trackBPC-> GetMomDir();
    if( beam_pid==Beam_Kaon ){
      h2 = (TH2F*)rtFile-> Get("BLC2BPC_diff"), h2-> Fill((BLC2_track_pos-BPC_track_pos).X(), (BLC2_track_pos-BPC_track_pos).Y());
      h2 = (TH2F*)rtFile-> Get("BLC2BPC_dir_diff"), h2-> Fill((BLC2_mom_dir-BPC_mom_dir).X(), (BLC2_mom_dir-BPC_mom_dir).Y());
    } else if( beam_pid==Beam_Pion ){
      h2 = (TH2F*)rtFile-> Get("BLC2BPC_diff_pi"), h2-> Fill((BLC2_track_pos-BPC_track_pos).X(), (BLC2_track_pos-BPC_track_pos).Y());
      h2 = (TH2F*)rtFile-> Get("BLC2BPC_dir_diff_pi"), h2-> Fill((BLC2_mom_dir-BPC_mom_dir).X(), (BLC2_mom_dir-BPC_mom_dir).Y());
    }

    BLC2BPC_diff = BLC2_track_pos-BPC_track_pos;
    BLC2BPC_dir_diff = BLC2_mom_dir-BPC_mom_dir;
  }

  //***** off BLC2 & BPC matching *****//
  //  BLC2BPC_flag = true;
  if( BLC2BPC_diff.X()>BLC2BPC_X_MIN && BLC2BPC_X_MAX>BLC2BPC_diff.X() && BLC2BPC_diff.Y()>BLC2BPC_Y_MIN && BLC2BPC_Y_MAX>BLC2BPC_diff.Y() &&
      BLC2BPC_dir_diff.X()>BLC2BPC_dX_MIN && BLC2BPC_dX_MAX>BLC2BPC_dir_diff.X() && BLC2BPC_dir_diff.Y()>BLC2BPC_dY_MIN && BLC2BPC_dY_MAX>BLC2BPC_dir_diff.Y() ){
    BLC2BPC_flag = true;
  }  

  // if( bltrackMan->ntrackBLC1a()==1 ) fillBLC1a(bltrackMan->trackBLC1a(0));
  // if( bltrackMan->ntrackBLC1b()==1 ) fillBLC1b(bltrackMan->trackBLC1b(0));
  // if( bltrackMan->ntrackBLC2a()==1 ) fillBLC2a(bltrackMan->trackBLC2a(0));
  // if( bltrackMan->ntrackBLC2b()==1 ) fillBLC2b(bltrackMan->trackBLC2b(0));
  // if( trackBPC ){
  //   fillBPC(trackBPC);
  // }
  TVector3 BPC_FF;

  if( trackBPC && BLC2BPC_flag ){
    BPC_FF = trackBPC-> GetPosatZ(-3);

    if( beam_pid==Beam_Kaon ) h1= (TH1F*)gFile->Get("BeamProf"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
    else if( beam_pid==Beam_Pion ) h1= (TH1F*)gFile->Get("BeamProf_pi"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
    
    if( header->IsTrig(Trig_KCDH2) ){
      h1= (TH1F*)gFile->Get("BeamProf_CDH2"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
    }

    if( header->IsTrig(Trig_Kf) ){
      h1= (TH1F*)gFile->Get("BeamProf_Kf"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
    }
    if( beam_pid==Beam_Kaon){
      if( GeomTools::GetID(BPC_FF)==CID_Fiducial ) beam_fiducial_flag=true;
    }
  }

  if( beam_pid==Beam_Kaon ){
    if( BPC_flag ){
      h1_ER->Fill(3);
      if( header->IsTrig(Trig_Kf) ) h1_Kf->Fill(3);
      if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(3);
      if( BLC1_flag ){
	h1_ER->Fill(4);
	if( header->IsTrig(Trig_Kf) ) h1_Kf->Fill(4);
	if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(4);
	if( BLC2_flag ){
	  h1_ER->Fill(5);
	  if( header->IsTrig(Trig_Kf) ) h1_Kf->Fill(5);
	  if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(5);
	  if( BLC2BPC_flag ){
	    h1_ER->Fill(6);
	    if( header->IsTrig(Trig_Kf) ) h1_Kf->Fill(6);
	    if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(6);
	    if( beam_mom_flag ){
	      h1_ER->Fill(7);
	      if( header->IsTrig(Trig_Kf) ) h1_Kf->Fill(7);
	      if( header->IsTrig(Trig_Neutral) ) h1_N->Fill(7);
	      if( beam_fiducial_flag ){
		if( header->IsTrig(Trig_Kf) ) h1_Kf->Fill(8);
	      }
	    }
	  }
	}
      }
    }
  }

  if( BLC2_flag ){
    TVector3 posAC = trackBLC2-> GetPosatZ(-87.2);
    h2 = (TH2F*)rtFile-> Get("hitposAC"), h2-> Fill(posAC.X(), posAC.Y());
    if( header-> IsTrig(Trig_Kf) )  h2 = (TH2F*)rtFile-> Get("hitposAC_Kf"), h2-> Fill(posAC.X(), posAC.Y());
    if( header-> IsTrig(Trig_Kaon) )  h2 = (TH2F*)rtFile-> Get("hitposAC_Kaon"), h2-> Fill(posAC.X(), posAC.Y());
    if( header-> IsTrig(Trig_Pion) )  h2 = (TH2F*)rtFile-> Get("hitposAC_Pion"), h2-> Fill(posAC.X(), posAC.Y());
  }

  if( beam_pid==Beam_Kaon ){
    if( BLC2_flag && BPC_flag && BLC2BPC_flag ){
      if( bltrackMan-> ntrackBLC1a()==1 ){
	LocalTrack *track = bltrackMan->trackBLC1a(0);
	double dx = track->dx(), dy = track->dy();
	double x, y, z=15.;
	track-> XYLocalPosatZ(x, y, z);
	double time = track->GetTrackTime();
	// if( -0.001<dx && dx<0.001 && -0.001<dy && dy<0.001 && track->chi2all()<10. && -10.<time && time<10 ){
	//   if( -8.<x && x<8. && -.8<y && y<8. ){
	//     h1 = (TH1F*)gFile-> Get("BLC1b_eff"), h1->Fill(0);
	//     if( bltrackMan->ntrackBLC1b()>0 ) h1->Fill(1);
	//   }
	// }
      }
      if( bltrackMan-> ntrackBLC1b()==1 ){
	LocalTrack *track = bltrackMan->trackBLC1b(0);
	double dx = track->dx(), dy = track->dy();
	double x, y, z=-15.;
	track-> XYLocalPosatZ(x, y, z);
	double time = track->GetTrackTime();
	// if( -0.001<dx && dx<0.001 && -0.001<dy && dy<0.001 && track->chi2all()<10. && -10.<time && time<10 ){
	//   if( -8.<x && x<8. && -.8<y && y<8. ){
	//     h1 = (TH1F*)gFile-> Get("BLC1a_eff"), h1->Fill(0);
	//     if( bltrackMan->ntrackBLC1a()>0 ) h1->Fill(1);
	//   }
	// }
      }
    }

    if( BLC1_flag && BPC_flag ){
      if( bltrackMan-> ntrackBLC2a()==1 ){
	LocalTrack *track = bltrackMan->trackBLC2a(0);
	double dx = track->dx(), dy = track->dy();
	double x, y, z=-15.;
	track-> XYLocalPosatZ(x, y, z);
	double time = track->GetTrackTime();
	// if( -0.001<dx && dx<0.001 && -0.001<dy && dy<0.001 && track->chi2all()<10. && -10.<time && time<10 ){
	//   if( -8.<x && x<8. && -.8<y && y<8. ){
	//     h1 = (TH1F*)gFile-> Get("BLC2b_eff"), h1->Fill(0);
	//     if( bltrackMan->ntrackBLC2b()>0 ) h1->Fill(1);
	//   }
	// }
      }
      if( bltrackMan-> ntrackBLC2b()==1 ){
	LocalTrack *track = bltrackMan->trackBLC2b(0);
	double dx = track->dx(), dy = track->dy();
	double x, y, z=15.;
	track-> XYLocalPosatZ(x, y, z);
	double time = track->GetTrackTime();
	// if( -0.001<dx && dx<0.001 && -0.001<dy && dy<0.001 && track->chi2all()<10. && -10.<time && time<10 ){
	//   if( -8.<x && x<8. && -.8<y && y<8. ){
	//     h1 = (TH1F*)gFile-> Get("BLC2a_eff"), h1->Fill(0);
	//     if( bltrackMan->ntrackBLC2a()>0 ) h1->Fill(1);
	//   }
	// }
      }
    }
  }

  if( BLC1_flag && beam_mom_flag && BLC2_flag ){
    if( bltrackMan->ntrackBLC2()==1 ){
      LocalTrack *track = bltrackMan-> trackBLC2(0);
      TVector3 dir = track-> GetMomDir();
      TVector3 pos = track-> GetPosatZ(-20);
      double time = track-> GetTrackTime();
      // if( -0.001<dir.X() && dir.X()<0.001 && -0.001<dir.Y() && dir.Y()<0.001 && track->chi2all()<10 && -5.<time && time<5 ){
      // 	if( beam_pid==Beam_Kaon ){
      // 	  h1 = (TH1F*)gFile-> Get("BPC_eff_k"), h1->Fill(0);
      // 	  if( bltrackMan->ntrackBPC()>0 ) h1->Fill(1);
      // 	}
      // 	else if( beam_pid==Beam_Pion ){
      // 	  h1 = (TH1F*)gFile-> Get("BPC_eff_pi"), h1->Fill(0);
      // 	  if( bltrackMan->ntrackBPC()>0 ) h1->Fill(1);
      // 	}
      // }
    }
  }

  //******************************************************************************************************************************************//
  //*** return BLDC condition                                                                                                              ***//
  //******************************************************************************************************************************************//
  if( !BLC1_flag || !BLC2_flag || !BPC_flag || !BLC2BPC_flag || !beam_mom_flag ){
    Clear();
    return true;
  }

  //*******************************************//
  //*** for CDC Tracking Efficeincy Trigger ***//
  //*******************************************//
  int nCDH=0, nIH=0;
  HodoscopeLikeHit *IH_hit=0, *CDH_hit=0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->CDH(i);
    int seg=hit->seg();
    double time = hit-> ctmean();
    double dE = hit->emean();
    if( hit->CheckRange() ){
      CDH_hit = hit;
      nCDH++;
    }
  }
  for( int i=0; i<cdsMan->nIH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan-> IH(i);
    double time = hit->tu();
    double dE = hit->eu();
    int adc = hit->adcu();
    int tdc = hit->tdcu();
    if( -1<tdc && tdc<4095 ){
      IH_hit = hit;
      nIH++;
    }
  }

  bool CDC_trig_flag=false;
  if( nIH==1 && nCDH==1 && beam_fiducial_flag && beam_pid==Beam_Kaon ){
    if( 0<IH_hit->tu() && IH_hit->tu()<10 && IH_hit->eu()>IH_Thre && 0<CDH_hit->ctmean() && CDH_hit->ctmean()<20 && CDH_hit->emean()>CDH_Thre ){
      TVector3 CDHpos = CDH_hit->pos();
      TVector3 IHpos = IH_hit->pos();
      if( fabs(CDHpos.Phi()-IHpos.Phi())<20.*TwoPi/360. || fabs(CDHpos.Phi()-TwoPi-IHpos.Phi())<20.*TwoPi/360. || fabs(CDHpos.Phi()+TwoPi-IHpos.Phi())<20.*TwoPi/360 ){
	// if( IH_hit->seg()==11 ){
	// 	cout<<"  CDH seg"<<CDH_hit->seg()<<" : "<<CDHpos.Phi()<<"  IH : "<<IHpos.Phi()<<endl;
	// }
	CDC_trig_flag=true;
	h1 = (TH1F*)gFile-> Get("CDC_trig"), h1-> Fill(CDH_hit->seg());
	h1 = (TH1F*)gFile-> Get("CDC_IH_trig"), h1-> Fill(IH_hit->seg());
      }
    }
  }

  //**************************//
  //*** CDC raw data Fill  ***//
  //**************************//
  for( int l=1; l<=NumOfCDCLayers; l++ ){
    for( int i=0; i<cdsMan->nCDC(l); i++ ){
      CDCHit *hit = cdsMan->CDC(l, i);
      int wire=hit->wire();
      h1 = (TH1F*)gFile-> Get(Form("hitpatCDC_%d", l)), h1-> Fill(wire);
    }
  }

  //****************************************//
  //*** CDC Tracking File Event Matching ***//
  //****************************************//
  cdcTree-> GetEntry(CDC_Event_Num);
  if( CDC_Event_Num>cdcTree->GetEntries() ){
    Clear();
    std::cout<<"  CDC File Reach End"<<std::endl;
    return true;
  }
  while( cdcHeader->ev()<Event_Number ){
    CDC_Event_Num++;
    if( CDC_Event_Num>cdcTree->GetEntries() ){
      Clear();
      return true;
    }
    cdcTree-> GetEntry(CDC_Event_Num);
  }
  if( cdcHeader->ev()>Event_Number ){
    Clear();
    return true;
  }
  if( cdcHeader->ev()!=Event_Number ){
    std::cout<<"  !!! CDC File Event matching miss !!!"<<std::endl;
    Clear();
    return false;
  }

  //*** dumping ******//
  // std::cout<<" CDC File match Event Number : "<<Event_Number<<std::endl;
  // std::cout<<"> nT0 : "<<nT0<<std::endl;
  // std::cout<<"> nT0 : "<<nT0<<std::endl;
  // std::cout<<"> Kaon : "<<std::boolalpha<<kaon_flag<<std::endl;
  // std::cout<<"> BLC1 : "<<std::boolalpha<<BLC1_flag<<"  BLC2 : "<<std::boolalpha<<BLC2_flag<<std::endl;
  // std::cout<<"> Beam Mom : "<<std::boolalpha<<beam_mom_flag<<std::endl;
  // std::string str; std::cin>>str; if( str=="q" ) exit(0);

  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *track= cdstrackMan->GoodTrack(i);
    for( int i=0; i<track->nCDHHit(); i++ ){
      HodoscopeLikeHit *hit = track->CDHHit(cdsMan, i);
      int seg= hit-> seg();
      double dE = hit->emean();
      h1 = (TH1F*)gFile-> Get(Form("CDH_dE_wCDC_%d", seg)), h1-> Fill(dE);
    }
  
    if( track->Chi()<30. ){
      track-> SearchHodoHit(cdsMan, confMan);
      TVector3 vtxIH  = track->IHVertex();
      TVector3 vtxCDH = track->CDHVertex();
      // std::cout<<" IH  dphi : "<<ihmaxdphi<<"  vtx("<<vtxIH.x()<<", "<<vtxIH.y()<<", "<<vtxIH.z()<<")"<<std::endl;
      // std::cout<<" CDH dphi : "<<cdhmaxdphi<<"  vtx("<<vtxCDH.x()<<", "<<vtxCDH.y()<<", "<<vtxCDH.z()<<")"<<std::endl;
      bool single=true;
      for( int lay=1; lay<=NumOfCDCLayers; lay++ ){ if( track->nTrackHit(lay)!=1 ) false; }
      //      cout<<" CDC singl : "<<boolalpha<<single<<endl;

      double dis;
      TVector3 hpos, lpos;
      bool vtx_flag = TrackTools::CalcLineHelixVertex( trackBPC, track, lpos, hpos, dis);
      if( dis>3.0 ) vtx_flag=false;
      bool fiducial_flag = false;
      if( GeomTools::GetID(lpos)==CID_Fiducial ) fiducial_flag=true;
      //      cout<<" vertex : "<<boolalpha<<vtx_flag<<"  dis="<<dis<<endl;
    
      if( vtx_flag && single ){
	//	cout<<" phi CDH : "<<phiCDH<<endl
	if( header->IsTrig(Trig_Kf) ){
	  int trig_seg=0;
	  for( int seg=1; seg<=36; seg++ ){
	    TVector3 pos;
	    confMan->GetGeomMapManager()->GetPos(CID_CDH, seg, pos);
	    double dphi=pos.Phi()-vtxCDH.Phi();
	    if( dphi>TMath::Pi()       ) dphi -= TMath::Pi()*2;
	    else if( dphi<-TMath::Pi() ) dphi += TMath::Pi()*2;
	    dphi *= TMath::RadToDeg();

	    if( fabs(dphi)<4.5 && fabs(vtxCDH.Z())<389.5 ){
	      trig_seg=seg;
	    }
	  }

	  if( trig_seg>0 ){
	    h1 = (TH1F*)rtFile-> Get("CDH_trig"), h1-> Fill(trig_seg);
	    bool flag=false;
	    for( int j=0; j<cdsMan->nCDH(); j++ ){
	      HodoscopeLikeHit *hit = cdsMan->CDH(j);
	      if( hit->seg()!=trig_seg ) continue;
	      if( hit->CheckRange() ) flag=true;
	    }
	    if( flag ) h1 = (TH1F*)rtFile-> Get("CDH_eff"), h1-> Fill(trig_seg);
	    // cout<<"  CDH trigger event  seg : "<<trig_seg<<" hit:"<<boolalpha<<flag<<endl;
	  }
	}

	if( fiducial_flag ){	
	  int trig_seg=0;
	  for( int seg=1; seg<=24; seg++ ){
	    TVector3 pos;
	    confMan->GetGeomMapManager()->GetGPos(CID_IH, seg, pos);
	    double dphi=pos.Phi()-vtxIH.Phi();
	    if( dphi>TMath::Pi()       ) dphi -= TMath::Pi()*2;
	    else if( dphi<-TMath::Pi() ) dphi += TMath::Pi()*2;
	    dphi *= TMath::RadToDeg();
	    
	    if( fabs(dphi)<7.0 ) trig_seg=seg;
	  }
	  if( trig_seg>0 ){
	    h1 = (TH1F*)rtFile-> Get("IH_trig"), h1-> Fill(trig_seg);
	    bool flag=false;
	    for( int j=0; j<cdsMan->nIH(); j++ ){
	      HodoscopeLikeHit *hit = cdsMan->IH(j);
	      if( hit->seg()!=trig_seg ) continue;
	      TVector3 dpos=hit->pos();
	      if( 0<hit->tdcu() && hit->tdcu()<4000 ) flag=true;
	    }
	    if( flag ) h1 = (TH1F*)rtFile-> Get("IH_eff"), h1-> Fill(trig_seg);
	    //	  cout<<"  IH trigger event  seg : "<<trig_seg<<" hit:"<<boolalpha<<flag<<endl;
	  }
	}
      }

      double diff_CDH = DBL_MAX;
      for( int i=0; i<cdsMan-> nCDH(); i++ ){
	HodoscopeLikeHit *hit = cdsMan->CDH(i);
	if( hit-> CheckRange() ){
	  if( fabs(vtxCDH.DeltaPhi(hit->pos()))<diff_CDH ){
	    diff_CDH = fabs(vtxCDH.DeltaPhi(hit->pos()));
	    CDH_hit = hit;
	  }
	}
      }

      double diff_IH = DBL_MAX;
      for( int i=0; i<cdsMan-> nIH(); i++ ){
	HodoscopeLikeHit *hit = cdsMan->IH(i);
	if( -1<hit-> tdcu() && hit->tdcu()<4095 ){
	  if( fabs(vtxIH.DeltaPhi(hit->pos()))<diff_IH ){
	    diff_IH = fabs(vtxIH.DeltaPhi(hit->pos()));
	    IH_hit = hit;
	  }
	}
      }

      if( CDH_hit ){
	h1 = (TH1F*)gFile-> Get(Form("CDH_deg_diff_%d", CDH_hit->seg())), h1-> Fill(360.*vtxCDH.DeltaPhi(CDH_hit->pos())/TwoPi);
      }

      if( IH_hit ){
	h1 = (TH1F*)gFile-> Get(Form("IH_deg_diff_%d", IH_hit->seg())), h1-> Fill(360.*vtxIH.DeltaPhi(IH_hit->pos())/TwoPi);
	// if( IH_hit->seg()==13 ){
	//   if( 360.*vtxIH.DeltaPhi(IH_hit->pos())/TwoPi>18. ){
	//     TVector3 pos = IH_hit->pos();
	//     std::cout<<" ev : "<<Event_Number<<std::endl;
	//     std::cout<<"    deg diff : "<<vtxIH.DeltaPhi(IH_hit->pos())/TwoPi<<std::endl;
	//     std::cout<<"  vtx IH : ("<<vtxIH.x()<<", "<<vtxIH.y()<<", "<<vtxIH.z()<<")"<<std::endl;
	//     std::cout<<"  IH pos : ("<<pos.x()<<", "<<pos.y()<<", "<<pos.z()<<")"<<std::endl;
	//   }
	// }
      }
    }
  }

  histMan-> setT0time(T0time);
  histMan-> setD5mom(D5mom);
  histMan-> setTrackBPC(trackBPC);
  histMan-> ana(header);

  if( CDC_trig_flag ){
    bool eff_flag=false, eff_flag2=false;
    for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
      CDSTrack *track = cdstrackMan->GoodTrack(i);
      bool IH_flag=false, CDH_flag=false;
      TVector3 vtxIH=track->IHVertex();
      TVector3 vtxCDH=track->CDHVertex();
      TVector3 posIH=IH_hit->pos();
      TVector3 posCDH=CDH_hit->pos();

      if( fabs(vtxIH.Phi()-posIH.Phi())<15.*TwoPi/360. || fabs(vtxIH.Phi()-TwoPi-posIH.Phi())<15.*TwoPi/360. ||  fabs(vtxIH.Phi()+TwoPi-posIH.Phi())<15.*TwoPi/360. ) IH_flag=true;
      if( fabs(vtxCDH.Phi()-posCDH.Phi())<10.*TwoPi/360. || fabs(vtxCDH.Phi()-TwoPi-posCDH.Phi())<10.*TwoPi/360. ||  fabs(vtxCDH.Phi()+TwoPi-posCDH.Phi())<10.*TwoPi/360. ) CDH_flag=true;

      // for( int j=0; j<track->nCDHHit(); j++ ){
      // 	HodoscopeLikeHit *hit = track->CDHHit(cdsMan, j);
      // 	if( hit->seg()==CDH_hit->seg() ) CDH_flag = true;
      // }
      // for( int j=0; j<track->nIHHit(); j++ ){
      // 	HodoscopeLikeHit *hit = track->IHHit(cdsMan, j);
      // 	if( hit->seg()==IH_hit->seg() ) IH_flag = true;
      // }

      // if( IH_hit->seg()==11 ){
      // 	cout<<" vtxIH phi  : "<<vtxIH.Phi()<<"  posIH phi : "<<posIH.Phi()<<"  flag : "<<boolalpha<<IH_flag<<endl;
      // 	cout<<" vtxCDH phi : "<<vtxCDH.Phi()<<"  posCDH phi : "<<posCDH.Phi()<<"  flag : "<<boolalpha<<CDH_flag<<endl;
      // }

      if( CDH_flag && IH_flag ){
	eff_flag = true;
	if( track->Chi()<30. ){
	  eff_flag2 = true;
	}
      }
      if( eff_flag  ) h1 = (TH1F*)gFile-> Get("CDC_eff"), h1-> Fill(CDH_hit->seg());
      if( eff_flag2 ) h1 = (TH1F*)gFile-> Get("CDC_eff2"), h1-> Fill(CDH_hit->seg());

      if( eff_flag  ) h1 = (TH1F*)gFile-> Get("CDC_IH_eff"), h1-> Fill(IH_hit->seg());
      if( eff_flag2 ) h1 = (TH1F*)gFile-> Get("CDC_IH_eff2"), h1-> Fill(IH_hit->seg());
    }
  }

  histMan-> fill(header);

  Clear();
  return true;
}

void EventAnalysisHistwMC::Clear()
{
  header->Clear();
  histMan-> Clear();
  beamSpec-> Clear();
}

void EventAnalysisHistwMC::InitializeHistogram()
{
  rtFile-> cd();
  for( int i=0; i<20; i++ ){
    new TH1F(Form("trig_tdc_%d", i),     Form("Trigger TDC %d", i), 4095, 0.5, 4095.5);
    new TH1F(Form("trig_tdc_%d_eff", i), Form("Trigger TDC %d", i), 4095, 0.5, 4095.5);
  }

  new TH1F("TrigEff_CDH2",   "K#timesCDH2/f Trig Efficeincy",       6, -0.5, 5.5);
  new TH1F("TrigEff_N",   "   K#timesCDH1 Neutral Trig Efficeincy", 6, -0.5, 5.5);
  new TH1F("Trigger",   "Trigger",       20, -0.5, 19.5);
  new TH1F("TrigMode",  "Trigger Mode",  20, -0.5, 19.5);
  new TH1F("TrigMode1", "Trigger Mode",  20, -0.5, 19.5);
  new TH1F("TrigMode2", "Trigger Mode",  20, -0.5, 19.5);

  new TH2F("hitposBLC1", "BLC1 hit position", 500, -25, 25, 500, -25, 25);
  new TH2F("hitposBLC2", "BLC2 hit position", 500, -25, 25, 500, -25, 25);

  new TH1F("EventReduction", "Event Reduction", 21, -1.5, 19.5);
  new TH1F("Kf_Reduction", "K/f Event Reduction", 10, -0.5, 9.5);
  new TH1F("N_Reduction", "Neutral Event Reduction", 15, -0.5, 14.5);

  new TH1F("check_hist", "histogram for check", 100, -0.5, 99.5);

  new TH1F("nT0", "nT0", 6, -0.5, 5.5);

  new TH1F("hitpatBHD", "BHD hit pattern", 20, 0.5, 20.5);
  new TH1F("nBHD", "nBHD", 21, -0.5, 21.5);
  new TH1F("nBHD_woOut", "nBHD", 21, -0.5, 21.5);
  new TH1F("BHDT0",    "BHD-T0", 5000, 0.0, 50);
  new TH1F("BHDT0_Kf", "BHD-T0", 5000, 0.0, 50);
  new TH1F("BHDT0_K",  "BHD-T0", 5000, 0.0, 50);
  new TH1F("BHDT0_pi", "BHD-T0", 5000, 0.0, 50);

  new TH2F("hitposAC",    "AC hit pos", 500, -25, 25, 500, -25, 25);
  new TH2F("hitposAC_Kf", "AC hit pos", 500, -25, 25, 500, -25, 25);
  new TH2F("hitposAC_Kaon", "AC hit pos", 500, -25, 25, 500, -25, 25);
  new TH2F("hitposAC_Pion", "AC hit pos", 500, -25, 25, 500, -25, 25);

  new TH1F("ntrackBLC1", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_1", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_2", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("BLC1_time", "BLC1 track time", 10000, -200, 800);
  new TH1F("BLC1_chi2", "BLC1 chi-square", 1000, 0.0, 100);
  new TH1F("ntrackBLC1_pi", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_1_pi", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC1_2_pi", "BLC1 n track", 10, -0.5, 9.5);
  new TH1F("BLC1_time_pi", "BLC1 track time", 10000, -200, 800);
  new TH1F("BLC1_chi2_pi", "BLC1 chi-square", 1000, 0.0, 100);

  new TH1F("ntrackBLC2", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_1", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_2", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("BLC2_time", "BLC2 track time", 10000, -200, 800);
  new TH1F("BLC2_chi2", "BLC2 chi-square", 1000, 0.0, 100);
  new TH1F("ntrackBLC2_pi", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_1_pi", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("ntrackBLC2_2_pi", "BLC2 n track", 10, -0.5, 9.5);
  new TH1F("BLC2_time_pi", "BLC2 track time", 10000, -200, 800);
  new TH1F("BLC2_chi2_pi", "BLC2 chi-square", 1000, 0.0, 100);

  new TH1F("ntrackBPC", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_1", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_2", "BPC n track", 10, -0.5, 9.5);
  new TH1F("BPC_time", "BPC track time", 10000, -200, 800);
  new TH1F("BPC_chi2", "BPC chi-square", 1000, 0.0, 100);
  new TH1F("ntrackBPC_pi", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_1_pi", "BPC n track", 10, -0.5, 9.5);
  new TH1F("ntrackBPC_2_pi", "BPC n track", 10, -0.5, 9.5);
  new TH1F("BPC_time_pi", "BPC track time", 10000, -200, 800);
  new TH1F("BPC_chi2_pi", "BPC chi-square", 1000, 0.0, 100);

  new TH2F("BLC2BPC_diff", "BLC2 BPC pos diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff", "BLC2 BPC dir diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);
  new TH2F("BLC2BPC_diff_pi", "BLC2 BPC pos diff", 1000, -10, 10, 1000, -10, 10);
  new TH2F("BLC2BPC_dir_diff_pi", "BLC2 BPC dir diff", 1000, -0.1, 0.1, 1000, -0.1, 0.1);

  new TH1F("D5_mom",  "Beam Momentum by D5", 1000,  0.5, 1.5);
  for( int seg=1; seg<=20; seg++ )  new TH1F(Form("D5_mom_BHD%d", seg),  "Beam Momentum by D5", 1000,  0.5, 1.5);
  new TH1F("D5_chi2", "D5 chisquare", 1000,  0.0, 100);
  new TH1F("D5_mom_pi",  "Beam Momentum by D5", 1000,  0.5, 1.5);
  new TH1F("D5_chi2_pi", "D5 chisquare", 1000,  0.0, 100);

  new TH2F("BeamProf",    "Beam Profile at FF",        1000, -30, 30, 1000, -30, 30);
  new TH2F("BeamProf_pi",    "Beam Profile at FF",        1000, -30, 30, 1000, -30, 30);
  new TH2F("BeamProf_Kf", "Beam Profile at FF w/ K/f", 1000, -30, 30, 1000, -30, 30);
  new TH2F("BeamProf_CDH2", "Beam Profile at FF w/ CDH2", 1000, -30, 30, 1000, -30, 30);

  //  initHistBLDC();
}

void EventAnalysisHistwMC::Finalize()
{
  std::cout<<"EventAnalysisHistwMC Finish   Event Number : "<<Event_Number<<"  Event Number on Spill : "<<Event_Number_onSpill<<std::endl;
  delete cdcFile;
  delete blMan;
  delete cdsMan;
  delete header;

  delete bltrackMan;
  delete beamSpec;

  histMan->finit();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisHistwMC *event = new EventAnalysisHistwMC();
  return (EventTemp*)event;
}
