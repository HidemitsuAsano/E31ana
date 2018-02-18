//#define LOG

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "DataContainer.h"
#include "AnalysisMan.h"

#include <TNtuple.h>

static const int CDSAnaMode = 6; // T0 1hit & BPC track & Forward hit

class EventAnalysisFull: public EventTemp
{
public:
  EventAnalysisFull();
  ~EventAnalysisFull();
private:
  std::ofstream ofs_log;
  time_t time0, now_time;

  TFile          *rtFile;

  ScalerMan      *scaMan;
  EventHeader    *header;

  CDSHitMan      *cdsMan;
  BeamLineHitMan *blMan;

  CDSTrackingMan   *cdsTrackMan;
  BeamLineTrackMan *blTrackMan;

  AnalysisMan    *anaMan;

public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysisFull::EventAnalysisFull()
  : EventTemp()
{
}

EventAnalysisFull::~EventAnalysisFull()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisFull::Initialize( ConfMan *conf )
{
  std::cout << " Enter EventAnalysisFull::Initialize " << std::endl;
  time(&time0);

  confMan = conf;  
#ifdef LOG
  struct tm *local = localtime(&time0);
  std::string tmp_str = confMan->GetOutFileName();
  tmp_str.erase(0,5);
  tmp_str.erase(tmp_str.size()-5,tmp_str.size());
  std::string log_name = "log/"+tmp_str+".log";

  ofs_log.open(log_name.c_str());
  ofs_log<<"########### Initialize ################################"<<std::endl;
  ofs_log<<"### Data : "<<local->tm_year+1900<<"/"<<local->tm_mon+1<<"/"<<local->tm_mday<<"  "<<local->tm_hour<<":"<<local->tm_min<<std::endl;
  ofs_log<<"### ConfFileName : "<<confMan->GetConfFileName()<<std::endl;
  ofs_log<<"### OutFileName  : "<<confMan->GetOutFileName()<<std::endl;
  ofs_log<<"#######################################################"<<std::endl;
  ofs_log<<"> Start !!!"<<std::endl;
#endif

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );

  scaMan      = new ScalerMan();
  header      = new EventHeader();
  cdsMan      = new CDSHitMan();
  blMan       = new BeamLineHitMan();
  cdsTrackMan = new CDSTrackingMan();
  blTrackMan  = new BeamLineTrackMan();

  anaMan = new AnalysisMan();
  anaMan-> SetDumpLevel(999);
  anaMan-> SetConfMan(confMan);
  anaMan-> SetScaMan(scaMan);
  anaMan-> SetHeader(header);
  anaMan-> SetCDSHitMan(cdsMan);
  anaMan-> SetBeamLineHitMan(blMan);
  anaMan-> SetCDSTrackingMan(cdsTrackMan);
  anaMan-> SetBeamLineTrackMan(blTrackMan);
  anaMan-> SetCDSAnaMode(6);

  anaMan-> Initialize();

  InitializeHistogram();
}

void EventAnalysisFull::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisFull::UAna( TKOHitCollection *tko )
{
  rtFile-> cd();

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%1000==0 ){
    time(&now_time);
    double diff_time = difftime(now_time, time0);
    int mm = (int)floor(diff_time/60);
    int ss = (int)(diff_time-mm*60);

    std::cout<<" Event# : "<<std::setw(6)<<Event_Number<<std::setw(16)
             <<"BlockEvent# : "<<std::setw(3)<<Block_Event_Number<<std::setw(9)
             <<"time : "<<mm<<" m "<<ss<<" s"<<std::endl;
#ifdef LOG
    ofs_log<<" Event# : "<<std::setw(6)<<Event_Number<<std::setw(16)
	   <<"BlockEvent# : "<<std::setw(3)<<Block_Event_Number<<std::setw(9)
	   <<"time : "<<mm<<" m "<<ss<<" s"<<std::endl;
#endif
  }

  TH1F *h1;
  TH2F *h2;
  TNtuple *tup;

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  anaMan-> Execute();

  // for user defined evnt 
  h1 = (TH1F*)gFile-> Get("UserEvent"), h1-> Fill(0.5);
  if( anaMan->nT0()==1 )                h1-> Fill(1.5);
  if( anaMan->nBHD()==1 )               h1-> Fill(2.5);
  if( anaMan->BeamPID()==Beam_Kaon )    h1-> Fill(3.5);
  if( anaMan->D5() )                    h1-> Fill(4.5);
  if( anaMan->BPC_ID()>=0 )             h1-> Fill(5.5);
  if( cdsTrackMan->nGoodTrack()>0 )     h1-> Fill(6.5);
  if( anaMan->GoodNC() )                h1-> Fill(7.5);
  if( anaMan->GoodPC() )                h1-> Fill(8.5);
  if( anaMan->GoodCVC() )               h1-> Fill(9.5);
  // for event reduction
  h1 = (TH1F*)gFile-> Get("EventReduction"), h1-> Fill(0.5);
  if( anaMan->nT0()==1 ){  h1-> Fill(1.5);
    if( anaMan->nBHD()==1 ){  h1-> Fill(2.5);
      if( anaMan->BeamPID()==Beam_Kaon ){  h1-> Fill(3.5);
	if( anaMan->D5() ){  h1-> Fill(4.5);
	  if( anaMan->BPC_ID()>=0 ){  h1-> Fill(5.5);
	    if( cdsTrackMan->nGoodTrack()>0 ){  h1-> Fill(6.5);
	      if( anaMan->GoodNC() ) h1-> Fill(7.5);
	      if( anaMan->GoodPC() ) h1-> Fill(8.5);
	      if( anaMan->GoodCVC() ) h1-> Fill(9.5);
	    }
	  }
	}
      }
    }
  }

  if( anaMan->D5() ) h1 = (TH1F*)gFile-> Get("D5mom"), h1-> Fill(anaMan->D5mom());
  if( anaMan->BHDT0() ) h1 = (TH1F*)gFile-> Get("BHDT0"), h1-> Fill(anaMan->BHDT0_TOF());

  if( anaMan->IsVertexPos() ){
    TVector3 VtxPos = anaMan-> GetVertexPos();
    tup = (TNtuple*)gFile-> Get("VtxPos"), tup-> Fill(VtxPos.x(),VtxPos.y(),VtxPos.z());
    if( !anaMan->VtxCut() ) tup = (TNtuple*)gFile-> Get("VtxPos_wcut"), tup-> Fill(VtxPos.x(),VtxPos.y(),VtxPos.z());
    else                    tup = (TNtuple*)gFile-> Get("VtxPos_out_cut"), tup-> Fill(VtxPos.x(),VtxPos.y(),VtxPos.z());
  }

  // for CDS 
  for( int i=0; i<anaMan->nGoodTrack(); i++ ){
    double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
    double beta = anaMan-> CDS_Beta(i);
    double mass2 = anaMan-> CDS_Mass2(i);

    h2 = (TH2F*)gFile-> Get("CDS_overbeta_mom"), h2-> Fill(1./beta, mom);
    h2 = (TH2F*)gFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
  }

  for( int i=0; i<anaMan->nCDSPiPi(); i++ ){
    h1 = (TH1F*)gFile-> Get("CDS_pipi"), h1-> Fill(anaMan->CDSPiPiIM(i)->GetIM());
  }

  for( int i=0; i<anaMan->nCDSPiP(); i++ ){
    h1 = (TH1F*)gFile-> Get("CDS_pip"), h1-> Fill(anaMan->CDSPiPIM(i)->GetIM());
  }

  for( int i=0; i<anaMan->nCDSKP(); i++ ){
    h1 = (TH1F*)gFile-> Get("CDS_kp"), h1-> Fill(anaMan->CDSKPIM(i)->GetIM());
  }
  // for CDS w Vertex Cut
  if( !anaMan->VtxCut() ){
    for( int i=0; i<anaMan->nGoodTrack(); i++ ){
      double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
      double beta = anaMan-> CDS_Beta(i);
      double mass2 = anaMan-> CDS_Mass2(i);

      h2 = (TH2F*)gFile-> Get("CDS_overbeta_mom"), h2-> Fill(1./beta, mom);
      h2 = (TH2F*)gFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
    }

    for( int i=0; i<anaMan->nCDSPiPi(); i++ ){
      h1 = (TH1F*)gFile-> Get("CDS_pipi"), h1-> Fill(anaMan->CDSPiPiIM(i)->GetIM());
    }

    for( int i=0; i<anaMan->nCDSPiP(); i++ ){
      h1 = (TH1F*)gFile-> Get("CDS_pip"), h1-> Fill(anaMan->CDSPiPIM(i)->GetIM());
    }

    for( int i=0; i<anaMan->nCDSKP(); i++ ){
      h1 = (TH1F*)gFile-> Get("CDS_kp"), h1-> Fill(anaMan->CDSKPIM(i)->GetIM());
    }
  }

  // for Forward Counter slect Kaon beam
  if( anaMan->BeamPID()==Beam_Kaon ){
    if( anaMan->GoodNC() ){
      h1 = (TH1F*)gFile-> Get("VtxNC1st"), h1-> Fill(anaMan->VtxNC1st_TOF());
      h1 = (TH1F*)gFile-> Get("NC1st_overbeta"), h1-> Fill(1./anaMan->NC1st_Beta());
      if( !anaMan->VtxCut() ){
	h1 = (TH1F*)gFile-> Get("VtxNC1st_wcut"), h1-> Fill(anaMan->VtxNC1st_TOF());
	h1 = (TH1F*)gFile-> Get("NC1st_overbeta_wcut"), h1-> Fill(1./anaMan->NC1st_Beta());
      }
      for( int i=0; i<anaMan->nNC(); i++ ){
	h1 = (TH1F*)gFile-> Get("VtxNC"), h1-> Fill(anaMan->VtxNC_TOF(i));
	h1 = (TH1F*)gFile-> Get("NC_overbeta"), h1-> Fill(1./anaMan->NC_Beta(i));
      }

      if( anaMan->Npip() ) h1 = (TH1F*)gFile-> Get("Npip_im"), h1-> Fill(anaMan->NPip_IM()->GetIM());
      if( anaMan->Npim() ) h1 = (TH1F*)gFile-> Get("Npim_im"), h1-> Fill(anaMan->NPim_IM()->GetIM());
      if( anaMan->Nk() )   h1 = (TH1F*)gFile-> Get("Nk_im"),   h1-> Fill(anaMan->NK_IM()->GetIM());
      if( !anaMan->VtxCut() ){
	if( anaMan->Npip() ) h1 = (TH1F*)gFile-> Get("Npip_im_wcut"), h1-> Fill(anaMan->NPip_IM()->GetIM());
	if( anaMan->Npim() ) h1 = (TH1F*)gFile-> Get("Npim_im_wcut"), h1-> Fill(anaMan->NPim_IM()->GetIM());
	if( anaMan->Nk() )   h1 = (TH1F*)gFile-> Get("Nk_im_wcut"),   h1-> Fill(anaMan->NK_IM()->GetIM());
      }
	
      if( anaMan->D5() ){
	if( anaMan->Neutron() ){
	  double mm = anaMan-> BeamN_MM()-> GetMM();
	  h1 = (TH1F*)gFile-> Get("KN_MM"), h1-> Fill(mm);
	  if( !anaMan->VtxCut() ) h1 = (TH1F*)gFile-> Get("KN_MM_wcut"), h1-> Fill(mm);
	}
      }
    }

    if( anaMan->ForwardCharged() ){
      double FCdE = anaMan-> GetFC()->emean();
      if( anaMan->GoodPC() ){
	h1 = (TH1F*)gFile-> Get("VtxPC"), h1-> Fill(anaMan->VtxFC_TOF());
	h2 = (TH2F*)gFile-> Get("PC_tof_mom"), h2-> Fill(anaMan->VtxFC_TOF(), anaMan->FCMom());
	h2 = (TH2F*)gFile-> Get("PC_tof_dE"),  h1-> Fill(anaMan->VtxFC_TOF(), FCdE);
	if( !anaMan->VtxCut() && !anaMan->FCADCCut() ){
	  h2 = (TH2F*)gFile-> Get("PC_tof_mom_wcut"), h2-> Fill(anaMan->VtxFC_TOF(), anaMan->FCMom());
	  h2 = (TH2F*)gFile-> Get("PC_tof_dE_wcut"),  h1-> Fill(anaMan->VtxFC_TOF(), FCdE);
	}
      }
      if( anaMan->GoodCVC() ){
	h1 = (TH1F*)gFile-> Get("VtxCVC"), h1-> Fill(anaMan->VtxFC_TOF());
	h2 = (TH2F*)gFile-> Get("CVC_tof_mom"), h2-> Fill(anaMan->VtxFC_TOF(), anaMan->FCMom());
	h2 = (TH2F*)gFile-> Get("CVC_tof_dE"),  h1-> Fill(anaMan->VtxFC_TOF(), FCdE);
	if( !anaMan->VtxCut() && !anaMan->FCADCCut() ) {
	  h2 = (TH2F*)gFile-> Get("CVC_tof_mom_wcut"), h2-> Fill(anaMan->VtxFC_TOF(), anaMan->FCMom());
	  h2 = (TH2F*)gFile-> Get("CVC_tof_dE_wcut"),  h1-> Fill(anaMan->VtxFC_TOF(), FCdE);
	}
      }

      if( anaMan->FC_PID()==CDS_Proton ){
	h2 = (TH2F*)gFile-> Get("P_mom_mom2"), h2-> Fill( anaMan->FCMom(), anaMan->FCMom2() );

	if( anaMan->Ppip() ) h1 = (TH1F*)gFile-> Get("Ppip_im"), h1-> Fill(anaMan->PPip_IM()->GetIM());
	if( anaMan->Ppim() ) h1 = (TH1F*)gFile-> Get("Ppim_im"), h1-> Fill(anaMan->PPim_IM()->GetIM());
	if( anaMan->Pk() )   h1 = (TH1F*)gFile-> Get("Pk_im"),   h1-> Fill(anaMan->PK_IM()->GetIM());

	if( !anaMan->VtxCut() && !anaMan->FCADCCut() ){
	  h2 = (TH2F*)gFile-> Get("P_mom_mom2_wcut"), h2-> Fill( anaMan->FCMom(), anaMan->FCMom2() );
	  if( anaMan->Ppip() ) h1 = (TH1F*)gFile-> Get("Ppip_im_wcut"), h1-> Fill(anaMan->PPip_IM()->GetIM());
	  if( anaMan->Ppim() ) h1 = (TH1F*)gFile-> Get("Ppim_im_wcut"), h1-> Fill(anaMan->PPim_IM()->GetIM());
	  if( anaMan->Pk() )   h1 = (TH1F*)gFile-> Get("Pk_im_wcut"),   h1-> Fill(anaMan->PK_IM()->GetIM());
	}
      }

      if( anaMan->FC_PID()==CDS_Deuteron ){
	h2 = (TH2F*)gFile-> Get("D_mom_mom2"), h2-> Fill( anaMan->FCMom(), anaMan->FCMom2() );

	if( anaMan->Dpip() ) h1 = (TH1F*)gFile-> Get("Dpip_im"), h1-> Fill(anaMan->DPip_IM()->GetIM());
	if( anaMan->Dpim() ) h1 = (TH1F*)gFile-> Get("Dpim_im"), h1-> Fill(anaMan->DPim_IM()->GetIM());
	if( anaMan->Dk() )   h1 = (TH1F*)gFile-> Get("Dk_im"),   h1-> Fill(anaMan->DK_IM()->GetIM());

	if( !anaMan->VtxCut() && !anaMan->FCADCCut() ){
	  h2 = (TH2F*)gFile-> Get("D_mom_mom2_wcut"), h2-> Fill( anaMan->FCMom(), anaMan->FCMom2() );

	  if( anaMan->Dpip() ) h1 = (TH1F*)gFile-> Get("Dpip_im_wcut"), h1-> Fill(anaMan->DPip_IM()->GetIM());
	  if( anaMan->Dpim() ) h1 = (TH1F*)gFile-> Get("Dpim_im_wcut"), h1-> Fill(anaMan->DPim_IM()->GetIM());
	  if( anaMan->Dk() )   h1 = (TH1F*)gFile-> Get("Dk_im_wcut"),   h1-> Fill(anaMan->DK_IM()->GetIM());
	}
      }

      if( anaMan->D5() ){
	if( anaMan->FC_PID()==CDS_Proton ){
	  double mm = anaMan-> BeamP_MM()-> GetMM();
	  h1 = (TH1F*)gFile-> Get("KP_MM"), h1-> Fill(mm);
	  if( !anaMan->VtxCut() && !anaMan->FCADCCut() ) h1 = (TH1F*)gFile-> Get("KP_MM_wcut"), h1-> Fill(mm);
	}
	if( anaMan->FC_PID()==CDS_Deuteron ){
	  double mm = anaMan-> BeamD_MM()-> GetMM();
	  h1 = (TH1F*)gFile-> Get("KD_MM"), h1-> Fill(mm);
	  if( !anaMan->VtxCut() && !anaMan->FCADCCut() ) h1 = (TH1F*)gFile-> Get("KD_MM_wcut"), h1-> Fill(mm);
	}
      }
    }
  }

  header-> Clear();
  blMan->  Clear();
  cdsMan-> Clear();
  cdsTrackMan-> Clear();
  blTrackMan->  Clear();

  anaMan-> Clear();

  return true;
}

void EventAnalysisFull::Finalize()
{
  std::cout << " Enter EventAnalysisFull::Finalize " << std::endl;
  rtFile-> cd();
  confMan->SaveCode();
  confMan->SaveParams();
  gFile->Write();
  gFile->Close();

#ifdef LOG
  time(&now_time);
  double diff_time = difftime(now_time, time0);
  int mm = (int)floor(diff_time/60);
  int ss = (int)(diff_time-mm*60);
  ofs_log<<"> Finish !!!"<<std::endl;
  ofs_log<<" Event# : "<<std::setw(6)<<Event_Number<<std::setw(16)
	 <<"BlockEvent# : "<<std::setw(3)<<Block_Event_Number<<std::setw(9)
	 <<"time : "<<mm<<" m "<<ss<<" s"<<std::endl;
#endif

  delete anaMan;

  delete blMan;
  delete cdsMan;
  delete header;
  delete cdsTrackMan;
  delete blTrackMan;
}

void EventAnalysisFull::InitializeHistogram()
{
  rtFile-> cd();
  std::cout<<"===== EventAnalysisFull::InitializeHistogram start ====="<<std::endl;
  new TH1F("UserEvent", "User Defined Event", 10, 0, 10);
  new TH1F("EventReduction", "User Event Reduction", 10, 0, 10);

  new TH1F("BHDT0", "BHD-To TOF", 50, 0, 50);
  new TH1F("D5mom", "beam momentum by D5", 1200, 0, 1.2);

  new TNtuple("VtxPos", "Vertex Position", "x:y:z");
  new TH2F("CDS_overbeta_mom", "CDS 1/#beta vs mom", 500, 0, 10, 300, -1.5, 1.5);
  new TH2F("CDS_mass2_mom", "CDS mass^{2} vs mom", 600, -1, 5, 300, -1.5, 1.5);
  new TH1F("CDS_pipi", "CDS #pi^{-} & #pi^{+} Invariant Mass", 1000, 0, 1.0);
  new TH1F("CDS_pip", "CDS #pi^{-} & proton Invariant Mass", 500, 1.0, 1.5);
  new TH1F("CDS_kp", "CDS K^{-} & proton Invariant Mass", 1000, 1.0, 2.0);

  new TH1F("VtxNC", "Vtx-NC TOF", 2000, 0, 500);
  new TH1F("VtxNC1st", "Vtx-NC1st TOF", 2000, 0, 500);
  new TH1F("NC1st_overbeta", "NC 1st hit 1/#beta", 1000, 0, 20);
  new TH1F("NC_overbeta", "NC 1/#beta", 1000, 0, 20);

  new TH1F("VtxPC", "Vtx-PC TOF", 2000, 0, 500);
  new TH2F("PC_tof_mom", "PC TOF vs momentum", 2000, 0, 500, 1500, 0, 1.5);
  new TH2F("PC_tof_dE", "PC TOF vs dE", 2000, 0, 500, 500, 0, 50);

  new TH1F("VtxCVC", "Vtx-CVC TOF", 2000, 0, 500);
  new TH2F("CVC_tof_mom", "CVC TOF vs momentum", 2000, 0, 500, 1500, 0, 1.5);
  new TH2F("CVC_tof_dE", "CVC TOF vs dE", 2000, 0, 500, 500, 0, 50);

  new TH2F("P_mom_mom2", "P momentum by angle vs tof", 150, 0, 1.5, 150, 0, 1.5);
  new TH2F("D_mom_mom2", "D momentum by angle vs tof", 150, 0, 1.5, 150, 0, 1.5);

  new TH1F("KN_MM", "^{3}He(K^{-}, n) Missing Mass", 500, 1.0, 3.0);
  new TH1F("KP_MM", "^{3}He(K^{-}, p) Missing Mass", 500, 1.0, 3.0);
  new TH1F("KD_MM", "^{3}He(K^{-}, d) Missing Mass", 500, 0.0, 2.0);

  new TH1F("Npip_im", "Neutron & #pi^{+} Invariant Mass", 500, 1.0 ,3.0);
  new TH1F("Npim_im", "Neutron & #pi^{-} Invariant Mass", 500, 1.0 ,3.0);
  new TH1F("Nk_im", "Neutron & K^{-} Invariant Mass", 500, 1.0 ,3.0);

  new TH1F("Ppip_im", "Forward P & #pi^{+} Invariant Mass", 2000, 1.0 ,3.0);
  new TH1F("Ppim_im", "Forward P & #pi^{-} Invariant Mass", 2000, 1.0 ,3.0);
  new TH1F("Pk_im", "Forward P & K^{-} Invariant Mass", 2000, 1.0 ,3.0);

  new TH1F("Dpip_im", "Forward D & #pi^{+} Invariant Mass", 2000, 1.0 ,3.0);
  new TH1F("Dpim_im", "Fprward D & #pi^{-} Invariant Mass", 2000, 1.0 ,3.0);
  new TH1F("Dk_im", "Forward D & K^{-} Invariant Mass", 2000, 1.0 ,3.0);

  std::cout<<"  === for with cut ==="<<std::endl;
  new TH2F("CDS_overbeta_mom_wcut", "CDS 1/#beta vs mom w cut", 500, 0, 10, 300, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_wcut", "CDS mass^{2} vs mom w cut", 600, -1, 5, 300, -1.5, 1.5);
  new TH1F("CDS_pipi_wcut", "CDS #pi^{-} & #pi^{+} Invariant Mass w cut", 1000, 0, 1.0);
  new TH1F("CDS_pip_wcut", "CDS #pi^{-} & proton Invariant Mass w cut", 500, 1.0, 1.5);
  new TH1F("CDS_kp_wcut", "CDS K^{-} & proton Invariant Mass w cut", 1000, 1.0, 2.0);

  new TNtuple("VtxPos_out_cut", "Vertex Position", "x:y:z");
  new TNtuple("VtxPos_wcut", "Vertex Position", "x:y:z");

  new TH1F("VtxNC1st_wcut", "Vtx-NC1st TOF w cut", 2000, 0, 500);
  new TH1F("NC1st_overbeta_wcut", "Vtx-NC1st 1/#beta w cut", 2000, 0, 500);

  new TH2F("PC_tof_mom_wcut", "PC TOF vs momentum w cut", 2000, 0, 500, 1500, 0, 1.5);
  new TH2F("PC_tof_dE_wcut", "PC TOF vs dE w cut", 2000, 0, 500, 500, 0, 50);

  new TH2F("CVC_tof_mom_wcut", "CVC TOF vs momentum w cut", 2000, 0, 500, 1500, 0, 1.5);
  new TH2F("CVC_tof_dE_wcut", "CVC TOF vs dE w cut", 2000, 0, 500, 500, 0, 50);

  new TH2F("P_mom_mom2_wcut", "P momentum by angle vs tof w cut", 150, 0, 1.5, 150, 0, 1.5);
  new TH2F("D_mom_mom2_wcut", "D momentum by angle vs tof w cut", 150, 0, 1.5, 150, 0, 1.5);

  new TH1F("KN_MM_wcut", "^{3}He(K^{-}, n) Missing Mass w cut", 500, 1.0, 3.0);
  new TH1F("KP_MM_wcut", "^{3}He(K^{-}, p) Missing Mass w cut", 500, 1.0, 3.0);
  new TH1F("KD_MM_wcut", "^{3}He(K^{-}, d) Missing Mass w cut", 500, 0.0, 2.0);

  new TH1F("Npip_im_wcut", "Neutron & #pi^{+} Invariant Mass w cut", 2000, 1.0 ,3.0);
  new TH1F("Npim_im_wcut", "Neutron & #pi^{-} Invariant Mass w cut", 2000, 1.0 ,3.0);
  new TH1F("Nk_im_wcut", "Neutron & K^{-} Invariant Mass w cut", 2000, 1.0 ,3.0);

  new TH1F("Ppip_im_wcut", "Forward P & #pi^{+} Invariant Mass w cut", 2000, 1.0 ,3.0);
  new TH1F("Ppim_im_wcut", "Forward P & #pi^{-} Invariant Mass w cut", 2000, 1.0 ,3.0);
  new TH1F("Pk_im_wcut", "Forward P & K^{-} Invariant Mass w cut", 2000, 1.0 ,3.0);

  new TH1F("Dpip_im_wcut", "Forward D & #pi^{+} Invariant Mass w cut", 2000, 1.0 ,3.0);
  new TH1F("Dpim_im_wcut", "Fprward D & #pi^{-} Invariant Mass w cut", 2000, 1.0 ,3.0);
  new TH1F("Dk_im_wcut", "Forward D & K^{-} Invariant Mass w cut", 2000, 1.0 ,3.0);

  std::cout<<"===== EventAnalysisFull::InitializeHistogram end ====="<<std::endl;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisFull *event = new EventAnalysisFull();
  return (EventTemp*)event;
}
