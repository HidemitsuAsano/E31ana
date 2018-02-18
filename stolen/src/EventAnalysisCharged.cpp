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

static const double PI = 6*asin(0.5);

class EventAnalysisCheck: public EventTemp
{
public:
  EventAnalysisCheck();
  ~EventAnalysisCheck();
private:
  std::ofstream   ofsLog;
  int NumOfGoodPC;
  int NumOfGoodPC_K;
  int NumOfPC_K_P;
  int NumOfPC_Pk;

  int NumOfPC_Ppim;
  int nFC_lambda1520;
  int nFC_lambda;
  int nFC_lambda_p;

  int NumOfPC_Ppip;

  int NumOfPC_K_D;
  int NumOfPC_D_pipi;
  int NumOfPC_D_pip;
  int NumOfPC_D_k;

  int NumOfGoodPC_Pi;
  int NumOfPC_Pi_P;
  int NumOfPC_Pi_D;

  int NumOfGoodCVC;
  int NumOfGoodCVC_K;
  int NumOfCVC_K_P;
  int NumOfCVC_K_D;

  int NumOfGoodCVC_Pi;
  int NumOfCVC_Pi_P;
  int NumOfCVC_Pi_D;


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

EventAnalysisCheck::EventAnalysisCheck()
  : EventTemp()
{
}

EventAnalysisCheck::~EventAnalysisCheck()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisCheck::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisCheck::Initialize " << std::endl;
#endif
  confMan = conf;

  std::string tmp_str = confMan->GetOutFileName();
  tmp_str.erase(0,5);
  tmp_str.erase(tmp_str.size()-5,tmp_str.size());
  std::string log_name = "log/"+tmp_str+".log";
  std::cout<<" Log File Name : "<<log_name<<std::endl;
  ofsLog.open(log_name.c_str());
  ofsLog<<"########### Initialize ################################"<<std::endl;
  ofsLog<<"### ConfFileName : "<<confMan->GetConfFileName()<<std::endl;
  ofsLog<<"### OutFileName  : "<<confMan->GetOutFileName()<<std::endl;
  ofsLog<<"#######################################################"<<std::endl;
  ofsLog<<"> Start !!!"<<std::endl;
  NumOfGoodPC     = 0;
  NumOfPC_Pk      = 0;
  NumOfPC_Ppim    = 0;
  NumOfPC_Ppip    = 0;

  nFC_lambda     = 0;
  nFC_lambda1520 = 0;
  nFC_lambda_p   = 0;

  NumOfGoodPC_K   = 0;
  NumOfGoodPC_Pi  = 0;
  NumOfPC_K_P = 0;
  NumOfPC_K_D = 0;
  NumOfPC_Pi_P = 0;
  NumOfPC_Pi_D = 0;

  NumOfPC_D_pipi = 0;
  NumOfPC_D_pip = 0;
  NumOfPC_D_k = 0;

  NumOfGoodCVC    = 0;
  NumOfGoodCVC_K  = 0;
  NumOfGoodCVC_Pi = 0;
  NumOfCVC_K_P = 0;
  NumOfCVC_K_D = 0;
  NumOfCVC_Pi_P = 0;
  NumOfCVC_Pi_D = 0;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );

  scaMan      = new ScalerMan();
  header      = new EventHeader();
  cdsMan      = new CDSHitMan();
  blMan       = new BeamLineHitMan();
  cdsTrackMan = new CDSTrackingMan();
  blTrackMan  = new BeamLineTrackMan();

  anaMan = new AnalysisMan();
  anaMan-> SetConfMan(confMan);
  anaMan-> SetScaMan(scaMan);
  anaMan-> SetHeader(header);
  anaMan-> SetCDSMan(cdsMan);
  anaMan-> SetBLMan(blMan);
  anaMan-> SetCDSTrackMan(cdsTrackMan);
  anaMan-> SetCDSTrackMode(99);
  anaMan-> SetBLTrackMan(blTrackMan);

  anaMan-> Initialize();

  InitializeHistogram();
}

void EventAnalysisCheck::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisCheck::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca<< std::endl;
  for( int i=0; i<nsca; i++ ){
   std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif
  scaMan->Clear();
}

bool EventAnalysisCheck::UAna( TKOHitCollection *tko )
{
  rtFile-> cd();
#if 0
  std::cout << " Enter EventAnalysisCheck::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;

  TH1F *h1;
  TH2F *h2;
    
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  header->SetBlockEventNumber(Block_Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  for( int i=0; i<20; i++ ){
    int val = header-> pattern(i);
    if( val>0 ){
      h1 = (TH1F*)gFile-> Get("Pattern"), h1-> Fill(i);
      h1 = (TH1F*)gFile-> Get(Form("Pattern_%d",i)), h1-> Fill(val);
    }
  }

  anaMan-> Ana();

  TH1F *h1_ev;
  h1_ev = (TH1F*)gFile-> Get("Event"), h1_ev-> Fill(0);
  TH1F *h1_er;
  h1_er = (TH1F*)gFile-> Get("EventReduction"), h1_er-> Fill(0);
  bool reduction_flag = true;
  if( anaMan-> nT0()==1 ){
    h1_ev-> Fill(1);
    h1_er-> Fill(1);
  }
  else reduction_flag = false;

  if( anaMan-> nBHD()==1 ){
    h1_ev-> Fill(2);
    if( reduction_flag ) h1_er-> Fill(2);
  }
  else reduction_flag = false;

  if( anaMan-> BeamPID()==Beam_Kaon ){
    h1_ev-> Fill(3);
    if( reduction_flag ) h1_er-> Fill(3);
  }
  else reduction_flag = false;

  if( anaMan-> D5() ){
    h1_ev-> Fill(4);
    if( reduction_flag ) h1_er-> Fill(4);
  }
  else reduction_flag = false;

  if( anaMan-> VtxBeam() ){
    h1_ev-> Fill(5);
    if( reduction_flag ) h1_er-> Fill(5);
  }
  else reduction_flag = false;

  if( anaMan-> ForwardCharged() ){
    h1_ev-> Fill(6);
    if( reduction_flag ) h1_er-> Fill(6);
  }
  else reduction_flag = false;

  if( anaMan-> GoodPC() ){
    h1_ev-> Fill(7);
    if( reduction_flag ) h1_er-> Fill(7);
  }
  else reduction_flag = false;

  if( anaMan-> GoodCVC() ){
    h1_ev-> Fill(8);
    if( reduction_flag ) h1_er-> Fill(8);
  }
  else reduction_flag = false;

  if( anaMan-> BHDT0() ){
    h1 = (TH1F*)gFile-> Get("BHDT0"), h1-> Fill(anaMan->BHDT0_TOF());
  }
  if( anaMan-> D5() ){
    h1 = (TH1F*)gFile-> Get("D5mom"), h1-> Fill(anaMan->D5mom());
  }

  if( anaMan-> GoodPC() || anaMan-> GoodCVC() ){
    h1 = (TH1F*)gFile-> Get("CDS_nGoodTrack"), h1-> Fill(anaMan->nGoodTrack());
    for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
      double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
      double beta = anaMan->VtxCDH_Beta(i);
      double mass2 = mom*mom*(1/(beta*beta)-1);
      double mass = sqrt(mass2);
      int cds_pid = anaMan-> CDSGoodTrack(i)-> PID();
      h1 = (TH1F*)gFile-> Get("CDS_PID"), h1->Fill(cds_pid);

      h2 = (TH2F*)gFile-> Get("CDS_overbeta_mom"), h2-> Fill(1./beta, mom);
      h2 = (TH2F*)gFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
      h2 = (TH2F*)gFile-> Get("CDS_mass_mom"), h2-> Fill(mass, mom);
    }
    if( anaMan->VtxBeam() ){
      TVector3 vtx_pos = anaMan->GetVertexBeam()->GetVertexPos_mean();
      h2 = (TH2F*)gFile-> Get("VtxPos_xy"), h2-> Fill(vtx_pos.x(), vtx_pos.y());
      h2 = (TH2F*)gFile-> Get("VtxPos_zx"), h2-> Fill(vtx_pos.z(), vtx_pos.x());
      h2 = (TH2F*)gFile-> Get("VtxPos_zy"), h2-> Fill(vtx_pos.z(), vtx_pos.y());
      if( anaMan->VtxCut() ){
	h2 = (TH2F*)gFile-> Get("VtxPos_xy_wCut"), h2-> Fill(vtx_pos.x(), vtx_pos.y());
	h2 = (TH2F*)gFile-> Get("VtxPos_zx_wCut"), h2-> Fill(vtx_pos.z(), vtx_pos.x());
	h2 = (TH2F*)gFile-> Get("VtxPos_zy_wCut"), h2-> Fill(vtx_pos.z(), vtx_pos.y());
      }
    }

    for( int i=0; i<anaMan->nT0NC(); i++ ){
      double T0NC_tof = anaMan-> T0NC_TOF(i);
      h1 = (TH1F*)gFile-> Get("T0NC"), h1-> Fill(T0NC_tof);
    }
    if( anaMan->T0CVC() ){
      double T0CVC_tof = anaMan-> T0CVC_TOF();
      h1 = (TH1F*)gFile-> Get("T0CVC"), h1-> Fill(T0CVC_tof);
    }
    if( anaMan->T0PC() ){
      double T0PC_tof = anaMan-> T0PC_TOF();
      h1 = (TH1F*)gFile-> Get("T0PC"), h1-> Fill(T0PC_tof);
    }

    for( int i=0; i<anaMan->nVtxNC(); i++ ){
      double VtxNC_tof = anaMan-> VtxNC_TOF(i);
      h1 = (TH1F*)gFile-> Get("VtxNC"), h1-> Fill(VtxNC_tof);
    }
    if( anaMan->VtxCVC() ){
      double VtxCVC_tof = anaMan-> VtxCVC_TOF();
      h1 = (TH1F*)gFile-> Get("VtxCVC"), h1-> Fill(VtxCVC_tof);
    }
    if( anaMan->VtxPC() ){
      double VtxPC_tof = anaMan-> VtxPC_TOF();
      h1 = (TH1F*)gFile-> Get("VtxPC"), h1-> Fill(VtxPC_tof);
    }

    if( anaMan-> ForwardCharged() ){
      double fc_mm = anaMan-> FC_MM();

      if( anaMan-> GoodPC() ){
	NumOfGoodPC++;
	int PCseg = anaMan->PC(0)->seg();
	double PCdE = anaMan-> PC(0)-> emean();

	if( anaMan->FCADC() ){
	  h2 = (TH2F*)gFile-> Get("PC_tof_angwADC"), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get("PC_tof_momwADC"), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get("PC_tof_flwADC"),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get("PC_mom_offsetwADC"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());
	  h2 = (TH2F*)gFile-> Get("PC_tof_dEwADC"),  h2-> Fill(anaMan->VtxPC_TOF(),PCdE);
	}
	if( anaMan->FCWVC() ){
	  h2 = (TH2F*)gFile-> Get("PC_tof_angwWVC"), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get("PC_tof_momwWVC"), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get("PC_tof_flwWVC"),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get("PC_mom_offsetwWVC"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());
	  h2 = (TH2F*)gFile-> Get("PC_tof_dEwWVC"),  h2-> Fill(anaMan->VtxPC_TOF(),PCdE);
	}

	h2 = (TH2F*)gFile-> Get("PC_tof_dE"),  h2-> Fill(anaMan->VtxPC_TOF(),PCdE);

	if( !anaMan->FCADC() ){
	  h2 = (TH2F*)gFile-> Get("PC_tof_ang"), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get("PC_tof_mom"), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get("PC_tof_fl"),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get("PC_mom_offset"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());
	  h2 = (TH2F*)gFile-> Get("PC_mass_mom"), h2-> Fill(anaMan->FCCalcMass(), anaMan->FC_Mom());

	  h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_ang",PCseg)), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_mom",PCseg)), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_fl",PCseg)),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get(Form("PC%d_mom_offset",PCseg)),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF());

	  if( anaMan-> FC_PID()==CDS_Proton ){
	    h2 = (TH2F*)gFile-> Get("PC_p_mom_mom2"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	    h2 = (TH2F*)gFile-> Get("PC_p_res_mom"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());
	  }
	  if( anaMan-> FC_PID()==CDS_Deuteron ){
	    h2 = (TH2F*)gFile-> Get("PC_d_mom_mom2"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	    h2 = (TH2F*)gFile-> Get("PC_d_res_mom"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());
	  }

	  if( anaMan->BeamPID()==Beam_Kaon ){
	    NumOfGoodPC_K++;
	    h2 = (TH2F*)gFile-> Get("PC_tof_ang_k"), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get("PC_tof_mom_k"), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get("PC_tof_fl_k"),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get("PC_mom_offset_k"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_ang_k",PCseg)), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_mom_k",PCseg)), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_fl_k",PCseg)),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get(Form("PC%d_mom_offset_k",PCseg)),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get("PC_mass_mom_k"), h2-> Fill(anaMan->FCCalcMass(), anaMan->FC_Mom());
	    if( anaMan-> FC_PID()==CDS_Proton ){
	      TVector3 vtx_pos = anaMan->GetVertexBeam()->GetVertexPos_mean();
	      h2 = (TH2F*)gFile-> Get("VtxPos_xy_PC_k_p"), h2-> Fill(vtx_pos.x(), vtx_pos.y());
	      h2 = (TH2F*)gFile-> Get("VtxPos_zx_PC_k_p"), h2-> Fill(vtx_pos.z(), vtx_pos.x());
	      h2 = (TH2F*)gFile-> Get("VtxPos_zy_PC_k_p"), h2-> Fill(vtx_pos.z(), vtx_pos.y());

	      NumOfPC_K_P++;
	      h2 = (TH2F*)gFile-> Get("PC_p_mom_mom2_k"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("PC_p_res_mom_k"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->FPk() ){
		double im = anaMan->FPk_IM();
	      //	      std::cout<<" \\(^O^)/ Forward Proton & CDS K- Invariant Mass : "<<im<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("IM_pk"), h1-> Fill(im);
		if( anaMan-> FCLambda1520() ){
		  nFC_lambda1520++;
		  std::cout<<" \\(^o^)/ Forward Proton & CDS K- Lambda(1520) Event "<<im<<"[GeV] !!!"<<std::endl;
		  std::cout<<"  > CDS low mass  hit : "<<anaMan->nCDSlow()<<std::endl;
		  std::cout<<"  > CDS pi+       hit : "<<anaMan->nCDSpip()<<std::endl;
		  std::cout<<"  > CDS K+        hit : "<<anaMan->nCDSkp()<<std::endl;
		  std::cout<<"  > CDS Proton    hit : "<<anaMan->nCDSp()<<std::endl;
		  std::cout<<"  > CDS Deuteron  hit : "<<anaMan->nCDSd()<<std::endl;
		  std::cout<<"  > CDS pi-       hit : "<<anaMan->nCDSpim()<<std::endl;
		  std::cout<<"  > CDS K-        hit : "<<anaMan->nCDSkm()<<std::endl;
		  std::cout<<"  > CDS high mass hit : "<<anaMan->nCDShigh()<<std::endl;

		  ofsLog<<" EventNumber : "<<Event_Number<<std::endl;
		  ofsLog<<" \\(^o^)/ Forward Proton & CDS K- Lambda(1520) Event "<<im<<"[GeV] !!!"<<std::endl;
		  ofsLog<<"  > CDS low mass  hit : "<<anaMan->nCDSlow()<<std::endl;
		  ofsLog<<"  > CDS pi+       hit : "<<anaMan->nCDSpip()<<std::endl;
		  ofsLog<<"  > CDS K+        hit : "<<anaMan->nCDSkp()<<std::endl;
		  ofsLog<<"  > CDS Proton    hit : "<<anaMan->nCDSp()<<std::endl;
		  ofsLog<<"  > CDS Deuteron  hit : "<<anaMan->nCDSd()<<std::endl;
		  ofsLog<<"  > CDS pi-       hit : "<<anaMan->nCDSpim()<<std::endl;
		  ofsLog<<"  > CDS K-        hit : "<<anaMan->nCDSkm()<<std::endl;
		  ofsLog<<"  > CDS high mass hit : "<<anaMan->nCDShigh()<<std::endl;
		  if( anaMan-> D5() ){
		    h1 = (TH1F*)gFile-> Get("PC_k_p_mm_lambda1520"), h1-> Fill(fc_mm);
		  }
		}
		NumOfPC_Pk++;
	      }

	      if( anaMan->FPpim() ){
		double im = anaMan->FPpim_IM();
	      //	      std::cout<<" \\(^O^)/ Forward Proton & CDS pi- Invariant Mass : "<<im<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("IM_ppim"), h1-> Fill(im);
		if( anaMan-> FCLambda() ){
		  nFC_lambda++;
		  std::cout<<" \\(^o^)/ Forward Proton & CDS pi- Lambda Event "<<im<<"[GeV] !!!"<<std::endl;
		  std::cout<<"  > CDS low mass  hit : "<<anaMan->nCDSlow()<<std::endl;
		  std::cout<<"  > CDS pi+       hit : "<<anaMan->nCDSpip()<<std::endl;
		  std::cout<<"  > CDS K+        hit : "<<anaMan->nCDSkp()<<std::endl;
		  std::cout<<"  > CDS Proton    hit : "<<anaMan->nCDSp()<<std::endl;
		  std::cout<<"  > CDS Deuteron  hit : "<<anaMan->nCDSd()<<std::endl;
		  std::cout<<"  > CDS pi-       hit : "<<anaMan->nCDSpim()<<std::endl;
		  std::cout<<"  > CDS K-        hit : "<<anaMan->nCDSkm()<<std::endl;
		  std::cout<<"  > CDS high mass hit : "<<anaMan->nCDShigh()<<std::endl;

		  ofsLog<<" EventNumber : "<<Event_Number<<std::endl;
		  ofsLog<<" \\(^o^)/ Forward Proton & CDS pi- Lambda Event "<<im<<"[GeV] !!!"<<std::endl;
		  ofsLog<<"  > CDS low mass  hit : "<<anaMan->nCDSlow()<<std::endl;
		  ofsLog<<"  > CDS pi+       hit : "<<anaMan->nCDSpip()<<std::endl;
		  ofsLog<<"  > CDS K+        hit : "<<anaMan->nCDSkp()<<std::endl;
		  ofsLog<<"  > CDS Proton    hit : "<<anaMan->nCDSp()<<std::endl;
		  ofsLog<<"  > CDS Deuteron  hit : "<<anaMan->nCDSd()<<std::endl;
		  ofsLog<<"  > CDS pi-       hit : "<<anaMan->nCDSpim()<<std::endl;
		  ofsLog<<"  > CDS K-        hit : "<<anaMan->nCDSkm()<<std::endl;
		  ofsLog<<"  > CDS high mass hit : "<<anaMan->nCDShigh()<<std::endl;

		  if( anaMan-> D5() ){
		    h1 = (TH1F*)gFile-> Get("PC_k_p_mm_lambda"), h1-> Fill(fc_mm);
		  }
		  if( anaMan->nCDSp()==1 ){
		    double im2 = anaMan-> FCLambdaP_IM();
		    std::cout<<" \\(^O^)/ Forward Proton & CDS pi- Lambda & CDS p Event "<<std::endl;
		    std::cout<<"          Lambda & CDS p Invariant Mass : "<<im2<<"[GeV]"<<std::endl;
		    ofsLog<<" \\(^O^)/ Forward Proton & CDS pi- Lambda & CDS p Event "<<std::endl;
		    ofsLog<<"          Lambda & CDS p Invariant Mass : "<<im2<<"[GeV]"<<std::endl;

		    h1 = (TH1F*)gFile-> Get("IM_lambda_p"), h1-> Fill(im2);
		    nFC_lambda_p++;
		  }
		}
		NumOfPC_Ppim++;
	      }

	      if( anaMan->FPpip() ){
		double im = anaMan->FPpip_IM();
	      //	      std::cout<<" \\(^O^)/ Forward Proton & CDS pi+ Invariant Mass : "<<im<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("IM_ppip"), h1-> Fill(im);
		NumOfPC_Ppip++;
	      }

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ PC hit K(3He, p) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("PC_k_p_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_PC_k_p"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_PC_k_p"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_PC_k_p"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	    if( anaMan-> FC_PID()==CDS_Deuteron ){
	      TVector3 vtx_pos = anaMan->GetVertexBeam()->GetVertexPos_mean();
	      h2 = (TH2F*)gFile-> Get("VtxPos_xy_PC_k_d"), h2-> Fill(vtx_pos.x(), vtx_pos.y());
	      h2 = (TH2F*)gFile-> Get("VtxPos_zx_PC_k_d"), h2-> Fill(vtx_pos.z(), vtx_pos.x());
	      h2 = (TH2F*)gFile-> Get("VtxPos_zy_PC_k_d"), h2-> Fill(vtx_pos.z(), vtx_pos.y());

	      NumOfPC_K_D++;
	      h2 = (TH2F*)gFile-> Get("PC_d_mom_mom2_k"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("PC_d_res_mom_k"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
		std::cout<<" \\(^O^)/ PC hit K(3He, d) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		std::cout<<"  > CDS low mass  hit : "<<anaMan->nCDSlow()<<std::endl;
		std::cout<<"  > CDS pi+       hit : "<<anaMan->nCDSpip()<<std::endl;
		std::cout<<"  > CDS K+        hit : "<<anaMan->nCDSkp()<<std::endl;
		std::cout<<"  > CDS Proton    hit : "<<anaMan->nCDSp()<<std::endl;
		std::cout<<"  > CDS Deuteron  hit : "<<anaMan->nCDSd()<<std::endl;
		std::cout<<"  > CDS pi-       hit : "<<anaMan->nCDSpim()<<std::endl;
		std::cout<<"  > CDS K-        hit : "<<anaMan->nCDSkm()<<std::endl;
		std::cout<<"  > CDS high mass hit : "<<anaMan->nCDShigh()<<std::endl;

		ofsLog<<" EventNumber : "<<Event_Number<<std::endl;
		ofsLog<<" \\(^O^)/ PC hit K(3He, d) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		ofsLog<<"  > CDS low mass  hit : "<<anaMan->nCDSlow()<<std::endl;
		ofsLog<<"  > CDS pi+       hit : "<<anaMan->nCDSpip()<<std::endl;
		ofsLog<<"  > CDS K+        hit : "<<anaMan->nCDSkp()<<std::endl;
		ofsLog<<"  > CDS Proton    hit : "<<anaMan->nCDSp()<<std::endl;
		ofsLog<<"  > CDS Deuteron  hit : "<<anaMan->nCDSd()<<std::endl;
		ofsLog<<"  > CDS pi-       hit : "<<anaMan->nCDSpim()<<std::endl;
		ofsLog<<"  > CDS K-        hit : "<<anaMan->nCDSkm()<<std::endl;
		ofsLog<<"  > CDS high mass hit : "<<anaMan->nCDShigh()<<std::endl;

		h1 = (TH1F*)gFile-> Get("PC_k_d_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_PC_k_d"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_PC_k_d"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_PC_k_d"), h2-> Fill(1./beta, mom);
		}
		if( anaMan->nCDSpim()==1 && anaMan->nCDSpip()==1 ){
		  std::cout<<" \\(^o^)/ K(3He, d)X pi+ pi- tag event : "<<Event_Number<<std::endl;
		  std::cout<<"   Missing mass = "<<fc_mm<<"[GeV]"<<std::endl;

		  ofsLog<<" \\(^o^)/ K(3He, d)X pi+ pi- tag event : "<<Event_Number<<std::endl;
		  ofsLog<<"   Missing mass = "<<fc_mm<<"[GeV]"<<std::endl;
		  h1 = (TH1F*)gFile-> Get("PC_k_d_mm_pipi"), h1-> Fill(fc_mm);
		  NumOfPC_D_pipi++;
		}
		if( anaMan->nCDSpim()==1 && anaMan->nCDSp()==1 ){
		  std::cout<<" \\(^o^)/ K(3He, d)X p pi- tag event : "<<Event_Number<<std::endl;
		  std::cout<<"   Missing mass = "<<fc_mm<<"[GeV]"<<std::endl;

		  ofsLog<<" \\(^o^)/ K(3He, d)X p pi- tag event : "<<Event_Number<<std::endl;
		  ofsLog<<"   Missing mass = "<<fc_mm<<"[GeV]"<<std::endl;
		  h1 = (TH1F*)gFile-> Get("PC_k_d_mm_ppi"), h1-> Fill(fc_mm);
		  NumOfPC_D_pip++;
		}
		if( anaMan->nCDSkm()==1 ){
		  std::cout<<" \\(^o^)/ K(3He, d)X k- tag event : "<<Event_Number<<std::endl;
		  std::cout<<"   Missing mass = "<<fc_mm<<"[GeV]"<<std::endl;

		  ofsLog<<" \\(^o^)/ K(3He, d)X k- tag event : "<<Event_Number<<std::endl;
		  ofsLog<<"   Missing mass = "<<fc_mm<<"[GeV]"<<std::endl;
		  h1 = (TH1F*)gFile-> Get("PC_k_d_mm_k"), h1-> Fill(fc_mm);
		  NumOfPC_D_k++;
		}
	      }
	      if( anaMan->FDpim() ){
		double im = anaMan->FDpim_IM();
		std::cout<<" \\(^O^)/ Forward Deuteron & CDS pi- Invariant Mass : "<<im<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("IM_dpim"), h1-> Fill(im);
	      }
	      if( anaMan->FDpip() ){
		double im = anaMan->FDpip_IM();
		std::cout<<" \\(^O^)/ Forward Deuteron & CDS pi+ Invariant Mass : "<<im<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("IM_dpip"), h1-> Fill(im);
	      }
	      if( anaMan->FDk() ){
		double im = anaMan->FDk_IM();
		std::cout<<" \\(^O^)/ Forward Deuteron & CDS K- Invariant Mass : "<<im<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("IM_dk"), h1-> Fill(im);
	      }
	    }
	  }

	  else if( anaMan->BeamPID()==Beam_Pion ){
	    NumOfGoodPC_Pi++;
	    h2 = (TH2F*)gFile-> Get("PC_tof_ang_pi"), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get("PC_tof_mom_pi"), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get("PC_tof_fl_pi"),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get("PC_mom_offset_pi"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_ang_pi",PCseg)), h2-> Fill(anaMan->VtxPC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_mom_pi",PCseg)), h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get(Form("PC%d_tof_fl_pi",PCseg)),  h2-> Fill(anaMan->VtxPC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get(Form("PC%d_mom_offset_pi",PCseg)),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxPC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get("PC_mass_mom_pi"), h2-> Fill(anaMan->FCCalcMass(), anaMan->FC_Mom());
	    if( anaMan-> FC_PID()==CDS_Proton ){
	      NumOfPC_Pi_P++;
	      h2 = (TH2F*)gFile-> Get("PC_p_mom_mom2_pi"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("PC_p_res_mom_pi"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ PC hit pi(3He, p) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("PC_pi_p_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_PC_pi_p"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_PC_pi_p"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_PC_pi_p"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	    if( anaMan-> FC_PID()==CDS_Deuteron ){
	      NumOfPC_Pi_D++;
	      h2 = (TH2F*)gFile-> Get("PC_d_mom_mom2_pi"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("PC_d_res_mom_pi"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ PC hit pi(3He, d) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("PC_pi_d_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_PC_pi_d"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_PC_pi_d"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_PC_pi_d"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	  }
	}
      }
      if( anaMan-> GoodCVC() ){
	int CVCseg = anaMan->CVC(0)->seg();
	double CVCdE = anaMan->CVC(0)->seg();
	if( anaMan->FCADC() ){
	  h2 = (TH2F*)gFile-> Get("CVC_tof_angwADC"), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get("CVC_tof_momwADC"), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get("CVC_tof_flwADC"),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get("CVC_mom_offsetwADC"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());
	  h2 = (TH2F*)gFile-> Get("CVC_tof_dEwADC"),  h2-> Fill(anaMan->VtxCVC_TOF(), CVCdE);
	}
	if( anaMan->FCWVC() ){
	  h2 = (TH2F*)gFile-> Get("CVC_tof_angwWVC"), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get("CVC_tof_momwWVC"), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get("CVC_tof_flwWVC"),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get("CVC_mom_offsetwWVC"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());
	  h2 = (TH2F*)gFile-> Get("CVC_tof_dEwWVC"),  h2-> Fill(anaMan->VtxCVC_TOF(), CVCdE);
	}

	h2 = (TH2F*)gFile-> Get("CVC_tof_dE"),  h2-> Fill(anaMan->VtxCVC_TOF(), CVCdE);

	if( !anaMan->FCADC() ){
	  h2 = (TH2F*)gFile-> Get("CVC_tof_ang"), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get("CVC_tof_mom"), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get("CVC_tof_fl"),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get("CVC_mom_offset"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());
	  h2 = (TH2F*)gFile-> Get("CVC_mass_mom"), h2-> Fill(anaMan->FCCalcMass(), anaMan->FC_Mom());

	  h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_ang",CVCseg)), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	  h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_mom",CVCseg)), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	  h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_fl",CVCseg)),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	  h2 = (TH2F*)gFile-> Get(Form("CVC%d_mom_offset",CVCseg)),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());

	  if( anaMan-> FC_PID()==CDS_Proton ){
	    h2 = (TH2F*)gFile-> Get("CVC_p_mom_mom2"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	    h2 = (TH2F*)gFile-> Get("CVC_p_res_mom"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());
	  }
	  if( anaMan-> FC_PID()==CDS_Deuteron ){
	    h2 = (TH2F*)gFile-> Get("CVC_d_mom_mom2"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	    h2 = (TH2F*)gFile-> Get("CVC_d_res_mom"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());
	  }

	  if( anaMan->BeamPID()==Beam_Kaon ){
	    NumOfGoodCVC_K++;
	    h2 = (TH2F*)gFile-> Get("CVC_tof_ang_k"), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get("CVC_tof_mom_k"), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get("CVC_tof_fl_k"),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get("CVC_mom_offset_k"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_ang_k",CVCseg)), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_mom_k",CVCseg)), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_fl_k",CVCseg)),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_mom_offset_k",CVCseg)),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get("CVC_mass_mom_k"), h2-> Fill(anaMan->FCCalcMass(), anaMan->FC_Mom());
	    if( anaMan-> FC_PID()==CDS_Proton ){
	      NumOfCVC_K_P++;
	      h2 = (TH2F*)gFile-> Get("CVC_p_mom_mom2_k"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("CVC_p_res_mom_k"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ CVC hit K(3He, p) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("CVC_k_p_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_CVC_k_p"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_CVC_k_p"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_CVC_k_p"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	    if( anaMan-> FC_PID()==CDS_Deuteron ){
	      NumOfCVC_K_D++;
	      h2 = (TH2F*)gFile-> Get("CVC_d_mom_mom2_k"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("CVC_d_res_mom_k"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ CVC hit K(3He, d) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("CVC_k_d_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_CVC_k_d"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_CVC_k_d"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_PC_k_d"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	  }
	  else if( anaMan->BeamPID()==Beam_Pion ){
	    NumOfGoodCVC_Pi++;
	    h2 = (TH2F*)gFile-> Get("CVC_tof_ang_pi"), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get("CVC_tof_mom_pi"), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get("CVC_tof_fl_pi"),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get("CVC_mom_offset_pi"),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_ang_pi",CVCseg)), h2-> Fill(anaMan->VtxCVC_TOF(), 180.*anaMan->BendingAngle()/PI);
	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_mom_pi",CVCseg)), h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_Mom());
	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_tof_fl_pi",CVCseg)),  h2-> Fill(anaMan->VtxCVC_TOF(), anaMan->FC_FL());
	    h2 = (TH2F*)gFile-> Get(Form("CVC%d_mom_offset_pi",CVCseg)),  h2-> Fill(anaMan->FC_Mom(), anaMan->VtxCVC_TOF()-anaMan->FCCalcTOF());

	    h2 = (TH2F*)gFile-> Get("CVC_mass_mom_pi"), h2-> Fill(anaMan->FCCalcMass(), anaMan->FC_Mom());
	    if( anaMan-> FC_PID()==CDS_Proton ){
	      NumOfCVC_Pi_P++;
	      h2 = (TH2F*)gFile-> Get("CVC_p_mom_mom2_pi"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("CVC_p_res_mom_pi"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ CVC hit pi(3He, p) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("CVC_pi_p_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_CVC_pi_p"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_CVC_pi_p"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_CVC_pi_p"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	    if( anaMan-> FC_PID()==CDS_Deuteron ){
	      NumOfCVC_Pi_D++;
	      h2 = (TH2F*)gFile-> Get("CVC_d_mom_mom2_pi"), h2-> Fill(anaMan->FC_Mom(), anaMan->FC_Mom2());
	      h2 = (TH2F*)gFile-> Get("CVC_d_res_mom_pi"), h2-> Fill(anaMan->FC_Mom2()-anaMan->FC_Mom(), anaMan->FC_Mom());

	      if( anaMan->D5() ){
	      //	      std::cout<<" \(^O^)/ CVC hit pi(3He, d) Missing Mass "<<fc_mm<<"[GeV]"<<std::endl;
		h1 = (TH1F*)gFile-> Get("CVC_pi_d_mm"), h1-> Fill(fc_mm);

		h1 = (TH1F*)gFile-> Get("nCDS_CVC_pi_d"), h1-> Fill(anaMan->nGoodTrack());
		for( int i=0; i<anaMan->nGoodTrack2(); i++ ){
		  double mom = anaMan-> CDSGoodTrack(i)-> Momentum();
		  double beta = anaMan->VtxCDH_Beta(i);
		  double mass2 = mom*mom*(1/(beta*beta)-1);
		  double mass = sqrt(mass2);

		  int cds_pid = anaMan->CDSGoodTrack(i)-> PID();
		  h1 = (TH1F*)gFile-> Get("CDSpid_CVC_pi_d"), h1->Fill(cds_pid);
		  h2 = (TH2F*)gFile-> Get("CDS_CVC_pi_d"), h2-> Fill(1./beta, mom);
		}
	      }
	    }
	  }
	  NumOfGoodCVC++;
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

void EventAnalysisCheck::Finalize()
{
  rtFile-> cd();
  std::cout << " Enter EventAnalysisCheck::Finalize " << std::endl;
  ofsLog<<"##### Finalize #############################"<<std::endl;
  ofsLog<<" All Event : "<<Event_Number<<std::endl;
  ofsLog<<"### PC hit Event ###########################"<<std::endl;
  ofsLog<<" Good PC hit Event  : "<<NumOfGoodPC<<"   "<<100.*NumOfGoodPC/Event_Number<<"[%]"<<std::endl;
  ofsLog<<"   Good PC hit Event w k beam  : "<<NumOfGoodPC_K<<std::endl;
  ofsLog<<"       K-(3He, p) reaction : "<<NumOfPC_K_P<<std::endl;
  ofsLog<<"          Forward p + CDS K-  : "<<NumOfPC_Pk<<std::endl;
  ofsLog<<"           > Forward p + CDS K- Lambda(1520) : "<<nFC_lambda1520<<std::endl;
  ofsLog<<"          Forward p + CDS pi- : "<<NumOfPC_Ppim<<std::endl;
  ofsLog<<"           > Forward p + CDS pi- Lambda      : "<<nFC_lambda<<std::endl;
  ofsLog<<"           -> Forward p + CDS pi- Lambda & CDS p : "<<nFC_lambda_p<<std::endl;
  ofsLog<<"          Forward p + CDS pi+ : "<<NumOfPC_Ppim<<std::endl;
  ofsLog<<"       K-(3He, d) reaction : "<<NumOfPC_K_D<<std::endl;
  ofsLog<<"   Num Of Good PC Event by pi beam : "<<NumOfGoodPC_Pi<<std::endl;
  ofsLog<<"       pi-(3He, p) reaction : "<<NumOfPC_Pi_P<<std::endl;
  ofsLog<<"       pi-(3He, d) reaction : "<<NumOfPC_Pi_D<<std::endl;
  ofsLog<<"### CVC hit Event ##########################"<<std::endl;
  ofsLog<<" Num Of Good CVC Event : "<<NumOfGoodCVC<<"   "<<100.*NumOfGoodCVC/Event_Number<<"[%]"<<std::endl;
  ofsLog<<"   Num Of Good CVC Event by k beam  : "<<NumOfGoodCVC_K<<std::endl;
  ofsLog<<"       K-(3He, p) reaction : "<<NumOfCVC_K_P<<std::endl;
  ofsLog<<"       K-(3He, d) reaction : "<<NumOfCVC_K_D<<std::endl;
  ofsLog<<"   Num Of Good CVC Event by pi beam : "<<NumOfGoodCVC_Pi<<std::endl;
  ofsLog<<"       pi-(3He, p) reaction : "<<NumOfCVC_Pi_P<<std::endl;
  ofsLog<<"       pi-(3He, d) reaction : "<<NumOfCVC_Pi_D<<std::endl;
  ofsLog<<"### foward Charged w CDS tagging event #####"<<std::endl;
  ofsLog<<" K-(3He, d) pi+ pi- tag event : "<<NumOfPC_D_pipi<<std::endl;
  ofsLog<<" K-(3He, d) p   pi- tag event : "<<NumOfPC_D_pip<<std::endl;
  ofsLog<<" K-(3He, d) k-      tag event : "<<NumOfPC_D_k<<std::endl;
  ofsLog<<"##### Finish ###############################"<<std::endl;

  gFile->Write();
  gFile->Close();

  delete anaMan;

  delete blMan;
  delete cdsMan;
  delete header;
  delete cdsTrackMan;
  delete blTrackMan;
}

void EventAnalysisCheck::InitializeHistogram()
{
  rtFile-> cd();
  std::cout<<" EventAnalysisCheck::InitializeHistogram start ... "<<std::endl;
  std::cout<<"  for Check etc. "<<std::endl;
  new TH1F("Event", "Event", 10, 0, 10);
  new TH1F("EventReduction", "Event Reduction", 10, 0, 10);

  new TH1F("Pattern", "Pattern", 21, 0, 21);
  for( int i=0; i<20; i++ ){
    new TH1F(Form("Pattern_%d",i), Form("Pattern_%d",i), 4095, 0, 4095);
  }

  new TH1F("BHDT0", "BHD-T0 TOF", 1000, 0, 50);
  new TH1F("D5mom", "beam momentum by D5", 1200, 0, 1.2);

  std::cout<<"  for CDS Check"<<std::endl;
  new TH1F("CDS_nGoodTrack", "CDS Num Of Good Track", 10, 0, 10);
  new TH1F("CDS_PID", "CDS PID", 7, 0, 7);

  new TH2F("CDS_overbeta_mom", "CDS 1/#beta vs mom", 500, 0, 10, 500, -1.2, 1.2);
  new TH2F("CDS_mass_mom", "CDS mass vs mom", 500, -1, 10, 500, -1.2, 1.2);
  new TH2F("CDS_mass2_mom", "CDS mass^{2} vs mom", 500, -1, 10, 500, -1.2, 1.2);

  new TH2F("VtxPos_xy", "VertexPos XY", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zx", "VertexPos ZX", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zy", "VertexPos ZY", 500, -20, 20, 500, -20, 20);

  new TH2F("VtxPos_xy_wCut", "VertexPos XY w Cut", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zx_wCut", "VertexPos ZX w Cut", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zy_wCut", "VertexPos ZY w Cut", 500, -20, 20, 500, -20, 20);

  new TH2F("VtxPos_xy_PC_k_p", "VertexPos XY", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zx_PC_k_p", "VertexPos ZX", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zy_PC_k_p", "VertexPos ZY", 500, -20, 20, 500, -20, 20);

  new TH2F("VtxPos_xy_PC_k_d", "VertexPos XY", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zx_PC_k_d", "VertexPos ZX", 500, -20, 20, 500, -20, 20);
  new TH2F("VtxPos_zy_PC_k_d", "VertexPos ZY", 500, -20, 20, 500, -20, 20);

  new TH1F("VtxNC",  "Vtx-NC TOF",  2000, -10, 190);
  new TH1F("VtxCVC", "Vtx-CVC TOF", 3000, -50, 250);
  new TH1F("VtxPC",  "Vtx-PC TOF",  3000, -50, 250);

  new TH1F("T0NC",  "Vtx-NC TOF",  2000, -10, 190);
  new TH1F("T0CVC", "Vtx-CVC TOF", 3000, -50, 250);
  new TH1F("T0PC",  "Vtx-PC TOF",  3000, -50, 250);

  std::cout<<"  for Forward Charged"<<std::endl;
  new TH1F("FC_PID", "Forward Charged PID", 5, 0, 5);

  new TH2F("CVC_mass_mom",    "CVC_mass vs mom",          500, 0, 5, 200, 0, 2.0);
  new TH2F("CVC_mass_mom_k",  "CVC_mass vs mom  k beam",  500, 0, 5, 200, 0, 2.0);
  new TH2F("CVC_mass_mom_pi", "CVC_mass vs mom  pi beam", 500, 0, 5, 200, 0, 2.0);

  new TH2F("CVC_p_mom_mom2",    "CVC P hit mom by ang vs mom by TOF",          200, 0, 2.0, 200, 0, 2.0);
  new TH2F("CVC_p_mom_mom2_k",  "CVC P hit mom by ang vs mom by TOF  k beam",  200, 0, 2.0, 200, 0, 2.0);
  new TH2F("CVC_p_mom_mom2_pi", "CVC P hit mom by ang vs mom by TOF  pi beam", 200, 0, 2.0, 200, 0, 2.0);

  new TH2F("CVC_d_mom_mom2",    "CVC D hit mom by ang vs mom by TOF",          200, 0, 2.0, 200, 0, 2.0);
  new TH2F("CVC_d_mom_mom2_k",  "CVC D hit mom by ang vs mom by TOF  k beam",  200, 0, 2.0, 200, 0, 2.0);
  new TH2F("CVC_d_mom_mom2_pi", "CVC D hit mom by ang vs mom by TOF  pi beam", 200, 0, 2.0, 200, 0, 2.0);

  new TH2F("CVC_p_res_mom",    "CVC P hit res vs mom ",         200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("CVC_p_res_mom_k",  "CVC P hit res vs mom  k beam",  200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("CVC_p_res_mom_pi", "CVC P hit res vs mom  pi beam", 200, -0.2, 0.2, 200, 0, 2.0);

  new TH2F("CVC_d_res_mom",    "CVC D hit res vs mom ",         200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("CVC_d_res_mom_k",  "CVC D hit res vs mom  k beam",  200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("CVC_d_res_mom_pi", "CVC D hit res vs mom  pi beam", 200, -0.2, 0.2, 200, 0, 2.0);

  new TH2F("CVC_tof_ang", "T0-CVC TOF vs bending angle", 400, -10, 190, 350, -5, 35);
  new TH2F("CVC_tof_mom", "T0-CVC TOF vs momentum", 400, -10, 190, 500, 0, 2.0);
  new TH2F("CVC_tof_fl",  "T0-CVC TOF vs Fligth length", 400, -10, 190, 400, 1200, 1700);
  new TH2F("CVC_mom_offset",  "CVC mom vs Calc Offset", 500, 0, 2.0, 1000, -100, 100);
  new TH2F("CVC_tof_dE", "T0-CVC tof vs dE", 400, -10, 190, 500, 0, 50);

  new TH2F("CVC_tof_ang_k", "T0-CVC TOF vs bending angle  k beam", 400, -10, 190, 350, -5, 35);
  new TH2F("CVC_tof_mom_k", "T0-CVC TOF vs momentum  k beam", 400, -10, 190, 500, 0, 2.0);
  new TH2F("CVC_tof_fl_k",  "T0-CVC TOF vs Fligth length  k beam", 400, -10, 190, 400, 1200, 1700);
  new TH2F("CVC_mom_offset_k",  "CVC mom vs Calc Offset  k beam", 500, 0, 2.0, 1000, -100, 100);

  new TH2F("CVC_tof_ang_pi", "T0-CVC TOF vs bending angle  #pi beam", 400, -10, 190, 350, -5, 35);
  new TH2F("CVC_tof_mom_pi", "T0-CVC TOF vs momentum  #pi beam", 400, -10, 190, 500, 0, 2.0);
  new TH2F("CVC_tof_fl_pi",  "T0-CVC TOF vs Fligth length  #pi beam", 400, -10, 190, 400, 1200, 1700);
  new TH2F("CVC_mom_offset_pi",  "CVC mom vs Calc Offset  #pi beam", 500, 0, 2.0, 1000, -100, 100);

  for( int seg=1; seg<=NumOfCVCSegments; seg++ ){
    new TH2F(Form("CVC%d_tof_ang",seg), Form("T0-CVC seg%d TOF vs bending angle",seg), 400, -10, 190, 350, -5, 35);
    new TH2F(Form("CVC%d_tof_mom",seg), Form("T0-CVC seg%d TOF vs momentum",seg), 400, -10, 190, 500, 0, 2.0);
    new TH2F(Form("CVC%d_tof_fl",seg),  Form("T0-CVC seg%d TOF vs Fligth length",seg), 400, -10, 190, 400, 1200, 1700);
    new TH2F(Form("CVC%d_mom_offset",seg),  Form("CVC seg%d mom vs Calc Offset",seg), 500, 0, 2.0, 1000, -100, 100);

    new TH2F(Form("CVC%d_tof_ang_k",seg), Form("T0-CVC seg%d TOF vs bending angle  k beam",seg), 400, -10, 190, 350, -5, 35);
    new TH2F(Form("CVC%d_tof_mom_k",seg), Form("T0-CVC seg%d TOF vs momentum  k beam",seg), 400, -10, 190, 500, 0, 2.0);
    new TH2F(Form("CVC%d_tof_fl_k",seg),  Form("T0-CVC seg%d TOF vs Fligth length  k beam",seg), 400, -10, 190, 400, 1200, 1700);
    new TH2F(Form("CVC%d_mom_offset_k",seg),  Form("CVC seg%d mom vs Calc Offset  k beam",seg), 500, 0, 2.0, 1000, -100, 100);

    new TH2F(Form("CVC%d_tof_ang_pi",seg), Form("T0-CVC seg%d TOF vs bending angle  #pi beam",seg), 400, -10, 190, 350, -5, 35);
    new TH2F(Form("CVC%d_tof_mom_pi",seg), Form("T0-CVC seg%d TOF vs momentum  #pi beam",seg), 400, -10, 190, 500, 0, 2.0);
    new TH2F(Form("CVC%d_tof_fl_pi",seg),  Form("T0-CVC seg%d TOF vs Fligth length  #pi beam",seg), 400, -10, 190, 400, 1200, 1700);
    new TH2F(Form("CVC%d_mom_offset_pi",seg),  Form("CVC seg%d mom vs Calc Offset  #pi beam",seg), 500, 0, 2.0, 1000, -100, 100);
  }

  new TH2F("PC_mass_mom",    "PC_mass vs mom",          500, 0, 5, 200, 0, 2.0);
  new TH2F("PC_mass_mom_k",  "PC_mass vs mom  k beam",  500, 0, 5, 200, 0, 2.0);
  new TH2F("PC_mass_mom_pi", "PC_mass vs mom  pi beam", 500, 0, 5, 200, 0, 2.0);

  new TH2F("PC_p_mom_mom2",    "PC P hit mom by ang vs mom by TOF",          200, 0, 2.0, 200, 0, 2.0);
  new TH2F("PC_p_mom_mom2_k",  "PC P hit mom by ang vs mom by TOF  k beam",  200, 0, 2.0, 200, 0, 2.0);
  new TH2F("PC_p_mom_mom2_pi", "PC P hit mom by ang vs mom by TOF  pi beam", 200, 0, 2.0, 200, 0, 2.0);

  new TH2F("PC_d_mom_mom2",    "PC D hit mom by ang vs mom by TOF",          200, 0, 2.0, 200, 0, 2.0);
  new TH2F("PC_d_mom_mom2_k",  "PC D hit mom by ang vs mom by TOF  k beam",  200, 0, 2.0, 200, 0, 2.0);
  new TH2F("PC_d_mom_mom2_pi", "PC D hit mom by ang vs mom by TOF  pi beam", 200, 0, 2.0, 200, 0, 2.0);

  new TH2F("PC_p_res_mom",    "PC P hit res vs mom ",         200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("PC_p_res_mom_k",  "PC P hit res vs mom  k beam",  200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("PC_p_res_mom_pi", "PC P hit res vs mom  pi beam", 200, -0.2, 0.2, 200, 0, 2.0);

  new TH2F("PC_d_res_mom",    "PC D hit res vs mom ",         200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("PC_d_res_mom_k",  "PC D hit res vs mom  k beam",  200, -0.2, 0.2, 200, 0, 2.0);
  new TH2F("PC_d_res_mom_pi", "PC D hit res vs mom  pi beam", 200, -0.2, 0.2, 200, 0, 2.0);

  new TH2F("PC_tof_ang", "T0-PC TOF vs bending angle", 400, -10, 190, 350, -5, 35);
  new TH2F("PC_tof_mom", "T0-PC TOF vs momentum", 400, -10, 190, 500, 0, 2.0);
  new TH2F("PC_tof_fl",  "T0-PC TOF vs Fligth length", 400, -10, 190, 400, 1200, 1700);
  new TH2F("PC_mom_offset",  "PC mom vs Calc Offset", 500, 0, 2.0, 1000, -100, 100);
  new TH2F("PC_tof_dE", "T0-PC tof vs dE", 400, -10, 190, 500, 0, 50);

  new TH2F("PC_tof_ang_k", "T0-PC TOF vs bending angle  k beam", 400, -10, 190, 350, -5, 35);
  new TH2F("PC_tof_mom_k", "T0-PC TOF vs momentum  k beam", 400, -10, 190, 500, 0, 2.0);
  new TH2F("PC_tof_fl_k",  "T0-PC TOF vs Fligth length  k beam", 400, -10, 190, 400, 1200, 1700);
  new TH2F("PC_mom_offset_k",  "PC mom vs Calc Offset  k beam", 500, 0, 2.0, 1000, -100, 100);

  new TH2F("PC_tof_ang_pi", "T0-PC TOF vs bending angle  #pi beam", 400, -10, 190, 350, -5, 35);
  new TH2F("PC_tof_mom_pi", "T0-PC TOF vs momentum  #pi beam", 400, -10, 190, 500, 0, 2.0);
  new TH2F("PC_tof_fl_pi",  "T0-PC TOF vs Fligth length  #pi beam", 400, -10, 190, 400, 1200, 1700);
  new TH2F("PC_mom_offset_pi",  "PC mom vs Calc Offset  #pi beam", 500, 0, 2.0, 1000, -100, 100);
  for( int seg=1; seg<=NumOfPCSegments; seg++ ){
    new TH2F(Form("PC%d_tof_ang",seg), Form("T0-PC seg%d TOF vs bending angle",seg), 400, -10, 190, 350, -5, 35);
    new TH2F(Form("PC%d_tof_mom",seg), Form("T0-PC seg%d TOF vs momentum",seg), 400, -10, 190, 500, 0, 2.0);
    new TH2F(Form("PC%d_tof_fl",seg),  Form("T0-PC seg%d TOF vs Fligth length",seg), 400, -10, 190, 400, 1200, 1700);
    new TH2F(Form("PC%d_mom_offset",seg),  Form("PC seg%d mom vs Calc Offset",seg), 500, 0, 2.0, 1000, -100, 100);

    new TH2F(Form("PC%d_tof_ang_k",seg), Form("T0-PC seg%d TOF vs bending angle  k beam",seg), 400, -10, 190, 350, -5, 35);
    new TH2F(Form("PC%d_tof_mom_k",seg), Form("T0-PC seg%d TOF vs momentum  k beam",seg), 400, -10, 190, 500, 0, 2.0);
    new TH2F(Form("PC%d_tof_fl_k",seg),  Form("T0-PC seg%d TOF vs Fligth length  k beam",seg), 400, -10, 190, 400, 1200, 1700);
    new TH2F(Form("PC%d_mom_offset_k",seg),  Form("PC seg%d mom vs Calc Offset  k beam",seg), 500, 0, 2.0, 1000, -100, 100);

    new TH2F(Form("PC%d_tof_ang_pi",seg), Form("T0-PC seg%d TOF vs bending angle  #pi beam",seg), 400, -10, 190, 350, -5, 35);
    new TH2F(Form("PC%d_tof_mom_pi",seg), Form("T0-PC seg%d TOF vs momentum  #pi beam",seg), 400, -10, 190, 500, 0, 2.0);
    new TH2F(Form("PC%d_tof_fl_pi",seg),  Form("T0-PC seg%d TOF vs Fligth length  #pi beam",seg), 400, -10, 190, 400, 1200, 1700);
    new TH2F(Form("PC%d_mom_offset_pi",seg),  Form("PC seg%d mom vs Calc Offset  #pi beam",seg), 500, 0, 2.0, 1000, -100, 100);
  }

  std::cout<<"  for Forward Charge Missing Mass "<<std::endl;
  new TH1F("nCDS_PC_k_p", "ntrack of CDS  K^-(^{3}He, p)", 5, 0, 5);
  new TH1F("nCDS_PC_k_d", "ntrack of CDS  K^-(^{3}He, d)", 5, 0, 5);
  new TH1F("nCDS_PC_pi_p", "ntrack of CDS  #pi^-(^{3}He, p)", 5, 0, 5);
  new TH1F("nCDS_PC_pi_d", "ntrack of CDS  #pi^-(^{3}He, d)", 5, 0, 5);

  new TH1F("nCDS_CVC_k_p", "ntrack of CDS  K^-(^{3}He, p)", 5, 0, 5);
  new TH1F("nCDS_CVC_k_d", "ntrack of CDS  K^-(^{3}He, d)", 5, 0, 5);
  new TH1F("nCDS_CVC_pi_p", "ntrack of CDS  #pi^-(^{3}He, p)", 5, 0, 5);
  new TH1F("nCDS_CVC_pi_d", "ntrack of CDS  #pi^-(^{3}He, d)", 5, 0, 5);

  new TH1F("CDSpid_PC_k_p", "CDS PID  K^-(^{3}He, p)", 8, -1, 7);
  new TH1F("CDSpid_PC_k_d", "CDS PID  K^-(^{3}He, d)", 8, -1, 7);
  new TH1F("CDSpid_PC_pi_p", "CDS PID  #pi^-(^{3}He, p)", 8, -1, 7);
  new TH1F("CDSpid_PC_pi_d", "CDS PID  #pi^-(^{3}He, d)", 8, -1, 7);

  new TH1F("CDSpid_CVC_k_p", "CDS PID  K^-(^{3}He, p)", 8, -1, 7);
  new TH1F("CDSpid_CVC_k_d", "CDS PID  K^-(^{3}He, d)", 8, -1, 7);
  new TH1F("CDSpid_CVC_pi_p", "CDS PID  #pi^-(^{3}He, p)", 8, -1, 7);
  new TH1F("CDSpid_CVC_pi_d", "CDS PID  #pi^-(^{3}He, d)", 8, -1, 7);

  new TH2F("CDS_PC_k_p", "CDS 1/#beta vs mom  K^-(^{3}He, p)", 500, 0, 5, 300, -1.5, 1.5);
  new TH2F("CDS_PC_k_d", "CDS 1/#beta vs mom  K^-(^{3}He, d)", 500, 0, 5, 300, -1.5, 1.5);
  new TH2F("CDS_PC_pi_p", "CDS 1/#beta vs mom  #pi^-(^{3}He, p)", 500, 0, 5, 300, -1.5, 1.5);
  new TH2F("CDS_PC_pi_d", "CDS 1/#beta vs mom  #pi^-(^{3}He, d)", 500, 0, 5, 300, -1.5, 1.5);

  new TH2F("CDS_CVC_k_p", "CDS 1/#beta vs mom  K^-(^{3}He, p)", 500, 0, 5, 300, -1.5, 1.5);
  new TH2F("CDS_CVC_k_d", "CDS 1/#beta vs mom  K^-(^{3}He, d)", 500, 0, 5, 300, -1.5, 1.5);
  new TH2F("CDS_CVC_pi_p", "CDS 1/#beta vs mom  #pi^-(^{3}He, p)", 500, 0, 5, 300, -1.5, 1.5);
  new TH2F("CDS_CVC_pi_d", "CDS 1/#beta vs mom  #pi^-(^{3}He, d)", 500, 0, 5, 300, -1.5, 1.5);

  new TH1F("PC_k_p_mm", "K^{-}(^{3}He, p) Missing Mass", 500, 0, 5);
  new TH1F("PC_k_d_mm", "K^{-}(^{3}He, d) Missing Mass", 500, 0, 5);
  new TH1F("PC_pi_p_mm", "#pi^{-}(^{3}He, p) Missing Mass", 500, 0, 5);
  new TH1F("PC_pi_d_mm", "#pi^{-}(^{3}He, d) Missing Mass", 500, 0, 5);

  new TH1F("PC_k_p_mm_lambda1520", "K^{-}(^{3}He, d) Missing Mass w #Lambda(1520)", 500, 0, 5);
  new TH1F("PC_k_p_mm_lambda", "K^{-}(^{3}He, d) Missing Mass w #Lambda", 500, 0, 5);

  new TH1F("PC_k_d_mm_pipi", "K^{-}(^{3}He, d) Missing Mass w CDS #pi^{+} #pi^{-} tag", 500, 0, 5);
  new TH1F("PC_k_d_mm_ppi",  "K^{-}(^{3}He, d) Missing Mass w CDS p #pi^{-} tag", 500, 0, 5);
  new TH1F("PC_k_d_mm_k",    "K^{-}(^{3}He, d) Missing Mass w CDS K^{-} tag", 500, 0, 5);

  new TH1F("CVC_k_p_mm", "K^{-}(^{3}He, p) Missing Mass", 500, 0, 5);
  new TH1F("CVC_k_d_mm", "K^{-}(^{3}He, d) Missing Mass", 500, 0, 5);
  new TH1F("CVC_pi_p_mm", "#pi^{-}(^{3}He, p) Missing Mass", 500, 0, 5);
  new TH1F("CVC_pi_d_mm", "#pi^{-}(^{3}He, d) Missing Mass", 500, 0, 5);

  new TH1F("IM_pk", "foward p K^{-} Invariat Mass", 1000, 1.4, 2.4);
  new TH1F("IM_lambda_p", "foward #Lambda & CDS p Invariat Mass", 2000, 2.0, 4.0);
  new TH1F("IM_ppim", "foward p #pi^{-} Invariat Mass", 1400, 1.0, 2.4);
  new TH1F("IM_ppip", "foward p #pi^{+} Invariat Mass", 1400, 1.0, 2.4);

  new TH1F("IM_dk", "foward p K^{-} Invariat Mass", 2000, 1.4, 3.4);
  new TH1F("IM_dpim", "foward p #pi^{-} Invariat Mass", 2400, 1.0, 3.4);
  new TH1F("IM_dpip", "foward p #pi^{+} Invariat Mass", 2400, 1.0, 3.4);

  new TH2F("PC_tof_angwADC", "T0-PC TOF vs bending angle", 400, -10, 190, 350, -5, 35);
  new TH2F("PC_tof_momwADC", "T0-PC TOF vs momentum", 400, -10, 190, 500, 0, 2.0);
  new TH2F("PC_tof_flwADC",  "T0-PC TOF vs Fligth length", 400, -10, 190, 400, 1200, 1700);
  new TH2F("PC_mom_offsetwADC",  "PC mom vs Calc Offset", 500, 0, 2.0, 1000, -100, 100);
  new TH2F("PC_tof_dEwADC", "T0-PC tof vs dE", 400, -10, 190, 500, 0, 50);

  new TH2F("PC_tof_angwWVC", "T0-PC TOF vs bending angle", 400, -10, 190, 350, -5, 35);
  new TH2F("PC_tof_momwWVC", "T0-PC TOF vs momentum", 400, -10, 190, 500, 0, 2.0);
  new TH2F("PC_tof_flwWVC",  "T0-PC TOF vs Fligth length", 400, -10, 190, 400, 1200, 1700);
  new TH2F("PC_mom_offsetwWVC",  "PC mom vs Calc Offset", 500, 0, 2.0, 1000, -100, 100);
  new TH2F("PC_tof_dEwWVC", "T0-PC tof vs dE", 400, -10, 190, 500, 0, 50);

  new TH2F("CVC_tof_angwADC", "T0-CVC TOF vs bending angle", 400, -10, 190, 350, -5, 35);
  new TH2F("CVC_tof_momwADC", "T0-CVC TOF vs momentum", 400, -10, 190, 500, 0, 2.0);
  new TH2F("CVC_tof_flwADC",  "T0-CVC TOF vs Fligth length", 400, -10, 190, 400, 1200, 1700);
  new TH2F("CVC_mom_offsetwADC",  "CVC mom vs Calc Offset", 500, 0, 2.0, 1000, -100, 100);
  new TH2F("CVC_tof_dEwADC", "T0-CVC tof vs dE", 400, -10, 190, 500, 0, 50);

  new TH2F("CVC_tof_angwWVC", "T0-CVC TOF vs bending angle", 400, -10, 190, 350, -5, 35);
  new TH2F("CVC_tof_momwWVC", "T0-CVC TOF vs momentum", 400, -10, 190, 500, 0, 2.0);
  new TH2F("CVC_tof_flwWVC",  "T0-CVC TOF vs Fligth length", 400, -10, 190, 400, 1200, 1700);
  new TH2F("CVC_mom_offsetwWVC",  "CVC mom vs Calc Offset", 500, 0, 2.0, 1000, -100, 100);
  new TH2F("CVC_tof_dEwWVC", "T0-CVC tof vs dE", 400, -10, 190, 500, 0, 50);

  std::cout<<" EventAnalysisCheck::InitializeHistogram end ... "<<std::endl;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisCheck *event = new EventAnalysisCheck();
  return (EventTemp*)event;
}
