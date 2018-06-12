#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "BeamSpectrometer.h"

#include "MyTools.h"
#include "TTree.h"
#include "NpipiData.h"

static const double K0_mean = 0.497519;
static const double K0_sigma = 0.00565;
static const double K0_MAX = K0_mean+2.*K0_sigma;
static const double K0_MIN = K0_mean-2.*K0_sigma;

static const double Sm_mean = 1.19694;
static const double Sm_sigma = 0.00499;
static const double Sm_MAX = Sm_mean+2.*Sm_sigma;
static const double Sm_MIN = Sm_mean-2.*Sm_sigma;

static const double Sp_mean = 1.1889;
static const double Sp_sigma = 0.00384;
static const double Sp_MAX = Sp_mean+2.*Sp_sigma;
static const double Sp_MIN = Sp_mean-2.*Sp_sigma;

static const TLorentzVector tgt_lmom = TLorentzVector(0., 0., 0., dMass);

class EventAnalysisNpipi: public EventTemp
{
public:
  EventAnalysisNpipi();
  ~EventAnalysisNpipi();

private:
  int Event_Number_onSpill;
  int CDC_Event_Num;
  TFile *cdcFile;
  TTree *cdcTree;
  CDSTrackingMan *cdstrackMan;
  EventHeader *cdcHeader;

  TFile *rtFile;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackingMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;

  TTree *evTree;
  NpipiData *data;

  //*** for Event Check ***//
  std::ofstream fOfs;
  int fNumNpipi;

public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Clear();
  void Finalize();
};

EventAnalysisNpipi::EventAnalysisNpipi()
  : EventTemp()
{
  std::cout<<"EventAnalysisNpipi::Constractor Call"<<std::endl;
}

EventAnalysisNpipi::~EventAnalysisNpipi()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisNpipi::Initialize( ConfMan *conf )
{
  std::cout<<"EventAnalysisNpipi::Initialization START"<<std::endl;
  Event_Number_onSpill = 0;
  CDC_Event_Num = 0;
  confMan = conf;
  cdcFile = new TFile((TString)confMan->GetCDSTrackFileName());
  cdstrackMan = new CDSTrackingMan();
  cdcHeader = new EventHeader();
  cdcTree = (TTree*)cdcFile-> Get("EventTree");
  cdcTree-> SetBranchAddress("EventHeader", &cdcHeader);
  cdcTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  evTree = new TTree("EventTree", "EventTree");

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  scaMan = new ScalerMan();

  bltrackMan = new BeamLineTrackMan();
  beamSpec = new BeamSpectrometer(conf);

  data = new NpipiData();

  evTree-> Branch("EventHeader", &header);
  evTree-> Branch("NpipiData", &data);

  fOfs.open("evcheck.txt");
  fNumNpipi=0;

  InitializeHistogram();
  std::cout<<"EventAnalysisNpipi::Initialization FINISH"<<std::endl;
}

void EventAnalysisNpipi::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisNpipi::UAna( TKOHitCollection *tko )
{
  TH1F *h1;
  TH2F *h2;

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
    
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) ){
    Clear();
    return true;
  }
  Event_Number_onSpill++;

  TH1F *h1_Kf = (TH1F*)gFile-> Get("Kf_Reduction");
  TH1F *h1_N  = (TH1F*)gFile-> Get("N_Reduction");

  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(0);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(0);

  //*** Check T0 1hit ***//
  int nT0=0;
  double T0time = DBL_MIN;
  HodoscopeLikeHit *T0hit=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      T0time = blMan->T0(i)->ctmean();
      T0hit = blMan->T0(i);
    }
  }
  if( nT0!=1 ){
    Clear();
    return true;
  }
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(1);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(1);

  //*** Check BHD-T0 TOF ***//
  int beam_pid=Beam_Other;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      double tof = T0time-blMan->BHD(i)->ctmean();
      h1 = (TH1F*)gFile-> Get("BHDT0"), h1-> Fill(tof);
      if( header->IsTrig(Trig_Kaon) ) h1 = (TH1F*)gFile-> Get("BHDT0_K"), h1-> Fill(tof);
      if( header->IsTrig(Trig_Pion) ) h1 = (TH1F*)gFile-> Get("BHDT0_pi"), h1-> Fill(tof);
      if( 27<tof && tof<31 ) beam_pid=Beam_Kaon;
      if( beam_pid!=Beam_Kaon && 24<tof && tof<27 ) beam_pid=Beam_Pion;
    }
  }
  if( beam_pid!=Beam_Kaon ){
    Clear();
    return true;
  }
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(2);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(2);

  bltrackMan-> DoTracking(blMan, confMan, true, true);

  int ntrackBLC1=0;
  LocalTrack *trackBLC1=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    double tracktime = bltrackMan->trackBLC1(i)->GetTrackTime();
    h1 = (TH1F*)gFile-> Get("BLC1_time"), h1-> Fill(tracktime);
    if( -20<tracktime && tracktime<20 ){
      ntrackBLC1++;
      trackBLC1 = bltrackMan->trackBLC1(i);
    }
  }
  if( ntrackBLC1!=1 ){
    Clear();
    return true;
  }
  h1 = (TH1F*)gFile-> Get("BLC1_chi2"), h1-> Fill(trackBLC1->chi2all());
  if( trackBLC1->chi2all()>30 ){
    Clear();
    return true;
  }
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(3);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(3);

  int ntrackBLC2=0;
  LocalTrack *trackBLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    double tracktime = bltrackMan->trackBLC2(i)->GetTrackTime();
    h1 = (TH1F*)gFile-> Get("BLC2_time"), h1-> Fill(tracktime);
    if( -20<tracktime && tracktime<20 ){
      ntrackBLC2++;
      trackBLC2 = bltrackMan->trackBLC2(i);
    }
  }
  if( ntrackBLC2!=1 ){
    Clear();
    return true;
  }
  h1 = (TH1F*)gFile-> Get("BLC2_chi2"), h1-> Fill(trackBLC2->chi2all());
  if( trackBLC2->chi2all()>30 ){
    Clear();
    return true;
  }
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(4);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(4);

  beamSpec->TMinuitFit(trackBLC1, trackBLC2, confMan);
  h1 = (TH1F*)gFile-> Get("D5_chi2"), h1-> Fill(beamSpec->chisquare());
  if( beamSpec->chisquare()>30 ){
    Clear();
    return true;
  }
  double D5mom = beamSpec-> mom();
  h1 = (TH1F*)gFile-> Get("D5_mom"), h1->Fill(beamSpec->mom());
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(5);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(5);  

  int ntrackBPC=0;
  LocalTrack *trackBPC=0;
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    double tracktime = bltrackMan->trackBPC(i)->GetTrackTime();
    h1 = (TH1F*)gFile-> Get("BPC_time"), h1-> Fill(tracktime);
    if( -20<tracktime && tracktime<20 ){
      ntrackBPC++;
      trackBPC = bltrackMan->trackBPC(i);
    }
  }
  if( ntrackBPC!=1 ){
    Clear();
    return true;
  }
  h1 = (TH1F*)gFile-> Get("BPC_chi2"), h1-> Fill(trackBPC->chi2all());
  if( trackBPC->chi2all()>30 ){
    Clear();
    return true;
  }
  if( header->IsTrig(Trig_Kf)      ) h1_Kf-> Fill(6);
  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(6);

  TVector3 BPC_FF = trackBPC-> GetPosatZ(-4);
  h1= (TH1F*)gFile->Get("BeamProf"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
  if( header->IsTrig(Trig_Kf) ){
    h1= (TH1F*)gFile->Get("BeamProf_Kf"), h1->Fill(BPC_FF.X(), BPC_FF.Y());
    if( GeomTools::GetID(BPC_FF)==CID_Fiducial ) h1_Kf-> Fill(7);
  }

  //****************************************//
  //*** CDC Tracking File Event Matching ***//
  //****************************************//
  if( CDC_Event_Num>cdcTree->GetEntries() ){
    Clear();
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

  TVector3 T0pos = trackBPC-> GetPosatZ(-110.5);
  std::vector<CDSTrack*> trackPim;
  std::vector<CDSTrack*> trackPip;
  std::vector<CDSTrack*> trackP;
  std::vector<CDSTrack*> trackKm;
  TVector3 good_vtxb;
  TVector3 good_vtxcds;
  double DCA=DBL_MAX;
  bool cdsflag = false;

  //*** for Event Check ***//
  double pim_mass2, pip_mass2;
  double pim_mom0, pip_mom0;
  double pim_CDHtime ,pip_CDHtime;
  for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
    CDSTrack *cdc = cdstrackMan->GoodTrack(i);
    double CDHtime, dis;
    int CDHseg;
    if( !cdc-> GetCDHHit(cdsMan, CDHseg, CDHtime) ) continue;
    double par[5];
    TVector3 vtx;
    cdc-> GetParameters(CID_CDC, par, vtx);
    double mom = cdc->Momentum();
    TVector3 vtxb, vtxcds;
    TrackTools::CalcLineHelixVertex(trackBPC, cdc, vtxb, vtxcds, dis);
    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(T0pos, vtxb, D5mom, parMass[beam_pid], beam_out, beam_tof);
    TVector3 vtxCDH = cdc-> CDHVertex();
    double cdc_dis = MathTools::CalcHelixArc(par, vtxCDH, vtxcds);
    double beta = cdc_dis/(CDHtime-T0time)/(100.*Const);
    double mass2 = mom*mom*(1./(beta*beta)-1);
      
    h1 = (TH1F*)gFile-> Get("T0CDH_TOF"), h1-> Fill(CDHtime-T0time);
    double calc_beta, tofvtxcdc;
    if( !TrackTools::FindMass2(cdc, trackBPC, CDHtime-T0time, D5mom, beam_pid, calc_beta, mass2, tofvtxcdc) ) continue;
    h2 = (TH2F*)gFile->Get("CDS_overbeta_mom"), h2-> Fill(1./calc_beta, mom);
     
    int pid=TrackTools::PID(mom, mass2);
    if( MyTools::IsElectron(calc_beta, mom) ){
      pid=CDS_Other;
      h2 = (TH2F*)gFile->Get("CDS_overbeta_mom_e"), h2-> Fill(1./calc_beta, mom);
    }
    cdc-> SetPID(pid);

    if( GeomTools::GetID(vtxb)==CID_Fiducial ){
      h2 = (TH2F*)gFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);
      if( mom<0 ) h1 = (TH1F*)gFile->Get("CDS_mass2_minus"), h1-> Fill(mass2);
      else h1 = (TH1F*)gFile->Get("CDS_mass2_plus"), h1-> Fill(mass2);
      if( pid==CDS_PiMinus ){
	pim_mass2 = mass2;
	pim_mom0 = mom;
	pim_CDHtime = CDHtime;
	h1 = (TH1F*)gFile->Get("CDS_mass2_pim"), h1-> Fill(mass2);
	h2 = (TH2F*)gFile->Get("CDS_mass2_mom_pim"), h2-> Fill(mass2, mom);
      }
      if( pid==CDS_PiPlus  ){
	pip_mass2 = mass2;
	pip_mom0 = mom;
	pip_CDHtime = CDHtime;
	h1 = (TH1F*)gFile->Get("CDS_mass2_pip"), h1-> Fill(mass2);
	h2 = (TH2F*)gFile->Get("CDS_mass2_mom_pip"), h2-> Fill(mass2, mom);
      }
      if( pid==CDS_Kaon    ){
	h1 = (TH1F*)gFile->Get("CDS_mass2_K"), h1-> Fill(mass2);
	h2 = (TH2F*)gFile->Get("CDS_mass2_mom_K"), h2-> Fill(mass2, mom);
      }
      if( pid==CDS_Proton  ){
	h1 = (TH1F*)gFile->Get("CDS_mass2_p"), h1-> Fill(mass2);
	h2 = (TH2F*)gFile->Get("CDS_mass2_mom_p"), h2-> Fill(mass2, mom);
      }
    }

    double cds_tof, tmpl;
    cdc-> CalcVertexTimeLength(T0pos, trackBPC->GetMomDir(), cdsMass[pid], vtxb, vtxcds, cds_tof, tmpl, true);
    double dis2;
    TrackTools::CalcLineHelixVertex(trackBPC, cdc, vtxb, vtxcds, dis2);
    if( dis<DCA ){
      DCA = dis2;
      good_vtxb = vtxb;
      good_vtxcds = vtxcds;
      cdsflag = true;
    }
    if( pid==CDS_PiMinus ) trackPim.push_back(cdc);
    else if( pid==CDS_Kaon ) trackKm.push_back(cdc);
    else if( pid==CDS_PiPlus ) trackPip.push_back(cdc);
    else if( pid==CDS_Proton ) trackP.push_back(cdc);
  }
  if( !cdsflag ){
    Clear();
    return true;
  }
  h2 = (TH2F*)gFile-> Get("Vtx_XY"), h2-> Fill(good_vtxb.X(), good_vtxb.Y());
  h2 = (TH2F*)gFile-> Get("Vtx_ZX"), h2-> Fill(good_vtxb.Z(), good_vtxb.X());
  h2 = (TH2F*)gFile-> Get("Vtx_ZY"), h2-> Fill(good_vtxb.Z(), good_vtxb.Y());
  if( GeomTools::GetID(good_vtxb)==CID_Fiducial ){
    h2 = (TH2F*)gFile-> Get("Vtx_XY_wtar"), h2-> Fill(good_vtxb.X(), good_vtxb.Y());
    h2 = (TH2F*)gFile-> Get("Vtx_ZX_wtar"), h2-> Fill(good_vtxb.Z(), good_vtxb.X());
    h2 = (TH2F*)gFile-> Get("Vtx_ZY_wtar"), h2-> Fill(good_vtxb.Z(), good_vtxb.Y());
  }

  if( header->IsTrig(Trig_Neutral) ) h1_N-> Fill(7);
  bool vtx_flag = false;
  if( GeomTools::GetID(good_vtxb)==CID_Fiducial ) vtx_flag = true;
  if( GeomTools::GetID(good_vtxb)==CID_Fiducial && header->IsTrig(Trig_Neutral) ) h1_N-> Fill(8);

  bool K0_flag = false;
  for( int i=0; i<trackPim.size(); i++ ){
    for( int j=0; j<trackPip.size(); j++ ){
      TVector3 vtx_pim, vtx_pip;
      if( TrackTools::Calc2HelixVertex(trackPim[i], trackPip[j], vtx_pim, vtx_pip) ){
        TVector3 pim_mom, pip_mom;
        if( trackPim[i]-> GetMomentum(vtx_pim, pim_mom, true, true) &&
            trackPip[j]-> GetMomentum(vtx_pip, pip_mom, true, true) ){
          TLorentzVector pim_lmom, pip_lmom;
          pim_lmom.SetVectM(pim_mom, piMass);
          pip_lmom.SetVectM(pip_mom, piMass);
          TVector3 vtx_mean = 0.5*(vtx_pim+vtx_pip);
          TVector3 mom_sum = pim_mom+pip_mom;

          TVector3 vtxb, vtxcds;
          double dltmp, dca;
	  MathTools::LineToLine(vtx_mean, mom_sum.Unit(), T0pos, trackBPC->GetMomDir(), dltmp, dca, vtxcds, vtxb);

          double beam_out, beam_tof;
	  ELossTools::CalcElossBeamTGeo(T0pos, vtxb, D5mom, kpMass, beam_out, beam_tof);
          TVector3 beam_mom = trackBPC->GetMomDir();
          beam_mom.SetMag(beam_out);
          TLorentzVector beam_lmom;
	  beam_lmom.SetVectM(beam_mom, kpMass);

	  if( K0_MIN<(pim_lmom+pip_lmom).M() && (pim_lmom+pip_lmom).M()<K0_MAX ) K0_flag=true;
          if( GeomTools::GetID(vtxb)==CID_Fiducial ){
            h1 = (TH1F*)rtFile-> Get("CDS_IM_pipi"), h1-> Fill((pim_lmom+pip_lmom).M());
          }
        }
      }
    }
  }

  std::vector<HodoscopeLikeHit*> NC_hit[8];
  for( int i=0; i<blMan->nNC(); i++ ){
    if( blMan->NC(i)->CheckRange() ){
      int seg = blMan->NC(i)-> seg();
      int lay = (seg-1)/14;
      NC_hit[lay].push_back(blMan->NC(i));
    }
  }

  int nBVC=0;
  for( int i=0; i<blMan->nBVC(); i++ ){
    if( blMan->BVC(i)->CheckRange() ) nBVC++;
  }

  int nCVC=0;
  for( int i=0; i<blMan->nCVC(); i++ ){
    if( blMan->CVC(i)->CheckRange() ) nCVC++;
  }

  bool NCflag = false;
  TVector3 nc_pos;
  int NCseg=-1;
  double NC_hitpos;
  double NCtime = DBL_MAX;
  double NCdE = 0.0;
  for( int lay=0; lay<8; lay++ ){
    for( int i=0; i<NC_hit[lay].size(); i++ ){
      if( NC_hit[lay][i]->emean()>8.0 ){
	if( NCdE<NC_hit[lay][i]->emean() ){
	  NCtime = NC_hit[lay][i]->ctmean();
	  NCseg  = NC_hit[lay][i]->seg();
	  NCdE   = NC_hit[lay][i]->emean();
	  NC_hitpos = NC_hit[lay][i]->hitpos();
	  NCflag = true;
	}
      }
    }
    if( NCflag ){
      confMan-> GetGeomMapManager()-> GetGPos(CID_NC, NCseg, nc_pos);
      nc_pos.SetY(NC_hitpos);
      break;
    }
  }

  if( !NCflag ){
    Clear();
    return true;
  }

  double fl = (nc_pos-good_vtxb).Mag();
  double beam_out, beam_tof;
  ELossTools::CalcElossBeamTGeo(T0pos, good_vtxb, D5mom, parMass[beam_pid], beam_out, beam_tof);
  double NC_beta = fl/((NCtime-T0time-beam_tof)*100.*Const);

  bool Sp_flag = false;
  bool Sm_flag = false;
  if( vtx_flag ){
    h1 = (TH1F*)gFile->Get("NC_hitpos"), h1->Fill(NC_hitpos);

    h1 = (TH1F*)gFile->Get("NC_overbeta"), h1-> Fill(1/NC_beta);
    h2 = (TH2F*)gFile->Get("NC_overbeta_dE"), h2-> Fill(1/NC_beta, NCdE);
    if( NCdE>8. ){
      h1 = (TH1F*)gFile->Get("NC_overbeta_8MeVee"), h1-> Fill(1/NC_beta);
      if( NC_beta<0.9 ){
	double nc_mom = nMass*NC_beta/sqrt(1-NC_beta*NC_beta);
	TVector3 n_mom = (nc_pos-good_vtxb);
	n_mom.SetMag(nc_mom);
	TLorentzVector n_lmom;
	n_lmom.SetVectM(n_mom, nMass);

	TLorentzVector fn_lmom;
	fn_lmom.SetVectM(n_mom, nMass);
	
	double beam_out, beam_tof;
	ELossTools::CalcElossBeamTGeo(T0pos, good_vtxb, D5mom, parMass[beam_pid], beam_out, beam_tof);
	TVector3 beam_mom = trackBPC->GetMomDir();
	beam_mom.SetMag(beam_out);
	TLorentzVector beam_lmom;
	beam_lmom.SetVectM(beam_mom, parMass[beam_pid]);

	TLorentzVector kn_lmom = beam_lmom+tgt_lmom-n_lmom;

	for( int i=0; i<trackPim.size(); i++ ){
	  TVector3 vtxb, vtxcds;
	  double dis;
	  if( !TrackTools::CalcLineHelixVertex(trackBPC, trackPim[i], vtxb, vtxcds, dis) ) continue;

	  TVector3 p;
	  if( !trackPim[i]->GetMomentum(vtxcds, p, true, true) ) continue;
	  TLorentzVector lmom;
	  lmom.SetVectM(p, piMass);

	  if( GeomTools::GetID(good_vtxb)==CID_Fiducial ) h1 = (TH1F*)gFile-> Get("Npim_IM"), h1-> Fill((fn_lmom+lmom).M());
	  if( Sm_MIN<(fn_lmom+lmom).M() && (fn_lmom+lmom).M()<Sm_MAX ) Sm_flag = true;
	}

	for( int i=0; i<trackPip.size(); i++ ){
	  TVector3 vtxb, vtxcds;
	  double dis;
	  if( !TrackTools::CalcLineHelixVertex(trackBPC, trackPip[i], vtxb, vtxcds, dis) ) continue;

	  TVector3 p;
	  if( !trackPip[i]->GetMomentum(vtxcds, p, true, true) ) continue;
	  TLorentzVector lmom;
	  lmom.SetVectM(p, piMass);

	  if( GeomTools::GetID(good_vtxb)==CID_Fiducial ) h1 = (TH1F*)gFile-> Get("Npip_IM"), h1-> Fill((fn_lmom+lmom).M());
	  if( Sp_MIN<(fn_lmom+lmom).M() && (fn_lmom+lmom).M()<Sp_MAX ) Sp_flag = true;
	}

	if( GeomTools::GetID(good_vtxb)==CID_Fiducial ){
	  h1 = (TH1F*)gFile-> Get("KN_MM"), h1-> Fill(kn_lmom.M());
	  if( trackPim.size()>0 ) h1 = (TH1F*)gFile-> Get("KN_MM_pim"), h1-> Fill(kn_lmom.M());
	  if( trackPip.size()>0 ) h1 = (TH1F*)gFile-> Get("KN_MM_pip"), h1-> Fill(kn_lmom.M());
	  if( trackKm.size()>0 )  h1 = (TH1F*)gFile-> Get("KN_MM_km"), h1-> Fill(kn_lmom.M());
	  if( trackP.size()>0 )   h1 = (TH1F*)gFile-> Get("KN_MM_p"), h1-> Fill(kn_lmom.M());

	  if( K0_flag ) h1 = (TH1F*)gFile-> Get("KN_MM_K0"), h1-> Fill(kn_lmom.M());
	  if( Sp_flag ) h1 = (TH1F*)gFile-> Get("KN_MM_Sp"), h1-> Fill(kn_lmom.M());
	  if( Sm_flag ) h1 = (TH1F*)gFile-> Get("KN_MM_Sm"), h1-> Fill(kn_lmom.M());
	  if( Sm_flag || Sp_flag ) h1 = (TH1F*)gFile-> Get("KN_MM_S"), h1-> Fill(kn_lmom.M());	
	}
      }
    }
  }

  //***** for n pi+ pi- Event *****//
  if( nBVC==0 && nCVC==0 && header->IsTrig(Trig_Neutral) ){
    if( NC_beta<0.9 && NCdE>8.){
      if( vtx_flag ) h1_N-> Fill(9);
      if( trackPim.size()==1 && trackPip.size()==1 ){
	if( vtx_flag ) h1_N-> Fill(10);
	TVector3 vtxPim, vtxPip;
	if( !TrackTools::Calc2HelixVertex(trackPim[0], trackPip[0], vtxPim, vtxPip) ){
	  std::cout<<"  !!! TrackTools::Calc2HelixVertex false !!!"<<std::endl;
	  Clear();
	  return true;
	}
	TVector3 pim_mom, pip_mom;
	if( !trackPim[0]->GetMomentum(vtxPim, pim_mom, true, true) ){
	  std::cout<<"  !!! pi- track::GetMomentum false !!!"<<std::endl;
	  Clear();
	  return true;
	}
	if( !trackPip[0]->GetMomentum(vtxPip, pip_mom, true, true) ){
	  std::cout<<"  !!! pi+ track::GetMomentum false !!!"<<std::endl;
	  Clear();
	  return true;
	}
	
	TVector3 vtx_mean = 0.5*(vtxPim+vtxPip);
	TVector3 mom_sum = pim_mom+pip_mom;
	TVector3 vtx_cds, vtxb;
	double dltmp, dca;

	MathTools::LineToLine(vtx_mean, mom_sum.Unit(), T0pos, trackBPC->GetMomDir(), dltmp, dca, vtx_cds, vtxb);
	ELossTools::CalcElossBeamTGeo(T0pos, vtxb, D5mom, parMass[beam_pid], beam_out, beam_tof);
	
	TVector3 vtxBeamPim, vtxBeamPip, vtxPimBeam, vtxPipBeam;
	TLorentzVector fn_lmom_pim, fn_lmom_pip;

	double dis_pim, dis_pip;
	TrackTools::CalcLineHelixVertex(trackBPC, trackPim[0], vtxBeamPim, vtxPimBeam, dis_pim);
	TrackTools::CalcLineHelixVertex(trackBPC, trackPip[0], vtxBeamPip, vtxPipBeam, dis_pip);
	TVector3 pim_mom2, pip_mom2;
	if( !trackPim[0]->GetMomentum(vtxPimBeam, pim_mom2, true, true) ){
	  std::cout<<"  !!! pi- track::GetMomentum return false 222 !!!"<<std::endl;
	  return true;
	}
	if( !trackPip[0]->GetMomentum(vtxPipBeam, pip_mom2, true, true) ){
	  std::cout<<"  !!! pi+ track::GetMomentum return false 222 !!!"<<std::endl;
	  return true;
	}
	TLorentzVector pim_lmom2, pip_lmom2;
	pim_lmom2.SetVectM(pim_mom2, piMass);
	pip_lmom2.SetVectM(pip_mom2, piMass);
	
	//*** add for n pi+ pi- ***//
	bool good_pim_flag;
	TVector3 beam_mom = trackBPC->GetMomDir();
	TLorentzVector beam_lmom;
	double nc_mom;
	TVector3 n_mom;
	TLorentzVector n_lmom;
	double NC_beta_c;
	if( (vtxBeamPim-vtxPimBeam).Mag()<(vtxPipBeam-vtxBeamPip).Mag() ){
	  ELossTools::CalcElossBeamTGeo(T0pos, vtxBeamPim, D5mom, parMass[beam_pid], beam_out, beam_tof);

	  NC_beta_c = (nc_pos-vtxBeamPim).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	  nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	  n_mom = (nc_pos-vtxBeamPim);
	  n_mom.SetMag(nc_mom);
	  n_lmom.SetVectM(n_mom, nMass);
	  beam_mom.SetMag(beam_out);
	  beam_lmom.SetVectM(beam_mom, kpMass);

	  good_pim_flag=true;
	}
	else{
	  ELossTools::CalcElossBeamTGeo(T0pos, vtxBeamPip, D5mom, parMass[beam_pid], beam_out, beam_tof);

	  NC_beta_c = (nc_pos-vtxBeamPip).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	  nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	  n_mom = (nc_pos-vtxBeamPip);
	  n_mom.SetMag(nc_mom);
	  n_lmom.SetVectM(n_mom, nMass);
	  beam_mom.SetMag(beam_out);
	  beam_lmom.SetVectM(beam_mom, kpMass);

	  good_pim_flag=false;
	}

	TLorentzVector pim_lmom, pip_lmom;
	pim_lmom.SetVectM(pim_mom, piMass);
	pip_lmom.SetVectM(pip_mom, piMass);
	
	data-> setBeamLmom(beam_lmom);
	data-> setFNLmom(n_lmom);
	data-> setPimLmom(pim_lmom);
	data-> setPipLmom(pip_lmom);
	data-> setVtxCDS(vtx_cds);
	data-> setVtxBeam(vtxb);
	data-> setVtxPim(vtxPim);
	data-> setVtxPip(vtxPip);
	
	double NC_beta_pim = (nc_pos-vtxPimBeam).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	double nc_mom_pim = nMass*NC_beta_c/sqrt(1-NC_beta_pim*NC_beta_pim);
	TVector3 n_mom_pim = (nc_pos-vtxPimBeam);
	n_mom_pim.SetMag(nc_mom_pim);
	TLorentzVector n_lmom_pim;
	n_lmom_pim.SetVectM(n_mom_pim, nMass);
	
	double NC_beta_pip = (nc_pos-vtxPipBeam).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	double nc_mom_pip = nMass*NC_beta_c/sqrt(1-NC_beta_pip*NC_beta_pip);
	TVector3 n_mom_pip = (nc_pos-vtxPipBeam);
	n_mom_pip.SetMag(nc_mom_pip);
	TLorentzVector n_lmom_pip;
	n_lmom_pip.SetVectM(n_mom_pip, nMass);
	
	double NC_beta_pimb = (nc_pos-vtxBeamPim).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	double nc_mom_pimb = nMass*NC_beta_c/sqrt(1-NC_beta_pim*NC_beta_pim);
	TVector3 n_mom_pimb = (nc_pos-vtxBeamPim);
	n_mom_pimb.SetMag(nc_mom_pim);
	TLorentzVector n_lmom_pim_beam;
	n_lmom_pim_beam.SetVectM(n_mom_pimb, nMass);
	
	double NC_beta_pipb = (nc_pos-vtxBeamPip).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	double nc_mom_pipb = nMass*NC_beta_c/sqrt(1-NC_beta_pip*NC_beta_pip);
	TVector3 n_mom_pipb = (nc_pos-vtxBeamPip);
	n_mom_pipb.SetMag(nc_mom_pip);
	TLorentzVector n_lmom_pip_beam;
	n_lmom_pip_beam.SetVectM(n_mom_pipb, nMass);
	
	data-> setFNLmomPim(n_lmom_pim);
	data-> setFNLmomPip(n_lmom_pip);
	data-> setFNLmomPimB(n_lmom_pim_beam);
	data-> setFNLmomPipB(n_lmom_pip_beam);
	data-> setPimLmom2(pim_lmom2);
	data-> setPipLmom2(pip_lmom2);
	data-> setVtxPimBeam(vtxPimBeam);
	data-> setVtxPipBeam(vtxPipBeam);
	data-> setVtxBeamPim(vtxBeamPim);
	data-> setVtxBeamPip(vtxBeamPip);

	//*** for Event Check ***// 
	bool fiducial_flag = false;
	TVector3 good_vtx;
	if( dis_pim<dis_pip ){
	  if( GeomTools::GetID(vtxBeamPim)==CID_Fiducial ){
	    fiducial_flag = true;
	    good_vtx = vtxBeamPim;
	  }
	}
	else{
	  if( GeomTools::GetID(vtxBeamPip)==CID_Fiducial ){
	    fiducial_flag = true;
	    good_vtx = vtxBeamPip;
	  }
	}
	if( fiducial_flag ){
	  TVector3 momdir = trackBPC-> GetMomDir();
	  TVector3 T0pos = trackBPC-> GetPosatZ(-110.5);
	  fOfs<<"===== Event Number : "<<Event_Number<<"====="<<std::endl;
	  fOfs<<"> T0  time : "<<T0time<<std::endl;
	  fOfs<<">> T0 raw data   seg "<<T0hit->seg()<<std::endl;
	  fOfs<<">> adc up : "<<T0hit->adcu()<<"  down : "<<T0hit->adcd()<<std::endl;
	  fOfs<<">> tdc up : "<<T0hit->tdcu()<<"  down : "<<T0hit->tdcd()<<std::endl;
	  fOfs<<">> T0 ctu : "<<T0hit->ctu()<<"   ctd : "<<T0hit->ctd()<<std::endl;
	  fOfs<<">> T0 dE_u : "<<T0hit->eu()<<"   dE_d : "<<T0hit->ed()<<std::endl;
	  fOfs<<"> D5 mom : "<<D5mom<<std::endl;
	  fOfs<<"> pi- mass2 : "<<pim_mass2<<" mom : "<<pim_mom0<<" CDH time : "<<pim_CDHtime<<std::endl;
	  fOfs<<">     vtx beam : ("<<vtxBeamPim.X()<<", "<<vtxBeamPim.Y()<<", "<<vtxBeamPim.Z()<<")"<<std::endl;
	  fOfs<<">     vtx CDS  : ("<<vtxPimBeam.X()<<", "<<vtxPimBeam.Y()<<", "<<vtxPimBeam.Z()<<")"<<std::endl;
	  fOfs<<"> pi+ mass2 : "<<pip_mass2<<" mom : "<<pip_mom0<<" CDH time : "<<pip_CDHtime<<std::endl;
	  fOfs<<">     vtx beam : ("<<vtxBeamPip.X()<<", "<<vtxBeamPip.Y()<<", "<<vtxBeamPip.Z()<<")"<<std::endl;
	  fOfs<<">     vtx CDS  : ("<<vtxPipBeam.X()<<", "<<vtxPipBeam.Y()<<", "<<vtxPipBeam.Z()<<")"<<std::endl;
	  fOfs<<"> Mom Dir ("<<momdir.X()<<", "<<momdir.Y()<<", "<<momdir.Z()<<")"<<std::endl;
	  fOfs<<"> T0 pos ("<<T0pos.X()<<", "<<T0pos.Y()<<", "<<T0pos.Z()<<")"<<std::endl;
	  fOfs<<"> NC hit pos : "<<nc_pos.X()<<", "<<nc_pos.Y()<<", "<<nc_pos.z()<<std::endl;
	  if( good_pim_flag ) fOfs<<"> Start pos by pi-"<<std::endl;
	  else  fOfs<<"> Start pos by pi+"<<std::endl;
	  fOfs<<"> T0-vtx tof : "<<beam_tof<<std::endl;
	  fOfs<<"> NC time : "<<NCtime<<std::endl;
	  fOfs<<"======================================"<<std::endl;
	  fNumNpipi++;
	}
      }
      
      evTree->Fill();
    }
  }
  else if( NCdE>8. && vtx_flag ) h1_N-> Fill(11);
  else if( vtx_flag ) h1_N-> Fill(12);
  Clear();
  return true;
}

void EventAnalysisNpipi::Clear()
{
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  bltrackMan-> Clear();
  beamSpec-> Clear();
  data-> clear();
}

void EventAnalysisNpipi::InitializeHistogram()
{
  rtFile-> cd();
  new TH1F("Kf_Reduction", "K/f Event Reduction", 10, -0.5, 9.5);
  new TH1F("N_Reduction", "Neutral Event Reduction", 15, -0.5, 14.5);

  new TH1F("BHDT0",    "BHD-T0", 1000, 0.0, 100);
  new TH1F("BHDT0_K",  "BHD-T0", 1000, 0.0, 100);
  new TH1F("BHDT0_pi", "BHD-T0", 1000, 0.0, 100);

  new TH1F("BLC1_time", "BLC1 track time", 10000, -200, 800);
  new TH1F("BLC1_chi2", "BLC1 chi-square", 1000, 0.0, 100);

  new TH1F("BLC2_time", "BLC2 track time", 10000, -200, 800);
  new TH1F("BLC2_chi2", "BLC2 chi-square", 1000, 0.0, 100);

  new TH1F("BPC_time", "BPC track time", 10000, -200, 800);
  new TH1F("BPC_chi2", "BPC chi-square", 1000, 0.0, 100);

  new TH1F("D5_mom",  "Beam Momentum by D5", 1000,  0.5, 1.5);
  new TH1F("D5_chi2", "D5 chisquare", 1000,  0.0, 100);

  new TH2F("BeamProf",    "Beam Profile at FF",        1000, -30, 30, 1000, -30, 30);
  new TH2F("BeamProf_Kf", "Beam Profile at FF w/ K/f", 1000, -30, 30, 1000, -30, 30);

  new TH1F("T0CDH_TOF", "T0-CDH TOF", 1100, -10, 100);
  new TH2F("CDS_overbeta_mom",   "CDS mass2 vs mom", 1000, 0.0, 10.0, 1500, -1.5, 1.5);
  new TH2F("CDS_overbeta_mom_e", "CDS mass2 vs mom", 1000, 0.0, 10.0, 1500, -1.5, 1.5);
  new TH2F("CDS_mass2_mom",      "CDS mass2 vs mom", 5000, -0.5, 9.5, 1500, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_pim",  "CDS mass2 vs mom", 5000, -0.5, 9.5, 1500, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_pip",  "CDS mass2 vs mom", 5000, -0.5, 9.5, 1500, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_K",    "CDS mass2 vs mom", 5000, -0.5, 9.5, 1500, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_p",    "CDS mass2 vs mom", 5000, -0.5, 9.5, 1500, -1.5, 1.5);
  new TH1F("CDS_mass2_minus", "CDS mass2 mom<0", 5000, -0.5, 9.5);
  new TH1F("CDS_mass2_plus",  "CDS mass2 mom>0", 5000, -0.5, 9.5);
  new TH1F("CDS_mass2_pim",   "CDS mass2 mom>0", 5000, -0.5, 9.5);
  new TH1F("CDS_mass2_pip",   "CDS mass2 mom>0", 5000, -0.5, 9.5);
  new TH1F("CDS_mass2_K",     "CDS mass2 mom>0", 5000, -0.5, 9.5);
  new TH1F("CDS_mass2_p",     "CDS mass2 mom>0", 5000, -0.5, 9.5);

  new TH2F("Vtx_XY", "Vertex XY Plane", 1000, -30, 30, 1000, -30, 30);
  new TH2F("Vtx_ZX", "Vertex ZX Plane", 1000, -30, 30, 1000, -30, 30);
  new TH2F("Vtx_ZY", "Vertex ZY Plane", 1000, -30, 30, 1000, -30, 30);

  new TH2F("Vtx_XY_wtar", "Vertex XY Plane", 1000, -30, 30, 1000, -30, 30);
  new TH2F("Vtx_ZX_wtar", "Vertex ZX Plane", 1000, -30, 30, 1000, -30, 30);
  new TH2F("Vtx_ZY_wtar", "Vertex ZY Plane", 1000, -30, 30, 1000, -30, 30);

  new TH1F("CDS_IM_pipi", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);
  new TH1F("CDS_IM_ppim", "CDS p #pi^{-} IM", 1000, 0.0, 1.0);
  new TH1F("Npim_IM", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM", "n #pi^{+} IM", 1000, 1.0, 2.0);

  new TH1F("KN_MM",    "d(K^{-}, n)",   1000, 1.0, 2.0);
  new TH1F("KN_MM_pim",  "d(K^{-}, n)", 1000, 1.0, 2.0);
  new TH1F("KN_MM_pip",  "d(K^{-}, n)", 1000, 1.0, 2.0);
  new TH1F("KN_MM_km",  "d(K^{-}, n)",  1000, 1.0, 2.0);
  new TH1F("KN_MM_p",  "d(K^{-}, n)",   1000, 1.0, 2.0);
  new TH1F("KN_MM_d",  "d(K^{-}, n)",   1000, 1.0, 2.0);

  new TH1F("KN_MM_K0", "d(K^{-}, n)",   1000, 1.0, 2.0);
  new TH1F("KN_MM_Sp", "d(K^{-}, n)",   1000, 1.0, 2.0);
  new TH1F("KN_MM_Sm", "d(K^{-}, n)",   1000, 1.0, 2.0);
  new TH1F("KN_MM_S",  "d(K^{-}, n)",   1000, 1.0, 2.0);

  new TH1F("NC_hitpos",          "NC hit pos", 10000, -500, 500);
  new TH1F("NC_overbeta",        "NC 1/#beta", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_8MeVee", "NC 1/#beta", 10000, 0.0, 10.0);
  new TH2F("NC_overbeta_dE",     "NC 1/#beta vs dE", 10000, 0.0, 10.0, 1000, 0.0, 100.);
}

void EventAnalysisNpipi::Finalize()
{
  std::cout<<"EventAnalysisNpipi Finish   Event Number : "<<Event_Number<<"  Event Number on Spill : "<<Event_Number_onSpill<<std::endl;
  gFile->Write();
  confMan-> SaveParams();
  confMan-> SaveCode();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;

  delete bltrackMan;
  delete beamSpec;

  //*** for Event Check ***//
  fOfs<<" NumNpipi/Analysis Event : "<<fNumNpipi<<"/"<<Event_Number<<std::endl;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisNpipi *event = new EventAnalysisNpipi();
  return (EventTemp*)event;
}
