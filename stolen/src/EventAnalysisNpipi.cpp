#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "HistManBeamAna.h"
#include "CDSTrackingMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "BeamSpectrometer.h"

#include "TTree.h"
#include "NpipiData.h"

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

public:
  void Initialize( ConfMan *conf );
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

  new TH2F("CDS_mass2_mom",      "CDS mass2 vs mom", 5000, -0.5, 9.5, 1500, -1.5, 1.5);
  new TH1F("NC_overbeta",        "NC 1/#beta", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_8MeVee", "NC 1/#beta", 10000, 0.0, 10.0);

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





  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

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
  //*** only Neutral Trigger Analysis ***//
  if( header->IsTrig(Trig_Neutral) ){
    int nT0=0;
    double T0time = DBL_MIN;
    for( int i=0; i<blMan->nT0(); i++ ){
      if( blMan->T0(i)->CheckRange() ){
	nT0++;
	T0time = blMan->T0(i)->ctmean();
      }
    }
    if( nT0!=1 ){
      Clear();
      return true;
    }

    int beam_pid=Beam_Other;
    for( int i=0; i<blMan->nBHD(); i++ ){
      if( blMan->BHD(i)->CheckRange() ){
	double tof = T0time-blMan->BHD(i)->ctmean();
	if( 27<tof && tof<31 ) beam_pid=Beam_Kaon;
	if( beam_pid!=Beam_Kaon && 24<tof && tof<27 ) beam_pid=Beam_Pion;
      }
    }
    if( beam_pid!=Beam_Kaon ){
      Clear();
      return true;
    }

    bltrackMan-> DoTracking(blMan, confMan, true, true);

    int ntrackBLC1=0;
    LocalTrack *trackBLC1=0;
    for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
      double tracktime = bltrackMan->trackBLC1(i)->GetTrackTime();
      if( -20<tracktime && tracktime<20 ){
	ntrackBLC1++;
	trackBLC1 = bltrackMan->trackBLC1(i);
      }
    }
    if( ntrackBLC1!=1 ){
      Clear();
      return true;
    }
    if( trackBLC1->chi2all()>30 ){
      Clear();
      return true;
    }

    int ntrackBLC2=0;
    LocalTrack *trackBLC2=0;
    for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
      double tracktime = bltrackMan->trackBLC2(i)->GetTrackTime();
      if( -20<tracktime && tracktime<20 ){
	ntrackBLC2++;
	trackBLC2 = bltrackMan->trackBLC2(i);
      }
    }
    if( ntrackBLC2!=1 ){
      Clear();
      return true;
    }
    if( trackBLC2->chi2all()>30 ){
      Clear();
      return true;
    }

    beamSpec->TMinuitFit(trackBLC1, trackBLC2, confMan);
    if( beamSpec->chisquare()>30 ){
      Clear();
      return true;
    }
    double D5mom = beamSpec-> mom();

    int ntrackBPC=0;
    LocalTrack *trackBPC=0;
    for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
      double tracktime = bltrackMan->trackBPC(i)->GetTrackTime();
      if( -20<tracktime && tracktime<20 ){
	ntrackBPC++;
	trackBPC = bltrackMan->trackBPC(i);
      }
    }
    if( ntrackBPC!=1 ){
      Clear();
      return true;
    }
    if( trackBPC->chi2all()>30 ){
      Clear();
      return true;
    }

    TVector3 T0pos = trackBPC-> GetPosatZ(-110.5);
    std::vector<CDSTrack*> trackPim;
    std::vector<CDSTrack*> trackPip;
    TVector3 good_vtxb;
    TVector3 good_vtxcds;
    double DCA=DBL_MAX;
    bool cdsflag = false;
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

      double beta_calc, tofvtxcdc;
      std::cout<<"input mass2="<<mass2<<std::endl;
      if( !TrackTools::FindMass2C(cdc, trackBPC, CDHtime-T0time, D5mom, beam_pid, beta, mass2, tofvtxcdc) ) continue;
      std::cout<<"mass2="<<mass2<<std::endl;
      h2 = (TH2F*)gFile-> Get("CDS_mass2_mom"), h2-> Fill(mass2, mom);

      int pid=TrackTools::PID(mom, mass2);
      cdc-> SetPID(pid);

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
      else if( pid==CDS_PiPlus ) trackPip.push_back(cdc);
    }
    if( !cdsflag ){
      Clear();
      return true;
    }

    std::vector<HodoscopeLikeHit*> NC_hit[8];
    for( int i=0; i<blMan->nNC(); i++ ){
      if( blMan->NC(i)->CheckRange() ){
	int seg = blMan->NC(i)-> seg();
	int lay = (seg-1)/14;
	NC_hit[lay].push_back(blMan->NC(i));
      }
    }

    bool NCflag = false;
    TVector3 nc_pos;
    int NCseg=-1;
    double NCtime = DBL_MAX;
    double NCdE = 0.0;
    for( int lay=0; lay<8; lay++ ){
      if( NC_hit[lay].size()>0 ){
	for( int i=0; i<NC_hit[lay].size(); i++ ){
	  if( NCtime>NC_hit[lay][i]->ctmean() ){
	    NCtime = NC_hit[lay][i]->ctmean();
	    NCseg  = NC_hit[lay][i]->seg();
	    NCdE   = NC_hit[lay][i]->emean();
	    NCflag = true;
	  }
	}
	confMan-> GetGeomMapManager()-> GetGPos(CID_NC, NCseg, nc_pos);
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
    h1 = (TH1F*)gFile->Get("NC_overbeta"), h1-> Fill(1/NC_beta);
    if( NCdE>8. ) h1 = (TH1F*)gFile->Get("NC_overbeta_8MeVee"), h1-> Fill(1/NC_beta);

    if( NC_beta<0.9 && NCdE>8.){
      if( trackPim.size()==1 && trackPip.size()==1 ){
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
	TVector3 beam_mom = trackBPC->GetMomDir();
	beam_mom.SetMag(beam_out);
	TLorentzVector beam_lmom;
	beam_lmom.SetVectM(beam_mom, parMass[beam_pid]);

	double NC_beta_c = (nc_pos-vtxb).Mag()/((NCtime-T0time-beam_tof)*100.*Const);
	double nc_mom = nMass*NC_beta_c/sqrt(1-NC_beta_c*NC_beta_c);
	TVector3 n_mom = (nc_pos-vtxb);
	n_mom.SetMag(nc_mom);
	TLorentzVector n_lmom;
	n_lmom.SetVectM(n_mom, nMass);
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

	evTree->Fill();
      }
    }
  }
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
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisNpipi *event = new EventAnalysisNpipi();
  return (EventTemp*)event;
}
