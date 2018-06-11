// H. Asano Mar. 15th 2018
// This code is for the drift time correlation calibration for beamline drift chambers.
// 
// input: raw data, GainMapBL, X-Tmap.
// output: root files which contain histograms of correlatation of drift time of each chamebr with several cuts

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
#include "GlobalVariables.h"
#include "Tools.h"

#include <TParameter.h>
#include <TROOT.h>

class EventAnalysisDCTimeCorr: public EventTemp
{
public:
  EventAnalysisDCTimeCorr();
  ~EventAnalysisDCTimeCorr();
private:
  TFile *rtFile;
  //CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  EventHeader *header;
  ScalerMan *scaMan;
  double scainit[40];
  double scaend[40];
  bool INIT;
  bool SCAOVERFLOW[40];
  Time_t Time;
  long t_init;
  int spillinit;
  int spillfini;

  int Verbosity_;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  bool EndOfAnEvent(bool flag=true);
  void UTime( int time ){ Time = time; };
};

EventAnalysisDCTimeCorr::EventAnalysisDCTimeCorr()
  : EventTemp()
{
}

EventAnalysisDCTimeCorr::~EventAnalysisDCTimeCorr()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisDCTimeCorr::Initialize( ConfMan *conf )
{
  Verbosity_ = 0;

  gROOT->SetStyle("Plain");
#if 0
  std::cout << " Enter EventAnalysisDCTimeCorr::Initialize " << std::endl;
#endif
  confMan = conf;
  rtFile = new TFile( conf->GetOutFileName().c_str() , "recreate" );
  INIT=true;
  spillinit=spillfini=-1;
  for(int i=0;i<40;i++){
    SCAOVERFLOW[i]=false;
    scaend[i]=0;
  }
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!! object of EventHeader is NULL " << std::endl; return; }
  //cdsMan = new CDSHitMan();
  //if( cdsMan==NULL ){ std::cerr << "!!!! object of cdsMan is NULL" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!! object of BeamLineHitMan is NULL" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!! object of ScalerMan is NULL" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan==NULL ){ std::cerr << "!!!! object of bltrackMan is NULL" << std::endl; return; }

}

void EventAnalysisDCTimeCorr::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisDCTimeCorr::USca " << std::endl;
#endif
  rtFile->cd();
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  for( int i=0; i<8; i++ )
    std::cout<<sca[i]<<"\t";
  std::cout<<std::endl;
#endif
  //  std::cout<<"nsca: "<<nsca<<std::endl;
  if(INIT&&sca[0]<5){
    t_init=Time;    
    INIT=false;
    for( int i=0; i<scaMan->nsca(); i++ )   scainit[i] = scaMan->sca(i)->val();
  }
  if(INIT){
    scaMan->Clear();
    return;
  }
  for( int i=0; i<scaMan->nsca(); i++ ){
    //   int val = scaMan->sca(i)->val()-scaend[i];
    //   TString name = scaMan->sca(i)->name();
    //   Tools::Fill2D(name,Block_Event_Number,val);
    if(scaend[i]>9.9e+07&&scaend[i]>scaMan->sca(i)->val()) SCAOVERFLOW[i]=true;
    scaend[i] = scaMan->sca(i)->val();
    if(SCAOVERFLOW[i]) scaend[i]+=1.0e+08;
    Tools::H1(scaMan->sca(i)->name(),Block_Event_Number,3000,-0.5,2999.5,scaend[i]);
  }
  Tools::H1("Time",Block_Event_Number,3000,-0.5,2999.5,Time-t_init);
  scaMan->Clear();
}


bool EventAnalysisDCTimeCorr::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return EndOfAnEvent();
    if( status==2 ) return EndOfAnEvent(false);}
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number <<"  " <<Time<<"  "<<scaend[10]<<std::endl;

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  DetectorList *dlist=DetectorList::GetInstance();
//#########

  header->Convert( tko, confMan );
  if(spillinit<0) spillinit=header->spillt(0);
  spillfini=header->spillt(0);
  blMan->Convert( tko, confMan );
  //cdsMan->Convert( tko, confMan );
  rtFile->cd();

  // Trigger Pattern
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      Tools::H1("Pattern",i,21,-0.5,20.5);
      Tools::H1(Form("TPattern%d",i),val,4196,-0.5,4195.5);
    }
  }


  for(int i=0;i<20;i++)
    if(header->trigmode(i))  Tools::H1("TrigMode",i, 21,-0.5,20.5); 

  int DAQFLAG=0;
  for(int i=0;i<20;i++)
    if(header->trigmode2(i)) 
      {
	Tools::H1("DAQMode",i,21,-0.5,20.5); 
	DAQFLAG++;
      }
  if(DAQFLAG==0)  Tools::H1("DAQMode",18,21,-0.5,20.5); 
  if(DAQFLAG==1)  Tools::H1("DAQMode",19,21,-0.5,20.5); 
  if(DAQFLAG>1)   Tools::H1("DAQMode",20,21,-0.5,20.5); 
  
  if(header->IsTrig(Trig_Cosmic)) return EndOfAnEvent();

  
  //Asano memo:
  //DoTracking( BeamLineHitMan *blMan, ConfMan *conf, const bool &SEMILOCAL=false, const bool &TIMING=false );
  //flag SEMILOCAL -> if true, string option "linear" is added (default is false).
  //flag TIMING ->if false, string option "notiming" is added (default). In addion to this, FDC is added "noslope". If "notiming" is declared in BeamLineTrackingMan, the RMS of track time (= drift time corrected by pair plane calibration) < 20 (15 ?) is required.
  
  bltrackMan-> DoTracking(blMan, confMan);


  // BeamLine Chamber
  // # of layers
  //  BLC1a,b ( 8 layers UU'VV'UU'VV')
  //  BLC2a,b ( 8 layers UU'VV'UU'VV')
  //  BPC     ( 8 layers XX'YY'XX'YY')
  //  FDC     ( 6 layers VV'XX'UU')
  //
  
  //prepare BVC info. for FDC1 study
  int nBVC=0;
  int BVCseg=-1;

  for( int i=0; i<blMan->nBVC();i++){
    HodoscopeLikeHit *hit=blMan->BVC(i);
    if(hit->CheckRange()){
      nBVC++;
      BVCseg=hit->seg();
    }
  }
  Tools::H1( Form("MulBVC") , nBVC ,8,0.5,8.5, "multiplicity");
  
  const int ndc=6;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		  CID_BLC2a,CID_BLC2b,
		  CID_BPC,
		  CID_FDC1
  };
  bool all_1hit[ndc]={true,true,true,true,true,true};//each layer has 1 hit 
  int dchit[8][ndc]={{0}};
  
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    //no cuts

    for( int layer=1; layer<=nlays; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( cid, layer );
      int tmpmul=blMan->nBLDC(cid,layer);
      Tools::H1( Form("Mul%s_%d",name,layer),tmpmul, nwire+1, -0.5, nwire+0.5,"Multiplicity" ); 
      if(tmpmul!=1) all_1hit[idc]=false;
      dchit[layer-1][idc]=tmpmul;
      for( int i=0; i<tmpmul; i++ ){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        int wire = hit->wire();
        int tdc = hit->tdc();
        Tools::H1( Form("T%s_%d",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
        Tools::H1( Form("T%s_%d_%d",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
        Tools::H1( Form("HitPat%s_%d",name,layer) , wire , nwire, 0.5, nwire+0.5,"wire#");
  
        double dt = hit->dt();
        Tools::H1( Form("d%s_%d",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
        Tools::H1( Form("d%s_%d_%d",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
        
        //FDC1 study
        if(cid == CID_FDC1){
          //with trigger cuts
          if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) ){
            Tools::H1( Form("d%s_%d_beam_kf_trig",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            Tools::H1( Form("d%s_%d_%d_beam_kf_trig",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          }
           
          //with BVC multiplicity cut
          if(nBVC==1){
            Tools::H1( Form("d%s_%d_BVC1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            Tools::H1( Form("d%s_%d_%d_BVC1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          }
          //both
          if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) && (nBVC==1)  ){
            Tools::H1( Form("d%s_%d_BVC1_beam_kf_trig",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            Tools::H1( Form("d%s_%d_%d_BVC1_beam_kf_trig",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          }

        }

        if(tmpmul==1){
          Tools::H1( Form("T%s_%d_mul1",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
          Tools::H1( Form("T%s_%d_%d_mul1",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
          Tools::H1( Form("d%s_%d_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          Tools::H1( Form("d%s_%d_%d_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          //FDC1 study
          if(cid == CID_FDC1){
            //with trigger cuts
            if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) ){
              Tools::H1( Form("d%s_%d_beam_kf_trig_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_beam_kf_trig_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }

            //with BVC multiplicity cut
            if(nBVC==1){
              Tools::H1( Form("d%s_%d_BVC1_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_BVC1_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }
            //both
            if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) && (nBVC==1)  ){
              Tools::H1( Form("d%s_%d_BVC1_beam_kf_trig_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_BVC1_beam_kf_trig_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }
          }
        }


	
        if(layer%2==1)
        {
          for( int i2=0; i2< (blMan->nBLDC(cid,layer+1)); i2++ ){
            ChamberLikeHit *hit2 = blMan->BLDC(cid,layer+1,i2);
            int wire2 = hit2->wire();
            //	      int tdc2 = hit2->tdc();
            Tools::H2( Form("WirePat%s_%d",name,layer), wire,wire2, nwire, 0.5, nwire+0.5, nwire, 0.5, nwire+0.5 );
            Tools::SetXTitleH2( Form("WirePat%s_%d",name,layer), Form("Layer%d wire#",layer));
            Tools::SetYTitleH2( Form("WirePat%s_%d",name,layer), Form("Layer%d wire#",layer+1));
          }
        }
      }//for mul
    }//for layer
    
    //requiring other layers' hits
    for( int layer=1; layer<=nlays; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( cid, layer );
      int tmpmul=blMan->nBLDC(cid,layer);
      int OtherNhitLayer=0;
      for(int ilr=1;ilr<=8;ilr++){
        if(ilr!=layer && dchit[ilr-1][idc]!=0) OtherNhitLayer++;
      }
      
      for( int i=0; i<tmpmul; i++ ){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        int wire = hit->wire();
        int tdc = hit->tdc();
        double dt = hit->dt();
        if(all_1hit[idc]){  
          Tools::H1( Form("T%s_%d_all1hit",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
          Tools::H1( Form("T%s_%d_%d_all1hit",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
          Tools::H1( Form("HitPat%s_%d_all1hit",name,layer) , wire , nwire, 0.5, nwire+0.5,"wire#");
  
          Tools::H1( Form("d%s_%d_all1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          Tools::H1( Form("d%s_%d_%d_all1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          //FDC1 study
          if(cid == CID_FDC1){
            //with trigger cuts
            if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) ){
              Tools::H1( Form("d%s_%d_beam_kf_trig_all1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_beam_kf_trig_all1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }

            //with BVC multiplicity cut
            if(nBVC==1){
              Tools::H1( Form("d%s_%d_BVC1_all1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_BVC1_all1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }
            //both
            if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) && (nBVC==1)  ){
              Tools::H1( Form("d%s_%d_BVC1_beam_kf_trig_all1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_BVC1_beam_kf_trig_all1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }
          }
	      }

        if(OtherNhitLayer>0){
          Tools::H1( Form("T%s_%d_other1hit",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
          Tools::H1( Form("T%s_%d_%d_other1hit",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
          Tools::H1( Form("HitPat%s_%d_other1hit",name,layer) , wire , nwire, 0.5, nwire+0.5,"wire#");
  
          Tools::H1( Form("d%s_%d_other1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          Tools::H1( Form("d%s_%d_%d_other1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          if(cid == CID_FDC1){
            //with trigger cuts
            if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) ){
              Tools::H1( Form("d%s_%d_beam_kf_trig_other1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_beam_kf_trig_other1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }

            //with BVC multiplicity cut
            if(nBVC==1){
              Tools::H1( Form("d%s_%d_BVC1_other1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_BVC1_other1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }
            //both
            if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) && (nBVC==1)  ){
              Tools::H1( Form("d%s_%d_BVC1_beam_kf_trig_other1hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              Tools::H1( Form("d%s_%d_%d_BVC1_beam_kf_trig_other1hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            }
          }

          if(tmpmul==1){
            Tools::H1( Form("T%s_%d_other1hit_mul1",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
            Tools::H1( Form("T%s_%d_%d_other1hit_mul1",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
            Tools::H1( Form("d%s_%d_other1hit_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            Tools::H1( Form("d%s_%d_%d_other1hit_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          
            if(cid == CID_FDC1){
              //with trigger cuts
              if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) ){
                Tools::H1( Form("d%s_%d_beam_kf_trig_other1hit_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
                Tools::H1( Form("d%s_%d_%d_beam_kf_trig_other1hit_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              }

              //with BVC multiplicity cut
              if(nBVC==1){
                Tools::H1( Form("d%s_%d_BVC1_other1hit_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
                Tools::H1( Form("d%s_%d_%d_BVC1_other1hit_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              }
              //both
              if( (header->trigmode(Mode_Beam) || header->trigmode(Mode_Kf)) && (nBVC==1)  ){
                Tools::H1( Form("d%s_%d_BVC1_beam_kf_trig_other1hit_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
                Tools::H1( Form("d%s_%d_%d_BVC1_beam_kf_trig_other1hit_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
              }
            }
          }
        }
        
        if(OtherNhitLayer>3){
          Tools::H1( Form("T%s_%d_other4hit",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
          Tools::H1( Form("T%s_%d_%d_other4hit",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
          Tools::H1( Form("d%s_%d_other4hit",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          Tools::H1( Form("d%s_%d_%d_other4hit",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          if(tmpmul==1){
            Tools::H1( Form("T%s_%d_other4hit_mul1",name,layer) , tdc ,2048,-0.5,2047.5, "TDC Ch.");
            Tools::H1( Form("T%s_%d_%d_other4hit_mul1",name,layer,wire) , tdc  ,2048,-0.5,2047.5,"TDC Ch.");
            Tools::H1( Form("d%s_%d_other4hit_mul1",name,layer) , dt , 6000, -1500,1500, "drift time [nsec.]" );
            Tools::H1( Form("d%s_%d_%d_other4hit_mul1",name,layer,wire) , dt , 6000, -1500,1500, "drift time [nsec.]" );
          }
        }
      }//for mul
    }//for layer
  }//for idc
  
  //fill # of reconstructed track.
  for(int idc=0;idc<6;idc++){
    const char* name= dlist->GetName(dccid[idc]).data();
    int ntrack = bltrackMan->ntrackBLDC(dccid[idc]);
    Tools::H1( Form("nTrack%s",name),ntrack,20,0,20,"# of tracks");
  }

  if(Verbosity_){
    std::cout << "Event:" << Event_Number << std::endl;
  }
  for(int idc=0;idc<6;idc++){
    const char* name=dlist->GetName(dccid[idc]).data();
    int ntrack = bltrackMan->ntrackBLDC(dccid[idc]);
    if(Verbosity_){
      std::cout << name << "  ntrack:" << ntrack << std::endl;
    }
    //requires only 1 track.
    if(ntrack!=1) continue;
    for(int itrk=0;itrk<ntrack;itrk++){
      LocalTrack *ltrack = bltrackMan->trackBLDC(dccid[idc],itrk);

      //hit loop
      int nhit = ltrack->nhit();
      for(int ihit=0;ihit<nhit;ihit++){
        ChamberLikeHit *chit = ltrack->hit(ihit);
        Tools::H1( Form("%s_dt_%d_track",name, chit->layer() ), chit->dt(),1000,-200,800,"drift time [nsec.]");
        Tools::H1( Form("%s_dt_%d_%d_track",name, chit->layer(), chit->wire()), chit->dt(),1000,-200,800,"drift time [nsec.]");
        Tools::H1( Form("%s_res_%d_track",name, chit->layer() ), chit->resl(),500,-0.1,0.1,"hit residual [cm]");
        Tools::H1( Form("%s_res_%d_%d_track",name, chit->layer(), chit->wire()), chit->resl(),500,-0.1,0.1,"hit residual [cm]");
      }

      //cluster loop
      for(int ixy=0;ixy<2;ixy++){
        int ncls = ltrack->ncluster(ixy);
        for(int icls=0;icls<ncls;icls++){
          BLDCCluster *cluster = ltrack->cluster(ixy,icls);
          if( cluster->nhit()!=2 ) continue;// necessary for pair analysis 
          int layer = ltrack->hit(0)->layer();
          double timesub = cluster->GetTimeSub();
          double timemean = cluster->GetTimeMean();
          double ctimemean = cluster->GetCTime();//corrected time mean (?)
          Tools::H2( Form("%s%d_TimeSubVSTimeMean_track",name,layer),timesub,timemean,500,-250,250,500,-100,400, "Time Sub. [nsec]", "Time mean [nsec.]");
          Tools::H2( Form("%s%d_TimeSubVSTCimeMean_track",name,layer),timesub,ctimemean,500,-250,250,500,-100,400, "Time Sub. [nsec]", "Corrected Time mean [nsec.]" );
        }
      } //ixy
    }//for itrk
  }//for idc


  return EndOfAnEvent();  
}

void EventAnalysisDCTimeCorr::Finalize()
{
  std::cout << " Enter EventAnalysisDCTimeCorr::Finalize " << std::endl;
  for(int i=0;i<40;i++){
    Tools::H1("Scaler",i,41,-0.5,40.5,scaend[i]-scainit[i]);
    std::cout<<i<<"  "<<scaend[i]<<std::endl;
  }
  TParameter<long> inittime("Time_init",t_init);
  inittime.Write();
  std::cout<<inittime.GetVal()<<std::endl;
  TParameter<long> finitime("Time_fini",Time);
  finitime.Write();
  std::cout<<finitime.GetVal()<<std::endl;
  TParameter<int> initspill("Spill_init",spillinit);
  initspill.Write();
  TParameter<int> finispill("Spill_fini",spillfini);
  finispill.Write();

  rtFile->cd();
  gFile->Write();
  gFile->Close();
  delete blMan;
  //delete cdsMan;
  delete bltrackMan;
  delete header;
}

bool EventAnalysisDCTimeCorr::EndOfAnEvent(bool flag){
  header->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  return flag;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisDCTimeCorr *event = new EventAnalysisDCTimeCorr();
  return (EventTemp*)event;
}
