#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h> 

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "SimDataMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "MathTools.h"
#include "Particle.h"
#include "Tools.h"
#include "TrackTools.h"

//<<<<<<< UserSimDatG4.cpp
//#######Need knucl lib file!! Set your dir!! ###
//#include "KnuclRootData.h"
//#include "/home/had/thashi/g4src/knucl4.10/include/KnuclRootData.h"
#include "/home/yamaga/geant/knucl4/include/KnuclRootData.h"
//#######################################
#define DEBUG 0
#define TRACKING 1

const double ADC_TH_BVC = 0.1; // MeV
const double ADC_TH_CDH =-1.0; // MeV
const double ADC_TH_CVC = 0.1; // MeV
const double ADC_TH_NC = -0.1; // MeV
//=======
#include "KnuclRootData.h"
//>>>>>>> 1.17

int main( int argc, char **argv )
{
  gSystem->Load(".lib/KnuclRootData_cc.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");

  if(argc != 4 ){
    std::cout << "Plese set Conffile Outputfile Inputfile  "<< std::endl;
    std::cout << argv[0] << " <conf file> <out file> <in file> " << std::endl;
    return 0;
  }

  std::string confFile,outFile,inFile;
//<<<<<<< UserSimDatG4.cpp
  if( argv[1][0] != '!' ) confFile = argv[1];
  if( argv[2][0] != '!' ) outFile = argv[2];
  if( argv[3][0] != '!' ) inFile = argv[3];

  //#######Need knucl lib file!! Set your dir!! ###
  //  gSystem->Load("~/work/ana/geant/knucl3/libKnuclRootData.so");
  gSystem->Load("/home/yamaga/geant/knucl4/build/KnuclRootData_cc.so");
  //#######################################
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("./src/lib/libAll.so");
//=======
  if( argv[1][0] != '!' ){
    confFile = argv[1];
    if( confFile.find_last_of(".conf")==std::string::npos){
      std::cout<<"ConfFile should be *.conf"<<std::endl;
      return 0;
    }
  }
  if( argv[2][0] != '!' ){
    outFile = argv[2];
    if( outFile.find_last_of(".root")==std::string::npos ){
      std::cout<<"Oututfile should be *.root"<<std::endl;
      return 0;
    }
  }
  if( argv[3][0] != '!' ){
    inFile = argv[3];
    if( inFile.find_last_of(".root")==std::string::npos ){
      std::cout<<"Inputfile should be *.root"<<std::endl;
      return 0;
    }
  }
//>>>>>>> 1.17

  TFile *f = new TFile( inFile.c_str() );
//<<<<<<< UserSimDatG4.cpp
  TTree *evtree = (TTree*)f->Get("tree");

  /*** conf file for new parameter ***/
//=======
//>>>>>>> 1.17
  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();
  SimDataMan *simMan = new SimDataMan(f, conf);
  CDSHitMan *cdsMan = new CDSHitMan();
  BeamLineHitMan *blMan = new BeamLineHitMan();

  TFile *of = new TFile( outFile.c_str(), "recreate" ); 

  int evn = simMan-> nEvent();
  if( conf->GetStopEvNum()>0 ){
    if( evn>conf->GetStopEvNum() ) evn=conf->GetStopEvNum();
  }

//<<<<<<< UserSimDatG4.cpp
  //  Header* header = 0;
  DetectorData *detectorData = 0;
  MCData* mcData = 0;
  ReactionData* reacData = 0;

  /*** declaration of classes ***/  
  //  evtree->SetBranchAddress( "Header", &header );
  evtree->SetBranchAddress( "DetectorData", &detectorData );
  evtree->SetBranchAddress( "MCData", &mcData );
  evtree->SetBranchAddress( "ReactionData", &reacData );

  /*** make output file & make tree ***/

  TFile *of = new TFile(outFile.c_str(),"recreate");

  TTree *otree = new TTree( "EventTree", "EventTree");
  //--- copy input to output ---//  
  //otree->Branch( "header", &header );
  //otree->Branch( "detectorData", &detectorData );
  //  otree->Branch( "MCData", &mcData );
  //  otree->Branch( "ReactionData", &reacData );
  //--- attach new branch on output ---//

  CDSHitMan* cdsMan = 0;
  BeamLineHitMan* blMan = 0;
  CDSTrackingMan* cdsTrackMan = 0;
  BeamLineTrackMan* blTrackMan = new BeamLineTrackMan();
  otree->Branch( "CDSHitMan", &cdsMan );
  otree->Branch( "BeamLineHitMan", &blMan );
  otree->Branch( "CDSTrackingMan", &cdsTrackMan );
  //otree->Branch( "BeamLineTrackMan", &blTrackMan );

  int nev = evtree->GetEntries();
  int fev = 0;
  int tev = 0;
  //  nev = 10000;
  std::cerr<<"# of events = "<<nev<<std::endl;

  int t0=time(0);
  for(int iev = 0; iev<nev; iev++ ){
    evtree->GetEvent(iev);
    SimDataMan *simMan = new SimDataMan(conf);

    int status = conf->CheckEvNum( iev, 1 );
    if(status){ 
      if( status==1 ) continue;
      gFile->Write();
      gFile->Close();      
      if( status==2 ) return 0;
    }    
    if( iev%200==0 )
      std::cout << " Event : " << iev << " Time : " << time(0)-t0 << " Fill# "<< fev<<"  track "<<tev<<  std::endl;

    //###########################//
    //### check detector hits ###//
    //###########################//
    // check hit size
    int nT0 = 0;
    int nCDH = 0;
    int nBLC2a = 0;
    int nBLC2b = 0;
    int nCDC = 0;
    int nNC = 0;
    int nCVC = 0;
    int nPC = 0;
    int nBVC = 0;
    double time0 = 0;
    for( int i=0; i<detectorData->detectorHitSize(); i++ ){
      DetectorHit *detectorHit = detectorData->detectorHit(i);
      int detectorID = detectorHit->detectorID();
      int parentID = detectorHit->parentID();
      int pdgCode  = detectorHit->pdg();      
      if( detectorID==CID_T0 && parentID==0 && pdgCode==321){
        nT0++;
        time0 = detectorHit->tdc();
      }
      else if( detectorID==CID_CDH ) nCDH++;
      else if( detectorID==CID_BLC2a ) nBLC2a++;
      else if( detectorID==CID_BLC2b ) nBLC2b++;
      else if( detectorID==CID_CDC ) nCDC++;
      else if( detectorID==CID_NC && detectorHit->adc()>ADC_TH_NC ) nNC++;
      else if( detectorID==CID_CVC && detectorHit->adc()>ADC_TH_CVC ) nCVC++;
      else if( detectorID==CID_PC ) nPC++;
      else if( detectorID==CID_BVC && detectorHit->adc()>ADC_TH_BVC ) nBVC++;
    }
#if DEBUG
    std::cerr<<"###### event number  = "<<iev<<std::endl;
    std::cerr<<"# of T0 hits   = "<<nT0<<std::endl;
    std::cerr<<"# of CDH hits  = "<<nCDH<<std::endl;
    std::cerr<<"# of BLC1 hits = "<<nBLC2a<<std::endl;
    std::cerr<<"# of BLC2 hits = "<<nBLC2b<<std::endl;
    std::cerr<<"# of CDC hits  = "<<nCDC<<std::endl;
    std::cerr<<"# of NC hits  = "<<nNC<<std::endl;
    std::cerr<<"# of CVC hits  = "<<nCVC<<std::endl;
    std::cerr<<"# of PC hits  = "<<nPC<<std::endl;
    std::cerr<<"# of BVC hits  = "<<nBVC<<std::endl;
#endif

    //if(!nT0==1) continue;
    //if(nCDC==0) continue;
    //if(nCDH==0) continue;
    //if(Ncdh<2) continue;
#if 0
    if ( !nT0 ) { continue;
#else
      if(nT0==0)
      {
        //	simMan->SetT0Hit( 3, 0, 0, 0 ); // 20130629
        time0=0;
      }
#endif     
      //##########################//
      //### fill detector hits ###//
      //##########################//
      DetectorHit *detectorHit = 0;
      double val0=0;
      double val1=0;
      double val2=0;
      std::map<int,int> CDHpdg;
      std::map<int,int> CDHparent;
      std::map<int,int> CDHid;
      std::map<int,bool> CDHin;
      for(int i=0; i<detectorData->detectorHitSize(); i++ ){
        detectorHit = detectorData->detectorHit(i);
        //--- T0 ---//
        if(detectorHit->detectorID()==CID_T0){
          Int_t seg = detectorHit->channelID()+1;
          val0=detectorHit->tdc()-time0;
          val1=detectorHit->adc();
          val2=detectorHit->pos().y()/10.;
          //	 if(val1==0) val1=3.;//MeV
          simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
        }
        //--- CDH ---//
        if(detectorHit->detectorID()==CID_CDH&& detectorHit->adc()>ADC_TH_CDH){
          Int_t seg = detectorHit->channelID()+1;
          val0=detectorHit->tdc()+time0;
          val1=detectorHit->adc();
          val2=detectorHit->pos().z()/10.+gRandom->Gaus(0,2);
          CDHpdg[seg]=detectorHit->pdg();
          CDHparent[seg]=detectorHit->parentID();
          CDHid[seg]=detectorHit->trackID();
          int itr;
          for(itr=0;itr<mcData->trackSize();itr++){
            int tmpid=mcData->track(itr)->trackID();
            if(tmpid==CDHid[seg]) break;
          }
          if(mcData->track(itr)->vertex().Perp()<151.) CDHin[seg]=true;
          else CDHin[seg]=false;
          if(val1==0) continue;//MeV
          simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
        }
        if(detectorHit->detectorID()==CID_IH||
            detectorHit->detectorID()==CID_PC
          ){	
          Int_t seg = detectorHit->channelID()+1;
          val0=detectorHit->tdc()+time0;
          val1=detectorHit->adc();
          val2=detectorHit->pos().z()/10.;
          if(val1==0) continue;//MeV
          simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
        }
        if(detectorHit->detectorID()==CID_CDC ||
            detectorHit->detectorID()==CID_BPC ||
            detectorHit->detectorID()==CID_FDC1 ||
            detectorHit->detectorID()==CID_BLC2a ||
            detectorHit->detectorID()==CID_BLC2b ){
          Int_t layer = detectorHit->layerID()+1;
          Int_t wire = detectorHit->channelID()+1;
          val0=fabs(detectorHit->dx()/10.);
#if 0
          if(detectorHit->detectorID()==CID_CDC){
            TVector3 tmppos=detectorHit->pos()*0.1;
            TVector3 tmpdir=detectorHit->momentum().Unit();
            TVector3 wpos,wdir;
            if(conf->GetCDCWireMapManager()->GetWirePosDir(layer,wire,wpos,wdir)){
              TVector3 tmp1,tmp2,tmp3;
              double dist1,dist2;
              MathTools::LineToLine(tmppos,tmpdir,wpos,wdir,0.,dist1,tmp1,tmp2);
              // //	    MathTools::PointToLine(tmppos,wpos,wdir,dist,tmp1);
              conf->GetCDCWireMapManager()->GetWirePosDir(layer,wire,wpos,wdir);
              MathTools::LineToLine(tmppos,tmpdir,wpos,wdir,0.,dist2,tmp2,tmp3);
              std::cout<<"================layer  "<<layer<<"  ,knucl  "<<val0<<"  , k18ana  "<<dist1<<"  "<<dist2<<std::endl;      
              tmppos.Print();
              tmp1.Print();
              tmp2.Print();
              //	    val0=dist;
            }
          }
#endif
          val1=0;//time offset
          val2=detectorHit->momentum().Pt();
          simMan->SetChamberLikeHit( detectorHit->detectorID(), layer, wire, val0, val1,val2 );
        }
        //--- NC ---//
        if(detectorHit->detectorID()==CID_NC && detectorHit->adc()>ADC_TH_NC){
          Int_t seg = detectorHit->channelID()+1;
          val0=detectorHit->tdc()+time0;
          val1=detectorHit->adc();
          val2=detectorHit->pos().y()/10.;
          simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
        }
        //--- CVC ---//
        if(detectorHit->detectorID()==CID_CVC && detectorHit->adc()>ADC_TH_CVC){
          Int_t seg = detectorHit->channelID()+1;
          val0=detectorHit->tdc()+time0;
          val1=detectorHit->adc();
          val2=detectorHit->pos().z()/10.;
          simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
        }
        //--- BVC ---//
        if(detectorHit->detectorID()==CID_BVC && detectorHit->adc()>ADC_TH_BVC){
          Int_t seg = detectorHit->channelID()+1;
          val0=detectorHit->tdc()+time0;
          val1=detectorHit->adc();
          val2=detectorHit->pos().y()/10.;
          simMan->SetHodoscopeLikeHit(detectorHit->detectorID(), seg, val0, val1, val2 );
        }
      } // for(int i=0; i<detectorData->detectorHitSize(); i++ ){

      cdsMan = simMan->GetCDSHit();
      cdsTrackMan->Execute( cdsMan, conf );
      blMan = simMan->GetBeamLineHit();
      //blTrackMan->DoTracking( blMan, conf );

#if DEBUG
      std::cout<<"######Filled data#############"<<std::endl;
      for(int i=0; i<blMan->nT0(); i++){
        HodoscopeLikeHit *hod = blMan->T0(i);
        std::cerr<<"T0 : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
      }
      for(int i=0; i<cdsMan->nCDH(); i++){
        HodoscopeLikeHit *hod = cdsMan->CDH(i);
        std::cerr<<"CDH : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
      }
      for(int i=1; i<=NumOfBLCLayers; i++){
        for(int j=0; j<blMan->nBLC2a(i); j++){
          ChamberLikeHit *chm = blMan->BLC2a(i, j);
          std::cerr<<"BLC2a : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<std::endl;
        }
      }
      for(int i=1; i<=NumOfBLCLayers; i++){
        for(int j=0; j<blMan->nBLC2b(i); j++){
          ChamberLikeHit *chm = blMan->BLC2b(i, j);
          std::cerr<<"BLC2b : "<<chm->hid()<<" "<<chm->layer()<<" "<<chm->wire()<<" "<<chm->dl()<<std::endl;
        }
      }
      for(int i=1; i<=NumOfCDCLayers; i++){
        for(int j=0; j<cdsMan->nCDC(i); j++){
          CDCHit *cdc = cdsMan->CDC(i, j);
          std::cerr<<"CDC : "<<cdc->hid()<<" "<<cdc->layer()<<" "<<cdc->wire()<<" "<<cdc->dl()<<std::endl;
        }
      }
      for(int i=0; i<blMan->nNC(); i++){
        HodoscopeLikeHit *hod = blMan->NC(i);
        std::cerr<<"NC : "<<hod->hid()<<" "<<hod->seg()<<" ("<<hod->seg()/16+1<<","<<hod->seg()%16<<") "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
      }
      for(int i=0; i<blMan->nCVC(); i++){
        HodoscopeLikeHit *hod = blMan->CVC(i);
        std::cerr<<"CVC : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
      }
      for(int i=0; i<blMan->nPC(); i++){
        HodoscopeLikeHit *hod = blMan->PC(i);
        std::cerr<<"PC : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
      }
      for(int i=0; i<blMan->nBVC(); i++){
        HodoscopeLikeHit *hod = blMan->BVC(i);
        std::cerr<<"BVC : "<<hod->hid()<<" "<<hod->seg()<<" "<<hod->ctmean()<<" "<<hod->emean()<<std::endl;
      }
#endif

#if DEBUG
      std::cerr<<"ev-num = "<<iev<<" cdc-track-num = "<<cdsTrackMan->nGoodTrack()<<" bpc-track-num = "<<blTrackMan->ntrackBPC()<<std::endl;
#endif
      otree->Fill();

      int reacID=reacData->ReactionID();
      //    header->SetSimReacID(reacID);

      bool T0OK=false;
      int t0seg=-1;
      for(int i=0;i<blMan->nT0();i++)
        if(TMath::Abs(blMan->T0(i)->ctmean())<0.1){
          T0OK=true;
          t0seg=blMan->T0(i)->seg();
        }
      if(T0OK&&blTrackMan->ntrackBPC()==1){
        //-----------generated data analysis
        double beam_mom;
        for(int itr=0;itr<mcData->trackSize();itr++){
          int parid=mcData->track(itr)->parentTrackID();
          int pdgid=mcData->track(itr)->pdgID();
          TVector3 mom=mcData->track(itr)->momentum();
          if(pdgid==321&&parid==0){
            beam_mom=mom.Mag()*0.001;
          }      
        }

        //##################################
        // Beam Analysis
        //##################################
        LocalTrack *bpctrack=blTrackMan->trackBLDC(CID_BPC,0);
        TVector3 t0gpos;
        conf->GetGeomMapManager()->GetPos(CID_T0,0,t0gpos);
        TVector3 bpcpos=bpctrack->GetPosatZ(0);		  
        TVector3 t0pos=bpctrack->GetPosatZ(t0gpos.Z()-0.5);		  
        TVector3 bpcdir=bpctrack->GetMomDir().Unit();

        pBeam *pbeam=new pBeam();
        pbeam->SetPID(Beam_Kaon);
        pbeam->SetT0seg(t0seg);
        pbeam->SetT0Time(0.);
        pbeam->SetT0Pos(t0pos);
        pbeam->SetBPCPos(bpcpos);
        pbeam->SetBPCDir(bpcdir);
        pbeam->SetMomentum(beam_mom);
        //##################################
        // CDS 1 track analysis
        //##################################
        int ngood=cdsTrackMan->nGoodTrack();
        for( int it1=0; it1<ngood; it1++ ){      
          //      std::cout<<it1<<"  /  "<<ngood<<std::endl;
          int trackID=cdsTrackMan->GoodTrackID(it1);
          CDSTrack *track1 = cdsTrackMan->Track(trackID);
          Tools::H1("CDCreducedChi2",track1->Chi(),1000,0,100);
          pCDS* cds=TrackTools::CalcSingleAll(pbeam,track1,cdsMan,true);
          if(!cds) continue;
          cds->SetTrackID(trackID);
          Tools::H1("Mass2"       ,cds->mass(),  500, -0.5, 1.5 );
          // mass2 should be added in pCDS
          double mom=track1->Momentum();
          Tools::H2("overbeta_mom",1./cds->beta(),cds->mom(),  300, 0, 30, 200,-2.,2. );	
          Tools::H2("CDC_mass2_mom", cds->mass2(), mom ,300,-1,11, 300,-1.5,1.5);
          Tools::H2("CDC_mass2_mom_zoom", cds->mass2(), mom ,300,-0.1,0.5, 200,-1.,1.);
          double beta_calc,calc_mass2,tofvtxcdc,tof;
          int seg;
          track1->GetCDHHit(cdsMan,seg,tof);

          TrackTools::FindMass2(track1,bpctrack,tof,beam_mom,Beam_Kaon,beta_calc,calc_mass2,tofvtxcdc);
          Tools::H2("coverbeta_mom" ,1./beta_calc ,mom  ,200,0,10,300,-1.5,1.5);	
          Tools::H2("CDC_cmass_mom" , sqrt(calc_mass2), mom ,400, 0, 4, 300,-1.5,1.5);
          Tools::H2("CDC_cmass2_mom", calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
          Tools::H2("CDC_cmass2_mom_zoom", calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          Tools::H2(Form("CDC_cmass2_mom_id%d",reacID), calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
          Tools::H2(Form("CDC_cmass2_mom_zoom_id%d",reacID), calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          if(CDHparent[seg]==0){
            Tools::H2("CDC_cmass2_mom_primary", calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
            Tools::H2("CDC_cmass2_mom_zoom_primary", calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          }
          if(CDHin[seg]){
            Tools::H2("CDC_cmass2_mom_in", calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
            Tools::H2("CDC_cmass2_mom_zoom_in", calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          }
          if(CDHpdg[seg]==-321){
            Tools::H2("CDC_cmass2_mom_kaon", calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
            Tools::H2("CDC_cmass2_mom_zoom_kaon", calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          }
          if(CDHpdg[seg]==211){
            Tools::H2("CDC_cmass2_mom_zoom_pip", calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
            Tools::H2(Form("CDC_cmass2_mom_zoom_pip_id%d",reacID), calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          }
          if(CDHpdg[seg]==-211){
            Tools::H2("CDC_cmass2_mom_zoom_pim", calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
            Tools::H2(Form("CDC_cmass2_mom_zoom_pim_id%d",reacID), calc_mass2, mom ,300,-0.1,0.5, 200,-1.,1.);
          }
          if(CDHpdg[seg]==2212){
            Tools::H2("CDC_cmass2_mom_proton", calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
            Tools::H2(Form("CDC_cmass2_mom_proton_id%d",reacID), calc_mass2, mom ,300,-1,11, 300,-1.5,1.5);
          }
          //	particle->AddCDS(*cds);
          delete cds;
        }//1track      
        delete pbeam;
      }
      fev++;
      tev+=cdsTrackMan->nGoodTrack();
      //######clear###########
      cdsMan->Clear();
      blMan->Clear(); 
      cdsTrackMan->Clear();
      blTrackMan->Clear(); 
      simMan->Clear();

    } // for(int iev = 0; iev<nev; iev++ ){

    std::cout<<"filled events = "<<fev<<std::endl;
    std::cout<<"tracks = "<<tev<<std::endl;

    of->Write();
    of->Close();
//=======
  std::cout<<" Read "<<simMan-> nEvent()<<" Event read"<<std::endl;
  for( int iev=0; iev<evn; iev++ ){
    std::cout<<" Event Number : "<<iev<<std::endl;
    simMan-> Get(iev, cdsMan, blMan);
    simMan-> PrintEvent();
    simMan-> PrintHit();

    std::cout<<" Please input any char"<<std::endl;
    std::cout<<" q: quit"<<std::endl;
    char in;
    std::cin>>in;
    if( in=='q' ) break;

    cdsMan-> Clear();
    blMan-> Clear();
    simMan-> Clear();
  }
 
  of->Write();
  of->Close();
//>>>>>>> 1.17

    return 0;
    }
