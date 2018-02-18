#include <iostream>
#include <string>
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
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TLorentzVector.h>

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
//#include "Display3D.h"

//#######Need knucl lib file!! Set your dir!! ###
#include "/home/yamaga/geant/knucl4/include/KnuclRootData.h"
//#######################################
#define SIM 0

int main( int argc, char **argv )
{
  std::cout << " argc:" << argc << std::endl;
  if(argc != 4 ){
    std::cout << "Plese set Conffile Outputfile Inputfile  "<< std::endl;
    std::cout << argv[0] << " <conf file> <out file> <in file> " << std::endl;
    return 0;
  }
  for(int i=0; i<argc; i++ ){
    std::cout << "  " << argv[i] << std::endl;
  }

  std::string confFile,outFile,inFile;
  if( argv[1][0] != '!' ) confFile = argv[1];
  if( argv[2][0] != '!' ) outFile = argv[2];
  if( argv[3][0] != '!' ) inFile = argv[3];

  int runnum=99999;
  if(argc>1) runnum=atoi(argv[1]);


  gSystem->Load("libPhysics.so");
  gSystem->Load( "./lib/libAll.so" );
  //#######Need knucl lib file!! Set your dir!! ###
  //gSystem->Load("~/geant/knucl4/build/libKnuclRootData.so");
  gSystem->Load("~/geant/knucl4/build/KnuclRootData_cc.so");
  //#######################################

  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();

  TFile *f;
  TTree *evtree;
  f = new TFile( inFile.c_str());
  evtree = (TTree*)f->Get( "EventTree" );

  std::cout<<"histogram initailization start"<<std::endl;
  TFile *fout;
  fout = new TFile( outFile.c_str(), "recreate" );
  fout->cd();
  {  
    new TH1F( "nTrack", "nTrack", 200, 0, 100 );
    new TH1F( "chi", "chi", 300, 0, 30 );
    new TH1F( "mom", "mom", 300, 0, 3. );
    new TH1F( "costheta", "costheta", 200, -1., 1. );
    new TH2F( "pt_gene_ana", "pt_gene_ana", 200, 0, 2., 200, 0, 2. );
    new TH2F( "pt_diff", "pt_diff", 200, 0, 2., 300, -0.15, 0.15 );
    new TH2F( "mom_gene_ana", "mom", 300, 0, 3. , 300, 0, 3. );
    new TH2F( "mom_diff", "mom_diff", 300, 0, 3. , 300, -0.15, 0.15 );
    new TH2F( "costheta_gene_ana","costheta", 200, -1., 1. , 200, -1., 1. );
    new TH2F( "costheta_diff","costheta_diff", 200, -1., 1. , 200, -0.1, 0.1 );
    new TH2F( "angle_diff","angle_diff", 200, -1., 1. , 200, -0.1, 0.1 );
    new TH2F( "vertex_z_sim_ana","vertex_z_sim_ana", 300, -150,150. , 300, -150., 150 );
    new TH2F( "vertex_z_diff","vertex_z_diff", 300, -150,150. , 400, -20., 20 );
    new TH2F( "vertex_xy","vertex_xy", 400, -20.,20. , 400, -20., 20 );
    new TH2F( "vertex_xy_sim","vertex_xy_sim", 400, -20.,20. , 400, -20., 20 );
  }

  std::cout<<"histogram initailization end"<<std::endl;
  TH1F *h1;

  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  BeamLineTrackMan *blTrackMan = 0;
  CDSTrackingMan *cdsTrackMan = 0;
  EventHeader *head = 0;
  MCData *mcdata=0;

  std::cout << "###" << std::endl;
  
  evtree->SetBranchAddress( "CDSHitMan",      &cdsMan );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan  );
  evtree->SetBranchAddress( "CDSTrackingMan", &cdsTrackMan );
  evtree->SetBranchAddress( "BeamLineTrackMan", &blTrackMan  );
  evtree->SetBranchAddress( "MCData",         &mcdata );

  CDSTrackingMan *trackMan = new CDSTrackingMan();

  int nev = evtree->GetEntries();
  std::cout << " AllEvent : " << nev << std::endl;
  std::cout << "|0%                  |50%                |100% "     
    <<   std::endl;
  int moniter=0;  
  for( int iev=0; iev<nev; iev++ ){
    if( iev==0 )
      std::cout << "|";
    if( iev%(nev/40)==0 )
      std::cout << "*"<<std::flush;
    if( iev%(100)==0 )
    {
      if(moniter==0) std::cout << "\\\b"<<std::flush;
      else if(moniter==1) std::cout << "-\b"<<std::flush;
      else if(moniter==2) std::cout << "/\b"<<std::flush;
      else if(moniter==3) std::cout << "|\b"<<std::flush;
      moniter++;
      if(moniter==4) moniter=0;
    }
    if( iev==nev-1 )
      std::cout << "| fin!" << std::endl;

    int status = conf->CheckEvNum( iev, 1 );
    if(status){ 
      gFile->Write();
      gFile->Close();
      if( status==1 ) return true;
      if( status==2 ) return false;
    }

    evtree->GetEvent(iev);

    // for(int i=1; i<=NumOfCDCLayers; i++){
    //   for(int j=0; j<cdsMan->nCDC(i); j++){
    // 	CDCHit *cdc = cdsMan->CDC(i, j);
    // 	std::cerr<<"CDC : "<<cdc->hid()<<" "<<cdc->layer()<<" "<<cdc->wire()<<" "<<cdc->dl()<<std::endl;
    //   }
    // }
    int ntr=cdsTrackMan->nGoodTrack();
    h1= (TH1F*)gFile->Get("nTrack"); h1->Fill(ntr);    
    for( int it1=0; it1<cdsTrackMan->nGoodTrack(); it1++ ){      
      CDSTrack *track1 = cdsTrackMan->Track(cdsTrackMan->GoodTrackID(it1));
      //      track1->HelixFitting(cdsMan);
      track1->Calc(conf);
      h1= (TH1F*)gFile->Get("chi"); h1->Fill(track1->Chi());
      h1= (TH1F*)gFile->Get("mom"); h1->Fill(track1->Momentum());
      //      std::cout<<track1->Momentum()<<std::endl;
      double par[5];
      track1->GetParameters(par);
      // for(int i=0;i<5;i++)
      // 	std::cout<<par[i]<<"\t";
      // std::cout<<std::endl;
      //      std::cout<<par[4]<<"\t"<<TMath::ATan(par[4])<<std::endl;
      double cos=-TMath::Cos(TMath::ATan(par[4])+TMath::Pi()/2.);
      h1= (TH1F*)gFile->Get("costheta"); h1->Fill( cos );

      TVector3 vertexhana;
      TVector3 vertexlana;
      double dis;
      MathTools::LineToHelix(TVector3(0.,0.,0.),TVector3(0.,0.,1.),par,vertexlana,vertexhana,dis);
      vertexhana *= 10;
      for(int itr=0;itr<mcdata->trackSize();itr++){
        int pdgid=mcdata->track(itr)->pdgID();
        int pid=mcdata->track(itr)->parentTrackID();
        if(pid==0&&pdgid==-211){
          TH2F* h2;
          TVector3 mom=mcdata->track(itr)->momentum()*0.001;	  
          TVector3 vertexsim=mcdata->track(itr)->vertex();
          h2=(TH2F*)gFile->Get("mom_gene_ana");
          h2->Fill(mom.Mag(),track1->Momentum());
          h2=(TH2F*)gFile->Get("mom_diff");
          h2->Fill(mom.Mag(),track1->Momentum()-mom.Mag());
          h2=(TH2F*)gFile->Get("pt_gene_ana");
          h2->Fill(mom.Pt(),track1->Pt());
          h2=(TH2F*)gFile->Get("pt_diff");
          h2->Fill(mom.Pt(),track1->Pt()-mom.Pt());
          h2=(TH2F*)gFile->Get("costheta_gene_ana");
          h2->Fill(mom.CosTheta(),cos);
          h2=(TH2F*)gFile->Get("costheta_diff");
          h2->Fill(mom.CosTheta(),cos-mom.CosTheta());
          h2=(TH2F*)gFile->Get("angle_diff");
          h2->Fill(mom.CosTheta(),TMath::ACos(cos)-TMath::ACos(mom.CosTheta()));

          h2=(TH2F*)gFile->Get("vertex_xy");
          h2->Fill(vertexhana.X(),vertexhana.Y());
          h2=(TH2F*)gFile->Get("vertex_xy_sim");
          h2->Fill(vertexsim.X(),vertexsim.Y());
          h2=(TH2F*)gFile->Get("vertex_z_sim_ana");
          h2->Fill(vertexsim.Z(),vertexhana.Z());
          h2=(TH2F*)gFile->Get("vertex_z_diff");
          h2->Fill(vertexsim.Z(),vertexhana.Z()-vertexsim.Z());
        }
      }      
    }//1track      
  } // iev

  gFile->Write();
  gFile->Close();

  return 0;
}

