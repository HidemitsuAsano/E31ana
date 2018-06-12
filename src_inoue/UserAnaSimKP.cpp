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
#include "HodoClusterMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "Tools.h"
#include "TrackTools.h"
#include "Particle.h"
#include "TDatabasePDG.h"

//#######Need knucl lib file!! Set your dir!! ###
#include "KnuclRootData.h"
//#######################################
#define SIM 0
const char* PCName[5]={"Kaon","Pion","Proton","Deuteron","Gamma"};
const char* cdspar[]={"Piplus","Proton","Deuteron","Triton","Helium3","Piminus","Kaon","Lambda","K0"};

void FillAll(Particle *particle, TLorentzVector L_cds, TVector3 vtx, const char* cname);

bool PROTONARM;
bool NEUTRAL;
bool genOnNC;
bool genOnNC2;
TVector3 vtx_gen;
TLorentzVector L_mm3hen_gen;

int main( int argc, char **argv )
{
  std::cout << " argc:" << argc << std::endl;
  if(argc != 5 ){
    std::cout << "Plese set Conffile Outputfile Inputfile  "<< std::endl;
    return 0;
  }
  
  for(int i=0; i<argc; i++ ){
    std::cout << "  " << argv[i] << std::endl;
  }

  std::string confFile,outFile,inFile,inFile2;
  if( argv[1][0] != '!' ) confFile = argv[1];
  if( argv[2][0] != '!' ) outFile = argv[2];
  if( argv[3][0] != '!' ) inFile = argv[3];
  if( argv[4][0] != '!' ) inFile2 = argv[4];
 
  int runnum=99999;
  if(argc>1) runnum=atoi(argv[1]);
  
  gSystem->Load("libPhysics.so");
  gSystem->Load( "./lib/libAll.so" );
  //#######Need knucl lib file!! Set your dir!! ###
  //  gSystem->Load("~/work/geant/knucl3/libKnuclRootData.so");
  //#######################################

  ConfMan *conf = new ConfMan(confFile);
  conf->Initialize();

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");

  TFile *f;
  TTree *evtree;
  f = new TFile( inFile.c_str());
  evtree = (TTree*)f->Get( "EventTree" );
  double missncenter=0.94;

  TH1F *h1;
  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  BeamLineTrackMan *blTrackMan = new BeamLineTrackMan();
  CDSTrackingMan *cdsTrackMan = 0;
  EventHeader *head = 0;

  MCData *mcData=new MCData();
  ReactionData *reacData=0;
  DetectorData *detectorData=0;

  TFile *f2 = new TFile( inFile2.c_str() );
  TTree *evtree2 = (TTree*)f2->Get("tree");

  evtree->SetBranchAddress( "CDSHitMan",      &cdsMan );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan  );
  evtree->SetBranchAddress( "CDSTrackingMan", &cdsTrackMan );

  evtree2->SetBranchAddress( "ReactionData", &reacData  );
  evtree2->SetBranchAddress( "MCData",         &mcData );
  evtree2->SetBranchAddress( "DetectorData", &detectorData );

  std::cout<<"histogram initailization start"<<std::endl;
  TFile *fout;
  fout = new TFile( outFile.c_str(), "recreate" );

  Particle *particle   = new Particle();
  int nev = evtree->GetEntries();
  //  if(nev>10000) nev=10000;
  std::cout << " AllEvent : " << nev << std::endl;
  std::cout << "|0%                  |50%                |100% "     
	    <<   std::endl;
  int moniter=0;  
  int FermiFlag = 0;
  TString initname[2];
  TString name[6];

  for( int iev=0; iev<nev; iev++ ){
    particle->Clear();
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
      if( status==1 ) continue;
      gFile->Write();
      gFile->Close();
      if( status==2 ) return false;
    }    
    evtree->GetEvent(iev);
    evtree2->GetEvent(iev);

    // Check T0 hit information
    bool T0OK=false;
    for(int i=0;i<blMan->nT0();i++)
      if(blMan->T0(i)->ctmean()<0.1) T0OK=true;
    if(!T0OK) continue;

    //
    // Check wether the primary neutron passes through NC or not
    //
    genOnNC=false;
    for( int i=0; i<detectorData->detectorHitSize(); i++ ){
      DetectorHit *detectorHit = detectorData->detectorHit(i);
      int detectorID = detectorHit->detectorID();
      int parentID = detectorHit->parentID();
      int pdgCode  = detectorHit->pdg();      
      if(detectorID==CID_NC){
	if(pdgCode==2112&&parentID==0)
	  {
	    genOnNC=true; 
	    break;
	  }
      }
    }

    //-----------generated data analysis
    //p1: kaon, p2: proton
    double mass[2]={kpMass,pMass};
    int pid[2]={-321,2212};
    const char* name[2]={"K","FP"};
    TLorentzVector L_gen_p1;
    TLorentzVector L_gen_p2;
    TLorentzVector L_gen_beam;
    double beam_mom;
    double trid_k0s=-1;
    bool KCHARGE=false;
    bool KONCDH=false;
    for(int itr=0;itr<mcData->trackSize();itr++){
      int parid=mcData->track(itr)->parentTrackID();
      int pdgid=mcData->track(itr)->pdgID();
      TVector3 mom=mcData->track(itr)->momentum();
      if(pdgid==pid[0]&&parid==0){
	L_gen_p1.SetVectM(mom*0.001,mass[0]);
	KCHARGE=true;
	double fl=mcData->track(itr)->FlightLength();
	if(fl>0) KONCDH=true;
      }
      if(pdgid==-311&&parid==0){
	L_gen_p1.SetVectM(mom*0.001,k0Mass);

      }
      if(pdgid==pid[1]&&parid==0)  L_gen_p2.SetVectM(mom*0.001,mass[1]);
      if(pdgid==2112&&parid==0)    L_gen_p2.SetVectM(mom*0.001,nMass);
      if(pdgid==321&&parid==0){
	L_gen_beam.SetVectM(-mom*0.001,kpMass);
	beam_mom=mom.Mag()*0.001; 
	vtx_gen=mcData->track(itr)->vertex()*0.1;
      }      
      if(pdgid==130&&parid==0){
	L_gen_p1.SetVectM(mom*0.001,k0Mass);
      }      
    }
    if(GeomTools::GetID(vtx_gen)!=CID_Fiducial) continue;
    
    // detector data analysis 
    blTrackMan->Clear();
    blTrackMan->DoTracking( blMan, conf );
    cdsTrackMan->Calc( cdsMan, conf );  
    int ntr=cdsTrackMan->nGoodTrack();
    Tools::H1("nTrackCDC" ,ntr                     , 10,0,10);    
    Tools::H1("nTrackBPC" ,blTrackMan->ntrackBPC() , 10,0,10);    
    Tools::H1("nTrackFDC1",blTrackMan->ntrackFDC1(), 10,0,10);    
    Tools::H1("nHitIH"    ,cdsMan->nIH()           , 10,0,10);    
    Tools::H1("nHitCDH"   ,cdsMan->nCDH()          , 10,0,10);    
    Tools::H1("nHitNC"    ,blMan->nNC()            , 10,0,10);    
    Tools::H1("nHitPC"    ,blMan->nPC()            , 10,0,10);    
    Tools::H1("nHitCVC"   ,blMan->nCVC()           , 10,0,10);    
    Tools::H1("nHitBVC"   ,blMan->nBVC()           , 10,0,10);    
    int ncdh=0;
    //    std::cout<<"------------------------"<<std::endl;
    for(int i=0;i<cdsMan->nCDH();i++){
      HodoscopeLikeHit *hit=cdsMan->CDH(i);
      double time=hit->ctmean();
      double edep=hit->emean();
      //      std::cout<<i<<"  "<<hit->seg()<<"  "<<time<<"  "<<edep<<std::endl;
      if(time<50&&edep>1) ncdh++;
    }
    Tools::H1("nHitCDHcut"   ,ncdh          , 10,0,10);    
    int nih=0;
    for(int i=0;i<cdsMan->nIH();i++){
      HodoscopeLikeHit *hit=cdsMan->IH(i);
      double time=hit->ctmean();
      double edep=hit->emean();
      if(time<50&&edep>0.2) nih++;
    }
    Tools::H1("nHitIHcut"   ,nih          , 10,0,10);    

    TLorentzVector L_gen_tar=L_gen_p1+L_gen_p2-L_gen_beam;
    TLorentzVector L_target; L_target.SetVectM( TVector3(0,0,0), ThreeHeMass);
    TLorentzVector L_proton; L_proton.SetVectM( TVector3(0,0,0), pMass);
    
    TVector3 p2gendir=L_gen_p2.Vect().Unit();
    TVector3 p2genpos=vtx_gen+p2gendir*((1470-vtx_gen.Z())/p2gendir.Z());
    TVector3 vboost=(L_gen_beam+L_gen_tar).BoostVector();
    TVector3 vboost2=(L_gen_beam+L_proton).BoostVector();
    TVector3 vboost3=(L_gen_beam+L_target).BoostVector();

    TLorentzVector L_gencm_tar=L_gen_tar;
    L_gencm_tar.Boost(-vboost);

    TLorentzVector L_p1_cm=L_gen_p1;
    L_p1_cm.Boost(-vboost);
    TLorentzVector L_p1_cm2=L_gen_p1;
    L_p1_cm2.Boost(-vboost2);
    TLorentzVector L_p1_cm3=L_gen_p1;
    L_p1_cm3.Boost(-vboost3);


    //    if(TMath::Abs(p2genpos.X())<160&&TMath::Abs(p2genpos.Y())<75) genOnNC=true;
    Tools::H1( Form("%sR_gen",name[0]), vtx_gen.Perp() ,100,0,5.);
    Tools::H1( Form("%sCosThetaLab_gen",name[0]), L_gen_p1.Vect().CosTheta() ,100,-1,1.);
    Tools::H1( Form("%sCosThetaCM_gen" ,name[0]), L_p1_cm.Vect().CosTheta() ,100,-1,1.);
    Tools::H1( Form("%sCosThetaCM2_gen",name[0]), L_p1_cm2.Vect().CosTheta() ,100,-1,1.);
    Tools::H1( Form("%sCosThetaCM3_gen",name[0]), L_p1_cm3.Vect().CosTheta() ,100,-1,1.);
    Tools::H1( Form("FermiMomGene")             , L_gen_tar.Vect().Mag() ,100,0,1.);
    Tools::H1( Form("FermiAngleLab")            , L_gen_tar.Vect().CosTheta() ,100,-1,1.);
    Tools::H1( Form("FermiAngleCM")             , L_gencm_tar.Vect().CosTheta() ,100,-1,1.);
    Tools::H2( Form("MissPos%s_gen",name[0]), p2genpos.X(),p2genpos.Y(),100,-500,500,60,-300,300);
    if(blMan->nNC()>0){
      Tools::H2( Form("MissPos%sifNC_gen",name[0]), p2genpos.X(),p2genpos.Y() ,100,-500,500,60,-300,300);
    }

    if(blTrackMan->ntrackBPC()!=1) continue;
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
    pbeam->SetT0seg(3);
    pbeam->SetT0Time(0.);
    pbeam->SetT0Pos(t0pos);
    //    pbeam->SetBHDseg(segBHD);
    //    pbeam->SetBHDT0TOF(tofBHDT0);    
    pbeam->SetBPCPos(bpcpos);
    pbeam->SetBPCDir(bpcdir);
    //    pbeam->SetVertex(vertex);
    pbeam->SetMomentum(beam_mom*(1+gRandom->Gaus(0,0.002)));
    particle->AddBeam(*pbeam);
    delete pbeam;

    //##################################
    // CDS 1 track analysis
    //##################################
    int ngood=cdsTrackMan->nGoodTrack();
    for( int it1=0; it1<ngood; it1++ ){      
      //      std::cout<<it1<<"  /  "<<ngood<<std::endl;
      int trackID=cdsTrackMan->GoodTrackID(it1);
      CDSTrack *track1 = cdsTrackMan->Track(trackID);
      //      track1->HelixFitting(cdsMan);
      //      track1->Calc(conf);      
      Tools::H1("CDCreducedChi2",track1->Chi(),1000,0,100);
      pCDS* cds=TrackTools::CalcSingleAll(particle->beam(0),track1,cdsMan,true);
      if(!cds) continue;

      cds->SetTrackID(trackID);
      Tools::H1("Mass2"       ,cds->mass(),  500, -0.5, 1.5 );
      // mass2 should be added in pCDS
      Tools::H2("overbeta_mom",1./cds->beta(),cds->mom(),  300, 0, 30, 200,-2.,2. );	
      particle->AddCDS(*cds);
      delete cds;
    }//1track      

    //##################################
    // CDS 2 track analysis
    //##################################
    for( int id1=0; id1<particle->nCDS(); id1++ ){
      for( int id2=id1+1; id2<particle->nCDS(); id2++ ){
	pCDS* cds1=particle->cdsi(id1);
	pCDS* cds2=particle->cdsi(id2);
	if(!cds1||!cds2) continue;
	pCDS *pro=TrackTools::Calc2HelixAll(particle->beam(0),cds1,cds2,cdsTrackMan,true);
	if(!pro) continue;
	pro->SetDaughterID1(id1);
	pro->SetDaughterID2(id2);	
	particle->AddProduct(*pro);  
	delete pro;
      }
    }

    //###############################
    //####Proton Arm analysis #######
    //###############################  
    PROTONARM=false;
    pPC* ppc;
    if(blMan->nPC()==1){
      ppc=new pPC();
      HodoscopeLikeHit *hit = blMan->PC(0);
      ppc->SetSegment(hit->seg());
      ppc->SetTime(hit->ctmean());
      TVector3 pos;
      conf->GetGeomMapManager()->GetGPos(CID_PC,hit->seg(),pos);
      ppc->SetHitPosition(pos);
      int ntrfdc1=blTrackMan->ntrackBLDC(CID_FDC1);
      //      std::cout<<ntrfdc1<<std::endl;
      if(ntrfdc1==1&&blMan->nBVC()==1){
	TVector3 gfdcpos,gfdcdir,fdc1pos;
	conf->GetBLDCWireMapManager()->GetGParam(CID_FDC1,gfdcpos,gfdcdir);
	fdc1pos=blTrackMan->trackFDC1(0)->GetPosatZ(gfdcpos.Z());
	ppc->SetFDC1Pos(fdc1pos);
	particle->AddPC(*ppc);
	//	if(blMan->nBVC()>0)
	PROTONARM=true;
	  //	std::cout<<"!!!"<<std::endl;
      }
      delete ppc;
    }
    if(PROTONARM) ppc=particle->pc(0);

    //##################################
    //####Neutron Counter analysis #####
    //##################################  
    NEUTRAL=false;
    HodoClusterMan* clMan = new HodoClusterMan();
    pNC *pnc;
    clMan->Clustering(blMan,conf);
    if(blMan->nCVC()==0&&blMan->nBVC()==0&&clMan->ncluster(CID_NC)>0){
      int seg;
      double time;
      TVector3 pos;     
      if(clMan->ncfirsttime(-1,seg,time,pos)){
	conf->GetGeomMapManager()->GetGPos(CID_NC,seg,pos);
	pnc=new pNC(seg,time,pos);
	NEUTRAL=true;
	particle->AddNC(*pnc);
	pnc->CalcMom(particle->beam(0),TVector3(0,0,0));
	int cid=CID_NC;
	delete pnc;
      }
    }
    if(NEUTRAL) pnc=particle->nc(0);
    genOnNC2=false;
    if(NEUTRAL&&blMan->nBVC()==0){
      if(pnc->mom().Mag()>1.0){
	genOnNC2=true;
	//	break;
      }
    }

    //##################################
    // Fill Histogram
    //##################################
    TLorentzVector L_deuteron; L_deuteron.SetVectM( TVector3(0,0,0), dMass);
    L_mm3hen_gen=L_target+L_gen_beam-L_gen_p2;
    TLorentzVector L_beam=particle->beam(0)->GetLorentzVector();
    Tools::H1("nCDS",particle->nCDS(),10,-0.5,9.5);

    //##################################
    // CDS 1 track
    //##################################
    bool FNEUTRON=false;
    bool FPROTON=false;
    bool MISSN=false;
    bool SPEC=false;
    bool CDSK=false;
    bool TWOCHARGE=false;
    for(int i=0;i<particle->nCDS();i++){
      pCDS* cds=particle->cdsi(i);
      TLorentzVector L_cds=cds->GetLorentzVector();
      int pid=cds->pid();
      TVector3 vtx=cds->vbeam();
      TLorentzVector L_cdscm=L_cds;
      L_cdscm.Boost(-vboost2);

      const char *cname=cdspar[pid];
      if(NEUTRAL)     pnc->CalcMom(particle->beam(0),vtx);
      if(PROTONARM){
	ppc->CalcMom(particle->beam(0),vtx);
	//	std::cout<<L_gen_p2.Vect().Mag()<<"  "<<ppc->mom().Mag()<<std::endl;
      }
      FillAll(particle,L_cds,vtx,cname);      

      if(NEUTRAL&&pnc->isneutron()&&pnc->mom().Mag()>1.0){
	FNEUTRON=true;
	TLorentzVector L_p=pnc->GetLorentzVector();
	TLorentzVector L_x=L_target+L_beam-L_p-L_cds;
	if(L_x.M()<1.95&&L_x.M()>1.85){	    
	  SPEC=true;
	}
      }
      if(PROTONARM&&ppc->pid()==F_Proton&&ppc->mom().Mag()>1.0){
	FPROTON=true;
	TLorentzVector L_p=ppc->GetLorentzVector();
	TLorentzVector L_x=L_target+L_beam-L_p-L_cds;
	if(L_x.M()<1.95&&L_x.M()>1.85){
	  SPEC=true;
	}
      }
      if(cds->pid()==CDS_Kaon&&KCHARGE&&KONCDH) CDSK=true;
    }
    //##################################
    // CDS 2 track
    //##################################

    for(int i=0;i<particle->nProduct();i++){
      pCDS *pro=particle->product(i);
      TString tmpname;
      if(pro->comb()==(int)( pow(2,CDS_PiMinus)+pow(2,CDS_PiPlus) )) tmpname="PiPi";
      else if(pro->comb()==(int)( pow(2,CDS_PiMinus)+pow(2,CDS_Proton) )) tmpname="PiP";
      else continue;
      const char* cname=tmpname.Data();
      TLorentzVector L_cds=pro->GetLorentzVector();	
      TLorentzVector L_cdscm=L_cds;
      L_cdscm.Boost(-vboost2);
      TVector3 vtx=pro->vbeam();
      bool K0FLAG=false;
      if(NEUTRAL&&pnc->isneutron()&&pnc->mom().Mag()>1.0){
	FNEUTRON=true;
	TLorentzVector L_p=pnc->GetLorentzVector();
	TLorentzVector L_x=L_target+L_beam-L_p-L_cds;
	if(L_x.M()<1.95&&L_x.M()>1.85){	    
	  SPEC=true;
	}
      }
      if(PROTONARM&&ppc->pid()==F_Proton&&ppc->mom().Mag()>1.0){
	FPROTON=true;
	TLorentzVector L_p=ppc->GetLorentzVector();
	TLorentzVector L_x=L_target+L_beam-L_p-L_cds;
	if(L_x.M()<1.95&&L_x.M()>1.85){
	  SPEC=true;
	}
      }
      if(!strcmp(cname,"PiPi") && L_cds.M()>k0Mass-0.012 && L_cds.M()<(k0Mass+0.012) ){
	K0FLAG=true;
	CDSK=true;
      }
      if(K0FLAG&&!KCHARGE){
	TVector3 nmom=L_gen_p2.Vect();
	TLorentzVector L_tmpn;
	double missn=(L_target+L_gen_beam-L_gen_p2).M();
	double missk0n=(L_target+L_gen_beam-L_gen_p1-L_gen_p2).M();
	for(int ishift=-10;ishift<10;ishift++){
	  L_tmpn.SetVectM(nmom.Unit()*(nmom.Mag()+0.001*ishift),nMass);
	  double mm3hek0n=(L_target+L_gen_beam-L_tmpn-L_gen_p1).M()-missk0n;
	  double tmpmm=(L_target+L_gen_beam-L_tmpn).M()-missn;
	  mm3hek0n*=1000;
	  tmpmm*=1000;
	  Tools::H2(Form("Mm3hek0n_Mm3Hen_shift%d",ishift),mm3hek0n,tmpmm, 200,-10,10,200,-10,10);
	  Tools::H2(Form("NMom_Mm3Hen_shift%d",ishift),nmom.Mag(),tmpmm,150,0,1.5,200,-10,10);
	  Tools::H2(Form("Mm3Hen_Mm3Hen_shift%d",ishift),missn,tmpmm,200,2.2,2.7,200,-10,10);
	  Tools::H2(Form("NMom_Mm3hek0n_shift%d",ishift),nmom.Mag(),mm3hek0n,150,0,1.5,200,-10,10);
	}

	TLorentzVector L_missn=L_proton+L_beam-L_cds;
	TVector3 diru=L_missn.Vect().Unit();
	TVector3 ncpos=vtx+diru*((1470-vtx.Z())/diru.Z());
	//	ncpos.Print();
	double xxx=L_missn.M();
	//	if(xxx>0.8&&xxx<1.08){
	if(xxx>0.88&&xxx<1.00){
	  MISSN=true;
	  for(int i=0;i<10;i++){
	    if(( TMath::Abs(ncpos.X())<(20+20*i) )&&( TMath::Abs(ncpos.Y())<(10+10*i) )){
	      Tools::H2( Form("MissPosK0_gen_cut%d" ,i),  p2genpos.X(), p2genpos.Y() ,100,-500,500,60,-300,300);
	      if(genOnNC)	   
		Tools::H2( Form("MissPosK0ifNC_gen_cut%d" ,i),  p2genpos.X(), p2genpos.Y() ,100,-200,200,60,-100,100);
	      if(genOnNC2)	   
		Tools::H2( Form("MissPosK0ifNC_gen2_cut%d" ,i),  p2genpos.X(), p2genpos.Y() ,100,-200,200,60,-100,100);
	    }
	  }
	}
	if(particle->nCDS()==2&&nih==2&&ncdh==2) TWOCHARGE=true;
      }
      FillAll(particle,L_cds,vtx,cname);      
    } // 2track

    static const int ncut=9;
    bool CUT[ncut];
    TString cutname[ncut]={"_gen2",
			   "_wN",
			   "_wP",
			   "_CDSK",
			   "_CDSK_1track",
			   "_CDSKwP",
			   "_CDSKwN",
			   "_CDSK_2charge",
			   "_CDSKwNonNC"			   
    };

    CUT[0]=true;
    CUT[1]=FNEUTRON;
    CUT[2]=FPROTON;
    CUT[3]=CDSK&&MISSN;
    CUT[4]=CDSK&&MISSN&&particle->nCDS()==1;
    CUT[5]=CDSK&&MISSN&&FPROTON&&SPEC;
    CUT[6]=CDSK&&MISSN&&FNEUTRON&&SPEC;
    CUT[7]=CDSK&&MISSN&&TWOCHARGE;
    CUT[8]=CDSK&&MISSN&&FNEUTRON&&SPEC&&genOnNC;
    for(int icut=0;icut<ncut;icut++){
      if(!CUT[icut]) continue;
      Tools::H1( Form("%sR%s",name[0],cutname[icut].Data()), vtx_gen.Perp() ,100,0,5.);
      Tools::H1( Form("%sCosThetaLab%s",name[0],cutname[icut].Data()), L_gen_p1.Vect().CosTheta() ,100,-1,1.);
      Tools::H1( Form("%sCosThetaCM%s" ,name[0],cutname[icut].Data()), L_p1_cm.Vect().CosTheta() ,100,-1,1.);
      Tools::H1( Form("%sCosThetaCM2%s",name[0],cutname[icut].Data()), L_p1_cm2.Vect().CosTheta() ,100,-1,1.);
      Tools::H1( Form("%sCosThetaCM3%s",name[0],cutname[icut].Data()), L_p1_cm3.Vect().CosTheta() ,100,-1,1.);
      //      Tools::H1( Form("%sCosThetaCM4%s",name[0],cutname[icut].Data()), L_cdscm.Vect().CosTheta() ,100,-1,1.);
    }
  } // iev
  gFile->Write();
  gFile->Close();  
  return 0;
}

void FillAll(Particle *particle, TLorentzVector L_cds, TVector3 vtx, const char* cname){

  TLorentzVector L_target; L_target.SetVectM( TVector3(0,0,0), ThreeHeMass);
  TLorentzVector L_proton; L_proton.SetVectM( TVector3(0,0,0), pMass);
  TLorentzVector L_deuteron; L_deuteron.SetVectM( TVector3(0,0,0), dMass);
  TLorentzVector L_beam=particle->beam(0)->GetLorentzVector();
  TLorentzVector L_mm3hecds= L_beam + L_target - L_cds;
  TLorentzVector L_mmdcds= L_beam + L_deuteron - L_cds;
  TLorentzVector L_mmpcds= L_beam + L_proton - L_cds;
  double mom=L_cds.Vect().Mag();
  double angle=TMath::Cos(L_cds.Vect().Angle(particle->beam(0)->bpcdir()));
  
  Tools::H1(Form("Inv%s",cname),L_cds.M(),1000,0,2);
  Tools::H1(Form("Mm3He%s",cname),L_mm3hecds.M(),1000,0,4);
  Tools::H1(Form("Mmd%s"  ,cname),L_mmdcds.M()  ,1000,0,3);
  Tools::H1(Form("Mmp%s"  ,cname),L_mmpcds.M()  ,1000,0,2);
  Tools::H2(Form("Mom_CosLab_%s",cname),TMath::Abs(mom),angle,150, 0, 1.5 ,200, -1., 1.);
  
  Tools::H2(Form("Vertex_XY_%s",cname),vtx.X(),vtx.Y(), 200,  -5,  5 ,200, -5, 5);
  Tools::H2(Form("Vertex_ZX_%s",cname),vtx.Z(),vtx.X(), 200, -10, 10 ,200, -5, 5);
  Tools::H2(Form("Vertex_ZY_%s",cname),vtx.Z(),vtx.Y(), 200, -10, 10 ,200, -5, 5);
  Tools::H1(Form("Vertex_Xdiff_%s",cname),(vtx-vtx_gen).X(), 200,  -5,  5);
  Tools::H1(Form("Vertex_Ydiff_%s",cname),(vtx-vtx_gen).Y(), 200,  -5,  5);
  Tools::H1(Form("Vertex_Zdiff_%s",cname),(vtx-vtx_gen).Z(), 200, - 5,  5);

  if(PROTONARM){
    pPC* ppc=particle->pc(0);
    int pcpid=ppc->pid();
    const char* pname=PCName[pcpid];
    TLorentzVector L_p=ppc->GetLorentzVector();
    TLorentzVector L_invcds=L_cds + L_p;
    TLorentzVector L_mm3he=L_beam + L_target - L_p;
    TLorentzVector L_mm3hecds=L_beam + L_target - L_p - L_cds;
    Tools::H2( Form("BetaAnglePC%s"    ,pname), 1./ppc->beta(), ppc->angle(), 200,0.5,4.5,100,0,1.);  
    Tools::H1( Form("Mmp%swPC%s"     ,cname,pname), L_mmpcds.M(),1000,0,2);  
    Tools::H1( Form("Mm3HePC%s"            ,pname), L_mm3he.M() ,2000,1,3);
    Tools::H1( Form("Mm3HePC%s_%stag",pname,cname), L_mm3he.M() ,2000,1,3);
    Tools::H1( Form("Mm3HePC%sgen"            ,pname), L_mm3hen_gen.M() ,2000,1,3);
    Tools::H1( Form("Mm3HePC%sgen_%stag",pname,cname), L_mm3hen_gen.M() ,2000,1,3);
    Tools::H1( Form("Mmp%swPC%s"     ,cname,pname), L_mmpcds.M(),1000,0,2);
    Tools::H2( Form("Mom_CosLab_%swPC%s",cname,pname),
	       TMath::Abs(mom),angle,150, 0, 1.5 ,200, -1., 1.);
    Tools::H1( Form("InvPC%s_%s"        ,pname,cname),
	       L_invcds.M(),2000,0,4);
    Tools::H2( Form("InvMm3HePC%s_%s"   ,pname,cname),
	       L_invcds.M(),L_mm3hecds.M(),200,1,3,200,0.5,1.5);
    Tools::H1( Form("Mm3He%sPC%s"       ,cname,pname),
	       L_mm3hecds.M(),2000,0,4);      
    Tools::H2( Form("Mm3He%sPC%s_Mm3HePC%s",cname,pname,pname),
	       L_mm3hecds.M(),L_mm3he.M(),200,1,3,200,2,3);      
  }

  if(NEUTRAL){
    pNC* pnc=particle->nc(0); 
    if(pnc->isneutron()){
      TLorentzVector L_n=pnc->GetLorentzVector();
      TLorentzVector L_invcds=L_cds + L_n;
      TLorentzVector L_mm3he=L_beam + L_target - L_n;
      TLorentzVector L_mm3hecds=L_beam + L_target - L_n - L_cds;
      Tools::H1( Form("Mmp%swN"        ,cname), L_mmpcds.M() ,1000,0,2);   
      Tools::H1( Form("Mm3HeN"               ), L_mm3he.M()  ,2000,1,3);
      Tools::H1( Form("MomN"                 ), pnc->mom().Mag() ,2000,0,2);
      Tools::H1( Form("TOFN"                 ), pnc->time()   ,3000,0,300);
      Tools::H1( Form("Mm3HeN_%stag"   ,cname), L_mm3he.M()  ,2000,1,3);
      Tools::H1( Form("Mm3HeNgen"               ), L_mm3hen_gen.M()  ,2000,1,3);
      Tools::H1( Form("Mm3HeNgen_%stag"   ,cname), L_mm3hen_gen.M()  ,2000,1,3);
      if(genOnNC){
	Tools::H1( Form("Mm3HeNonNC"               ), L_mm3he.M()  ,2000,1,3);
	Tools::H1( Form("Mm3HeNonNC_%stag"   ,cname), L_mm3he.M()  ,2000,1,3);
	Tools::H1( Form("Mm3HeNgenonNC"               ), L_mm3hen_gen.M()  ,2000,1,3);
	Tools::H1( Form("Mm3HeNgenonNC_%stag"   ,cname), L_mm3hen_gen.M()  ,2000,1,3);
      }
      Tools::H1( Form("Mmp%swN"        ,cname), L_mmpcds.M() ,1000,0,2);
      Tools::H2( Form("Mom_CosLab_%swN",cname),
		 TMath::Abs(mom),angle,150, 0, 1.5 ,200, -1., 1.);
      Tools::H1( Form("InvN_%s"        ,cname), L_invcds.M(),2000,0,4);
      Tools::H2( Form("InvMm3HeN_%s"   ,cname),
		 L_invcds.M(),L_mm3hecds.M(),200,1,3,200,0.5,1.5);
      Tools::H1( Form("Mm3He%sN"       ,cname), L_mm3hecds.M(),2000,0,4);      
      Tools::H2( Form("Mm3He%sN_Mm3HeN",cname), L_mm3hecds.M(),L_mm3he.M(),200,1,3,200,2,3);      	  

      pnc->CalcMom(particle->beam(0),vtx_gen);
      L_n=pnc->GetLorentzVector();
      L_mm3he=L_beam + L_target - L_n;
      Tools::H1( Form("Mm3HeNvtxgen"               ), L_mm3he.M()  ,2000,1,3);
      Tools::H1( Form("MomNvtxgen"                 ), pnc->mom().Mag() ,2000,0,2);
      Tools::H1( Form("TOFNvtxgen"                 ), pnc->time()   ,3000,0,300);
      Tools::H1( Form("Mm3HeNvtxgen_%stag"   ,cname), L_mm3he.M()  ,2000,1,3);
    }
  }
}
