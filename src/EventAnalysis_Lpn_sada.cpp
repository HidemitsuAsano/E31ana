#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include <TLorentzVector.h>
#include "TrackTools.h"
#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "Tools.h"
#include "ELossTools.h"


class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TFile *rtFile2;
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *evTree;

  Float_t tlogprob1;
  Float_t tlogprob2;

  Float_t twoimprob1;
  Float_t twoimprob2;


  Float_t tvtxp_x;
  Float_t tvtxp_y;
  Float_t tvtxp_z;
  Float_t tvtxp_dis;

  Float_t tvtxL_x;
  Float_t tvtxL_y;
  Float_t tvtxL_z;
  Float_t tvtxL_dis;
  Float_t tvtxLp_dis;

  Float_t tvtxpip_x;
  Float_t tvtxpip_y;
  Float_t tvtxpip_z;
  Float_t tvtxpip_dis;

  Float_t tvtxR_x;
  Float_t tvtxR_y;
  Float_t tvtxR_z;
  Float_t tvtxR_dis;

  TLorentzVector tL_L;
  TLorentzVector tL_p;
  TLorentzVector tL_n;
  TLorentzVector tL_K;


  Float_t tim_ppi;
  Float_t tim_Lp;
  Float_t tmm_Lp;
  Float_t tim_Ln;
  Float_t tcosL,tcosp,tcosn;
  Float_t tcosOALp,tcosOALn,tcosOApn;
  Float_t tcosLlab,tcosplab,tcosnlab;
  Float_t tmomL,tmomp,tmomn;
  Float_t tmomLlab,tmomplab,tmomnlab;
  Float_t tmompLlab,tmompiLlab;
  Float_t tcospLlab,tcospiLlab;
  Float_t tcosLkpp,tcospkpp;
  Float_t tTL,tTp,tTn,ttQ;
  Float_t tqtransn;
  Float_t tqtransnCM;
  Float_t tbeam;
  Float_t tRinFid;


  TTree *LpTree;
  TTree *vtxRTree;

  Float_t tvtx_x;
  Float_t tvtx_y;
  Float_t tvtx_z;
  Float_t tvtxb_x;
  Float_t tvtxb_y;
  Float_t tvtxb_z;

  TTree *vtxTree;


  TTree *LKtree;
  Float_t tflagK0;
  Float_t tdca_pipi_K0;
  Float_t tdca_K0k_K0;
  Float_t timpipi_K0;
  Float_t tflagL;
  Float_t tdca_ppi_L;
  Float_t tdca_Lk_L;
  Float_t timppi_L;


  TTree *scaTree;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  CDSTrackingMan *trackMan;
  BeamLineTrackMan *bltrackMan;
  EventHeader *header;
  EventHeader *header2;
  ScalerMan *scaMan;
  int t0,t1;

  int ncheckev;
  int checkrun[26];
  int checkev[26];

  double scainit[40];
  double scaend[40];
  bool INIT;
  bool SCAOVERFLOW[40];
  Time_t Time;
  long t_init;


  int AllGoodTrack;
  int nTrack;
  int CDC_Event_Number;
  std::ofstream ofs;
  
public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

const int MaxTreeSize = 1000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  INIT=true;
  //  spillinit=spillfini=-1;
  for(int i=0;i<40;i++){
    SCAOVERFLOW[i]=false;
    scaend[i]=0;
  }

  //Conf file open
  confMan = conf;
  cdcFile = new TFile(confMan->GetCDSTrackFileName().c_str());
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<confMan->GetCDSTrackFileName().c_str()<< "  !!!"<<std::endl;
    exit(false);
  }

  
  gSystem->Load("libPhysics.so");


  //Getting CDSTracking info. from CDCfile
  cdcTree=(TTree*)cdcFile->Get("EventTree");
  header2=0;
  trackMan=0;
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan );
  cdcTree->SetBranchAddress( "EventHeader" ,&header2);  


  //Making outputfile and tree
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();
  InitializeHistogram();

  std::string outfile2=confMan->GetOutFileName();
  outfile2.insert(outfile2.size()-5,"_tr");
  std::cout<<"tr file "<<outfile2<<std::endl;
  rtFile2 = new TFile( outfile2.c_str(), "recreate" );
  rtFile2->cd();

  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  evTree->Branch( "CDSTrackingMan", &trackMan );

  blMan = new BeamLineHitMan();
  evTree->Branch( "BeamLineHitMan", &blMan );
  bltrackMan = new BeamLineTrackMan();
  evTree->Branch( "BeamLineTrackMan", &bltrackMan );

  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  rtFile->cd();

  t0=clock();
  AllGoodTrack=0;
  nTrack=0;
  CDC_Event_Number=0;





}



void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::USca " << std::endl;
#endif
  rtFile->cd();
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  //  std::cout<<Block_Event_Number<<"\t";
  for( int i=0; i<11; i++ )
    std::cout<<sca[i]<<"\t";
  std::cout<<std::endl;
#endif
  //  std::cout<<"nsca: "<<nsca<<std::endl;
  TH1F *h1;
  if(INIT&&sca[0]<5){
    t_init=Time;    
    for( int i=0; i<scaMan->nsca(); i++ ){
      scainit[i] = scaMan->sca(i)->val();

      std::cout<<scainit[i]<<"\t";

    }
    INIT=false;
  }
  //  std::cout<<std::endl;
  if(INIT){
    scaMan->Clear();
    return;
  }
  for( int i=0; i<scaMan->nsca(); i++ ){
    //  double val = scaMan->sca(i)->val()-scaend[i];
    TString name = scaMan->sca(i)->name();
    //    Tools::Fill2D(name,Block_Event_Number,val);
    if(scaend[i]>9.9e+07&&scaend[i]>scaMan->sca(i)->val()) SCAOVERFLOW[i]=true;
    scaend[i] = scaMan->sca(i)->val();
    if(SCAOVERFLOW[i]) scaend[i]+=1.0e+08;

  }
  h1 = (TH1F*)gFile->Get("Time"); h1->Fill(Block_Event_Number,Time-t_init);
  scaMan->Clear();
}




bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif


  Tools::Fill1D(Form("EventCheck"), 1);//EventC 1 
  Event_Number++;


  //Checking event number consistentcy with CDC tracking file
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) {
      if(CDC_Event_Number>=cdcTree->GetEntries()) return false;
      cdcTree->GetEntry(CDC_Event_Number);

      if(header2->ev()<Event_Number) return true;
      else   if(header2->ev()==Event_Number)
	{
	  CDC_Event_Number++;    	  
	  return true;
	}
    }
    if( status==2 ) {return false;} }


  if( Event_Number%1000==1)
    {
      t1=clock();
      std::cout <<"Run "<<confMan->GetRunNumber()<< " Event# : " << Event_Number <<" "<<CDC_Event_Number<< "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0)*1.0e-6 << std::endl;
    }


  if(CDC_Event_Number>=cdcTree->GetEntries()) return false;
  cdcTree->GetEntry(CDC_Event_Number);

  if(header2->ev()!=Event_Number) return true;



  CDC_Event_Number++;
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  //Converting data=>tree and re-calc. of trackMan of CDS.  
  header->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  trackMan->Calc( cdsMan, confMan,true);  

  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();
  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;


  if( nGoodTrack<1 ){
    header->Clear();
    cdsMan->Clear();
    blMan->Clear();
    bltrackMan->Clear();
    return true;
  }



  //###beam track####
  bltrackMan->DoTracking(blMan,confMan,true,true);


  for(int i=0;i<bltrackMan->ntrackBLC1();i++)
    {
      LocalTrack *blc1=bltrackMan->trackBLC1(i);
    }
  for(int i=0;i<bltrackMan->ntrackBLC2();i++)
    {
      LocalTrack *blc2=bltrackMan->trackBLC2(i);
    }
  for(int i=0;i<bltrackMan->ntrackBPC();i++)
    {
      LocalTrack *bpc=bltrackMan->trackBPC(i);
    }


  //###T0 ##BHD##

  int numT0=0;
  int segT0=0;
  double euT0=0;
  double edT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    HodoscopeLikeHit *hit=blMan->T0(i);
    if( hit->CheckRange() )
      {
	numT0++;
	segT0=hit->seg();
	euT0=hit->eu();
	edT0=hit->ed();
      }
  }
  
  //select T0=1hit  
  if( numT0!=1  /*|| segT0!=3*/) 
    { 
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;      
    }


  //Beam PID by TOF of T0-BHD
  TVector3 vtxT0;
  double ctmT0=0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      ctmT0 = blMan->T0(i)->ctmean();
      confMan->GetGeomMapManager()->GetGPos(CID_T0, blMan->T0(i)->seg(),vtxT0);
    }
  }
	
  int nBHD=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ) nBHD++;
  }
  //       if( nBHD!=1 ) continue;
  double ctmBHD;
  int pid_beam=3;//0:pi 1:K 3:else      
  for( int i=0; i<blMan->nBHD(); i++ ){
    if( blMan->BHD(i)->CheckRange() ){
      ctmBHD = blMan->BHD(i)->ctmean();
      double tofBHDT0 = ctmT0-ctmBHD;
      Tools::Fill1D(Form("tofT0BHD"),tofBHDT0);      

      if(header->kaon()&&tofBHDT0>28.176&&tofBHDT0<29.400)  pid_beam=Beam_Kaon; 
      else if(header->pion()&& tofBHDT0<28.176&&tofBHDT0>24.0) pid_beam=Beam_Pion;
    }
  }
  if(pid_beam==3) 
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }



  //BPC selection
  int nbpc=0;
  int bpcid=-1;
  double chibpc=999;
  for(int i=0;i<bltrackMan->ntrackBPC();i++)
    {

      if(bltrackMan->trackBPC(i)->CheckRange(-30,100)){
	nbpc++;
	bpcid=i;
	chibpc=bltrackMan->trackBPC(i)->chi2all();
      }
    }


  if(nbpc!=1 ) 
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }
  LocalTrack *bpctrack=bltrackMan->trackBPC(bpcid);    

  if(!(bpctrack->CheckRange(-10,10)) || bpctrack->chi2all()>10) 
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }


  for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ){

    trackMan->CalcVertex_beam(trackMan->GoodTrackID(it1),bltrackMan,confMan);

  }

  tflagL=0;tflagK0=0;

  //vector for PIDcon
  std::vector <int> pip_ID;
  std::vector <int> pim_ID;
  std::vector <int> km_ID;
  std::vector <int> p_ID;
  std::vector <int> d_ID;

  std::vector <int> spip_ID;
  std::vector <int> spim_ID;
  std::vector <int> skm_ID;
  std::vector <int> sp_ID;


  std::vector <int> vCDHseg;


  int gamma_CDH=0;
  int nue_CDH=0;
  const double ELV=14.0;// [cm/ns]
  TVector3 tar_cen;tar_cen.SetXYZ(-0.475,0,1.4);
  double time_beam=999;
  bool flagbmom=false,flagany=false,flag_L=false,flag_K0s=false,flag_defK0s=false,flagtar=false;
  TLorentzVector L3_beambf,L3_beam,L3_target,L3_targetP;
  TLorentzVector L3_beambfCM,L3_beamCM,L3_targetCM,L3_targetPCM;
  TLorentzVector L3_dou[3];
  TLorentzVector Lpipdef,Lpimdef;

  TVector3 vtx_react;
  int p_id_L=-1;
  //######beam mom###########

  int nblc1=0;
  int blc1id=-1;
  int nblc2=0;
  int blc2id=-1;

  //Beamline chamber timing selection
  for(int i=0;i<bltrackMan->ntrackBLC1();i++)
    {
      LocalTrack *blc1=bltrackMan->trackBLC1(i);
      if(blc1->CheckRange(-30,100)){
	nblc1++;
	if(blc1->CheckRange(-5,5) &&  bltrackMan->trackBLC1(i)->chi2all()<10)   blc1id=i;
      }	
    }
  for(int i=0;i<bltrackMan->ntrackBLC2();i++)
    {
      LocalTrack *blc2=bltrackMan->trackBLC2(i);
      if(blc2->CheckRange(-30,100)){
	nblc2++;

	if(blc2->CheckRange(-5,5) && bltrackMan->trackBLC2(i)->chi2all()<10) 	blc2id=i;
      }	
    }

  if(nblc1!=1 && blc1id==-1 && nblc2!=1 && blc2id==-1 )
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }


  //BLC2-BPC position consistency
  bool fblc2bpc=false;
  for(int ii=0;ii<bltrackMan->ntrackBLC2();ii++)
    {
      if(ii!=blc2id) continue;
      LocalTrack *blc2=bltrackMan->trackBLC2(ii);
      double xblc2bpc[2],  yblc2bpc[2];	 
      double xmom[2],  ymom[2];	 
      bpctrack->XYPosatZ(-74.85,xblc2bpc[0],yblc2bpc[0]);
      bpctrack->XYPosatZ(-19.7,xmom[0],ymom[0]);
      blc2->XYPosatZ(-74.85,xblc2bpc[1],yblc2bpc[1]);
      blc2->XYPosatZ(-130,xmom[1],ymom[1]);
      double dxdz[2],dydz[2];
      dxdz[0]=(xmom[0]-xblc2bpc[0])/(-19.7+74.85);
      dxdz[1]=(xblc2bpc[1]-xmom[1])/(130-74.85);
      dydz[0]=(ymom[0]-yblc2bpc[0])/(-19.7+74.85);
      dydz[1]=(yblc2bpc[1]-ymom[1])/(130-74.85);

      if((xblc2bpc[1]-xblc2bpc[0])<-0.795 || (xblc2bpc[1]-xblc2bpc[0]) >0.822) fblc2bpc=false;
      else if( (yblc2bpc[1]-yblc2bpc[0])<-0.865|| (yblc2bpc[1]-yblc2bpc[0]) >0.871) fblc2bpc=false;
      else if( (dxdz[1]-dxdz[0])<-0.0240 || (dxdz[1]-dxdz[0]) >0.0250) fblc2bpc=false;//check after!!
      else if((dydz[1]-dydz[0]) <-0.02481 ||(dydz[1]-dydz[0]) >0.02489) fblc2bpc=false;//check after!!
      else fblc2bpc=true;
    }

  if(!fblc2bpc)
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }




  // Calc beam momentum
  double bchi=999;
  for(int ii=0;ii<bltrackMan->ntrackBLC2();ii++)
    {
      for(int i=0;i<bltrackMan->ntrackBLC1();i++)
	{
	  if(nblc2!=1 || ii!=blc2id) continue;
	  if(nblc1!=1 || i!=blc1id) continue;
	  LocalTrack *blc1=bltrackMan->trackBLC1(i);
	  LocalTrack *blc2=bltrackMan->trackBLC2(ii);
	  TVector3 blc2t0= blc2->GetPosatZ(-110);
	  TVector3 bpct0= bpctrack->GetPosatZ(-110);

	  BeamSpectrometer *beamsp=new BeamSpectrometer(confMan);
	  beamsp->TMinuitFit(blc1,blc2,confMan);
	  double beammom=beamsp->mom();
	  double bchitmp=beamsp->chisquare(); 

	  if(pid_beam==Beam_Kaon && bchitmp<bchi)
	    {
	      bchi=bchitmp;
	      double x1,y1,x2,y2;
	      double z1=0,z2=20;
	      bpctrack->XYPosatZ(z1,x1,y1);		  
	      bpctrack->XYPosatZ(z2,x2,y2);
	      TVector3 lp;lp.SetXYZ(x1,y1,z1);
	      TVector3 ls;ls.SetXYZ(x2-x1,y2-y1,z2-z1);ls=ls.Unit();	       
	      TVector3 Pp_beam=beammom*ls; 
	      TVector3 Pp_target;Pp_target.SetXYZ(0,0,0); 
	      //	      std::cout<<"bmom "<<beammom<<std::endl;


	      L3_beambf.SetVectM(Pp_beam , kpMass );
	      L3_target.SetVectM( Pp_target, ThreeHeMass);
	      L3_targetP.SetVectM( Pp_target, pMass);
	      L3_beam=L3_beambf;
	      TVector3 boost=(L3_target+L3_beam).BoostVector();
	      L3_beambfCM=L3_beam;
	      L3_targetCM=L3_target;
	      L3_targetPCM=L3_targetP;
	      L3_beambfCM.Boost(-1*boost);
	      L3_targetCM.Boost(-1*boost);
	      L3_targetPCM.Boost(-1*boost);

	      if(bchi<20)    flagbmom=true;

	    }
	  delete beamsp;
	  
	}
    }

  if(!flagbmom)
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }


  //Vertex study
  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
    double param[5];
    track->GetParameters(param);
    double mom = track->Momentum();
    TVector3 vtxb1, vtxb2,vtxb;
    track->GetVertex(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),vtxb1,vtxb2);

    vtxb=(vtxb1+vtxb2)*0.5;

  }

  //Fiducial cut
  bool   ftarget=false;
  double z=0;
  double x,y;
  bpctrack->XYPosatZ(z,x,y);
  TVector3 vtx_tar;
  vtx_tar.SetXYZ(x,y,z);
  if(GeomTools::GetID(vtx_tar)==CID_Fiducial)
    {
      ftarget=true;
    } 

  if(!ftarget) 
    {
      header->Clear();
      cdsMan->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      return true;
    }




  //PID of CDS##
  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
    double param[5];
    track->GetParameters(param);
    double mom = track->Momentum();
    TVector3 vtxb1, vtxb2,vtxb;

    track->GetVertex(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),vtxb1,vtxb2);
    track->SetPID(-1);
    vtxb=(vtxb1+vtxb2)*0.5;
    if(GeomTools::GetID(vtxb)==CID_Fiducial) 
      {
	flagtar=true; 
	vtx_react=vtxb;
      }

    if(!track->CDHFlag() ) continue;


    double beta=-1.;
    double tof=999.;
    double mass2=-999.;
    int CDHseg=-1;
    for(int icdh=0;icdh<track->nCDHHit();icdh++){
      HodoscopeLikeHit *cdhhit=track->CDHHit(cdsMan,icdh);
      double tmptof    = cdhhit->ctmean()-ctmT0;      

      if(tmptof<tof||tof==999.){
	tof=tmptof;
	CDHseg=cdhhit->seg();
      }
    }
    bool CDHflag=true;
    for(int m=0;m<(int)vCDHseg.size();m++)
      {
	if(CDHseg== vCDHseg[m] ) {CDHflag=false;}
      }
    if(!CDHflag) {

      continue;}  
    vCDHseg.push_back(CDHseg);
    double tmptof,beta_calc;
    TrackTools::CalcBetaMass2(vtxb, bpctrack, track, confMan, pid_beam, (double)L3_beam.Vect().Mag(),tof,beta_calc,mass2);

    //Retiming of CDC track by CDH info.
    track->Retiming(cdsMan,confMan,beta_calc,true);
    for(int m=0;m<5;m++) track->HelixFitting(cdsMan);
    track->Calc(confMan);

    //finalize PID
    if(  !TrackTools::FindMass2C(track,bpctrack,tof,L3_beam.Vect().Mag(),pid_beam,beta_calc,mass2,tmptof) ) continue;

    double mass=sqrt(fabs(mass2));

    int pid=TrackTools::PIDcorr(mom,mass2);      
    track->SetPID(pid);



    //Calc energy loss
    double tmpl;
    TVector3 vtx_beam,vtx_cds;
    if(    !track->CalcVertexTimeLength(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),track->Mass(),vtx_beam,vtx_cds,tmptof,tmpl,true) ) std::cout<<"Eloss false!"<<std::endl;
    //    track->CalcELoss();
    if(pid==CDS_PiMinus)
      {
	pim_ID.push_back(trackMan->GoodTrackID(it));
      }
    else if(pid==CDS_PiPlus)
      {
	pip_ID.push_back(trackMan->GoodTrackID(it));
      }
    else if(pid==CDS_Proton)
      {
	p_ID.push_back(trackMan->GoodTrackID(it));
      }
    else if(pid==CDS_Deuteron)
      {
	d_ID.push_back(trackMan->GoodTrackID(it));
      }
    else if(pid==CDS_Kaon)
      {
	km_ID.push_back(trackMan->GoodTrackID(it));
      }
   

  }// cdc itrack
  //#end of PID


  //nCDC
  int nCDCin=0,nCDC=0;
  for(int layer=1;layer<=15;layer++)
    {
      if(layer<=3) nCDCin+=cdsMan->nCDC(layer);
      nCDC+=cdsMan->nCDC(layer);

    }




  //###short PID####
  // for( int it1=0; it1<trackMan->nShortTrack(); it1++ ){

  //   CDSTrack *track1 = trackMan->ShortTrack(it1);
  //   //    track1->Calc(confMan);
  //   if(!track1->IHFlag()) continue;
  //   HodoscopeLikeHit *ih=0;
  //   for(int l=0;l<track1->nIHHit();l++) ih=track1->IHHit(cdsMan,l);
  //   if(!(ih->seg()>=1 && ih->seg()<=24) ) continue;
  //   TVector3 pos;
  //   //    double len,wid,th,lv;
  //   //    confMan->GetGeomMapManager()->GetParam(CID_IH, ih->seg(),pos,dir,len,wid,th,lv);
  //   confMan->GetGeomMapManager()->GetGPos(CID_IH, ih->seg(),pos);

  //   double phi=pos.Phi();
  //   TVector3 momih;
  //   track1->GetMomentum(track1->IHVertex(),momih);
  //   double mom=track1->Momentum();		    
  //   double deih=ih->eu();
  //   double cospt=cos(phi-momih.Phi());
  //   double deih2=deih*sin(momih.Theta())*cospt;
  //   Tools::Fill2D("IHPID",deih,mom);
  //   Tools::Fill2D("IHPID2",deih2,mom);
  //   //mom=0.7*sqrt(1/de)*exp(-0.5*de)
  //   double thr=0.7*sqrt(1/deih)*exp(-0.5*deih);
  //   double mass;
  //   if(fabs(mom)<thr ) 
  //     {
  // 	mass=piMass;
  // 	if(mom<0)
  // 	  {spim_ID.push_back(it1);
  // 	    track1->SetPID(1);}
  // 	else 
  // 	  {  
  // 	    spip_ID.push_back(it1);
  // 	    track1->SetPID(2);}
  //     }
  //   else  if(fabs(mom)>=thr ) 	 
  //     {
  // 	    mass=pMass;
  // 	    if(mom<0)
  // 	      { 
  // 		track1->SetPID(6);
  // 	      }
  // 	    else
  // 	      {
  // 		sp_ID.push_back(it1);
  // 		track1->SetPID(5);}
  //     }


  // }
  //##end short PID


          

  bool chargedhit=false;
  for(int i=0;i<blMan->nTOF();i++)
    {
      HodoscopeLikeHit *hit=blMan->TOF(i);
      if(hit->CheckRange() )
	{
	  chargedhit=true;
	}
    }
  for(int i=0;i<blMan->nBVC();i++)
    {
      HodoscopeLikeHit *hit=blMan->BVC(i);
      if(hit->CheckRange() )
	{
	  chargedhit=true;
	}
    }

  int nIH=0;
  for(int i=0;i<cdsMan->nIH();i++)
    {
      HodoscopeLikeHit *hit=cdsMan->IH(i);
      if(hit->tdcu()>0 && hit->eu()>0.2)
	{
	  nIH++;
	}
    }

    


  //Lp n_miss event study         
  
  if(flagbmom&&p_ID.size()==2&& pim_ID.size()==1 &&trackMan->nGoodTrack()==3 && !chargedhit  )
    {
      int npL=-1;
      double dca_ppi[2]={999,999};
      double dca_p[2]={999,999};
      double dca_Lp[2]={999,999};
      double dca_L[2]={999,999};
      double dca_reac[2]={999,999};
      double dis_Lct[2]={0,0};
      double dca_sum[2]={999,999};
      bool LinFid[2]={false};      
      bool PinFid[2]={false};
      bool RinFid[2]={false};
      TVector3 vtx_piL[2],vtx_pL[2],vtx_L[2],vtx_p[2];
      TVector3 vtxreac[2];
      TLorentzVector L_pL,L_piL,L_p,L_n,L_L;
      double pt_pL[2],pt_piL[2],pt_p[2];

      double im_ppi[2],mm_ppipi[2];
  
      //calc dca_ppi & vertex
      for(int ip=0;ip<(int)p_ID.size();ip++)
	{
	  for(int npi=0;npi<(int)pim_ID.size();npi++)
	    {
	  
	      CDSTrack *track_p=trackMan->Track(p_ID[ip]);
	      CDSTrack *track_pim=trackMan->Track(pim_ID[npi]);
	      int p2=-1;
	      if(ip==0) p2=1;
	      else if(ip==1) p2=0;
	      CDSTrack *track_p2=trackMan->Track(p_ID[p2]);


	      double dis_vtx_ppi=0;
	      TrackVertex vertex;
	      TVector3 vtx_pi, vtx_pd, vtx;
	      if( !TrackTools::Calc2HelixVertex(track_p,track_pim,vtx_pd,vtx_pi) )
		{std::cout<<"cvtx!! "<<p_ID[ip]<<" "<<pim_ID[npi]<<std::endl;continue;}

	      vtx=(vtx_pd+vtx_pi)*0.5;
	      dis_vtx_ppi=(vtx_pd-vtx_pi).Mag();
	      dca_ppi[ip]=dis_vtx_ppi;
	      vtx_piL[ip]=vtx_pi;
	      vtx_pL[ip]=vtx_pd;

	      pt_pL[ip]=track_p->Pt();
	      pt_piL[ip]=track_pim->Pt();
	      pt_p[ip]=track_p2->Pt();


	      TVector3 Pp_p,Pp_pim,Pp_ppi;
	      if( !track_p->GetMomentum(vtx_pL[ip],Pp_p,true,true) )
		{std::cout<<"cP p!!"<<std::endl;continue;}
	      if( !track_pim->GetMomentum(vtx_piL[ip],Pp_pim,true,true) ) 
		{std::cout<<"cP pim lpn!!"<<std::endl;continue;}
	      TLorentzVector tL_p1,tL_p2,tL_pi;	      
	      tL_p1.SetVectM(Pp_p,pMass);
	      tL_pi.SetVectM(Pp_pim,piMass);

	      im_ppi[ip]=(tL_p1+tL_pi).M();

	      Pp_ppi=Pp_p+Pp_pim;
	      double x1,y1,x2,y2;
	      double z1=0,z2=20;
	      bpctrack->XYPosatZ(z1,x1,y1);		  
	      bpctrack->XYPosatZ(z2,x2,y2);
	      TVector3 lp;lp.SetXYZ(x1,y1,z1);
	      TVector3 ls;ls.SetXYZ(x2-x1,y2-y1,z2-z1);ls=ls.Unit();
	      
	      double dist_L,dltmp=0;TVector3 vtx_Lt,nest;
	      MathTools::LineToLine( vtx,Pp_ppi.Unit(),lp, ls,dltmp,dist_L,vtx_Lt,nest );

	      dca_L[ip]=dist_L;
	      vtx_L[ip]=vtx_Lt;
	      dis_Lct[ip]=(vtx_Lt-vtx).Mag();
	      
	      if(GeomTools::GetID(vtx_Lt)==CID_Fiducial)
		{
		  LinFid[ip]=true;
		}
	     
	      TVector3 vtxb1,vtxb2;
	      track_p2->GetVertex(bpctrack->GetPosatZ(0),bpctrack->GetMomDir(),vtxb1,vtx_p[ip]);
	      dca_p[ip]=(vtx_p[ip]-vtxb1).Mag();
	      dca_Lp[ip]=(vtx_p[ip]-vtx_L[ip]).Mag();
	      TVector3 Pp_p2;


	      vtxreac[ip]=0.5*((vtxb1)+(nest) );
	      dca_reac[ip]=((vtxb1)-(nest)).Mag();

	      if(GeomTools::GetID(vtx_p[ip])==CID_Fiducial)
		{
		  PinFid[ip]=true;
		}

	      if(GeomTools::GetID(vtxreac[ip])==CID_Fiducial)
		{
		  RinFid[ip]=true;
		}


	      if( !track_p2->GetMomentum(vtx_p[ip],Pp_p2,true,true) ) 
		{std::cout<<"cP p lpn!!"<<std::endl;continue;}
	      tL_p2.SetVectM(Pp_p2,pMass);



	      //Eloss of Beam mom
	      double beamtof,momout;
	      double z_pos=-120;
	      
	      ELossTools::CalcElossBeamTGeo(bpctrack->GetPosatZ(z_pos), vtxreac[ip],L3_beambf.Vect().Mag() , kpMass, momout,beamtof);
	      
	      L3_beam.SetVectM( momout*L3_beambf.Vect().Unit(),kpMass);
	      
	      TVector3 boost=(L3_target+L3_beam).BoostVector();
	      L3_beamCM=L3_beam;
	      L3_targetCM=L3_target;
	      L3_targetPCM=L3_targetP;
	      L3_beamCM.Boost(-1*boost);
	      L3_targetCM.Boost(-1*boost);
	      L3_targetPCM.Boost(-1*boost);

	      mm_ppipi[ip]=(L3_target+L3_beam-tL_p1-tL_pi-tL_p2).M();
		
	    }

	}


      //Vertex check
      TF1 *f_dca=new TF1("f_dca","[0]*exp([1]*x)+[2]*exp([3]*(x-[4])**[5])");
      TF1 *f_imL=new TF1("f_imL","gaus(0)");
      f_imL->SetParameters(1,1.1158,0.0019);
      double dcaparam_L[5][6]={

//#hx mcdca_ppi_Lpn
      {0.0111829,  -0.386551,  1.86103,  -1.65383,  -0.488821,  1.50069},  
//#hx mcdca_L_Lpn
      {1.38649,  -6.8509,  -31683.6,  -25.5,  -3.52,  21.9},  
//#hx mcdca_p_Lpn
      {1.43316,  -8.27214,  -31500.1,  -25.5,  -5.1463,  21.9},  
//#hx mcdca_Lp_Lpn
      {2.79166,  -2.44547,  -7.82249,  -6.82579,  -0.225766,  1.29841},  
//#hx mcdca_sum_Lpn
      {5.82406,  -1.5735,  -8.50142,  -0.542443,  -0.89754,  2.91009},  
      };



      double probsum_sin[2]={99,99};
      double probsum_imsin[2]={99,99};
      double prob_L[2][5]={};
      //      double prob_p[2][5]={};
      double prob_ML[2]={};


      //Vertex check
      for(int ip=0;ip<2;ip++)//for cut
	{
	  dca_sum[ip]=dca_ppi[ip]+dca_L[ip]+dca_p[ip]+dca_Lp[ip];

	  prob_ML[ip]=f_imL->Eval(im_ppi[ip]);
	  for(int id=0;id<5;id++)
	    {
	      f_dca->SetParameters(dcaparam_L[id][0],dcaparam_L[id][1],dcaparam_L[id][2],dcaparam_L[id][3],dcaparam_L[id][4],dcaparam_L[id][5]);
	      double dcatmp=0;
	      if(id==0) dcatmp=dca_ppi[ip];
	      else if(id==1) dcatmp=dca_L[ip];
	      else if(id==2) dcatmp=dca_p[ip];
	      else if(id==3) dcatmp=dca_Lp[ip];
	      else if(id==4) dcatmp=dca_sum[ip];

	      prob_L[ip][id]=f_dca->Eval(dcatmp);
	    }

	  probsum_sin[ip]=1;
	  for(int id=0;id<4;id++) probsum_sin[ip]*=prob_L[ip][id];

	  probsum_imsin[ip]=1;
	  for(int id=0;id<4;id++) probsum_imsin[ip]*=prob_L[ip][id];
	  probsum_imsin[ip]*=prob_ML[ip];

	  if(isnan(probsum_sin[ip]) ) probsum_sin[ip]=0;
	  if(isnan(probsum_imsin[ip]) ) probsum_imsin[ip]=0;

	}//ip

      double logprob[2];
      if(probsum_imsin[0]<1e-50) logprob[0]=50-1;
      else   logprob[0]=-1*log10(probsum_imsin[0]);

      if(probsum_imsin[1]<1e-50) logprob[1]=50-1;
      else   logprob[1]=-1*log10(probsum_imsin[1]);

      double woimlogprob[2];
      if(probsum_sin[0]<1e-50) woimlogprob[0]=50-1;
      else   woimlogprob[0]=-1*log10(probsum_sin[0]);

      if(probsum_sin[1]<1e-50) woimlogprob[1]=50-1;
      else   woimlogprob[1]=-1*log10(probsum_sin[1]);

              
      if(probsum_imsin[0] > probsum_imsin[1]) npL=0;
      else npL=1;


      if(logprob[npL]<2.5)
	{	  
	  vtxRTree->Fill();
	}
      //Asano memo 
      //bug?
      if(!RinFid[npL]) npL=-1;


      if(RinFid[npL])
	{
	  tlogprob1=logprob[npL];  tlogprob2=logprob[1-npL];
	  twoimprob1=woimlogprob[npL];  twoimprob2=woimlogprob[1-npL];
	}


      bool flagbprob=true;
      if(logprob[0]>2.5 && logprob[1]>2.5  ) 
	{
	  flagbprob=false;
	}

      if( !RinFid[npL] )  flagbprob=false;

      if( (logprob[0]<5 || logprob[1]<5) && npL!=-1) 
	{
	  rtFile2->cd();
	  evTree->Fill();
	  rtFile->cd();
	}

      //after defining  pi-p pair of  Lambda.
      for(int n=0;n<2;n++){
	if( n!=npL) {continue;}
	
	  CDSTrack *track_pL=trackMan->Track(p_ID[npL]);
	  CDSTrack *track_p=trackMan->Track(p_ID[1-npL]);
	  CDSTrack *track_piL=trackMan->Track(pim_ID[0]);
	  
	  TVector3 P_pL,P_piL,P_p;
	  if( !track_pL->GetMomentum(vtx_pL[npL],P_pL,true,true) )
	    {std::cout<<"cP pL lpn!!"<<std::endl;continue;}
	  if( !track_piL->GetMomentum(vtx_piL[npL],P_piL,true,true) ) 
	    {std::cout<<"cP piL lpn!!"<<std::endl;continue;}
	  if( !track_p->GetMomentum(vtx_p[npL],P_p,true,true) ) 
	    {std::cout<<"cP p lpn!!"<<std::endl;continue;}

	  vtx_react=(vtx_p[npL]+vtx_L[npL])*0.5;	  
	  L_pL.SetVectM(P_pL,pMass);
	  L_piL.SetVectM(P_piL,piMass);
	  L_p.SetVectM(P_p,pMass);

	  //Eloss of Beam mom
	  double beamtof,momout;
	  double z_pos=-120;
	  
	  ELossTools::CalcElossBeamTGeo(bpctrack->GetPosatZ(z_pos), vtx_react,L3_beambf.Vect().Mag() , kpMass, momout,beamtof);
	  
	  L3_beam.SetVectM( momout*L3_beambf.Vect().Unit(),kpMass);
	  tL_K=L3_beam;
	  TVector3 boost=(L3_target+L3_beam).BoostVector();
	  L3_beamCM=L3_beam;
	  L3_targetCM=L3_target;
	  L3_targetPCM=L3_targetP;
	  L3_beamCM.Boost(-1*boost);
	  L3_targetCM.Boost(-1*boost);
	  L3_targetPCM.Boost(-1*boost);

	  double mm_lp=(L3_target+L3_beam-L_pL-L_piL-L_p).M();

	  TVector3 P_missn=(L3_target+L3_beam-L_pL-L_piL-L_p).Vect();
	  L_L=(L_pL+L_piL);
	  double im_ppi_t=(L_pL+L_piL).M();



	  TLorentzVector L_Llab,L_plab,L_nlab;
	  L_Llab=L_L;
	  L_plab=L_p;
	  tL_L=L_L;
	  tL_p=L_p;

	  //	  std::cout<<"tL_L "<<tL_L.M()<<std::endl;

	  TVector3 boost1=(L3_target+L3_beam).BoostVector();	    
	  L_L.Boost(-1*boost1);
	  L_p.Boost(-1*boost1);
	  double cosOAlp=L_L.Vect().Dot(L_p.Vect())/L_L.Vect().Mag()/L_p.Vect().Mag();
	  
	  double cosL=L_L.Vect().Dot(L3_beamCM.Vect())/L_L.Vect().Mag()/L3_beamCM.Vect().Mag();
	  double cosp=L_p.Vect().Dot(L3_beamCM.Vect())/L_p.Vect().Mag()/L3_beamCM.Vect().Mag();

	  double coslabL=L_Llab.Vect().Dot(L3_beam.Vect())/L_Llab.Vect().Mag()/L3_beam.Vect().Mag();
	  double coslabp=L_plab.Vect().Dot(L3_beam.Vect())/L_plab.Vect().Mag()/L3_beam.Vect().Mag();

	  double coslabpL=L_pL.Vect().Dot(L3_beam.Vect())/L_pL.Vect().Mag()/L3_beam.Vect().Mag();
	  double coslabpiL=L_piL.Vect().Dot(L3_beam.Vect())/L_piL.Vect().Mag()/L3_beam.Vect().Mag();



	  L_n.SetVectM(P_missn,nMass);
	  L_nlab=L_n;
	  tL_n=L_n;
	  
		  //daliz
	  L_n.Boost(-1*boost1);
	  
	  double Tp=L_p.E()-pMass;
	  double TL=L_L.E()-lMass;
	  double Tn=L_n.E()-nMass;
	  double tQ=(L3_beamCM+L3_targetCM).M()-lMass-pMass-nMass;
	  
	  
	  double cosOAln=L_L.Vect().Dot(L_n.Vect())/L_L.Vect().Mag()/L_n.Vect().Mag();
	  double cosOApn=L_p.Vect().Dot(L_n.Vect())/L_p.Vect().Mag()/L_n.Vect().Mag();
	  double cosn=L_n.Vect().Dot(L3_beamCM.Vect())/L_n.Vect().Mag()/L3_beamCM.Vect().Mag();
	  
	  double coslabn=L_nlab.Vect().Dot(L3_beam.Vect())/L_nlab.Vect().Mag()/L3_beam.Vect().Mag();
	  
	  
	  //	  double qtransn=pow(L3_beam.Vect().Mag(),2)+pow(L_nlab.Vect().Mag(),2)-2*L3_beam.Vect().Mag()*L_nlab.Vect().Mag()*coslabn;
	  double qtransn=(L_nlab.Vect()- L3_beam.Vect()).Mag();
	  double qtransnCM=(L_n.Vect()- L3_beamCM.Vect()).Mag();
	  

	  //Kpp Flame
	  TVector3 boost_kpp=(L_Llab+L_plab).BoostVector();	    
	  TLorentzVector L_Lkpp=L_Llab;
	  TLorentzVector L_pkpp=L_plab;
	  TLorentzVector L_kpp=(L_Llab+L_plab);
	  L_Lkpp.Boost(-1*boost_kpp);
	  L_pkpp.Boost(-1*boost_kpp);
	  
	  
	  double cosLkpp=L_Lkpp.Vect().Dot(L_kpp.Vect())/L_Lkpp.Vect().Mag()/L_kpp.Vect().Mag();
	  double cospkpp=L_pkpp.Vect().Dot(L_kpp.Vect())/L_pkpp.Vect().Mag()/L_kpp.Vect().Mag();
	  
	  


	  if(flagbprob)
	    {

	      Tools::Fill1D(Form("logprobim_Lpn"),-log10(prob_ML[n]) );	      
	      Tools::Fill2D(Form("nIM_Mm_Lp"),(L_L+L_p).M(),mm_lp);	      
	      Tools::Fill1D("ncosL_Lp", cosL);
	      Tools::Fill1D("ncosp_Lp", cosp);	      
	      Tools::Fill1D("ncoslabL_Lp", coslabL);
	      Tools::Fill1D("ncoslabp_Lp", coslabp);	      
	      Tools::Fill1D("nmomL_Lp", L_L.Vect().Mag());
	      Tools::Fill1D("nmomp_Lp", L_p.Vect().Mag());
	      Tools::Fill1D("ncosOA_Lp", cosOAlp);


	      if(mm_lp>0.84 &&mm_lp<1.04)//select missn
		{
	      
		  

		  Tools::Fill1D(Form("logprobim_sum_Lpn"),logprob[npL] );
		  
		  Tools::Fill1D("nIM_Lp_Lpn",(L_L+L_p).M());		  
		  Tools::Fill1D("nIM_Ln_Lpn",(L_L+L_n).M());		  
		  Tools::Fill1D("nOAlp_Lpn", cosOAlp);
		  Tools::Fill1D("nOAln_Lpn", cosOAln);
		  Tools::Fill1D("nOApn_Lpn", cosOApn);
		  Tools::Fill1D("ncosL_Lpn", cosL);
		  Tools::Fill1D("ncosp_Lpn", cosp);
		  Tools::Fill1D("ncosn_Lpn", cosn);
		  
		  Tools::Fill1D("ncoslabL_Lpn", coslabL);
		  Tools::Fill1D("ncoslabp_Lpn", coslabp);
		  Tools::Fill1D("ncoslabn_Lpn", coslabn);
	      
		  Tools::Fill1D(Form("ncosLkpp"),cosLkpp );		  
		  Tools::Fill1D(Form("ncospkpp"),cospkpp );		  
		  
		  Tools::Fill1D("nmomL_Lpn", L_L.Vect().Mag());
		  Tools::Fill1D("nmomp_Lpn", L_p.Vect().Mag());
		  Tools::Fill1D("nmomn_Lpn", L_n.Vect().Mag());		  
		  

		  Tools::Fill2D("nIM_Mm_Lpn", (L_L+L_p).M(),mm_lp);
		  Tools::Fill2D("nIM_OAlp_Lpn", (L_L+L_p).M(),cosOAlp);
		  Tools::Fill2D("nIM_OAln_Lpn", (L_L+L_p).M(),cosOAln);
		  Tools::Fill2D("nIM_OApn_Lpn", (L_L+L_p).M(),cosOApn);
		  Tools::Fill2D("nIM_cosL_Lpn", (L_L+L_p).M(),cosL);
		  Tools::Fill2D("nIM_cosp_Lpn", (L_L+L_p).M(),cosp);
		  Tools::Fill2D("nIM_cosn_Lpn", (L_L+L_p).M(),cosn);
		  Tools::Fill2D("nIM_coslabL_Lpn", (L_L+L_p).M(),coslabL);
		  Tools::Fill2D("nIM_coslabp_Lpn", (L_L+L_p).M(),coslabp);
		  Tools::Fill2D("nIM_coslabn_Lpn", (L_L+L_p).M(),coslabn);
		  
		  Tools::Fill2D("nIM_qtransn_Lpn", (L_L+L_p).M(),qtransn);
	      
		  Tools::Fill2D("nIM_momL_Lpn", (L_L+L_p).M(),L_L.Vect().Mag() );
		  Tools::Fill2D("nIM_momp_Lpn", (L_L+L_p).M(),L_p.Vect().Mag());
		  Tools::Fill2D("nIM_momn_Lpn", (L_L+L_p).M(),L_n.Vect().Mag());
		  




		  Tools::Fill1D("nmomL_Lpn_lab", L_Llab.Vect().Mag());
		  Tools::Fill1D("nmomp_Lpn_lab", L_plab.Vect().Mag());
		  Tools::Fill1D("nmomn_Lpn_lab", L_nlab.Vect().Mag());		  
		  
		  Tools::Fill2D("ndalitz_Lpn", (Tp-Tn)/sqrt(3)/tQ,TL/tQ );
		  Tools::Fill2D("ndalitz2_Lpn", (Tp-TL)/sqrt(3)/tQ,Tn/tQ );
		  Tools::Fill2D("ndalitzIM_Lpn", pow((L_L+L_p).M(),2),pow((L_L+L_n).M(),2) );
		  Tools::Fill2D("nIM_Tn_Lpn",(L_L+L_p).M(),Tn);
	     		  
		}	//missing neutron  
	    }
      }// for nLp    
 
    }//Lpn 
  


    


  header->Clear();
  cdsMan->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  rtFile->cd();
  for(int i=0;i<40;i++){
    TH1F* h1 = (TH1F*)gFile->Get("Scaler"); 
    std::cout<<i<<"  "<<scaend[i]<<std::endl;
    h1->Fill(i,scaend[i]-scainit[i]);
  }


  //  confMan->SaveCDSParam();
  gFile->Write();
  gFile->Close();
  cdcFile->Close();

  delete cdsMan;
  delete header;
}



void EventAnalysis::InitializeHistogram()
{

  //  new TTree("Lptree","tree for Lpn event");

  Tools::newTH1F(Form("EventCheck"),20,0,20);// 

  new TH1F("Time","Time",3000,-0.5,2999.5);
  new TH1F("Scaler","Scaler",41,-0.5,40.5);



  Tools::newTH1F(Form("tofT0BHD"),1100,-5,50);


  //Lpn new

  Tools::newTH2F(Form("nIM_Mm_Lp"),1000,2,3,1400,0.4,1.8);		  
  Tools::newTH1F(Form("ncosL_Lp"),800, -1,1);
  Tools::newTH1F(Form("ncosp_Lp"),800, -1,1);
  
  Tools::newTH1F(Form("ncoslabL_Lp"),800, -1,1);
  Tools::newTH1F(Form("ncoslabp_Lp"),800, -1,1);
  
  Tools::newTH1F(Form("ncosLkpp"),800, -1,1);
  Tools::newTH1F(Form("ncospkpp"),800, -1,1);

  Tools::newTH1F(Form("nmomL_Lp"),700, 0,1.4);
  Tools::newTH1F(Form("nmomp_Lp"),700, 0,1.4);
  Tools::newTH1F(Form("ncosOA_Lp"),800, -1,1);
  
  Tools::newTH1F(Form("nIM_Lp_Lpn"),1000,2,3);		  
  Tools::newTH1F(Form("nIM_Ln_Lpn"),1000,2,3);		  
  Tools::newTH1F(Form("nOAlp_Lpn"),800, -1,1);
  Tools::newTH1F(Form("nOAln_Lpn"),800, -1,1);
  Tools::newTH1F(Form("nOApn_Lpn"),800, -1,1);
  Tools::newTH1F(Form("ncosL_Lpn"),800, -1,1);
  Tools::newTH1F(Form("ncosp_Lpn"),800, -1,1);
  Tools::newTH1F(Form("ncosn_Lpn"),800, -1,1);
  Tools::newTH1F(Form("ncoslabL_Lpn"),800, -1,1);
  Tools::newTH1F(Form("ncoslabp_Lpn"),800, -1,1);
  Tools::newTH1F(Form("ncoslabn_Lpn"),800, -1,1);
  
  Tools::newTH1F(Form("nmomL_Lpn"),700, 0,1.4);
  Tools::newTH1F(Form("nmomp_Lpn"),700, 0,1.4);
  Tools::newTH1F(Form("nmomn_Lpn"),700, 0,1.4);
  
  Tools::newTH2F(Form("nIM_Mm_Lpn"),1000,2,3,300, 0.8,1.1);
  Tools::newTH2F(Form("nIM_OAlp_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_OAln_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_OApn_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_cosL_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_cosp_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_cosn_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_coslabL_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_coslabp_Lpn"),1000,2,3,800, -1,1);
  Tools::newTH2F(Form("nIM_coslabn_Lpn"),1000,2,3,800, -1,1);
  
  Tools::newTH2F(Form("nIM_qtransn_Lpn"),1000,2,3,800, -0.3,1.3);
  
  Tools::newTH2F(Form("nIM_momL_Lpn"),1000,2,3,700, 0,1.4);
  Tools::newTH2F(Form("nIM_momp_Lpn"),1000,2,3,700, 0,1.4);
  Tools::newTH2F(Form("nIM_momn_Lpn"),1000,2,3,700, 0,1.4);
  
  
  Tools::newTH1F(Form("nmomL_Lpn_lab"),700, 0,1.4);
  Tools::newTH1F(Form("nmomp_Lpn_lab"),700, 0,1.4);
  Tools::newTH1F(Form("nmomn_Lpn_lab"),700, 0,1.4);
  
  Tools::newTH2F(Form("ndalitz_Lpn"),480, -0.6,0.6,400,0,1 );
  Tools::newTH2F(Form("ndalitz2_Lpn"), 480, -0.6,0.6,400,0,1 );
  Tools::newTH2F(Form("ndalitz_n_Lpn"),480, -0.6,0.6,400,0,1 );
  Tools::newTH2F(Form("ndalitz2_n_Lpn"), 480, -0.6,0.6,400,0,1 );
  Tools::newTH2F(Form("nIM_Tn_Lpn"),1200,0,3,400,0,1);
  
    



}
EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
