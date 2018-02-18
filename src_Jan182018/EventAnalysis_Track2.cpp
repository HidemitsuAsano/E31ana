#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  int counter[5];
public:
  void Initialize( ConfMan *conf );
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

const int MaxTreeSize = 1900000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  trackMan = new CDSTrackingMan();  
  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSTrackingMan", &trackMan );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );

  rtFile->cd();
  new TH1F( "chi", "chi", 2000, 0, 1000 );
  new TH1F( "ctCDH", "ctCDH", 500, -10, 10 );
  new TH1F( "emeanCDH", "emeanCDH", 500, -1, 49 );
  new TH1F( "ctT0", "ctT0", 500, -10, 10 );
  new TH1F( "TOF", "TOF", 500, -50, 50 );
  new TH1F( "Rho", "Rho", 500, -500, 500 );
  new TH1F( "Mom", "Mom", 500, -1.5, 1.5 );
  new TH2F( "PID1", "PID1", 300, -11, 19, 300, -1.5, 1.5 );
  new TH2F( "PID2", "PID2", 300, -11, 19, 500, -1, 49 );

  new TH2F( "chinum", "chinu", 1000, 0, 10000, 1000, 0, 5000 );

  counter[0]=0;
  counter[1]=0;
  counter[2]=0;
  counter[3]=0;
  counter[4]=0;

}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
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
  scaTree->Fill();
  scaMan->Clear();
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%1000==0 )
    {
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
    std::cout << " counter[0]# : " << counter[0] 
	      << " [1]# : " << counter[1] 
	      << " [2]# : " << counter[2] 
	      << " [3]# : " << counter[3] 
	      << " [4]# : " << counter[4] 
	      << std::endl;
    }
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  int nT0 = 0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
    }
  }
  //std::cout << " nT0:" << nT0 << std::endl;
  //if( nT0!=1 ) return true;
  //std::cout << " OK" << std::endl;

  /***********************/
  /**** Tracking **********/
  /***********************/

  /*                */
  /* event analysis */
  /*                */

  //######## Time Offset ##############// 
  double toffset=0;
  cdsMan->Calc(confMan,toffset);
  //###################################//


  trackMan->Clear();
  //  int counter[4]={0};
  int numCDC;
  for(int tracking=0;tracking<1;tracking++){


    //###### event cut######
    int numCDH=0;
    for( int i=0; i<cdsMan->nCDH(); i++ )
      {
	if( 0<cdsMan->CDH(i)->tdcu() && cdsMan->CDH(i)->tdcu()<4000
	    && 0<cdsMan->CDH(i)->tdcd() && cdsMan->CDH(i)->tdcd()<4000)
	  numCDH++;
      }

    numCDC=0;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	int tdc = cdsMan->CDC(layer,i)->tdc();
	if( 400<tdc && tdc<900){
	  numCDC++;
// 	  std::cout << "2 numCDC :" << numCDC << std::endl;
// 	  std::cout << " layer:" << layer << " tdc:" << tdc << std::endl;
	}
      }      
    }
    
    if(  !(10<numCDC&&numCDC<70) ){
      //continue;
      header->Clear();
      blMan->Clear();
      cdsMan->Clear();
      trackMan->Clear();
      return true;
    }


    //    trackMan->SetKilledLayer(8);
    trackMan->SearchCluster(cdsMan);
    trackMan->FindingTrack(cdsMan);

    //std::cout<<"tracknum= "<<trackMan->nTrack()<<std::endl;

    if(trackMan->nTrack()>4) {continue;}
      
    counter[0]+=trackMan->nTrack();
    //####### each fit###############
    
    double arho[2],ax_c[2],ay_c[2],aPt[2];
    double aparam[2][5];
    double param[5];
    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);

	double a=track->A();
	double b=track->B();
	double c=track->C();
 
	double x_mid=track->MidX();
	double y_mid=track->MidY();

	double ap=b/a;
	double bp=-1;
	double cp=y_mid-ap*x_mid;

	double jitter=track->Jitter();
	//	double theta =track->Theta();
	if(fabs(jitter)<1e-100 ) continue;
	double rho=fabs( (30*30+jitter*jitter)/(2*jitter) );
	//if(theta>180 || theta<0 ) rho=-1*rho;

	double xp=rho/sqrt(1+b*b);

	double x_c;
	if(jitter>0 && y_mid>0 )       x_c=x_mid-xp; 
	else if(jitter<0 && y_mid>0 )  x_c=x_mid+xp; 
	else if(jitter>0 && y_mid<0 )  x_c=x_mid-xp; 
	else if(jitter<0 && y_mid<0 )  x_c=x_mid+xp; 
	double y_c=b*x_c+cp;

	double d_c=fabs( sqrt(x_c*x_c+y_c*y_c) );
	double cos_c=-x_c/d_c;
	double sin_c=-y_c/d_c;
	
	double x_o=x_c+rho*cos_c;
	double y_o=y_c+rho*sin_c;
	param[0]=sqrt(x_o*x_o+y_o*y_o );
	//double x_o=-a*c/(b*b+a*a);
	//double y_o=-b*c/(b*b+a*a);
	//param[0]=fabs(c)/sqrt(a*a+b*b);
	
	rho=fabs(rho);
	if(rho<d_c) param[0]=param[0];
	else if(rho>d_c) param[0]=-param[0];
       
	if(cos_c>0 && sin_c>0)   param[1]=atan(sin_c/cos_c);
	else if(cos_c<0 && sin_c>0)   param[1]=TMath::Pi()+atan(sin_c/cos_c);
	else if(cos_c<0 && sin_c<0)   param[1]=TMath::Pi()+atan(sin_c/cos_c);
	else if(cos_c>0 && sin_c<0)   param[1]=2*TMath::Pi()+atan(sin_c/cos_c);

	if(jitter>0 && y_mid>0 ) 
	  {param[2]=-1./rho;param[0]=-1*param[0];param[1]=param[1];}
	else if(jitter<0 && y_mid>0 ) 
	  {param[2]=1./rho;param[0]=param[0];param[1]-=TMath::Pi();}
	else if(jitter>0 && y_mid<0 ) 
	  {param[2]=1./rho;param[0]=param[0];param[1]-=TMath::Pi();}
	else if(jitter<0 && y_mid<0 ) 
	  {param[2]=-1./rho;param[0]=-1*param[0];param[1]=param[1];}
	param[3]=0;
	param[4]=0;

	counter[1]++;
	track->SetParameters(param);
       	if(!track->FirstTrackFitting() ) continue;
	//track->SecondTrackFitting();
	bool flag=true;
	for(int trial=0;trial<10;trial++) { if(!track->TrackFitting() ) flag=false;}
	if(!flag) continue;
	
	counter[2]++;
	double dis=999;
	int CDHnum=-1;
	for( int ii=0; ii<cdsMan->nCDH(); ii++ )
	  {
	    track->GetParameters(param);
	    double arho=1./fabs(param[2]);
	    double ax_c=(param[0]+1./param[2] )*cos(param[1]);
	    double ay_c=(param[0]+1./param[2] )*sin(param[1]);
	    double cdh_x=cdsMan->CDH(ii)->x();
	    double cdh_y=cdsMan->CDH(ii)->y();
	    if(y_mid*cdh_y<0 ) continue;
	    double distmp,xest,yest;
	    track->PointToCircle(cdh_x,cdh_y,arho,ax_c,ay_c,distmp,xest,yest);
	    if(dis==999){CDHnum=ii;dis=distmp; }
	    else if(distmp<dis){CDHnum=ii;dis=distmp; }
	   
	  }
	if( 0<cdsMan->CDH(CDHnum)->tdcu() && cdsMan->CDH(CDHnum)->tdcu()<4000
	    && 0<cdsMan->CDH(CDHnum)->tdcd() && cdsMan->CDH(CDHnum)->tdcd()<4000)
	  track->SetCDHHit(*cdsMan->CDH(CDHnum) );

	
	if(! track->FirstHelixFitting() ) continue;
	if( track->Chi()>80 ) continue;

	std::cout<<"first "<<track->Chi();

	if(! track->SecondHelixFitting() ) continue;
	std::cout<<" second "<<track->Chi();
	//	if(Event_Number%10==0) std::cout<<" second "<<track->Chi();
	counter[3]++;
	track->SetHitPos();
		
	for(int trial=0;trial<5;trial++) { if(!track->HelixFitting() ) {flag=false;break;}}
	if(!flag) continue;
	//counter[4]++;

	
	track->DeleteExtraHit();
	
	//	if(! track->FinalHelixFitting() ) continue;
	///	std::cout<<" final "<<track->Chi()<<std::endl;

	int numt=0;
	double Chi=999,Chitmp=999;
	bool hflag1=true;	
	while(  numt<20 )

	  {
	    Chi=Chitmp;
	    for(int trial=0;trial<10;trial++)
	      { if(!track->HelixFitting() ){hflag1=false; break;}  }
	    //    std::cout<<"Chi2/dof= "<<track->Chi()<<" / "<<track->Dof()<<std::endl;
	    Chitmp=track->Chi();
	    if( (Chi-Chitmp)<0.1 ) break;
	    numt++;
	  }
	if(!hflag1) {std::cout<<"skip event"<<std::endl;continue;}
	//	if( Chitmp>2000 ) {std::cout<<"skip2 event"<<std::endl;continue;}
	counter[4]++;
	std::cout<<" helix "<<track->Chi()<<std::endl;
	//	if(Event_Number%10==0)  std::cout<<" helix "<<track->Chi()<<std::endl;
	/*
	track->GetParameters(aparam[i]);
	arho[i]=1./aparam[i][2];
	ax_c[i]=(aparam[i][0]+1./aparam[i][2])*cos(aparam[i][1]);
	ay_c[i]=(aparam[i][0]+1./aparam[i][2])*sin(aparam[i][1]);
	//    std::cout<<"rho= "<<rho<<" x_c= "<<x_c<<" y_c= "<<y_c<<std::endl;
	aPt[i]=0.3*0.5*arho[i]/100.;
	//	std::cout<<"Pt"<<i<<"= "<<aPt[i]<<"GeV/c"<<std::endl;
	*/
#if 0
	std::cout<<"aparam"<<i<<" ";
    for(int n=0;n<5;n++) std::cout<<aparam[i][n]<<" ";
    std::cout<<std::endl;
#endif 
        
    }
    
    }       
    /***********************/
    /**********************/

  rtFile->cd();
  TH1F *h1;
  TH2F *h2;

  //std::cout << " nTrack:" << trackMan->nTrack() << std::endl;

  double nGoodTrack=0;
  for( int i=0; i<trackMan->nTrack(); i++ ){
    double chi = trackMan->Track(i)->Chi();
    int numsecond = trackMan->Track(i)->GetSecondComb();
    //std::cout << " i:" << i << " chi:" << chi << std::endl;
    h1 = (TH1F*)rtFile->Get("chi"); h1->Fill( chi );
    h2 = (TH2F*)rtFile->Get("chinum"); h2->Fill( chi, numsecond );
    if( chi<10 ) nGoodTrack++;
  }
  if( nGoodTrack<1 ){
    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }
  
  rtFile->cd();
  for( int i=0; i<trackMan->nTrack(); i++ ){
    CDSTrack *track = trackMan->Track(i);
    double param[5];
    track->GetParameters(param);
    double rho = 1./param[2];
    double p = 0.003*0.5*rho;
    HodoscopeLikeHit *cdh = track->CDHHit();
    double ctCDH = cdh->ctmean();
    double emeanCDH = cdh->emean();
  rtFile->cd();
    h1 = (TH1F*)rtFile->Get("ctCDH"); h1->Fill( ctCDH );
    h1 = (TH1F*)rtFile->Get("emeanCDH"); h1->Fill( emeanCDH );
//     std::cout << " track:" << track << " rho:" << rho << " ctCDH:" << ctCDH << " emeanCDH:" << emeanCDH
// 	      << " p:" << p
// 	      << std::endl;

    for( int j=0; j<blMan->nT0(); j++ ){
      if( blMan->T0(j)->CheckRange() ){
	double ctT0 = blMan->T0(j)->ctmean();
	double tofT0CDH = ctCDH-ctT0;
  rtFile->cd();
// 	std::cout << " j:" << j << " ctT0:" << ctT0 << " tofT0CDH:" << tofT0CDH << std::endl;
	h1 = (TH1F*)rtFile->Get("ctT0"); h1->Fill( ctT0 );
	h1 = (TH1F*)rtFile->Get("TOF"); h1->Fill( tofT0CDH );
	h1 = (TH1F*)rtFile->Get("Rho"); h1->Fill( rho );
	h1 = (TH1F*)rtFile->Get("Mom"); h1->Fill( p );
	h2 = (TH2F*)rtFile->Get("PID1"); h2->Fill( tofT0CDH, p );
	h2 = (TH2F*)rtFile->Get("PID2"); h2->Fill( tofT0CDH, emeanCDH );
      }
    }
    

   
  }


  evTree->Fill();

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete trackMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
