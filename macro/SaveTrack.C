#include <fstream.h>
#include <iostream.h>
#include <vector.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"

#define DISPLAY 0
#define DISPLAY2 0

class ConfMan;
class Display;
class CDSHitMan;
class CDCHit;
class CDSTrackingMan;
class CDSTrack;
class CircleFit;
class HelixFit;
void SaveTrack()
{

  /*** load library ***/
  gSystem->Load("lib/libAll.so");



  /*** assign input file & call tree ***/
  TFile *f = new TFile("/w/e15/sada/tree/1736.root");
  TTree *evtree = (TTree*)f->Get("CDSTree");


  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/analyzertmp.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */




#if DISPLAY
  Display *disp = new Display();
  TCanvas *c_disp = new TCanvas( "c_disp", "c_disp", 0, 0, 600, 600 );
  c_disp->Divide(2,2);

    disp->SetCDSFrameXY(-70,70,-70,70);
    disp->SetCDSFrameYZ(-70,70,-70,70);
#endif

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader", &head );
  

  CDSTrackingMan *trackMan=new CDSTrackingMan();

  TFile *of=new TFile("tree/tmptr1.root","recreate");
  TTree *otree = new TTree( "CDSTrackTree", "CDSTrackTree" );
  otree->Branch( "TKOHitCol", &tko );
  otree->Branch( "CDSHitMan",  &cdsMan );
  otree->Branch( "EventHeader", &head );
  otree->Branch( "CDSTrack", &trackMan );


  /*                */
  /* event analysis */
  /*                */
  int nev = evtree->GetEntries();
  for( int iev=0; iev<nev; iev++ ){
    trackMan->Clear();
    evtree->GetEvent(iev);

    if( iev%100==0 )
      std::cout << " Event : " << iev << std::endl;

    //    if(iev%1000==0 ) c_disp->Update();    
    if(iev>5000 ) break;    
    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    //    cdsMan->Calc(conf);


    //###### event cut######
    int numCDH=0;
    for( int i=0; i<cdsMan->nCDH(); i++ )
      {
	if( 0<cdsMan->CDH(i)->tdcu() && cdsMan->CDH(i)->tdcu()<4000
	    && 0<cdsMan->CDH(i)->tdcd() && cdsMan->CDH(i)->tdcd()<4000)
	  numCDH++;
      }

    int numCDC=0;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ )
      {

	for( int i=0; i<cdsMan->nCDC(layer); i++ )
	  {
	   
	    int tdc = cdsMan->CDC(layer,i)->tdc();
	    if( 400<tdc && tdc<900) numCDC++;
	    
	  }      
      }
    

    if( !(10< numCDC &&numCDC<40 && numCDH==2 ) ) {continue;}


#if DISPLAY2
    TVirtualPad *pad;
    pad = c_disp->GetPad(1);
    disp->DrawCDSFrameXY( pad );
    disp->DrawSegmentsXY( pad, conf, CID_CDH );
    disp->DrawCDCLayersXY( pad, conf );
    disp->DrawCDSHitXY( pad, conf, cdsMan, CID_CDH );
    disp->DrawCDSHitXY( pad, conf, cdsMan, CID_CDC );
    
    pad = c_disp->GetPad(2);
    disp->DrawCDSFrameYZ( pad );
    disp->DrawCDCLayersYZ( pad, conf );
    disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDH );
    disp->DrawCDSHitYZ( pad, conf, cdsMan, CID_CDC );

    c_disp->Update();

#endif    


    //    trackMan->SetKilledLayer(8);
    trackMan->SearchCluster(cdsMan);
    trackMan->FindingTrack(cdsMan);

    //   std::cout<<"tracknum= "<<trackMan->nTrack()<<std::endl;

    if(trackMan->nTrack()!=2) {continue;}
      
    double dtheta; 
    double a[2],b[2],c[2];
    double ap[2],bp[2],cp[2];
    double x_mid[2],y_mid[2];
    
   
    dtheta=fabs( trackMan->Track(0)->Theta() - trackMan->Track(1)->Theta() );
    if(dtheta<0.01 ) {continue;}
    
    for(int i=0;i<2;i++)
      {
	a[i]=trackMan->Track(i)->A();
	b[i]=trackMan->Track(i)->B();
	c[i]=trackMan->Track(i)->C();
 
	x_mid[i]=trackMan->Track(i)->MidX();
	y_mid[i]=trackMan->Track(i)->MidY();

	ap[i]=b[i];
	bp[i]=-1;
	cp[i]=y_mid[i]-ap[i]*x_mid[i];
#if 0
	std::cout<<"x_mid["<<i<< "] "<<x_mid[i];  
	std::cout<<" y_mid["<<i<< "] "<<y_mid[i];  
	std::cout<<" ap["<<i<< "] "<<ap[i];  
	std::cout<<" bp["<<i<< "] "<<cp[i];  
	std::cout<<" cp["<<i<< "] "<<cp[i]<<std::endl; 
 	std::cout<<" inphi["<<i<< "] "<<trackMan->Track(i)->TrackHit(1,0)->phi()<<std::endl;  
#endif
      }

    double x_c,y_c,rho,rhotmp,sign;
    x_c=-(cp[0]-cp[1])/(ap[0]-ap[1]);
    y_c=ap[0]*x_c+cp[0];

    rho=sqrt( (x_c-x_mid[0])*(x_c-x_mid[0])+(y_c-y_mid[0])*(y_c-y_mid[0]) );
    rhotmp=sqrt( (x_c-x_mid[1])*(x_c-x_mid[1])+(y_c-y_mid[1])*(y_c-y_mid[1]) );

    if(fabs(rho-rhotmp)>10. || rho<1e-5)
      {
	//std::cout<<"Don't match rho!!"<<std::endl;
	{continue;}
      }

    //    if( x_c>x_mid[0] && x_c>x_mid[1] ) sign=1;
    //    else sign =-1;
    //  std::cout<<"rho ="<<rho<<" x_c y_c=  "<<x_c<<" "<<y_c<<std::endl;

    double param[5];
    double d_c=fabs( sqrt(x_c*x_c+y_c*y_c) );
    double cos_c=-x_c/d_c;
    double sin_c=-y_c/d_c;

    double x_o=x_c+rho*cos_c;
    double y_o=y_c+rho*sin_c;

    param[0]=sqrt(x_o*x_o+y_o*y_o );
    if(x_o>0 && y_o>0)   param[1]=atan(y_o/x_o);
    else if(x_o<0 && y_o>0)   param[1]=TMath::Pi()+atan(y_o/x_o);
    else if(x_o<0 && y_o<0)   param[1]=TMath::Pi()+atan(y_o/x_o);
    else if(x_o>0 && y_o<0)   param[1]=2*(TMath::Pi() )+atan(y_o/x_o);

    if(d_c>rho) param[2]=1./rho;
    else if(d_c<rho) param[2]=-1./rho;

    param[3]=-999;
    param[4]=-999;

#if 0
    std::cout<<"param ";
    for(int i=0;i<5;i++) std::cout<<param[i]<<" ";
    std::cout<<std::endl;
#endif
    //########## 2 Track Fitting ########## 

    CDSTrack *tracktmp=new CDSTrack();

    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);
	for(int numlayer=0;numlayer<15;numlayer++)
	  {
	    for(int n=0;n<track->nTrackHit(numlayer+1);n++)
	      {
		tracktmp->AddHit( *(track->TrackHit(numlayer+1,n)) );
	      }
	  }
      }

    tracktmp->SetParameters(param);
    tracktmp->FirstTrackFitting();
    for(int trial=0;trial<10;trial++)  tracktmp->TrackFitting();
    tracktmp->FirstHelixFitting();
    
    bool hflag1=true;
    for(int trial=0;trial<30;trial++)
      { if(!tracktmp->HelixFitting() ){hflag1=false; break;}  }
    if(!hflag1) {std::cout<<"skip event"<<std::endl;continue;}
    //  std::cout<<"Chi2/dof= "<<tracktmp->Chi()<<" / "<<tracktmp->Dof()<<std::endl;

    rho=fabs(1./param[2]);
    x_c=(param[0]+1./param[2])*cos(param[1]);
    y_c=(param[0]+1./param[2])*sin(param[1]);
    double Pt=0.3*0.5*rho/100.;
    // std::cout<<"Pt= "<<Pt<<"GeV/c"<<std::endl;
    // if(Pt>3 || Pt<0.05) {continue;}

    tracktmp->GetParameters(param);
    delete tracktmp;

    //####### each fit###############
    
    double arho[2],ax_c[2],ay_c[2],aPt[2];
    double aparam[2][5];
    for(int i=0;i<trackMan->nTrack();i++)
      {
	CDSTrack *track=trackMan->Track(i);
	track->SetParameters(param);
	track->SetHitPos();
	track->DeleteExtraHit();

	int numt=0;
	double Chi=999,Chitmp=999;
	bool hflag1=true;	
	while(  numt<10 )
	  {
	    Chi=Chitmp;
	    for(int trial=0;trial<30;trial++)
	      { if(!track->HelixFitting() ){hflag1=false; break;}  }
	    //    std::cout<<"Chi2/dof= "<<track->Chi()<<" / "<<track->Dof()<<std::endl;
	    Chitmp=track->Chi();
	    if( (Chi-Chitmp)<0.3 ) break;
	    numt++;
	  }
	if(!hflag1) {std::cout<<"skip event"<<std::endl;continue;}
	if( Chitmp>2000 ) {std::cout<<"skip2 event"<<std::endl;continue;}

	track->GetParameters(aparam[i]);
	arho[i]=fabs(1./aparam[i][2]);
	ax_c[i]=(aparam[i][0]+1./aparam[i][2])*cos(aparam[i][1]);
	ay_c[i]=(aparam[i][0]+1./aparam[i][2])*sin(aparam[i][1]);
	//    std::cout<<"rho= "<<rho<<" x_c= "<<x_c<<" y_c= "<<y_c<<std::endl;
	aPt[i]=0.3*0.5*arho[i]/100.;
	//	std::cout<<"Pt"<<i<<"= "<<aPt[i]<<"GeV/c"<<std::endl;

#if 0
	std::cout<<"aparam"<<i<<" ";
    for(int n=0;n<5;n++) std::cout<<aparam[i][n]<<" ";
    std::cout<<std::endl;
#endif 
    
      }
    otree->Fill();
    trackMan->Clear();
   
#if DISPLAY2 
    pad = c_disp->GetPad(1);
    pad->cd();
    TArc arc1,arc2,arc3;


    arc1.SetFillStyle(0);
    arc1.SetLineColor(2);

    arc2.SetFillStyle(0);
    arc2.SetLineColor(3);

    arc3.SetFillStyle(0);
    arc3.SetLineColor(4);

    arc1.DrawArc(x_c,y_c,rho,0,360);
    arc2.DrawArc(ax_c[0],ay_c[0],arho[0],0,360);
    arc3.DrawArc(ax_c[1],ay_c[1],arho[1],0,360);
      
    TF1 *func_y[3];
    func_y[0]=new TF1("func_y0","([0])*sin([1])+[2]*(sin([1])-sin( [1]+([3]-x)/([4]*[2]) ) )",param[3]-1./param[2]*param[4]*(-TMath::Pi()/2. ),param[3]-1./param[2]*param[4]*(TMath::Pi()/2. ) );
    func_y[1]=new TF1("func_y1","([0])*sin([1])+[2]*(sin([1])-sin( [1]+([3]-x)/([4]*[2]) ) )",aparam[0][3]-aparam[0][2]*aparam[0][4]*(-TMath::Pi()/2. ),aparam[0][3]-aparam[0][2]*aparam[0][4]*(TMath::Pi()/2. ) );
    func_y[2]=new TF1("func_y2","([0])*sin([1])+[2]*(sin([1])-sin( [1]+([3]-x)/([4]*[2]) ) )",aparam[1][3]-aparam[1][2]*aparam[1][4]*(-TMath::Pi()/2. ),aparam[1][3]-aparam[1][2]*aparam[1][4]*(TMath::Pi()/2. ) );

    func_y[0]->SetParameters(param[0],param[1],1./param[2],param[3],param[4]);
    func_y[1]->SetParameters(aparam[0][0],aparam[0][1],aparam[0][2],aparam[0][3],aparam[0][4]);
    func_y[2]->SetParameters(aparam[1][0],aparam[1][1],aparam[1][2],aparam[1][3],aparam[1][4]);
    func_y[0]->SetLineColor(2);
    func_y[1]->SetLineColor(3);
    func_y[2]->SetLineColor(4);
    c_disp->cd(2);
    func_y[0]->Draw("LPSAME");
    func_y[1]->Draw("LPSAME");
    func_y[2]->Draw("LPSAME");

    c_disp->Update();    


    bool status = disp->Wait();
    if( !status ) return;
    c_disp->Update();    
#endif



  }

    of->Write();
    of->Close();

}

