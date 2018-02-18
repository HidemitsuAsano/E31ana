const double MassThPPi = 0.8; // GeV/c2
const double TOFOffset = 3; // nsec
const double piMass = 0.138;
const double pMass = 0.934;
#include "TMinuit.h"
#include "TVector3.h"
void drawimpipi()
{
  
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
 

  int nrun = sizeof(run)/sizeof(int);
  std::cout << " nrun:" << nrun << std::endl;
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("IMpipi");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("IMpipi");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("Invariant Mass of #pi#pi [GeV/c^2]");
  h1->SetAxisRange( 0.25, 0.6, "X" );
  h1->Draw();

  c1->Print("impipi.ps");
}

void drawimppi()
{
  
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
  
  //  int run[] = {3060,3066,3067};
  int nrun = sizeof(run)/sizeof(int);
  std::cout << " nrun:" << nrun << std::endl;
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("IMppi");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("IMppi");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("Invariant Mass of p#pi [GeV/c^2]");
  h1->SetAxisRange( 1.05, 1.25, "X" );
  h1->Draw();
  h1->SetName("IMppisum");
  c1->Print("imppi.ps");
}

void drawimpipi_nopid()
{
  
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
 

  int nrun = sizeof(run)/sizeof(int);
  std::cout << " nrun:" << nrun << std::endl;
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("IMpipi_nopid");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("IMpipi_nopid");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("Invariant Mass of #pi#pi [GeV/c^2] no PID");
  h1->SetAxisRange( 0.25, 0.6, "X" );
  h1->Draw();

  c1->Print("impipi_nopid.ps");
}


void drawphi()
{
  
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
 

  int nrun = sizeof(run)/sizeof(int);
  std::cout << " nrun:" << nrun << std::endl;
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("phi");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("phi");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("Phi dis of track");
  //h1->SetAxisRange( 0.25, 0.6, "X" );
  h1->Draw();

  c1->Print("phi.ps");
}


void drawcosOA()
{
  
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
 

  int nrun = sizeof(run)/sizeof(int);
  std::cout << " nrun:" << nrun << std::endl;
  TFile *f;
  TH1F *h1,*h2;
  TH1F *hpipi1,*hpipi2;
  TH1F *h2dppi1,*h2dppi2;
  TH1F *h2dpipi1,*h2dpipi;

  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("cosOAppi");
  h2dppi1 = (TH1F*)f->Get("IMcosOAppi");
  hpipi1 = (TH1F*)f->Get("cosOApipi");
  h2dpipi1 = (TH1F*)f->Get("IMcosOApipi");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("cosOAppi");
    h2dppi2 = (TH1F*)f->Get("IMcosOAppi");
    hpipi2 = (TH1F*)f->Get("cosOApipi");
    h2dpipi2 = (TH1F*)f->Get("IMcosOApipi");

    h1->Add( h1, h2, 1, 1 );
    h2dppi1->Add( h2dppi1, h2dppi2, 1, 1 );
    hpipi1->Add( hpipi1, hpipi2, 1, 1 );
    h2dpipi1->Add( h2dpipi1, h2dpipi2, 1, 1 );

  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  // h1->SetXTitle("Invariant Mass of #pi#pi [GeV/c^2] no PID");
  //  h1->SetAxisRange( 0.25, 0.6, "X" );
  h1->Draw();

  c1->Print("cosOAppi.ps");


  TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
  // h1->SetXTitle("Invariant Mass of #pi#pi [GeV/c^2] no PID");
  //  h1->SetAxisRange( 0.25, 0.6, "X" );
  h2dppi1->Draw("colz");

  c2->Print("IMcosOAppi.ps");


  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  // h1->SetXTitle("Invariant Mass of #pi#pi [GeV/c^2] no PID");
  //  h1->SetAxisRange( 0.25, 0.6, "X" );
  hpipi1->Draw();

  c3->Print("cosOApipi.ps");


  TCanvas *c2 = new TCanvas( "c4", "c4", 800, 600 );
  // h1->SetXTitle("Invariant Mass of #pi#pi [GeV/c^2] no PID");
  //  h1->SetAxisRange( 0.25, 0.6, "X" );
  h2dpipi1->Draw("colz");

  c4->Print("IMcosOApipi.ps");
}



void drawmom()
{
  //  int run[] = {3060,3066,3067};
  int run[] = {3038};
  int nrun = sizeof(run)/sizeof(int);
  TFile *f;
  TH1F *h1;
  TH1F *h_mom, *h_momp, *h_mompip, *h_mompim;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h_mom = (TH1F*)f->Get("mom");
  h_momp = (TH1F*)f->Get("momp");
  h_mompip = (TH1F*)f->Get("mompip");
  h_mompim = (TH1F*)f->Get("mompim");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h1 = (TH1F*)f->Get("mom"); h_mom->Add( h_mom, h1, 1, 1 );
    h1 = (TH1F*)f->Get("momp"); h_mom->Add( h_momp, h1, 1, 1 );
    h1 = (TH1F*)f->Get("mompip"); h_mom->Add( h_mompip, h1, 1, 1 );
    h1 = (TH1F*)f->Get("mompim"); h_mom->Add( h_mompim, h1, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h_mom->SetXTitle("Momentum [GeV/c]");
  h_mom->Draw();
  h_momp->SetLineColor(2);
  h_momp->Draw("same");
  h_mompip->SetLineColor(3);
  h_mompip->Draw("same");
  h_mompim->SetLineColor(4);
  h_mompim->Draw("same");

  c1->Print("mom.ps");
}



void drawvrx()
{
  //  int run[] = {3060,3061,3065,3066,3067};
  //  int run[] = {3060,3061,3065,3066,3067,3068,3069,3070,3071,3072,
  //       3073,3074,3075,3076,3077,3078,3079};
  /*
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
  */
  int run[] = {3029,3030,3031,3032,3033,3034};
  int nrun = sizeof(run)/sizeof(int);
  TFile *f;
  TH1F *h1;
  TH1F *h_vx, *h_vy, *h_xz;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );

  h_vx = (TH1F*)f->Get("vx");
  h_vy = (TH1F*)f->Get("vy");
  h_vz = (TH1F*)f->Get("vz");

  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h1 = (TH1F*)f->Get("vx"); h_vx->Add( h_vx, h1, 1, 1 );
    h1 = (TH1F*)f->Get("vy"); h_vy->Add( h_vy, h1, 1, 1 );
    h1 = (TH1F*)f->Get("vz"); h_vz->Add( h_vz, h1, 1, 1 );

  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );
  h_vx->SetXTitle("Vertex Point [cm]");
  h_vx->SetLineColor(2);
  h_vx->Draw();
  c1->Print("vrx.ps(");

  h_vy->SetXTitle("Vertex Point [cm]");
  h_vy->SetLineColor(2);
  h_vy->Draw("");
  c1->Print("vrx.ps");

  h_vz->SetXTitle("Vertex Point [cm]");
  h_vz->SetLineColor(2);
  h_vz->Draw("");
  h_vz->SetName("vzsum");
  c1->Print("vrx.ps)");
}


void drawpid()
{
  /*
  int run[] = {3057,3058,3059,
               3060,3061,     3063,3064,3065,3066,3067,3068,3069,
	       3070,3071,3072,3073,3074,3075,3076,3077,3078,3079,
	       3080,3081,3082,    ,3084,3085,3086,3087,3088,3089.
	       3090,3091,3092,3093,3094,3095,3096,3097,3098,3099,
	       3100,3101,3102,3103,3104,3105,3106};
  */
  int run[] = {3057,3058,3059
	       //               3063,3064,3065,

	       /* 3068,3069,3070*/};

  int nrun = sizeof(run)/sizeof(int);
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("pid");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("pid");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("TOF between CDH and T0 [ns]");
  h1->SetYTitle("Momentum [GeV/c]");
  h1->Draw("colz");

  TF1 *fun1 = new TF1("fun1",Form("%lf/sqrt( pow((x+%lf),2)-1)",MassThPPi,TOFOffset),-2,14);
  fun1->SetNpx(999);
  fun1->Draw("same");
  TF1 *fun2 = new TF1("fun2",Form("-%lf/sqrt( pow((x+%lf),2)-1)",MassThPPi,TOFOffset),-2,14);
  fun2->SetNpx(999);
  fun2->Draw("same");

  c1->Print("pid.ps");
}

void drawvzppi()
{
  int run[] = {3060,3066,3067};
  int nrun = sizeof(run)/sizeof(int);
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_new_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("vzppi");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_new_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("vzppi");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("TOF between CDH and T0 [ns]");
  h1->SetYTitle("Momentum [GeV/c]");
  h1->Draw();

  c1->Print("vzppi.ps");
}

void drawvppi()
{
  int run[] = {3060,3066,3067};
  int nrun = sizeof(run)/sizeof(int);
  TFile *f;
  TH1F *h1,*h2;
  f = new TFile( Form("./root/trout_re2_%d.root",run[0]) );
  h1 = (TH1F*)f->Get("vppi");
  for( int i=1; i<nrun; i++ ){
    f = new TFile( Form("./root/trout_re2_%d.root",run[i]) );
    h2 = (TH1F*)f->Get("vppi");
    h1->Add( h1, h2, 1, 1 );
  }

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  h1->SetXTitle("TOF between CDH and T0 [ns]");
  h1->SetYTitle("Momentum [GeV/c]");
  h1->Draw();

  c1->Print("vppi.ps");
}

void track(int runnum=3029)
{
  
  //runnum = 3060;
  //  int run[] = {1,2,3,4,5,6,7,8,9};
  //   int run[] = {1,2,3,4,5,6};
  //int run[] = {1};
  int run[] = {0};
  int nrun = sizeof(run)/sizeof(int);
  
  gSystem->Load("libPhysics.so");
  gSystem->Load( "./lib/libAll.so" );
  
  ConfMan *conf = new ConfMan( "conf/Oct2010/analyzer.conf" );

  MathTools *mtool = new MathTools();

  TFile *f[10];
  TTree *evtree[10];
  for( int i=0; i<nrun; i++ ){
    //f[i] = new TFile( Form( "./root/track_%d_%d.root", runnum,run[i] ) );
    //       f[i] = new TFile( Form( "./root/track_re_%d.root", runnum ) );
    f[i] = new TFile( Form( "./root/track_new_%d.root", runnum ) );
    evtree[i] = (TTree*)f[i]->Get( "EventTree" );
  }

  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan = 0;
  
  TFile *fout = new TFile( Form("./root/trout_new_%d.root",runnum), "recreate" );
  //TFile *fout = new TFile( "./root/trout.root", "recreate" );
  
  new TH1F( "nTrack", "nTrack", 200, 0, 100 );
  new TH1F( "chi", "chi", 200, 0, 100 );
  new TH1F( "tofBHDT0", "tof BHD and T0", 200, -5, 15 );
  new TH1F( "mom", "mom", 500, -1, 1 );
  new TH2F( "pid", "pid", 200, -5, 15, 300, -1.5, 1.5 );
  new TH1F( "phi", "phi of track", 200, 0, 360 );
  new TH1F( "momp", "mom proton", 500, -1, 1 );
  new TH1F( "mompip", "mom piplus", 500, -1, 1 );
  new TH1F( "mompim", "mom piminus", 500, -1, 1 );
  new TH1F( "vx", "vx", 500, -25, 25 );
  new TH1F( "vy", "vy", 500, -25, 25 );
  new TH1F( "vz", "vz", 1000, -50, 50 );
  new TH2F( "vxy", "vxy", 500, -25, 25 , 500, -25, 25 );
  new TH1F( "vrx_dis", "vrx_dis", 2000, 0, 10 );

  new TH1F( "vxppi", "vxppi", 500, -25, 25 );
  new TH1F( "vyppi", "vyppi", 500, -25, 25 );
  new TH1F( "vzppi", "vzppi", 1000, -50, 50 );
  new TH3F( "vppi",  "vppi", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "IMppi", "IM proton piminus", 500, 1.0, 2 );
  new TH1F( "IMppi_m", "IM proton piminus p*0.99", 500, 1.0, 2 );
  new TH1F( "IMppi_p", "IM proton piminus p*1.01", 500, 1.0, 2 );
  new TH1F( "IMppiKselected", "IM proton piminus", 500, 1.0, 2 );
  new TH1F( "cosOAppi", "cosOA between proton and #pi^-", 200, -1, 1 );
  new TH2F( "IMcosOAppi", "IM and cosOA for ppi", 200, 1, 2, 200, -1, 1);
  new TH1F( "costppi", "cost of ppi system", 200, -1, 1 );
  new TH1F( "phippi", "phi of ppi system", 200, -180, 180 );
  new TH2F( "IMcostppi", "IM and cost of ppi system", 500, 1, 2, 200, -1, 1 );

  new TH1F( "vxpipi", "vxpipi", 500, -25, 25 );
  new TH1F( "vypipi", "vypipi", 500, -25, 25 );
  new TH1F( "vzpipi", "vzpipi", 1000, -50, 50 );



  new TH3F( "vpipi",  "vpipi", 200, -25, 25, 200, -25, 25, 200, -25, 25 );
  new TH1F( "IMpipi", "IM piplus piminus", 500, 0, 1.0 );
  new TH1F( "IMpipi_m", "IM piplus piminus p*0.99", 500, 0, 1.0 );
  new TH1F( "IMpipi_p", "IM piplus piminus p*1.01", 500, 0, 1.0 );
  new TH1F( "IMpipi_cut", "IM piplus piminus cut OA", 500, 0, 1.0 );
  new TH1F( "IMpipi_nopid", "IM piplus piminus nopid", 500, 0, 1.0 );
  new TH1F( "cosOApipi", "cosOA between #pi^+ and #pi^-", 200, -1, 1 );
  new TH2F( "IMcosOApipi", "IM and cosOA for pipi", 500, 0, 1, 200, -1, 1);
  new TH1F( "costpipi", "cost of pipi system", 200, -1, 1);
  new TH1F( "phipipi", "phi of pipi system", 200, -180, 180 );
  new TH2F( "IMcostpipi", "IM and cost of pipi system", 500, 0, 1, 200, -1, 1 );
  

  TH1F *h1;
  TH2F *h2;
  TH3F *h3;
  for( int irun=0; irun<nrun; irun++ ){
    evtree[irun]->SetBranchAddress( "CDSHitMan", &cdsMan );
    evtree[irun]->SetBranchAddress( "BeamLineHitMan", &blMan );
    evtree[irun]->SetBranchAddress( "EventHeader", &head );
    evtree[irun]->SetBranchAddress( "CDSTrackingMan", &trackMan );
    int nev = evtree[irun]->GetEntries();
    std::cout << " i:" << i << " entry:" << nev << std::endl;

    for( int iev=0; iev<nev; iev++ ){
      if( iev%2000 == 0 )
	std::cout << " Event : " << iev << std::endl;
      evtree[irun]->GetEvent(iev);
      //cdsMan->Calc(conf);

      //cdsMan->CheckContainerSize();
      //blMan->CheckContainerSize();
      
      int nT0=0;
      for( int i=0; i<blMan->nT0(); i++ ){
	if( blMan->T0(i)->CheckRange() ) nT0++;
      }
      if( nT0!=1 ) continue;
      double ctmT0;
      for( int i=0; i<blMan->nT0(); i++ ){
	if( blMan->T0(i)->CheckRange() ){
	  ctmT0 = blMan->T0(i)->ctmean();
	}
      }

    int GoodTrack=0;
      for( int it=0; it<trackMan->nTrack(); it++ ){
	CDSTrack *track = trackMan->Track(it);
	double chi = track->Chi();
	if(chi<20) GoodTrack++;
      }
            
      int nBHD=0;
       for( int i=0; i<blMan->nBHD(); i++ ){
	 if( blMan->BHD(i)->CheckRange() ) nBHD++;
       }
       //       if( nBHD!=1 ) continue;
       double ctmBHD;
       for( int i=0; i<blMan->nBHD(); i++ ){
	 if( blMan->BHD(i)->CheckRange() ){
 	  ctmBHD = blMan->BHD(i)->ctmean();
 	}
       }
       double tofBHDT0 = ctmT0-ctmBHD;
       h1 = (TH1F*)gFile->Get("tofBHDT0"); h1->Fill( tofBHDT0 );
       int pid_beam;//0:pi 1:K 3:else

       if(-2<tofBHDT0 && tofBHDT0<2) pid_beam=0;
       else if(2<tofBHDT0 && tofBHDT0<5) pid_beam=1;
       else pid_beam=3;
       //if(pid_beam!=0) continue;
      
      //std::cout << " nTrack:" << trackMan->nTrack() << std::endl;
      h1 = (TH1F*)gFile->Get("nTrack"); h1->Fill( trackMan->nGoodTrack() );
      for( int it=0; it<trackMan->nGoodTrack(); it++ ){
	CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
	HodoscopeLikeHit *cdh = track->CDHHit();
	double tof = cdh->ctmean()-ctmT0;
	double chi = track->Chi();
	double param[5];
	track->GetParameters(param);
	if(param[2]==0 ||param[4]==0   ) continue;
	double drho=param[0], phi0=param[1], rho=1./param[2], dz=param[3], tlam=param[4];
	double mom = 0.003*0.5*rho;
	mom=mom*sqrt(1+param[4]*param[4]);
	h1 = (TH1F*)gFile->Get("chi"); h1->Fill( chi );
	h1 = (TH1F*)gFile->Get("mom"); h1->Fill( mom );
	h2 = (TH2F*)gFile->Get("pid"); h2->Fill( tof, mom );
	int ptype = PID( tof, mom );
	if( ptype==1 ){ h1 = (TH1F*)gFile->Get("mompip"); h1->Fill( mom ); }
	if( ptype==2 ){ h1 = (TH1F*)gFile->Get("mompim"); h1->Fill( mom ); }
	if( ptype==3 ){ h1 = (TH1F*)gFile->Get("momp"); h1->Fill( mom ); }

	double phi=track->Theta();
	h1 = (TH1F*)gFile->Get("phi"); h1->Fill( phi);
      }

      for( int it1=0; it1<trackMan->nGoodTrack(); it1++ ){

	CDSTrack *track1 = trackMan->Track(trackMan->GoodTrackID(it1));
	HodoscopeLikeHit *cdh1 = track1->CDHHit();
	double tof1 = cdh1->ctmean()-ctmT0;
	double param1[5];
	track1->GetParameters(param1);
	if(param1[2]==0 ||param1[4]==0   ) continue;
	double drho1=param1[0], phi01=param1[1], rho1=1./param1[2], dz1=param1[3], tlam1=param1[4];
	double mom1 = 0.003*0.5*rho1;
	mom1=mom1*sqrt(1+param1[4]*param1[4]);
	int ptype1 = PID( tof1, mom1 );
	for( int it2=it1+1; it2<trackMan->nGoodTrack(); it2++ ){
	  CDSTrack *track2 = trackMan->Track(trackMan->GoodTrackID(it2));
	  HodoscopeLikeHit *cdh2 = track2->CDHHit();
	  double tof2 = cdh2->ctmean()-ctmT0;
	  double param2[5];
	  track2->GetParameters(param2);

	  if(param2[2]==0 ||param2[4]==0   ) continue;

	  double drho2=param2[0], phi02=param2[1], rho2=1./param2[2], dz2=param2[3], tlam2=param2[4];
	  double mom2 = 0.003*0.5*rho2;
	  mom2=mom2*sqrt(1+param[4]*param2[4]);
	  int ptype2 = PID( tof2, mom2 );
	  
	  if(trackMan->nGoodTrack()==2)
	    {
	      TVector3 vtx1, vtx2, vtx;
	      //    mtool->CalcHelixDCA( param1, param2, vtx1, vtx2, vtx );
	      vtx=trackMan->GetVertex(trackMan->GoodTrackID(it1),trackMan->GoodTrackID(it2));
	      h1 = (TH1F*)gFile->Get("vx"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vy"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vxy"); h1->Fill( vtx.X(),vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vz"); h1->Fill( vtx.Z() );

	      TVector3 nest;
	      double dis;
	      //    cout<<"test1 "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<endl;
	      track1->PointToHelix(vtx,param1,nest,dis);
	      // cout<<"test2 "<<dis<<endl;
	      h1 = (TH1F*)gFile->Get("vrx_dis"); h1->Fill( dis );

	    }



	  if( (ptype1!=3 || ptype1!=0)  && (ptype2!=0 ||ptype2!=3) ){// no use PID
	    CDSTrack *tracktmp; 
	    HodoscopeLikeHit *cdhtmp;
	    double paramtmp[5];
	    
	    double mass1, mass2;
	    mass1 = piMass; mass2 = piMass; 
	    TVector3 vtx1, vtx2, vtx;

	    //	    mtool->CalcHelixDCA( param1, param2, vtx1, vtx2, vtx );
	    vtx=trackMan->GetVertex(trackMan->GoodTrackID(it1),trackMan->GoodTrackID(it2));
	    TVector3 Pp1 = mtool->CalcHelixMom(  param1, vtx.Z() );
	    TVector3 Pp2 = mtool->CalcHelixMom(  param2, vtx.Z() );

	    TVector3 Pp = Pp1 + Pp2;
	    TLorentzVector L_p1; L_p1.SetVectM( Pp1, mass1 );
	    TLorentzVector L_p2; L_p2.SetVectM( Pp2, mass2);

	    double cosOA = (Pp1.Dot(Pp2))/Pp1.Mag()/Pp2.Mag();
	    double im = (L_p1+L_p2).M();
	      h1 = (TH1F*)gFile->Get("IMpipi_nopid"); h1->Fill( im );

	  }
	  

	  if( ((ptype1==3&&ptype2==2) || (ptype1==2&&ptype2==3)) ||// select p pi-
		((ptype1==1&&ptype2==2) || (ptype1==2&&ptype2==1)) ){// select pi+ pi-
	    CDSTrack *tracktmp; 
	    HodoscopeLikeHit *cdhtmp;
	    double paramtmp[5];
	    
	    if( (ptype1==3&&ptype2==2) || (ptype1==1&&ptype2==2) ){}
	    else{
	      tracktmp = track1; track1 = track2; track2 = tracktmp;
	      cdhtmp = cdh1; cdh1 = cdh2; cdh2 = cdhtmp;
	      for( int i=0; i<5; i++ ){ paramtmp[i] = param1[i]; param1[i] = param2[i]; param2[i] = paramtmp[i]; }
	    }
	    
	    double mass1, mass2;
	    if( ptype1==3 ){ mass1 = pMass;  mass2 = piMass; }
	    else           { mass1 = piMass; mass2 = piMass; }
	    TVector3 vtx1, vtx2, vtx;

	    mtool->CalcHelixDCA( param1, param2, vtx1, vtx2, vtx );
	    //  vtx=trackMan->GetVertex(trackMan->GoodTrackID(it1),trackMan->GoodTrackID(it2));
	    TVector3 Pp1 = mtool->CalcHelixMom(  param1, vtx1.Z() );
	    TVector3 Pp2 = mtool->CalcHelixMom(  param2, vtx2.Z() );

	    TVector3 Pp = Pp1 + Pp2;
	    TLorentzVector L_p1; L_p1.SetVectM( Pp1, mass1 );
	    TLorentzVector L_p2; L_p2.SetVectM( Pp2, mass2);

	    TLorentzVector L_p1_p; L_p1_p.SetVectM( Pp1*1.01, mass1 );
	    TLorentzVector L_p2_p; L_p2_p.SetVectM( Pp2*1.01, mass2);

	    TLorentzVector L_p1_m; L_p1_m.SetVectM( Pp1*0.99, mass1 );
	    TLorentzVector L_p2_m; L_p2_m.SetVectM( Pp2*0.99, mass2);

	    double cosOA = (Pp1.Dot(Pp2))/Pp1.Mag()/Pp2.Mag();
	    double im = (L_p1+L_p2).M();
	    double im_p = (L_p1_p+L_p2_p).M();
	    double im_m = (L_p1_m+L_p2_m).M();
	    if( ptype1==3 ){ // ppi-
	      h1 = (TH1F*)gFile->Get("vxppi"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vyppi"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vzppi"); h1->Fill( vtx.Z() );
	      h3 = (TH3F*)gFile->Get("vppi");  h3->Fill( vtx.X(), vtx.Y(), vtx.Z() );
	      h1 = (TH1F*)gFile->Get("cosOAppi"); h1->Fill( cosOA );
	      h1 = (TH1F*)gFile->Get("IMppi"); h1->Fill( im );
	      h1 = (TH1F*)gFile->Get("IMppi_m"); h1->Fill( im_m );
	      h1 = (TH1F*)gFile->Get("IMppi_p"); h1->Fill( im_p );
	      h2 = (TH2F*)gFile->Get("IMcosOAppi"); h2->Fill( im, cosOA );
 	      h1 = (TH1F*)gFile->Get("costppi"); h1->Fill( Pp.CosTheta() );
 	      h1 = (TH1F*)gFile->Get("phippi"); h1->Fill( Pp.Phi()*TMath::RadToDeg() );
 	      h2 = (TH2F*)gFile->Get("IMcostppi"); h2->Fill( im, Pp.CosTheta() );
	      //   std::cout<<"Phi "<<Pp.Phi()*TMath::RadToDeg()<<std::endl;
// 	      if( 2<tofBHDT0 ){
// 		h1 = (TH1F*)gFile->Get("IMppiKselected"); h1->Fill( im );
// 	      }
	    }
	    else if( ptype1==1){ // pi+pi-
	      
	      h1 = (TH1F*)gFile->Get("vxpipi"); h1->Fill( vtx.X() );
	      h1 = (TH1F*)gFile->Get("vypipi"); h1->Fill( vtx.Y() );
	      h1 = (TH1F*)gFile->Get("vzpipi"); h1->Fill( vtx.Z() );
	      h3 = (TH3F*)gFile->Get("vpipi");  h3->Fill( vtx.X(), vtx.Y(), vtx.Z() );
	      h1 = (TH1F*)gFile->Get("cosOApipi"); h1->Fill( cosOA );
	      h1 = (TH1F*)gFile->Get("IMpipi"); h1->Fill( im );
	      h1 = (TH1F*)gFile->Get("IMpipi_m"); h1->Fill( im_m );
	      h1 = (TH1F*)gFile->Get("IMpipi_p"); h1->Fill( im_p );

	      h2 = (TH2F*)gFile->Get("IMcosOApipi"); h2->Fill( im, cosOA );
	      h1 = (TH1F*)gFile->Get("costpipi"); h1->Fill( Pp.CosTheta() );
	      h1 = (TH1F*)gFile->Get("phipipi"); h1->Fill( Pp.Phi()*TMath::RadToDeg() );
 	      h2 = (TH2F*)gFile->Get("IMcostpipi"); h2->Fill( im, Pp.CosTheta() );
	      if(  (cosOA < 0.6 && cosOA > -0.95)) h1 = (TH1F*)gFile->Get("IMpipi_cut"); h1->Fill( im );
	    }

	  }
	}
      }
    }
  }
  
  gFile->Write();
  gFile->Close();
  
}

int PID( double tof, double mom ) // 0:unknown, 1:pi+, 2;pi-, 3:proton
{
  //TF1 *fun1 = new TF1("fun1",Form("%lf/sqrt( pow((x+%lf),2)-1)",MassThPPi,TOFOffset),-2,14);
  int val=0;
  if( tof<TOFOffset ){
  }
  if( 0<mom ){
    if( mom < MassThPPi/sqrt(pow(tof+TOFOffset,2)-1) ){
      val = 1;
    }
    else{
      val = 3;
    }
  }
  else{
    if( fabs(mom) < MassThPPi/sqrt(pow(tof+TOFOffset,2)-1) ){
      val = 2;
    }
    else{
      val = 0;
    }
  }

  //std::cout << " tof:" << tof << " mom:" << mom << " type:" << val << std::endl;
  return val;
}

double CalcHelixDCA2(double *par1, double *par2,
		     TVector3& vtx1, TVector3& vtx2, TVector3& vtx)
{
  if(par1[2]==0 || par1[4]==0 ||par2[2]==0 || par2[4]==0 ) 
    {vtx1.SetXYZ(-999,-999,-999);vtx2.SetXYZ(-999,-999,-999);vtx.SetXYZ(-999,-999,-999); return -999;}  
  const int npoint = 100;
  const double initz = 0; //cm
  double region = 60.0; //+/-cm
  double nowz = initz;

  double x1,y1,z1,x2,y2,z2,x,y,z;
  double x1c,y1c,z1c,x2c,y2c,z2c,xc,yc,zc;
  double dl;
  double dlmax = 99999;
  for(int i1=0; i1<npoint+1; i1++){
    for(int i2=0; i2<npoint+1; i2++){
      z1 = nowz+(2*region/npoint)*(-npoint/2+i1);
      z2 = nowz+(2*region/npoint)*(-npoint/2+i2);
      CalcHelixPos2(par1,z1,x1,y1);
      CalcHelixPos2(par2,z2,x2,y2);
      x = (x1+x2)/2.; y = (y1+y2)/2.; z = (z1+z2)/2.; 
      dl = sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) );
      if( dl<dlmax ){
	x1c=x1; y1c=y1; z1c=z1; x2c=x2; y2c=y2; z2c=z2; xc=x; yc=y; zc=z;
	dlmax = dl;
      }
    }
  }
  vtx1.SetXYZ(x1c,y1c,z1c);
  vtx2.SetXYZ(x2c,y2c,z2c);
  vtx.SetXYZ(xc,yc,zc);
  return dlmax;
}

void CalcHelixPos2( double *par, double z, double &x, double &y)
{
  if(par[2]==0 || par[4]==0  ) 
  {x=-999;y=-999;return; }  
  double phi2 = par[2]/par[4]*(par[3]-z);
  double x    = (par[0]+1./par[2] )*cos(par[1])-1./par[2]*cos(par[1]+phi2);
  double y    = (par[0]+1./par[2] )*sin(par[1])-1./par[2]*sin(par[1]+phi2);
  //return TVector3(x,y,z);
}


double CalcHelixDCA(const double par1[5], const double par2[5],
		    TVector3& vtx1, TVector3& vtx2, TVector3& vtx)
{
  if(par1[2]==0 || par1[4]==0 ||par2[2]==0 || par2[4]==0  ) 
  {vtx1.SetXYZ(-999,-999,-999);vtx2.SetXYZ(-999,-999,-999);vtx.SetXYZ(-999,-999,-999);return -999; }  
  const int NUM = 2;
  const int npoint = 100;
  const double initz = 0; //cm
  double region = 60.0; //+/-cm
  TVector3 pos1[npoint+1];
  TVector3 pos2[npoint+1];
  int num = 1;
  TVector3 now_vtx1;
  TVector3 now_vtx2;
  TVector3 now_vtx;
  double nowz = initz;
  double minl = 999.0;
  while(num<=NUM){
    //cerr<<"---"<<num<<endl;
    for(int i=0; i<npoint+1; i++){
      double z=nowz+(2*region/npoint)*(-npoint/2+i);
      //cerr<<z<<endl;
      pos1[i] = CalcHelixPos(par1,z);
      pos2[i] = CalcHelixPos(par2,z);
    }
    for(int i=0; i<npoint+1; i++){
      for(int j=0; j<npoint+1; j++){
	TVector3 diff = pos1[i]-pos2[j];
        double l = diff.Mag();
        if(l<minl){
          minl = l;
          now_vtx1 = pos1[i];
          now_vtx2 = pos2[j];
          now_vtx = now_vtx1+now_vtx2;
	  now_vtx *=0.5;
          nowz = now_vtx.z();
        }
      }
    }
    num++;
    region = 2*region/10;
  }
  if(now_vtx1.x()==now_vtx1.y()==now_vtx1.z()==0
     ||now_vtx2.x()==now_vtx2.y()==now_vtx2.z()==0
     ||now_vtx.x()==now_vtx.y()==now_vtx.z()==0)
  {vtx1.SetXYZ(-999,-999,-999);vtx2.SetXYZ(-999,-999,-999);vtx.SetXYZ(-999,-999,-999);return -999; }  
  vtx1 = now_vtx1;
  vtx2 = now_vtx2;
  vtx = now_vtx;
  return minl;
}

TVector3 CalcHelixPos(const double par[5], double z)
{
  if(par[2]==0 || par[4]==0  ) 
  {return TVector3(-999,-999,-999); }  
  double phi2 = par[2]/par[4]*(par[3]-z);
  double x    = (par[0]+1./par[2] )*cos(par[1])-1./par[2]*cos(par[1]+phi2);
  double y    = (par[0]+1./par[2] )*sin(par[1])-1./par[2]*sin(par[1]+phi2);
  return TVector3(x,y,z);
}


TVector3 CalcHelixMom(const double par[5], double z)
{
  if(par[2]==0 || par[4]==0  ) 
  {return TVector3(-999,-999,-999); }  
  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = -1*0.5; //T, "-1" is needed.
  double phi2 = (par[3]-z)*par[2]/par[4];
  double pt = fabs(1/par[2])*(Const*dMagneticField)/100.;
  double px = pt*(-1*sin(par[1]+phi2));
  double py = pt*(cos(par[1]+phi2));
  double pz = pt*(par[4]);
  return TVector3(px,py,pz);
}

double CalcInvMass(const TVector3 & piMom, const TVector3& pMom)
{
  double piE = sqrt(piMass*piMass+piMom.Mag2());
  double pE  = sqrt(pMass*pMass+pMom.Mag2());
  double E   = piE+pE;
  TVector3 P = piMom+pMom;
  double M   = E*E-P.Mag2();
  if(M>0) return sqrt(M);
  else    return 0;
}
