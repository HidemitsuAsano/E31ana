static int sEvNum = -1;

void disp( int evnum = -1 )
{
  gSystem->Load("./lib/libAll.so");
  sEvNum++;

  if( 0<evnum ) sEvNum = evnum;
 
  TFile *f = new TFile("./root/tr3086.root");
  TTree *evtree = (TTree*)f->Get("EventTree");
    
  ConfMan *conf = new ConfMan( "conf/Oct2010/analyzer.conf" );
  conf->Initialize();

  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  CDSTrackingMan *trackMan = 0;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader",  &head );
  evtree->SetBranchAddress( "CDSTrackingMan",  &trackMan );
  
  int nev = evtree->GetEntries();
  if( nev<=sEvNum ) return;

  std::cout << " Event# : " << sEvNum << std::endl;
  evtree->GetEvent(sEvNum);

  //cdsMan->Calc(conf);

  if( trackMan->nTrack()<2 ) return;

  double nGoodTrack=0;
  for( int i=0; i<trackMan->nTrack(); i++ ){
    double chi = trackMan->Track(i)->Chi();
    if( chi<10 ) nGoodTrack++;
  }
  if( nGoodTrack<2 ) return;

  trackMan->CalcVertex();

  // first display
  TCanvas *c1 = new TCanvas( "c1", "c1", 200, 200, 500, 500 );
  Display3D *disp = new Display3D(c1);
  //disp->SetCDSSolenoidBody(conf,180,40,40);
  //disp->SetCDSSolenoidCap(conf,0xb,40);
  disp->SetTargetOct2010(conf,50);
  disp->SetCDH(conf,70);
  //disp->SetCDCFrame(conf);

  disp->SetCDHHit(conf,cdsMan);
  //disp->SetCDCHitWire(conf,cdsMan);
  disp->SetCDCHit(conf,cdsMan);
  disp->SetCDCTrackHit(conf,trackMan);
  disp->SetCDCTrack(conf,trackMan);
  disp->SetCDCVertexPoint(conf,trackMan);
  //disp->PrintObj();

  disp->Draw("ogl");
}
