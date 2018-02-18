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
  //evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader",  &head );
  
  int nev = evtree->GetEntries();
  if( nev<=sEvNum ) return;

  std::cout << " Event# : " << sEvNum << std::endl;
  evtree->GetEvent(sEvNum);

  cdsMan->Calc(conf);

  // first display
  TCanvas *c1 = new TCanvas( "c1", "c1", 200, 200, 500, 500 );
  Display3D *disp = new Display3D(c1);
  disp->SetCDSSolenoidBody(conf,180,40,40);
  disp->SetCDSSolenoidCap(conf,0xb,40);
  disp->SetTarget(conf,50);
  disp->SetCDH(conf,70);
  disp->SetCDCFrame(conf);

  disp->SetCDHHit(conf,cdsMan);
  disp->SetCDCHitWire(conf,cdsMan);
  disp->SetCDCHit(conf,cdsMan);
  //disp->PrintObj();
  disp->Draw("ogl");

  // second display
  TCanvas *c2 = new TCanvas( "c2", "c2", 200, 200, 500, 500 );
  Display3D *disp2 = new Display3D(c2);
  disp2->SetCDH(conf,70);
  disp2->SetCDHHit(conf,cdsMan);
  disp2->Draw("ogl");
}
