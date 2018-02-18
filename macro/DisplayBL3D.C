static int sEvNum = -1;

void disp( int evnum = -1 )
{
  gSystem->Load("./lib/libAll.so");
  sEvNum++;

  if( 0<evnum ) sEvNum = evnum;
 
  TFile *f = new TFile("./root/1164.root");
  TTree *evtree = (TTree*)f->Get("EventTree");
    
  ConfMan *conf = new ConfMan( "conf/analyzer.conf" );
  conf->Initialize();

  TKOHitCollection *tko = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
  evtree->SetBranchAddress( "EventHeader",  &head );
  
  int nev = evtree->GetEntries();
  if( nev<=sEvNum ) return;

  std::cout << " Event# : " << sEvNum << std::endl;
  evtree->GetEvent(sEvNum);

  blMan->Calc(conf);

  // set
  TCanvas *c1 = new TCanvas( "c1", "c1", 200, 200, 500, 500 );
  Display3D *disp1 = new Display3D(c1);
  disp1->SetCounterFrame( conf, CID_TOFstop, true );
  disp1->SetCounterHit( conf, blMan, CID_TOFstop, true );
//   disp1->SetCounter( conf, CID_PA );
//   disp1->SetCounter( conf, CID_BHD );
  //disp->PrintObj();
  disp1->Draw("ogl");

  TCanvas *c2 = new TCanvas( "c2", "c2", 200, 200, 500, 500 );
  Display3D *disp2 = new Display3D(c2);
  disp2->SetCounterFrame( conf, CID_T0, false );
  disp2->SetCounterFrame( conf, CID_B1, false );
  disp2->SetCounterFrame( conf, CID_B2, false );
  disp2->SetCounterHit( conf, blMan, CID_T0, false );
  disp2->SetCounterHit( conf, blMan, CID_B1, false );
  disp2->SetCounterHit( conf, blMan, CID_B2, false );
  disp2->Draw("ogl");

  TCanvas *c3 = new TCanvas( "c3", "c3", 200, 200, 500, 500 );
  Display3D *disp3 = new Display3D(c3);
  disp3->SetCounterFrame( conf, CID_PA, true );
  disp3->SetCounterHit( conf, blMan, CID_PA, true );
  disp3->Draw("ogl");

  TCanvas *c4 = new TCanvas( "c4", "c4", 200, 200, 500, 500 );
  Display3D *disp4 = new Display3D(c4);
  disp4->SetCounterFrame( conf, CID_BHD, true );
  disp4->SetCounterHit( conf, blMan, CID_BHD, true );
  disp4->Draw("ogl");
}
