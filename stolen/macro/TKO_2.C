void TKO_2()
{
    /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  //TFile *f = new TFile("tmp.root");
  TFile *f = new TFile("tmp.root");
  TTree *evtree = (TTree*)f->Get("EventTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/analyzer.conf");
  conf->Initialize();

  /*** assign output file ***/
  TFile *of = new TFile("root/tkoview.root","recreate");

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;

  cout << "check1" << endl;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  cout << "check2" << endl;
  

  of->cd();
  TH1F *h1;
  /*                */
  /* event analysis */
  /*                */
  for( int cr=0; cr<1; cr++ ){
    for( int sl=0; sl<8; sl++ ){
      for( int ch=0; ch<=31; ch++ ){
	new TH1F( Form( "c%dn%da%d", cr, sl,ch), Form( "c%dn%da%d", cr ,sl,ch), 4096 ,0, 4095 );
      }
    }
  }

  int nev = evtree->GetEntries(); 
  std:: cout << "AllEvent" <<nev << endl;
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    for( int i=0; i<tko->entries(); i++ ){
      TKOHit *hit = tko->hit(i);
      h1 = 0;
      h1 = (TH1F*)gFile->Get( Form( "c%dn%da%d", hit->cr(),hit->sl(),hit->ch()) );
      if( h1!=0 )  h1->Fill( hit->data()-1 );
  
      std::cout << hit->cr() << "  " << hit->sl() << "  " << hit->ch() << "  " << hit->data() << std::endl;
    }
  }

  gFile->Write();
  gFile->Close();
}
