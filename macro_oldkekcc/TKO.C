void TKO()
{
  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("tmp.root");
  TTree *evtree = (TTree*)f->Get("Tree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/analyzer.conf");
  conf->Initialize();

  /*** assign output file ***/
  TFile *of = new TFile("root/out.root","recreate");

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  
  of->cd();
  TH1F *h1;

  /*                */
  /* event analysis */
  /*                */
  for( int cr=11; cr<=14; cr++ ){
    for( int sl=1; sl<=22; sl++ ){
      for( int ch=0; ch<=31; ch++ ){
	new TH1F( Form( "c%dn%da%d", cr,sl,ch), Form( "c%dn%da%d", cr,sl,ch), 4096, 0, 4096 );
      }
    }
  }

  int nev = evtree->GetEntries();
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    for( int i=0; i<tko->entries(); i++ ){
      TKOHit *hit = tko->hit(i);
      h1 = 0;
      h1 = (TH1F*)gFile->Get( Form( "c%dn%da%d", hit->cr(),hit->sl(),hit->ch()) );
      if( h1!=0 )  h1->Fill( hit->data() );

      std::cout << hit->cr() << "  " << hit->sl() << "  " << hit->ch() << "  " << hit->data() << std::endl;

    }    
  }

  gFile->Write();
  gFile->Close();
}
