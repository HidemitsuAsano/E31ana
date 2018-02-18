void TOF()
{
  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("root/1164.root");
  TTree *evtree = (TTree*)f->Get("EventTree");
  TTree *scatree = (TTree*)f->Get("ScalerTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/730.conf");
  conf->Initialize();

  /*** assign output file ***/
  //TFile *of = new TFile("root/out.root","recreate");
  TFile *of = new TFile("root/out1164.root","recreate");

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  ScalerMan *scaMan = 0;
  //evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
  evtree->SetBranchAddress( "EventHeader",  &head );
  scatree->SetBranchAddress( "ScalerMan", &scaMan );
  
  of->cd();
  TH1F *h1;

  /*                 */
  /* scaler analysis */
  /*                 */
  new TH1F( "sca", "sca", 50, 0, 50 );
  int nev = scatree->GetEntries();
  std::cout << " # Scaler read : " << nev << std::endl;
  for( int iev=0; iev<nev; iev++ ){
    scatree->GetEvent(iev);
    for( int i=0; i<scaMan->nsca(); i++ ){
      int val = scaMan->sca(i)->val();
      TString name = scaMan->sca(i)->name();
      h1=(TH1F*)gFile->Get("sca"); h1->Fill(i,val);
    }
  }

  /*                */
  /* event analysis */
  /*                */
  new TH1F( "Pattern", "Pattern", 10, 0, 10 );
  for( int i=1; i<=NumOfTOFstopSegments; i++ ){
    new TH1F( Form( "TOF%dtdcu", i ), Form( "TOF%dtdcu", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dtdcd", i ), Form( "TOF%dtdcd", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dadcu", i ), Form( "TOF%dadcu", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dadcd", i ), Form( "TOF%dadcd", i ), 4000, 0, 4000 );
    new TH1F( Form( "TOF%dtimeu", i ), Form( "TOF%dtimeu", i ), 400, -20, 20 );
    new TH1F( Form( "TOF%dtimed", i ), Form( "TOF%dtimed", i ), 400, -20, 20 );
    new TH1F( Form( "TOF%dtmean", i ), Form( "TOF%dtmean", i ), 400, -20, 20 );
    new TH1F( Form( "TOF%dtsub",  i ), Form( "TOF%dtsub",  i ), 400, -20, 20 );
  }

  int nev = evtree->GetEntries();
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    blMan->Calc(conf);

    /*** after here, we only fill histograms. ***/
    for( int i=0; i<20; i++ ){
      int val = head->pattern(i);
      if( 0<val ){ h1=(TH1F*)gFile->Get("Pattern"); h1->Fill(i); }
    }

    of->cd();
    for( int i=0; i<blMan->nTOF(); i++ ){
      int seg  = blMan->TOF(i)->seg();
      int tdcu = blMan->TOF(i)->tdcu();
      int tdcd = blMan->TOF(i)->tdcd();
      int adcu = blMan->TOF(i)->adcu();
      int adcd = blMan->TOF(i)->adcd();
      h1=(TH1F*)gFile->Get( Form("TOF%dtdcu",seg ) ); h1->Fill(tdcu);
      h1=(TH1F*)gFile->Get( Form("TOF%dtdcd",seg ) ); h1->Fill(tdcd);
      h1=(TH1F*)gFile->Get( Form("TOF%dadcu",seg ) ); h1->Fill(adcu);
      h1=(TH1F*)gFile->Get( Form("TOF%dadcd",seg ) ); h1->Fill(adcd);

      if( !blMan->TOF(i)->CheckRange() ) continue;
      double tu = blMan->TOF(i)->tu();
      double td = blMan->TOF(i)->td(); 
      double tmean = blMan->TOF(i)->tmean(); 
      double tsub  = blMan->TOF(i)->tsub(); 
      h1=(TH1F*)gFile->Get( Form("TOF%dtimeu",seg ) ); h1->Fill(tu);
      h1=(TH1F*)gFile->Get( Form("TOF%dtimed",seg ) ); h1->Fill(td);
      h1=(TH1F*)gFile->Get( Form("TOF%dtmean",seg ) ); h1->Fill(tmean);
      h1=(TH1F*)gFile->Get( Form("TOF%dtsub", seg ) ); h1->Fill(tsub);
    }

  }

  gFile->Write();
  gFile->Close();
}
