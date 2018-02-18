void CDS()
{
  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("root/1367.root");
  TTree *evtree = (TTree*)f->Get("CDSTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/730.conf");
  conf->Initialize();

  /*** assign output file ***/
  //TFile *of = new TFile("root/out.root","recreate");
  TFile *of = new TFile("root/out1367.root","recreate");

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  /*** declaration of classes ***/
  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader",  &head );
  
  of->cd();
  TH1F *h1;

  /*                */
  /* event analysis */
  /*                */
  new TH1F( "Pattern", "Pattern", 10, 0, 10 );
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    int nwire = (int)360/(conf->GetCDCWireMapManager()->dphi(layer));
    new TH1F( Form("CDCHitPattern%d",layer), Form("CDCHitPattern%d",layer), nwire, 0, nwire );
    new TH1F( Form("CDCTDC%d",layer), Form("CDCTDC%d",layer), 1500, 0, 1500 );
  }

  new TH1F( "CDHHitPattern", "CDHHitPattern", 40, 0, 40 );
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    new TH1F( Form("CDHADCU%d",seg), Form("CDHADCU%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHADCD%d",seg), Form("CDHADCD%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHTDCU%d",seg), Form("CDHTDCU%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHTDCD%d",seg), Form("CDHTDCD%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHADCUWT%d",seg), Form("CDHADCUWT%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHADCDWT%d",seg), Form("CDHADCDWT%d",seg), 2000, 0, 2000 );    
  }

  int nev = evtree->GetEntries();
  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);

    if( iev%5000==0 )
      std::cout << " Event : " << iev << std::endl;

    /*** if some parameters are newed, you can re-calc all quantities as below, ***/
    cdsMan->Calc(conf);

    /*** after here, we only fill histograms. ***/
    for( int i=0; i<20; i++ ){
      int val = head->pattern(i);
      if( 0<val ){ h1=(TH1F*)gFile->Get("Pattern"); h1->Fill(i); }
    }

    for( int i=0; i<tko->entries(); i++ ){
      TKOHit *hit = tko->hit(i);
      //std::cout << hit->cr() << " " << hit->sl() << " " << hit->ch() << " " << hit->data() << std::endl;
    }    
    
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      //std::cout << "Layer" << layer << "  #ev=" << cdsMan->nCDC(layer) << std::endl;
      for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	int wire = cdsMan->CDC(layer,i)->wire();
	int tdc = cdsMan->CDC(layer,i)->tdc();
	//std::cout << "   wire=" << wire << " tdc=" << tdc << std::endl;
	h1 = (TH1F*)gFile->Get( Form("CDCHitPattern%d",layer) ); h1->Fill(wire);
	h1 = (TH1F*)gFile->Get( Form("CDCTDC%d",layer) ); h1->Fill(tdc);
      }      
    }

    for( int i=0; i<cdsMan->nCDH(); i++ ){
      int seg  = cdsMan->CDH(i)->seg();
      int tdcu = cdsMan->CDH(i)->tdcu();
      int tdcd = cdsMan->CDH(i)->tdcd();
      int adcu = cdsMan->CDH(i)->adcu();
      int adcd = cdsMan->CDH(i)->adcd();
//       std::cout << " seg:" << seg << " tdcu:" << tdcu << " tdcd:" << tdcd
// 		<< " adcu:" << adcu << " adcd:" << adcd << std::endl;
      h1 = (TH1F*)gFile->Get( Form("CDHADCU%d",seg) ); h1->Fill(adcu);
      h1 = (TH1F*)gFile->Get( Form("CDHADCD%d",seg) ); h1->Fill(adcd);
      h1 = (TH1F*)gFile->Get( Form("CDHTDCU%d",seg) ); h1->Fill(tdcu);
      h1 = (TH1F*)gFile->Get( Form("CDHTDCD%d",seg) ); h1->Fill(tdcd);
      if( 0<tdcu && tdcu<4096 ){
	h1 = (TH1F*)gFile->Get( Form("CDHADCUWT%d",seg) ); h1->Fill(adcu);
      }
      if( 0<tdcd && tdcd<4096 ){
	h1 = (TH1F*)gFile->Get( Form("CDHADCDWT%d",seg) ); h1->Fill(adcd);
      }
      if( cdsMan->CDH(i)->CheckRange() ){
	h1 = (TH1F*)gFile->Get("CDHHitPattern"); h1->Fill(seg);
      }
    }

  }

  gFile->Write();
  gFile->Close();
}
