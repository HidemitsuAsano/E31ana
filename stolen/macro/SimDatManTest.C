void SimDatManTest()
{
  gSystem->Load("./lib/libAll.so");

  TFile *f = new TFile("./root/1367.root");
  TTree *evtree = (TTree*)f->Get("CDSTree");

  ConfMan *conf = new ConfMan( "conf/analyzer.conf" );
  conf->Initialize();

  TKOHitCollection *tko = 0;
  CDSHitMan *cdsMan = 0;
  EventHeader *head = 0;
  //evtree->SetBranchAddress( "TKOHitCol", &tko );
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader", &head );

  SimDataMan *simMan = new SimDataMan(conf);

  std::cout << " start " << std::endl;
  int nev = evtree->GetEntries();
  //for( int iev=0; iev<nev; iev++ ){
  for( int iev=0; iev<3; iev++ ){
    evtree->GetEvent(iev);

    //if( iev%5000==0 )
    std::cout << std::endl << " Event : " << iev << std::endl;

    cdsMan->Calc(conf);

    // from data
    std::cout << " nCDH : " << cdsMan->nCDH() << std::endl;
    for( int i=0; i<cdsMan->nCDH(); i++ ){
      HodoscopeLikeHit *cdh = cdsMan->CDH(i);
      int seg = cdh->seg();
      double time = cdh->ctmean();
      double dene = cdh->emean();
      double hitpos = cdh->hitpos();
//       std::cout << " seg:" << seg
// 		<< " crau:" << cdh->crau() << " slau:" << cdh->slau() << " chau:" << cdh->chau()
// 		<< " crad:" << cdh->crad() << " slad:" << cdh->slad() << " chad:" << cdh->chad()
// 		<< " crtu:" << cdh->crtu() << " sltu:" << cdh->sltu() << " chtu:" << cdh->chtu()
// 		<< " crtd:" << cdh->crtd() << " sltd:" << cdh->sltd() << " chtd:" << cdh->chtd()
// 		<< std::endl;
//       std::cout << " seg:" << seg
// 		<< " hitpos:" << cdh->hitpos()
// 		<< " lv:" << cdh->lv()
// 		<< " ctmean:" << cdh->ctmean() << " ctsub:" << cdh->ctsub()
// 		<< " ctu:" << cdh->ctu() << " ctd:" << cdh->ctd()
// 		<< " eu:" << cdh->eu() << " ed:" << cdh->ed() << std::endl;
//       std::cout << " seg:" << seg << " ctmean:" << cdh->ctmean() << " dene:" << cdh->emean()
// 		<< " hitpos:" << cdh->hitpos() << std::endl;
//       std::cout << " seg:" << seg << " adcu:" << cdh->adcu() << " adcd:" << cdh->adcd()
// 		<< " tdcu:" << cdh->tdcu() << " tdcd:" << cdh->tdcd()
// 		<< std::endl;
      simMan->SetCDHHit( seg, time, dene, hitpos ); // <-- create new CDH in CDS
    }

    // from data
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      std::cout << " nCDC(layer" << layer << ") :" << cdsMan->nCDC(layer) << std::endl;
      for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	CDCHit *hit = cdsMan->CDC(layer,i);
	int wire = hit->wire();
	double dl = hit->dl();

	std::cout << " seg:" << seg << " tdc:" << hit->tdc() << std::endl;
	std::cout << " seg:" << seg << " dl:" << dl << " dt:" << hit->dt() << std::endl;

	simMan->SetCDCHit( layer, wire, dl, 0 );
      }
    }


    // created cds
    CDSHitMan *cdssim = simMan->GetCDSHit();
    cdssim->CheckContainerSize();
    std::cout << " nCDH(2):" << cdssim->nCDH() << std::endl;
    for( int i=0; i<cdssim->nCDH(); i++ ){
      HodoscopeLikeHit *cdh = cdssim->CDH(i);
      int seg = cdh->seg();
      double time = cdh->ctmean();
      double dene = cdh->emean();
      double hitpos = cdh->hitpos();
//       std::cout << " seg:" << seg
// 		<< " crau:" << cdh->crau() << " slau:" << cdh->slau() << " chau:" << cdh->chau()
// 		<< " crad:" << cdh->crad() << " slad:" << cdh->slad() << " chad:" << cdh->chad()
// 		<< " crtu:" << cdh->crtu() << " sltu:" << cdh->sltu() << " chtu:" << cdh->chtu()
// 		<< " crtd:" << cdh->crtd() << " sltd:" << cdh->sltd() << " chtd:" << cdh->chtd()
// 		<< std::endl;
//       std::cout << " seg:" << seg
// 		<< " hitpos:" << cdh->hitpos()
// 		<< " lv:" << cdh->lv()
// 		<< " ctmean:" << cdh->ctmean() << " ctsub:" << cdh->ctsub()
// 		<< " ctu:" << cdh->ctu() << " ctd:" << cdh->ctd()
// 		<< " eu:" << cdh->eu() << " ed:" << cdh->ed() << std::endl;
//       std::cout << " seg:" << seg << " ctmean:" << cdh->ctmean() << " dene:" << cdh->emean()
// 		<< " hitpos:" << cdh->hitpos() << std::endl;
//       std::cout << " seg:" << seg << " adcu:" << cdh->adcu() << " adcd:" << cdh->adcd()
// 		<< " tdcu:" << cdh->tdcu() << " tdcd:" << cdh->tdcd()
// 		<< std::endl;
    }

    // from data
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      std::cout << " nCDC(layer" << layer << ") :" << cdssim->nCDC(layer) << std::endl;
      for( int i=0; i<cdssim->nCDC(layer); i++ ){
	CDCHit *hit = cdssim->CDC(layer,i);
	int wire = hit->wire();
	double dl = hit->dl();

	std::cout << " seg:" << seg << " tdc:" << hit->tdc() << std::endl;
	std::cout << " seg:" << seg << " dl:" << dl << " dt:" << hit->dt() << std::endl;
      }
    }



    simMan->Clear();
  }

}
