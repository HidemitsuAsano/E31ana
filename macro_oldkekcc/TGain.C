//memo 
//modified from original macro by H.Asano
//last update Jan. 24th 2018

void TGain(const char* filename="a.root", int cr=2,int sl=7)
{
  std::cout << "Macro TGain.C " << std::endl;
  std::cout << "input file name: " << filename << std::endl;
  std::cout << "Crate : " << cr << std::endl;
  std::cout << "Slot  : " << sl << std::endl;

  
  const int NTKOch = 31;
  /*** load library ***/
  gSystem->Load("lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile(filename);

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run78/analyzer.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TH1F *h1;
  TGraphErrors *gra;
  TLine line; line.SetLineColor(2);
  int npeak;
  double peakpos[100]; // position of peaks (rough, before fit)
  double fitmean[100]; // position of peaks (by fit)
  double fitsigma[100]; // sigma of peaks (by fit)
  double tdiff[100];
  for( int i=0; i<100; i++ ) tdiff[i] = 10.*i; // time calibrator was set to period=10nsec.

  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );
  c1->Divide(4,4);
  TCanvas *c2 = new TCanvas( "c2", "c2", 600, 600 );
  c2->Divide(4,4);

  for( int ch=0; ch<=NTKOch; ch++ ){
    npeak = 0;
    c1->cd( ch+1 );
    h1 = 0;
    h1 = (TH1F*)f->Get( Form( "TKOc%dn%da%d", cr,sl,ch ) );  
    if(h1==0){ std::cout << " !!! " << std::endl; continue; } // just check
    h1->Draw();
    int nbinx = h1->GetNbinsX();
    for( int i=5; i<nbinx-5; i++ ){ // search peak positions roughly (both edge of histogram are avoided)
      int content = h1->GetBinContent(i);
      if( 500<content ){ // at peak, content should be large enough.
	peakpos[npeak] = h1->GetBinCenter(i);
	npeak++;
	line.DrawLine( peakpos[npeak], 0, peakpos[npeak], h1->GetMaximum() );
	i += 20; // peak width is maybe smaller than 20 bins, and distance between peaks are maybe larger than 20bins.
      }
    }
    
    for( int i=0; i<npeak; i++ ){ // fit each peak
      h1->Fit( "gaus","0q","", peakpos[i]-10., peakpos[i]+10. );
      fitmean[i] = h1->GetFunction( "gaus" )->GetParameter(1);
      //fitsigma[i]= h1->GetFunction( "gaus" )->GetParameter(2);
      fitsigma[i]= h1->GetFunction( "gaus" )->GetParError(1);
    }
    
    c2->cd( ch+1 );
    gra = new TGraphErrors( npeak, tdiff, fitmean, 0, fitsigma );
    gra->Fit( "pol1" ); // linear fit
    gra->Draw("alpw");
    
    double gain = 1./gra->GetFunction("pol1")->GetParameter(1);
    conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, 0. ); // set a new parameter into map
  }
  
  ofstream ofs("tmp.param");
  //if( cr==Crate_PDC || cr==Crate_BLC || cr==Crate_BLHodo ){
  if( cr==Crate_BLC1 || cr==Crate_BLC2 || cr==Crate_BLHodo ){
    conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  }
  if( cr==Crate_CDC1 || cr==Crate_CDC2 || cr==Crate_CDC3 || cr==Crate_CDSHodo ){
    conf->GetGainMapManager()->PrintMapBL(ofs); // print map for CDS
  }
  //   conf->GetGainMapManager()->PrintMapHeader(ofs); // print header
  //   conf->GetGainMapManager()->PrintMap(0,ofs); // print map for crate 0
  //   conf->GetGainMapManager()->PrintMap(1,ofs); // print map for crate 1
  //   conf->GetGainMapManager()->PrintMap(2,ofs); // print map for crate 2
  ofs.close();

  //gFile->Write();
  //gFile->Close();
}
