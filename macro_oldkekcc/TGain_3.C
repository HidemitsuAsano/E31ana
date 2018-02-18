//from Kawasaki's macro

void TGain_3()
{


  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("src/lib/libAll.so");



  //SET//
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int BL =1;
  int CDS=0;

  int cr   =  1;       /// crate number



  std::string confname          = "conf/Run74/analyzer.conf";
          int runnum            = 1;
  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan(confname.c_str(), runnum);    /// param (Base)
  conf->Initialize();


  for(int k=1; k<=20; k++){
//if(k==-17 || k==-21)continue;
  if(k!=15)continue;

  int sl   =  k;       /// slot number

  int ch_s   =  0;     
  int ch_e   = 31;

  int d      =  6;     /// TCanvas divide (d x d)

  //broken ch
  int ch_b1= -16;
  int ch_b2= -12;
  int ch_b3= -13;

  int period = 40; ///[ns]

  int cut_s =   5; 
  int cut_e =   5; 


  std::string infilename        = "/group/had/knucl/e15/shinngo/Run74pre/tcal_0020.root";
  //std::string confname        = "conf/Run74/analyzer.conf";
  //        int runnum          = 1;
  std::string peakfit_name      = Form("TDCCalib/Run74pre/tcal_peakfit_cr%dsl%d_%d-%d.pdf",     cr,sl,ch_s,ch_e);       ///
  std::string linearfit_name    = Form("TDCCalib/Run74pre/tcal_linearfit_cr%dsl%d_%d-%d.pdf",   cr,sl,ch_s,ch_e);       ///
  std::string resi_name         = Form("TDCCalib/Run74pre/tcal_resi_cr%dsl%d_%d-%d.pdf",        cr,sl,ch_s,ch_e);       ///
  std::string peakfitsigma_name = Form("TDCCalib/Run74pre/tcal_peakfitsigma_cr%dsl%d_%d-%d.pdf",cr,sl,ch_s,ch_e);       ///

  std::string outfilename       = "tmp.param";                                 /// out param file name
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /*** assign input file & call tree ***/
  TFile *f = new TFile(infilename.c_str());  /// use root file

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TH1F *h1;
  TGraphErrors *gra;
  TGraphErrors *gra2;
  TGraphErrors *gra3;


  TLine line; line.SetLineColor(2);
  int npeak;
  double peakpos[100]; // position of peaks (rough, before fit)
  double fitmean[100]; // position of peaks (by fit)
  double fitsigma[100]; // sigma of peaks (by fit)
  double tdiff[100];
  double tdiff_err[100];

  double resi[100];
  for(int i=0; i<100; i++){
  resi[i]=0;
  }


  for( int i=0; i<100; i++ ) tdiff[i] = period*i;  /// time calibrator was set to period = 10nsec.

  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );
  c1->Divide(d,d);                           
  TCanvas *c2 = new TCanvas( "c2", "c2", 600, 600 );
  c2->Divide(d,d);                            
  TCanvas *c3 = new TCanvas( "c3", "c3", 600, 600 );
  c3->Divide(d,d);                            
  TCanvas *c4 = new TCanvas( "c4", "c4", 600, 600 );
  c4->Divide(d,d);                            


  


  for( int ch=ch_s; ch<=ch_e; ch++ ){                /// set ch Number
  if(!(ch==ch_b1 || ch==ch_b2 || ch==ch_b3)){
    npeak = 0;
    c1->cd( ch-ch_s+1 );                             //
    h1 = 0;
    h1 = (TH1F*)f->Get( Form( "TKOc%ds%da%d", cr,sl,ch ) );  
    if(h1==0){ std::cout << " !!! " << std::endl; continue; } // just check
    h1->GetXaxis()->SetRangeUser(0, 3000); 
    h1->Draw();
    int nbinx = h1->GetNbinsX();

    for( int i=cut_s; i<nbinx-cut_e; i++ ){ // search peak positions roughly (both edge of histogram are avoided)
      int content = h1->GetBinContent(i);
      if( 10<content ){ // at peak, content should be large enough.
	peakpos[npeak] = h1->GetBinCenter(i);
	
	line.DrawLine( peakpos[npeak], 0, peakpos[npeak], h1->GetMaximum() );
        npeak++;
        i += 20; // peak width is maybe smaller than 20 bins, and distance between peaks are maybe larger than 20bins.

	  }
    }
    
    for( int i=0; i<npeak; i++ ){ // fit each peak

    do    
      { h1->Fit( "gaus","0q","", peakpos[i]-10., peakpos[i]+10. );
      }
    
    while(h1->GetFunction( "gaus" )->GetParError(1)>10);
  
  
      fitmean[i] = h1->GetFunction( "gaus" )->GetParameter(1);
      //fitsigma[i]= h1->GetFunction( "gaus" )->GetParameter(2);
      fitsigma[i]= h1->GetFunction( "gaus" )->GetParError(1);

   if(fitmean[i]>4000)fitmean[i]=0;
    }
    


    gra = new TGraphErrors( npeak, tdiff, fitmean, 0, fitsigma);    
    gra->Fit("pol1"); // linear fit



    //Calc Residual************************************************
    double para[2];
    pol1->GetParameters(para);

    for( int i=0; i<npeak; i++ ){
    double y = para[1]*tdiff[i] +para[0];
    resi[i]  = y - fitmean[i];
    tdiff_err[i]=period*0.5;
    
    }


    //Remove large resi point in graph****************************

    int n=0;
    for( int i=0; i<npeak; i++ ){
      if(fabs(resi[i])>1){
	gra->RemovePoint(i-n);
      n++;
      }
    }


    c2->cd( ch-ch_s+1 );                           //
    //gra->Fit("pol1"); // linear fit
    gra->Draw("alpw");


    //Re-Calc Residual************************************************
    
    /*
    pol1->GetParameters(para);

    for( int i=0; i<npeak; i++ ){
    double y = para[1]*tdiff[i] +para[0];
    resi[i]  = y - fitmean[i];
    tdiff_err[i]=period*0.5;
    
    }
    */

    gra2 = new TGraphErrors( npeak, tdiff, resi, tdiff_err, 0);
    gra3 = new TGraphErrors( npeak, tdiff, fitsigma, tdiff_err, 0);

    
    n=0; 
    for( int i=0; i<npeak; i++ ){
      if(fabs(resi[i])>1){
      gra2->RemovePoint(i-n);
      gra3->RemovePoint(i-n);
      n++;
      }
     }
    


    c3->cd( ch-ch_s+1 );                           //
    gra2->Draw("AP");
    
    c4->cd( ch-ch_s+1 );                           //
    gra3->Draw("AP");


    gra->Fit("pol1"); // linear fit    
    double gain = 1./gra->GetFunction("pol1")->GetParameter(1);

    ////////////////////////////////////////////////////////////
    double gain_err = 1./gra->GetFunction("pol1")->GetParError(1);
    ////////////////////////////////////////////////////////////

    conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, 0. );   // set a new parameter into map
    


  }
  }
  


  c1->SaveAs(peakfit_name.c_str());          /// picture name 1
  c2->SaveAs(linearfit_name.c_str());        /// picture name 2
  c3->SaveAs(resi_name.c_str());             /// picture name 3
  c4->SaveAs(peakfitsigma_name.c_str());     /// picture name 4

  }

  ofstream ofs(outfilename.c_str());          

  if(BL)conf->GetGainMapManager()->PrintMapBL(ofs); 
  if(CDS)conf->GetGainMapManager()->PrintMapCDS(ofs); 

  ofs.close();



}
