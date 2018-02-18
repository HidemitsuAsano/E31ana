TString pdfname="tgain.pdf";
gErrorIgnoreLevel=1001;
void TGain(){
  //  TGaintofs();
  //  TGainRun46pre();
  //  TGainRun47pre();
  if(0){
    int confnum[]={60,65,71,76,81,86,91,100,106,110,115,120,125,130,135,139,145,151,155};
    int nconf=sizeof(confnum)/sizeof(int);
    for(int i=3;i<nconf;i++){
      TGainRun49cprechamber(Form("conf/Run49c/20130618/analyzer%d.conf",confnum[i]),"pass1","pass2");
    }
  }
  TGainRun49cpre();
  //  TGainRun49cprechamber("conf/Run49c/analyzer.conf","","");
  //  TGainRun47chamber();
  //  TGainRun46chamber();
  //TGainTDCTest();
    //  TGainRun35();
  //  TGainRun40();
  //  TGainRun43after();
}
void Calc(TFile *f,ConfMan *conf,int cr,int sl, double interval=10.,int init=0, bool POL2=false)
{
  TH1F *h1;
  TGraphErrors *gra;
  TGraphErrors *gra2;
  TLine line; line.SetLineColor(2);
  int npeak;
  double peakpos[100]; // position of peaks (rough, before fit)
  double fitmean[100]; // position of peaks (by fit)
  double fitsigma[100]; // sigma of peaks (by fit)
  double tdiff[100];
  double linearity[100];
  double linsigma[100];
  for( int i=0; i<100; i++ ) tdiff[i] = interval*i; // time calibrator was set to period=10nsec.

  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );
  c1->Divide(4,4);
  TCanvas *c2 = new TCanvas( "c2", "c2", 600, 600 );
  c2->Divide(4,4);
  TCanvas *c3 = new TCanvas( "c3", "c3", 600, 600 );
  c3->Divide(4,4);

  for( int tmpch=0; tmpch<=15; tmpch++ ){
    npeak = 0;
    int ch=tmpch+init;
    c1->cd( tmpch+1 );
    h1 = 0;
    h1 = (TH1F*)f->Get( Form( "TKOc%dn%da%d", cr,sl,ch) );  
    if(h1==0){ std::cout << " !!! " << cr <<" "<<sl<<" "<<ch<<std::endl; continue; } // just check
    if(h1->Integral(h1->FindBin(10),h1->FindBin(4000))<200) { std::cout << " !!!! " << cr <<" "<<sl<<" "<<ch<<std::endl; continue; } 
    h1->Draw();
    int nbinx = h1->GetNbinsX();
    for( int i=h1->FindBin(300); i<h1->FindBin(4000); i++ ){ // search peak positions roughly (both edge of histogram are avoided)
      int content = h1->GetBinContent(i);
      if( 200<content ){ // at peak, content should be large enough.
	peakpos[npeak] = h1->GetBinCenter(i);
	line.DrawLine( peakpos[npeak], 0, peakpos[npeak], h1->GetMaximum() );
	npeak++;
	i += 20; // peak width is maybe smaller than 20 bins, and distance between peaks are maybe larger than 20bins.
      }
    }
    if(npeak<5) { std::cout << " !!!!! " << cr <<" "<<sl<<" "<<ch<<std::endl; continue; }
    for( int i=0; i<npeak; i++ ){ // fit each peak
      h1->Fit( "gaus","L0q","", peakpos[i]-10., peakpos[i]+10. );
      fitmean[i] = h1->GetFunction( "gaus" )->GetParameter(1);
      //fitsigma[i]= h1->GetFunction( "gaus" )->GetParameter(2);
      fitsigma[i]= h1->GetFunction( "gaus" )->GetParError(1);
    }
    
    c2->cd( tmpch+1 );
    gra = new TGraphErrors( npeak-1 , fitmean, tdiff, fitsigma, 0 );
    gra->SetTitle(Form( "TKOc%dn%da%d", cr,sl,ch ));
    TString func="pol1";
    gra->Fit( "pol1" ,"q"); // linear fit
    if(POL2){
      func="pol2";      
      TF1 f2(func,func,0,4000);
      double par[3];
      gra->GetFunction("pol1")->GetParameters(par);
      par[2]=0.;
      f2.SetParameters(par);
      gra->Fit( &f2 ,"q"); // linear fit
    }
    gra->SetLineStyle(2);
    gra->SetMarkerStyle(20);
    gra->Draw("alpw");

    c3->cd( tmpch+1 );    
    TF1* f1=gra->GetFunction(func);
    f1->SetLineColor(2);
    f1->Draw("same");

    for( int i=0; i<npeak; i++ ){ // fit each peak
      linearity[i]=f1->Eval(fitmean[i])-tdiff[i];
      linsigma[i]=f1->Eval(fitmean[i]+fitsigma[i])-f1->Eval(fitmean[i]);
    }
    gra2 = new TGraphErrors( npeak-1 , fitmean, linearity, 0, linsigma );
    gra2->SetTitle(Form( "TKOc%dn%da%d", cr,sl,ch ));
    gra2->SetLineStyle(2);
    gra2->SetMarkerStyle(20);
    gra2->Draw("alpw");

    double p1 = f1->GetParameter(1);
    double p2 = f1->GetParameter(2);
    double offs,temp1;
    conf->GetGainMapManager()->GetParam( cr,sl,ch, temp1, offs );
    conf->GetGainMapManager()->SetParam( cr,sl,ch, offs, p1, p2 ); // set a new parameter into map
    std::cout<<"cr,sl,ch,gain before->after\t"<<cr<<"\t"<<sl<<"\t"<<ch<<"\t"<<p1<<" <- "<<temp1<<std::endl;
  }  
  c1->Print(pdfname);
  c2->Print(pdfname);
  c3->Print(pdfname);
}

void CalcSingle(TFile *f,ConfMan *conf,int cr,int sl,int ch, double interval=10.,int upto=15)
{
  TH1F *h1;
  TGraphErrors *gra;
  TLine line; line.SetLineColor(2);
  int npeak;
  double peakpos[100]; // position of peaks (rough, before fit)
  double fitmean[100]; // position of peaks (by fit)
  double fitsigma[100]; // sigma of peaks (by fit)
  double tdiff[100];
  for( int i=0; i<100; i++ ) tdiff[i] = interval*i; // time calibrator was set to period=10nsec.
  if(!c1) c1=new TCanvas();
  c1->Clear();
  c1->Divide(2,1);
  c1->cd(1);
  npeak = 0;
  h1 = 0;
  h1 = (TH1F*)f->Get( Form( "TKOc%dn%da%d", cr,sl,ch ) );  
  if(h1==0){ std::cout << " !!! " << cr <<" "<<sl<<" "<<ch<<std::endl;return; } // just check
  if(h1->Integral(h1->FindBin(10),h1->FindBin(4000)<100)) return; 
  h1->Draw();
  //  c1->Print("temp.pdf");
  int nbinx = h1->GetNbinsX();
  for( int i=h1->FindBin(300); i<h1->FindBin(4000); i++ ){ // search peak positions roughly (both edge of histogram are avoided)
    int content = h1->GetBinContent(i);
    if( 200<content ){ // at peak, content should be large enough.
      peakpos[npeak] = h1->GetBinCenter(i);
      line.DrawLine( peakpos[npeak], 0, peakpos[npeak], h1->GetMaximum() );
      npeak++;
      i += 20; // peak width is maybe smaller than 20 bins, and distance between peaks are maybe larger than 20bins.
    }
  }
  for( int i=0; i<npeak; i++ ){ // fit each peak
    h1->Fit( "gaus","L0q","", peakpos[i]-10., peakpos[i]+10. );
    fitmean[i] = h1->GetFunction( "gaus" )->GetParameter(1);
    //fitsigma[i]= h1->GetFunction( "gaus" )->GetParameter(2);
    fitsigma[i]= h1->GetFunction( "gaus" )->GetParError(1);
  }
  c1->cd(2);
  gra = new TGraphErrors( npeak-1 , tdiff, fitmean, 0, fitsigma );
  gra->Fit( "pol1" ,"q","",tdiff[0]+interval/2.,tdiff[npeak-1]-interval/2.); // linear fit
  gra->SetLineStyle(2);
  gra->SetMarkerStyle(20);
  gra->Draw("alpw");
  
  double gain = 1./gra->GetFunction("pol1")->GetParameter(1);
  double offs = gain * fitmean[0];
  //  conf->GetGainMapManager()->GetParam( cr,sl,ch, temp1, offs );
  //  conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, offs ); // set a new parameter into map
  conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, offs ); // set a new parameter into map
  std::cout<<"cr,sl,ch,gain\t"<<cr<<"\t"<<sl<<"\t"<<ch<<"\t"<<gain<<std::endl;
  c1->Print(pdfname);
}

void CalcChamber(TFile *f,ConfMan *conf,int cr,int sl, double interval=80,bool POL2=false)
{
  TH1F *h1;
  TGraphErrors *gra;
  TLine line; line.SetLineColor(2);
  int npeak;
  double peakpos[100]; // position of peaks (rough, before fit)
  double fitmean[100]; // position of peaks (by fit)
  double fitsigma[100]; // sigma of peaks (by fit)
  double tdiff[100];
  double linearity[100];
  double linsigma[100];
  for( int i=0; i<100; i++ ) tdiff[i] = interval*i; // time calibrator was set to period=10nsec.

  TCanvas *c1 = new TCanvas( "c1", "c1", 600, 600 );
  c1->Divide(4,4);
  TCanvas *c2 = new TCanvas( "c2", "c2", 600, 600 );
  c2->Divide(4,4);
  TCanvas *c3 = new TCanvas( "c3", "c3", 600, 600 );
  c3->Divide(4,4);
  int ipad=0;
  for( int ch=0; ch<=31; ch++ ){
    npeak = 0;
    ipad++;
    c1->cd( ipad );
    h1 = 0;
    h1 = (TH1F*)f->Get( Form( "TKOc%dn%da%d", cr,sl,ch ) );  
    if(h1==0){ std::cout << " !!! " << cr <<" "<<sl<<" "<<ch<<std::endl; continue; } // just check
    if(h1->GetEntries()<100) continue; 
    h1->Draw();
    int nbinx = h1->GetNbinsX();
    for( int i=h1->FindBin(300); i<h1->FindBin(4000); i++ ){ // search peak positions roughly (both edge of histogram are avoided)
      int content = h1->GetBinContent(i);
      if( 100<content ){ // at peak, content should be large enough.
	peakpos[npeak] = h1->GetBinCenter(i);
	line.DrawLine( peakpos[npeak], 0, peakpos[npeak], h1->GetMaximum() );
	npeak++;
	i += 20; // peak width is maybe smaller than 20 bins, and distance between peaks are maybe larger than 20bins.
      }
    }
    if(npeak<5){ std::cout << " !!! " << cr <<" "<<sl<<" "<<ch<<std::endl; continue; }
    for( int i=0; i<npeak; i++ ){ // fit each peak
      h1->Fit( "gaus","L0q","", peakpos[i]-10., peakpos[i]+10. );
      fitmean[i] = h1->GetFunction( "gaus" )->GetParameter(1);
      //fitsigma[i]= h1->GetFunction( "gaus" )->GetParameter(2);
      fitsigma[i]= h1->GetFunction( "gaus" )->GetParError(1);
    }
    h1->GetXaxis()->SetRangeUser(peakpos[0]-100,peakpos[npeak-1]+100);
    
    c2->cd( ipad );
    gra = new TGraphErrors( npeak-1, fitmean, tdiff, fitsigma, 0 );
    gra->SetTitle(Form( "TKOc%dn%da%d", cr,sl,ch ));

    TString func="pol1";
    gra->Fit( "pol1" ,"q"); // linear fit
    if(POL2){
      func="pol2";
      TF1 f2(func,func,0,4000);
      double par[3];
      gra->GetFunction("pol1")->GetParameters(par);
      par[2]=0.;
      f2.SetParameters(par);
      gra->Fit( &f2 ,"q"); // linear fit
    }
    gra->Draw("alpw");
    double offs,temp1;
    //    double gain = -1./func->GetParameter(1);

    c3->cd( ipad );    
    TF1* f1=gra->GetFunction(func);
    f1->SetLineColor(2);
    f1->Draw("same");

    for( int i=0; i<npeak; i++ ){ // fit each peak
      linearity[i]=f1->Eval(fitmean[i])-tdiff[i];
      linsigma[i]=f1->Eval(fitmean[i]+fitsigma[i])-f1->Eval(fitmean[i]);
    }

    gra2 = new TGraphErrors( npeak-1 , fitmean, linearity, 0, linsigma );
    gra2->SetTitle(Form( "TKOc%dn%da%d", cr,sl,ch ));
    gra2->SetLineStyle(2);
    gra2->SetMarkerStyle(20);
    gra2->Draw("alpw");


    double p1 =-f1->GetParameter(1);
    double p2 =-f1->GetParameter(2);
    if(cr==0||(cr==1&&sl>8&&sl<17)){
      p1*=-1;
      p2*=-1;
    }
    double offs,temp1;
    conf->GetGainMapManager()->GetParam( cr,sl,ch, temp1, offs );
    conf->GetGainMapManager()->SetParam( cr,sl,ch, offs, p1, p2 ); // set a new parameter into map
    std::cout<<"cr,sl,ch,gain before -> after"<<cr<<"\t"<<sl<<"\t"<<ch<<"\t"<<temp1<<" -> "<<p1<<std::endl;
    if(ch==15||ch==31){
      ipad=0;
      c1->Print(pdfname);
      c2->Print(pdfname);
      c3->Print(pdfname);
    }
  }  
}
TGainRun43before(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  /*** conf file for new parameters ***/
  //  ConfMan *conf = new ConfMan("conf/Apr2012/NCcalib20120413.conf");
  //  ConfMan *conf = new ConfMan("conf/Cosmic2008/anaT0.conf");
  ConfMan *conf = new ConfMan("conf/Run43/analyzer.conf");
  conf->Initialize();
  /*  
  int run[]={67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84};
  const nrun=18;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,
		     2,2,
		     2,2,
		     2,2,
		     2,2,
		     6,6,
		     6,6,
		     6,6,
		     7,7,
		     7,7,
		     7,7,
		     7,7,
		     7,7,
		     8,8,
		     8,8,
		     8,8,
		     8,8,
		     8,8};
  int sl[nrun][ncr]={1,2,
		     3,4,
		     7,9,
		     11,12,
		     13,14,
		     5,6,
		     7,8,
		     9,9,
		     1,2,
		     3,4,
		     5,6,
		     7,8,
		     9,10,
		     1,2,
		     3,4,
		     5,6,
		     7,8,
		     10,11,
		     };
  */
  /*
  const int nrun=3;
  int run[nrun]={157,158,159};
  //  static const int nrun = sizeof(run)/sizeof(int);
  const int ncr=2;
  int cr[nrun][ncr]={8,8,
		     8,8,
		     8,8};
  int sl[nrun][ncr]={1,2,
		     3,4,
		     5,5};
  */
  /*
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("./root/NCcalib/NCr000%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr]);
  }
  */
  //  TFile *f=new TFile("/s/e17/hashimoto/k18ana/cosmic2008/run0396.root");
  //  TFile *f=new TFile("/s/e17/hashimoto/k18ana/tko1460.root");
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run43//tko14%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr]);
  }
  int cr2[5]={0,1,3,4,5};
  int nsl[5]={16,16,20,20,20};
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run43//tko14%d.root",85));
    for(int icr=0;icr<5;icr++)
      for(int isl=1;isl<=nsl[icr];isl++)
	CalcChamber(f,conf,cr2[icr],isl);
    
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmp2.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();
}

TGainRun43after(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan("conf/Run43/analyzer.conf");
  conf->Initialize();

  int run[]={192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213};
  const nrun=22;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,     2,2,     2,2,     2,2,	     2,2,
		     7,7,     7,7,     7,7,     7,7,	     7,7,
		     8,8,     8,8,     8,8,     8,8,	     8,8,
		     6,6,     6,6,     6,6,
		     9,9,     9,9,     9,9,     9,9
  };
  int sl[nrun][ncr]={1,2,     3,4,     7,9,     11,12,	     13,14,
		     1,2,     3,4,     5,6,     7,8,	     9,10,
		     1,2,     3,4,     5,6,     7,8,	     10,11,
		     5,6,     7,8,     9,9,
		     1,2,     4,5,     6,9,     10,10
		     };
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run43//evanatko%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr]);
  }
  
  int cr2[5]={0,1,3,4,5};
  int nsl[5]={16,16,20,20,20};
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run43/evanatko%d.root",217));
  for(int icr=0;icr<5;icr++)
    for(int isl=1;isl<=nsl[icr];isl++)
      CalcChamber(f,conf,cr2[icr],isl,40.);
  f=new TFile(Form("/s/e17/hashimoto/k18ana/run43/evanatko%d.root",218));
  for(int isl=17;isl<=22;isl++)
    CalcChamber(f,conf,1,isl,40);
  
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmp2.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();
}

TGainRun45pre(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan("conf/Run45/analyzer.conf");
  conf->Initialize();

  int run[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  const nrun=19;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,     2,2,     2,2, 
		     6,6,     6,6,     6,6,     6,6,
		     2,2,     2,2,    
		     7,7,     7,7,     7,7,     7,7,     7,7,     
		     8,8,     8,8,     8,8,     8,8,     8,8
		     //		     9,9,     9,9,     9,9,     9,9
  };
  int sl[nrun][ncr]={2,3,     4,5,     7,7,
		     1,2,     5,6,     7,8,     9,9,
		     11,12,    13,14,
		     1,2,     4,5,     7,8,     9,11,	    13,14,
		     1,2,     4,5,     7,8,     10,11,      12,13
		     //		     1,2,     4,5,     6,9,     10,10
  };
  int interval[nrun]={10,10,10,
		      10,10,10,10,
		      20,20,
		      20,20,20,20,20,
		      20,20,20,20,20
  };

  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run45/evanatko%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr],interval[irun]);
  }
  
  int cr2[5]={0,1,3,4,5};
  int nsl[5]={16,22,20,20,20};
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run45/evanatko%d.root",20));
  //  for(int icr=0;icr<5;icr++)
  //    for(int isl=1;isl<=nsl[icr];isl++)
  //      CalcChamber(f,conf,cr2[icr],isl,40.);
  //  f=new TFile(Form("/s/e17/hashimoto/k18ana/run43/evanatko%d.root",218));
  //  for(int isl=17;isl<=22;isl++)
  //    CalcChamber(f,conf,1,isl,40);
  
  ofstream ofs("tmpbl.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmpcds.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();
}

TGainRun47pre(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan("conf/Run47/analyzer.conf");
  conf->Initialize();

  int run[]={1,2,3,4,5,6,7,
	     8,9,10,11,12,13,14,
	     15,20,17,18,19,
	     21,22,23,24,25,26,
	     27
	     //	     28,29,30,31,32,33
  };
  const nrun=26;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,     2,2,     2,2,     2,2,     2,2,
		     2,2,     2,2,  
		     6,6,     6,6,     6,6,     6,6,
		     6,6,     6,6,     6,6,  
		     7,7,     7,7,     7,7,     7,7,     7,7,     
		     8,8,     8,8,     8,8,     8,8,     8,8, 8,8,
		     6,6
		     //		     0,0,     0,0,     0,0,     0,0,     0,0, 0,0
		     //		     9,9,     9,9,     9,9,     9,9
  };
  int sl[nrun][ncr]={1,1,     2,3,     4,4,     7,7,     9,9,    
		     11,12,    13,14,
		     1,2,     3,3,      5,6,     7,8, 
		     9,9,     19,19,   20,20,
		     1,2,     4,5,     7,8,     9,11,	    13,14,
		     1,2,     4,5,     7,10,     11,12,      13,13,   8,8,
		     21,21		     
		     //		     1,2,     4,5,     6,9,     10,10
  };
  int interval[nrun]={40,20,20,10,40,20,20,
		      20,20,10,10,10,40,10,
		      20,20,20,20,20,
		      20,20,20,20,20,20,
		      40
  };
  
  c1=new TCanvas("gain","gain",600,600);
  c1->Print(pdfname+"[");
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run47pre/tko%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr],interval[irun]);
  }
  
  c1->Print(pdfname+"]");
  ofstream ofs("tmpbl.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmpcds.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();
}

TGainRun49cpre(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  // log book 13 p130- p121
  //  ConfMan *conf = new ConfMan("conf/Run49cpre/analyzer.conf");
  ConfMan *conf = new ConfMan("conf/Run49c/analyzer.conf",100);
  conf->Initialize();
  TCanvas *c1;

  if(0){
    int run[]={15,16,17,18,19,
	       20,21,22,23,24,
	       25,26,27,28,29,
	       30,31,32,33,34,
	       35,36,37,38,39,
	       40,41
    };
    const int nrun=sizeof(run)/sizeof(int);
    const int ncr=2;
    int cr[nrun][ncr]={2,2,     2,2,     2,6,     2,2,     2,2,
		       6,6,     2,6,     6,6,     6,6,     6,6, 
		       7,7,     7,7,     7,7,     7,7,     7,7,     
		       8,8,     8,8,     8,8,     8,8,     8,8,
		       0,0,     0,0,     0,0,     0,0,     0,0,
		       0,0,     6,6
		       //		     9,9,     9,9,     9,9,     9,9
    };
    int sl[nrun][ncr]={1,9,     2,3,     4,1,     11,12,     13,14,    
		       2,3,     7,5,     6,11,     8,9,      19,20, 
		       1,2,     4,5,     7,8,     9,11,	    13,14,
		       1,2,     4,5,     7,8,     10,11,      12,13,
		       1,1,     1,1,     2,2,     2,2,        3,3,
		       5,5,     19,21
		       //		     1,2,     4,5,     6,9,     10,10
    };
    int interval[nrun]={20,10,10,10,10,
			10,10,10,10,10,
			10,10,10,10,10,
			10,10,10,10,10,
			10,10,10,10,10,
			20,20
    };
    
    int init[nrun][ncr]={0,0, 0,0, 0,0, 0,0, 0,0,
			 0,0, 0,0, 0,0, 0,0, 0,0,
			 0,0, 0,0, 0,0, 0,0, 0,0,
			 0,0, 0,0, 0,0, 0,0, 0,0,
			 0,16, 32,48, 0,16, 32,48, 0,0,
			 0,16, 0,0		       
    };
    pdfname=Form("tgain_run49cpre_counter.pdf");
    c1=new TCanvas("gain","gain",600,600);
    c1->Print(pdfname+"[");
    for(int irun=0;irun<nrun;irun++){
      TFile *f=new TFile(Form("~/k18ana/root/run49cpre/tko%d.root",run[irun]));
      for(int icr=0;icr<ncr;icr++)
	Calc(f,conf,cr[irun][icr],sl[irun][icr],interval[irun],init[irun][icr],true);
    }
    c1=new TCanvas();
    c1->Print(pdfname+"]");
  }

  // FDC
  if(1){
    pdfname=Form("tgain_run49cpre_fdc1_pol1.pdf");
    int cr=9;
    int run[]={14,2,3,4,5,6,7,8,9,10,11,12};
    const int nrun=12;
    int sl[nrun]={9,10,12,13,14,
		  16,17,18,19,20,
		  21,22};
    c1=new TCanvas();
    c1->Print(pdfname+"[");
    for(int irun=0;irun<nrun;irun++){
      TFile *f=new TFile(Form("~/k18ana/root/run49cpre/tko%d.root",run[irun]));
      CalcChamber(f,conf,cr,sl[irun],40.);
    }
    c1=new TCanvas();
    c1->Print(pdfname+"]");
  }
  ofstream ofs("tmpbl.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmpcds.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();
  ofs2.close();
}

TGainRun49cprechamber(TString confname,TString before,TString after){
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan(confname.Data());
  conf->Initialize();

  pdfname=Form("tgain_run49cpre_chamber.pdf");
  int cr2[5]={0,1,3,4,5};
  int slini[5]={7,1,1,1,1};
  int slfin[5]={22,22,20,20,20};
  TFile *f=new TFile(Form("~/k18ana/root/run49cpre/tko%d.root",13));
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  for(int icr=0;icr<1;icr++)
    for(int isl=slini[icr];isl<=slfin[icr];isl++)
      CalcChamber(f,conf,cr2[icr],isl,40.,true);
  if(0){
    f=new TFile(Form("root/run49cpre/tko%d.root",42));
    for(int icr=2;icr<5;icr++)
      for(int isl=slini[icr];isl<=slfin[icr];isl++)
	CalcChamber(f,conf,cr2[icr],isl,40.,true);
  }
  std::cout<<"end"<<std::endl;
  c1=new TCanvas();
  c1->Print(pdfname+"]");
  std::cout<<"printed"<<std::endl;
  if(0){
    TString blname=conf->GetGainMapManager()->GetFileNameBL();
    blname.ReplaceAll(before,after);
    ofstream ofs(blname);
    conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
    TString cdsname=conf->GetGainMapManager()->GetFileNameCDS();
    cdsname.ReplaceAll(before,after);
    ofstream ofs2("tmpcds.param");
    conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
    ofs.close();
    ofs2.close();
  }
}

TGainRun47chamber(){
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan("conf/Run47/analyzer.conf");
  conf->Initialize();

  pdfname=Form("tgain_run47_chamber.pdf");
  int cr2[5]={0,1,3,4,5};
  int slini[5]={7,1,1,1,1};
  int slfin[5]={22,22,20,20,20};
  TFile *f=new TFile(Form("root/run47/tko%d.root",1));
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  for(int icr=0;icr<5;icr++)
    for(int isl=slini[icr];isl<=slfin[icr];isl++)
      CalcChamber(f,conf,cr2[icr],isl,40.);
  std::cout<<"end"<<std::endl;
  c1=new TCanvas();
  c1->Print(pdfname+"]");
  std::cout<<"printed"<<std::endl;
  ofstream ofs("tmpbl_run47_chamber.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmpcds_run47_chamber.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();

}

TGainRun46chamber(){
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan("conf/Run46/analyzer.conf");
  conf->Initialize();

  pdfname=Form("tgain_run46_chamber.pdf");
  int cr2[5]={0,1,3,4,5};
  int slini[5]={7,1,1,1,1};
  int slfin[5]={22,22,20,20,20};
  TFile *f=new TFile(Form("root/run46/tko%d.root",0));
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  for(int icr=0;icr<5;icr++)
    for(int isl=slini[icr];isl<=slfin[icr];isl++)
      CalcChamber(f,conf,cr2[icr],isl,40.);
  std::cout<<"end"<<std::endl;
  c1=new TCanvas();
  c1->Print(pdfname+"]");
  std::cout<<"printed"<<std::endl;
  ofstream ofs("tmpbl_run46_chamber.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmpcds_run46_chamber.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();

}
TGainRun46pre(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  ConfMan *conf = new ConfMan("conf/Run46/analyzer_upto6.conf");
  conf->Initialize();

  //  int run[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  int run[]={0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21};
  const nrun=21;
  const int ncr=2;
  int cr[nrun][ncr]={9,9,     9,9,     9,9,     9,9,     9,9,     
		     9,9,     9,9,     9,9,     9,9,     9,9,   
		     9,9,     9,9,
		     2,2,     2,2,     2,2,     2,2,     2,2,     2,2,
		     6,6,     6,6,     6,6
  };
  int sl[nrun][ncr]={10,10,   11,11,   12,12,   13,13,   14,14,
		     16,16,   17,17,   18,18,   19,19,   20,20,
		     21,21,   22,22,
		     2,3,     4,4,     5,5,     7,7,     11,12,   13,14,
		     1,2,     19,19,   20,20
  };
  int interval[nrun]={40,  40,  40,  40,  40,
		      40,  40,  40,  40,  40,
		      40,  40,
		      20,  20,   40,   10,   20,   20,
		      20,  40,   10
  };
  for(int irun=0;irun<12;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run46/evanatko%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr],interval[irun],31);
  }
  
  int cr2[5]={0,1,3,4,5};
  int slini[5]={7,1,1,1,1};
  int slfin[5]={22,22,20,20,20};
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run46/evanatko%d.root",22));
  //  for(int icr=0;icr<5;icr++)
  //    for(int isl=slini[icr];isl<=slfin[icr];isl++)
  //      CalcChamber(f,conf,cr2[icr],isl,40.);
  //  f=new TFile(Form("/s/e17/hashimoto/k18ana/run43/evanatko%d.root",218));
  //  for(int isl=17;isl<=22;isl++)
  //    CalcChamber(f,conf,1,isl,40);
  
  ofstream ofs("tmpbl.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs2("tmpcds.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs.close();
}

TGainRun40(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run40/analyzer.conf");
  conf->Initialize();
  int run[]={51,52,53,54,55,56};
  const nrun=6;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,
		     2,2,
		     2,2,
		     2,2,
		     2,2,
		     2,2
  };
  int sl[nrun][ncr]={1,2,
		     3,4,
		     5,7,
		     9,10,
		     11,12,
		     13,14};
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run40/evanatko14%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr]);
  }
  /*
  int cr2[2]={0,1};
  int nsl[2]={16,16};
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run40/evanatko14%d.root",50));
  for(int icr=0;icr<2;icr++)
    for(int isl=1;isl<=nsl[icr];isl++)
      CalcChamber(f,conf,cr2[icr],isl);
  */
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  /*
  ofstream ofs2("tmp2.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs2.close();
  */

}


TGainRun35(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run35/analyzer.conf");
  conf->Initialize();
  int run[]={1287,1288,1289,1290,1291,1292};
  const nrun=6;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,
		     2,2,
		     2,2,
		     2,2,
		     2,2,
		     2,2
  };
  int sl[nrun][ncr]={1,2,
		     3,4,
		     5,7,
		     9,10,
		     11,12,
		     13,14};
  int interval[nrun][ncr]={10,10,
			   10,10,
			   10,10,
			   10,10,
			   20,20,
			   20,20};
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run35/evanatko%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr],interval[irun][icr]);
  }
  /*
  int cr2[2]={0,1};
  int nsl[2]={16,16};
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run40/evanatko14%d.root",50));
  for(int icr=0;icr<2;icr++)
    for(int isl=1;isl<=nsl[icr];isl++)
      CalcChamber(f,conf,cr2[icr],isl);
  */
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  /*
  ofstream ofs2("tmp2.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs2); // print map for BeamLine
  ofs2.close();
  */

}
TGaintofs(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/tofs/analyzer.conf");
  conf->Initialize();
  int run[]={544,545};
  const nrun=2;
  const int ncr=2;
  int cr[nrun][ncr]={2,2,
		     2,2
  };
  int sl[nrun][ncr]={11,12,
		     13,14};
  int interval[nrun][ncr]={10,10,
			   10,10};
  for(int irun=0;irun<nrun;irun++){
    TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/tofs/evanatko%d.root",run[irun]));
    for(int icr=0;icr<ncr;icr++)
      Calc(f,conf,cr[irun][icr],sl[irun][icr],interval[irun][icr]);
  }
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
}

TGainTDCTest(){
  /*** load library ***/
  gSystem->Load("lib/libAll.so");
  
  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run47pre/analyzer.conf");
  conf->Initialize();
  int run=10;
  int interval=10;
  TFile *f=new TFile(Form("/s/e17/hashimoto/k18ana/run47pre/shift%d.root",run));
  int cr=2;
  int ch=15;
  c1=new TCanvas("gain","gain",800,400);
  c1->Print(pdfname+"[");
  CalcSingle(f,conf,cr,2,ch,interval);
  CalcSingle(f,conf,cr,3,ch,interval);
  CalcSingle(f,conf,cr,4,14,interval);
  CalcSingle(f,conf,cr,5,14,interval);

  for(int i=0;i<8;i++)
    CalcSingle(f,conf,cr,7,i,interval);
  CalcSingle(f,conf,cr,11,ch,interval);
  CalcSingle(f,conf,cr,12,ch,interval);
  CalcSingle(f,conf,cr,13,ch,interval);
  CalcSingle(f,conf,cr,14,ch,interval);
  
  cr=3;
  CalcSingle(f,conf,cr,1,ch,interval);
  CalcSingle(f,conf,cr,2,ch,interval);
  CalcSingle(f,conf,cr,4,ch,interval);
  CalcSingle(f,conf,cr,5,ch,interval);
  CalcSingle(f,conf,cr,7,ch,interval);
  CalcSingle(f,conf,cr,8,ch,interval);
  CalcSingle(f,conf,cr,9,ch,interval);
  CalcSingle(f,conf,cr,11,ch,interval);
  CalcSingle(f,conf,cr,13,ch,interval);
  CalcSingle(f,conf,cr,14,ch,interval);
  cr=0;
  CalcSingle(f,conf,cr,1,ch,interval);
  CalcSingle(f,conf,cr,2,ch,interval);
  CalcSingle(f,conf,cr,4,ch,interval);
  CalcSingle(f,conf,cr,5,ch,interval);
  CalcSingle(f,conf,cr,7,ch,interval);
  CalcSingle(f,conf,cr,8,ch,interval);
  c1->Print(pdfname+"]");

  ofstream ofs("tmpbl.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  ofs.open("tmpcds.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs); // print map for BeamLine
  ofs.close();
}
