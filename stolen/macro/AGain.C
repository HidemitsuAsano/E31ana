TString pdfname="again.pdf";
double edep[8]={0.587,1.957,1,1,1.957,5.87,5.87,10.};
int cid[8]={CID_BHD,CID_T0,CID_DEF,CID_BPD,CID_BVC,CID_CVC,CID_PC,CID_NC};
TCanvas *c1;
TH1F *h1;
TString psname="again.pdf";


void AGain(){
  gROOT->LoadMacro("~/k18ana/macro/FitTools.C");
  //  gROOT->LoadMacro("~/k18ana/macro/DrawRawHistNew.C");

  CVCGain("20131029_1",60);
  CVCGain("20131029_1",70);
  CVCGain("20131029_1",119);
  CVCGain("20131029_1",141);

}
void CVCGain(TString pre,int run){
  /*** assign input file & call tree ***/
  TFile *f = new TFile("~/k18ana/root/run49c/evanaraw43.root");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run49c/analyzer.conf",run);
  conf->Initialize();
  c1=new TCanvas("c1","c1",800,600);
  psname="cvcgain.pdf";
  c1->Print(psname+"[");
  calc(f,conf,5,1,1);
  calc(f,conf,4,1,1);
  c1->Print(psname+"]");
  ofstream ofs(Form("GainMapBL_%s_run%d.param",pre.Data(),run));
  conf->GetGainMapManager()->PrintMapBL(ofs);
  ofs.close();
}


void calc(TFile* f,ConfMan* conf, int ihodo, int iud,int seg,int xmax=500,int mipmin=0,double wid= 20){
  int cr,sl,ch;
  //  std::cout<< Form( "AwoT%s%s%d",name[ic].Data(),ud[iud].Data(),seg) << std::endl;
  std::cout << "Drawing Histgram for "<< hodoname[ihodo] << "  "<<adcmax[ihodo]<<std::endl;    
  c1= new TCanvas(hodoname[ihodo],hodoname[ihodo],800, 600 );
  c1->Divide(divx[ihodo],divy[ihodo]);
  int ipad=1;
  char* ud[2]={"U","D"};
  int nud=2;
  bool CALC=true;
  if(hodoname[ihodo]=="IH") nud=1;
  for(int iud=0;iud<2;iud++){
    ipad=1;
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      //    c1->Clear("D");
      gPad->SetLogy();
      
      h1 = (TH1F*)f->Get( Form("A%s%s%d",hodoname[ihodo],ud[iud],seg) );
      h1->Rebin(rebin[ihodo]);
      h1->SetAxisRange(0,adcmax[ihodo],"X");
      h1->Draw();
      double max=h1->GetMaximum();
      h1 = (TH1F*)f->Get( Form("AwoT%s%s%d",hodoname[ihodo],ud[iud],seg) );
      h1->Rebin(rebin[ihodo]);
      h1->SetLineColor(4);
      h1->SetLineWidth(0.5);
      h1->Draw("same");
      double par[4];
      double ped,mip;
      if( CALC&&SearchPedPeak(h1,par) ){
	TLine line;
	line.SetLineColor(4);
	line.SetLineWidth(1);
	line.SetLineStyle(2);
	line.DrawLine(par[1],0.1,par[1],max);
	std::cout<<"Ped: "<<hodoname[ihodo]<<"  "<<seg<<"  "<<iud;
	for(int ipar=0;ipar<4;ipar++){
	  std::cout<<"  "<<par[ipar];
	}
	std::cout<<std::endl;
	ped=par[1];
      }
      h1 = (TH1F*)f->Get( Form("AwT%s%s%d",hodoname[ihodo],ud[iud],seg) );
      h1->SetLineColor(2);
      h1->Rebin(rebin[ihodo]);
      h1->DrawCopy("same");
      
      //      h1->Rebin(5);
      h1->GetXaxis()->SetRangeUser(par[1]+adccut[ihodo],4000);
      if( CALC&&SearchMIPPeak(h1,par,100) ){
	TLine line;
	line.SetLineColor(2);
	line.SetLineWidth(1);
	line.SetLineStyle(2);
	line.DrawLine(par[1],0.1,par[1],max);
	std::cout<<"MIP: "<<hodoname[ihodo]<<"  "<<seg<<"  "<<iud;
	for(int ipar=0;ipar<4;ipar++){
	  std::cout<<"  "<<par[ipar];
	}
	std::cout<<std::endl;
	mip=par[1];
      }
      double gain=edep[ihodo]/(mip-ped);
      
      conf->GetCounterMapManager()->GetCNA(cid[ihodo],seg,0,iud,cr,sl,ch);
      conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, ped ); // set a new parameter into map
      // std::cout<<name[ic]<<"\t"<<ud[iud]<<seg<<"\t"<<ped<<"\t"<<mip<<"\t"<<edep[ic]<<"\t"<<gain<<std::endl;
      
      if(seg==nseg[ihodo]){
	for(int ipad=seg+1;ipad<=maxpad[ihodo];ipad++){
	  c1->cd(ipad);
	  gPad->Clear();
	  c1->Update();
	}
	c1->cd();
	tex.DrawTextNDC(0,0,text);
	c1->Print(psname);
      }	  
      else if(ipad<maxpad[ihodo]) ipad++;
      else{
	c1->Update();
	c1->cd();
	tex.DrawTextNDC(0,0,text);
	c1->Print(psname);
	c1->Clear();
	c1->Divide(divx[ihodo],divy[ihodo]);
	ipad=1;
      }
    }    
  }


  // h1=(TH1F*)f->Get( Form( "AwoT%s%s%d",name[ic].Data(),ud[iud].Data(),seg) );
  // //  h1=(TH1F*)f->Get( Form("%sADCwoT%d%s",name[ic].Data(),seg,ud[iud].Data()) );
  // if(h1->GetMaximum()<50){
  //   return;
  //   /*
  //   conf->GetCounterMapManager()->GetCNA(cid[ic],nSeg[ic]/2-1,0,iud,cr,sl,ch);
  //   conf->GetGainMapManager()->GetParam( cr,sl,ch, gain, mip );
  //   mip=ped+edep[ic]/gain;
  //   */
  // }
  // //  h1->Rebin(xmax/125);
  // h1->GetXaxis()->SetRangeUser(0,xmax);
  // double ped,mip,gain;

  // ped=h1->GetBinCenter(h1->GetMaximumBin());
  // //  std::cout<<ped<<std::endl;
  // h1->Fit("gaus","LI0q","",ped-10,ped+10);
  // ped=h1->GetFunction("gaus")->GetParameter(1);
  // h1->Draw();
  // h1->GetFunction("gaus")->Draw("same");
  // int ymax=  h1->GetMaximum();
  // h1=(TH1F*)f->Get( Form( "AwT%s%s%d",name[ic].Data(),ud[iud].Data(),seg) );
  // //  h1=(TH1F*)f->Get(Form( "%sADCwT%d%s",name[ic].Data(),seg,ud[iud].Data()) );
  // //  h1->Rebin(xmax/125*4);
  // h1->GetXaxis()->SetRangeUser(mipmin,xmax);
  // mip=h1->GetBinCenter(h1->GetMaximumBin());
  // TLine line2; line2.SetLineColor(4);
  // line2.DrawLine( mip, 1, mip, ymax );
  // h1->Fit("gaus","LI0q","",mip-wid,mip+wid);
  // mip=h1->GetFunction("gaus")->GetParameter(1);

  // h1->Draw("same");
  // h1->GetFunction("gaus")->Draw("same");

  // TLine line; line.SetLineColor(2);
  // line.SetLineWidth(2);
  // line.DrawLine( ped, 1, ped, ymax );
  // line.DrawLine( mip, 1, mip, ymax );

  // gain=edep[ic]/(mip-ped);

  // conf->GetCounterMapManager()->GetCNA(cid[ic],seg,0,iud,cr,sl,ch);
  // conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, ped ); // set a new parameter into map
  // std::cout<<name[ic]<<"\t"<<ud[iud]<<seg<<"\t"<<ped<<"\t"<<mip<<"\t"<<edep[ic]<<"\t"<<gain<<std::endl;
}

void AGainBPD(TFile *f, ConfMan* conf) //for bhd,t0,tofstop
{
  conf->Initialize();
  int ic=5;
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   800,800);
    c1->Divide(6,4);
    int ipad=1;
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(ipad);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg,2000,0,200);
      //void calc(TFile* f,ConfMan* conf, int ic, int iud,int seg,int xmax=500,int mipmin=0,double wid= 20){
      if(ipad==24){
	c1->Print(pdfname);
	ipad=1;
      }	else
	ipad++;
    }
  }
  
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  c1->Print(pdfname+"]");
}

void AGainBVC(TFile *f, ConfMan* conf) //for bhd,t0,tofstop
{
  conf->Initialize();
  int ic=6;
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   800,800);
    c1->Divide(5,2);
    int ipad=1;
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(ipad);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg,2000,0,200);
      //void calc(TFile* f,ConfMan* conf, int ic, int iud,int seg,int xmax=500,int mipmin=0,double wid= 20){
      if(ipad==8){
	c1->Print(pdfname);
	ipad=1;
      }	else
	ipad++;
    }
  }
  
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  c1->Print(pdfname+"]");
}

void AGainCVC(TFile *f, ConfMan* conf) //for bhd,t0,tofstop
{
  conf->Initialize();
  int ic=4;
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   800,800);
    c1->Divide(6,3);
    int ipad=1;
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(ipad);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg,2000,0,200);
      //void calc(TFile* f,ConfMan* conf, int ic, int iud,int seg,int xmax=500,int mipmin=0,double wid= 20){
      if(ipad==18||seg==nSeg[ic]){
	c1->Print(pdfname);
	ipad=1;
      }	else
	ipad++;
    }
  }
  
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  c1->Print(pdfname+"]");
}
void AGainBL(TFile *f, ConfMan* conf) //for bhd,t0,tofstop
{
  /*** load library ***/
  gSystem->Load("~/lib/libAll.so");

  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TH1F *h1;
  TGraphErrors *gra;

  double ped; 
  double mip; 
  double gain;

  int ic=0;
  TCanvas *c1=new TCanvas();
  c1->Print(pdfname+"[");
  // for(int iud=0;iud<2;iud++){
  //   c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
  // 		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
  // 		   800,800);
  //   c1->Divide(4,5);
  //   for(int seg=1;seg<=nSeg[ic];seg++){      
  //     c1->cd(seg);
  //     gPad->SetLogy();
  //     calc(f,conf,ic,iud,seg);
  //   }
  //   c1->Print(pdfname);
  // }

  ic=1;
  c1=new TCanvas(Form("%s ADC",name[ic].Data()),
		 Form("%s ADC",name[ic].Data()),
		 600,600);
  c1->Divide(3,4);
  for(int iud=0;iud<2;iud++){
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(seg+iud*6);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg);
    }
  }
  c1->Print(pdfname);

  ic=2;
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   600,600);
    c1->Divide(6,3);
    int ipad=1;
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(ipad);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg,1000,200,100);
      ipad++;
      if(ipad==19){
	c1->Print(pdfname);
	ipad=1;
      }
    }
    c1->Print(pdfname);
  }

  ic=3;
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   600,600);
    c1->Divide(4,4);
    int ipad=1;
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(ipad);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg,1000,200,100);
      ipad++;
      if(ipad==17){
	c1->Print(pdfname);
	ipad=1;
      }
    }
  }

  ic=4;
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s ADC",name[ic].Data(),ud[iud].Data()),
		   600,600);
    c1->Divide(8,4);
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(seg);
      gPad->SetLogy();
      calc(f,conf,ic,iud,seg,1000,200,100);
    }
    c1->Print(pdfname);
  }

  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
  c1->Print(pdfname+"]");
  //gFile->Write();
  //gFile->Close();
}


void CalcOffset(TFile *f,ConfMan *conf,int cr,int sl){
  double offs,gain,temp1;
  TH1F* h1;
  for( int ch=0; ch<=31; ch++ ){
    h1 = 0;
    h1 = (TH1F*)f->Get( Form( "TKOc%dn%da%d", cr,sl,ch ) );  
    if(h1==0){ std::cout << " !!! " << cr <<" "<<sl<<" "<<ch<<std::endl; continue; } // just check
    if(h1->GetEntries()<100) continue; 
    int max=h1->GetBinCenter(h1->GetMaximumBin());
    h1->Fit("gaus","0q","L",max-20,max+20);
    offs=h1->GetFunction("gaus")->GetParameter(1);
    conf->GetGainMapManager()->GetParam( cr,sl,ch, gain, temp1 );
    conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, offs ); // set a new parameter into map
    std::cout<<"cr,sl,ch,offset\t"<<cr<<"\t"<<sl<<"\t"<<ch<<"\t"<<offs<<std::endl;
  }  
}
void AOffset(){
  /*** load library ***/
  gSystem->Load("~/lib/libAll.so");

  /*** assign input file & call tree ***/
  TFile *f = new TFile("~/data/k18ana/run43/tko1485.root");
  if(!f) return;
  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run43/analyzer.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */
  int cr=2;
  for(int sl=14;sl<22;sl++)
    CalcOffset(f,conf,cr,sl);
  int cr=6;
  for(int sl=1;sl<4;sl++)
    CalcOffset(f,conf,cr,sl);
  int cr=7;
  for(int sl=16;sl<22;sl++)
    CalcOffset(f,conf,cr,sl);
  int cr=8;
  for(int sl=13;sl<18;sl++)
    CalcOffset(f,conf,cr,sl);

  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofstream ofs("tmp2.param");
  conf->GetGainMapManager()->PrintMapCDS(ofs); // print map for BeamLine
  ofs.close();
}
void AGainTOF() //for bhd,t0,tofstop
{
  /*** load library ***/
  gSystem->Load("~/lib/libAll.so");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/tofs/analyzer.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TH1F *h1;
  TGraphErrors *gra;
  TCanvas *c1;

  double ped; 
  double mip; 
  double gain;

  int runinit = 546;
  int seginit[8]={1,5,9,13,17,21,25,29};
  c1=new TCanvas();  
  c1->Print("again.pdf[");
  int ic=2;
  /*** assign input file & call tree ***/
  for(int irun=0;irun<8;irun++){
    TFile *f = new TFile(Form("~/data/k18ana/tofs/ncpc2%d.root",runinit+irun));
    c1=new TCanvas(Form("%s ADC run %d",name[ic].Data(),runinit+irun),
		   Form("%s ADC run %d",name[ic].Data(),runinit+irun)
		   ,1);
    c1->Divide(4,2);  
    for(int iud=0;iud<2;iud++){    
      for(int seg=1;seg<=4;seg++){      
	c1->cd(seg+4*iud);
	gPad->SetLogy();
	calc(f,conf,ic,iud,seginit[irun]+seg-1,1000);
      }
    }
    c1->Print("again.pdf");
  }
  c1->Print("again.pdf]");
  ofstream ofs("tmp.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
}
