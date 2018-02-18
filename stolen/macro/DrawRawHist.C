//for run47/49a/49c
#define STDCOUT 0
char* prefig="rawfig/";
char* predat="dat/";
char* folder_print="rawfig/";
gErrorIgnoreLevel=1001;
TText tex;
TString text;
const int nhodo=15;
char* hodoname[nhodo]={"BHD","T0",      "DEF","BPD",
		       "BVC", "CVC",     "PC", "NC",
		       "BD","LB",     "CDH","IH",
		       "WVC","HVC1","HVC2"
};
int nseg[nhodo]={20,5,	 8,70,
		 8,34,	 27,112,
		 5,2,	 36,24,
		 2,4,4
};
int divx[nhodo]={5,3,4,5,
		 4,6,6,4,
		 1,1,6,4,
		 1,2,2
};
int divy[nhodo]={4,2,2,7,
		 2,3,3,4,
		 1,1,3,6,
		 1,2,2
};
int maxpad[nhodo]={20,5,8,35,
		   8,18,18,16,
		   1,1,18,24,
		   1,4,4
};
#if 1
int adcmax[nhodo]={300,300,  1000,1000,
		   1500,1500,   1500,2000,
		   1000,1000,   1000,2000,
		   1000,1000,   500
};
#else
int adcmax[nhodo]={4100,4100,  4100,4100,
		   4100,4100,   4100,4100,
		   4100,4100,   4100,4100,
		   4100,4100,   410
};
#endif
double adccut[nhodo]={10,10,10,10,
		      10, 150,   150,150,
		      10,10,10,10,
		      10,10,10
};

int rebin[nhodo]={4,4,16,32,
		  8,8,8,8,
		  8,8,4,16,
		  4,4,4
};
int tdcrebin=8;
double acmin=0;
double acmax=16000;

const int ndc=7;
int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		CID_BLC2a,CID_BLC2b,
		CID_BPC,
		CID_FDC1,
		CID_CDC
};
char* dcname[ndc]={"BLC1a","BLC1b",
		   "BLC2a","BLC2b",
		   "BPC",
		   "FDC1",
		   "CDC"
};

int divdcx[ndc]={4,4,
		 4,4,
		 4,
		 3,
		 5
};
int divdcy[ndc]={2,2,
		 2,2,
		 2,
		 2,
		 3
};

int NumOfLayers[ndc]={8,8,8,8,8,6,15};
TCanvas* c1;
TH1F *h1;

void DrawHodoscope(TFile* f,int ihodo, TString psname, char* mrrun, int run,ofstream &ofs,bool CALC=true){
  char *name = hodoname[ihodo];
  std::cout << "Drawing Histgram for "<< name << std::endl;    
  c1= new TCanvas(name,name,800, 600 );
  TString psname2=Form("%s/Run%s/run%d/%s_adc_run%d.pdf",prefig,mrrun,run,name,run);
  c1->Print(psname2+"["); 
  c1->Divide(divx[ihodo],divy[ihodo]);
  int ipad=1;
  char* ud[2]={"U","D"};
  int nud=2;
  double max=0;
  if(name=="IH") nud=1;
  for(int iud=0;iud<nud;iud++){
    ipad=1;
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      //    c1->Clear("D");
      gPad->SetLogy();
      h1=(TH1F*)DrawHistR(f, Form("A%s%s%d",name,ud[iud],seg),rebin[ihodo],0,adcmax[ihodo],"hist" );
      if(h1)	max=h1->GetMaximum();
      h1=(TH1F*)DrawHistR(f, Form("AwoT%s%s%d",name,ud[iud],seg),rebin[ihodo],0,adcmax[ihodo],"histsame",4,1,0.5);
      double par[4];
      if( h1&&CALC&&SearchPedPeak(h1,par) ){
	DrawLineV(par[1],0.1,max,4,1,2);
	ofs<<"Ped: "<<name<<"  "<<seg<<"  "<<iud;
	for(int ipar=0;ipar<4;ipar++){
	  ofs<<"  "<<par[ipar];
	}
	ofs<<std::endl;
      }
      h1=(TH1F*)DrawHistR(f, Form("AwT%s%s%d",name,ud[iud],seg),rebin[ihodo],0,adcmax[ihodo],"histsame",2);
      if(h1){
	h1->GetXaxis()->SetRangeUser(par[1]+adccut[ihodo],4000);
	if( CALC&&SearchMIPPeak(h1,par) ){
	  DrawLineV(par[1],0.1,max,2,1,2);
	  ofs<<"MIP: "<<name<<"  "<<seg<<"  "<<iud;
	  for(int ipar=0;ipar<4;ipar++) ofs<<"  "<<par[ipar];
	  ofs<<std::endl;
	}  
      }
      if(seg==nseg[ihodo]){
	for(int ipad=seg+1;ipad<=maxpad[ihodo];ipad++){
	  c1->cd(ipad);
	  gPad->Clear();
	  c1->Update();
	}
	c1->cd();
	Print(c1,psname,text);
	Print(c1,psname2,text);
      }	  
      else if(ipad<maxpad[ihodo]) ipad++;
      else{
	c1->Update();
	c1->cd();
	Print(c1,psname,text);
	Print(c1,psname2,text);
	c1->Clear();
	c1->Divide(divx[ihodo],divy[ihodo]);
	ipad=1;
      }
    }    
  }
  c1->Print(psname2+"]");
  psname2=Form("%s/Run%s/run%d/%s_tdc_run%d.pdf",prefig,mrrun,run,name,run);
  c1->Print(psname2+"[");

  ipad=1;
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    gPad->SetLogy();    
    h1=(TH1F*)DrawHistR(f, Form("T%sU%d",name,seg),tdcrebin,10,4000,"hist");
    if(h1&&nud==2){
      h1=(TH1F*)DrawHistR(f, Form("T%sD%d",name,seg),tdcrebin,10,4000,"histsame",2 );
    }
    if(seg==nseg[ihodo]){
      for(int ipad=seg+1;ipad<=maxpad[ihodo];ipad++){
	c1->cd(ipad);
	gPad->Clear();
      }
      c1->Update();
      c1->cd();
      Print(c1,psname,text);
      Print(c1,psname2,text);
    }	  
    else if(ipad<maxpad[ihodo]) ipad++;
    else{
      c1->Update();
      c1->cd();
      Print(c1,psname,text);
      Print(c1,psname2,text);
      c1->Clear();
      c1->Divide(divx[ihodo],divy[ihodo]);
      ipad=1;
    }
  }
  c1->Print(psname2+"]");

  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  h1=(TH1F*)DrawHist(f, Form("HitPat%s",name) );
  if(h1)   h1->SetMinimum(0);
  c1->cd(2);
  gPad->SetLogy();
  h1=(TH1F*)DrawHist(f, Form("Mul%s",name) );
  c1->Update();
  c1->cd();
  psname2=Form("%s/Run%s/run%d/%s_pat_run%d.pdf",prefig,mrrun,run,name,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);

  if(ihodo==1){
    c1=new TCanvas(Form("%s_pile",name),Form("%s_pile",name),800,600);
    psname2=Form("%s/Run%s/run%d/%s_pile_run%d.pdf",prefig,mrrun,run,name,run);
    c1->Print(psname2+"[");
    
    c1->Divide(5,3);
    for(int iud=0;iud<2;iud++){
      for(int seg=1;seg<=5;seg++){
	c1->cd(seg);
	gPad->SetLogy();
	DrawHistR(f, Form("T%spre%s%d",name,ud[iud],seg),tdcrebin,10,4000,"hist" );
      }
      for(int seg=1;seg<=5;seg++){
	c1->cd(seg+5);
	gPad->SetLogy();
	DrawHistR(f, Form("T%s%s%d",name,ud[iud],seg),1,10,4000,"hist" );
      }
      for(int seg=1;seg<=5;seg++){
	c1->cd(seg+10);
	gPad->SetLogy();
	DrawHistR(f, Form("T%spost%s%d",name,ud[iud],seg),tdcrebin,10,4000,"hist" );
      }
      c1->cd();
      Print(c1,psname,text);
      Print(c1,psname2,text);
    }
    c1->Print(psname2+"]");
  }
}


void DrawAC(TFile *f,TString psname, char* mrrun, int run){
  c1 = new TCanvas( "AC", "AC", 800, 600 );
  std::cout << "Drawing Histgram for AC" << std::endl;    
  c1->Divide(2,3);
  for( int seg=1; seg<=4; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    DrawHistR(f,Form("AAC%d",seg),1,0,4100,"");
    DrawHistR(f,Form("AwTAC%d",seg),1,0,4100,"same",2);
  }  
  c1->cd(5);
  gPad->SetLogy();
  DrawHistR(f,Form("AACSUM"),1,acmin,acmax);
  DrawHistR(f,Form("AwTACSUM"),1,acmin,acmax,"same",2);

  c1->cd(6);  
  DrawHistR(f,Form("TAC1"),tdcrebin,10,4000);
  gPad->SetLogy();

  c1->Update();
  c1->cd();
  TString psname2=Form("%s/Run%s/run%d/AC_run%d.pdf",prefig,mrrun,run,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);
}

DrawDC(TFile* f,int idc, TString psname, char* mrrun, int run){
  char* name=dcname[idc];
  std::cout << "Drawing Histgram for "<< name << std::endl;    

  c1 = new TCanvas( name,name, 800, 600 );
  c1->Divide(divdcx[idc],divdcy[idc]);
  for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
    c1->cd(layer);
    DrawHistR(f,Form("T%s_%d",name,layer), tdcrebin, 10, 2000 );
  }
  c1->Update();
  c1->cd();
  TString psname2=Form("%s/Run%s/run%d/%s_tdc_run%d.pdf",prefig,mrrun,run,name,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);

  for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
    c1->cd(layer);
    h1=(TH1F*)DrawHist(f,Form("HitPat%s_%d",name,layer) );
    if(h1) h1->SetMinimum(0);
  }
  c1->Update();
  c1->cd();
  psname2=Form("%s/Run%s/run%d/%s_pat_run%d.pdf",prefig,mrrun,run,name,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);

  for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
    c1->cd(layer);
    gPad->Clear();
    gPad->SetLogy();
    DrawHist(f, Form("Mul%s_%d",name,layer) );
  }

  c1->Update();
  c1->cd();
  psname2=Form("%s/Run%s/run%d/%s_mul_run%d.pdf",prefig,mrrun,run,name,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);
}

DrawBDBAR(TFile *f,TString psname, char* mrrun, int run){
  c1 = new TCanvas( "Others", "Others", 800, 600 );
  int ipad=1;
  c1->Divide(5,3);
  for(int ihodo=8;ihodo<10;ihodo++){
    char* name=hodoname[ihodo];
    std::cout << "Drawing Histgram for "<< name << std::endl;    
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      gPad->SetLogy();
      DrawHistR(f, Form("A%sU%d",name,seg), rebin[ihodo], 0, adcmax[ihodo]);
      DrawHistR(f, Form("AwT%sU%d",name,seg), rebin[ihodo], 0, adcmax[ihodo], "same",2);
      ipad++;
      c1->cd(ipad);
      gPad->SetLogy();
      DrawHistR(f, Form("A%sD%d",name,seg), rebin[ihodo], 0, adcmax[ihodo]);
      DrawHistR(f, Form("AwT%sD%d",name,seg), rebin[ihodo], 0, adcmax[ihodo], "same",2);
      ipad++;
    }
  }
  c1->cd();
  TString psname2=Form("%s/Run%s/run%d/BDBAR_adc_run%d.pdf",prefig,mrrun,run,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);

  int ipad=1;
  c1->cd();
  c1->Clear();
  c1->Divide(5,3);
  for(int ihodo=8;ihodo<10;ihodo++){
    char* name=hodoname[ihodo];
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      gPad->SetLogy();
      h1=(TH1F*)DrawHistR(f, Form("T%sU%d",name,seg), tdcrebin, 10, 4000);
      TString option="";
      if(h1) option="same";
      DrawHistR(f, Form("T%sD%d",name,seg), tdcrebin, 10, 4000, "same",2);
      ipad++;
    }
  }
  for(int ihodo=8;ihodo<10;ihodo++){
    char* name=hodoname[ihodo];
    c1->cd(9+(ihodo-7)*2);
    DrawHist(f, Form("HitPat%s",name) );
    c1->cd(10+(ihodo-7)*2);
    gPad->SetLogy();
    DrawHist(f, Form("Mul%s",name) );
  }
  c1->cd();
  psname2=Form("%s/Run%s/run%d/BDBAR_tdc_run%d.pdf",prefig,mrrun,run,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);
  psname2=Form("%s/Run%s/run%d/BDBAR_pat_run%d.pdf",prefig,mrrun,run,run);
  Print(c1,psname2,text);
}
void DrawWVC(TFile *f,TString psname, char* mrrun, int run){
  c1 = new TCanvas( "WVC", "WVC", 800, 600 );
  int ipad=1;
  c1->Divide(2,2);
  int ihodo=12;
  char *name=hodoname[ihodo];
  std::cout << "Drawing Histogram for "<< name << std::endl;    
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    gPad->SetLogy();
    DrawHistR(f, Form("A%sU%d",name,seg), rebin[ihodo], 0, adcmax[ihodo]);
    DrawHistR(f, Form("AwT%sU%d",name,seg), rebin[ihodo], 0, adcmax[ihodo], "same",2);
    ipad++;
  }
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    gPad->SetLogy();
    DrawHistR(f, Form("A%sU%d",name,seg), tdcrebin, 10, 4000);
    ipad++;
  }
  c1->cd();
  TString psname2=Form("%s/Run%s/run%d/WVC_adc_run%d.pdf",prefig,mrrun,run,run);
  TString psname3=Form("%s/Run%s/run%d/WVC_tdc_run%d.pdf",prefig,mrrun,run,run);
  c1->Print(psname2);
  TString psname4=Form("%s/Run%s/run%d/WVC_pat_run%d.pdf",prefig,mrrun,run,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);
  Print(c1,psname3,text);
  Print(c1,psname4,text);
}


void PrintBeamline( TFile *f, TString psname, char* mrrun, int run,ofstream &ofs)
{
  //  for(int ihodo=0;ihodo<1;ihodo++){
  for(int ihodo=0;ihodo<4;ihodo++){
    DrawHodoscope(f,ihodo,psname,mrrun,run,ofs);
  }
  for(int ihodo=13;ihodo<15;ihodo++){
    DrawHodoscope(f,ihodo,psname,mrrun,run,ofs,false);
  }  
  DrawAC(f,psname,mrrun,run);
  for(int idc=0;idc<ndc-1;idc++){
    DrawDC(f,idc,psname,mrrun,run);
  }
}


void PrintCDS( TFile *f, TString psname, char* mrrun, int run,ofstream &ofs)
{
  for(int ihodo=10;ihodo<12;ihodo++){
    DrawHodoscope(f,ihodo,psname,mrrun,run,ofs,false);
  }
  DrawDC(f,ndc-1,psname,mrrun,run);
}

void PrintForwardTOF( TFile *f, TString psname, char* mrrun, int run,ofstream &ofs)
{
  for(int ihodo=4;ihodo<8;ihodo++){
    DrawHodoscope(f,ihodo,psname,mrrun,run,ofs);
  } 
  c1=new TCanvas("c2","c2",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  TH2F* h2;
  h2=DrawHist(f,"HitPatNC2D","colz");
  if(h2) h2->SetMinimum(0);
  c1->cd(2);
  h2=DrawHist(f,"HitPatNC2DwC","colz");
  if(h2) h2->SetMinimum(0);
  c1->cd(3);
  h2=DrawHist(f,"HitPatNC2DwoC","colz");
  if(h2) h2->SetMinimum(0);
  c1->cd(4);
  h2=DrawHist(f,"HitPatCharged","colz");
  if(h2) h2->SetMinimum(0);
  c1->cd();
  TString psname2=Form("%s/Run%s/run%d/NC_2d_run%d.pdf",prefig,mrrun,run,run);
  Print(c1,psname,text);
  Print(c1,psname2,text);
  
  DrawBDBAR(f,psname,mrrun,run);
  DrawWVC(f,psname,mrrun,run);
}

  
void PrintTKO(TString filename="tmp.root", double amin=0,double amax=1000)
{
  f = new TFile( filename );
  c1 = new TCanvas( "c1", "c1", 800, 600 );
  TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  c1->Divide(5,4);

  for( int c=3; c<=5; c++ ){
    for( int n=1; n<=20; n++ ){ 
      c1->cd(n);
      h1 = (TH1F*)f->Get( Form("TKOc%ds%d",c,n) );
      h1->SetAxisRange(amin,amax,"Y");
      h1->SetFillColor(2);
      h1->Draw();
      if(c==3 && n==20)   c1->Print("tko.ps(");
      else if(c==4 && n==20)   c1->Print("tko.ps");
      else if(c==5 && n==20)   c1->Print("tko.ps)");
    }
  }
}
void PrintStart(TFile *f, TString psname, char* mrrun, int run,ofstream &ofs)
{
  c1 = new TCanvas();
  c1->Divide(5,3);
  for(int seg=1;seg<=15;seg++){
    c1->cd(seg);
    h1 = (TH1F*)f->Get( Form("start%d",seg) );
    double par[4];
    if( SearchPedPeak(h1,par) ){
      h1->GetXaxis()->SetRangeUser(par[1]-10,par[1]+10);
      double max=h1->GetMaximum();
      h1->Draw();
      TLine line;
      line.SetLineColor(4);
      line.SetLineWidth(1);
      line.SetLineStyle(2);
      line.DrawLine(par[1],0.1,par[1],max);
      ofs<<"Start: "<<seg;
      for(int ipar=0;ipar<4;ipar++){
	ofs<<"  "<<par[ipar];
      }
      ofs<<std::endl;
    }    
  }
  c1->Print(psname);
}

void PrintAll(TString filename,char* mrrun,int runnum, double amin=0,double amax=1000)
{
  TFile *f=new TFile(filename);
  if(!f.IsOpen()){
    std::cout<<filename<<" does not exist !!!"<<std::endl;
    return;
  }
  text=filename;
#if STDCOUT
  ofstream &ofs=std::cout;
#else
  if(gSystem->Exec(Form("test -d %s",predat))){
    std::cout<<predat<<" does not exist !!!"<<std::endl;
    return;
  }
  if(gSystem->Exec(Form("test -d %s/Run%s",predat,mrrun))){
    gSystem->Exec(Form("mkdir %s/Run%s",predat,mrrun));
  }
  TString fname=Form("%s/Run%s/gain_run%d.dat",predat,mrrun,runnum);
  std::cout<<"trying to open "<<fname<<std::endl;
  ofstream ofs(fname);
  if(!ofs.is_open()){
    std::cout<<"cannot open "<<fname<<" !!!"<<std::endl; 
    return;
  }
  std::cout<<fname<<" opend."<<std::endl;
#endif
  if(gSystem->Exec(Form("test -d %s/Run%s",folder_print,mrrun))){
    gSystem->Exec(Form("mkdir %s/Run%s",folder_print,mrrun));
  }
  if(gSystem->Exec(Form("test -d %s/Run%s",prefig,mrrun))){
    gSystem->Exec(Form("mkdir %s/Run%s",prefig,mrrun));
  }
  if(gSystem->Exec(Form("test -d %s/Run%s/run%d",prefig,mrrun,runnum))){
    gSystem->Exec(Form("mkdir %s/Run%s/run%d",prefig,mrrun,runnum));
  }

  TString psname;
  c1=new TCanvas("c1","c1",800,600);
  psname=Form("%s/Run%s/beam_run%s_%d.pdf",folder_print,mrrun,mrrun,runnum);
  c1->Print(psname+"[");
  PrintBeamline(f,psname,mrrun,runnum,ofs);
  c1->Print(psname+"]");
  psname=Form("%s/Run%s/cds_run%s_%d.pdf",folder_print,mrrun,mrrun,runnum);
  c1->Print(psname+"[");
  PrintCDS(f,psname,mrrun,runnum,ofs);
  c1->Print(psname+"]");
  psname=Form("%s/Run%s/ftof_run%s_%d.pdf",folder_print,mrrun,mrrun,runnum);
  c1->Print(psname+"[");
  PrintForwardTOF(f,psname,mrrun,runnum,ofs);
  PrintStart(f,psname,mrrun,runnum,ofs);
  c1->Print(psname+"]");
  f->Close();
#if STDCOUT
  return;
#else
  ofs.close();
#endif
}

void DrawRawHist(char* mrrun,int runnum){
  //  PrintAll("/s/e17/hashimoto/k18ana/run45/evanaraw183.root"); 
  gSystem->Load("macro/HistTools_C.so");
  gSystem->Load("macro/FitTools_C.so");
  PrintAll(Form("root/run%s/evanaraw%d.root",mrrun,runnum),mrrun,runnum); 
  //  PrintAll("raw.root"); 
}
void DrawRawHistNew(){
  gSystem->Load("macro/HistTools_C.so");
  gSystem->Load("macro/FitTools_C.so");
  PrintAll(Form("root/run49c/test70.root"),"49c",999); 
}
