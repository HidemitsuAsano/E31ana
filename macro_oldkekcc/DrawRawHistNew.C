const int nhodo=15;
char* hodoname[nhodo]={"BHD","T0",
		       "DEF","BPD",
		       "BVC", "CVC",
		       "PC", "NC",
		       "Beamdump","Longbar",
		       "IH","CDH",
		       "WVC","HVC1","HVC2"
};
int nseg[nhodo]={20,5,
		 8,70,
		 8,34,
		 27,112,
		 5,2,
		 4,36,
		 2,4,4
};
int divx[nhodo]={5,3,4,5,4,
		 6,6,4,1,1,
		 4,6,
		 1,2,2
};
int divy[nhodo]={4,2,2,7,2,
		 3,3,4,1,1,
		 6,3,
		 1,2,2
};
int maxpad[nhodo]={20,5,8,35,8,
		   18,18,16,1,1,
		   24,18,1,4,4
};
int adcmax[nhodo]={500,1000,
		   1000,4000,
		   1000,2000,
		   2000,2000,
		   2000,2000,
		   1000,1000,
		   1000,1000,
		   1000
};

const int ndc=6;
int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		CID_BLC2a,CID_BLC2b,
		CID_BPC,
		CID_FDC1
};
char* dcname[ndc]={"BLC1a","BLC1b",
		   "BLC2a","BLC2b",
		   "BPC",
		   "FDC1"
};

//int NumOfLayers[ndc]={NumOfBLCLayers,NumOfBLCLayers,
//		      NumOfBLCLayers,NumOfBLCLayers,
//		      NumOfBPCLayers,
//		      NumOfFDC1Layers
//};

int NumOfLayers[ndc]={8,8,8,8,8,6};
TCanvas* c1;
TFile *f;
TH1F *h1;

void DrawHodoscope(int ihodo, TString psname){
  std::cout << "Drawing Histgram for "<< hodoname[ihodo] << std::endl;    
  c1= new TCanvas(hodoname[ihodo],hodoname[ihodo],800, 600 );
  c1->Divide(divx[ihodo],divy[ihodo]);
  int ipad=1;
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    //    c1->Clear("D");
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("A%sU%d",hodoname[ihodo],seg) );
    h1->SetAxisRange(0,adcmax[ihodo],"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwT%sU%d",hodoname[ihodo],seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
    if(seg==nseg[ihodo]){
      for(int ipad=seg+1;ipad<=maxpad[ihodo];ipad++){
	c1->cd(ipad);
	gPad->Clear();
	c1->Update();
      }
      c1->Print(psname);
    }	  
    else if(ipad<maxpad[ihodo]) ipad++;
    else{
      c1->Update();
      c1->Print(psname);
      c1->cd();
      c1->Clear();
      c1->Divide(divx[ihodo],divy[ihodo]);
      ipad=1;
    }
    
  }
  ipad=1;
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("A%sD%d",hodoname[ihodo],seg) );
    h1->SetAxisRange(0,adcmax[ihodo],"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwT%sD%d",hodoname[ihodo],seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
    if(seg==nseg[ihodo]){
      for(int ipad=seg+1;ipad<=maxpad[ihodo];ipad++){
	c1->cd(ipad);
	gPad->Clear();
	c1->Update();
      }
      c1->Print(psname);
    }	  
    else if(ipad<maxpad[ihodo]) ipad++;
    else{
      c1->Update();
      c1->Print(psname);
      c1->cd();
      c1->Clear();
      c1->Divide(divx[ihodo],divy[ihodo]);
      ipad=1;
    }
  }
  ipad=1;
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("T%sU%d",hodoname[ihodo],seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("T%sD%d",hodoname[ihodo],seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
    if(seg==nseg[ihodo]){
      for(int ipad=seg+1;ipad<=maxpad[ihodo];ipad++){
	c1->cd(ipad);
	gPad->Clear();
      }
      c1->Update();
      c1->Print(psname);
    }	  
    else if(ipad<maxpad[ihodo]) ipad++;
    else{
      c1->Update();
      c1->Print(psname);
      c1->cd();
      c1->Clear();
      c1->Divide(divx[ihodo],divy[ihodo]);
      ipad=1;
    }
  }
  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  h1 = (TH1F*)f->Get( Form("HitPat%s",hodoname[ihodo]) );
  h1->SetMinimum(0);
  h1->Draw();
  c1->cd(2);
  gPad->SetLogy();
  h1 = (TH1F*)f->Get( Form("Mul%s",hodoname[ihodo]) );
  h1->Draw();
  c1->Update();
  c1->Print(psname);
  if(ihodo==1){
    c1=new TCanvas("c1","c1",800,600);
    c1->Divide(5,3);
    for(int seg=1;seg<=5;seg++){
      c1->cd(seg);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%spreU%d",hodoname[ihodo],seg) );
      h1->Draw();
    }
    for(int seg=1;seg<=5;seg++){
      c1->cd(seg+5);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%sU%d",hodoname[ihodo],seg) );
      h1->Draw();
    }
    for(int seg=1;seg<=5;seg++){
      c1->cd(seg+10);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%spostU%d",hodoname[ihodo],seg) );
      h1->Draw();
    }
    c1->Print(psname);
    for(int seg=1;seg<=5;seg++){
      c1->cd(seg);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%spreD%d",hodoname[ihodo],seg) );
      h1->Draw();
    }
    for(int seg=1;seg<=5;seg++){
      c1->cd(seg+5);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%sD%d",hodoname[ihodo],seg) );
      h1->Draw();
    }
    for(int seg=1;seg<=5;seg++){
      c1->cd(seg+10);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%spostD%d",hodoname[ihodo],seg) );
      h1->Draw();
    }
    c1->Print(psname);
  }
}

void PrintBeamline( TString filename= "tmp.root", TString psname="beam.ps",
		    double amin=0, double amax=1000 )
{
  f = new TFile( filename );
  c1=new TCanvas("c1","c1",800,600);
  c1->Print(psname+"[");
  for(int ihodo=0;ihodo<4;ihodo++){
    DrawHodoscope(ihodo,psname);
  }
  for(int ihodo=13;ihodo<15;ihodo++){
    DrawHodoscope(ihodo,psname);
  }
  
  c1 = new TCanvas( "AC", "AC", 800, 600 );
  std::cout << "Drawing Histgram for AC" << std::endl;    
  c1->Divide(2,3);
  for( int seg=1; seg<=4; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("AAC%d",seg) );
    h1->SetAxisRange(0,500,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTAC%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }  
  c1->cd(5);
  gPad->SetLogy();
  
  h1 = (TH1F*)f->Get( Form("AACSUM") );
  h1->SetAxisRange(amin,amax,"X");
  h1->Draw();
  
  h1 = (TH1F*)f->Get( Form("AwTACSUM") );
  h1->SetLineColor(2);
  h1->Draw("same");

  c1->cd(6);  
  h1 = (TH1F*)f->Get( Form("TAC1") );
  h1->Draw();
  h1->SetLineColor(1);
  gPad->SetLogy();
  c1->Update();
  c1->Print(psname);

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Drawing Histgram for "<< dcname[idc] << std::endl;    
    c1 = new TCanvas( dcname[idc],dcname[idc], 800, 600 );
    c1->Divide(4,4);
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      c1->cd(layer);
      h1 = (TH1F*)f->Get( Form("T%s_%d",dcname[idc],layer) );
      h1->Draw();
    }
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      c1->cd(layer+8);
      h1 = (TH1F*)f->Get( Form("HitPat%s_%d",dcname[idc],layer) );
      if(h1->GetEntries()>10) h1->SetMinimum(0);
      h1->Draw();
    }
    c1->Update();
    c1->Print(psname);
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      c1->cd(layer+8);
      gPad->Clear();
      c1->cd(layer);
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("Mul%s_%d",dcname[idc],layer) );
      h1->Draw();
    }
      c1->Update();
    c1->Print(psname);
  }

  c1->Print(psname+"]");
}


void PrintCDS( TString filename = "tmp.root", TString psname="cds.ps", 
	       double amin=0, double amax=1000 )
{
  f = new TFile( filename );
  c1=new TCanvas("c1","c1",800,600);
  c1->Print(psname+"[");
  for(int ihodo=11;ihodo<12;ihodo++){
    DrawHodoscope(ihodo,psname);
  }

  c1 = new TCanvas( "IH", "IH", 800, 600 );
  int ipad=1;
  c1->Divide(6,4);
  int ihodo=10;
  std::cout << "Drawing Histgram for "<< hodoname[ihodo] << std::endl;    
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("A%sU%d",hodoname[ihodo],seg) );
    h1->SetAxisRange(0,adcmax[ihodo],"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwT%sU%d",hodoname[ihodo],seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
    ipad++;
  }
  ipad=1;
  for( int seg=1; seg<=nseg[ihodo]; seg++ ){
    c1->cd(ipad);
    c1->Clear("D");
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("T%sU%d",hodoname[ihodo],seg) );
    h1->SetAxisRange(10,4000,"X");
    h1->Draw();
    ipad++;
  }
  //    c1->Clear("D");
  c1 = new TCanvas( "IH", "IH", 800, 600 );
  c1->Divide(1,2);
  ipad++;
  c1->cd(ipad);
  h1=(TH1F*)f->Get( Form("HitPat%s",hodoname[ihodo]) );
  //    h1->SetMinimum(0);
  h1->Draw();
  c1->cd(ipad);
  ipad++;
  gPad->SetLogy();
  h1 = (TH1F*)f->Get( Form("Mul%s",hodoname[ihodo]) );
  h1->Draw();
  c1->Print(psname);

  c1 = new TCanvas( "cdc", "cdc", 800, 600 );
  c1->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c1->cd(layer);
    h1 = (TH1F*)f->Get( Form("TCDC%d",layer) );
    h1->Draw();
  }
  c1->Update();
  c1->Print(psname);
  for( int layer=1; layer<=15; layer++ ){
    c1->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatCDC%d",layer) );
    h1->SetFillColor(2);
    h1->SetMinimum(0);
    h1->Draw();
  }
  c1->Update();
  c1->Print(psname);
  for( int layer=1; layer<=15; layer++ ){
    c1->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulCDC%d",layer) );
    h1->Draw();
  }
  c1->Update();
  c1->Print(psname);
  c1->Print(psname+"]");
}



void PrintForwardTOF( TString filename = "tmp.root", TString psname="ftof.ps",
		      double amin=0, double amax=1000 )
{
  f = new TFile( filename );
  c1=new TCanvas("c1","c1",800,600);
  c1->Print(psname+"[");
  for(int ihodo=4;ihodo<8;ihodo++){
    DrawHodoscope(ihodo,psname);
  } 
  c1=new TCanvas("c2","c2",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  HitPatNC2D->SetMinimum(0);
  HitPatNC2D->Draw("colz");
  c1->cd(2);
  HitPatNC2DwC->SetMinimum(0);
  HitPatNC2DwC->Draw("colz");
  c1->cd(3);
  HitPatNC2DwoC->SetMinimum(0);
  HitPatNC2DwoC->Draw("colz");
  c1->cd(4);
  HitPatCharged->SetMinimum(0);
  HitPatCharged->Draw();

  c1->Print(psname);

  c1 = new TCanvas( "Others", "Others", 800, 600 );
  int ipad=1;
  c1->Divide(5,3);
  for(int ihodo=8;ihodo<10;ihodo++){
    std::cout << "Drawing Histgram for "<< hodoname[ihodo] << std::endl;    
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      //      c1->Clear("D");
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("A%sU%d",hodoname[ihodo],seg) );
      h1->SetAxisRange(0,adcmax[ihodo],"X");
      h1->Draw();
      h1 = (TH1F*)f->Get( Form("AwT%sU%d",hodoname[ihodo],seg) );
      h1->SetLineColor(2);
      h1->Draw("same");
      ipad++;
      if(hodoname[ihodo]=="IH") continue;
      c1->cd(ipad);
      //      c1->Clear("D");
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("A%sD%d",hodoname[ihodo],seg) );
      h1->SetAxisRange(0,adcmax[ihodo],"X");
      h1->Draw();
      h1 = (TH1F*)f->Get( Form("AwT%sD%d",hodoname[ihodo],seg) );
      h1->SetLineColor(2);
      h1->Draw("same");
      ipad++;
    }
  }
  c1->Print(psname);
  int ipad=1;
  c1->cd();
  c1->Clear();
  c1->Divide(5,3);
  for(int ihodo=8;ihodo<10;ihodo++){
    //    std::cout << "Drawing Histgram for "<< hodoname[ihodo] << std::endl;    
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      //      c1->Clear("D");
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%sU%d",hodoname[ihodo],seg) );
      h1->SetAxisRange(10,4000,"X");
      h1->Draw();
      h1 = (TH1F*)f->Get( Form("T%sD%d",hodoname[ihodo],seg) );
      h1->SetLineColor(2);
      h1->Draw("same");
      ipad++;
    }
  }
  for(int ihodo=8;ihodo<10;ihodo++){
    c1->cd(9+(ihodo-7)*2);
    //    c1->Clear("D");
    h1=(TH1F*)f->Get( Form("HitPat%s",hodoname[ihodo]) );
    //    h1->SetMinimum(0);
    h1->Draw();
    c1->cd(10+(ihodo-7)*2);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("Mul%s",hodoname[ihodo]) );
    h1->Draw();
  }
  c1->Print(psname);
  c1 = new TCanvas( "Others2", "Others2", 800, 600 );
  int ipad=1;
  c1->Divide(2,2);
  for(int ihodo=12;ihodo<13;ihodo++){
    std::cout << "Drawing Histgram for "<< hodoname[ihodo] << std::endl;    
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      //      c1->Clear("D");
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("A%sU%d",hodoname[ihodo],seg) );
      h1->SetAxisRange(0,adcmax[ihodo],"X");
      h1->Draw();
      h1 = (TH1F*)f->Get( Form("AwT%sU%d",hodoname[ihodo],seg) );
      h1->SetLineColor(2);
      h1->Draw("same");
      ipad++;
    }
  }
  for(int ihodo=12;ihodo<13;ihodo++){
    //    std::cout << "Drawing Histgram for "<< hodoname[ihodo] << std::endl;    
    for( int seg=1; seg<=nseg[ihodo]; seg++ ){
      c1->cd(ipad);
      //      c1->Clear("D");
      gPad->SetLogy();
      h1 = (TH1F*)f->Get( Form("T%sU%d",hodoname[ihodo],seg) );
      h1->SetAxisRange(10,4000,"X");
      h1->Draw();
      ipad++;
    }
  }
  c1->Print(psname);
  c1->Print(psname+"]");
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


void PrintAll(TString filename,int runnum, double amin=0,double amax=1000)
{
  PrintCDS(filename,Form("cds_run%d.pdf",runnum),amin,amax);
  PrintBeamline(filename,Form("beam_run%d.pdf",runnum),amin,amax);
  PrintForwardTOF(filename,Form("ftof_run%d.pdf",runnum),amin,amax);
}

void DrawRawHistNew(int runnum){
  //  PrintAll("/s/e17/hashimoto/k18ana/run45/evanaraw183.root"); 

  PrintAll(Form("/s/e17/hashimoto/k18ana/run47/evanaraw%d.root",runnum),runnum); 
  //  PrintAll("raw.root"); 
}
