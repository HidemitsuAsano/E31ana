void All( TString filename = "tmp.root")
{
  CDH(filename);
  CDC(filename);
  BHD(filename);
  PA(filename);
  T0(filename);
  TOF(filename);
  LC(filename);
  PDC(filename);
  BLC(filename);
}

void CDH( TString filename = "tmp.root", double amin=0, double amax=1000 )
{
  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
  c1->Divide(6,3);
  c2->Divide(6,3);
  for( int seg=1; seg<=36; seg++ ){
    if( seg<=18) c1->cd(seg);
    else c2->cd(seg-18);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ACDHU%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTCDHU%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  TCanvas *c4 = new TCanvas( "c4", "c4", 800, 600 );
  c3->Divide(6,3);
  c4->Divide(6,3);
  for( int seg=1; seg<=36; seg++ ){
    if( seg<=18) c3->cd(seg);
    else c4->cd(seg-18);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ACDHD%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTCDHD%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  TCanvas *c6 = new TCanvas( "c6", "c6", 800, 600 );
  c5->Divide(6,3);
  c6->Divide(6,3);
  for( int seg=1; seg<=36; seg++ ){
    if( seg<=18) c5->cd(seg);
    else c6->cd(seg-18);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TCDHU%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TCDHD%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(1,2);
  c7->cd(1);
  HitPatCDH->SetMinimum(0);
  HitPatCDH->Draw();
  c7->cd(2);
  gPad->SetLogy();
  MulCDH->Draw();

  c1->Print("cdh.ps(");
  c2->Print("cdh.ps");
  c3->Print("cdh.ps");
  c4->Print("cdh.ps");
  c5->Print("cdh.ps");
  c6->Print("cdh.ps");
  c7->Print("cdh.ps)");
}

void CDC( TString filename = "tmp.root" )
{
  TFile *f = new TFile( filename );

  TH1F *h1;
  TCanvas *c6 = new TCanvas( "c6", "c6", 800, 600 );
  c6->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c6->cd(layer);
    h1 = (TH1F*)f->Get( Form("TCDC%d",layer) );
    h1->Draw();
  }
  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c7->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatCDC%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c8 = new TCanvas( "c8", "c8", 800, 600 );
  c8->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c8->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulCDC%d",layer) );
    h1->Draw();
  }

  c6->Print("cdc.ps(");
  c7->Print("cdc.ps");
  c8->Print("cdc.ps)");

}

void BHD( TString filename = "tmp.root", double amin=0, double amax=500 )
{
  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  c1->Divide(5,4);
  for( int seg=1; seg<=20; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ABHDU%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTBHDU%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  c3->Divide(5,4);
  for( int seg=1; seg<=20; seg++ ){
    c3->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ABHDD%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTBHDD%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  c5->Divide(5,4);
  for( int seg=1; seg<=20; seg++ ){
    c5->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TBHDU%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TBHDD%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(1,2);
  c7->cd(1);
  HitPatBHD->SetMinimum(0);
  HitPatBHD->Draw();
  c7->cd(2);
  gPad->SetLogy();
  MulBHD->Draw();

  c1->Print("bhd.ps(");
  c3->Print("bhd.ps");
  c5->Print("bhd.ps");
  c7->Print("bhd.ps)");
}

void PA( TString filename = "tmp.root", double amin=0, double amax=500 )
{
  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  c1->Divide(4,2);
  for( int seg=1; seg<=8; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("APAU%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTPAU%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  c3->Divide(4,2);
  for( int seg=1; seg<=8; seg++ ){
    c3->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("APAD%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTPAD%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  c5->Divide(4,2);
  for( int seg=1; seg<=8; seg++ ){
    c5->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TPAU%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TPAD%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(1,2);
  c7->cd(1);
  HitPatPA->SetMinimum(0);
  HitPatPA->Draw();
  c7->cd(2);
  gPad->SetLogy();
  MulPA->Draw();

  c1->Print("pa.ps(");
  c3->Print("pa.ps");
  c5->Print("pa.ps");
  c7->Print("pa.ps)");
}

void T0( TString filename = "tmp.root", double amin=0, double amax=500 )
{
  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  c1->Divide(3,2);
  for( int seg=1; seg<=5; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("AT0U%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTT0U%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  c3->Divide(3,2);
  for( int seg=1; seg<=5; seg++ ){
    c3->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("AT0D%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTT0D%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  c5->Divide(3,2);
  for( int seg=1; seg<=5; seg++ ){
    c5->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TT0U%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TT0D%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(1,2);
  c7->cd(1);
  HitPatT0->SetMinimum(0);
  HitPatT0->Draw();
  c7->cd(2);
  gPad->SetLogy();
  MulT0->Draw();

  c1->Print("t0.ps(");
  c3->Print("t0.ps");
  c5->Print("t0.ps");
  c7->Print("t0.ps)");
}

void TOF( TString filename = "tmp.root", double amin=0, double amax=1000 )
{
  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  c1->Divide(8,4);
  for( int seg=1; seg<=32; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ATOFU%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTTOFU%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  c3->Divide(8,4);
  for( int seg=1; seg<=32; seg++ ){
    c3->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ATOFD%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTTOFD%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  c5->Divide(8,4);
  for( int seg=1; seg<=32; seg++ ){
    c5->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TTOFU%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TTOFD%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(1,2);
  c7->cd(1);
  HitPatTOF->SetMinimum(0);
  HitPatTOF->Draw();
  c7->cd(2);
  gPad->SetLogy();
  MulTOF->Draw();

  c1->Print("tof.ps(");
  c3->Print("tof.ps");
  c5->Print("tof.ps");
  c7->Print("tof.ps)");
}

void PDC( TString filename = "tmp.root" )
{
  TFile *f = new TFile( filename );

  TH1F *h1;
  TCanvas *c6 = new TCanvas( "c6", "c6", 800, 600 );
  c6->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c6->cd(layer);
    h1 = (TH1F*)f->Get( Form("TPDC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c6->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("TPDC2_%d",layer) );
    h1->Draw();
  }
  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c7->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatPDC1_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c7->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("HitPatPDC2_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c8 = new TCanvas( "c8", "c8", 800, 600 );
  c8->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c8->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulPDC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c8->cd(layer+8);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulPDC2_%d",layer) );
    h1->Draw();
  }

  c6->Print("pdc.ps(");
  c7->Print("pdc.ps");
  c8->Print("pdc.ps)");

}

void BLC( TString filename = "tmp.root" )
{
  TFile *f = new TFile( filename );

  TH1F *h1;
  TCanvas *c6 = new TCanvas( "c6", "c6", 800, 600 );
  c6->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c6->cd(layer);
    h1 = (TH1F*)f->Get( Form("TBLC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c6->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("TBLC2_%d",layer) );
    h1->Draw();
  }
  TCanvas *c7 = new TCanvas( "c7", "c7", 800, 600 );
  c7->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c7->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatBLC1_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c7->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("HitPatBLC2_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c8 = new TCanvas( "c8", "c8", 800, 600 );
  c8->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c8->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulBLC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c8->cd(layer+8);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulBLC2_%d",layer) );
    h1->Draw();
  }

  c6->Print("blc.ps(");
  c7->Print("blc.ps");
  c8->Print("blc.ps)");

}

void LC( TString filename = "tmp.root", double amin=0, double amax=1000 )
{
  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  c1->Divide(2,2);
  for( int seg=1; seg<=2; seg++ ){
    c1->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ALC2U%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTLC2U%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }
  for( int seg=1; seg<=2; seg++ ){
    c1->cd(seg+2);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ALC2D%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTLC2D%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

//   TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
//   c2->Divide(2,2);
//   for( int seg=1; seg<=2; seg++ ){
//     c2->cd(seg);
//     gPad->SetLogy();
//     h1 = (TH1F*)f->Get( Form("TLC2U%d",seg) );
//     h1->SetAxisRange(amin,amax,"X");
//     h1->Draw();
//   }
//   for( int seg=1; seg<=2; seg++ ){
//     c2->cd(seg+2);
//     gPad->SetLogy();
//     h1 = (TH1F*)f->Get( Form("TLC2D%d",seg) );
//     h1->SetAxisRange(amin,amax,"X");
//     h1->Draw();
//   }
  c1->Print("lc.ps");

}

void PrintCDC()
{
}
