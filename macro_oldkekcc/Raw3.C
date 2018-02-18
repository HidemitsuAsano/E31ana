void All( TString filename = "tmp.root", double amin=0, double amax=1000 )
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
    //    h1->Rebin();
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTCDHU%d",seg) );
    //   h1->Rebin(4);
    h1->SetLineColor(2);
    // h1->Draw("");
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
    //h1->Rebin();
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTCDHD%d",seg) );
    //h1->Rebin(4);
    h1->SetLineColor(2);
    h1->Draw("same");
    //h1->Draw("");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  TCanvas *c6 = new TCanvas( "c6", "c6", 800, 600 );
  c5->Divide(6,3);
  c6->Divide(6,3);
  for( int seg=1; seg<=36; seg++ ){
    if( seg<=18) c5->cd(seg);
    else c6->cd(seg-18);
    //  gPad->SetLogy();
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

  c1->Print("all.ps(");
  c2->Print("all.ps");
  c3->Print("all.ps");
  c4->Print("all.ps");
  c5->Print("all.ps");
  c6->Print("all.ps");
  c7->Print("all.ps");

  TCanvas *c8 = new TCanvas( "c8", "c8", 800, 600 );
  c8->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c8->cd(layer);
    h1 = (TH1F*)f->Get( Form("TCDC%d",layer) );
    h1->Draw();
  }
  TCanvas *c9 = new TCanvas( "c9", "c9", 800, 600 );
  c9->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c9->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatCDC%d",layer) );
    h1->SetFillColor(2);
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c10 = new TCanvas( "c10", "c10", 800, 600 );
  c10->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c10->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulCDC%d",layer) );
    h1->Draw();
  }

  c8->Print("all.ps");
  c9->Print("all.ps");
  c10->Print("all.ps");


  TCanvas *c11 = new TCanvas( "c11", "c11", 800, 600 );
  c11->Divide(5,4);
  for( int seg=1; seg<=20; seg++ ){
    c11->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ABHDU%d",seg) );
    h1->SetAxisRange(0,500,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTBHDU%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c12 = new TCanvas( "c12", "c12", 800, 600 );
  c12->Divide(5,4);
  for( int seg=1; seg<=20; seg++ ){
    c12->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ABHDD%d",seg) );
    h1->SetAxisRange(0,500,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTBHDD%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c13 = new TCanvas( "c13", "c13", 800, 600 );
  c13->Divide(5,4);
  for( int seg=1; seg<=20; seg++ ){
    c13->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TBHDU%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TBHDD%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c14 = new TCanvas( "c14", "c14", 800, 600 );
  c14->Divide(1,2);
  c14->cd(1);
  HitPatBHD->SetMinimum(0);
  HitPatBHD->Draw();
  c14->cd(2);
  gPad->SetLogy();
  MulBHD->Draw();

  c11->Print("all.ps");
  c12->Print("all.ps");
  c13->Print("all.ps");
  c14->Print("all.ps");

//   TCanvas *c15 = new TCanvas( "c15", "c15", 800, 600 );
//   c15->Divide(4,2);
//   for( int seg=1; seg<=8; seg++ ){
//     c15->cd(seg);
//     gPad->SetLogy();
//     h1 = (TH1F*)f->Get( Form("APAU%d",seg) );
//     h1->SetAxisRange(amin,amax,"X");
//     h1->Draw();
//     h1 = (TH1F*)f->Get( Form("AwTPAU%d",seg) );
//     h1->SetLineColor(2);
//     h1->Draw("same");
//   }

//   TCanvas *c16 = new TCanvas( "c16", "c16", 800, 600 );
//   c16->Divide(4,2);
//   for( int seg=1; seg<=8; seg++ ){
//     c16->cd(seg);
//     gPad->SetLogy();
//     h1 = (TH1F*)f->Get( Form("APAD%d",seg) );
//     h1->SetAxisRange(amin,amax,"X");
//     h1->Draw();
//     h1 = (TH1F*)f->Get( Form("AwTPAD%d",seg) );
//     h1->SetLineColor(2);
//     h1->Draw("same");
//   }

//   TCanvas *c17 = new TCanvas( "c17", "c17", 800, 600 );
//   c17->Divide(4,2);
//   for( int seg=1; seg<=8; seg++ ){
//     c17->cd(seg);
//     gPad->SetLogy();
//     h1 = (TH1F*)f->Get( Form("TPAU%d",seg) );
//     h1->SetAxisRange(10,4000);
//     h1->Draw();
//     h1 = (TH1F*)f->Get( Form("TPAD%d",seg) );
//     h1->SetLineColor(2);
//     h1->SetAxisRange(10,4000);
//     h1->Draw("same");
//   }

//   TCanvas *c18 = new TCanvas( "c18", "c18", 800, 600 );
//   c18->Divide(1,2);
//   c18->cd(1);
//   HitPatPA->SetMinimum(0);
//   HitPatPA->Draw();
//   c18->cd(2);
//   gPad->SetLogy();
//   MulPA->Draw();

//   c15->Print("all.ps");
//   c16->Print("all.ps");
//   c17->Print("all.ps");
//   c18->Print("all.ps");

  TCanvas *c19 = new TCanvas( "c19", "c19", 800, 600 );
  c19->Divide(3,2);
  for( int seg=1; seg<=5; seg++ ){
    c19->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("AT0U%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTT0U%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c20 = new TCanvas( "c20", "c20", 800, 600 );
  c20->Divide(3,2);
  for( int seg=1; seg<=5; seg++ ){
    c20->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("AT0D%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTT0D%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c21 = new TCanvas( "c21", "c21", 800, 600 );
  c21->Divide(3,2);
  for( int seg=1; seg<=5; seg++ ){
    c21->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TT0U%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TT0D%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c22 = new TCanvas( "c22", "c22", 800, 600 );
  c22->Divide(1,2);
  c22->cd(1);
  HitPatT0->SetMinimum(0);
  HitPatT0->Draw();
  c22->cd(2);
  gPad->SetLogy();
  MulT0->Draw();

  c19->Print("all.ps");
  c20->Print("all.ps");
  c21->Print("all.ps");
  c22->Print("all.ps");

  TCanvas *c23 = new TCanvas( "c23", "c23", 800, 600 );
  c23->Divide(8,4);
  for( int seg=1; seg<=32; seg++ ){
    c23->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ATOFU%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTTOFU%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c24 = new TCanvas( "c24", "c24", 800, 600 );
  c24->Divide(8,4);
  for( int seg=1; seg<=32; seg++ ){
    c24->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ATOFD%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTTOFD%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }

  TCanvas *c25 = new TCanvas( "c25", "c25", 800, 600 );
  c25->Divide(8,4);
  for( int seg=1; seg<=32; seg++ ){
    c25->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("TTOFU%d",seg) );
    h1->SetAxisRange(10,4000);
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TTOFD%d",seg) );
    h1->SetLineColor(2);
    h1->SetAxisRange(10,4000);
    h1->Draw("same");
  }

  TCanvas *c26 = new TCanvas( "c26", "c26", 800, 600 );
  c26->Divide(1,2);
  c26->cd(1);
  HitPatTOF->SetMinimum(0);
  HitPatTOF->Draw();
  c26->cd(2);
  gPad->SetLogy();
  MulTOF->Draw();

  c23->Print("all.ps");
  c24->Print("all.ps");
  c25->Print("all.ps");
  c26->Print("all.ps");

  TCanvas *c27 = new TCanvas( "c27", "c27", 800, 600 );
  c27->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c27->cd(layer);
    h1 = (TH1F*)f->Get( Form("TPDC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c27->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("TPDC2_%d",layer) );
    h1->Draw();
  }
  TCanvas *c28 = new TCanvas( "c28", "c28", 800, 600 );
  c28->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c28->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatPDC1_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c28->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("HitPatPDC2_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c29 = new TCanvas( "c29", "c29", 800, 600 );
  c29->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c29->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulPDC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c29->cd(layer+8);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulPDC2_%d",layer) );
    h1->Draw();
  }

  c27->Print("all.ps");
  c28->Print("all.ps");
  c29->Print("all.ps");

  TCanvas *c30 = new TCanvas( "c30", "c30", 800, 600 );
  c30->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c30->cd(layer);
    h1 = (TH1F*)f->Get( Form("TBLC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c30->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("TBLC2_%d",layer) );
    h1->Draw();
  }
  TCanvas *c31 = new TCanvas( "c31", "c31", 800, 600 );
  c31->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c31->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatBLC1_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c31->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("HitPatBLC2_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c32 = new TCanvas( "c32", "c32", 800, 600 );
  c32->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c32->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulBLC1_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c32->cd(layer+8);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulBLC2_%d",layer) );
    h1->Draw();
  }

  c30->Print("all.ps");
  c31->Print("all.ps");
  c32->Print("all.ps");

  //BPC
  TCanvas *c40 = new TCanvas( "c40", "c40", 800, 600 );
  c40->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c40->cd(layer);
    h1 = (TH1F*)f->Get( Form("TBPC_%d",layer) );
    h1->Draw();
  }
  for( int layer=1; layer<=8; layer++ ){
    c40->cd(layer+8);
    h1 = (TH1F*)f->Get( Form("HitPatBPC_%d",layer) );
    h1->SetMinimum(0);
    h1->Draw();
  }

  TCanvas *c41 = new TCanvas( "c41", "c41", 800, 600 );
  c41->Divide(4,4);
  for( int layer=1; layer<=8; layer++ ){
    c41->cd(layer);
    //    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulBPC_%d",layer) );
    h1->Draw();
  }
  c40->Print("all.ps");
  c41->Print("all.ps");


  TCanvas *c33 = new TCanvas( "c33", "c33", 800, 600 );
  c33->Divide(2,3);
  for( int seg=1; seg<=2; seg++ ){
    c33->cd(seg);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ALC2U%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTLC2U%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }
  for( int seg=1; seg<=2; seg++ ){
    c33->cd(seg+2);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ALC2D%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTLC2D%d",seg) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }
  
  for( int seg=1; seg<=2; seg++ ){
    c33->cd(seg+4);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ALC2SUM%d",seg) );
    h1->SetAxisRange(amin,amax,"X");
    h1->Draw();
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
  c33->Print("all.ps)");

}

void PrintCDS( TString filename = "tmp.root", double amin=0, double amax=1000 )
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
    //    h1->Rebin();
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTCDHU%d",seg) );
    //   h1->Rebin(4);
    h1->SetLineColor(2);
    // h1->Draw("");
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
    //h1->Rebin();
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTCDHD%d",seg) );
    //h1->Rebin(4);
    h1->SetLineColor(2);
    h1->Draw("same");
    //h1->Draw("");
  }

  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  TCanvas *c6 = new TCanvas( "c6", "c6", 800, 600 );
  c5->Divide(6,3);
  c6->Divide(6,3);
  for( int seg=1; seg<=36; seg++ ){
    if( seg<=18) c5->cd(seg);
    else c6->cd(seg-18);
    //  gPad->SetLogy();
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

  c1->Print("cds.ps(");
  c2->Print("cds.ps");
  c3->Print("cds.ps");
  c4->Print("cds.ps");
  c5->Print("cds.ps");
  c6->Print("cds.ps");
  c7->Print("cds.ps");

  TCanvas *c8 = new TCanvas( "c8", "c8", 800, 600 );
  c8->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c8->cd(layer);
    h1 = (TH1F*)f->Get( Form("TCDC%d",layer) );
    h1->Draw();
  }
  TCanvas *c9 = new TCanvas( "c9", "c9", 800, 600 );
  c9->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c9->cd(layer);
    h1 = (TH1F*)f->Get( Form("HitPatCDC%d",layer) );
    h1->SetFillColor(2);
    h1->SetMinimum(0);
    h1->Draw();
  }
  TCanvas *c10 = new TCanvas( "c10", "c10", 800, 600 );
  c10->Divide(5,3);
  for( int layer=1; layer<=15; layer++ ){
    c10->cd(layer);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("MulCDC%d",layer) );
    h1->Draw();
  }

  c8->Print("cds.ps");
  c9->Print("cds.ps");
  c10->Print("cds.ps");

}

void PrintTKO(TString filename="tmp.root", double amin=0,double amax=1000)
{

  TFile *f = new TFile( filename );

  TH1F *h1;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
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
