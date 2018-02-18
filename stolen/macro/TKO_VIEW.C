void TKO_VIEW()
{
  TString filename = "root/tkoview.root";
  TFile *f = new TFile( filename );
  int i = 0;
  TH1F *h1,*h2;
  int cr = 0;
  int sl = 1; 
  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 700 );
  c1->Divide(4,4);
  for( int ch = 0; ch <= 15;ch ++ ){i++; 
  //for( int ch = 16; ch <= 31;ch ++ ){i++; 
    c1->cd(i);
    h1 = (TH1F*)f->Get( Form("c%dn%da%d",cr,sl,ch) );
    //        gPad->SetLogy();
    //h1->fit("gaus",1000,2000);
    h1->Draw();
    //    h2 = (TH1F*)f->Get( Form("c%dn%da%d_tri",cr,sl,ch) );
    //h2->SetLineColor(2);
    //h2->Draw("same");
  
  }
 c1->Print("tmp.pdf");
 if(1){
   TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
   c2->Divide(2,2);
   i=0;
    for( int sl = 1; sl < 5;sl ++ ){i++; 
    //for( int ch = 16; ch <= 31;ch ++ ){i++; 
    c2->cd(i);
    h2 = (TH1F*)f->Get( Form("c%dn%d",cr,sl) );
    //     gPad->SetLogy();
    h2->Draw();
  } 
 c2->Print("tmp_2.pdf");

 }  
}
