void NCRaw( TString filename = "root/NCcalib/NCr00034.root", double amin=0, double amax=1000 )
{
  //Layer 
  int Layer[2]={1,2};
  int sseg[2][4]={1,2,3,4,1,2,3,4};

  int seg[2][4];
  for(int i=0;i<4;i++)
    {
      seg[0][i]=(Layer[0]-1)*16+sseg[0][i];
      seg[1][i]=(Layer[1]-1)*16+sseg[1][i];
    }
  TFile *f = new TFile( filename );

  TH1F *h1;
  TH2F *h2;

  TCanvas *c1 = new TCanvas( "c1", "c1", 800, 600 );
  TCanvas *c2 = new TCanvas( "c2", "c2", 800, 600 );
  c1->Divide(4,3);
  c2->Divide(4,3);
  for( int iset=1;iset<=2; iset++ ){
  for( int iseg=0; iseg<4; iseg++ ){
   if(iset==1) c1->cd(iseg+1);
   else if(iset==2) c2->cd(iseg+1);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ANCU%d",seg[iset-1][iseg]) );
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTNCU%d",seg[iset-1][iseg]) );
    h1->SetLineColor(2);
    h1->Draw("same");
    if(iset==1) c1->cd(iseg+5);
    else if(iset==2) c2->cd(iseg+5);
    gPad->SetLogy();
    h1 = (TH1F*)f->Get( Form("ANCD%d",seg[iset-1][iseg]) );
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("AwTNCD%d",seg[iset-1][iseg]) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }
  if(iset==1) c1->cd(9);
  else if(iset==2) c2->cd(9);
  gPad->SetLogy();
  h1 = (TH1F*)f->Get( Form("APCU%d",1+(iset-1)*2) );
  h1->SetAxisRange(0,2000,"X");
  h1->Draw();
  h1 = (TH1F*)f->Get( Form("AwTPCU%d",1+(iset-1)*2) );
  h1->SetLineColor(2);
  h1->Draw("same");
  if(iset==1) c1->cd(10);
  else if(iset==2) c2->cd(10);
  gPad->SetLogy();
  h1 = (TH1F*)f->Get( Form("APCU%d",2+(iset-1)*2) );
  h1->SetAxisRange(0,2000,"X");
  h1->Draw();
  h1 = (TH1F*)f->Get( Form("AwTPCU%d",2+(iset-1)*2) );
  h1->SetLineColor(2);
  h1->Draw("same");
  }

  TCanvas *c3 = new TCanvas( "c3", "c3", 800, 600 );
  c3->Divide(4,3);
  for( int iset=1;iset<=2; iset++ ){
  for( int iseg=0; iseg<4; iseg++ ){
   if(iset==1) c3->cd(iseg+1);
   else if(iset==2) c3->cd(iseg+5);
    h1 = (TH1F*)f->Get( Form("TNCU%d",seg[iset-1][iseg]) );
    h1->SetAxisRange(0,2000,"X");
    h1->Draw();
    h1 = (TH1F*)f->Get( Form("TNCD%d",seg[iset-1][iseg]) );
    h1->SetLineColor(2);
    h1->Draw("same");
  }
  if(iset==1) c3->cd(9);
  else if(iset==2) c3->cd(10);
  h1 = (TH1F*)f->Get( Form("TPCU%d",1+(iset-1)*2) );
  h1->SetAxisRange(0,2000,"X");
  h1->Draw();
  h1 = (TH1F*)f->Get( Form("TPCU%d",2+(iset-1)*2) );
  h1->SetLineColor(2);
  h1->Draw("same");
  if(iset==1) c3->cd(11);
  else if(iset==2) c3->cd(12);
  h1 = (TH1F*)f->Get( Form("TPCU%d",5+(iset-1)) );
  h1->SetAxisRange(0,2000,"X");
  h1->Draw();

  }



  TCanvas *c4 = new TCanvas( "c4", "c4", 800, 600 );
  TCanvas *c5 = new TCanvas( "c5", "c5", 800, 600 );
  c4->Divide(4,3);
  c5->Divide(4,3);
  for( int iset=1;iset<=2; iset++ ){
  for( int iseg=0; iseg<4; iseg++ ){
   if(iset==1) c4->cd(iseg+1);
   else if(iset==2) c5->cd(iseg+1);
    h2 = (TH2F*)f->Get( Form("ATNCU%d",seg[iset-1][iseg]) );
    //    h2->SetAxisRange(0,2000,"X");
    h2->Draw("colz");
    if(iset==1) c4->cd(iseg+5);
    else if(iset==2) c5->cd(iseg+5);
    h2 = (TH2F*)f->Get( Form("ATNCD%d",seg[iset-1][iseg]) );
    //    h2->SetAxisRange(0,2000,"X");
    h2->Draw("colz");
  }
  if(iset==1) c4->cd(9);
  else if(iset==2) c5->cd(9);

  h2 = (TH2F*)f->Get( Form("ATPCU%d",1+(iset-1)*2) );
  h2->Draw("colz");

  if(iset==1) c4->cd(10);
  else if(iset==2) c5->cd(10);
  h2 = (TH2F*)f->Get( Form("ATPCU%d",2+(iset-1)*2) );
  h2->Draw("colz");
  }


}
