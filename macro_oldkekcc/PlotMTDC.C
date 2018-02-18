{
  TFile *f=new TFile("/s/e17/hashimoto/k18ana/run47/mtdc52.root");
  TString pdfname="mtdc.pdf";

  TCanvas *c1;
  TH1F* h1;
  c1=new TCanvas();
  c1->Divide(4,4);
  c1->Print(pdfname+"[");
  int ipad=1;
  for(int ii=0;ii<2;ii++){
    for(int jj=0;jj<32;jj++){
      h1=(TH1F*)f->Get(Form("MTDC%dch%d",ii,jj));
      c1->cd(ipad);
      h1->Draw();
      ipad++;
      if(ipad>16){
	c1->Print(pdfname);
	ipad=1;
      }
    }
  }
  int ipad=1;
  for(int ii=0;ii<2;ii++){
    for(int jj=0;jj<32;jj++){
      h1=(TH1F*)f->Get(Form("nhitMTDC%dch%d",ii,jj));
      c1->cd(ipad);
      h1->Draw();
      ipad++;
      if(ipad>16){
	c1->Print(pdfname);
	ipad=1;
      }
    }
  }
  c1=new TCanvas();
  c1->Divide(3,2);
  ipad=1;
  TH2F* h2;
  for(int seg=1;seg<=5;seg++){
    h2=(TH2F*)f->Get(Form("TKO_MTDC_T0t%d",seg));
    c1->cd(ipad);
    h2->Draw("colz");
    ipad++;
  }
  c1->Print(pdfname);
  ipad=1;
  for(int seg=1;seg<=5;seg++){
    h2=(TH2F*)f->Get(Form("TKO_MTDC_T0b%d",seg));
    c1->cd(ipad);
    h2->Draw("colz");
    ipad++;
  }
  c1->Print(pdfname);
    
  c1->Print(pdfname+"]");
}
