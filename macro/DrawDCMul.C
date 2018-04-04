void DrawDCMul(const char* filename){
  TFile *f1 = new TFile(filename,"READ");
  if(!f1) return;
  
  char dcname[256][6]={"BLC1a","BLC1b","BLC2a","BLC2b","BPC","FDC1"};
  int nlayer[6]      ={      8,      8,      8,      8,    8,    6};
  gStyle->SetOptStat("em");
  TCanvas *c[6];
  for(int idc=0;idc<6;idc++){
    c[idc] = new TCanvas(Form("c%s",dcname[idc]),Form("c%s",dcname[idc]),800,800);
    c[idc]->Divide(4,2);

    for(int ilr=0;ilr<nlayer[idc];ilr++){
      char hname[256];
      sprintf(hname,"Mul%s_%d",dcname[idc],ilr+1);
      //std::cout << hname << std::endl;
      c[idc]->cd(ilr+1);
      TH1F *htemp = (TH1F*)f1->Get(hname);
      htemp->Draw();
    }
  }


  TCanvas *ccdc = new TCanvas("cCDC","cCDC",800,800);
  ccdc->Divide(4,4);
  for(int ilr=0;ilr<15;ilr++){
    char hname[256];
    sprintf(hname,"MulCDC_%d",ilr+1);
    ccdc->cd(ilr+1);
    TH1F *htemp = (TH1F*)f1->Get(hname);
    htemp->Draw();
  }


}
