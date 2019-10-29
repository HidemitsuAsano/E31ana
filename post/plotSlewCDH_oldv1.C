void plotSlewCDH(const char* filename){
  TFile *file = new TFile(filename,"READ");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(6,6);
  TH2F* dECDH_T03_CDHU_PIM[36];
  
  gSystem->Load("./lib/libAll.so");
  
  ConfMan *conf = new ConfMan("/gpfs/home/had/hiasano/ana/k18ana/conf/Run78/analyzer_kwsk.conf");
  
  conf->SetRunNumber(800);
  conf->Initialize();
  double egatemin = 0.5;
  double egatemax = 49.5;
  TF1 *slewfunc[36*2]; 
  for(int iseg=0;iseg<36;iseg++){
    c1->cd(iseg+1);
    dECDH_T03_CDHU_PIM[iseg] = (TH2F*) file->Get(Form("dECDH_T03_CDHU%d_PIM",iseg+1));
    dECDH_T03_CDHU_PIM[iseg]->Draw("colz");
    double meanx = dECDH_T03_CDHU_PIM[iseg]->GetMean(1);
    double meany = dECDH_T03_CDHU_PIM[iseg]->GetMean(2);
    //slewfunc[iseg] = new TF1( "slewfunc", "[0]/sqrt(x)", egatemin, egatemax );
    slewfunc[iseg] = new TF1( "slewfunc", "[0]/sqrt(x)+[1]+[2]*x", egatemin, egatemax );
    slewfunc[iseg]->FixParameter(0,0.0);
    dECDH_T03_CDHU_PIM[iseg]->Fit("slewfunc","q","",egatemin,egatemax);
    //slewfunc[iseg]->Draw("same");
    double p0 = -slewfunc[iseg]->GetParameter(0);
    double p1 = slewfunc[iseg]->GetParameter(1);
    double p2 = -slewfunc[iseg]->GetParameter(2);
    
    const int UP=0;
    const int type=2;
    std::cout << "mean x y " << meanx << "  " << meany << std::endl;
    std::cout << "p1" << p1 << std::endl;
    std::cout << (1-p2)*meanx << std::endl;
    double par[3]={p0,(1-p2)*meanx-12.63-0.5,p2};
    conf->GetSlewingMapManager()->SetParam( CID_CDH, iseg+1, UP, 1, type, 3, par );
  }
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(6,6);
  TH2F* dECDH_T03_CDHD_PIM[36];
  for(int iseg=0;iseg<36;iseg++){
    c2->cd(iseg+1);
    dECDH_T03_CDHD_PIM[iseg] = (TH2F*) file->Get(Form("dECDH_T03_CDHD%d_PIM",iseg+1));
    dECDH_T03_CDHD_PIM[iseg]->Draw("colz");
    double meanx = dECDH_T03_CDHD_PIM[iseg]->GetMean(1);
    double meany = dECDH_T03_CDHD_PIM[iseg]->GetMean(2);
    slewfunc[iseg+36] = new TF1( "slewfunc", "[0]/sqrt(x)+[1]+[2]*x", egatemin, egatemax );
    slewfunc[iseg+36]->FixParameter(0,0.0);
    dECDH_T03_CDHD_PIM[iseg]->Fit("slewfunc","","",egatemin,egatemax);
    //slewfunc[iseg+36]->Draw("same");
    double p0 = -slewfunc[iseg+36]->GetParameter(0);
    double p1 = slewfunc[iseg+36]->GetParameter(1);
    double p2 = -slewfunc[iseg+36]->GetParameter(2);
    
    const int DOWN=1;
    const int type=2;
    double par[3]={p0,(1-p2)*meanx-12.63-0.5,p2};
    conf->GetSlewingMapManager()->SetParam( CID_CDH, iseg+1, DOWN, 1, type, 3, par );
  }
  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_CDH,ofs);
  ofs.close();
}

/*
//correct ctmean by meas_tof-calc_tof
//calc_tof is derived from CDC track momentum
void plotSlewCDH_2ndpros(const *filename){
  TFile *file = new TFile(filename,"READ");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(6,6);
  TH2F* dECDH_T03_CDHU_PIM[36];
  
  gSystem->Load("./lib/libAll.so");
  
  ConfMan *conf = new ConfMan("/gpfs/home/had/hiasano/ana/k18ana/conf/Run78/analyzer_kwsk.conf");

  conf->SetRunNumber(800);
  conf->Initialize();

  for(int iseg=0;iseg<36;iseg++){
   TH2F*  



  }


}*/
