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
  double egatemax = 28.0;
  TF1 *slewfunc[36*2]; 
  for(int iseg=0;iseg<36;iseg++){
    c1->cd(iseg+1);
    dECDH_T03_CDH_PIM[iseg] = (TH2F*) file->Get(Form("dECDH_T03_CDHU%d_PIM",iseg+1));
    dECDH_T03_CDH_PIM[iseg]->Draw("colz");
    double meanx = dECDH_T03_CDH_PI[iseg]->GetMean(1);
    double meany = dECDH_T03_CDH_PI[iseg]->GetMean(2);
    //slewfunc[iseg] = new TF1( "slewfunc", "[0]/sqrt(x)", egatemin, egatemax );
    //slewfunc[iseg] = new TF1( "slewfunc", "[0]/sqrt(x)+[1]+[2]*x", egatemin, egatemax );
    slewfunc[iseg] = new TF1( "slewfunc", "[0]/sqrt(x)+[1]+[2]/x", egatemin, egatemax );
    //slewfunc[iseg]->FixParameter(0,0.0);
    //slewfunc[iseg]->SetParLimits(2,-0.01,0.01);
    dECDH_T03_CDH_PI[iseg]->ProfileX()->Fit("slewfunc","","",egatemin,egatemax);
    //slewfunc[iseg]->Draw("same");
    double p0 = -slewfunc[iseg]->GetParameter(0);
    double p1 = slewfunc[iseg]->GetParameter(1);
    double p2 = -slewfunc[iseg]->GetParameter(2);
    
    const int UP=0;
    const int DOWN=1;
    const int type=2;
    std::cout << "mean x y " << meanx << "  " << meany << std::endl;
    std::cout << "p1" << p1 << std::endl;
    std::cout << (1-p2)*meanx << std::endl;
    //double par[3]={p0,(1-p2)*meanx-7,p2};
    double par[3]={p0,(1-p2)*meanx-7,p2};
    conf->GetSlewingMapManager()->SetParam( CID_CDH, iseg+1, UP, 0, type, 3, par );
    conf->GetSlewingMapManager()->SetParam( CID_CDH, iseg+1, DOWN, 0, type, 3, par );
  }
  c1->Print("a.png");
  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_CDH,ofs);
  ofs.close();
}


//correct ctmean by meas_tof-calc_tof
//calc_tof is derived from CDC track momentum
void plotSlewCDH_2ndpros(const char *filename){
  TFile *file = new TFile(filename,"READ");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(6,6);
  TH2F* dECDH_T03_CDHU_PIM[36];
  
  gSystem->Load("./lib/libAll.so");
  
  ConfMan *conf = new ConfMan("/gpfs/home/had/hiasano/ana/k18ana/conf/Run78/analyzer_kwsk.conf");

  conf->SetRunNumber(800);
  conf->Initialize();

  for(int iseg=0;iseg<36;iseg++){
    TH2F* htemp = (TH2F*)file->Get(Form("CDH%d_mom_TOF_pi",iseg+1));
    double mean = htemp->GetMean();
    const int UP=0;
    const int DOWN=1;
    const int type=2;
    double par[3]={0,-mean,0};
    conf->GetSlewingMapManager()->SetParam( CID_CDH, iseg+1, UP, 3, type, 3, par );
    conf->GetSlewingMapManager()->SetParam( CID_CDH, iseg+1, DOWN, 3, type, 3, par );
  }

  ofstream ofs("tmp.dat");
  conf->GetSlewingMapManager()->PrintMap(CID_CDH,ofs);
  ofs.close();

}
