const double d_mass  = 1.87561;
const double K_mass    = 0.493677;
const double n_mass    = 0.939565;
const double pK = 1.00; //default value of simulation
const double EK = sqrt(K_mass*K_mass+pK*pK);
double func_EM(double *x, double *par);
double func_cq(double *x, double *par);
double func_C(double *x, double *par);
double func_S2(double *x, double *par);
double funcos(double q, double M);


const double costhetacutCMhi = 1.0; 
const double costhetacutCMlo = -1.0; 

void CalcLpimCSFromFit()
{
  TFile *flpim = TFile::Open("CSLpimFit.root","READ");

  const double br_s1385ToLambdapi = 0.87;
  const double br_s1385TopiSigma = 0.117;
  const double br_s1385TopiSigma_err = 0.015;
  const double IsospinCGFactor = 2.0;  

  //get Lpim
  TH2F* CS_lpim_sum = (TH2F*)flpim->Get("CS_sum");
  TH2F* CS_lpim_fit = (TH2F*)flpim->Get("Func");

  
  TH2F* CS_lpim_qcut[3];//iq
  for(int iq=0;iq<3;iq++){
    CS_lpim_qcut[iq] = (TH2F*) CS_lpim_fit->Clone(Form("CS_lpim_qcut%d",iq));
    CS_lpim_qcut[iq]->SetMinimum(0);
    CS_lpim_qcut[iq]->SetFillColor(0);
  }
  CS_lpim_qcut[0]->GetYaxis()->SetRangeUser(0,0.65);//total
  CS_lpim_qcut[1]->GetYaxis()->SetRangeUser(0,0.35);//
  CS_lpim_qcut[2]->GetYaxis()->SetRangeUser(0.35,0.65);
  double binwidthq = CS_lpim_qcut[0]->GetYaxis()->GetBinWidth(1)*1000.0; 
  std::cout << "binq width " << binwidthq  << std::endl;
  
  TH1D* CS_S1385_ToSp[3][3];//0:sysdown,1:center,2:sysup
  TH1D* CS_S1385_ToSm[3][3];//0:sysdown,1:center,2:sysup
  TH1D* CS_S1385_ToSpSm[3][3];//0:sysdown,1:center,2:sysup
 
  const int sysdef=1;
   
  //assume C.S. Sigma(1385)- ~ Sigma(1385)0
  for(int iq=0;iq<3;iq++){
    CS_lpim_qcut[iq]->Scale(br_s1385TopiSigma/2.0/br_s1385ToLambdapi*binwidthq/IsospinCGFactor);
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSp[iq][isys] = (TH1D*)CS_lpim_qcut[iq]->ProjectionX(Form("CS_S1385_ToSp%d_sys%d",iq,isys));
      CS_S1385_ToSm[iq][isys] = (TH1D*)CS_lpim_qcut[iq]->ProjectionX(Form("CS_S1385_ToSm%d_sys%d",iq,isys));
      CS_S1385_ToSpSm[iq][isys] = (TH1D*)CS_lpim_qcut[iq]->ProjectionX(Form("CS_S1385_ToSpSm%d_sys%d",iq,isys));
      CS_S1385_ToSpSm[iq][isys]->Scale(2.0);
    }
  }

  //qlow
  CS_S1385_ToSp[1][0]->Scale(0.0);//sys down
  CS_S1385_ToSp[1][2]->Scale(3);//sys up
  CS_S1385_ToSm[1][0]->Scale(0.0);//sys down
  CS_S1385_ToSm[1][2]->Scale(3);//sys up
  CS_S1385_ToSpSm[1][0]->Scale(0.0);//sys down
  CS_S1385_ToSpSm[1][2]->Scale(3);//sys up
  //qhi 
  CS_S1385_ToSp[2][0]->Scale(0.7);//sys down
  CS_S1385_ToSp[2][2]->Scale(1.3);//sys up
  CS_S1385_ToSm[2][0]->Scale(0.7);//sys down
  CS_S1385_ToSm[2][2]->Scale(1.3);//sys up
  CS_S1385_ToSpSm[2][0]->Scale(0.7);//sys down
  CS_S1385_ToSpSm[2][2]->Scale(1.3);//sys up
  
  //reset qall hist, because sys. err. are different btw qlow and qhi
  CS_S1385_ToSp[0][0]->Reset();
  CS_S1385_ToSp[0][2]->Reset();
  CS_S1385_ToSm[0][0]->Reset();
  CS_S1385_ToSm[0][2]->Reset();
  CS_S1385_ToSpSm[0][0]->Reset();
  CS_S1385_ToSpSm[0][2]->Reset();
  CS_S1385_ToSp[0][0]->Add(CS_S1385_ToSp[1][0]);
  CS_S1385_ToSp[0][0]->Add(CS_S1385_ToSp[2][0]);
  CS_S1385_ToSm[0][0]->Add(CS_S1385_ToSm[1][0]);
  CS_S1385_ToSm[0][0]->Add(CS_S1385_ToSm[2][0]);
  CS_S1385_ToSpSm[0][0]->Add(CS_S1385_ToSpSm[1][0]);
  CS_S1385_ToSpSm[0][0]->Add(CS_S1385_ToSpSm[2][0]);
  CS_S1385_ToSp[0][2]->Add(CS_S1385_ToSp[1][2]);
  CS_S1385_ToSp[0][2]->Add(CS_S1385_ToSp[2][2]);
  CS_S1385_ToSm[0][2]->Add(CS_S1385_ToSm[1][2]);
  CS_S1385_ToSm[0][2]->Add(CS_S1385_ToSm[2][2]);
  CS_S1385_ToSpSm[0][2]->Add(CS_S1385_ToSpSm[1][2]);
  CS_S1385_ToSpSm[0][2]->Add(CS_S1385_ToSpSm[2][2]);

  TH2F *CS_S1385_ToSp_coscut[3][3];//q-cut,isys
  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSp_coscut[iq][isys] = new TH2F(Form("CS_S1385_ToSp_coscut%d_%d",iq,isys),Form("CS_S1385_ToSp_coscut%d_%d",iq,isys),60,1.2,2.1,30,0,1.5);
    }
  }
  TH2F *CS_S1385_coscut_ref = new TH2F("CS_S1385_coscut_ref","CS_S1385_coscut_ref",60,1.2,2.1,30,0,1.5);
  
  TH2F *CS_S1385_ToSm_coscut[3][3];//q,isys
  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSm_coscut[iq][isys] = new TH2F(Form("CS_S1385_ToSm_coscut%d_%d",iq,isys),Form("CS_S1385_ToSm_coscut%d_%d",iq,isys),60,1.2,2.1,30,0,1.5);
    }
  }
  
  TH2F *CS_S1385_ToSpSm_coscut[3][3];//q,isys
  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){  
      CS_S1385_ToSpSm_coscut[iq][isys] = new TH2F(Form("CS_S1385_ToSpSm_coscut%d_%d",iq,isys),Form("CS_S1385_ToSpSm_coscut%d_%d",iq,isys),60,1.2,2.1,30,0,1.5);
    }
  }
  

  const int nrand = 1000000;
  for(int i=0;i<nrand;i++){
    double m,q;
    CS_lpim_qcut[0]->GetRandom2(m,q);//qall
    double cosCM = funcos(q,m);
    if(costhetacutCMlo < cosCM && cosCM <costhetacutCMhi){ 
      if(q<=0.35){
        for(int isys=0;isys<3;isys++){
          CS_S1385_ToSp_coscut[1][isys]->Fill(m,q);//qlow, isys
          CS_S1385_ToSm_coscut[1][isys]->Fill(m,q);//qlow, isys
        }
      }else if(0.35<q && q<0.65){
        for(int isys=0;isys<3;isys++){
          CS_S1385_ToSp_coscut[2][isys]->Fill(m,q);//qhi, sysdef
          CS_S1385_ToSm_coscut[2][isys]->Fill(m,q);//qhi, sysdef
        }
      }
    }
    CS_S1385_coscut_ref->Fill(m,q);
  }
  double scaleF = CS_S1385_coscut_ref->Integral();
  double scaleR = CS_lpim_qcut[0]->Integral();
  for(int isys=0;isys<3;isys++){
    CS_S1385_ToSp_coscut[1][isys]->Scale(scaleR/scaleF);
    CS_S1385_ToSp_coscut[2][isys]->Scale(scaleR/scaleF);
    CS_S1385_ToSm_coscut[1][isys]->Scale(scaleR/scaleF);
    CS_S1385_ToSm_coscut[2][isys]->Scale(scaleR/scaleF);
  }

  
  CS_S1385_ToSm_coscut[0][sysdef]->Add(CS_S1385_ToSm_coscut[1][sysdef]);
  CS_S1385_ToSm_coscut[0][sysdef]->Add(CS_S1385_ToSm_coscut[2][sysdef]);

  CS_S1385_ToSpSm_coscut[1][sysdef]->Add(CS_S1385_ToSp_coscut[1][sysdef]);
  CS_S1385_ToSpSm_coscut[1][sysdef]->Add(CS_S1385_ToSm_coscut[1][sysdef]);
  CS_S1385_ToSpSm_coscut[2][sysdef]->Add(CS_S1385_ToSp_coscut[2][sysdef]);
  CS_S1385_ToSpSm_coscut[2][sysdef]->Add(CS_S1385_ToSm_coscut[2][sysdef]);
  CS_S1385_ToSpSm_coscut[0][sysdef]->Add(CS_S1385_ToSp_coscut[1][sysdef]);
  CS_S1385_ToSpSm_coscut[0][sysdef]->Add(CS_S1385_ToSp_coscut[2][sysdef]);
  CS_S1385_ToSpSm_coscut[0][sysdef]->Add(CS_S1385_ToSm_coscut[1][sysdef]);
  CS_S1385_ToSpSm_coscut[0][sysdef]->Add(CS_S1385_ToSm_coscut[2][sysdef]);

  TCanvas *ctestcoscutSp_d = new TCanvas("ctestcoscutSp_d","ctestcoscutSp_d");
  //sys down,qall sysdown
  const int sysdown = 0;
  CS_S1385_ToSp_coscut[1][sysdown]->Scale(0.5);
  CS_S1385_ToSp_coscut[2][sysdown]->Scale(0.9);
  CS_S1385_ToSp_coscut[0][sysdown]->Add(CS_S1385_ToSp_coscut[1][sysdown]);
  CS_S1385_ToSp_coscut[0][sysdown]->Add(CS_S1385_ToSp_coscut[2][sysdown]);
  CS_S1385_ToSp_coscut[0][sysdown]->Draw("colz");

  CS_S1385_ToSm_coscut[1][sysdown]->Scale(0.5);
  CS_S1385_ToSm_coscut[2][sysdown]->Scale(0.9);
  CS_S1385_ToSm_coscut[0][sysdown]->Add(CS_S1385_ToSm_coscut[1][sysdown]);
  CS_S1385_ToSm_coscut[0][sysdown]->Add(CS_S1385_ToSm_coscut[2][sysdown]);
  
  CS_S1385_ToSpSm_coscut[1][sysdown]->Add(CS_S1385_ToSp_coscut[1][sysdown]);
  CS_S1385_ToSpSm_coscut[1][sysdown]->Add(CS_S1385_ToSm_coscut[1][sysdown]);
  CS_S1385_ToSpSm_coscut[2][sysdown]->Add(CS_S1385_ToSp_coscut[2][sysdown]);
  CS_S1385_ToSpSm_coscut[2][sysdown]->Add(CS_S1385_ToSm_coscut[2][sysdown]);
  CS_S1385_ToSpSm_coscut[0][sysdown]->Add(CS_S1385_ToSp_coscut[1][sysdown]);
  CS_S1385_ToSpSm_coscut[0][sysdown]->Add(CS_S1385_ToSp_coscut[2][sysdown]);
  CS_S1385_ToSpSm_coscut[0][sysdown]->Add(CS_S1385_ToSm_coscut[1][sysdown]);
  CS_S1385_ToSpSm_coscut[0][sysdown]->Add(CS_S1385_ToSm_coscut[2][sysdown]);


  TCanvas *ctestcoscutSp_u = new TCanvas("ctestcoscutSp_u","ctestcoscutSp_u");
  //sys up, qall sysup
  const int sysup = 2;
  CS_S1385_ToSp_coscut[1][sysup]->Scale(2);
  CS_S1385_ToSp_coscut[2][sysup]->Scale(1.1);
  CS_S1385_ToSp_coscut[0][sysup]->Add(CS_S1385_ToSp_coscut[1][sysup]);
  CS_S1385_ToSp_coscut[0][sysup]->Add(CS_S1385_ToSp_coscut[2][sysup]);
  CS_S1385_ToSp_coscut[0][sysup]->Draw("colz");
  
  CS_S1385_ToSm_coscut[1][sysup]->Scale(2);
  CS_S1385_ToSm_coscut[2][sysup]->Scale(1.1);
  CS_S1385_ToSm_coscut[0][sysup]->Add(CS_S1385_ToSm_coscut[1][sysup]);
  CS_S1385_ToSm_coscut[0][sysup]->Add(CS_S1385_ToSm_coscut[2][sysup]);

  CS_S1385_ToSpSm_coscut[1][sysup]->Add(CS_S1385_ToSp_coscut[1][sysup]);
  CS_S1385_ToSpSm_coscut[1][sysup]->Add(CS_S1385_ToSm_coscut[1][sysup]);
  CS_S1385_ToSpSm_coscut[2][sysup]->Add(CS_S1385_ToSp_coscut[2][sysup]);
  CS_S1385_ToSpSm_coscut[2][sysup]->Add(CS_S1385_ToSm_coscut[2][sysup]);
  CS_S1385_ToSpSm_coscut[0][sysup]->Add(CS_S1385_ToSp_coscut[1][sysup]);
  CS_S1385_ToSpSm_coscut[0][sysup]->Add(CS_S1385_ToSp_coscut[2][sysup]);
  CS_S1385_ToSpSm_coscut[0][sysup]->Add(CS_S1385_ToSm_coscut[1][sysup]);
  CS_S1385_ToSpSm_coscut[0][sysup]->Add(CS_S1385_ToSm_coscut[2][sysup]);


  TCanvas *ctestcoscutSp_qlo = new TCanvas("ctestcoscutSp_qlo","ctestcoscutSp_qlo");
  CS_S1385_ToSp_coscut[1][1]->Draw("colz");
  
  TCanvas *ctestcoscutSp_qhi = new TCanvas("ctestcoscutSp_qhi","ctestcoscutSp_qhi");
  CS_S1385_ToSp_coscut[2][1]->Draw("colz");

  TCanvas *ctestcoscutSp_ref = new TCanvas("ctestcoscutSp_ref","ctestcoscutSp_ref");
  CS_S1385_coscut_ref->Draw("colz");

  TCanvas *ctestcoscutSp = new TCanvas("ctestcoscutSp","ctestcoscutSp");
  CS_S1385_ToSp_coscut[0][sysdef]->Add(CS_S1385_ToSp_coscut[1][sysdef]);
  CS_S1385_ToSp_coscut[0][sysdef]->Add(CS_S1385_ToSp_coscut[2][sysdef]);
  double maxref = CS_S1385_coscut_ref->GetMaximum();
  //CS_S1385_ToSp_coscut[0][sysdef]->SetMaximum(maxref);
  CS_S1385_ToSp_coscut[0][sysdef]->Draw("colz");

  TCanvas *clpim = new TCanvas("clpim","clpim",1000,800);
  clpim->cd();
  //CS_lpim_sum->SetMaximum(maxref);
  CS_lpim_sum->Draw("colz");
  //CS_lpim_fit->Draw("boxsame");



  TFile *flpim_calc = new TFile("CSLpimFit_calc.root","RECREATE");
  flpim_calc->cd();
  CS_lpim_sum->Write();
  CS_lpim_fit->Write();
  for(int iq=0;iq<3;iq++){
    for(int isys=0;isys<3;isys++){
      CS_S1385_ToSp_coscut[iq][isys]->Write();
      CS_S1385_ToSm_coscut[iq][isys]->Write();
      CS_S1385_ToSpSm_coscut[iq][isys]->Write();
      CS_S1385_ToSp[iq][isys]->Write();
      CS_S1385_ToSm[iq][isys]->Write();
      CS_S1385_ToSpSm[iq][isys]->Write();
    }
  }



}

double funcos(double q, double M)
{
  double x[1];
  x[0]=q;
  double par[1];
  par[0]=M;
  double costhetaCM = func_C(x,par)/TMath::Sqrt(func_C(x,par)*func_C(x,par)+func_S2(x,par));
  return costhetaCM;

}

double func_EM(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = TMath::Sqrt(par[0]*par[0]+q*q);
  return f;
}

double func_cq(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = (n_mass*n_mass+pK*pK+q*q-TMath::Power(EK+d_mass-func_EM(x,par),2))/(2*q*pK);
  return f;
}

double func_C(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = (EK+d_mass)*(pK-q*func_cq(x,par))-pK*(EK+d_mass-func_EM(x,par));
  return f;
}


double func_S2(double *x, double *par)
// x      = q
// par[0] = m
{
  double q = x[0];
  double f = q*q*(1-func_cq(x,par)*func_cq(x,par))*(d_mass*d_mass+2*d_mass*EK+K_mass*K_mass);
  return f;
}
