#include "../post/weightfunc.h"

//woK0
TF1 *fweight_q_v300 = NULL;
TF1 *fweight_q_v302 = NULL;
TF1 *fweight_q_v304 = NULL;
TF1 *fweight_q_v309 = NULL;
TF1 *fweight_q_v311 = NULL;
TF1 *fweight_q_v342 = NULL;
TF1 *fweight_q_v344 = NULL;
TF1 *fweight_MMnmiss_v301 = NULL;
TF1 *fweight_MMnmiss_v310 = NULL;
TF1 *fweight_MMnmiss_v315 = NULL;
TF1 *fweight_MMnmiss_v343 = NULL;
TF1 *fweight_nmom_v303 = NULL;
TF1 *fweight_nmom_v305 = NULL;
TF1 *fweight_nmom_v317 = NULL;
TF1 *fweight_nmom_v341 = NULL;
TF1 *fweight_nmom_v345 = NULL;
TF1 *fweight_IMnpip_v307 = NULL;
TF1 *fweight_IMnpip_v314 = NULL;
TF1 *fweight_IMnpip_v327 = NULL;
TF1 *fweight_IMnpim_v308 = NULL;
TF1 *fweight_IMnpim_v313 = NULL;
TF1 *fweight_IMnpim_v328 = NULL;
TF1 *fweight_IMpippim_v306 = NULL;
TF1 *fweight_IMpippim_v312 = NULL;
TF1 *fweight_IMpippim_v329 = NULL;


//wK0
TF1 *fweight_q_wK0_v308 = NULL;
TF1 *fweight_q_wK0_v310 = NULL;
TF1 *fweight_q_wK0_v312 = NULL;
TF1 *fweight_q_wK0_v314 = NULL;
TF1 *fweight_q_wK0_v318 = NULL;
TF1 *fweight_q_wK0_v320 = NULL;
TF1 *fweight_q_wK0_v325 = NULL;
TF1 *fweight_q_wK0_v328 = NULL;
TF1 *fweight_q_wK0_v338 = NULL;
TF1 *fweight_q_wK0_v342 = NULL;
TF1 *fweight_q_wK0_v344 = NULL;
TF1 *fweight_MMnmiss_wK0_v309 = NULL;
TF1 *fweight_MMnmiss_wK0_v315 = NULL;
TF1 *fweight_MMnmiss_wK0_v319 = NULL;
TF1 *fweight_MMnmiss_wK0_v324 = NULL;
TF1 *fweight_MMnmiss_wK0_v327 = NULL;
TF1 *fweight_MMnmiss_wK0_v332 = NULL;
TF1 *fweight_MMnmiss_wK0_v334 = NULL;
TF1 *fweight_MMnmiss_wK0_v337 = NULL;
TF1 *fweight_MMnmiss_wK0_v340 = NULL;
TF1 *fweight_MMnmiss_wK0_v343 = NULL;
TF1 *fweight_nmom_wK0_v311 = NULL;
TF1 *fweight_nmom_wK0_v313 = NULL;
TF1 *fweight_nmom_wK0_v326 = NULL;
TF1 *fweight_nmom_wK0_v333 = NULL;
TF1 *fweight_nmom_wK0_v339 = NULL;
TF1 *fweight_nmom_wK0_v341 = NULL;
TF1 *fweight_nmom_wK0_v345 = NULL;
TF1 *fweight_IMnpip_wK0_v316 = NULL;
TF1 *fweight_IMnpip_wK0_v321 = NULL;
TF1 *fweight_IMnpip_wK0_v323 = NULL;
TF1 *fweight_IMnpip_wK0_v330 = NULL;
TF1 *fweight_IMnpip_wK0_v335 = NULL;
TF1 *fweight_IMnpim_wK0_v317 = NULL;
TF1 *fweight_IMnpim_wK0_v322 = NULL;
TF1 *fweight_IMnpim_wK0_v329 = NULL;
TF1 *fweight_IMnpim_wK0_v331 = NULL;
TF1 *fweight_IMnpim_wK0_v336 = NULL;



Double_t func_qmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return fweight_q_v300->Eval(xx)*
         fweight_q_v302->Eval(xx)*
         fweight_q_v304->Eval(xx)*
         fweight_q_v309->Eval(xx)*
         fweight_q_v311->Eval(xx)*
         fweight_q_v342->Eval(xx)*
         fweight_q_v344->Eval(xx);
}

Double_t func_MMmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_MMnmiss_v301->Eval(xx)*
  fweight_MMnmiss_v310->Eval(xx)*
  fweight_MMnmiss_v315->Eval(xx)*
  fweight_MMnmiss_v343->Eval(xx);
}

Double_t func_nmommul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_nmom_v303->Eval(xx)*
  fweight_nmom_v305->Eval(xx)*
  fweight_nmom_v317->Eval(xx)*
  fweight_nmom_v341->Eval(xx)*
  fweight_nmom_v345->Eval(xx);
}

Double_t func_IMnpipmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_IMnpip_v307->Eval(xx)*
  fweight_IMnpip_v314->Eval(xx)*
  fweight_IMnpip_v327->Eval(xx);
  //fweight_IMnpip_v346->Eval(xx);
}

Double_t func_IMnpipmul_s(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<1.92){
    return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))*(1./(1.0+exp((x[0]-1.25)/par[11])))+
           (par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))*(1.0-1./(1.0+exp((x[0]-1.25)/par[11])));
  }else{
    return 1.;
  }
}

Double_t param_IMnpip_s[12]={
-327.854,
748.394,
-561.008,
138.25,
-3277.16,
12041.2,
-18359.5,
14891.1,
-6782.67,
1645.87,
-166.272,
0.01
};



Double_t func_IMnpimmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_IMnpim_v308->Eval(xx)*
  fweight_IMnpim_v313->Eval(xx)*
  fweight_IMnpim_v328->Eval(xx);
}


Double_t func_IMpippimmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return
  fweight_IMpippim_v306->Eval(xx)*
  fweight_IMpippim_v312->Eval(xx)*
  fweight_IMpippim_v329->Eval(xx);
}


Double_t func_q_wK0mul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_q_wK0_v308->Eval(xx)*
  fweight_q_wK0_v310->Eval(xx)*  
  fweight_q_wK0_v312->Eval(xx)*  
  fweight_q_wK0_v314->Eval(xx)*  
  fweight_q_wK0_v318->Eval(xx)*  
  fweight_q_wK0_v320->Eval(xx)*  
  fweight_q_wK0_v325->Eval(xx)*  
  fweight_q_wK0_v328->Eval(xx)*  
  fweight_q_wK0_v338->Eval(xx)*  
  fweight_q_wK0_v342->Eval(xx)*  
  fweight_q_wK0_v344->Eval(xx);  
}

Double_t func_MMnmiss_wK0mul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_MMnmiss_wK0_v309->Eval(xx)*
  fweight_MMnmiss_wK0_v315->Eval(xx)* 
  fweight_MMnmiss_wK0_v319->Eval(xx)* 
  fweight_MMnmiss_wK0_v324->Eval(xx)* 
  fweight_MMnmiss_wK0_v327->Eval(xx)* 
  fweight_MMnmiss_wK0_v332->Eval(xx)* 
  fweight_MMnmiss_wK0_v334->Eval(xx)* 
  fweight_MMnmiss_wK0_v337->Eval(xx)* 
  fweight_MMnmiss_wK0_v340->Eval(xx)* 
  fweight_MMnmiss_wK0_v343->Eval(xx); 
}


Double_t func_nmom_wK0mul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_nmom_wK0_v311->Eval(xx)*
  fweight_nmom_wK0_v313->Eval(xx)*    
  fweight_nmom_wK0_v326->Eval(xx)*    
  fweight_nmom_wK0_v333->Eval(xx)*    
  fweight_nmom_wK0_v339->Eval(xx)*    
  fweight_nmom_wK0_v341->Eval(xx)*    
  fweight_nmom_wK0_v345->Eval(xx);
}


Double_t func_IMnpip_wK0mul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_IMnpip_wK0_v316->Eval(xx)*
  fweight_IMnpip_wK0_v321->Eval(xx)*
  fweight_IMnpip_wK0_v323->Eval(xx)*
  fweight_IMnpip_wK0_v330->Eval(xx)*
  fweight_IMnpip_wK0_v335->Eval(xx);
}



Double_t func_IMnpim_wK0mul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_IMnpim_wK0_v317->Eval(xx)*
  fweight_IMnpim_wK0_v322->Eval(xx)*
  fweight_IMnpim_wK0_v329->Eval(xx)*
  fweight_IMnpim_wK0_v331->Eval(xx)*
  fweight_IMnpim_wK0_v336->Eval(xx);
}








void FakeMCWeight()
{
  //TFile *forg = TFile::Open("comp_fakedata_out_org.root","READ");
  //TFile *f = TFile::Open("comp_fakedata_out_v346.root","READ");

  TCanvas *c_woK0_func = new TCanvas("c_woK0_func","c_woK0_func",1800,1000);
  c_woK0_func->Divide(3,2);
  
  //q
  c_woK0_func->cd(1);
  
  fweight_q_v300 = new TF1("fweight_q_v300",func_q,0,1.5,8);
  fweight_q_v300->SetParameters(param_q);

  fweight_q_v302 = new TF1("fweight_q_v302",func_q,0,1.5,8);
  fweight_q_v302->SetParameters(param_q_corr);
  
  fweight_q_v304 = new TF1("fweight_q_v304",func_q,0,1.5,8);
  fweight_q_v304->SetParameters(param_q_corr2);
  
  fweight_q_v309 = new TF1("fweight_q_v309",func_q,0,1.5,8);
  fweight_q_v309->SetParameters(param_q_corr3);
  
  fweight_q_v311 = new TF1("fweight_q_v311",func_q,0,1.5,8);
  fweight_q_v311->SetParameters(param_q_corr4);
  
  fweight_q_v316 = new TF1("fweight_q_v316",func_q,0,1.5,8);
  fweight_q_v316->SetParameters(param_q_corr5);
  
  fweight_q_v342 = new TF1("fweight_q_v342",func_q,0,1.5,8);
  fweight_q_v342->SetParameters(param_q_corr6);
  
  fweight_q_v344 = new TF1("fweight_q_v344",func_q,0,1.5,8);
  fweight_q_v344->SetParameters(param_q_corr7);
  //q
  TF1* f_qmul = new TF1("q",func_qmul,0,1.38,8*8);
  f_qmul->SetTitle("");
  f_qmul->GetXaxis()->SetTitle("q [GeV/c]");
  f_qmul->GetXaxis()->CenterTitle();
  f_qmul->Draw("c");
  
  //TCanvas *cqmul = new TCanvas("cqmu;","cqmul");
  //cqmul->cd();
  //TH1D* h_qmul = (TH1D*)f_qmul->GetHistogram();


  //MMnmiss
  c_woK0_func->cd(2);
  fweight_MMnmiss_v301 = new TF1("fweight_MMnmiss_v301",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_v301->SetParameters(param_MMnmiss);
  
  fweight_MMnmiss_v310 = new TF1("fweight_MMnmiss_v310",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_v310->SetParameters(param_MMnmiss_corr);

  fweight_MMnmiss_v315 = new TF1("fweight_MMnmiss_v315",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_v315->SetParameters(param_MMnmiss_corr2);
  
  fweight_MMnmiss_v343 = new TF1("fweight_MMnmiss_v343",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_v343->SetParameters(param_MMnmiss_corr3);
  
  TF1* f_MMnmissmul = new TF1("MissMass",func_MMmul,0,1.5,18*4);
  f_MMnmissmul->SetNpx(1000);
  f_MMnmissmul->SetTitle("");
  f_MMnmissmul->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  f_MMnmissmul->GetXaxis()->CenterTitle();
  f_MMnmissmul->Draw("c");
  
  //nmom
  c_woK0_func->cd(3);

  fweight_nmom_v303 = new TF1("fweight_nmom_v303",func_nmom,0,1.0,9);
  fweight_nmom_v303->SetParameters(param_nmom);
  fweight_nmom_v305 = new TF1("fweight_nmom_v305",func_nmom,0,1.0,9);
  fweight_nmom_v305->SetParameters(param_nmom_v305);
  fweight_nmom_v317 = new TF1("fweight_nmom_v317",func_nmom,0,1.0,9);
  fweight_nmom_v317->SetParameters(param_nmom_v317);
  fweight_nmom_v341 = new TF1("fweight_nmom_v341",func_nmom,0,1.0,9);
  fweight_nmom_v341->SetParameters(param_nmom_v341);
  fweight_nmom_v345 = new TF1("fweight_nmom_v345",func_nmom,0,1.0,9);
  fweight_nmom_v345->SetParameters(param_nmom_v345);

  TF1* f_nmommul = new TF1("f_nmommul",func_nmommul,0,1.0,9*5);
  f_nmommul->SetNpx(1000);
  f_nmommul->SetTitle("");
  f_nmommul->GetXaxis()->SetTitle("n_{CDS} mom. [GeV/c^{2}]");
  f_nmommul->GetXaxis()->CenterTitle();
  f_nmommul->Draw();
  
  //IMnpip
  c_woK0_func->cd(4);
  
  fweight_IMnpip_v307 = new TF1("fweight_IMnpip_v307",func_IMnpip,1,2.0,9);
  fweight_IMnpip_v307->SetParameters(param_IMnpip);
  fweight_IMnpip_v314 = new TF1("fweight_IMnpip_v314",func_IMnpip,1,2.0,9);
  fweight_IMnpip_v314->SetParameters(param_IMnpip_corr);
  fweight_IMnpip_v327 = new TF1("fweight_IMnpip_v327",func_IMnpip,1,2.0,9);
  fweight_IMnpip_v327->SetParameters(param_IMnpip_corr2);
  //fweight_IMnpip_v346 = new TF1("fweight_IMnpip_v346",func_IMnpip_corr,1,2.0,5);
  //fweight_IMnpip_v346->SetParameters(param_IMnpip_corr3);
  
  TF1* f_IMnpipmul = new TF1("f_IMnpipmul",func_IMnpipmul,1.08,2.0,9*3);
  f_IMnpipmul->SetTitle("");
  f_IMnpipmul->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  f_IMnpipmul->GetXaxis()->CenterTitle();
  //f_IMnpipmul->SetNpx(5000);
  //f_IMnpipmul->Draw("c");
  //f_IMnpipmul->SetLineColor(1);
  //TCanvas *c1 = new TCanvas("c1","c1");
  //TH1D* h_IMnpipmul = (TH1D*)f_IMnpipmul->GetHistogram();
  //h_IMnpipmul->Smooth();
  //h_IMnpipmul->SetLineColor(3);
  //h_IMnpipmul->Draw("H");
  
  TF1* f_IMnpipmul_s = new TF1("f_IMnpipmul_s",func_IMnpipmul_s,1.06,2.0,12);
  f_IMnpipmul_s->SetParameters(param_IMnpip_s);
  f_IMnpipmul_s->SetTitle("");
  f_IMnpipmul_s->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  f_IMnpipmul_s->GetXaxis()->CenterTitle();
  /*
  f_IMnpipmul_s->SetParameter(4,-3277.38);
  f_IMnpipmul_s->SetParameter(5,12041.4);
  f_IMnpipmul_s->SetParameter(6,-18359.4);
  f_IMnpipmul_s->SetParameter(7,14891.1);
  f_IMnpipmul_s->SetParameter(8,-6782.71);
  f_IMnpipmul_s->SetParameter(9,1645.85);
  f_IMnpipmul_s->SetParameter(10,-166.262);
  f_IMnpipmul_s->SetParameter(11,1.0);
  f_IMnpipmul_s->SetParLimits(11,0.01,1.0);
  h_IMnpipmul->Fit(f_IMnpipmul_s,"","",1.06,1.92);
  */
  f_IMnpipmul_s->SetLineColor(2);
  f_IMnpipmul_s->Draw("c");

  //fweight_IMnpip_v307s = new TF1("fweight_IMnpip_v307s",func_IMnpip_s,1,2.0,9);
  //fweight_IMnpip_v307s->SetParameters(param_IMnpip);
  //fweight_IMnpip_v314s = new TF1("fweight_IMnpip_v314s",func_IMnpip_s,1,2.0,9);
  //fweight_IMnpip_v314s->SetParameters(param_IMnpip_corr);
  //fweight_IMnpip_v327s = new TF1("fweight_IMnpip_v327s",func_IMnpip_s,1,2.0,9);
  //fweight_IMnpip_v327s->SetParameters(param_IMnpip_corr2);
  


  c_woK0_func->cd(5);
  
  fweight_IMnpim_v308 = new TF1("fweight_IMnpim_v308",func_IMnpim,1,2.0,8);
  fweight_IMnpim_v308->SetParameters(param_IMnpim);
  fweight_IMnpim_v313 = new TF1("fweight_IMnpim_v313",func_IMnpim,1,2.0,8);
  fweight_IMnpim_v313->SetParameters(param_IMnpim_corr);
  fweight_IMnpim_v328 = new TF1("fweight_IMnpim_v328",func_IMnpim,1,2.0,8);
  fweight_IMnpim_v328->SetParameters(param_IMnpim_corr2);

  TF1* f_IMnpimmul = new TF1("f_IMnpimmul",func_IMnpimmul,1.08,2.0,8*3);
  f_IMnpimmul->SetTitle("");
  f_IMnpimmul->GetXaxis()->SetTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  f_IMnpimmul->GetXaxis()->CenterTitle();
  f_IMnpimmul->Draw("c");

  //IMpippim
  c_woK0_func->cd(6);
  
  fweight_IMpippim_v306 = new TF1("fweight_IMpippim_v306",func_IMpippim,0,1.0,7);
  fweight_IMpippim_v306->SetParameters(param_IMpippim);
  fweight_IMpippim_v312 = new TF1("fweight_IMpippim_v312",func_IMpippim_corr,0,1.0,12);
  fweight_IMpippim_v312->SetParameters(param_IMpippim_corr);
  fweight_IMpippim_v329 = new TF1("fweight_IMpippim_v329",func_IMpippim_corr,0,1.0,12);
  fweight_IMpippim_v329->SetParameters(param_IMpippim_corr2);
  
  TF1 *f_IMpippimmul = new TF1("f_IMpippimmul",func_IMpippimmul,0,1.0,7+12*2);
  f_IMpippimmul->SetTitle("");
  f_IMpippimmul->GetXaxis()->SetTitle("IM(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  f_IMpippimmul->GetXaxis()->CenterTitle();
  f_IMpippimmul->Draw("c");

  
  TCanvas *c_wK0_func = new TCanvas("c_wK0_func","c_wK0_func",1800,1000);
  c_wK0_func->Divide(3,2);
  
  //q
  c_wK0_func->cd(1);
  
  fweight_q_wK0_v308 = new TF1("fweight_q_wK0_v308",func_q,0,1.5,8);
  fweight_q_wK0_v308->SetParameters(param_q_wK0);
  fweight_q_wK0_v310 = new TF1("fweight_q_wK0_v310",func_q,0,1.5,8);
  fweight_q_wK0_v310->SetParameters(param_q_wK0_corr);
  fweight_q_wK0_v312 = new TF1("fweight_q_wK0_v312",func_q,0,1.5,8);
  fweight_q_wK0_v312->SetParameters(param_q_wK0_corr2);
  fweight_q_wK0_v314 = new TF1("fweight_q_wK0_v314",func_q,0,1.5,8);
  fweight_q_wK0_v314->SetParameters(param_q_wK0_corr3);
  fweight_q_wK0_v318 = new TF1("fweight_q_wK0_v318",func_q,0,1.5,8);
  fweight_q_wK0_v318->SetParameters(param_q_wK0_corr4);
  fweight_q_wK0_v320 = new TF1("fweight_q_wK0_v320",func_q,0,1.5,8);
  fweight_q_wK0_v320->SetParameters(param_q_wK0_corr5);
  fweight_q_wK0_v325 = new TF1("fweight_q_wK0_v325",func_q,0,1.5,8);
  fweight_q_wK0_v325->SetParameters(param_q_wK0_corr6);
  fweight_q_wK0_v328 = new TF1("fweight_q_wK0_v328",func_q,0,1.5,8);
  fweight_q_wK0_v328->SetParameters(param_q_wK0_corr7);
  fweight_q_wK0_v338 = new TF1("fweight_q_wK0_v338",func_q_wK0,0,1.5,9);
  fweight_q_wK0_v338->SetParameters(param_q_wK0_corr8);
  fweight_q_wK0_v342 = new TF1("fweight_q_wK0_v342",func_q_wK0,0,1.5,9);
  fweight_q_wK0_v342->SetParameters(param_q_wK0_corr9);
  fweight_q_wK0_v344 = new TF1("fweight_q_wK0_v344",func_q_wK0,0,1.5,9);
  fweight_q_wK0_v344->SetParameters(param_q_wK0_corr10);


  //q
  TF1* f_qmul_wK0 = new TF1("q_wK0",func_q_wK0mul,0,1.36,8*8+9*3);
  f_qmul_wK0->SetTitle("");
  f_qmul_wK0->GetXaxis()->SetTitle("q [GeV/c]");
  f_qmul_wK0->GetXaxis()->CenterTitle();
  f_qmul_wK0->Draw("c");


  //MMnmiss
  c_wK0_func->cd(2);
  
  fweight_MMnmiss_wK0_v309 = new TF1("fweight_MMnmiss_wK0_v309",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v309->SetParameters(param_MMnmiss_wK0);
  fweight_MMnmiss_wK0_v315 = new TF1("fweight_MMnmiss_wK0_v315",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v315->SetParameters(param_MMnmiss_wK0_corr);
  fweight_MMnmiss_wK0_v319 = new TF1("fweight_MMnmiss_wK0_v319",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v319->SetParameters(param_MMnmiss_wK0_corr2);
  fweight_MMnmiss_wK0_v324 = new TF1("fweight_MMnmiss_wK0_v324",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v324->SetParameters(param_MMnmiss_wK0_corr3);
  fweight_MMnmiss_wK0_v327 = new TF1("fweight_MMnmiss_wK0_v327",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v327->SetParameters(param_MMnmiss_wK0_corr4);
  fweight_MMnmiss_wK0_v332 = new TF1("fweight_MMnmiss_wK0_v332",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v332->SetParameters(param_MMnmiss_wK0_corr5);
  fweight_MMnmiss_wK0_v334 = new TF1("fweight_MMnmiss_wK0_v334",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v334->SetParameters(param_MMnmiss_wK0_corr6);
  fweight_MMnmiss_wK0_v337 = new TF1("fweight_MMnmiss_wK0_v337",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v337->SetParameters(param_MMnmiss_wK0_corr7);
  fweight_MMnmiss_wK0_v340 = new TF1("fweight_MMnmiss_wK0_v340",func_MMnmiss_wK0,0,1.5,18);
  fweight_MMnmiss_wK0_v340->SetParameters(param_MMnmiss_wK0_corr8);
  fweight_MMnmiss_wK0_v343 = new TF1("fweight_MMnmiss_wK0_v343",func_MMnmiss,0,1.5,18);
  fweight_MMnmiss_wK0_v343->SetParameters(param_MMnmiss_wK0_corr9);


  TF1* f_MMnmiss_wK0mul = new TF1("MissMass_wK0",func_MMnmiss_wK0mul,0,1.3,18*10);
  f_MMnmiss_wK0mul->SetNpx(50);
  f_MMnmiss_wK0mul->SetTitle("");
  f_MMnmiss_wK0mul->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  f_MMnmiss_wK0mul->GetXaxis()->CenterTitle();
  f_MMnmiss_wK0mul->Draw("c");
  
  //nmom
  c_wK0_func->cd(3);
  
  fweight_nmom_wK0_v311 = new TF1("fweight_nmom_wK0_v311",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v311->SetParameters(param_nmom_wK0);
  fweight_nmom_wK0_v313 = new TF1("fweight_nmom_wK0_v313",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v313->SetParameters(param_nmom_wK0_v313);
  fweight_nmom_wK0_v326 = new TF1("fweight_nmom_wK0_v326",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v326->SetParameters(param_nmom_wK0_v326);
  fweight_nmom_wK0_v333 = new TF1("fweight_nmom_wK0_v333",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v333->SetParameters(param_nmom_wK0_v333);
  fweight_nmom_wK0_v339 = new TF1("fweight_nmom_wK0_v339",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v339->SetParameters(param_nmom_wK0_v339);
  fweight_nmom_wK0_v341 = new TF1("fweight_nmom_wK0_v341",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v341->SetParameters(param_nmom_wK0_v341);
  fweight_nmom_wK0_v345 = new TF1("fweight_nmom_wK0_v345",func_nmom,0,1.0,9);
  fweight_nmom_wK0_v345->SetParameters(param_nmom_wK0_v345);

  TF1* f_nmom_wK0mul = new TF1("f_nmom_wK0",func_nmom_wK0mul,0,1.0,9*7);
  f_nmom_wK0mul->SetNpx(1000);
  f_nmom_wK0mul->SetTitle("");
  f_nmom_wK0mul->GetXaxis()->SetTitle("n_{CDS} mom. [GeV/c^{2}]");
  f_nmom_wK0mul->GetXaxis()->CenterTitle();
  f_nmom_wK0mul->Draw();
  
  //IMnpip
  c_wK0_func->cd(4);
  
  fweight_IMnpip_wK0_v316 = new TF1("fweight_IMnpip_wK0_v316",func_IMnpip,1,2.0,9);
  fweight_IMnpip_wK0_v316->SetParameters(param_IMnpip_wK0);
  fweight_IMnpip_wK0_v321 = new TF1("fweight_IMnpip_wK0_v321",func_IMnpip,1,2.0,9);
  fweight_IMnpip_wK0_v321->SetParameters(param_IMnpip_wK0_corr);
  fweight_IMnpip_wK0_v323 = new TF1("fweight_IMnpip_wK0_v323",func_IMnpip_wK0_corr,1,2.0,13);
  fweight_IMnpip_wK0_v323->SetParameters(param_IMnpip_wK0_corr2);
  fweight_IMnpip_wK0_v330 = new TF1("fweight_IMnpip_wK0_v330",func_IMnpip_wK0_corr2,1,2.0,12);
  fweight_IMnpip_wK0_v330->SetParameters(param_IMnpip_wK0_corr3);
  fweight_IMnpip_wK0_v335 = new TF1("fweight_IMnpip_wK0_v335",func_IMnpip_wK0_corr2,1,2.0,12);
  fweight_IMnpip_wK0_v335->SetParameters(param_IMnpip_wK0_corr4);


  TF1* f_IMnpip_wK0mul = new TF1("f_IMnpipmul",func_IMnpip_wK0mul,1.08,2.0,9*2+13+12*2);
  f_IMnpip_wK0mul->SetTitle("");
  f_IMnpip_wK0mul->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  f_IMnpip_wK0mul->GetXaxis()->CenterTitle();
  f_IMnpip_wK0mul->Draw("c");
  

  c_wK0_func->cd(5);
  
  fweight_IMnpim_wK0_v317 = new TF1("fweight_IMnpim_wK0_v317",func_IMnpim_wK0,1,2.0,9);
  fweight_IMnpim_wK0_v317->SetParameters(param_IMnpim_wK0);
  fweight_IMnpim_wK0_v322 = new TF1("fweight_IMnpim_wK0_v322",func_IMnpim_wK0,1,2.0,9);
  fweight_IMnpim_wK0_v322->SetParameters(param_IMnpim_wK0_corr);
  fweight_IMnpim_wK0_v329 = new TF1("fweight_IMnpim_wK0_v329",func_IMnpim_wK0,1,2.0,9);
  fweight_IMnpim_wK0_v329->SetParameters(param_IMnpim_wK0_corr2);
  fweight_IMnpim_wK0_v331 = new TF1("fweight_IMnpim_wK0_v331",func_IMnpim_wK0,1,2.0,9);
  fweight_IMnpim_wK0_v331->SetParameters(param_IMnpim_wK0_corr3);
  fweight_IMnpim_wK0_v336 = new TF1("fweight_IMnpim_wK0_v336",func_IMnpim_wK0,1,2.0,9);
  fweight_IMnpim_wK0_v336->SetParameters(param_IMnpim_wK0_corr4);
  
  
  TF1* f_IMnpim_wK0mul = new TF1("f_IMnpim_wK0mul",func_IMnpim_wK0mul,1.08,2.0,9*5);
  f_IMnpim_wK0mul->SetTitle("");
  f_IMnpim_wK0mul->GetXaxis()->SetTitle("IM(n#pi^{-}) [GeV/c^{2}]");
  f_IMnpim_wK0mul->GetXaxis()->CenterTitle();
  f_IMnpim_wK0mul->Draw("c");

  //IMpippim (N/A)
  c_wK0_func->cd(6);
  
  

  std::ofstream os;
  os.open("param_corr.txt");
  os << "IMnpip" << endl;
  for(int i=0;i<f_IMnpipmul_s->GetNpar();i++){
    os << f_IMnpipmul_s->GetParameter(i) << ",";
    os << endl;
  }

};
