#include "../post/weightfunc.h"
#include "../post/anacuts.h"

//woK0
TF1 *fweight_q_v348 = NULL;
TF1 *fweight_MMnmiss_v352 = NULL;
TF1 *fweight_nmom_v303 = NULL;
TF1 *fweight_nmom_v305 = NULL;
TF1 *fweight_nmom_v317 = NULL;
TF1 *fweight_nmom_v341 = NULL;
TF1 *fweight_nmom_v345 = NULL;
TF1 *fweight_IMnpip_v346 = NULL;
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
  return fweight_q_v348->Eval(xx);
}

Double_t func_MMmul(Double_t *x,Double_t *par)
{
  const Double_t xx=x[0];
  return 
  fweight_MMnmiss_v352->Eval(xx);
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
  fweight_IMnpip_v346->Eval(xx);
}


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



Double_t func_MMnmiss_mod(Double_t *x,Double_t *par)
{
   //connection using woods-saxon
   if(0.0<x[0] && x[0]<1.5){
   //  return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
   //    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
   //    +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0))*(1./(1.0+exp((x[0]-1.08)/par[10])))
   //    +(1.0-1./(1.0+exp((x[0]-1.08)/par[10])))*(par[11]*exp(-0.5*pow(((x[0]-par[12])/par[13]),2.0)))*(1./(1.0+exp((x[0]-1.18)/par[14])))
   //    +(1.0-1./(1.0+exp((x[0]-1.18)/par[14])))*
   //     (par[15]+par[16]*x[0]+par[17]*pow(x[0],2.0)+par[18]*pow(x[0],3.0)+par[19]*pow(x[0],4.0));
   //  return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
   //    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
   //    +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0))*(1./(1.0+exp((x[0]-1.08)/par[9])))
   //    +(1.0-1./(1.0+exp((x[0]-1.08)/par[9])))*(par[10]*exp(-0.5*pow(((x[0]-par[11])/par[12]),2.0)))*(1./(1.0+exp((x[0]-1.18)/par[13])))
   //    +(1.0-1./(1.0+exp((x[0]-1.18)/par[13])))*
   //     (par[14]+par[15]*x[0]+par[16]*pow(x[0],2.0)+par[17]*pow(x[0],3.0)+par[18]*pow(x[0],4.0));
     return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0))*(1./(1.0+exp((x[0]-1.08)/par[8])))
       +(1.0-1./(1.0+exp((x[0]-1.08)/par[8])))*(par[9]*exp(-0.5*pow(((x[0]-par[10])/par[11]),2.0)))*(1./(1.0+exp((x[0]-1.17)/par[12])))
       +(1.0-1./(1.0+exp((x[0]-1.17)/par[12])))*
        (par[13]+par[14]*x[0]+par[15]*pow(x[0],2.0)+par[16]*pow(x[0],3.0)+par[17]*pow(x[0],4.0));
   }else{
     return 1.0;
   }
}





void FakeMCWeight()
{
  //TFile *forg = TFile::Open("comp_fakedata_out_org.root","READ");
  //TFile *f = TFile::Open("comp_fakedata_out_v346.root","READ");

  TCanvas *c_woK0_func = new TCanvas("c_woK0_func","c_woK0_func",1800,1000);
  c_woK0_func->Divide(3,2);
  
  //q
  c_woK0_func->cd(1);
  fweight_q_v348 = new TF1("fweight_q_v348",func_q,0,1.5,8);
  fweight_q_v348->SetParameters(param_q_mul);
  //q
  TF1* f_qmul = new TF1("f_qmul",func_qmul,0,1.5,8);
  
  f_qmul->SetTitle("");
  f_qmul->GetXaxis()->SetTitle("q [GeV/c]");
  f_qmul->GetXaxis()->CenterTitle();
  f_qmul->Draw("c");
  std::cout << __LINE__ << std::endl;
  //TCanvas *cqmul = new TCanvas("cqmul","cqmul");
  //cqmul->cd();
  //std::cout << __LINE__ << std::endl;
  //TH1D* h_qmul = (TH1D*)f_qmul->GetHistogram();
  //h_qmul->Draw();
  //std::cout << __LINE__ << std::endl;
  //h_qmul->Fit(f_qmul,"","",0,1.5);
 // std::cout << __LINE__ << std::endl;

  //MMnmiss
  c_woK0_func->cd(2);
  fweight_MMnmiss_v352 = new TF1("fweight_MMnmiss_v352",func_MMnmiss_mod,0,1.5,18);
  fweight_MMnmiss_v352->SetParameters(param_MMnmiss_mod);
  
  TF1* f_MMnmissmul = new TF1("MissMass",func_MMmul,0,1.5,18);
  //f_MMnmissmul->SetNpx(1000);
  f_MMnmissmul->SetTitle("");
  f_MMnmissmul->GetXaxis()->SetTitle("Miss. Mass [GeV/c^{2}]");
  f_MMnmissmul->GetXaxis()->CenterTitle();
  f_MMnmissmul->SetMinimum(0);
  f_MMnmissmul->Draw("c");
  
  TBox *box_neutron = new TBox(anacuts::neutron_MIN,0,anacuts::neutron_MAX,3);
  box_neutron->SetFillColor(4);
  box_neutron->SetFillStyle(3002);
  box_neutron->Draw();
  /*
  TCanvas *cMMnmiss_mod = new TCanvas("cMMnmiss_mod","cMMnmiss_mod");
  cMMnmiss_mod->cd();
  f_MMnmissmul->SetLineColor(3);
  TH1D* h_MMnmiss = (TH1D*)f_MMnmissmul->GetHistogram();
  h_MMnmiss->SetLineColor(3);
  h_MMnmiss->SetMarkerColor(3);
  h_MMnmiss->Draw("H");
  box_neutron->Draw();
  
  TF1 *f_MMnmiss_mod = new TF1("f_MMnmiss_mod",func_MMnmiss_mod,0.0,1.5,18);
  double param_MMnmiss_mod[18];
  param_MMnmiss_mod[0]=-3.25369;                 
  param_MMnmiss_mod[1]=56.1105;                  
  param_MMnmiss_mod[2]=-380.608;                    
  param_MMnmiss_mod[3]=1358.72;                   
  param_MMnmiss_mod[4]=-2715.42;                    
  param_MMnmiss_mod[5]=3050.94;                   
  param_MMnmiss_mod[6]=-1788.38;                    
  param_MMnmiss_mod[7]=423.797;                   
  param_MMnmiss_mod[8]=0.010;//woods saxon         
  param_MMnmiss_mod[9]= 2.558;//gaus const        
  param_MMnmiss_mod[10]= 1.116;//gaus mean         
  param_MMnmiss_mod[11]= 8.45178e-02;//gaus sigma  
  param_MMnmiss_mod[12]=0.01;//woods saxon         
  param_MMnmiss_mod[13]=720.767;                   
  param_MMnmiss_mod[14]=-2156.71;                  
  param_MMnmiss_mod[15]=2433.82;                   
  param_MMnmiss_mod[16]=-1222.3;                   
  param_MMnmiss_mod[17]=229.81;                    
  //param_MMnmiss_mod[0]=0.00906858;
  //param_MMnmiss_mod[1]=-0.485126;
  //param_MMnmiss_mod[2]=28.4734;
  //param_MMnmiss_mod[3]=-257.647;
  //param_MMnmiss_mod[4]=1120.55;
  //param_MMnmiss_mod[5]=-2570.83;
  //param_MMnmiss_mod[6]=3196.52;
  //param_MMnmiss_mod[7]=-2027.73;
  //param_MMnmiss_mod[8]=513.021;
  //param_MMnmiss_mod[9]=0.010;//woods saxon
  //param_MMnmiss_mod[10]= 2.558;//gaus const      
  //param_MMnmiss_mod[11]= 1.116;//gaus mean       
  //param_MMnmiss_mod[12]= 8.45178e-02;//gaus sigma
  //param_MMnmiss_mod[13]=0.01;//woods saxon
  //param_MMnmiss_mod[14]=720.767; 
  //param_MMnmiss_mod[15]=-2156.71;
  //param_MMnmiss_mod[16]=2433.82; 
  //param_MMnmiss_mod[17]=-1222.3; 
  //param_MMnmiss_mod[18]=229.81;  


  //param_MMnmiss_mod[0]=0.00708634;
  //param_MMnmiss_mod[1]=-0.217851;
  //param_MMnmiss_mod[2]=20.9625;
  //param_MMnmiss_mod[3]=-178.022;
  //param_MMnmiss_mod[4]=708.546;
  //param_MMnmiss_mod[5]=-1399.05;
  //param_MMnmiss_mod[6]=1273.94;
  //param_MMnmiss_mod[7]=-217.497;
  //param_MMnmiss_mod[8]=-393.59;
  //param_MMnmiss_mod[9]=186.815;
  //param_MMnmiss_mod[10]=0.015;//woods saxon
  //param_MMnmiss_mod[11]= 2.558;//gaus const
  //param_MMnmiss_mod[12]= 1.116;//gaus mean
  //param_MMnmiss_mod[13]= 8.45178e-02;//gaus sigma 
  //param_MMnmiss_mod[14]=0.01;//woods saxon
  //param_MMnmiss_mod[15]=720.767;
  //param_MMnmiss_mod[16]=-2156.71;
  //param_MMnmiss_mod[17]=2433.82;
  //param_MMnmiss_mod[18]=-1222.3;
  //param_MMnmiss_mod[19]=229.81;
  f_MMnmiss_mod->SetParameters(param_MMnmiss_mod);
  f_MMnmiss_mod->FixParameter(8,0.015);
  f_MMnmiss_mod->FixParameter(12,0.010);
  f_MMnmiss_mod->Print("v");
  h_MMnmiss->Fit("f_MMnmiss_mod","","",0.2,1.5);
  f_MMnmiss_mod->Draw("same");
  
  */
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
   
  fweight_IMnpip_v346 = new TF1("fweight_IMnpip_v346",func_IMnpipmul_s,0,2.0,12);
  fweight_IMnpip_v346->SetParameters(param_IMnpip_s);
  TF1* f_IMnpipmul = new TF1("f_IMnpipmul",func_IMnpipmul,1.06,2.0,12);
  f_IMnpipmul->SetTitle("");
  f_IMnpipmul->GetXaxis()->SetTitle("IM(n#pi^{+}) [GeV/c^{2}]");
  f_IMnpipmul->GetXaxis()->CenterTitle();
  f_IMnpipmul->SetLineColor(2);
  f_IMnpipmul->Draw("c");
  TBox *box_sigmap = new TBox(anacuts::Sigmap_MIN,0,anacuts::Sigmap_MAX,3);
  box_sigmap->SetFillColor(4);
  box_sigmap->SetFillStyle(3002);
  box_sigmap->Draw();

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
  TBox *box_sigmam = new TBox(anacuts::Sigmam_MIN,0,anacuts::Sigmam_MAX,3);
  box_sigmam->SetFillColor(4);
  box_sigmam->SetFillStyle(3002);
  box_sigmam->Draw();

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
  TBox *box_K0 = new TBox(anacuts::pipi_MIN,0,anacuts::pipi_MAX,3);
  box_K0->SetFillColor(4);
  box_K0->SetFillStyle(3002);
  box_K0->Draw();

  
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
  box_neutron->Draw();
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
  box_sigmap->Draw();

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
  box_sigmam->Draw();
  //IMpippim (N/A)
  //c_wK0_func->cd(6);
  
  

   std::ofstream os;
   os.open("param_corr.txt");
   //os << "IMnpip" << endl;
  //for(int i=0;i<f_IMnpipmul_s->GetNpar();i++){
  //  os << f_IMnpipmul_s->GetParameter(i) << ",";
  //os << "MMnmiss " << endl;
  //for(int i=0;i<f_MMnmiss_mod->GetNpar();i++){
  //  os << std::setprecision(6);
  //  os << f_MMnmiss_mod->GetParameter(i) << ",";
  //  os << endl;
  //}
  os << "q " << endl;
  for(int i=0;i<f_qmul->GetNpar();i++){
  //  os << std::setprecision(6);
    os << f_qmul->GetParameter(i) << ",";
    os << endl;
  }

};
