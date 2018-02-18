#define RESET 0
TString hname="EdT";
TString psname=Form("blslew.pdf");
TString name[6]={"BHD","T0","CVC","NC","PC","BPD"};
int cid[6]={2,4,14,32,35,41};
int nSeg[6]={20,5,34,112,27,48};
double edep[6]={1.,1.957,5.87,10,5.87,1.};
//double tof[5]={25.93,53.88,50.51,53.88,50.51};
//double tof[5]={25.93,25.93.,57.,53.88,50.51};
//double tof[6]={25.93,57.,57.,53.88,50.51,5.};
//  double tof[5]={25.93,25.93,47.15,50.5,47.15};
double tof[6]={25.93,50.5,47.15,50.5,47.15,5};
TString ud[2]={"U","D"};
TH1F *h1;
TGraphErrors *gra;
TH2F* h2;
TCanvas *c1;
TFile *f1;
ConfMan* conf;

TFile *f[100];

char* aaa="49";
int run=1;
int run2=1;
char* prefix="hodo";

void Slewing()
{
  /*** load library ***/
  gSystem->Load("~/lib/libAll.so");

  /*** assign input file & call tree ***/
  //  f1 = new TFile("~/data/k18ana/run47/hodo57.root");
  f1 = new TFile("~/data/run49/hodo_20130910_1_1.root");
  
  /*** conf file for new parameters ***/
  conf = new ConfMan("conf/Run49/analyzer.conf");
  //  conf = new ConfMan("conf/Run40/analyzer.conf");
  //conf = new ConfMan("conf/Run35/analyzer.conf");
  conf->Initialize();

  for(int irun=run;irun<=run2;irun++)
    f[irun-run]=new TFile(Form("~/data/run%s/%s%d.root",aaa,prefix,irun));

#if RESET
  int ith=0, type=1, npar=2;
  double par[2]={0.,0.};
  for(int ic=0;ic<5;ic++)
    for(int iseg=1;iseg<=nSeg[ic];iseg++)
      for(int iud=0;iud<2;iud++)
	conf->GetSlewingMapManager()->SetParam( cid[ic], iseg, iud, ith, type, npar, par );	
#else
  std::cout<<"start!!!"<<std::endl;
  c1=new TCanvas();
  c1->Print(psname+"[");
  //  SlewBPD2(2);
  SlewBHD(2);
  //  SlewT0_BHD(3);
  //  SlewT0_BPD2(2);
  SlewT0_CVC(2);
  SlewCVC(2);
  //  SlewNC(1.);
  //  SlewT0_NC(3);
  //  SlewNC2(1.5);
  SlewPC(2);
  c1->Print(psname+"]");

#endif
  //  ofstream ofs(Form("blslew%d.dat",itai));
  ofstream ofs(Form("blslew.dat"));
  conf->GetSlewingMapManager()->PrintMapBL(ofs);
  ofs.close();
  ofs.open(Form("blslewoffs.param"));
  conf->GetGainMapManager()->PrintMapBL(ofs);
  ofs.close();

  //gFile->Write();
  //gFile->Close();

}
SlewBHD(double factor){
  int ic=0;
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s Slewing Correction",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s Slewing Correction",name[ic].Data(),ud[iud].Data()),
		   800,800);
    c1->Divide(4,5);
    for(int seg=1;seg<=nSeg[ic];seg++){           
      c1->cd(seg);
      h2=(TH2F*)f[0]->Get(Form("%sBHD%s%dT0",
			    hname.Data(),ud[iud].Data(),seg));
      //      h2->Draw();
      if(h2)
	h2->Rebin2D(2,2);
      calc(h2,conf,ic,iud,seg,factor);
    }
      c1->Print(psname);
  }
}


SlewT0_CVC(double factor){
  int ic=1;
  c1=new TCanvas(Form("%s Slewing Correction",name[ic].Data()),
		 Form("%s Slewing Correction",name[ic].Data()),
		 600,600);
  c1->Divide(3,4);
  for(int iud=0;iud<2;iud++){
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(seg+iud*6);
      h2=(TH2F*)f1->Get(Form("%sT0%s%dCVC",
			     hname.Data(),ud[iud].Data(),seg));
      h2->Draw();
      h2->Rebin2D(2,2);
      calc(h2,conf,ic,iud,seg,factor);
    }
  }
  c1->Print(psname);
}
SlewCVC(double factor){    
  int ic=2;
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s Slewing Correction",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s Slewing Correction",name[ic].Data(),ud[iud].Data()),
		   1500,750);
    c1->Divide(9,4);
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(seg);
      h2=(TH2F*)f1->Get(Form("%sCVC%s%dT0",
			    hname.Data(),ud[iud].Data(),seg));
      h2->Draw();
      h2->Rebin2D(4,2);
      calc(h2,conf,ic,iud,seg,factor,-1);
    }
    c1->Print(psname);
  }
}
SlewPC(double factor){  
  int ic=4;
  for(int iud=0;iud<2;iud++){
    c1=new TCanvas(Form("%s %s Slewing Correction",name[ic].Data(),ud[iud].Data()),
		   Form("%s %s Slewing Correction",name[ic].Data(),ud[iud].Data()),
		   1500,750);
    c1->Divide(7,4);
    for(int seg=1;seg<=nSeg[ic];seg++){      
      c1->cd(seg);
      h2=(TH2F*)f1->Get(Form("%sPC%s%dT0",
			    hname.Data(),ud[iud].Data(),seg));
      h2->Rebin2D(4,2);
      calc(h2,conf,ic,iud,seg,factor,-1);
    }
    c1->Print(psname);
  }
}
  


void calc(TH2F *h1,ConfMan* conf, int ic, int iud,int seg,double factor,int sign=1,double tofmean=-9999,double xmin=-9999,double xmax=-9999){
  if(!h1) return;
  TString ud[2]={"U","D"};
  double tempxmax[6]={1.5,5,15,30,15,3.};
  if(xmin<0) xmin=tempxmax[ic]*0.1;
  if(xmax<0) xmax=tempxmax[ic]*0.9;
  //  std::cout<<xmin<<"\t"<<xmax<<std::endl;
  double par[3]={0.,0.,0.};
  TH1D *h1_1;
  int ith=0, type=1, npar=2;
  if(tofmean<-999) tofmean=tof[ic];
  //  std::cout<<tofmean<<std::endl;
  //  h1->GetYaxis()->SetRangeUser(tofmean-3,tofmean+3);
  h1->GetXaxis()->SetRangeUser(0,tempxmax[ic]);
  TF1 *slewfunc = new TF1( "slewfunc", "[0]/sqrt(x)+[1]+[2]*x", 0.,tempxmax[ic]);
  //  TF1 *slewfunc = new TF1( "slewfunc", "[0]/sqrt(x)+[1]", 0, tempxmax[ic]);
  slewfunc->SetNpx(999);
  slewfunc->SetParameters(0.,0,0);
  //  slewfunc->FixParameter(0,-1);
  if(h1->GetEntries()<100){
    conf->GetSlewingMapManager()
      ->GetParam( cid[ic],nSeg[ic]/2 , iud, ith, type, npar, par );	
    conf->GetSlewingMapManager()
      ->SetParam( cid[ic], seg, iud, ith, type, npar, par );	
    slewfunc->SetParameters(par);
    h1->Draw("col");
  }else{
    TF1* gaus=new TF1("gaus","gaus(0)",-3,3);
    gaus->SetParameters(10,0,0.2);
    gaus->SetRange(-3,3);
    if(h1->GetEntries()>2000){
      h1->FitSlicesY(gaus,0,-1,3,"QNRG2S");
      //      h1->FitSlicesY(0,0,-1,3,"QNRG2S");
    }else {
      h1->FitSlicesY(gaus,0,-1,3,"QNRG3S");
      //      h1->FitSlicesY(0,0,-1,3,"QNRG3S");
    }
    TH1F* h3=(TH1F*)gFile->Get(Form("%s_1",(h1->GetName())));
    if(h1->GetEntries()>100){
      //      h3->Fit("slewfunc","0q","",xmin,xmax);
      h1->Fit("slewfunc","0q","",xmin,xmax);
    }else {
      slewfunc->FixParameter(2,0);
      h1->Fit("slewfunc","0q","",xmin-1.,xmax*0.8);
    }

    //       h1->ProfileX()->Fit("slewfunc","0q");    


    h1->Draw("col");
    h1->GetYaxis()->SetRangeUser(tofmean-2,tofmean+2);
    //    h1->ProfileX()->Draw("same");
    h3->Draw("same");
    double p0 = slewfunc->GetParameter(0);
    double p1 = slewfunc->GetParameter(1);
    double p2 = slewfunc->GetParameter(2);
    conf->GetSlewingMapManager()
      ->GetParam( cid[ic],seg , iud, ith, type, npar, par );	
    //    if(ic==0)   par[0]+=p0/factor;
    //    if(ic==1)   par[0]+=p0/factor;
    //    else        par[0]-=p0/factor;
    par[0]=par[0]+sign*p0/factor;
    par[1]=par[1]+sign*p1/factor;
    par[2]=par[2]+sign*p2/factor;
    type=2;
    npar=3;
    std::cout<<name[ic]<<seg<<ud[iud]<<" "<<par[0]<<" "<<par[1]<< " "<<par[2]<<std::endl;
    conf->GetSlewingMapManager()
      ->SetParam( cid[ic], seg, iud, ith, type, npar, par );

    int cr,sl,ch;
    double gain,tmpped;
    conf->GetCounterMapManager()->GetCNA(cid[ic],seg,1,iud,cr,sl,ch);
    conf->GetGainMapManager()->GetParam( cr,sl,ch, gain, tmpped );
    double ped = tmpped - sign*p1;
    conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, ped );
  }
  // 
  //  h1->ProfileX()->Draw("same");
  //  h1_1->Draw();
  slewfunc->SetLineColor(2);    
  slewfunc->Draw("same");
}

