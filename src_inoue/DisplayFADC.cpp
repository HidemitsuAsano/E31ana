#include <cmath>

#include "DisplayFADC.h"

ClassImp(DisplayFADC);

DisplayFADC::DisplayFADC() : TObject()
{
  frameFADC = 0;
  XMIN=XMAX=0;
}
bool DisplayFADC::Wait()
{
  char dispflag;
  char dispin[100]="";

  std::cout << " Input any word or return" << std::endl;
  std::cout << " (q:quit)" << std::endl;
  fgets(dispin,100,stdin);
  if(sscanf(dispin,"%c",&dispflag)==1){
    if( dispflag=='q' ){
      return false;
      //gSystem->Exit(1);
    }
  }
  return true;
}

void DisplayFADC::SetFADCFrame( double xmin, double xmax, double ymin, double ymax )
{
  XMIN=xmin;
  XMAX=xmax;
  frameFADC = new TH2F( "frameFADC", "", 100, xmin, xmax, 100, ymin, ymax );
  frameFADC->SetStats(0);
}

bool DisplayFADC::DrawFADCFrame( TVirtualPad *pad ){
  if( frameFADC!=0 ){
    pad->cd();
    frameFADC->SetTitle(Form("FADC wave forms of SDDs ;time (#mus);"));
    frameFADC->DrawCopy();
    return true;
  }
  return false;
}
bool DisplayFADC::DrawFADCFrameSep( TVirtualPad *pad ){
  pad->Divide(4,2);
  if( frameFADC!=0 ){
    for(int i=0;i<8;i++){
      pad->cd(i+1);
      frameFADC->SetTitle(Form("SDD%d;time (#mus);",i+1));
      frameFADC->DrawCopy();
      //      std::cout<<i<<std::endl;
    }
    return true;
  }
  return false;
}

bool DisplayFADC::DrawFADCHit( TVirtualPad *pad, ConfMan *conf, SDDHitMan *sdd)
{
  double val[64];
  double time[64];
  int nsamp,ph;
  int interval;
  TGraph* grgr=new TGraph();
  TLine* line=new TLine();
  grgr->SetMarkerStyle(20);
  line->SetLineColor(2);

  std::cout<<"start"<<std::endl;
  for(int i=0;i<sdd->nSDD();i++){
    grgr->SetMarkerColor(1);
    SDDHit* hit=sdd->SDD(i);
    nsamp=hit->nsample();
    interval=1;
    //    interval=hit->interval();
    for( int j=0; j<nsamp; j++ ){
      time[j]=interval*j;
      val[j]=hit->fadc(j);
    } 
    if(hit->tdcsdd()>700&&hit->tdcsdd()<800) grgr->SetMarkerColor(2);
    grgr->DrawGraph(nsamp,time,val,"p");
    TF1 *func=new TF1("func",Form("pol%d",hit->fpeaknpar()-1),0,nsamp);
    func->SetParameters(hit->fpeakpar());    
    func->DrawCopy("same");
    delete func;
    ph=hit->out();
    line->DrawLine(XMIN,ph,XMAX,ph);
  }
  delete line;
  delete grgr;
  return true;
}

bool DisplayFADC::DrawFADCHitSep( TVirtualPad *pad, ConfMan *conf, SDDHitMan *sdd)
{
  double val[64];
  double time[64];
  int nsamp,ph;
  TGraph* grgr=new TGraph();
  TLine* line=new TLine();

  grgr->SetMarkerStyle(20);
  line->SetLineColor(2);

  std::cout<<"start"<<std::endl;
  for(int i=0;i<sdd->nSDD();i++){
    SDDHit* hit=sdd->SDD(i);
    pad->cd(hit->seg());
    nsamp=hit->nsample();
    //    int interval=hit->interval();
    int interval=1;
    for( int j=0; j<nsamp; j++ ){
      time[j]=interval*j;
      val[j]=hit->fadc(j);
    } 
    grgr->SetMarkerColor(1);
    if(hit->tdcsdd()>700&&hit->tdcsdd()<800) grgr->SetMarkerColor(2);
    grgr->DrawGraph(nsamp,time,val,"p");
    TF1 *func=new TF1("func",
		      conf->GetFADCParamManager()->GetPeakFunc().c_str(),
		      conf->GetFADCParamManager()->GetLowerLimitPeak(),
		      conf->GetFADCParamManager()->GetUpperLimitPeak());
    func->SetParameters(hit->fpeakpar());    
    func->DrawCopy("same");
    delete func;
    func=new TF1("func","pol1",
		 conf->GetFADCParamManager()->GetLowerLimitBase(),
		 conf->GetFADCParamManager()->GetUpperLimitBase());
    func->SetParameters(hit->fbasepar());    
    func->DrawCopy("same");
    delete func;
    func=new TF1("func","expo",
		 conf->GetFADCParamManager()->GetLowerLimitPost(),
		 conf->GetFADCParamManager()->GetUpperLimitPost());
    func->SetParameters(hit->fpostpar());    
    func->DrawCopy("same");
    delete func;
    ph=hit->out();
    line->DrawLine(XMIN,ph,XMAX,ph);
  }
  delete grgr;
  delete line;
  return true;
}

