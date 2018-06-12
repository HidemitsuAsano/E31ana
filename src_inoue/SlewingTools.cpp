#include "SlewingTools.h"

using namespace std;

double SlewingTools::TimeReso()
{
  TH2F *h2 = (TH2F*)gFile->Get("eu_ctm");
  TH1F *h1 = projectionY(h2);
  //  CalibTools::FitGaus_2sigma(h1);
  h1-> Fit("gaus", "q", "q", h1->GetBinCenter(h1->GetMaximumBin())-2.*h1->GetRMS(), h1->GetBinCenter(h1->GetMaximumBin())+2.*h1->GetRMS());
  return h1->GetFunction("gaus")->GetParameter(2);
}

vector<double> SlewingTools::SetParam(int cid, int seg, int ud, int ith, int type, TH2F *h2, ConfMan *conf, bool add_old)
{
  vector<double> par;
  TF1 *func=Fit(h2, type);
  for( int i=0; i<func->GetNpar(); i++ ){
    par.push_back(-func->GetParameter(i));
  }

  if( add_old ){
    int npar;
    vector<double> old_par;
    if( conf-> GetSlewingMapManager()->GetParam(cid, seg, ud, ith, type, npar, old_par) ){
      if( npar!=(int)old_par.size() ){ cout<<"  !!! npar!=old_par.size() !!!"<<endl; exit(0); }
      if( npar!=(int)par.size() ){ cout<<"  !!! npar!=par.size() !!!"<<endl; exit(0); }

      for( int i=0; i<(int)par.size(); i++ ) par[i]+=old_par[i];
    }
  }
  conf-> GetSlewingMapManager()->SetParam(cid, seg, ud, ith, type, par.size(), par);
  return par;
}

vector<double> SlewingTools::SetParam(int cid, int seg, int ud, int ith, int type, TH2F *h2, ConfMan *conf, double weight, bool add_old)
{
  vector<double> par;
  TF1 *func=Fit(h2, type);
  for( int i=0; i<func->GetNpar(); i++ ){
    if( i==1 ){
      if( type!=5 ) par.push_back(-func->GetParameter(i));
      else par.push_back(-weight*func->GetParameter(i));
    }
    else par.push_back(-weight*func->GetParameter(i));
  }

  if( add_old ){
    int npar;
    vector<double> old_par;
    if( conf-> GetSlewingMapManager()->GetParam(cid, seg, ud, ith, type, npar, old_par) ){
      if( npar!=(int)old_par.size() ){ cout<<"  !!! npar!=old_par.size() !!!"<<endl; exit(0); }
      if( npar!=(int)par.size() ){ cout<<"  !!! npar!=par.size() !!!"<<endl; exit(0); }

      for( int i=0; i<(int)par.size(); i++ ) par[i]+=old_par[i];
    }
  }
  conf-> GetSlewingMapManager()->SetParam(cid, seg, ud, ith, type, par.size(), par);
  return par;
}

TF1* SlewingTools::Fit(TH2F *h2, int type)
{
  SetRange(h2);
  TProfile *prof=h2->ProfileX();
  TF1 *func=GetTF1(type);
  prof-> Fit(func, "q", "q");

  if( gPad ){
    h2-> Draw("colz");
    prof-> SetLineColor(kBlack);
    prof-> Draw("same");
    prof-> GetFunction(Form("slewing_%d", type))-> Draw("same");
  }
  return func;
}

TF1* SlewingTools::GetTF1(int type)
{
  if( gROOT->GetFunction(Form("slewing_%d", type)) ) return (TF1*)gROOT->GetFunction(Form("slewing_%d", type));
  return MakeTF1(type);
}

TF1 *SlewingTools::MakeTF1(int type)
{
  if( !gROOT->GetFunction(Form("slewing_%d", type)) ){
    if( type==1 ) return new TF1("slewing_1", "[0]/sqrt(x)+[1]");
    else if( type==2 ) return new TF1("slewing_2", "[0]/sqrt(x)+[1]+[2]*x");
    else if( type==3 ) return new TF1("slewing_3", "[0]/sqrt(x-[1])");
    else if( type==4 ) return new TF1("slewing_4", "[0]/sqrt(x)+[1]+[2]/x");
    else if( type==5 ) return new TF1("slewing_5", "[0]+[1]*x+[2]*x*x+[3]*exp([4]*x)");
    else{ cout<<"  !!! SlewingTools::MakeTF1 unknown type:"<<type<<" !!!"<<endl; return 0; }
  }
  else return (TF1*)gROOT->GetFunction(Form("slewing_%d", type));
}


void SlewingTools::MakeHist()
{
  cout<<"===== SlewingTools::MakeFile ====="<<endl;
  new TH2F("eu_tu", "eu_tu", 250, 0.0, 100, 10000, -100, 100);
  new TH2F("ed_td", "ed_td", 250, 0.0, 100, 10000, -100, 100);

  new TH2F("eu_ctu", "eu_ctu", 250, 0.0, 100, 10000, -50, 50);
  new TH2F("ed_ctd", "ed_ctd", 250, 0.0, 100, 10000, -50, 50);

  new TH2F("eu_ctm", "eu_ctm", 250, 0.0, 100, 10000, -50, 50);
  new TH2F("ed_ctm", "ed_ctm", 250, 0.0, 100, 10000, -50, 50);
}

void SlewingTools::Draw()
{
  TH2F *h2;
  TCanvas *c1 = new TCanvas("c1", "c1", 450, 600);
  c1-> Divide(2, 3);
  c1-> cd(1);

  h2 = (TH2F*)gFile-> Get("eu_tu");
  SetRange(h2);
  h2-> Draw("colz");

  c1-> cd(3);
  h2 = (TH2F*)gFile-> Get("eu_ctu");
  SetRange(h2);
  h2-> Draw("colz");

  c1-> cd(5);
  h2 = (TH2F*)gFile-> Get("eu_ctm");
  SetRange(h2);
  h2-> Draw("colz");

  c1-> cd(2);
  h2 = (TH2F*)gFile-> Get("ed_td");
  SetRange(h2);
  h2-> Draw("colz");

  c1-> cd(4);
  h2 = (TH2F*)gFile-> Get("ed_ctd");
  SetRange(h2);
  h2-> Draw("colz");

  c1-> cd(6);
  h2 = (TH2F*)gFile-> Get("ed_ctm");
  SetRange(h2);
  h2-> Draw("colz");

  c1-> Update();
  string str;
  cin>>str;
  delete c1;
}

void SlewingTools::Fill(int cid, int seg, ConfMan *conf, TNtuple *tup)
{
  TH2F *h2;
  float f_eu, f_ed, f_tu, f_td;
  tup-> SetBranchAddress("eu", &f_eu);
  tup-> SetBranchAddress("ed", &f_ed);
  tup-> SetBranchAddress("tu", &f_tu);
  tup-> SetBranchAddress("td", &f_td);

  //  cout<<"===== SlewingTools::Fill Entries : "<<tup->GetEntries()<<" ====="<<endl;
  for( int ev=0; ev<tup->GetEntries(); ev++ ){
    tup-> GetEntry(ev);
    double eu=f_eu, ed=f_ed, tu=f_tu, td=f_td;
    double ctu=conf-> GetSlewingMapManager()-> CalcCValue(cid, seg, 0, tu, eu);
    double ctd=conf-> GetSlewingMapManager()-> CalcCValue(cid, seg, 1, td, ed);
    double ctm=0.5*(ctu+ctd);

    h2 = (TH2F*)gFile->Get("eu_tu"), h2-> Fill(eu, tu);
    h2 = (TH2F*)gFile->Get("ed_td"), h2-> Fill(ed, td);
    h2 = (TH2F*)gFile->Get("eu_ctu"), h2-> Fill(eu, ctu);
    h2 = (TH2F*)gFile->Get("ed_ctd"), h2-> Fill(ed, ctd);
    h2 = (TH2F*)gFile->Get("eu_ctm"), h2-> Fill(eu, ctm);
    h2 = (TH2F*)gFile->Get("ed_ctm"), h2-> Fill(ed, ctm);
  }
}

void SlewingTools::Fill(int cid, int seg, ConfMan *conf, TNtuple *tup, double mean, double thre)
{
  TH2F *h2;
  float f_eu, f_ed, f_tu, f_td;
  tup-> SetBranchAddress("eu", &f_eu);
  tup-> SetBranchAddress("ed", &f_ed);
  tup-> SetBranchAddress("tu", &f_tu);
  tup-> SetBranchAddress("td", &f_td);

  //  cout<<"===== SlewingTools::Fill Entries : "<<tup->GetEntries()<<" ====="<<endl;
  for( int ev=0; ev<tup->GetEntries(); ev++ ){
    tup-> GetEntry(ev);
    double eu=f_eu, ed=f_ed, tu=f_tu, td=f_td;
    double ctu=conf-> GetSlewingMapManager()-> CalcCValue(cid, seg, 0, tu, eu);
    double ctd=conf-> GetSlewingMapManager()-> CalcCValue(cid, seg, 1, td, ed);
    double ctm=0.5*(ctu+ctd);

    h2 = (TH2F*)gFile->Get("eu_tu"), h2-> Fill(eu, tu);
    h2 = (TH2F*)gFile->Get("ed_td"), h2-> Fill(ed, td);
    h2 = (TH2F*)gFile->Get("eu_ctu"), h2-> Fill(eu, ctu);
    h2 = (TH2F*)gFile->Get("ed_ctd"), h2-> Fill(ed, ctd);
    if( mean-thre<ctm && ctm<mean+thre ){
      h2 = (TH2F*)gFile->Get("eu_ctm"), h2-> Fill(eu, ctm);
      h2 = (TH2F*)gFile->Get("ed_ctm"), h2-> Fill(ed, ctm);
    }
  }
}

bool SlewingTools::Reset()
{
  TH2F *h2;
  h2 = (TH2F*)gFile->Get("eu_tu");
  if( h2 ){
    h2-> Reset();
    h2-> GetXaxis()-> SetRangeUser(0, 100);
    h2-> GetYaxis()-> SetRangeUser(-100, 100);
  }
  else return false;

  h2 = (TH2F*)gFile->Get("ed_td");
  if( h2 ){
    h2-> Reset();
    h2-> GetXaxis()-> SetRangeUser(0, 100);
    h2-> GetYaxis()-> SetRangeUser(-100, 100);
  }
  else return false;

  h2 = (TH2F*)gFile->Get("eu_ctu");
  if( h2 ){
    h2-> Reset();
    h2-> GetXaxis()-> SetRangeUser(0, 100);
    h2-> GetYaxis()-> SetRangeUser(-50, 50);
  }
  else return false;

  h2 = (TH2F*)gFile->Get("ed_ctd");
  if( h2 ){
    h2-> Reset();
    h2-> GetXaxis()-> SetRangeUser(0, 100);
    h2-> GetYaxis()-> SetRangeUser(-50, 50);
  }
  else return false;

  h2 = (TH2F*)gFile->Get("eu_ctm");
  if( h2 ){
    h2-> Reset();
    h2-> GetXaxis()-> SetRangeUser(0, 100);
    h2-> GetYaxis()-> SetRangeUser(-50, 50);
  }
  else return false;

  h2 = (TH2F*)gFile->Get("ed_ctm");
  if( h2 ){
    h2-> Reset();
    h2-> GetXaxis()-> SetRangeUser(0, 100);
    h2-> GetYaxis()-> SetRangeUser(-50, 50);
  }
  else return false;

  return true;
}

void SlewingTools::SetRange(TH2F *h2)
{
  h2-> GetXaxis()-> SetRangeUser(0, 100);
  h2-> GetYaxis()-> SetRangeUser(-100, 100);

  TH1F *h1_py=projectionY(h2);
  //  CalibTools::FitGaus_2sigma(h1_py);
  h1_py-> Fit("gaus", "q", "q", h1_py->GetBinCenter(h1_py->GetMaximumBin())-2.*h1_py->GetRMS(), h1_py->GetBinCenter(h1_py->GetMaximumBin())+2.*h1_py->GetRMS());
  double mean=h1_py->GetFunction("gaus")->GetParameter(1);
  //  double sigma=h1_py->GetFunction("gaus")->GetParameter(2);

  TH1F *h1_px=projectionX(h2);
  h1_px->Fit("landau", "0q", "0q");
  double mpv=h1_px->GetFunction("landau")->GetParameter(1);
  h2->GetXaxis()-> SetRangeUser(0, 5*mpv);
  h2->GetYaxis()-> SetRangeUser(mean-5, mean+5);
} 

TH1F* SlewingTools::projectionX(TH2F *h2)
{
  //  int id=0;
  for( int i=0; i<10000; i++ ){
    TH1F *h1=(TH1F*)gFile-> Get(Form("%s_px_%d", h2->GetName(), i));
    if( !h1 ) return (TH1F*)h2->ProjectionX(Form("%s_px_%d", h2->GetName(),i));
  }
  return 0;
}

TH1F* SlewingTools::projectionY(TH2F *h2)
{
  //  int id=0;
  for( int i=0; i<10000; i++ ){
    TH1F *h1=(TH1F*)gFile-> Get(Form("%s_py_%d", h2->GetName(), i));
    if( !h1 ) return (TH1F*)h2->ProjectionY(Form("%s_py_%d", h2->GetName(),i));
  }
  return 0;
}
