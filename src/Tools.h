// Tools.h

#ifndef Tools_h
#define Tools_h 1

#include <string>
#include <iostream>
#include <cmath>

#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH2I.h>
#include <TFile.h>
#include "GlobalVariables.h"
#include <TMacro.h>
#include <TSystem.h>

namespace Tools
{
  inline void Fill1D(TString name, double val);
  inline void Fill2D(TString name, double val1,  double val2);
  inline void Fill3D(TString name, double val1,  double val2 ,double val3);
  inline void newTH1F(TString name,  int nbinx, double lbinx, double ubinx);
  inline void newTH1F(TString name, TString title, int nbin, double lbin, double ubin);
  inline void newTH2F(TString name, 
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny);
  inline void newTH3F(TString name, 
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny,
		      int nbinz, double lbinz, double ubinz);
  inline void H1(TString name,  double val,
		 int nbin, double lbin, double ubin);
  inline void H1(TString name,  double val,
		 int nbin, double lbin, double ubin, TString xtitle);
  inline void H1(TString name,  double val,
		 int nbin, double lbin, double ubin, double weight);

  inline void H2(TString name,  double val1, double val2,
		 int nbinx, double lbinx, double ubinx,
		 int nbiny, double lbiny, double ubiny);
  inline void H2(TString name,  double val1, double val2,
		 int nbinx, double lbinx, double ubinx,
		 int nbiny, double lbiny, double ubiny,TString xtitle, TString ytitle);
  inline void H2I(TString name,  double val1, double val2,
		 int nbinx, double lbinx, double ubinx,
		 int nbiny, double lbiny, double ubiny);
  inline void WriteFile(TString afile,TString direntry);
  inline void WriteFile(TString afile);
  inline void SaveFile(TString filename);
  
  inline void SetXTitleH1(TString histname,TString xtitle);
  inline void SetYTitleH1(TString histname,TString ytitle);
  inline void SetXTitleH2(TString histname,TString xtitle);
  inline void SetYTitleH2(TString histname,TString ytitle);
};

inline void Tools::WriteFile(TString afile,TString direntry){
  if(afile==DefaultFileName) return;
  //  std::cout<<afile<<std::endl;
  TMacro *m = new TMacro(afile);
  m->Write(direntry);
  delete m;
  return;
}

inline void Tools::WriteFile(TString afile){
  if(afile==DefaultFileName) return;
  //  std::cout<<afile<<std::endl;
  TMacro *m = new TMacro(afile);
  m->Write(afile);
  delete m;
  return;
}

inline void Tools::SaveFile(TString filename){
  std::cout<<filename<<std::endl;
  const char *dirname="file";
  char *slash = (char*)strrchr(dirname,'/');
  char *locdir;
  if (slash) locdir = slash+1;
  else       locdir = (char*)dirname;
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->GetDirectory(locdir);
  if(!adir) adir= savdir->mkdir(locdir);
  if(adir)  adir->cd();
  else{
    std::cout<<"!!!"<<std::endl;
    return;
  }
  //  gSystem->cd(locdir);
  //  void *dirp = gSystem->OpenDirectory(dirname);
  //  if (!dirp) return;
  //  TString afile = Form("%s/%s",dirname,filename.Data());
  Tools::WriteFile( filename);
  //  gSystem->FreeDirectory(dirp);
  savdir->cd();
  return;
}

inline void Tools::Fill1D(TString name, double val){
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(name.Data())) ){
    h1->Fill(val);
  }else{
    std::cout<<name<<" does not exist !!!!!"<<std::endl;
  }
}

inline void Tools::Fill2D(TString name, double val1, double val2){
  if( TH2F* h1 = dynamic_cast<TH2F *>(gFile->Get(name.Data())) ){
    h1->Fill(val1,val2);
  }else{
    std::cout<<name<<" does not exist !!!!!"<<std::endl;
  }
}

inline void Tools::Fill3D(TString name, double val1, double val2, double val3){
  if( TH3F* h1 = dynamic_cast<TH3F *>(gFile->Get(name.Data())) ){
    h1->Fill(val1,val2, val3);
  }else{
    std::cout<<name<<" does not exist !!!!!"<<std::endl;
  }
}

inline void Tools::newTH1F(TString name, 
			   int nbinx, double lbinx, double ubinx)
{
  new TH1F(name,name, nbinx,lbinx,ubinx);
}

inline void Tools::newTH1F(TString name, TString title,
			   int nbinx, double lbinx, double ubinx)
{
  new TH1F(name,title, nbinx,lbinx,ubinx);
}

inline void Tools::newTH2F(TString name, 
			   int nbinx, double lbinx, double ubinx,
			   int nbiny, double lbiny, double ubiny)
{
  new TH2F(name,name, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny); 
}

inline void Tools::newTH3F(TString name, 
			   int nbinx, double lbinx, double ubinx,
			   int nbiny, double lbiny, double ubiny,
         int nbinz, double lbinz, double ubinz)
{
  new TH3F(name,name, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny, nbinz,lbinz,ubinz); 
}


inline void Tools::H1(TString name, double val,
		      int nbinx, double lbinx, double ubinx)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(name.Data())) ){
    //    std::cout<<"exist "<<name<<std::endl;
    h1->Fill(val);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1 = new TH1F(name,name, nbinx,lbinx,ubinx);
    h1->Fill(val);
  }
}

inline void Tools::H1(TString name, double val,
		      int nbinx, double lbinx, double ubinx, TString xtitle)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(name.Data())) ){
    //    std::cout<<"exist "<<name<<std::endl;
    h1->Fill(val);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1 = new TH1F(name,name, nbinx,lbinx,ubinx);
    h1->SetXTitle(xtitle.Data());
    h1->Fill(val);
  }
}

inline void Tools::H1(TString name, double val,
		      int nbinx, double lbinx, double ubinx,
		      double weight)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(name.Data())) ){
    //    std::cout<<"exist "<<name<<std::endl;
    h1->Fill(val,weight);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1 = new TH1F(name,name, nbinx,lbinx,ubinx);
    h1->Fill(val,weight);
  }
}

inline void Tools::H2(TString name, double val1, double val2,
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny)
{
  if( TH2F* h1 = dynamic_cast<TH2F *>(gFile->Get(name.Data())) ){
    h1->Fill(val1,val2);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1= new TH2F(name,name, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny); 
    h1->Fill(val1,val2);
  }  
}

inline void Tools::H2(TString name, double val1, double val2,
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny,TString xtitle,TString ytitle)
{
  if( TH2F* h1 = dynamic_cast<TH2F *>(gFile->Get(name.Data())) ){
    h1->Fill(val1,val2);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1= new TH2F(name,name, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny); 
    h1->SetXTitle(xtitle.Data());
    h1->SetYTitle(ytitle.Data());
    h1->Fill(val1,val2);
  }  
}

inline void Tools::H2I(TString name, double val1, double val2,
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny)
{
  if( TH2I* h1 = dynamic_cast<TH2I *>(gFile->Get(name.Data())) ){
    h1->Fill(val1,val2);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1= new TH2I(name,name, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny); 
    h1->Fill(val1,val2);
  }  
}

inline void Tools::SetXTitleH1(TString histname,TString xtitle)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(histname.Data())) ){
    h1->SetXTitle(xtitle.Data());
  }else{
    std::cout << __FILE__ << " L. " << __LINE__  << " can not find histgram named   " << histname << std::endl;
  }
}

inline void Tools::SetYTitleH1(TString histname,TString xtitle)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(histname.Data())) ){
    h1->SetYTitle(xtitle.Data());
  }else{
    std::cout << __FILE__ << " L. " << __LINE__  << " can not find histgram named   " << histname << std::endl;
  }
}


inline void Tools::SetXTitleH2(TString histname,TString xtitle)
{
  if( TH2F* h2 = dynamic_cast<TH2F *>(gFile->Get(histname.Data())) ){
    h2->GetXaxis()->SetTitle(xtitle.Data());
  }else{
    std::cout << __FILE__ << " L. " << __LINE__  << " can not find histgram named   " << histname << std::endl;
  }
}

inline void Tools::SetYTitleH2(TString histname,TString ytitle)
{
  if( TH2F* h2 = dynamic_cast<TH2F *>(gFile->Get(histname.Data())) ){
    h2->GetYaxis()->SetTitle(ytitle.Data());
  }else{
    std::cout << __FILE__ << " L. " << __LINE__  << " can not find histgram named   " << histname << std::endl;
  }
}

#endif
