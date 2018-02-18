// Tools.h

#ifndef Tools_h
#define Tools_h 1

#include <string>
#include <iostream>
#include <cmath>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

namespace Tools
{
  inline void Fill1D(TString name,const double &val);
  inline void Fill2D(TString name,const double &val1, const double &val2);
  void newTH1F(TString name, const int &nbin,const double &lbin,const double &ubin);
  void newTH2F(TString name, 
	       const int &nbinx,const double &lbinx,const double &ubinx,
	       const int &nbiny,const double &lbiny,const double &ubiny);
};


inline void Tools::Fill1D(TString name,const double &val){
  TH1F* h1;
  //  gFile->GetObject(name.Data(),h1);
  h1=(TH1F*)gFile->Get(name.Data());
  if(h1){
    h1->Fill(val);
  }else{
    std::cout<<name<<" does not exist !!!!!"<<std::endl;
  }
}

inline void Tools::Fill2D(TString name,const double &val1,const double &val2){
  TH2F* h1;
  //  gFile->GetObject(name.Data(),h1);
  h1=(TH2F*)gFile->Get(name.Data());
  if(h1){
    h1->Fill(val1,val2);
  }else{
    std::cout<<name<<" does not exist !!!!!"<<std::endl;
  }
}

void Tools::newTH1F(TString name, 
			   const int &nbinx,const double &lbinx,const double &ubinx)
{
  new TH1F(name,name,
	   nbinx,lbinx,ubinx);
}

void Tools::newTH2F(TString name, 
			   const int &nbinx,const double &lbinx,const double &ubinx,
			   const int &nbiny,const double &lbiny,const double &ubiny)
{
  new TH2F(name,name,
	   nbinx,lbinx,ubinx,
	   nbiny,lbiny,ubiny); 
}

#endif
