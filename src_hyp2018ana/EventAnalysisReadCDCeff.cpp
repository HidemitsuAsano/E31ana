#include "EventAnalysisReadCDCeff.h"

using namespace std;

bool EventAnalysisReadCDCeff::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;
  if( !MyAnaTools::beamFiducial(anaInfo, confMan) ) return true;
  vector<HodoscopeLikeHit*> CDHhit=MyTools::getCDH(cdsMan);
  vector<HodoscopeLikeHit*> IHhit=MyTools::getIH(cdsMan);

  if( CDHhit.size()!=1 || IHhit.size()!= 1 ) return true;

  double phiCDH=CDHhit[0]->pos().Phi();
  double CDHdE=CDHhit[0]-> emean();
  double CDHtof=CDHhit[0]-> ctmean()-anaInfo->beam(0)->T0time();

  double phiIH=IHhit[0]->pos().Phi();
  double IHdE=IHhit[0]-> eu();
  double IHtof=IHhit[0]-> tu()-anaInfo->beam(0)->T0time();

  double dphi=phiCDH-phiIH;
  if( dphi>TMath::Pi() ){ dphi -= 2.*TMath::Pi(); }
  else if( dphi<-TMath::Pi() ){ dphi += 2.*TMath::Pi(); }
  // if( IHhit[0]->seg()==24 ){
  //   cout<<"CDH seg : "<<CDHhit[0]-> seg()<<endl;
  //   cout<<"IH seg : "<<IHhit[0]-> seg()<<endl;
  //   cout<<"phiCDH : "<<phiCDH*TMath::RadToDeg()<<endl;
  //   cout<<"phiIH : "<<phiIH*TMath::RadToDeg()<<endl;
  //   cout<<"dphi : "<<dphi*TMath::RadToDeg()<<endl;
  // }

  if( fabs(dphi*TMath::RadToDeg())>20 ) return true;
  //  if( IHhit[0]->seg()==24 ) cout<<" Trigger true"<<endl;
  MyHistTools::fillTH(Form("CDH_trig_IH%d", IHhit[0]->seg()), CDHhit[0]->seg());

  MyHistTools::fillTH("CDH_dE", CDHdE);
  if( 5.5<CDHdE && CDHdE<9.0 ) MyHistTools::fillTH("trigCDH_dE", CDHdE);
  MyHistTools::fillTH("CDH_tof", CDHtof);
  if( 5.0<CDHtof && CDHtof<9.5 ) MyHistTools::fillTH("trigCDH_tof", CDHtof);
  MyHistTools::fillTH("IH_dE", IHdE);
  if( 0.3<IHdE && IHdE<0.9 ) MyHistTools::fillTH("trigIH_dE", IHdE);

  MyHistTools::fillTH("IH_tof", IHtof);
  MyHistTools::fillTH(Form("IH_dE_%d",IHhit[0]->seg()), IHdE);
  MyHistTools::fillTH(Form("IH_tof_%d",IHhit[0]->seg()), IHtof);

  if( CDHdE<5.5 || 9<CDHdE ) return true;
  if( CDHtof<5  || 9.5<CDHtof ) return true;
  if( IHdE<0.3  || 0.9<IHdE ) return true;

  MyHistTools::fillTH("trigIH_tof", IHtof);

  //    cout<<"CDC trigger"<<endl;
  MyHistTools::fillTH("trigCDC", CDHhit[0]->seg());
  MyHistTools::fillTH("trigCDC_IH", IHhit[0]->seg());
  double cos=CDHhit[0]->hitpos()/54.4;
  MyHistTools::fillTH("trigCDC_cos", cos);

  if( cdsFileMatching() ){
    bool CDCeff=false;
    for( int i=0; i<cdstrackMan->nGoodTrack(); i++ ){
      CDSTrack *track = cdstrackMan->GoodTrack(i);

      bool CDH_flag=false;
      TVector3 vtxCDH=track->CDHVertex();
      double dphiCDH=vtxCDH.Phi()-phiCDH;
      if( dphiCDH>TMath::Pi() ) dphiCDH -= 2.*TMath::Pi();
      else if( dphiCDH<-TMath::Pi() ) dphiCDH += 2.*TMath::Pi();
      if( fabs(dphiCDH)<10.*2.*TMath::Pi()/360. ) CDH_flag=true;

      bool IH_flag=false;
      TVector3 vtxIH=track->IHVertex();
      double dphiIH=vtxIH.Phi()-phiIH;
      if( dphiIH>TMath::Pi() ) dphiIH -= 2.*TMath::Pi();
      else if( dphiIH<-TMath::Pi() ) dphiIH += 2.*TMath::Pi();
      if( fabs(dphiIH)<15.*2.*TMath::Pi()/360. ) IH_flag=true;

      if( IH_flag && CDH_flag ) CDCeff=true;
    }

    if( CDCeff ){
      //      cout<<"   CDC effective"<<endl;
      MyHistTools::fillTH("effCDC", CDHhit[0]->seg());
      MyHistTools::fillTH("effCDC_IH", IHhit[0]->seg());
      MyHistTools::fillTH("effCDC_cos", cos);
    }
  }
  return true;
}

void EventAnalysisReadCDCeff::InitializeHistogram()
{
  new TH1F("CDH_dE", "CDH dE", 1000, 0.0, 50);
  new TH1F("IH_dE",  "IH dE", 1000, -10, 10);
  new TH1F("CDH_tof", "CDH TOF", 2000, -100, 100);
  new TH1F("IH_tof",  "IH TOF", 2000, -100, 100);

  new TH1F("trigCDH_dE", "CDH dE", 1000, 0.0, 50);
  new TH1F("trigIH_dE",  "IH dE", 1000, -10, 10);
  new TH1F("trigCDH_tof", "CDH TOF", 2000, -100, 100);
  new TH1F("trigIH_tof",  "IH TOF", 2000, -100, 100);


  for( int seg=1; seg<=24; seg++ ){
    new TH1F(Form("IH_dE_%d", seg),  Form("IH dE seg%d",seg), 1000, -10, 10);
    new TH1F(Form("IH_tof_%d", seg), Form("IH TOF seg%d",seg), 2000, -100, 100);
    new TH1F(Form("CDH_trig_IH%d", seg), Form("CDH trig by IH seg%d", seg), 36, 0.5, 36.5);
  }

  new TH1F("trigCDC", "CDC trigger",   36, 0.5, 36.5);
  new TH1F("effCDC",  "CDC effective", 36, 0.5, 36.5);

  new TH1F("trigCDC_IH", "CDC trigger",   36, 0.5, 36.5);
  new TH1F("effCDC_IH",  "CDC effective", 36, 0.5, 36.5);

  new TH1F("trigCDC_cos", "CDC trigger",   36, 0.5, 36.5);
  new TH1F("effCDC_cos",  "CDC effective", 36, 0.5, 36.5);
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadCDCeff *event = new EventAnalysisReadCDCeff();
  return (EventTemp*)event;
}
