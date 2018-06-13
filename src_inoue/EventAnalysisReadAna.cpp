#include "EventAnalysisReadAna.h"

using namespace std;

bool EventAnalysisReadAna::UAna()
{
  //event by event QA 
  //beam info.
  int nBHD = blMan->nBHD();
  Tools::H1("h1_nBHD",nBHD,20,0,20,"BHD multiplicity");
  int nT0  = blMan->nT0();
  Tools::H1("h1_nT0",nT0,5,0,5,"T0 multiplicity");
  int nbeam = anaInfo->nBeam();
  Tools::H1("h1_nbeam",nbeam,10,0,10,"nbeam");
  int nBLC1 = anaInfo->beam(0)->nBLC1();//# of track
  Tools::H1("h1_nBLC1",nBLC1,10,0,10,"# of BLC1 tracks");
  int nBLC2 = anaInfo->beam(0)->nBLC2();//# of track
  Tools::H1("h1_nBLC2",nBLC2,10,0,10,"# of BLC2 tracks");
  int nBPC = anaInfo->beam(0)->nBPC();//# of track
  Tools::H1("h1_nBPC",nBPC,10,0,10,"# of BPC tracks");
  int nFDC1 = anaInfo->beam(0)->nFDC1();//# of track
  Tools::H1("h1_nFDC1",nFDC1,10,0,10,"# of FDC1 tracks");
  double BLC1time = anaInfo->beam(0)->BLC1time();//GetTrackTime
  Tools::H1("h1_BLC1time",BLC1time,600,-200,400,"BLC1 track time [nsec.]");
  double BLC2time = anaInfo->beam(0)->BLC2time();
  Tools::H1("h1_BLC2time",BLC2time,600,-200,400,"BLC2 track time [nsec.]");
  double BPCtime = anaInfo->beam(0)->BPCtime();
  Tools::H1("h1_BPCtime",BPCtime,600,-200,400,"BPC track time [nsec.]");
  double BLC1chi2 = anaInfo->beam(0)->BLC1chi2();
  Tools::H1("h1_BLC1chi2",BLC1chi2,100,0,100,"BLC1 chi2/ndf");
  double BLC2chi2 = anaInfo->beam(0)->BLC2chi2();
  Tools::H1("h1_BLC2chi2",BLC2chi2,100,0,100,"BLC2 chi2/ndf");
  double BPCchi2 = anaInfo->beam(0)->BPCchi2();
  Tools::H1("h1_BPCchi2",BPCchi2,100,0,100,"BPC chi2/ndf");
  double t0BHDTOF = anaInfo->beam(0)->tof();//BHD - T0 TOF
  Tools::H1("h1_T0BHDtof",t0BHDTOF,200,-100,100,"TOF (BHD - T0) [nsec.]");
  double D5specchi2 = anaInfo->beam(0)->D5chi2();//chi2/ndf
  Tools::H1("h1_D5chi2",D5specchi2,100,0,100,"beam reco. chi2/ndf");
  double D5mom = anaInfo->beam(0)->D5mom();
  Tools::H1("h1_D5mom",D5mom,100,0.8,1.2,"reco. Mom. [GeV/c]");
  bool IsconnectD5BHD = MyAnaTools::connectD5BHD(anaInfo);
  Tools::H1("h1_IsconnectD5BHD",IsconnectD5BHD,2,-0.5,1.5,"is connect D5BHD");
  bool IsconnectBLC2BPC = MyAnaTools::connectBLC2BPC(anaInfo);
  Tools::H1("h1_IsconnectBLC2BPC",IsconnectBLC2BPC,2,-0.5,1.5,"is connect BLC2BPC");
  //endQA
  
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    anaInfo-> ClearFC();
    ForwardChargeInfo fcInfo=MyTools::makeFC(confMan, blMan, anaInfo);
    //    fcInfo.dump();
    if( fcInfo.step()>1 ){
      anaInfo->AddCharge(fcInfo);
    }
    // for( int i=0; i<anaInfo->nCDS(); i++ ){
    //   CDSTrack *track=anaInfo->CDS(i)->track(confMan, cdsMan, cdstrackMan);
    //   // TVector3 vtxCDH=track->CDHVertex();
    //   // if( fabs(vtxCDH.Z())>40 ) anaInfo->CDS(i)->SetFlag(false);
    // }
    
    fillHistReadCDS(anaInfo, cdsMan, cdstrackMan);
    fillHistReadNC(header, blMan, anaInfo);
    fillHistReadFC(header, confMan, blMan, anaInfo);
    fillHistReadKNpim(header, blMan, anaInfo);
    fillHistReadKNpip(header, blMan, anaInfo);
    fillHistReadKNkm(header, blMan, anaInfo);
    static int first=0;
    if( MyAnaTools::targetNA()==1 ){
      if(first==0){
        std::cout << __FILE__ << " LHydrogen run is analyzed" << std::endl;
        first++;
      }
      fillHistReadKNpipi_H2(header, confMan, blMan, cdsMan, cdstrackMan, anaInfo);
    }
    if( MyAnaTools::targetNA()==2 ){
      if(first==0){
        std::cout << __FILE__ << " LDeuterium run is analyzed" << std::endl;
        first++;
      }
      fillHistReadKNpipi_D2(header, blMan, cdsMan, cdstrackMan, anaInfo);
      fillHistReadKPpimpim_D2(header, confMan, blMan, cdsMan, cdstrackMan, anaInfo);
      fillHistReadCDS_Lpim_D2(header, blMan, anaInfo);

      //      fillHistReadVtxDEF(header, confMan, anaInfo, blMan, cdsMan, cdstrackMan);
    }
    else{
      if(first==0){
        std::cout << __FILE__ << " target is not chosen" << std::endl;
        first++;
      }
      //      fillHistReadKNpipi(header, blMan, cdsMan, cdstrackMan, anaInfo);
    }
  }
  return true;
}

void EventAnalysisReadAna::InitializeHistogram()
{
  TTree *tree=new TTree("npipi_event", "npipi_event");
  tree-> Branch("AnaInfo", &anaInfo);

  TTree *tree2=new TTree("ppimpim_event", "ppimpim_event");
  tree2-> Branch("AnaInfo", &anaInfo);

  initHistReadCDS();
  initHistReadNC();
  initHistReadFC();
  initHistReadKNpim();
  initHistReadKNpip();
  initHistReadKNkm();
  if( MyAnaTools::targetNA()==1 ){
    initHistReadKNpipi_H2();
  }
  if( MyAnaTools::targetNA()==2 ){
    initHistReadKPpimpim_D2();
    initHistReadKNpipi_D2(header, anaInfo);
    initHistReadCDS_Lpim_D2();

    //    initHistReadVtxDEF();
  }
  else{
    //    initHistReadKNpipi();
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadAna *event = new EventAnalysisReadAna();
  return (EventTemp*)event;
}
