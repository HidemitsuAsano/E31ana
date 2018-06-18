#include "EventAnalysisChargedMode.h"

using namespace std;

bool EventAnalysisChargedMode::UAna()
{
  //event by event QA 
  //beam info.
  int nBHD = blMan->nBHD();//this returns # of segments, not hits
  int nBHDHit=0;
  for(int iseg=0;iseg<nBHD;iseg++){
    HodoscopeLikeHit *hit = blMan-> BHD(iseg);
    if( hit-> CheckRange() ) nBHDHit++;
  }
  Tools::H1("h1_nBHD",nBHDHit,20,0,20,"BHD multiplicity");
  int nT0  = blMan->nT0();
  int nT0Hit=0;
  for(int iseg=0;iseg<nT0;iseg++){
    HodoscopeLikeHit *hit = blMan-> T0(iseg);
    if( hit-> CheckRange() ) nT0Hit++;
  }
  Tools::H1("h1_nT0",nT0Hit,5,0,5,"T0 multiplicity");
  double t0BHDTOF = anaInfo->beam(0)->tof();//BHD - T0 TOF
  if(nT0Hit==1)Tools::H1("h1_T0BHDtof",t0BHDTOF,2000,-100,100,"TOF (BHD - T0) [nsec.]");
  
  if(nT0Hit!=1) return true;
  //end of QA before Kaon ID
  if(!MyAnaTools::isTOFKaon(t0BHDTOF)) return true;
  Tools::H1("h1_T0BHDtof_cut",t0BHDTOF,2000,-100,100,"TOF (BHD - T0) [nsec.]");
  
  int nBLC1 = anaInfo->beam(0)->nBLC1();//# of track
  Tools::H1("h1_nBLC1",nBLC1,10,0,10,"# of BLC1 tracks");
  double BLC1time = anaInfo->beam(0)->BLC1time();//GetTrackTime
  Tools::H1("h1_BLC1time",BLC1time,130,-30,100,"BLC1 track time [nsec.]");
  double BLC1chi2 = anaInfo->beam(0)->BLC1chi2();
  if(nBLC1==1 && (-10 < BLC1time) && (BLC1time < 10)){
    Tools::H1("h1_BLC1chi2",BLC1chi2,100,0,100,"BLC1 chi2/ndf");
  }

  int nBLC2 = anaInfo->beam(0)->nBLC2();//# of track
  double BLC2time = anaInfo->beam(0)->BLC2time();
  double BLC2chi2 = anaInfo->beam(0)->BLC2chi2();
  if(MyAnaTools::anaBLC1(anaInfo)){
    Tools::H1("h1_nBLC2",nBLC2,10,0,10,"# of BLC2 tracks");
    if(nBLC2==1){
      Tools::H1("h1_BLC2time",BLC2time,130,-30,100,"BLC2 track time [nsec.]");
    }
    if(nBLC2==1 && (-10 < BLC2time) && (BLC2time < 10)){
      Tools::H1("h1_BLC2chi2",BLC2chi2,100,0,100,"BLC2 chi2/ndf");
    }
  }
  //BLC1 and BLC2 selection
  if( !MyAnaTools::anaBLC1(anaInfo) ||  !MyAnaTools::anaBLC2(anaInfo)) return true;
  
  int nBPC = anaInfo->beam(0)->nBPC();//# of track
  double BPCtime = anaInfo->beam(0)->BPCtime();
  double BPCchi2 = anaInfo->beam(0)->BPCchi2();
  
  Tools::H1("h1_nBPC",nBPC,10,0,10,"# of BPC tracks");
  if(nBPC==1){
    Tools::H1("h1_BPCtime",BPCtime,130,-30,100,"BPC track time [nsec.]");
    if((-10 < BPCtime) && (BPCtime < 10)) Tools::H1("h1_BPCchi2",BPCchi2,100,0,100,"BPC chi2/ndf");
  }
  if(!MyAnaTools::anaBPC(anaInfo)) return true;

  int nFDC1 = anaInfo->beam(0)->nFDC1();//# of track
  Tools::H1("h1_nFDC1",nFDC1,10,0,10,"# of FDC1 tracks");
  double D5specchi2 = anaInfo->beam(0)->D5chi2();//chi2/ndf
  Tools::H1("h1_D5chi2",D5specchi2,100,0,100,"beam reco. chi2/ndf");

  double D5mom = anaInfo->beam(0)->D5mom();
  if(D5specchi2 < 30) Tools::H1("h1_D5mom",D5mom,100,0.8,1.2,"reco. Mom. [GeV/c]");
  
  //bool IsconnectD5BHD = MyAnaTools::connectD5BHD(anaInfo);
  //Tools::H1("h1_IsconnectD5BHD",IsconnectD5BHD,2,-0.5,1.5,"is connect D5BHD");
  bool IsconnectBLC2BPC = MyAnaTools::connectBLC2BPC(anaInfo);
  Tools::H1("h1_IsconnectBLC2BPC",IsconnectBLC2BPC,2,-0.5,1.5,"is connect BLC2BPC");
  
  //BLC2BPC matching cut QA
  BeamInfo *beam = anaInfo->beam(0);
  if( anaInfo->nBeam()==1 &&
      MyAnaTools::anaBLC1(anaInfo) &&
      MyAnaTools::anaBLC2(anaInfo) && 
      MyAnaTools::anaBPC(anaInfo) && 
      MyAnaTools::anaD5(anaInfo) ){

    BLDCTrackInfo infoBLC2=beam->BLC2(0);
    BLDCTrackInfo infoBPC=beam->BPC(0);

    double z = 0.5*(-130-20.3);//this number should be checked
    TVector3 dir_diff=infoBLC2.dir()-infoBPC.dir();
    TVector3 pos_diff=infoBLC2.GetPosatZ(z)-infoBPC.GetPosatZ(z);
    Tools::H2("h2_BLC2_BPC_posdiff",pos_diff.X(),pos_diff.Y(),1000,-2,2,1000,-2,2.);
    Tools::H2("h2_BLC2_BPC_dirdiff",dir_diff.X(),dir_diff.Y(),1000,-0.5,0.5,1000,-0.5,0.5);
  }

  //endQA beam 
  
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;
  
  //QA if good Beam
  

  
  //endQA

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

void EventAnalysisChargedMode::InitializeHistogram()
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
  EventAnalysisChargedMode *event = new EventAnalysisChargedMode();
  return (EventTemp*)event;
}
