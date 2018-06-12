#include "EventAnalysisReadBeam.h"

using namespace std;

bool EventAnalysisReadBeam::UAna()
{
  if( anaInfo->nBeam()!=1 ) return false;
  fillHistT0BVC(blMan, anaInfo);
  fillHistT0DEF(blMan, anaInfo);

  MyHistTools::fillTH("EventReduction", 0);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 0);
  //  cout<<" Reduction : "<<0<<endl;

  fillReadBeam_T0(header, blMan);
  if( MyTools::getHodo(blMan, CID_T0).size()==1 ){
    MyHistTools::fillTH("EventReduction", 1);
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 1);
    //    cout<<" Reduction : "<<1<<endl;
  }

  fillReadBeam_BHDT0(header, anaInfo);
  if( !MyAnaTools::isTOFKaon(anaInfo->beam(0)->tof()) ) return true;
  MyHistTools::fillTH("EventReduction", 2);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 2);
  //  cout<<" Reduction : "<<2<<endl;

  fillReadBeam_BLC1(header, confMan, anaInfo);
  if( !MyAnaTools::anaBLC1(anaInfo) ) return true;
  MyHistTools::fillTH("EventReduction", 3);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 3);
  //  cout<<" Reduction : "<<3<<endl;

  fillReadBeam_BLC2(header, confMan, anaInfo);
  if( !MyAnaTools::anaBLC2(anaInfo) ) return true;
  MyHistTools::fillTH("EventReduction", 4);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 4);
  //  cout<<" Reduction : "<<4<<endl;

  fillReadBeam_D5(header, anaInfo);
  if( !MyAnaTools::anaD5(anaInfo) ) return true;
  MyHistTools::fillTH("EventReduction", 5);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 5);
  //  cout<<" Reduction : "<<5<<endl;

  fillReadBeam_D5BHD(header, anaInfo);
  if( confMan->GetOutFileName().find("Run49c")==string::npos ){
    if( !MyAnaTools::connectD5BHD(anaInfo) ) return true;
  }
  MyHistTools::fillTH("EventReduction", 6);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 6);
  //  cout<<" Reduction : "<<6<<endl;

  fillReadBeam_BPC(header, confMan, anaInfo, blMan);
  if( !MyAnaTools::anaBPC(anaInfo) ) return true;
  MyHistTools::fillTH("EventReduction", 7);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 7);
  //  cout<<" Reduction : "<<7<<endl;

  fillReadBeam_BLC2BPC(header, anaInfo);
  if( !MyAnaTools::connectBLC2BPC(anaInfo) ) return true;
  MyHistTools::fillTH("EventReduction", 8);
  if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 8);
  //  cout<<" Reduction : "<<8<<endl;

  fillReadBeam_Profile(confMan, header, anaInfo);
  if( MyAnaTools::beamFiducial(anaInfo, confMan) ){
    MyHistTools::fillTH("EventReduction", 9);
    if( header->IsTrig(Trig_Kf) ) MyHistTools::fillTH("Kf_Reduction", 9);
    //    cout<<" Reduction : "<<9<<endl;
  }

  fillReadBeam_FDC1(header, confMan, anaInfo, blMan);

  if( cdsFileMatching() ){
  }

  return true;
}

void EventAnalysisReadBeam::InitializeHistogram()
{
  initHistT0BVC();
  initHistT0DEF();
  initHistReadBeam();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadBeam *event = new EventAnalysisReadBeam();
  return (EventTemp*)event;
}
