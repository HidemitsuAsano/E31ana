#include "EventAnalysisCalib2.h"

using namespace std;

void EventAnalysisCalib2::anaCDS()
{
  //  cout<<"Event Number : "<<Event_Number<<endl;
  fillCDC(cdstrackMan, cdsMan);
  BeamInfo beam=MyTools::makeBeamInfo(confMan, header, blMan, bltrackMan, beamSpec);

  if( beam.flag() ){
    fillBHDT0offset(&beam, blMan);
    //    cout<<"Event Number : "<<Event_Number<<endl;
    for( int i=0; i<cdstrackMan->nTrack(); i++ ){
      CDSTrack *track=cdstrackMan->Track(i);
      CDSInfo cds=MyTools::makeCDSInfo(i, cdsMan, cdstrackMan, &beam, bltrackMan);
      //      cds.dump();
      fillCDH(&beam, &cds, cdsMan);
      fillCDS(&cds, cdsMan, cdstrackMan);
    }
  }
}

void EventAnalysisCalib2::anaBL1()
{
  if( bltrackMan->ntrackBLC1a()==1 ){ fillBLC1a(bltrackMan->trackBLC1a(0)); }
  if( bltrackMan->ntrackBLC1b()==1 ){ fillBLC1b(bltrackMan->trackBLC1b(0)); }
  if( bltrackMan->ntrackBLC2a()==1 ){ fillBLC2a(bltrackMan->trackBLC2a(0)); }
  if( bltrackMan->ntrackBLC2b()==1 ){ fillBLC2b(bltrackMan->trackBLC2b(0)); }
  if( bltrackMan->ntrackBPC()==1   ){ fillBPC(bltrackMan->trackBPC(0)); }
  if( bltrackMan->ntrackFDC1()==1  ){ fillFDC1(bltrackMan->trackFDC1(0)); }

  fillReduction(blMan, bltrackMan, beamSpec, header);

  fillBPC(bltrackMan);
  fillBLC1(bltrackMan);
  fillBLC2(bltrackMan);
  fillBLC1(bltrackMan);
  fillBLC2BPC(bltrackMan);
  fillBeamSpectrometer(beamSpec);

  fillBeamProf(header, bltrackMan);
}

void EventAnalysisCalib2::anaBL0()
{
  fillT0_ADC(blMan);
  fillBHD_ADC(blMan);
  fillCVC_ADC(blMan);
  fillPC_ADC(blMan);
  fillNC_ADC(blMan);
  fillBPD_ADC(blMan);
  fillDEF_ADC(blMan);
  fillBVC_ADC(blMan);

  fillBLC1a(blMan);
  fillBLC1b(blMan);
  fillBLC2a(blMan);
  fillBLC2b(blMan);
  fillBPC(blMan);
  fillFDC1(blMan);

  fillCDH_ADC(cdsMan);
  fillCDC(cdsMan);

  fillBHDT0(header, blMan);
}

void EventAnalysisCalib2::InitializeHistogram()
{
  rtFile-> cd();
  initHistT0();
  initHistBHD();
  initHistCVC();
  initHistPC();
  initHistNC();
  initHistBPD();
  initHistDEF();
  initHistBVC();
  initHistCDH();

  initHistBLC1();
  initHistBLC1a();
  initHistBLC1b();

  initHistBLC2();
  initHistBLC2a();
  initHistBLC2b();

  initHistBPC();
  initHistFDC1();
  initHistCDC(confMan);

  initHistCDS();
  initHistBHDT0();

  initHistReduction();
  initHistBeamAna();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisCalib2 *event = new EventAnalysisCalib2();
  return (EventTemp*)event;
}
