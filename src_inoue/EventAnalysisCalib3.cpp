#include "EventAnalysisCalib3.h"

using namespace std;

void EventAnalysisCalib3::anaNC()
{
  ForwardNeutralInfo info0=MyTools::makeFN(blMan, anaInfo, 0.0);
  fillSlewingNC(blMan, anaInfo, &info0);

  ForwardNeutralInfo info=MyTools::makeFN(blMan, anaInfo, 8.0);
  if( info.pid()==F_Neutron ) anaInfo-> AddNeutral(info);
  else if( info.pid()==F_Gamma ){
    anaInfo-> AddNeutral(info);
    fillT0NC(blMan, anaInfo);
  }
}

void EventAnalysisCalib3::anaCDS()
{
  //  fillCDC(cdstrackMan, cdsMan);
  BeamInfo *beam=anaInfo->beam(0);

  if( beam->flag() ){
    for( int i=0; i<cdstrackMan->nTrack(); i++ ){
      CDSTrack *track=cdstrackMan->Track(i);
      CDSInfo cds=MyTools::makeCDSInfo(i, cdsMan, cdstrackMan, beam, bltrackMan);
      fillCDS(&cds, cdsMan, cdstrackMan);
      anaInfo-> AddCDS(cds);
    }
    if( anaInfo->minDCA() ){ beam-> SetVertex(anaInfo->minDCA()->vertexBeam()); }

    if( anaInfo->minDCA() && GeomTools::GetID(anaInfo->minDCA()->vertexBeam())==CID_Fiducial ){
      fillCDH(beam, anaInfo->minDCA(), cdsMan);
    }

    for( int i=0; i<cdstrackMan->nTrack(); i++ ){
      for( int j=i+1; j<cdstrackMan->nTrack(); j++ ){
	CDS2Info cds2=MyTools::makeCDS2Info(cdstrackMan, i, j, anaInfo);
	anaInfo->AddCDS2(cds2);
      }
    }
  }
}

void EventAnalysisCalib3::anaBL1()
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
  fillFDC1(header, blMan, bltrackMan);
  fillBLC2BPC(bltrackMan);
  fillBeamSpectrometer(beamSpec);

  fillBeamProf(header, bltrackMan);

  fillT0DEF(blMan, bltrackMan, anaInfo);
  fillT0BVC(blMan, bltrackMan, anaInfo);
}

void EventAnalysisCalib3::anaBL0()
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
  //  fillCDC(cdsMan);

  fillBHDT0(header, blMan);
}

void EventAnalysisCalib3::InitializeHistogram()
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
  //  initHistCDC(confMan);

  initHistCDS();
  initHistBHDT0();
  initHistT0DEF();
  initHistT0BVC();
  initHistT0NC();

  initHistReduction();
  initHistBeamAna();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisCalib3 *event = new EventAnalysisCalib3();
  return (EventTemp*)event;
}
