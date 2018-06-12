#include "MyHistMCChkCou.h"

void initHistMCChkCou(){
  new TH1F("T0_dE_MC",  "T0_dE_MC",  1000, 0.0, 100);

  new TH1F("BPD_dE_MC", "BPD_dE_MC", 1000, 0.0, 100);

  new TH1F("DEF_dE_MC", "DEF_dE_MC", 1000, 0.0, 100);

  new TH1F("CDH_dE_MC", "CDH_dE_MC", 1000, 0.0, 100);

  new TH1F("BVC_dE_MC", "BVC_dE_MC", 1000, 0.0, 100);

  new TH1F("CVC_dE_MC", "CVC_dE_MC", 1000, 0.0, 100);

  new TH1F("PC_dE_MC",  "PC_dE_MC",  1000, 0.0, 100);

  new TH1F("NC_dE_MC",  "NC_dE_MC",  1000, 0.0, 100);
}

void fillHistMCChkCou(EventHeaderMC *headerMC, DetectorData *detData, MCData *mcDatax){
  for( int i=0; i<detData->detectorHitSize(); i++ ){
    DetectorHit *hit=detData->detectorHit(i);
    if( hit->detectorID()==CID_T0  ) MyHistTools::fillTH("T0_dE_MC",  hit->adc());
    if( hit->detectorID()==CID_BPD ) MyHistTools::fillTH("BPD_dE_MC", hit->adc());
    if( hit->detectorID()==CID_DEF ) MyHistTools::fillTH("DEF_dE_MC", hit->adc());
    if( hit->detectorID()==CID_CDH ) MyHistTools::fillTH("CDH_dE_MC", hit->adc());
    if( hit->detectorID()==CID_BVC ) MyHistTools::fillTH("BVC_dE_MC", hit->adc());
    if( hit->detectorID()==CID_CVC ) MyHistTools::fillTH("CVC_dE_MC", hit->adc());
    if( hit->detectorID()==CID_PC  ) MyHistTools::fillTH("PC_dE_MC",  hit->adc());
    if( hit->detectorID()==CID_NC  ) MyHistTools::fillTH("NC_dE_MC", hit->adc());
  }
}

