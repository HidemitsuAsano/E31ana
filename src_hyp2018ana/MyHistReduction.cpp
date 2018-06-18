#include "MyHistReduction.h"

using namespace std;

void initHistReduction()
{
  new TH1F("Kf_Reduction", "K/f Event Reduction", 20, 0, 20);
}

void fillReduction(BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *spec, EventHeader *header)
{
  bool trigKf=header->IsTrig(Trig_Kf);
}
