#include "NppimData.h"

NppimData::NppimData()
{
  clear();
}

void NppimData::clear()
{
  fBeamLmom.SetXYZT(0.0, 0.0, 0.0, -999.);
  fFNLmom.SetXYZT(0.0, 0.0, 0.0, -999.);
  fPimLmom.SetXYZT(0.0, 0.0, 0.0, -999.);
  fPLmom.SetXYZT(0.0, 0.0, 0.0, -999.);

  fVtxBeam = DEFVECT;
  fVtxCDS = DEFVECT;
  fVtxPim = DEFVECT;
  fVtxP = DEFVECT;
}
