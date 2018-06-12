#include "NpipiData.h"

NpipiData::NpipiData()
{
  clear();
}

void NpipiData::clear()
{
  fBeamLmom.SetXYZT(0.0, 0.0, 0.0, -999.);
  fFNLmom.SetXYZT(0.0, 0.0, 0.0, -999.);
  fPimLmom.SetXYZT(0.0, 0.0, 0.0, -999.);
  fPipLmom.SetXYZT(0.0, 0.0, 0.0, -999.);

  fVtxBeam = DEFVECT;
  fVtxCDS = DEFVECT;
  fVtxPim = DEFVECT;
  fVtxPip = DEFVECT;
  
  fFNLmomPim.SetXYZT(0.0, 0.0, 0.0, -999.);
  fFNLmomPip.SetXYZT(0.0, 0.0, 0.0, -999.);
  fFNLmomPimB.SetXYZT(0.0, 0.0, 0.0, -999.);
  fFNLmomPipB.SetXYZT(0.0, 0.0, 0.0, -999.);
  fPimLmom2.SetXYZT(0.0, 0.0, 0.0, -999.);
  fPipLmom2.SetXYZT(0.0, 0.0, 0.0, -999.);
  fVtxBeamPim = DEFVECT;
  fVtxBeamPip = DEFVECT;
  fVtxPimBeam = DEFVECT;
  fVtxPipBeam = DEFVECT;
}
