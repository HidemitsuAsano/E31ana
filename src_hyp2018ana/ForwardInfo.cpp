#include "ForwardInfo.h"

ForwardInfo::ForwardInfo()
{
}

TLorentzVector ForwardInfo::lmom() const
{
  TLorentzVector lmom;
  lmom.SetVectM(fMomentum, pdgMass());
  return lmom;
}
