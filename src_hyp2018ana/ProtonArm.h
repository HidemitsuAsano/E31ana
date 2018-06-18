#ifndef PROTONARM_HH
#define PROTONARM_HH 1

#include "ConfMan.h"
#include "ChargeParticle.h"
#include "time.h"

#include "TPolyLine3D.h"

class ChargeParticle;

namespace ProtonArm{
  void Initialize(ConfMan *conf, double val=0.974854);  
  bool isUSWK_Field(const TVector3 &pos);
  TVector3 getField(const TVector3 &pos);

  //  bool nextStep_byRK4(const double &mass, const TVector3 &inpos, const TVector3 &inmom, TVector3 &outpos, TVector3 &outmom);

  ChargeParticle shootChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom, ConfMan *conf, std::vector<std::vector<double> > &trajectory);
  ChargeParticle shootChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom, ConfMan *conf, HodoscopeLikeHit *hit, std::vector<std::vector<double> > &trajectory);
  ChargeParticle shootChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom, ConfMan *conf);
  bool fit(const double &mass, TVector3 &vtx, TVector3 &mom, TVector3 &FDC1pos, HodoscopeLikeHit *fc_hit);
  bool fine_fit(const double &mass, TVector3 &vtx, TVector3 &mom, TVector3 &FDC1pos, HodoscopeLikeHit *fc_hit);

  bool fitBT(const double &mass, const double &mom, const TVector3 &BPCpos, const TVector3 &FDC1pos, ConfMan *conf);
  bool checkPreAna();

  TVector3 crossPointPCCVC(const TVector3 &pos, const TVector3 &dir);
  double disPCCVC(const TVector3 &pos, const TVector3 &dir);

  double momByAng(const double &angle);
}
#endif
