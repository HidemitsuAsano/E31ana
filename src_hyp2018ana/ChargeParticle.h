#ifndef CHARGEPARTICLE_HH
#define CHARGEPARTICLE_HH 1

#include "GeomTools.h"
#include "ELossTools.h"
#include "ProtonArm.h"
#include "TObject.h"
#include "ProtonArm.h"

#include "TCanvas.h"
#include "TLine.h"

class ChargeParticle : public TObject
{
 public:
  ChargeParticle(const double &mass, const TVector3 &pos, const TVector3 &mom);
  ~ChargeParticle(){};

 private:
  bool fFlag;
  double fMass;
  const double fInitMom;
  TVector3 fFDC1pos;
  double fDiffCounter;
  double fDeltaTime;
  double fFlightLength;
  double fTime;
  TVector3 fPos;
  TVector3 fMom;

 public:
  void SetFlag(const bool &flag){ fFlag=flag; }
  void SetFDC1pos(const TVector3 &pos){ fFDC1pos=pos; }
  void SetDiffCounter(const double &val){ fDiffCounter=val; }

  bool flag() const { return fFlag; };
  double mass() const { return fMass; };
  double initMom() const { return fInitMom; };
  TVector3 FDC1pos() const { return fFDC1pos; }
  double diffCounter() const { return fDiffCounter; }
  double deltaTime() const { return fDeltaTime; };
  double fl() const { return fFlightLength; };
  double time() const { return fTime; };
  TVector3 pos() const { return fPos; };
  TVector3 mom() const { return fMom; };
  TVector3 v() const { return (100.*Const/sqrt(fMom.Mag2()+fMass*fMass))*fMom; };

  bool calcELossBeamTGeo(const double &z);
  bool calcELossBeamTGeo(const double &z, const double &cds_field);
  bool nextStepRK4(const double &dt);

  void dump() const;

 private:
  double tofByInitMom(const double &fl) const { return fl*sqrt(fInitMom*fInitMom+fMass*fMass)/(100.*Const*fInitMom); }

  ClassDef(ChargeParticle, 1);
};
#endif
