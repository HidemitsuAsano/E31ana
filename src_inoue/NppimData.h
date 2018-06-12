#ifndef NppimData_h
#define NppimData_h 1

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "GlobalVariables.h"

class NppimData : public TObject
{
 public:
  NppimData();
  virtual ~NppimData(){};

 private:
  TLorentzVector fBeamLmom;
  TLorentzVector fFNLmom;
  TLorentzVector fPimLmom;
  TLorentzVector fPLmom;
  TVector3 fVtxBeam;
  TVector3 fVtxCDS;
  TVector3 fVtxPim;
  TVector3 fVtxP;

 public:
  void setBeamLmom(const TLorentzVector &lmom){ fBeamLmom=lmom; };
  void setFNLmom(const TLorentzVector &lmom){ fFNLmom=lmom; };
  void setPimLmom(const TLorentzVector &lmom){ fPimLmom=lmom; };
  void setPLmom(const TLorentzVector &lmom){ fPLmom=lmom; };
  void setVtxBeam(const TVector3 &vtx){ fVtxBeam=vtx; };
  void setVtxCDS(const TVector3 &vtx){ fVtxCDS=vtx; };
  void setVtxPim(const TVector3 &vtx){ fVtxPim=vtx; };
  void setVtxP(const TVector3 &vtx){ fVtxP=vtx; };

  TLorentzVector beamLmom() const { return fBeamLmom; };
  TLorentzVector fnLmom() const { return fFNLmom; };
  TLorentzVector pimLmom() const { return fPimLmom; };
  TLorentzVector pLmom() const { return fPLmom; };
  TVector3 vtxBeam() const { return fVtxBeam; };
  TVector3 vtxCDS() const { return fVtxCDS; };
  TVector3 vtxPim() const { return fVtxPim; };
  TVector3 vtxP() const { return fVtxP; };

  void clear();

  ClassDef(NppimData, 1);
};
#endif
