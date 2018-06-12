#ifndef NpipiData_h
#define NpipiData_h 1

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "GlobalVariables.h"

class NpipiData : public TObject
{
 public:
  NpipiData();
  virtual ~NpipiData(){};

 private:
  TLorentzVector fBeamLmom;
  TLorentzVector fFNLmom;
  TLorentzVector fPimLmom;
  TLorentzVector fPipLmom;
  TVector3 fVtxBeam;
  TVector3 fVtxCDS;
  TVector3 fVtxPim;
  TVector3 fVtxPip;
  //add 2015/6/14
  TLorentzVector fFNLmomPim;
  TLorentzVector fFNLmomPip;
  TLorentzVector fFNLmomPimB;
  TLorentzVector fFNLmomPipB;
  TLorentzVector fPimLmom2; //Calucleted by LineToHelixVertex;
  TLorentzVector fPipLmom2; //Calucleted by LineToHelixVertex;
  TVector3 fVtxPimBeam;
  TVector3 fVtxPipBeam;
  TVector3 fVtxBeamPim;
  TVector3 fVtxBeamPip;

 public:
  void setBeamLmom(const TLorentzVector &lmom){ fBeamLmom=lmom; };
  void setFNLmom(const TLorentzVector &lmom){ fFNLmom=lmom; };
  void setPimLmom(const TLorentzVector &lmom){ fPimLmom=lmom; };
  void setPipLmom(const TLorentzVector &lmom){ fPipLmom=lmom; };
  void setVtxBeam(const TVector3 &vtx){ fVtxBeam=vtx; };
  void setVtxCDS(const TVector3 &vtx){ fVtxCDS=vtx; };
  void setVtxPim(const TVector3 &vtx){ fVtxPim=vtx; };
  void setVtxPip(const TVector3 &vtx){ fVtxPip=vtx; };
  //add 2015/6/14
  void setFNLmomPim(const TLorentzVector &lmom){ fFNLmomPim=lmom; };
  void setFNLmomPip(const TLorentzVector &lmom){ fFNLmomPip=lmom; };
  void setFNLmomPimB(const TLorentzVector &lmom){ fFNLmomPimB=lmom; };
  void setFNLmomPipB(const TLorentzVector &lmom){ fFNLmomPipB=lmom; };
  void setPimLmom2(const TLorentzVector &lmom){ fPimLmom2=lmom; };
  void setPipLmom2(const TLorentzVector &lmom){ fPipLmom2=lmom; };
  void setVtxPimBeam(const TVector3 &vtx){ fVtxPimBeam=vtx; };
  void setVtxPipBeam(const TVector3 &vtx){ fVtxPipBeam=vtx; };
  void setVtxBeamPim(const TVector3 &vtx){ fVtxBeamPim=vtx; };
  void setVtxBeamPip(const TVector3 &vtx){ fVtxBeamPip=vtx; };

  TLorentzVector beamLmom() const { return fBeamLmom; };
  TLorentzVector fnLmom() const { return fFNLmom; };
  TLorentzVector pimLmom() const { return fPimLmom; };
  TLorentzVector pipLmom() const { return fPipLmom; };
  TVector3 vtxBeam() const { return fVtxBeam; };
  TVector3 vtxCDS() const { return fVtxCDS; };
  TVector3 vtxPim() const { return fVtxPim; };
  TVector3 vtxPip() const { return fVtxPip; };
  //add 2015/6/14
  TLorentzVector fnLmomPim() const { return fFNLmomPim; };
  TLorentzVector fnLmomPip() const { return fFNLmomPip; };
  TLorentzVector fnLmomPimB() const { return fFNLmomPimB; };
  TLorentzVector fnLmomPipB() const { return fFNLmomPipB; };
  TLorentzVector pimLmom2() const { return fPimLmom2; };
  TLorentzVector pipLmom2() const { return fPipLmom2; };
  TVector3 vtxPimBeam() const { return fVtxPimBeam; };
  TVector3 vtxPipBeam() const { return fVtxPipBeam; };
  TVector3 vtxBeamPim() const { return fVtxBeamPim; };
  TVector3 vtxBeamPip() const { return fVtxBeamPip; };

  void clear();

  ClassDef(NpipiData, 2);
};
#endif
