#ifndef KinFitMan_h
#define KinFitMan_h 1

#include "GlobalVariables.h"
#include "MyParam.h"
#include "SimDataReader.h"
#include "TMinuit.h"

class KinFitMan : public TObject
{
 public:
  KinFitMan();
  virtual ~KinFitMan(){ if( minuit ) delete minuit; };

 private:
  TMinuit *minuit;

  int fStatus;

  TLorentzVector fBeamLmom;
  TLorentzVector fTarLmom;
  int fPID1;
  TLorentzVector fLmom1;
  int fPID2;
  TLorentzVector fLmom2;

  double fMinMom;
  TLorentzVector fSpecLmom;
  TLorentzVector fReacLmom;

 public:
  double pipiMMSA_n_spec(const TLorentzVector &beam_lmom, const TLorentzVector &tar_lmom, const TLorentzVector &pip_lmom, const TLorentzVector &pim_lmom, TLorentzVector &n_lmom, TFile *f, SimDataReader *simReader=0);
  int getStatus() const { return fStatus; };
  double getMinMom() const { return fMinMom;};
  TLorentzVector getSpecLmom() const { return fSpecLmom; };
  TLorentzVector getReacLmom() const { return fReacLmom; };

  ClassDef(KinFitMan, 1);
};
#endif
