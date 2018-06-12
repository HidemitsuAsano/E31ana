#ifndef CDS2Data_h
#define CDS2Data_h 1

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "LocalTrack.h"

#include "GlobalVariables.h"
#include "TrackTools.h"

class CDS2Data : public TObject
{
 public:
  CDS2Data();
  CDS2Data(const int &id1, const int &id2);
  virtual ~CDS2Data(){};

 private:
  int fID1;
  int fID2;
  int fPID1;
  int fPID2;
  int fBeamPID;
  double fDCA;

  TVector3 fBeamMom;
  TVector3 fCDSMom1;
  TVector3 fCDSMom2;

  TVector3 fVtxBeam;
  TVector3 fVtxCDS;
  TVector3 fVtxCDS1;
  TVector3 fVtxCDS2;

 public:
  int id1() const { return fID1; };
  int id2() const { return fID2; };
  int pid1() const { return fPID1; };
  int pid2() const { return fPID2; };

  TLorentzVector lmom1() const{
    TLorentzVector lmom;
    lmom.SetVectM(fCDSMom1, cdsMass[fPID1]);
    return lmom;
  }
  TLorentzVector lmom2() const{
    TLorentzVector lmom;
    lmom.SetVectM(fCDSMom2, cdsMass[fPID2]);
    return lmom;
  }
  TLorentzVector beam_lmom() const {
    TLorentzVector lmom;
    lmom.SetVectM(fBeamMom, parMass[fBeamPID]);
    return lmom;
  }
  TLorentzVector lmom() const { return (lmom1()+lmom2()); };
  TLorentzVector lmom(const int &pid){
    if( pid==fPID1 ) return lmom1();
    else if( pid==fPID2 ) return lmom2();
    else{
      std::cout<<"  !!! CDS2Data::lmom("<<pid<<") not found !!!"<<std::endl;
      return TLorentzVector();
    }
  }

  double im() const { return (lmom1()+lmom2()).M(); };
  double mm(TLorentzVector &tgt) const { return (beam_lmom()+tgt-lmom1()-lmom2()).M(); };

  TVector3 vtxBeam() const { return fVtxBeam; };
  TVector3 vtxCDS()  const { return fVtxCDS; };
  TVector3 vtxCDS_mean()  const { return 0.5*(fVtxCDS1+fVtxCDS2); };
  TVector3 vtxCDS1() const { return fVtxCDS1; };
  TVector3 vtxCDS2() const { return fVtxCDS2; };
  TVector3 vtxCDS(const int &pid){
    if( pid==fPID1 ) return vtxCDS1();
    else if( pid==fPID2 ) return vtxCDS2();
    else{
      std::cout<<"  !!! CDS2Data::vtxCDS("<<pid<<") not found !!!"<<std::endl;
      return TVector3();
    }
  }

  // This method return status  0: succeed  1: false Calc2HelixVertex  2: false cds1->GetMomentum  3: cds2->GetMomentum
  int set(const int &beam_pid, const double &mom, LocalTrack *beam, CDSTrack *cds1, CDSTrack *cds2);

  ClassDef(CDS2Data, 1);
};
#endif
