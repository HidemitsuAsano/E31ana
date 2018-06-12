#ifndef ANAINFO_HH
#define ANAINFO_HH 1

#include "CDSTrackingMan.h"
#include "BeamLineTrackMan.h"
#include "BeamInfo.h"
#include "CDSInfo.h"
#include "CDS2Info.h"
#include "ForwardChargeInfo.h"
#include "ForwardNeutralInfo.h"

class AnaInfo : public TObject
{
 public:
  AnaInfo();
  virtual ~AnaInfo(){};

 private:
  std::vector<BeamInfo> fBeamInfoContainer;
  std::vector<CDSInfo> fCDSInfoContainer;
  std::vector<CDS2Info> fCDS2InfoContainer;
  std::vector<ForwardNeutralInfo> fFNInfoContainer;
  std::vector<ForwardChargeInfo> fFCInfoContainer;
  
 public:
  void AddBeam(BeamInfo &info){ fBeamInfoContainer.push_back(info); }
  void AddCDS(CDSInfo &info){ fCDSInfoContainer.push_back(info); }
  void AddCDS2(CDS2Info &info){ fCDS2InfoContainer.push_back(info); }
  void AddNeutral(ForwardNeutralInfo &info){ fFNInfoContainer.push_back(info); }
  void AddCharge(ForwardChargeInfo &info){ fFCInfoContainer.push_back(info); }

  int nBeam() const { return (int)fBeamInfoContainer.size(); }
  int nCDS() const { return (int)fCDSInfoContainer.size(); }
  int nCDS2() const { return (int)fCDS2InfoContainer.size(); }
  int nFCharge() const { return (int)fFCInfoContainer.size(); }
  int nFNeutral() const { return (int)fFNInfoContainer.size(); }

  BeamInfo* beam(const int &i=0){ return (0<=i && i<(int)fBeamInfoContainer.size() ) ? &fBeamInfoContainer[i] : 0; }
  CDSInfo* CDS(const int &i){ return (0<=i && i<(int)fCDSInfoContainer.size() ) ? &fCDSInfoContainer[i] : 0; }
  CDS2Info* CDS2(const int &i){ return (0<=i && i<(int)fCDS2InfoContainer.size() ) ? &fCDS2InfoContainer[i] : 0; }
  ForwardNeutralInfo* forwardNeutral(const int &i=0){ return (0<=i && i<(int)fFNInfoContainer.size() ) ? &fFNInfoContainer[i] : 0; }
  ForwardChargeInfo* forwardCharge(const int &i=0){ return (0<=i && i<(int)fFCInfoContainer.size() ) ? &fFCInfoContainer[i] : 0; }

  int nCDS(const int &pid);
  CDSInfo* CDS(const int &pid, const int &i);
  CDSInfo* CDSbyID(const int &id);

  int nCDS2(int pid1, int pid2);
  CDS2Info* CDS2(int pid1, int pid2, const int &i);

  CDSInfo *minDCA();

  void Clear();
  void ClearBeam(){ fBeamInfoContainer.clear(); }
  void ClearFN(){ fFNInfoContainer.clear(); }
  void ClearFC(){ fFCInfoContainer.clear(); }
  void ClearCDS(){ fCDSInfoContainer.clear(); fCDS2InfoContainer.clear(); }

  void dump();

  ClassDef(AnaInfo, 1);
};
#endif
