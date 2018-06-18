#ifndef BEAMINFO_HH
#define BEAMINFO_HH 1

#include "GlobalVariables.h"
#include "BeamSpectrometer.h"
#include "TLorentzVector.h"
#include "CDSTrackingMan.h"

class BLDCTrackInfo : public TObject
{
 public:
  BLDCTrackInfo();
  BLDCTrackInfo(int id, LocalTrack *track);
  virtual ~BLDCTrackInfo(){};

 private:
  int fID;
  double fTime;
  double fChi2;
  TVector3 fPos;
  TVector3 fDir;

 public:
  int id() const { return fID; };
  double time() const { return fTime; };
  double chi2() const { return fChi2; };
  TVector3 pos() const { return fPos; };
  TVector3 dir() const { return fDir; };
  TVector3 GetPosatZ(const double z) const { return fPos+fDir*((z-fPos.Z())/fDir.Z()); };

  ClassDef(BLDCTrackInfo, 1);
};

class BeamInfo : public TObject
{
 public:
  BeamInfo();
  virtual ~BeamInfo(){};

 private:
  bool fFlag;
  int fPID; 
  int fBHDseg;
  double fBHDtime;
  int fT0seg;
  double fT0time;
  double fD5mom;
  double fD5chisquare;
  int fBLC1id;
  int fBLC2id;
  int fBPCid;
  std::vector<BLDCTrackInfo> fBLC1info;
  std::vector<BLDCTrackInfo> fBLC2info;
  std::vector<BLDCTrackInfo> fBPCinfo;
  std::vector<BLDCTrackInfo> fFDC1info;
  TVector3 fT0pos;
  TVector3 fVertex;
  double fVertexMom;

  std::map<int, double> fBHDtimes;

  TVector3 fFDC1pos;
  TVector3 fFDC1dir;

 public:
  bool flag() const { return fFlag; };
  int pid() const { return fPID; };
  double mass() const { return (0<=fPID && fPID<(int)(sizeof(parMass)/sizeof(parMass[0]))) ? parMass[fPID] : 0.0; }
  int BHDseg() const { return fBHDseg; };
  int T0seg() const { return fT0seg; };
  double BHDtime() const { return fBHDtime; };
  double T0time() const { return fT0time; };
  double D5mom() const { return fD5mom; };

  int nBLC1() const { return fBLC1info.size(); };
  int nBLC2() const { return fBLC2info.size(); };
  int nBPC() const { return fBPCinfo.size(); };
  int nFDC1() const { return fFDC1info.size(); };
  BLDCTrackInfo BLC1(const int &i){ return fBLC1info[i]; };
  BLDCTrackInfo BLC2(const int &i){ return fBLC2info[i]; };
  BLDCTrackInfo BPC(const int &i){ return fBPCinfo[i]; };
  BLDCTrackInfo FDC1(const int &i){ return fFDC1info[i]; };

  double BLC1time() const { return fBLC1id<0 ? DEFAULTD : fBLC1info[fBLC1id].time(); };
  double BLC2time() const { return fBLC2id<0 ? DEFAULTD : fBLC2info[fBLC2id].time(); };
  double BPCtime() const { return fBPCid<0 ? DEFAULTD : fBPCinfo[fBPCid].time(); };

  double BLC1chi2() const { return fBLC1id<0 ? DEFAULTD : fBLC1info[fBLC1id].chi2(); };
  double BLC2chi2() const { return fBLC2id<0 ? DEFAULTD : fBLC2info[fBLC2id].chi2(); };
  double BPCchi2() const { return fBPCid<0 ? DEFAULTD : fBPCinfo[fBPCid].chi2(); };

  TVector3 BLC1dir() const { return fBLC1id<0 ? DEFVECT : fBLC1info[fBLC1id].dir(); };
  TVector3 BLC2dir() const { return fBLC2id<0 ? DEFVECT : fBLC2info[fBLC2id].dir(); };
  TVector3 BPCdir() const { return fBPCid<0 ? DEFVECT : fBPCinfo[fBPCid].dir(); };

  TVector3 BLC1pos() const { return fBLC1id<0 ? DEFVECT : fBLC1info[fBLC1id].pos(); };
  TVector3 BLC2pos() const { return fBLC2id<0 ? DEFVECT : fBLC2info[fBLC2id].pos(); };
  TVector3 BPCpos() const { return fBPCid<0 ? DEFVECT : fBPCinfo[fBPCid].pos(); };

  TVector3 T0pos() const { return fT0pos; };
  double D5chi2() const { return fD5chisquare; };
  TVector3 vertex() const { return fVertex; };
  double vertex_mom() const { return fVertexMom; };

  int nBHD() const { return fBHDtimes.size(); }
  bool getBHD(const int &i, int &seg, double &time);
  double BHDtime(const int &seg) const { return fBHDtimes.find(seg)!=fBHDtimes.end() ? fBHDtimes.at(seg) : 0.0; }

  double tof() const { return fT0time-fBHDtime; };
  TLorentzVector lmom() const;
  TLorentzVector lmom(const TVector3 &pos) const;
  HodoscopeLikeHit* T0(BeamLineHitMan *blMan) const;
  HodoscopeLikeHit* BHD(BeamLineHitMan *blMan) const;
  LocalTrack *BLC1(BeamLineTrackMan *bltrackMan) const { return bltrackMan->trackBLC1(fBLC2id); };
  LocalTrack *BLC2(BeamLineTrackMan *bltrackMan) const { return bltrackMan->trackBLC2(fBLC2id); };
  LocalTrack *BPC(BeamLineTrackMan *bltrackMan) const { return bltrackMan->trackBPC(fBPCid); };

  void SetBLMan(BeamLineHitMan *man);
  void SetBLDC(BeamLineTrackMan *bltrackMan);
  void SetFlag(bool flag){ fFlag=flag; };
  void SetPID(const int &id){ fPID=id; };
  void SetBLC1id(const int &id){ fBLC1id=id; };
  void SetBLC2id(const int &id){ fBLC2id=id; };
  void SetBPCid(const int &id, ConfMan *conf=0);
  void SetD5(BeamSpectrometer *beam){ fD5mom=beam->mom(); fD5chisquare=beam->chisquare(); };
  void SetT0(HodoscopeLikeHit *hit){ fT0seg=hit->seg(); fT0time=hit->ctmean(); };
  void SetBHD(HodoscopeLikeHit *hit){ fBHDseg=hit->seg(); fBHDtime=hit->ctmean(); };

  void SetVertex(const TVector3 &pos);
  void SetVtxMom(const double &mom){ fVertexMom=mom; };

  void SetD5mom(const double &mom){ fD5mom=mom; fD5chisquare=0.0; }; // for Monte Carlo Data
  void SetT0time(const double &time){ fT0time=time; }; // for Monte Carlo Data

  void dump() const;
  double calc_tof() const;

  ClassDef(BeamInfo, 1);
};


#endif
