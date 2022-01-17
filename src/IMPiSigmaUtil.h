#ifndef IMPISIGMAUTIL
#define IMPISIGMAUTIL 1

#include <iostream>
#include <cmath>
#include "ConfMan.h"
#include "HitMan.h"
#include "CDSHitMan.h"
#include "LocalTrack.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "CDSTrackingMan.h"
#include "IMPiSigmaAnaPar.h"
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "KnuclRootData.h"

namespace Util
{ 

  int GetCDHMul(CDSHitMan *cdsman,const int ntrack=0,const bool CDH3trgfired=true, const bool MCFlag=false);
  bool IsForwardCharge(BeamLineHitMan *blman);
  int GetCDHNeighboringNHits(const std::vector <int> &seg, const std::vector <int> &allhit, 
  const std::vector <int> &pippimseg,CDSHitMan *cdsman, bool MCFlag=false);
  
  int GetCDHTwoSegAwayNHits(const std::vector <int> &seg, const std::vector <int> &allhit, bool MCFlag=false);
  
  
  int GetNHitsCDCInner3Lay(CDSHitMan *cdsman);
  int GetNHitsCDCOuter(TVector3 PosCDH, CDSHitMan *cdsman, const double rangedeg=15.0);
  int GetNHitsCDCOuterNoAss(TVector3 PosCDH,CDSHitMan *cdsman, CDSTrackingMan *trackman, const double rangedeg=15.0);
  void AnaPipPimCDCCDH(const TVector3 PosCDH,const std::vector <int> &seg, const int pip_ID, const int pim_ID, CDSHitMan *cdsman,CDSTrackingMan *trackman);
  
  
  double AnaBeamSpec(ConfMan *confman,BeamLineTrackMan *bltrackman,const int blc1id, const int blc2id);
  int CDSChargedAna(const bool docdcretiming,
                    LocalTrack *bpctrack,
                    CDSHitMan *cdsman,
                    CDSTrackingMan *trackman,
                    ConfMan *confman,
                    BeamLineHitMan *blman,
                    const TLorentzVector beam,
                    const double ctmt0, 
                    std::vector <int> &cdhseg,
                    std::vector <int> &pimid,
                    std::vector <int> &pipid,
                    std::vector <int> &kmid,
                    std::vector <int> &protonid,
                    TVector3 &pim_projected,
                    TVector3 &pip_projected,
                    const bool MCFlag=false,
                    const int H2data=0
                    );
  double AnalyzeT0(BeamLineHitMan *blman, ConfMan *confman,int &t0seg);
  int BeamPID(EventHeader *header, const double ctmt0, BeamLineHitMan *blman, const int H2data=0); 
  int EveSelectBeamline(BeamLineTrackMan *bltrackman, 
                        CDSTrackingMan *trackman,
                        ConfMan *confman,
                        int &blc1id, 
                        int &blc2id, 
                        int &bpcid);

  void AnaCDHHitPos(const double meas_tof, const double beta_calc, 
                 LocalTrack *bpc,
                 const TLorentzVector LVec_beambf,
                 const double ctmt0,
                 CDSTrack *track,
                 CDSHitMan *cdsman,
                 ConfMan *confman,
                 BeamLineHitMan *blman,
                 const double correctedtof,
                 const int pid,
                 const bool MCFlag=false
                 );
   //MC only
   void AnaReactionData(ReactionData *reactionData);
   void AnaMcData(MCData *mcdata, 
                  DetectorData  *detdata,
                  CDSHitMan *cdsman,
                  ReactionData *reactionData,
                  double &ncanvtxr,
                  double &ncanvtxz,
                  int &ncangeneration
                  );
   void AnaMcData2(MCData *mcdata, 
                  DetectorData  *detdata,
                  const int CDHseg,
                  double &ncanvtxr,
                  double &ncanvtxz,
                  double &firstvtxr,
                  double &firstvtxz,
                  int &ncangeneration,
                  int &mcpattern,
                  int &nanc
                  );
   int ProcessNameToProcessID(const std::string &name);
   std::string ProcessIDToProcessName(const int &id);
   int CalcGeneration(MCData *mcdata,int parentid);
   double FillAncestryVertexR(MCData *mcdata,DetectorHit *dhit,double dE);
   double FillAncestryVertexZ(MCData *mcdata,DetectorHit *dhit,double dE);
   Track *FindTrackFromMcIndex(MCData *mcdata, int trackid);
   bool IsFromSigma(MCData *mcdata,DetectorHit *dhit);
   TLorentzVector *GetForwardNeutralLVec(BeamLineHitMan *blman,const TVector3 vtxpos,const double t0time,const double beamtof, const double thre);
   
   std::vector<std::vector<HodoscopeLikeHit*> > getNChits(BeamLineHitMan *blman);
   std::vector<HodoscopeLikeHit*>  getHodo(BeamLineHitMan *blman);

};

#endif
