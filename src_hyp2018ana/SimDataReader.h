#ifndef SimDataReader_h
#define SimDataReader_h 1

#include "GlobalVariables.h"
#include "KnuclRootData.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "CDSTrackingMan.h"
#include "MyParam.h"
#include "SimTools.h"

class SimDataReader
{
 public:
  SimDataReader();
  virtual ~SimDataReader(){};

 private:
  RunHeaderMC *runHeader;
  EventHeaderMC *evHeader;
  ReactionData *reacData;
  DetectorData *detData;
  MCData *mcData;

  bool fTrigger;
  TLorentzVector fL1405_lmom;
  int fBeamPDG;
  int fTarPDG;
  TLorentzVector fBeamLmom;
  TLorentzVector fTarLmom;

 public:
  RunHeaderMC *getRunHeader(){ return runHeader; };
  EventHeaderMC *getEvenHeader(){ return evHeader; };
  ReactionData *getReactionData(){ return reacData; };
  DetectorData *getDetectorData(){ return detData; };
  MCData *getMCData(){ return mcData; };

  void setTree2(TTree *tree);
  void setTree(TTree *tree);

  void get();
  void clear();

  void checkStatus();
  void printSize();
  void printReaction(const int &max_gen=-1);
  void initHist(TFile *f);
  bool fillNpipi(TFile *f, CDSHitMan *cdsMan, const TLorentzVector &beam_lmom, CDSTrack *pim, const TLorentzVector &pim_lmom, 
		 CDSTrack *pip, const TLorentzVector &pip_lmom, const int &NCseg, const TLorentzVector &n_lmom, const double &mm);
  void fill(TFile *f);
  bool isNC();
  bool isCDH2();

  std::string reacName();
  std::string parName(const int &pdg);
  Track* getMCTrack(const int &trackID);
  DetectorHit *getHit(const int &cid, const int &seg);
  std::vector<DetectorHit*> getHits(const int &cid, const int &seg);
  DetectorHit *getHit(const int &cid, const int &lay, const int &wire);
  Track* trace(const int &cid, const int &seg, const bool &print=false); // return Original Track
  Track* trace(CDSTrack *cds_track, CDSHitMan* cdsMan, const bool &print=false); // return Original Track
  Track* trace(DetectorHit *hit, const bool &print=false); // return Original Track
  Track* trace(Track *track, const bool &print=false); // return Original Track
  Track* traceAll(Track *track, const bool &print=false); // return Original Track

  Track* traceNC(const int &seg, const bool &print=false); // return Original Track

  int NChitNeutronStatus(const int &NCseg);
  bool isNeutronInElasticNC(const int &seg);
  //  bool isInitialNC(const int &cid, const int &seg);
  //  bool isSigmaDecayNC(const int &cid, const int &seg);
  //  bool isPiDecayNC(const int &cid, const int &seg);
  //  bool isLambdaDecayNC(const int &cid, const int &seg);
  //  bool isL1520DecayNC(const int &cid, const int &seg);

  bool isNeutronInElastic(const int &cid, const int &seg);
  bool isNeutronInElastic(DetectorHit *hit);
  bool isFromGamma(const int &cid, const int &seg);
  bool isFromGamma(DetectorHit *hit);
  bool isFromE(const int &cid, const int &seg);
  bool isFromE(DetectorHit *hit);

  bool isInitial(const int &cid, const int &seg);
  bool isInitial(DetectorHit *hit);
  bool isSigmaDecay(const int &cid, const int &seg);
  bool isSigmaDecay(DetectorHit *hit);
  bool isPiDecay(const int &cid, const int &seg);
  bool isPiDecay(DetectorHit *hit);
  bool isLambdaDecay(const int &cid, const int &seg);
  bool isLambdaDecay(DetectorHit *hit);
  bool isL1520Decay(const int &cid, const int &seg);
  bool isL1520Decay(DetectorHit *hit);
  bool isHyperon(const int &pdg);

  void fillMM_N(TFile *f);
  void fillMM_N_woK0(TFile *f);
  void fillMM_N_woSf(TFile *f);
  void fillMM_N_wo(TFile *f);
  void fillpipi(TFile *f);
  void fill_woK0(TFile *f);
  void fill_woSf(TFile *f);
  void fillCharged(TFile *f);
  void fillMM_N_woRcS(TFile *f);
  void fillSp(TFile *f);
  void fillSm(TFile *f);
  void fillKNpipiMM(TFile *f, const double &mm);
  void fillKNpipiMM_woAll(TFile *f, const double &mm);
  void fillKNpipiMM_wK0(TFile *f, const double &mm);
  void fillKNpipiMM_wSm(TFile *f, const double &mm);
  void fillKNpipiMM_wSp(TFile *f, const double &mm);
  void fillKN_MM(TFile *f, const double &mm);
  void fillKN_MM_wN(TFile *f, const double &mm);

  void fillLpim(TFile *f);
  void fillLpim_p(TFile *f);

  void fillppim_MC(TFile *f);
  void printReac(const int &parentID);
  void printTrackTrace(const int &trackID);
  std::string getProcessName(const int &processID);

  void initHist_knEl();
  void fillHist_knEl();
  void fillHist_knEl_data(const TLorentzVector &beam_lmom, const TLorentzVector &n_lmom, const TLorentzVector &pim_lmom, const TLorentzVector &pip_lmom,
			  const TVector3 &pim_vtx, const TVector3 &pip_vtx, 
			  const TVector3 &pim_b_vtx, const TVector3 &pip_b_vtx, const TVector3 &b_pim_vtx, const TVector3  &b_pip_vtx,
			  bool mm_n_flag, bool im_K0_flag, bool im_Sm_flag, bool im_Sp_flag);

  //  int checkCDS(CDSTrack *track, CDSHitMan *cdsMan, bool print=false); // return <0 : read error, 0 : no problem 
  int checkCDS(CDSTrack *track, CDSTrackingMan *cdstrackingMan, CDSHitMan *cdsMan, const double &mass2, const TVector3 &vtxCDS, const TVector3 &vtxBeam, bool print=false); // return <0 : read error, 0 : no problem 
  bool isSameHit(HodoscopeLikeHit *hit1, HodoscopeLikeHit *hit2);
  void printProcess(const int &trackID);
  bool printCDCTrack(CDSHitMan *cdsMan, CDSTrack *track);

  //  ClassDef(SimDataReader, 1);
};
#endif
