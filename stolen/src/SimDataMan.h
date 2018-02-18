#ifndef SimDataMan_h
#define SimDataMan_h 1

#include <vector>
#include <string>

#include "TObject.h"
#include "TTree.h"

#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "HodoscopeLikeHit.h"
#include "ChamberLikeHit.h"
#include "CDCHit.h"

#include "ComCrossSectionTable.hh"
#include "KnuclRootData.h"

class SimDataMan
{
 public:
  SimDataMan( TFile *file, ConfMan *conf, bool tracking_flag = false );
  virtual ~SimDataMan(){};

 private:
  SimDataMan();

 private:
  TFile *inFile;
  TFile *outFile;
  TTree *tree2;
  TTree *tree;

  ConfMan *confMan;

  RunHeaderMC   *runHeader;
  EventHeaderMC *eventHeader;
  DetectorData  *detectorData;
  MCData        *mcData;
  ReactionData  *reactionData;

  int NumOfEvent;

  std::vector<int> CDHHitPDG;
  int nIH;
  int nBPD;
  int nT0;
  int nDEF;
  int nBVC;
  int nNC;
  int nPC;
  int nCVC;

  int nCDC[NumOfCDCLayers];
  int nBLC2a[8];
  int nBLC2b[8];
  int nBPC[8];
  int nFDC1[6];

  double T0Time;
  double BeamMom;

 public:
  RunHeaderMC* GetRunHeader(){ return runHeader; };
  EventHeaderMC* GetEventHeader(){ return eventHeader; };
  DetectorData* GetDetectorData(){ return detectorData; };
  MCData* GetMCData(){ return mcData; };
  ReactionData* GetReactionData(){ return reactionData; };

  int nEvent() const { return NumOfEvent; };
  void GetRunInfo(){ tree2->GetEntry(0); };
  bool Get(const int &evn);
  bool Get(const int &evn, CDSHitMan *cdsMan, BeamLineHitMan *blMan);

  double T0time() const { return T0Time; };
  double beam_mom() const { return BeamMom; };

  void PrintHeader();
  void PrintEvent();
  void PrintTrack(Track *track);
  void PrintDetectorHit(DetectorHit *hit);
  void PrintHit();

  std::string PDGToName(const int &pdg_id) const;
  std::string CIDToName(const int &cid) const;

  void Clear();
  bool Check();

};

#endif

