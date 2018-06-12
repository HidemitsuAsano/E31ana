/* BeamLineTrackMan.h */

#ifndef BeamLineTrackMan_h
#define BeamLineTrackMan_h 1

#include <iostream>
#include <vector>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "BeamLineHitMan.h"
#include "BLDCClusterMan.h"
#include "LocalTrack.h"
#include "ChamberLikeHit.h"
#include "ConfMan.h"
#include "BLDCFittingParamMan.h"

typedef std::vector<LocalTrack> BeamLineTrackContainer;
class BeamLineTrackMan : public TObject
{
 private:
  BeamLineTrackContainer BLC1TrackContainer;
  BeamLineTrackContainer BLC1aTrackContainer;
  BeamLineTrackContainer BLC1bTrackContainer;
  BeamLineTrackContainer BLC2TrackContainer;
  BeamLineTrackContainer BLC2aTrackContainer;
  BeamLineTrackContainer BLC2bTrackContainer;
  BeamLineTrackContainer BPCTrackContainer;
  BeamLineTrackContainer FDC1TrackContainer;

  BeamLineTrackContainer* BLDCTrackContainer(const int &cid);
  int STATUS[10];

 public:
  BeamLineTrackMan();
  BeamLineTrackMan( const BeamLineTrackMan &right );
  ~BeamLineTrackMan();

 public:
  int status(const int &cid);
  int SetStatus(const int &cid, const int &sta);

  int ntrackBLDC(const int &cid,const double &mintime=-9999.,const double &maxtime=9999);
  LocalTrack *trackBLDC(const int &cid, const unsigned int &i);
  void SetTrack( const int &id, LocalTrack track );
  void SetTrackBLC1(  LocalTrack track ){ BLC1TrackContainer.push_back(track); }
  void SetTrackBLC1a( LocalTrack track ){ BLC1aTrackContainer.push_back(track); }
  void SetTrackBLC1b( LocalTrack track ){ BLC1bTrackContainer.push_back(track); }

  void SetTrackBLC2(  LocalTrack track ){ BLC2TrackContainer.push_back(track); }
  void SetTrackBLC2a( LocalTrack track ){ BLC2aTrackContainer.push_back(track); }
  void SetTrackBLC2b( LocalTrack track ){ BLC2bTrackContainer.push_back(track); }

  void SetTrackBPC(   LocalTrack track ){ BPCTrackContainer.push_back(track); }
  void SetTrackFDC1(  LocalTrack track ){ FDC1TrackContainer.push_back(track); }

  int ntrackBLC1()  const { return BLC1TrackContainer.size(); }
  int ntrackBLC1a() const { return BLC1aTrackContainer.size(); }
  int ntrackBLC1b() const { return BLC1bTrackContainer.size(); }

  int ntrackBLC2()  const { return BLC2TrackContainer.size(); }
  int ntrackBLC2a() const { return BLC2aTrackContainer.size(); }
  int ntrackBLC2b() const { return BLC2bTrackContainer.size(); }

  int ntrackBPC()  const { return BPCTrackContainer.size(); }

  int ntrackFDC1()  const { return FDC1TrackContainer.size(); }

  LocalTrack *trackBLC1( const int &i ){ return &(BLC1TrackContainer[i]); }
  LocalTrack *trackBLC1a(const int &i ){ return &(BLC1aTrackContainer[i]); }
  LocalTrack *trackBLC1b(const int &i ){ return &(BLC1bTrackContainer[i]); }

  LocalTrack *trackBLC2( const int &i ){ return &(BLC2TrackContainer[i]); }
  LocalTrack *trackBLC2a(const int &i ){ return &(BLC2aTrackContainer[i]); }
  LocalTrack *trackBLC2b(const int &i ){ return &(BLC2bTrackContainer[i]); }

  LocalTrack *trackBPC( const int &i ){ return &(BPCTrackContainer[i]); }
  LocalTrack *trackFDC1( const int &i ){ return &(FDC1TrackContainer[i]); }

  void DeleteTrackBLC1(  const int &i );
  void DeleteTrackBLC1a( const int &i );
  void DeleteTrackBLC1b( const int &i );
  void DeleteTrackBLC2(  const int &i );
  void DeleteTrackBLC2a( const int &i );
  void DeleteTrackBLC2b( const int &i );
  void DeleteTrackBPC(   const int &i );
  //  void DeleteTrackFDC1(  const int &i );

  void Clear();

 public:
  bool DoTracking( BeamLineHitMan *blMan, ConfMan *conf, const bool &SEMILOCAL=false, const bool &TIMING=false );  
  int LocalTracking( BeamLineHitMan *blMan, ConfMan *conf, const int &id, TString option="" );

  bool ConvertLocalToLinear(ConfMan* conf,const int &cid);

 private:
  bool Clustering( BeamLineHitMan *blMan, BLDCClusterMan *clMan, 
		   ConfMan *conf, const int &id, const int &id2=0 );
  void ConvLocalToGlobal(const int &id);

  ClassDef(BeamLineTrackMan,1);
};

#endif
