#ifndef DISPLAY_h
#define DISPLAY_h 1

#include <vector>
#include <iostream>
#include <fstream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include <TSystem.h>
#include <TVirtualPad.h>
#include <TH2F.h>
#include <TPolyLine.h>
#include <TArc.h>
#include <TLine.h>
#include <TMarker.h>


#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

class Display : public TObject
{
 public:
  Display();
  virtual ~Display() {};

 private:
  Display( Display & );

  TH2F *frameCDSXY;
  TH2F *frameCDSYZ;
  TH2F *frameCDSXZ;

  TH2F *frameBLXZ;
  TH2F *frameBLCXZ;
  TH2F *frameBLCYZ;

  TH2F *frameBLDCXZ;
  TH2F *frameBLDCYZ;

 public:
  bool Wait();
  bool DrawSegmentsXY( TVirtualPad *pad, ConfMan *conf, int cid );
  bool DrawSegmentsYZ( TVirtualPad *pad, ConfMan *conf, int cid );
  bool DrawSegmentsXZ( TVirtualPad *pad, ConfMan *conf, int cid );

  // local display for CDS
  void SetCDSFrameXY( double xmin=-80, double xmax=80, double ymin=-80, double ymax=80 );
  bool DrawCDSFrameXY( TVirtualPad *pad );
  bool DrawCDCLayersXY( TVirtualPad *pad, ConfMan *conf );
  bool DrawCDSHitXY( TVirtualPad *pad, ConfMan *conf, CDSHitMan *cds, int cid );
  bool DrawCDCTrackHitXY( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan );
  bool DrawCDCTrackXY( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan );

  void SetCDSFrameYZ( double ymin=-80, double ymax=80, double zmin=-80, double zmax=80 );
  bool DrawCDSFrameYZ( TVirtualPad *pad );
  bool DrawCDCLayersYZ( TVirtualPad *pad, ConfMan *conf );
  bool DrawCDSHitYZ( TVirtualPad *pad, ConfMan *conf, CDSHitMan *cds, int cid );
  bool DrawCDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan );
  bool DrawCDCTrackYZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan );

  void SetCDSFrameXZ( double xmin=-80, double xmax=80, double zmin=-80, double zmax=80 );
  bool DrawCDSFrameXZ( TVirtualPad *pad );
  bool DrawCDCLayersXZ( TVirtualPad *pad, ConfMan *conf );
  bool DrawCDSHitXZ( TVirtualPad *pad, ConfMan *conf, CDSHitMan *cds, int cid );
  bool DrawCDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan );
  bool DrawCDCTrackXZ( TVirtualPad *pad, ConfMan *conf, CDSTrackingMan *trackMan );

  // global display for BeamLine
  void SetBLFrameXZ( double xmin=-80, double xmax=80, double zmin=-80, double zmax=80 );
  bool DrawBLFrameXZ( TVirtualPad *pad );
  bool DrawBLHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid );
  bool DrawBLHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid );


  // display for Beamline DC
  void SetBLDCFrameXZ( double xmin=-80, double xmax=80, double zmin=-80, double zmax=80, char* title="" );
  void SetBLDCFrameYZ( double xmin=-80, double xmax=80, double zmin=-80, double zmax=80, char* title="" );

  bool DrawBLDCFrameXZ( TVirtualPad *pad );
  bool DrawBLDCFrameYZ( TVirtualPad *pad );
  bool DrawBLDCLayersXZ( TVirtualPad *pad, ConfMan *conf,int cid );
  bool DrawBLDCLayersYZ( TVirtualPad *pad, ConfMan *conf,int cid );
  bool DrawBLDCHit( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl, int cid,int xy );

  bool DrawBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);
  bool DrawBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);
  bool DrawTrackBLDCXZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);
  bool DrawTrackBLDCYZ( TVirtualPad *pad, ConfMan *conf, BeamLineHitMan *bl,
		   BeamLineTrackMan *track,int cid);

  bool DrawBLDCTrackHitXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid);
  bool DrawBLDCTrackHitYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid);

  bool DrawBLDCTrackXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid, const int &col=4);
  bool DrawBLDCTrackYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid, const int &col=4);

  bool DrawBLDCLTrackXZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid, const int &col=4);
  bool DrawBLDCLTrackYZ( TVirtualPad *pad, ConfMan *conf, BeamLineTrackMan *bl, const int &cid, const int &col=4);

  bool DrawBLC2TrackfromBLC1XZ( TVirtualPad *pad, ConfMan *conf, LinearTrack *track,const double &mom,const int &col);
  bool DrawBLC2TrackfromBLC1YZ( TVirtualPad *pad, ConfMan *conf, LinearTrack *track,const double &mom,const int &col);

  ClassDef( Display, 1 );
};
#endif
