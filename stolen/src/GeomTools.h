// GeomTools.h

#ifndef GeomTools_h
#define GeomTools_h 1

#include <string>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "GlobalVariables.h"
#include "ConfMan.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"

namespace GeomTools
{
  void MakeGeometry(ConfMan* conf);
  void MakeGeometryHeates(ConfMan* conf);
  void MakeGeometryPSI(ConfMan* conf);
  void ConstructHadronHall(ConfMan* conf);
  void ConstructTarget(ConfMan* conf);
  void ConstructPSICarbonTarget(Int_t CID,ConfMan *conf,TGeoVolume *mother=0);
  void ConstructHodoscopes(ConfMan* conf); 
  void ConstructChambers(ConfMan* conf);
  void ConstructMan();
  void ConstructMaterial();
  //  bool Check(const TVector3 &pos1,const TVector3 &pos2, TString &mat, int oldid);
  bool CrossBoundary(const TVector3 &pos1,const TVector3 &pos2,const double &step, TString &mat, double &tmpstep, int &newid);
  bool GetParam(const int &cid,const int &seg,ConfMan* conf, double *param);
  inline int GetID(const TVector3 &pos);
  inline int GetIDMat(const TVector3 &pos,TString &mat);
  //  bool SetInit(const TVector3& pos1,const TVector3 &pos2,TString &matname, int &id);
  bool StepToNextVolume(const TVector3 &in, const TVector3 &out,double &length, TString &newmat);
  bool HelixStepToNextVolume(const double param[5],const TVector3 &pos1,TVector3 &pos2, double &length, TString &mat, int &id);
  //  bool IsSameVolumeHelix(const double param[5], const TVector3 &pos1, const TVector3 &pos2,double margin=0.);
  bool IsSameVolume(const TVector3 &pos1, const TVector3 &pos2,double margin=0.);
  bool IsSameVolumeHelix(const double param[5],const TVector3 &pos1, const TVector3 &pos2,double margin=0.);
  double CalcLengthinFiducial(const TVector3 &pos1, const TVector3 &dir2);
  TGeoVolume* ConstructShape(Int_t CID,ConfMan *conf,TGeoVolume* mother=0);
  TGeoVolume* MakeShape(double *param, const TString &name, TGeoMedium* medium);
  TGeoCombiTrans *MakeTrans( double *param);
};

inline int GeomTools::GetID(const TVector3 &pos)
{
  //  if(pos.Perp()>200) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return 0;
  int id=node->GetNumber();
  id/=1000;
  return id;
}

inline int GeomTools::GetIDMat(const TVector3 &pos,TString &mat)
{
  //  if(pos.Mag()>100) return 0;
  TGeoNode *node=gGeoManager->FindNode(pos.X(),pos.Y(),pos.Z());
  if(!node) return 0;
  int id=node->GetNumber();
  id/=1000;
  mat = node->GetMedium()->GetName();
  return id;
}


#endif
