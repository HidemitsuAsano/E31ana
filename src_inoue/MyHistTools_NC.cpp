#include "MyHistTools.h"
#include "MyParam.h"

void MyHistTools::initNC()
{
  new TH1F("NC_overbeta", "NC 1/#beta", 5000, 0.0, 4.0);
  new TH1F("KN_MM", "KN_MM", 5000, 0.0, 5.0);
}

void MyHistTools::fillNC(AnaInfo *anaInfo, BeamLineHitMan *blMan)
{
  TH1F *h1;
  TH2F *h2;
  if( !MyTools::isFiducial(anaInfo) ) return;

  TLorentzVector beam_lmom=anaInfo->beam()->lmom(anaInfo->beam(0)->vertex());
  TLorentzVector tgt_lmom;
  DetectorList *dlist = DetectorList::GetInstance();
  if(dlist->GetMaterial(CID_Target)=="LHydrogen" ) tgt_lmom.SetVectM(TVector3(),pMass);
  if(dlist->GetMaterial(CID_Target)=="LDeuterium") tgt_lmom.SetVectM(TVector3(),dMass);
  if(dlist->GetMaterial(CID_Target)=="LHelium-3" ) tgt_lmom.SetVectM(TVector3(),ThreeHeMass);

  if( anaInfo->nFNeutral()==1 ){
    h1 = (TH1F*)gFile-> Get("NC_overbeta"), h1-> Fill(1./anaInfo->forwardNeutral(0)->beta());
    if( anaInfo->forwardNeutral(0)->pid()==F_Neutron ){
      TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
      TLorentzVector mm_lmom=beam_lmom+tgt_lmom-n_lmom;
      h1 = (TH1F*)gFile-> Get("KN_MM"), h1-> Fill(mm_lmom.M());
    }
  }
}
