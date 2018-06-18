#include "MyHistTools.h"
#include "MyParam.h"

using namespace std;

void MyHistTools::initLpim()
{
  new TH2F("CDS_ppim_ppim_IM", "CDS p pi^{-} IM vs p pi^{-} IM",  1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("Lpim_IM",  "#Lambda #pi^{-} IM",                      2000, 1.0, 3.0);
  new TH1F("KLpim_MM", "d(K^{-}, #Lambda #pi^{-})\"X\"",          2000, 0.0, 2.0);
  new TH2F("Lpim_IM_KLpim_MM", "#Lambda #pi^{-} IM vs d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 1.0, 3.0, 2000, 0.0, 2.0);

  new TH1F("KLpim_MM", "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_C", "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_C_wFDC1", "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_C_FC",    "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_C_fp",    "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);

  new TH2F("KLpim_MM_FC_mass2", "d(K^{-}, #Lambda #pi^{-})\"X\" vs forward mass2", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH1F("Lpim_IM_mm_p",          "#Lambda #pi^{-} IM", 2000, 1.0, 3.0);
  new TH1F("Lpim_IM_mm_p_CDH2", "#Lambda #pi^{-} IM", 2000, 1.0, 3.0);
  new TH1F("Lpim_IM_mm_p_C",    "#Lambda #pi^{-} IM", 2000, 1.0, 3.0);
  new TH1F("Lpim_IM_mm_p_C_wFDC1", "#Lambda #pi^{-} IM", 2000, 1.0, 3.0);
}

void MyHistTools::fillLpim(EventHeader *header, AnaInfo *anaInfo, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan)
{
  TH1F *h1;
  TH2F *h2;

  if( anaInfo->nBeam()!=1 ) return;
  if( !MyTools::isFiducial(anaInfo) ) return;
  if( cdstrackMan->nGoodTrack()!=3 ) return;
  if( anaInfo->nCDS(CDS_PiMinus)!=2 ) return;
  if( anaInfo->nCDS(CDS_Proton)!=1 ) return;
  BeamInfo *beam=anaInfo->beam(0);

  CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
  CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
  CDSInfo *p=anaInfo->CDS(CDS_Proton, 0);

  if( !pim0->flag() ) return;
  if( !pim1->flag() ) return;
  if( !p->flag() ) return;

  CDS2Info *ppim0=anaInfo->CDS2(CDS_Proton, CDS_PiMinus, 0);
  CDS2Info *ppim1=anaInfo->CDS2(CDS_Proton, CDS_PiMinus, 1);
  if( ppim1->trackID1()==pim0->trackID() || ppim1->trackID2()==pim0->trackID() ){
    swap(ppim1, ppim0);
  }
  if( !ppim0->flag() ) return;
  if( !ppim1->flag() ) return;

  double im0=ppim0->im();
  double im1=ppim1->im();
  if( 1.341<im0 && im0<1.342 ) ppim0->dump();
  if( 1.341<im1 && im1<1.342 ) ppim1->dump();


  if( ppim0->dca()>ppim1->dca() ){
    swap(ppim1, ppim0); swap(pim1, pim0);
  }
  h2 = (TH2F*)gFile-> Get("CDS_ppim_ppim_IM"), h2-> Fill(im0, im1);

  bool L_flag0=false, L_flag1=false;
  if( L_MIN<im0 && im0<L_MAX ) L_flag0=true;
  if( L_MIN<im1 && im1<L_MAX ) L_flag1=true;
  if( L_flag0 && L_flag1 ){
    double diff0 = fabs(im0-L_peak);
    double diff1 = fabs(im1-L_peak);
    if( diff0<diff1 ) L_flag1=false;
    else L_flag0=false;
  }
  if( !L_flag0 && !L_flag1 ) return;

  TLorentzVector l_lmom, pim_lmom;
  CDSInfo *pim=0;
  if( L_flag0 ){
    l_lmom=ppim0->lmom();
    pim=pim1;
    pim_lmom=pim1->lmom();
  }
  else{
    l_lmom=ppim1->lmom();
    pim=pim0;
    pim_lmom=pim0->lmom();
  }
  if( !MyTools::isFiducial(pim) ) return;

  TLorentzVector tgt_lmom; tgt_lmom.SetVectM(TVector3(0, 0, 0), dMass);
  double im=(l_lmom+pim_lmom).M();
  TLorentzVector beam_lmom=anaInfo->beam(0)->lmom(pim->vertexBeam());
  double mm=(beam_lmom+tgt_lmom-l_lmom-pim_lmom).M();

  h1 = (TH1F*)gFile-> Get("Lpim_IM"), h1-> Fill(im);
  h1 = (TH1F*)gFile-> Get("KLpim_MM"), h1-> Fill(mm);
  h2 = (TH2F*)gFile-> Get("Lpim_IM_KLpim_MM"), h2-> Fill(im, mm);
  if( header->IsTrig(Trig_Charged) ){
    h1 = (TH1F*)gFile-> Get("KLpim_MM_C"), h1-> Fill(mm);
    if( beam->nFDC1()==1 ){
      h1 = (TH1F*)gFile-> Get("KLpim_MM_C_wFDC1"), h1-> Fill(mm);
      if( anaInfo->nFCharge()==1 ){
	h1 = (TH1F*)gFile-> Get("KLpim_MM_C_FC"), h1-> Fill(mm);
	h2 = (TH2F*)gFile-> Get("KLpim_MM_FC_mass2"), h2-> Fill(mm, anaInfo->forwardCharge(0)->mass2());
	if( anaInfo->forwardCharge(0)->pid()==F_Proton ){
	  h1 = (TH1F*)gFile-> Get("KLpim_MM_C_fp"), h1-> Fill(mm);
	}
      }
    }
  }
  
  if( Lpim_P_MIN<mm && mm<Lpim_P_MAX ){
    h1 = (TH1F*)gFile->Get("Lpim_IM_mm_p"), h1-> Fill(im);
    if( header->trigmode2(Mode_KCDH2) ){
      h1 = (TH1F*)gFile->Get("Lpim_IM_mm_p_CDH2"), h1-> Fill(im);
    }
    if( header->IsTrig(Trig_Charged) ){
      h1 = (TH1F*)gFile->Get("Lpim_IM_mm_p_C"), h1-> Fill(im);
      if( beam->nFDC1()==1 ){
	h1 = (TH1F*)gFile->Get("Lpim_IM_mm_p_C_wFDC1"), h1-> Fill(im);
      }
    }
  }
}
