#include "MyHistTools.h"
#include "MyParam.h"

using namespace std;

void MyHistTools::initFC()
{
  new TH1F("FC_mass2", "FC mass^{2}", 5000, -0.5, 4.5);
  new TH2F("FC_mass2_mom", "FC mass^{2} vs mom", 5000, -0.5, 4.5, 2000, 0.0, 2.0);
  new TH2F("FC_mom_USWK_TOF", "mom by USWK vs mom by TOF", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH1F("KP_MM", "KP_MM", 5000, 0.0, 5.0);
  
  new TH1F("FC_hitpat_p", "CVC/PC proton hit pattern", 61, 0.5, 61.5);
  new TH1F("FC_hitpat_n", "CVC/PC proton hit pattern", 61, 0.5, 61.5);
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("CVC_%d_offset_p", seg), Form("CVC seg%d offsetp", seg), 1000, -25, 25);
    new TH1F(Form("CVC_%d_offset_n", seg), Form("CVC seg%d offsetp", seg), 1000, -25, 25);
    new TH2F(Form("FC_mass2_mom_CVC%d", seg), "FC mass^{2} vs mom", 5000, -0.5, 4.5, 2000, 0.0, 2.0);
  }
  for( int seg=1; seg<=27; seg++ ){
    new TH1F(Form("PC_%d_offset_p", seg), Form("PC seg%d offsetp", seg), 1000, -25, 25);
    new TH1F(Form("PC_%d_offset_n", seg), Form("PC seg%d offsetp", seg), 1000, -25, 25);
    new TH2F(Form("FC_mass2_mom_PC%d", seg), "FC mass^{2} vs mom", 5000, -0.5, 4.5, 2000, 0.0, 2.0);
  }
}

void MyHistTools::fillFC(ConfMan *conf, EventHeader *header, AnaInfo *anaInfo, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan)
{
  TH1F *h1;
  TH2F *h2;

  std::vector<HodoscopeLikeHit*> BVC_hits;
  for( int i=0; i<blMan->nBVC(); i++ ){ if( blMan->BVC(i)->CheckRange() ) BVC_hits.push_back(blMan->BVC(i)); }

  std::vector<HodoscopeLikeHit*> PC_hits;
  for( int i=0; i<blMan->nPC(); i++ ){ if( blMan->PC(i)->CheckRange() ) PC_hits.push_back(blMan->PC(i)); }

  std::vector<HodoscopeLikeHit*> CVC_hits;
  for( int i=0; i<blMan->nCVC(); i++ ){ if( blMan->CVC(i)->CheckRange() ) CVC_hits.push_back(blMan->CVC(i)); }

  BeamInfo *beam=anaInfo->beam(0);
  CDSInfo *minVtx=anaInfo->minDCA();
  // cout<<" nFDC : "<<beam->nFDC1()<<endl;
  // cout<<" nBVC : "<<BVC_hits.size()<<endl;
  //  if( beam->nFDC1()==0 && BVC_hits.size()==0 ){
  if( BVC_hits.size()==0 ){
    HodoscopeLikeHit *fc_hit=0;
    if( PC_hits.size()==1 && CVC_hits.size()==0 ) fc_hit=PC_hits[0];
    if( PC_hits.size()==0 && CVC_hits.size()==1 ) fc_hit=CVC_hits[0];
    if( fc_hit ){
      TVector3 pos;
      conf->GetGeomMapManager()->GetGPos(fc_hit->cid(), fc_hit->seg(), pos);
      double beam_out, beam_tof;
      ELossTools::CalcElossBeamTGeo(beam->T0pos(), minVtx->vertexBeam(), beam->D5mom(), beam->mass(), beam_out, beam_tof);
      double tof=fc_hit->ctmean()-beam->T0time()-beam_tof;
      double fl=(pos-minVtx->vertexBeam()).Mag();
      double calc_tof=fl/(100.*Const);
      double offset=tof-calc_tof;
      if( fc_hit->cid()==CID_CVC ){
	h1 = (TH1F*)gFile->Get("FC_hitpat_n"), h1->Fill(fc_hit->seg());
	h1 = (TH1F*)gFile->Get(Form("CVC_%d_offset_n", fc_hit->seg())), h1->Fill(offset);
      }
      else if( fc_hit->cid()==CID_PC ){
	h1 = (TH1F*)gFile->Get("FC_hitpat_n"), h1->Fill(34+fc_hit->seg());
	h1 = (TH1F*)gFile->Get(Form("PC_%d_offset_n", fc_hit->seg())), h1->Fill(offset);
      }
    }
  }

  if( anaInfo->nFCharge()==1 ){
    std::vector<HodoscopeLikeHit*> BVC_hits;
    for( int i=0; i<blMan->nBVC(); i++ ){ if( blMan->BVC(i)->CheckRange() ) BVC_hits.push_back(blMan->BVC(i)); }
    
    if( BVC_hits.size()>0 ){
      int seg=anaInfo->forwardCharge(0)->seg();
      double mass2=anaInfo->forwardCharge(0)->mass2();
      double mom=anaInfo->forwardCharge(0)->mom();
      h1 = (TH1F*)gFile-> Get("FC_mass2"), h1-> Fill(mass2);
      h2 = (TH2F*)gFile-> Get("FC_mass2_mom"), h2-> Fill(mass2, mom);
      if( anaInfo->forwardCharge(0)->isPC() ){
	h2 = (TH2F*)gFile-> Get(Form("FC_mass2_mom_PC%d", seg)), h2-> Fill(mass2, mom);
      }else{
	h2 = (TH2F*)gFile-> Get(Form("FC_mass2_mom_CVC%d", seg)), h2-> Fill(mass2, mom);
      }
      if( anaInfo->forwardCharge(0)->pid()==F_Proton ){
	double mom_by_tof=anaInfo->forwardCharge(0)->mom_by_tof();
	double offset=anaInfo->forwardCharge(0)->offset();
	h2 = (TH2F*)gFile-> Get("FC_mom_USWK_TOF"), h2-> Fill(mom, mom_by_tof);
	if( anaInfo->forwardCharge(0)->isPC() ){
	  h1 = (TH1F*)gFile->Get("FC_hitpat_p"), h1-> Fill(34+seg);
	  h1 = (TH1F*)gFile->Get(Form("PC_%d_offset_p", seg)), h1->Fill(offset);
	}else{
	  h1 = (TH1F*)gFile->Get("FC_hitpat_p"), h1-> Fill(seg);
	  h1 = (TH1F*)gFile->Get(Form("CVC_%d_offset_p", seg)), h1->Fill(offset);
	}
      }
    }
  }
}
