#include "MyHistReadCDS_Lpim.h"

using namespace std;

void initHistReadCDS_Lpim()
{
  new TH2F("CDS_IM_ppim_p_2pim_2D", "CDS p #pi^{-} IM vs p #pi^{-} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("CDS_IM_ppim_p_2pim_1",        "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_ppim_p_2pim_2",        "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_ppim_p_2pim_1_true",   "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_ppim_p_2pim_2_true",   "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_ppim_p_2pim_1_double", "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_ppim_p_2pim_2_double", "CDS p #pi^{-} IM", 1000, 1.0, 2.0);

  new TH1F("KLpim_MM",          "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_trigC",    "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_trigCDH3", "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_true_trigC",    "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KLpim_MM_true_trigCDH3", "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);

  new TH1F("Lpim_IM_trigC", "CDS #Lambda #pi^{-} IM", 2000, 0.0, 2.0);
  new TH1F("Lpim_IM_mmP_trigC", "CDS #Lambda #pi^{-} IM", 2000, 0.0, 2.0);
  new TH1F("Lpim_IM_mmP_trigCDH3", "CDS #Lambda #pi^{-} IM", 2000, 0.0, 2.0);

  new TH1F("Lpim_IM_mmP_trigC_whit", "CDS #Lambda #pi^{-} IM", 2000, 0.0, 2.0);
  new TH1F("Lpim_IM_mmP_trigC_wFP",  "CDS #Lambda #pi^{-} IM", 2000, 0.0, 2.0);

  new TH2F("Lpim_IM_KLpim_MM_trigC", "CDS #Lambda #pi^{-} IM vs d(K^{-}, #Lambda #pi^{-})\"X\"",    2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Lpim_IM_KLpim_MM_trigCDH3", "CDS #Lambda #pi^{-} IM vs d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH2F("Lpim_IM_p_cos_CM_trigC",    "CDS #Lambda #pi^{-} IM vs \"p\" cos#theta", 2000, 0.0, 2.0, 1000, -1.0, 1.0);
  new TH2F("Lpim_IM_p_cos_CM_trigCDH3", "CDS #Lambda #pi^{-} IM vs \"p\" cos#theta", 2000, 0.0, 2.0, 1000, -1.0, 1.0);

  new TH2F("Lpim_IM_KP_MM",      "CDS #Lambda #pi^{-} vs d(K^{-}, p)\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Lpim_IM_KP_MM_diff", "#Lambda #pi^{-} diff", 1000, -0.1, 1.0, 1000, 1.0, 2.0);
}

void fillHistReadCDS_Lpim(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;

  for( int i=0; i<anaInfo->nCDS2(CDS_PiMinus, CDS_Proton); i++ ){
    for( int j=1; j<anaInfo->nCDS2(CDS_PiMinus, CDS_Proton); j++ ){
      CDS2Info *ppim1=anaInfo->CDS2(CDS_PiMinus, CDS_Proton, i);
      CDS2Info *ppim2=anaInfo->CDS2(CDS_PiMinus, CDS_Proton, j);
      if( !ppim1-> flag() || !ppim2->flag() ) continue;
      int p_id1=ppim1->trackID1(); int pim_id1=ppim1->trackID2();
      if( ppim1->pid1()==CDS_PiMinus ) swap(p_id1, pim_id1);

      int p_id2=ppim2->trackID1(); int pim_id2=ppim2->trackID2();
      if( ppim2->pid1()==CDS_PiMinus ) swap(p_id2, pim_id2);

      if( pim_id1==pim_id2 ) continue;
      if( p_id1!=p_id2 ) continue;

      if( ppim1->displaced()>ppim2->displaced() ){
	swap(p_id1, p_id2); swap(pim_id1, pim_id2); swap(ppim1, ppim2);
      }
      if( GeomTools::GetID(ppim1->vertexBeam())!=CID_Fiducial || GeomTools::GetID(ppim2->vertexBeam())!=CID_Fiducial ) return;

      MyHistTools::fillTH("CDS_IM_ppim_p_2pim_2D", ppim1->im(), ppim2->im());
      MyHistTools::fillTH("CDS_IM_ppim_p_2pim_1", ppim1->im());
      MyHistTools::fillTH("CDS_IM_ppim_p_2pim_2", ppim2->im());
      bool L_flag1=false; bool L_flag2=false;
      if( L_MIN<ppim1->im() && ppim1->im()<L_MAX ){
	MyHistTools::fillTH("CDS_IM_ppim_p_2pim_1_true", ppim1->im());
	L_flag1=true;
      }
      if( L_MIN<ppim2->im() && ppim2->im()<L_MAX ){
	MyHistTools::fillTH("CDS_IM_ppim_p_2pim_2_true", ppim2->im());
	L_flag2=true;
      }

      if( L_flag1 && L_flag2 ){
	MyHistTools::fillTH("CDS_IM_ppim_p_2pim_1_double", ppim1->im());
	MyHistTools::fillTH("CDS_IM_ppim_p_2pim_2_double", ppim2->im());
      }

      CDS2Info *info_L=0; CDSInfo *info_pim=0;
      if( L_flag1 ){
	info_L=ppim1;
	info_pim=anaInfo->CDSbyID(pim_id2);
      }
      else if( L_flag2 ){
	info_L=ppim2;
	info_pim=anaInfo->CDSbyID(pim_id1);
      }

      if( info_L ){
	TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
	TLorentzVector target_lmom=MyAnaTools::target_lmom();

	TVector3 LabToCM=-(beam_lmom+target_lmom).BoostVector();
	TLorentzVector pim_lmom=info_pim->lmom();
	TLorentzVector L_lmom=info_L->lmom();
	TLorentzVector mm_p_lmom=beam_lmom+target_lmom-pim_lmom-L_lmom;
	TLorentzVector pim_lmom_CM=pim_lmom;
	TLorentzVector L_lmom_CM=L_lmom;
	TLorentzVector mm_p_lmom_CM=mm_p_lmom;
	pim_lmom_CM.Boost(LabToCM); L_lmom_CM.Boost(LabToCM); mm_p_lmom_CM.Boost(LabToCM);
	double p_cos_CM=mm_p_lmom_CM.CosTheta();

	double klpim_mm=(beam_lmom+target_lmom-info_L->lmom()-info_pim->lmom()).M();
	double lpim_im=(info_L->lmom()+info_pim->lmom()).M();
	MyHistTools::fillTH("KLpim_MM", klpim_mm);

	bool trigC=false;
	bool trigCDH3=false;
	if( header ){
	  if( header->trigmode2(Mode_KCDH1C) ) trigC=true;
	  if( header->trigmode2(Mode_KCDH3) ) trigCDH3=true;
	}
	else{
	  if( anaInfo->nFCharged()>0 ) trigC=true;
	  if( anaInfo->nCDS()>2 ) trigCDH3=true;
	}
	
	if( trigC ){
	  MyHistTools::fillTH("KLpim_MM_trigC", klpim_mm);
	  MyHistTools::fillTH("Lpim_IM_trigC", lpim_im);
	  MyHistTools::fillTH("Lpim_IM_KLpim_MM_trigC", lpim_im, klpim_mm);

	  if( 0.89<klpim_mm && klpim_mm<0.99 ){
	    MyHistTools::fillTH("KLpim_MM_true_trigC", klpim_mm);
	    MyHistTools::fillTH("Lpim_IM_p_cos_CM_trigC", lpim_im, p_cos_CM);
	    MyHistTools::fillTH("Lpim_IM_mmP_trigC", lpim_im);
	    
	    vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
	    vector<HodoscopeLikeHit*> PChits=MyTools::getHodo(blMan, CID_PC);
	    HodoscopeLikeHit* FChit=0;
	    if( CVChits.size()==1 && PChits.size()==0 ) FChit=CVChits[0];
	    if( PChits.size()==1 && CVChits.size()==0 ) FChit=PChits[0];
	    if( FChit ) MyHistTools::fillTH("Lpim_IM_mmP_trigC_whit", lpim_im);
	    if( anaInfo->nFCharge()==1 ){
	      ForwardChargeInfo *fcInfo = anaInfo->forwardCharge(0);
	      if( FC_P_MIN<fcInfo->mass2byRK() && fcInfo->mass2byRK()<FC_P_MAX ){
		TLorentzVector fp_lmom=fcInfo->lmom();
		double kp_mm=(beam_lmom+target_lmom-fp_lmom).M();
		
		MyHistTools::fillTH("Lpim_IM_mmP_trigC_wFP", lpim_im);
		MyHistTools::fillTH("Lpim_IM_KP_MM", lpim_im, kp_mm);
		MyHistTools::fillTH("Lpim_IM_KP_MM_diff", lpim_im-kp_mm, lpim_im);
	      }
	    }
	  }  
	}
	else if( trigCDH3 ){
	  MyHistTools::fillTH("KLpim_MM_trigCDH3", klpim_mm);
	  MyHistTools::fillTH("Lpim_IM_KLpim_MM_trigCDH3", lpim_im, klpim_mm);
	  
	  if( 0.89<klpim_mm && klpim_mm<0.99 ){
	    MyHistTools::fillTH("KLpim_MM_true_trigCDH3", klpim_mm);
	    MyHistTools::fillTH("Lpim_IM_p_cos_CM_trigCDH3", lpim_im, p_cos_CM);
	    MyHistTools::fillTH("Lpim_IM_mmP_trigCDH3", lpim_im);	  
	  }
	}
      }
    }
  }
}

