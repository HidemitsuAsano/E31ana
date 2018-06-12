#include "MyHistMC_pLpim.h"

using namespace std;

void initHistMC_pLpim()
{
  new TH2F("gen_Lpim_IM_p_cos_reac1600", "Generate Lpim_IM_p_cos_reac1600", 2000, 0.0, 2.0, 1000, -1.0, 1.0);
  new TH2F("acc_Lpim_IM_p_cos_reac1600", "Accepted Lpim_IM_p_cos_reac1600", 2000, 0.0, 2.0, 1000, -1.0, 1.0);

}


void fillHistMC_pLpim(DetectorData *detData, MCData* mcData, ReactionData *reacData, AnaInfo *anaInfo,
                      BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan)
{
  if( reacData->ReactionID()==1600 ){
    TLorentzVector beam_lmom_MC;
    TLorentzVector target_MC;
    for( int i=0; i<reacData->InitParticleSize(); i++ ){
      if( reacData->InitPDG(i)==-321 ) beam_lmom_MC=0.001*reacData->GetInitParticle(i);
      else target_MC=0.001*reacData->GetInitParticle(i);
    }

    TLorentzVector p_lmom_CM_MC;
    TLorentzVector L_lmom_CM_MC;
    TLorentzVector pim_lmom_CM_MC;
    for( int i=0; i<reacData->ParticleSize(); i++ ){
      if( reacData->PDG(i)==2212 ) p_lmom_CM_MC=0.001*reacData->GetCMParticle(i);
      if( reacData->PDG(i)==3122 ) L_lmom_CM_MC=0.001*reacData->GetCMParticle(i);
      if( reacData->PDG(i)==-211 ) pim_lmom_CM_MC=0.001*reacData->GetCMParticle(i);
    }

    double p_cos=p_lmom_CM_MC.CosTheta();
    double Lpim_im=(L_lmom_CM_MC+pim_lmom_CM_MC).M();

    if( anaInfo->nBeam()!=1 ) return;

    bool acc_flag=false;
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

	bool L_flag1=false; bool L_flag2=false;
	if( L_MIN<ppim1->im() && ppim1->im()<L_MAX ){
	  L_flag1=true;
	}
	if( L_MIN<ppim2->im() && ppim2->im()<L_MAX ){
	  L_flag2=true;
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

	  if( 0.89<klpim_mm && klpim_mm<0.99 ){
	    acc_flag=true;
	  }
	}
      }
    }
    MyHistTools::fillTH("gen_Lpim_IM_p_cos_reac1600", Lpim_im, p_cos);
    if( acc_flag ){
      MyHistTools::fillTH("acc_Lpim_IM_p_cos_reac1600", Lpim_im, p_cos);
    }
  }
}
