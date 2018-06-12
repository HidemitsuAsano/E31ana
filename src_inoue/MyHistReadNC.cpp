#include "MyHistReadNC.h"

using namespace std;

void initHistReadNC()
{
  new TH2F("NC_overbeta_dE", "NC 1/#beta vs dE", 5000, 0.0, 5.0, 200, 0.0, 200);
  new TH1F("NC_overbeta_gamma", "", 5000, 0.9, 1.1);
  new TH1F("NC_offset", "", 1000, -1.0, 1.0);
  for( int seg=1; seg<=112; seg++ ){
    new TH1F(Form("NC%d_overbeta_gamma",seg), "", 5000, 0.9, 1.1);
    new TH1F(Form("NC%d_offset",seg), "", 1000, -1.0, 1.0);
  }


  string target_mat=DetectorList::GetInstance()->GetMaterial(CID_Fiducial);
  if( target_mat=="LHydrogen" ){
    new TH1F("KN_MM", "p(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_pim", "p(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_pip", "p(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_km",  "p(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_p",   "p(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_d",   "p(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 

    new TH1F("KNpim_MM", "p(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F("KNpip_MM", "p(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0);
    new TH1F("KNkm_MM",  "p(K^{-}, n K^{-})\"X\"",   2000, 0.0, 2.0);
    new TH1F("KNp_MM2",  "p(K^{-}, n p)\"X\"",       2000, -1.0, 1.0);
    new TH1F("KNp_MM",  "p(K^{-}, n p)\"X\"",        2000, 0.0, 2.0);

    new TH2F("KN_MM_KNpim_MM", "p(K^{-}, n)\"X\" vs p(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
    new TH2F("KN_MM_KNpip_MM", "p(K^{-}, n)\"X\" vs p(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
    new TH2F("KN_MM_KNkm_MM",  "p(K^{-}, n)\"X\" vs p(K^{-}, n K^{-})\"X\"",   2000, 0.0, 2.0, 1000,  0.0, 1.0);
    new TH2F("KN_MM_KNp_MM2",  "p(K^{-}, n)\"X\" vs p(K^{-}, n p)\"X\"^{2}",   2000, 0.0, 2.0, 2000, -1.0, 1.0);
  }
  else if( target_mat=="LDeuterium" ){
    new TH1F("KN_MM", "d(K^{-}, n)\"X\"", 2000, 1.0, 3.0); 
    new TH1F("KN_MM_pim", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_pip", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_km",  "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_p",   "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 
    new TH1F("KN_MM_d",   "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0); 

    new TH1F("KNpim_MM", "d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0);
    new TH1F("KNpip_MM", "d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0);
    new TH1F("KNkm_MM",  "d(K^{-}, n K^{-})\"X\"",   2000, 0.0, 2.0);
    new TH1F("KNp_MM2",  "d(K^{-}, n p)\"X\"",       2000, -1.0, 1.0);
    new TH1F("KNp_MM",  "p(K^{-}, n p)\"X\"",        2000, 0.0, 2.0);

    new TH2F("KN_MM_KNpim_MM", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KN_MM_KNpip_MM", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KN_MM_KNkm_MM",  "d(K^{-}, n)\"X\" vs d(K^{-}, n K^{-})\"X\"",   2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KN_MM_KNp_MM2",  "d(K^{-}, n)\"X\" vs d(K^{-}, n p)\"X\"^{2}",   2000, 0.0, 2.0, 2000, -1.0, 1.0);
  }
  else if( target_mat=="LHelium-3" ){
    new TH1F("KN_MM", "^{3}He(K^{-}, n)\"X\"", 2000, 2.0, 4.0); 
    new TH1F("KN_MM_pim", "^{3}He(K^{-}, n)\"X\"", 2000, 2.0, 4.0); 
    new TH1F("KN_MM_pip", "^{3}He(K^{-}, n)\"X\"", 2000, 2.0, 4.0); 
    new TH1F("KN_MM_km",  "^{3}He(K^{-}, n)\"X\"", 2000, 2.0, 4.0); 
    new TH1F("KN_MM_p",   "^{3}He(K^{-}, n)\"X\"", 2000, 2.0, 4.0); 
    new TH1F("KN_MM_d",   "^{3}He(K^{-}, n)\"X\"", 2000, 2.0, 4.0); 

    new TH1F("KNpim_MM", "p(K^{-}, n #pi^{-})\"X\"", 2000, 1.0, 3.0);
    new TH1F("KNpip_MM", "p(K^{-}, n #pi^{+})\"X\"", 2000, 1.0, 3.0);
    new TH1F("KNkm_MM",  "p(K^{-}, n K^{-})\"X\"",   2000, 1.0, 3.0);
    new TH1F("KNp_MM2",  "p(K^{-}, n p)\"X\"",       2000, 0.0, 2.0);
    new TH1F("KNp_MM",  "p(K^{-}, n p)\"X\"",        2000, 0.0, 2.0);

    new TH2F("KN_MM_KNpim_MM", "^{3}He(K^{-}, n)\"X\" vs ^{3}He(K^{-}, n #pi^{-})\"X\"", 2000, 2.0, 4.0, 2000,  1.0, 3.0);
    new TH2F("KN_MM_KNpip_MM", "^{3}He(K^{-}, n)\"X\" vs ^{3}He(K^{-}, n #pi^{+})\"X\"", 2000, 2.0, 4.0, 2000,  1.0, 3.0);
    new TH2F("KN_MM_KNkm_MM",  "^{3}He(K^{-}, n)\"X\" vs ^{3}He(K^{-}, n K^{-})\"X\"",   2000, 2.0, 4.0, 2000,  1.0, 3.0);
    new TH2F("KN_MM_KNp_MM2",  "^{3}He(K^{-}, n)\"X\" vs ^{3}He(K^{-}, n p)\"X\"^{2}",   2000, 2.0, 4.0, 2000,  0.0, 2.0);
  }
  else{
    cout<<" !!! unknown target material "<<target_mat<<" !!! "<<endl;
    exit(0);
  }
}

void fillHistReadNC(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo)
{
  vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);

  ForwardNeutralInfo tmp_fnInfo=MyTools::makeFN(blMan, anaInfo, 0.0, false);
  if( tmp_fnInfo.pid()==F_Gamma || tmp_fnInfo.pid()==F_Neutron ){
    if( CVChits.size()==0 && BVChits.size()==0 ){
      MyHistTools::fillTH("NC_overbeta_dE", 1./tmp_fnInfo.beta(), tmp_fnInfo.NC(blMan)->emean());
    }
  }

  if( anaInfo->nFNeutral()!=1 ) return;
  if( header ){
    if( !header->IsTrig(Trig_Neutral) ) return;
  }

  if( CVChits.size()!=0 || BVChits.size()!=0 ) return;
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;

  BeamInfo *beam=anaInfo->beam(0);
  ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
  HodoscopeLikeHit *NChit=fnInfo->NC(blMan);

  if( fnInfo->pid()==F_Gamma ){
    double offset=fnInfo->time()-beam->T0time()-fnInfo->offset();
    MyHistTools::fillTH("NC_overbeta_gamma", 1./fnInfo->beta());
    MyHistTools::fillTH("NC_offset", offset);
    MyHistTools::fillTH(Form("NC%d_overbeta_gamma", NChit->seg()), 1./fnInfo->beta());
    MyHistTools::fillTH(Form("NC%d_offset", NChit->seg()), offset);
  }

  if( fnInfo->pid()==F_Neutron ){
    TLorentzVector beam_lmom=beam->lmom();
    TLorentzVector target_lmom=MyAnaTools::target_lmom();
    TLorentzVector fn_lmom=fnInfo->lmom();

    double mm=(beam_lmom+target_lmom-fn_lmom).M();
    MyHistTools::fillTH("KN_MM", mm);
    if( anaInfo->nCDS(CDS_PiMinus)>0 ){
      MyHistTools::fillTH("KN_MM_pim", mm);
    }
    else if( anaInfo->nCDS(CDS_PiPlus)>0 ){
      MyHistTools::fillTH("KN_MM_pip", mm);
    }
    else if( anaInfo->nCDS(CDS_Kaon)>0 ){
      MyHistTools::fillTH("KN_MM_km", mm);
    }
    else if( anaInfo->nCDS(CDS_Proton)>0 ){
      MyHistTools::fillTH("KN_MM_p", mm);
    }
    else if( anaInfo->nCDS(CDS_Deuteron)>0 ){
      MyHistTools::fillTH("KN_MM_d", mm);
    }

    for( int i=0; i<anaInfo->nCDS(CDS_PiMinus); i++ ){
      CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, i);
      TLorentzVector pim_lmom=pim->lmom();
      double mm_pim=(beam_lmom+target_lmom-fn_lmom-pim_lmom).M();
      MyHistTools::fillTH("KNpim_MM", mm_pim);
      MyHistTools::fillTH("KN_MM_KNpim_MM", mm, mm_pim);
    }

    for( int i=0; i<anaInfo->nCDS(CDS_PiPlus); i++ ){
      CDSInfo *pip=anaInfo->CDS(CDS_PiPlus, i);
      TLorentzVector pip_lmom=pip->lmom();
      double mm_pip=(beam_lmom+target_lmom-fn_lmom-pip_lmom).M();
      MyHistTools::fillTH("KNpip_MM", mm_pip);
      MyHistTools::fillTH("KN_MM_KNpip_MM", mm, mm_pip);
    }

    for( int i=0; i<anaInfo->nCDS(CDS_Kaon); i++ ){
      CDSInfo *km=anaInfo->CDS(CDS_Kaon, i);
      TLorentzVector km_lmom=km->lmom();
      double mm_km=(beam_lmom+target_lmom-fn_lmom-km_lmom).M();
      MyHistTools::fillTH("KNkm_MM", mm_km);
      MyHistTools::fillTH("KN_MM_KNkm_MM", mm, mm_km);
    }

    for( int i=0; i<anaInfo->nCDS(CDS_Proton); i++ ){
      CDSInfo *p=anaInfo->CDS(CDS_Proton, i);
      TLorentzVector p_lmom=p->lmom();
      double im=(fn_lmom+p_lmom).M();
      double mm_p=(beam_lmom+target_lmom-fn_lmom-p_lmom).M2();
      MyHistTools::fillTH("KNp_MM2", mm_p);
      MyHistTools::fillTH("KNp_MM", (beam_lmom+target_lmom-fn_lmom-p_lmom).M());
      MyHistTools::fillTH("KN_MM_KNp_MM2", mm, mm_p);
    }
  }
}
