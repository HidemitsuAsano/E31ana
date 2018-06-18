#include "MyHistReadFC.h"

using namespace std;

void initHistReadFC()
{
  new TH2F("FC_overbeta_dE", "FC 1/#beta vs dE", 5000, 0.0, 5.0, 2000, 0, 200);

  new TH1F("FC_mass2",    "Forward Charge mass^{2}", 5100, -0.1, 5.0);
  new TH1F("FC_mass2_p",  "Forward Charge mass^{2}", 5100, -0.1, 5.0);
  new TH1F("FC_mass2_pi", "Forward Charge mass^{2}", 5100, -0.1, 5.0);

  new TH2F("FC_mass2_mom_ang", "Forward Charge mass^{2} vs mom", 5100, -0.1, 5.0, 1500, 0.0, 1.5);
  new TH2F("FC_mass2_mom_RK",  "Forward Charge mass^{2} vs mom", 5100, -0.1, 5.0, 1500, 0.0, 1.5);

  new TH2F("FC_r_mom", "FC radius vs mom by TOF", 5000, 0.0, 5000, 1500, 0.0, 1.5);

  new TH2F("FC_mom_ang_TOF", "FC mom ang vs TOF", 1500, 0.0, 1.5, 1500, 0.0, 1.5);
  new TH2F("FC_mom_RK_TOF", "FC mom Runge-Kutta vs TOF", 1500, 0.0, 1.5, 1500, 0.0, 1.5);

  new TH2F("FC_diff_mom_ang", "FC mom difference vs mom", 1000, -0.5, 0.5, 1500, 0.0, 1.5);
  new TH2F("FC_diff_mom_RK",  "FC mom difference vs mom", 1000, -0.5, 0.5, 1500, 0.0, 1.5);

  new TH2F("CDS_mass2_mom_wFP", "CDS mass2 vs mom", 10000, -1.0, 9.0, 3000, -1.5, 1.5);

  new TH1F("Ppim_IM", "p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Ppim_IM_angle", "p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Ppim_IM_RK", "p #pi^{-} IM", 1000, 1.0, 2.0);

  new TH1F("Pkm_IM",  "p K^{-} IM",   1000, 1.0, 2.0);
  new TH1F("Pkm_IM_angle",  "p K^{-} IM",   1000, 1.0, 2.0);
  new TH1F("Pkm_IM_RK",  "p K^{-} IM",   1000, 1.0, 2.0);

  new TH1F("FC_offset_mmN",  "FC_offset_mmN", 500, -10, 10);

  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("CVC%d_offset_RK", seg), "", 1000, -10.0, 10.0);
  }
  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("PC%d_offset_RK", seg), "", 1000, -10.0, 10.0);
  }
  string target_mat=DetectorList::GetInstance()->GetMaterial(CID_Fiducial);
  if( target_mat=="LHydrogen" ){
    new TH1F("KP_MM", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_pim", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_pip", "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_km",  "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_p",   "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_d",   "p(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_ppim",   "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);

    new TH1F("KPkm_MM",   "p(K^{-}, p K^{-})\"X\"",    1000,  0.0, 1.0);
    new TH1F("KPkm_MM_angle", "d(K^{-}, p)\"X\"", 1000, 0.0, 1.0);
    new TH1F("KPkm_MM_RK",  "d(K^{-}, p)\"X\"", 1000, 0.0, 1.0);

    new TH1F("KPpim_MM",  "p(K^{-}, p #pi^{-})\"X\"",  1000,  0.0, 1.0);
    new TH2F("Ppim_IM_KPpim_MM", "p(K^{-}, p #pi^{-})\"X\"",  1000, 1.0, 2.0, 1000,  0.0, 1.0);
    new TH1F("KPpim_MM_pip",  "p(K^{-}, p #pi^{-})\"X\"",  1000,  0.0, 1.0);
    new TH1F("KPpip_MM",  "p(K^{-}, p #pi^{+})\"X\"",  1000,  0.0, 1.0);
    new TH1F("KPp_MM2",   "p(K^{-}, p p)\"X\"^{2}",    2000, -1.0, 1.0);

    new TH2F("KP_MM_KPpim_MM", "p(K^{-}, p)\"X\" vs p(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
    new TH2F("KP_MM_KPpip_MM", "p(K^{-}, p)\"X\" vs p(K^{-}, p #pi^{+})\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
    new TH2F("KP_MM_KPkm_MM",  "p(K^{-}, p)\"X\" vs p(K^{-}, p K^{-})\"X\"",   2000, 0.0, 2.0, 1000,  0.0, 1.0);
    new TH2F("KP_MM_KPp_MM2",   "p(K^{-}, p)\"X\" vs p(K^{-}, p p)\"X\"^{2}",   2000, 0.0, 2.0, 2000, -1.0, 1.0);

    new TH2F("KPkm_MM_FDC1x", "p(K^{-}, K^{-} p)\"X\" vs FDC1 pos x", 1000, -50, 50,   1000,  0.0, 1.0);
    new TH2F("KPkm_MM_FDC1y", "p(K^{-}, K^{-} p)\"X\" vs FDC1 pos y", 1000, -50, 50,   1000,  0.0, 1.0);
    new TH2F("KPkm_MM_dx",    "p(K^{-}, K^{-} p)\"X\" vs dx/dz",      1000, -0.2, 0.2, 1000,  0.0, 1.0);
    new TH2F("KPkm_MM_dy",    "p(K^{-}, K^{-} p)\"X\" vs dy/dz",      1000, -0.2, 0.2, 1000,  0.0, 1.0);
    new TH2F("KPkm_MM_FCseg", "p(K^{-}, K^{-} p)\"X\" vs CVC/PC seg", 51, 0.5, 51.5,   1000,  0.0, 1.0);
    new TH2F("KPkm_MM_ang",   "p(K^{-}, K^{-} p)\"X\" vs cos#alpha",  1000, 0.0, 0.5,  1000,  0.0, 1.0);
    new TH2F("KPkm_MM_dE",    "p(K^{-}, K^{-} p)\"X\" vs dE",         1000, 0.0, 100,  1000,  0.0, 1.0);

    new TH2F("KP_MM_cos",     "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 1000, 0.0, 1.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_pim", "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 1000, 0.0, 1.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_pip", "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 1000, 0.0, 1.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_km",  "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 1000, 0.0, 1.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_p",   "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 1000, 0.0, 1.0, 1000, 0.9, 1.0);

    new TH2F("KP_MM_Ppim_IM",  "p(K^{-}, p)\"X\" vs p #pi^{-} IM", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
    new TH2F("KP_MM_Pkm_IM",  "p(K^{-}, p)\"X\" vs p K^{-} IM", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  }
  else if( target_mat=="LDeuterium" ){
    new TH1F("KP_MM", "d(K^{-}, p)\"X\"", 2000, 1.0, 3.0);
    new TH1F("KP_MM_pim", "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_pip", "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_km",  "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KPkm_MM_angle", "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KPkm_MM_RK",  "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);

    new TH1F("KP_MM_p",   "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_d",   "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
    new TH1F("KP_MM_ppim",   "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);

    new TH1F("KPkm_MM",   "d(K^{-}, K^{-} p)\"X\"",  2000,  0.0, 2.0);
    new TH1F("KPpim_MM",  "d(K^{-}, K^{-} p)\"X\"",  2000,  0.0, 2.0);
    new TH2F("Ppim_IM_KPpim_MM", "p(K^{-}, p #pi^{-})\"X\"",  1000, 0.0, 1.0, 2000,  0.0, 2.0);
    new TH1F("KPpim_MM_pip",  "p(K^{-}, p #pi^{-})\"X\"",  2000,  0.0, 2.0);
    new TH1F("KPpip_MM",  "d(K^{-}, K^{-} p)\"X\"",  2000,  0.0, 2.0);
    new TH1F("KPp_MM2",   "d(K^{-}, K^{-} p)\"X\"^{2}",  2000, -1.0, 1.0);

    new TH2F("KPkm_MM_FDC1x", "d(K^{-}, K^{-} p)\"X\" vs FDC1 pos x", 1000, -50, 50,   2000,  0.0, 2.0);
    new TH2F("KPkm_MM_FDC1y", "d(K^{-}, K^{-} p)\"X\" vs FDC1 pos y", 1000, -50, 50,   2000,  0.0, 2.0);
    new TH2F("KPkm_MM_dx",    "d(K^{-}, K^{-} p)\"X\" vs dx/dz",      1000, -0.2, 0.2, 2000,  0.0, 2.0);
    new TH2F("KPkm_MM_dy",    "d(K^{-}, K^{-} p)\"X\" vs dy/dz",      1000, -0.2, 0.2, 2000,  0.0, 2.0);
    new TH2F("KPkm_MM_FCseg", "d(K^{-}, K^{-} p)\"X\" vs CVC/PC seg", 51, 0.5, 51.5,   2000,  0.0, 2.0);
    new TH2F("KPkm_MM_ang",   "d(K^{-}, K^{-} p)\"X\" vs cos#alpha",  1000, 0.0, 0.5,  2000,  0.0, 2.0);
    new TH2F("KPkm_MM_dE",    "d(K^{-}, K^{-} p)\"X\" vs dE",         1000, 0.0, 100,  2000,  0.0, 2.0);

    new TH2F("KP_MM_KPpim_MM", "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KP_MM_KPpip_MM", "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KP_MM_KPkm_MM",  "d(K^{-}, p)\"X\" vs d(K^{-}, p K^{-})\"X\"",   2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KP_MM_KPp_MM2",   "d(K^{-}, p)\"X\" vs d(K^{-}, p p)\"X\"^{2}",   2000, 0.0, 2.0, 2000, -1.0, 1.0);

    new TH2F("KP_MM_cos",     "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_pim", "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_pip", "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_km",  "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000, 0.9, 1.0);
    new TH2F("KP_MM_cos_p",   "d(K^{-}, p)\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000, 0.9, 1.0);

    new TH2F("KP_MM_Ppim_IM",  "d(K^{-}, p)\"X\" vs p #pi^{-} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
    new TH2F("KP_MM_Pkm_IM",  "d(K^{-}, p)\"X\" vs p K^{-} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  }
  else if( target_mat=="LHelium-3" ){
    new TH1F("KP_MM", "^{3}He(K^{-}, p)\"X\"", 2000, 2.0, 4.0);
    new TH1F("KP_MM_pim", "^{3}He(K^{-}, p)\"X\"", 2000, 2.0, 4.0);
    new TH1F("KP_MM_pip", "^{3}He(K^{-}, p)\"X\"", 2000, 2.0, 4.0);
    new TH1F("KP_MM_km",  "^{3}He(K^{-}, p)\"X\"", 2000, 2.0, 4.0);
    new TH1F("KP_MM_p",   "^{3}He(K^{-}, p)\"X\"", 2000, 2.0, 4.0);
    new TH1F("KP_MM_d",   "^{3}He(K^{-}, p)\"X\"", 2000, 2.0, 4.0);
    new TH1F("KP_MM_ppim",   "d(K^{-}, p)\"X\"",   2000, 2.0, 4.0);

    new TH1F("KPkm_MM",   "^{3}He(K^{-}, K^{-} p)\"X\"",     2000,  1.0, 3.0);
    new TH1F("KPpim_MM",  "^{3}He(K^{-}, K^{-} p)\"X\"",     2000,  1.0, 3.0);
    new TH2F("Ppim_IM_KPpim_MM", "p(K^{-}, p #pi^{-})\"X\"",  1000, 0.0, 1.0, 2000,  1.0, 3.0);
    new TH1F("KPpim_MM_pip",  "p(K^{-}, p #pi^{-})\"X\"",    2000,  1.0, 3.0);
    new TH1F("KPpip_MM",  "^{3}He(K^{-}, K^{-} p)\"X\"",     2000,  1.0, 3.0);
    new TH1F("KPp_MM2",   "^{3}He(K^{-}, K^{-} p)\"X\"^{2}", 2000, -1.0, 1.0);

    new TH2F("KPkm_MM_FDC1x", "d(K^{-}, K^{-} p)\"X\" vs FDC1 pos x", 1000, -50, 50,   2000,  1.0, 3.0);
    new TH2F("KPkm_MM_FDC1y", "d(K^{-}, K^{-} p)\"X\" vs FDC1 pos y", 1000, -50, 50,   2000,  1.0, 3.0);
    new TH2F("KPkm_MM_dx",    "d(K^{-}, K^{-} p)\"X\" vs dx/dz",      1000, -0.2, 0.2, 2000,  1.0, 3.0);
    new TH2F("KPkm_MM_dy",    "d(K^{-}, K^{-} p)\"X\" vs dy/dz",      1000, -0.2, 0.2, 2000,  1.0, 3.0);
    new TH2F("KPkm_MM_FCseg", "d(K^{-}, K^{-} p)\"X\" vs CVC/PC seg", 51, 0.5, 51.5,   2000,  1.0, 3.0);
    new TH2F("KPkm_MM_ang",   "d(K^{-}, K^{-} p)\"X\" vs cos#alpha",  1000, 0.0, 0.5,  2000,  1.0, 3.0);
    new TH2F("KPkm_MM_dE",    "d(K^{-}, K^{-} p)\"X\" vs dE",         1000, 0.0, 100,  2000,  1.0, 3.0);

    new TH2F("KP_MM_KPpim_MM", "^{3}He(K^{-}, p)\"X\" vs ^{3}He(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KP_MM_KPpip_MM", "^{3}He(K^{-}, p)\"X\" vs ^{3}He(K^{-}, p #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KP_MM_KPkm_MM",  "^{3}He(K^{-}, p)\"X\" vs ^{3}He(K^{-}, p K^{-})\"X\"",   2000, 0.0, 2.0, 2000,  0.0, 2.0);
    new TH2F("KP_MM_KPp_MM2",  "^{3}He(K^{-}, p)\"X\" vs ^{3}He(K^{-}, p p)\"X\"^{2}",   2000, 0.0, 2.0, 2000, -1.0, 1.0);

    new TH2F("KP_MM_Ppim_IM",  "^{3}He(K^{-}, p)\"X\" vs p #pi^{-} IM", 2000, 1.0, 3.0, 1000, 1.0, 2.0);
    new TH2F("KP_MM_Pkm_IM",  "^{3}He(K^{-}, p)\"X\" vs p K^{-} IM", 2000, 1.0, 3.0, 1000, 1.0, 2.0);
  }
  else{
    cout<<" !!! unknown target material "<<target_mat<<" !!! "<<endl;
    exit(0);
  }
}

void fillHistReadFC(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 ) return;
  if( !anaInfo->beam(0)->flag() ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
  if( anaInfo->nFCharge()!=1 ) return;

  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  if( BVChits.size()==0 ) return;

  ForwardChargeInfo *fcInfo = anaInfo->forwardCharge(0);
  HodoscopeLikeHit *hit=fcInfo->hodo(blMan);
  if( !hit ){
    fcInfo-> dump();
  }
  else{
    MyHistTools::fillTH("FC_overbeta_dE", 1./fcInfo->beta(), hit->emean());
  }


  MyHistTools::fillTH("FC_mass2_mom_ang", fcInfo->mass2byAng(), fcInfo->momByAng());
  MyHistTools::fillTH("FC_mass2_mom_RK", fcInfo->mass2byRK(), fcInfo->momByRK());

  MyHistTools::fillTH("FC_mass2", fcInfo->mass2byAng());
  if( 0.4<fcInfo->mass2byAng() && fcInfo->mass2byAng()<2.0 ){
    //    cout<<" mom by RK : "<<fcInfo->momByRK()<<"  by Ang : "<<fcInfo->momByAng()<<endl;
    MyHistTools::fillTH("FC_mass2_p", fcInfo->mass2byAng());

    MyHistTools::fillTH("FC_r_mom", fcInfo->radius(conf), fcInfo->momByTOF());
    MyHistTools::fillTH("FC_mom_ang_TOF", fcInfo->momByAng(), fcInfo->momByTOF());
    MyHistTools::fillTH("FC_mom_RK_TOF", fcInfo->momByRK(), fcInfo->momByTOF());

    MyHistTools::fillTH("FC_diff_mom_ang", fcInfo->momByAng()-fcInfo->momByTOF(), fcInfo->momByTOF());
    MyHistTools::fillTH("FC_diff_mom_RK", fcInfo->momByRK()-fcInfo->momByTOF(), fcInfo->momByTOF());

    int FCseg=fcInfo->seg();
    if( fcInfo->isPC() ) FCseg+=27;
    HodoscopeLikeHit *hit =fcInfo->hodo(blMan);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector target_lmom=MyAnaTools::target_lmom();
    TLorentzVector fc_lmom=fcInfo->lmom();
    TVector3 fc_mom=fc_lmom.Vect();
    TVector3 fc_mom_ang=fc_mom;
    TVector3 fc_mom_RK=fc_mom;
    fc_mom_ang.SetMag(fcInfo->momByAng());
    fc_mom_RK.SetMag(fcInfo->momByRK());
    TLorentzVector fc_lmom_ang, fc_lmom_RK;
    fc_lmom_ang.SetVectM(fc_mom_ang, pMass);
    fc_lmom_RK.SetVectM(fc_mom_RK, pMass);

    double cos_theta=beam_lmom.Vect().Dot(fc_lmom.Vect())/(beam_lmom.Vect().Mag()*fc_lmom.Vect().Mag());

    TVector3 FDC1pos=fcInfo->FDC1().pos();
    TVector3 vtx=fcInfo->vertex();
    TVector3 dir=fc_lmom.Vect().Unit();

    double beam_out, beam_tof;
    ELossTools::CalcElossBeamTGeo(anaInfo->beam(0)->T0pos(), anaInfo->minDCA()->vertexBeam(), 
				  anaInfo->beam(0)->D5mom(), anaInfo->beam(0)->mass(), beam_out, beam_tof);

    double tof=fcInfo->time()-anaInfo->beam(0)->T0time();
    double calc_offset=tof-beam_tof-fcInfo->calc_tof();
    //    cout<<" Calc offset : "<<calc_offset<<" = "<<tof<<" - "<<beam_tof<<" - "<<fcInfo->calc_tof()<<endl;
    if( fcInfo->isPC() ){
      MyHistTools::fillTH(Form("PC%d_offset_RK", hit->seg()), calc_offset);
    }
    else{
      MyHistTools::fillTH(Form("CVC%d_offset_RK", hit->seg()), calc_offset);      
    }

    double mm=(beam_lmom+target_lmom-fc_lmom).M();
    MyHistTools::fillTH("KP_MM", mm);
    MyHistTools::fillTH("KP_MM_cos", mm, cos_theta);
    if( anaInfo->nCDS(CDS_PiMinus)>0 ){
      MyHistTools::fillTH("KP_MM_pim", mm);
      MyHistTools::fillTH("KP_MM_cos_pim", mm, cos_theta);
      if( anaInfo->nCDS(CDS_Proton)>0 ){
	MyHistTools::fillTH("KP_MM_ppim", mm);
      }
    }
    if( anaInfo->nCDS(CDS_PiPlus)>0 ){
      MyHistTools::fillTH("KP_MM_pip", mm);
      MyHistTools::fillTH("KP_MM_cos_pip", mm, cos_theta);
    }
    if( anaInfo->nCDS(CDS_Kaon)>0 ){
      MyHistTools::fillTH("KP_MM_km", mm);
      MyHistTools::fillTH("KP_MM_cos_km", mm, cos_theta);
    }
    if( anaInfo->nCDS(CDS_Proton)>0 ){
      MyHistTools::fillTH("KP_MM_p", mm);
      MyHistTools::fillTH("KP_MM_cos_p", mm, cos_theta);
    }
    if( anaInfo->nCDS(CDS_Deuteron)>0 ){
      MyHistTools::fillTH("KP_MM_d", mm);
    }

    for( int i=0; i<anaInfo->nCDS(CDS_PiMinus); i++ ){
      CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, i);
      TLorentzVector pim_lmom=pim->lmom();
      double im=(fc_lmom+pim_lmom).M();
      double mm_pim=(beam_lmom+target_lmom-fc_lmom-pim_lmom).M();
      MyHistTools::fillTH("Ppim_IM", im);
      MyHistTools::fillTH("Ppim_IM_KPpim_MM", im, mm_pim);
      MyHistTools::fillTH("Ppim_IM_angle", (fc_lmom_ang+pim_lmom).M());
      MyHistTools::fillTH("Ppim_IM_RK", (fc_lmom_RK+pim_lmom).M());

      MyHistTools::fillTH("KPpim_MM", mm_pim);
      if( anaInfo->nCDS(CDS_PiPlus)>0 ) MyHistTools::fillTH("KPpim_MM_pip", mm_pim);
      MyHistTools::fillTH("KP_MM_KPpim_MM", mm, mm_pim);

      MyHistTools::fillTH("KP_MM_Ppim_IM", mm, (fc_lmom+pim_lmom).M());
    }

    for( int i=0; i<anaInfo->nCDS(CDS_PiPlus); i++ ){
      CDSInfo *pip=anaInfo->CDS(CDS_PiPlus, i);
      TLorentzVector pip_lmom=pip->lmom();
      double im=(fc_lmom+pip_lmom).M();
      double mm_pip=(beam_lmom+target_lmom-fc_lmom-pip_lmom).M();
      MyHistTools::fillTH("KPpip_MM", mm_pip);
      MyHistTools::fillTH("KP_MM_KPpip_MM", mm, mm_pip);
    }

    for( int i=0; i<anaInfo->nCDS(CDS_Kaon); i++ ){
      CDSInfo *km=anaInfo->CDS(CDS_Kaon, i);
      TLorentzVector km_lmom=km->lmom();

      double im=(fc_lmom+km_lmom).M();
      double mm_km=(beam_lmom+target_lmom-fc_lmom-km_lmom).M();
      double mm_km_ang=(beam_lmom+target_lmom-fc_lmom_ang-km_lmom).M();
      double mm_km_RK=(beam_lmom+target_lmom-fc_lmom_RK-km_lmom).M();
      MyHistTools::fillTH("Pkm_IM", im);
      MyHistTools::fillTH("Pkm_IM_angle", (fc_lmom_ang+km_lmom).M());
      MyHistTools::fillTH("Pkm_IM_RK", (fc_lmom_RK+km_lmom).M());

      MyHistTools::fillTH("KPkm_MM", mm_km);
      MyHistTools::fillTH("KPkm_MM_angle", mm_km_ang);
      MyHistTools::fillTH("KPkm_MM_RK", mm_km_RK);
      MyHistTools::fillTH("KP_MM_KPkm_MM", mm, mm_km);

      MyHistTools::fillTH("KP_MM_Pkm_IM", mm, im);

      MyHistTools::fillTH("KPkm_MM_FDC1x", FDC1pos.X(), mm_km);
      MyHistTools::fillTH("KPkm_MM_FDC1y", FDC1pos.Y(), mm_km);
      MyHistTools::fillTH("KPkm_MM_dx", dir.X()/dir.Z(), mm_km);
      MyHistTools::fillTH("KPkm_MM_dy", dir.Y()/dir.Z(), mm_km);
      MyHistTools::fillTH("KPkm_MM_FCseg", FCseg, mm_km);
      MyHistTools::fillTH("KPkm_MM_ang", fcInfo->angle(conf), mm_km);
      MyHistTools::fillTH("KPkm_MM_dE", hit->emean(), mm_km);

      if( 0.9<mm_km && mm_km<0.98 ){
	double offset=MyTools::calc_offset(nMass, beam_lmom+target_lmom-km_lmom, fc_lmom, anaInfo);
	MyHistTools::fillTH("FC_offset_mmN", offset);
      }
    }

    for( int i=0; i<anaInfo->nCDS(CDS_Proton); i++ ){
      CDSInfo *p=anaInfo->CDS(CDS_Proton, i);
      TLorentzVector p_lmom=p->lmom();
      double im=(fc_lmom+p_lmom).M();
      double mm_p=(beam_lmom+target_lmom-fc_lmom-p_lmom).M2();
      MyHistTools::fillTH("KPp_MM2", mm_p);
      MyHistTools::fillTH("KP_MM_KPp_MM2", mm, mm_p);
    }   


    for( int i=0; i<anaInfo->nCDS(); i++ ){
      if( !anaInfo->CDS(i)->flag() ) continue;
      if( GeomTools::GetID(anaInfo->CDS(i)->vertexBeam())!=CID_Fiducial ) continue;
      CDSInfo *cds=anaInfo->CDS(i);
      MyHistTools::fillTH("CDS_mass2_mom_wFP", cds->mass2(), cds->mom());
    }
  }
}
