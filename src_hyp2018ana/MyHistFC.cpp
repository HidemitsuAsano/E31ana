#include "MyHistFC.h"

using namespace std;

// static const double FC_P_MAX= 0.888345+2.*0.0454692;
// static const double FC_P_MIN= 0.888345-2.*0.0454692;

static const double FC_P_MAX= 0.833710+2.*0.0454692;
static const double FC_P_MIN= 0.833710-2.*0.0454692;

static const double FC_mmL_MIN= 1.11722-2*0.0170564;
static const double FC_mmL_MAX= 1.11722+2*0.0170564;
static const double FC_mmS0_MIN=1.19667-2*0.0118926;
static const double FC_mmS0_MAX=1.19667+2*0.0118926;
static const double FC_mmP_MIN=0.939025-2*0.0183338;
static const double FC_mmP_MAX=0.939025+2*0.0183338;

void initHistFC()
{
  new TH1F("FC_overbeta", "Forward Charged 1/#beta", 1000, -0.5, 9.5);
  new TH2F("FC_overbeta_dE", "Forward Charged 1/#beta vs dE", 1000, -0.5, 9.5, 100, 0, 100);

  new TH2F("FC_mass2_mom_ang", "Forward Charged mass2 vs mom", 1000, -0.5, 9.5, 1500, 0, 1.5);
  new TH2F("FC_mass2_mom_RK",  "Forward Charged mass2 vs mom", 1000, -0.5, 9.5, 1500, 0, 1.5);

  new TH2F("FC_mom_RK_TOF_p", "FC mom Runge-Kutta vs TOF", 1000, 0.5, 1.5,  1000, 0.5, 1.5);
  new TH2F("FC_mom_ang_TOF_p", "FC mom Runge-Kutta vs TOF", 1000, 0.5, 1.5, 1000, 0.5, 1.5);
  new TH2F("FC_ang_mom_TOF_p", "FC mom Runge-Kutta vs TOF", 1000, 0.0, 0.1, 1000, 0.5, 1.5);

  new TH1F("KP_MM",      "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_km",   "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_pim",  "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_pip",  "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_p",    "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_ppim", "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_2pim", "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_2pim_mmP",         "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_2pim_mmPgamma",    "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_2pim_mmP_wL",      "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KP_MM_2pim_mmPgamma_wS0", "d(K^{-}, p)\"X\"", 2000, 0.0, 2.0);

  new TH1F("KPpimpim_MM",      "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_mmL",  "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_mmS0", "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH2F("KPpim_KPpim_MM",          "d(K^{-}, p #pi^{-})\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 200, 1.0, 2.0, 200, 1.0, 2.0);
  new TH2F("KPpim_KPpim_MM_mmP",      "d(K^{-}, p #pi^{-})\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 200, 1.0, 2.0, 200, 1.0, 2.0);
  new TH2F("KPpim_KPpim_MM_mmPgamma", "d(K^{-}, p #pi^{-})\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 200, 1.0, 2.0, 200, 1.0, 2.0);
  new TH1F("KPpim_MM_mmP",      "d(K^{-}, p #pi^{-})\"X\"", 1000, 1.0, 2.0);
  new TH1F("KPpim_MM_mmPgamma", "d(K^{-}, p #pi^{-})\"X\"", 1000, 1.0, 2.0);

  new TH1F("Ppim_IM", "p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH2F("KPpim_MM_Ppim_IM", "d(K^{-}, p #pi^{-})\"X\" vs p #pi^{-} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpip_MM_Ppip_IM", "d(K^{-}, p #pi^{+})\"X\" vs p #pi^{+} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);

  new TH2F("KPkm_MM_Pkm_IM",  "d(K^{-}, p K^{-})\"X\" vs p K^{-} IM",  2000, 0.0, 2.0, 1000, 1.0, 2.0);
  for( int i=1; i<=61; i++ ){
    new TH1F(Form("KPkm_MM_%d",i), "d(K^{-}, p K^{-})\"X\"", 1000, 0.5, 1.5);
  }

  new TH1F("hitpatFC_mmN", "hitpattern FC", 61, 0.5, 61.5);
  new TH1F("hitpatFC_mmP", "hitpattern FC", 61, 0.5, 61.5);
  new TH1F("hitpatFC_gamma", "hitpattern FC", 61, 0.5, 61.5);

  new TH1F("FC_offset_mmN",  "FC_offset_mmN", 500, -10, 10);
  new TH1F("FC_offset_mmP",  "FC_offset_mmP", 500, -10, 10);
  new TH1F("FC_offset_mmL",  "FC_offset_mmN", 500, -10, 10);
  new TH1F("FC_offset_mmS0", "FC_offset_mmN", 500, -10, 10);

  for( int seg=1; seg<=34; seg++ ){
    new TH1F(Form("CVC%d_offset_gamma", seg), Form("CVC seg%d offset by #gamma", seg), 500, -10.0, 10.0);
    new TH1F(Form("CVC%d_offset_mmN",  seg), Form("CVC seg%d offset by \"n\"", seg), 500, -10.0, 10.0);
    new TH1F(Form("CVC%d_offset_mmP",  seg), Form("CVC seg%d offset by \"p\"", seg), 500, -10.0, 10.0);
    new TH1F(Form("CVC%d_offset_mmL",  seg), Form("CVC seg%d offset by \"p\"", seg), 500, -10.0, 10.0);
    new TH1F(Form("CVC%d_offset_mmS0", seg), Form("CVC seg%d offset by \"p\"", seg), 500, -10.0, 10.0);

    new TH2F(Form("CVC%d_eu_offset_gamma", seg), Form("CVC seg%d dE_{up} vs offset",seg), 500, 0, 100, 500, -10, 10);
    new TH2F(Form("CVC%d_ed_offset_gamma", seg), Form("CVC seg%d dE_{down} vs offset",seg), 500, 0, 100, 500, -10, 10);

    new TNtuple(Form("CVC%d_slewing_info_gamma", seg), Form("CVC seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");
  }
  for( int seg=1; seg<=27; seg++ ){
    new TH1F(Form("PC%d_offset_gamma", seg), Form("PC seg%d offset by #gamma", seg), 500, -10.0, 10.0);
    new TH1F(Form("PC%d_offset_mmN",  seg), Form("PC seg%d offset by \"n\"", seg), 500, -10.0, 10.0);
    new TH1F(Form("PC%d_offset_mmP",  seg), Form("PC seg%d offset by \"p\"", seg), 500, -10.0, 10.0);
    new TH1F(Form("PC%d_offset_mmL",  seg), Form("PC seg%d offset by \"p\"", seg), 500, -10.0, 10.0);
    new TH1F(Form("PC%d_offset_mmS0", seg), Form("PC seg%d offset by \"p\"", seg), 500, -10.0, 10.0);

    new TH2F(Form("PC%d_eu_offset_gamma", seg), Form("PC seg%d dE_{up} vs offset",seg), 500, 0, 100, 500, -10, 10);
    new TH2F(Form("PC%d_ed_offset_gamma", seg), Form("PC seg%d dE_{down} vs offset",seg), 500, 0, 100, 500, -10, 10);

    new TNtuple(Form("PC%d_slewing_info_gamma", seg), Form("PC seg%d slewing info", seg), "eu:ed:tu:td:ctu:ctd:ctm");
  }
}

void fillFC(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo, 
	    ForwardChargeInfo *fcInfo)
{
  if( header ){
    if( !header->IsTrig(Trig_Charged) ) return;
    if( !header->IsTrig(Trig_Kaon) ) return;
  }
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag () ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
  //  if( anaInfo->nFCharge()!=1 ) return;
  //  ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
  if( fcInfo->step()<2 ) return;

  BeamInfo *beam=anaInfo->beam(0);
  HodoscopeLikeHit *fc_hit=fcInfo->hodo(blMan);
  //  cout<<beam->pid()<<"  "<<beam->mass()<<endl;

  //  cout<<"mom : "<<fcInfo->momByAng()<<"  mass2 : "<<fcInfo->mass2byAng()<<endl;
  MyHistTools::fillTH("FC_overbeta", 1./fcInfo->beta());
  MyHistTools::fillTH("FC_overbeta_dE", 1./fcInfo->beta(), fc_hit->emean());
  MyHistTools::fillTH("FC_mass2_mom_ang", fcInfo->mass2byAng(), fcInfo->momByAng());
  MyHistTools::fillTH("FC_mass2_mom_RK", fcInfo->mass2byRK(), fcInfo->momByRK());

  if( FC_P_MIN<fcInfo->mass2byRK() && fcInfo->mass2byRK()<FC_P_MAX ){
    MyHistTools::fillTH("FC_mom_RK_TOF_p", fcInfo->momByRK(), fcInfo->momByTOF());
    MyHistTools::fillTH("FC_mom_ang_TOF_p", fcInfo->momByAng(), fcInfo->momByTOF());

    double mm=(beam->lmom()+D_LMOM-fcInfo->lmom()).M();
    MyHistTools::fillTH("KP_MM", mm);
    if( anaInfo->nCDS(CDS_PiMinus)>0 ) MyHistTools::fillTH("KP_MM_pim", mm);
    if( anaInfo->nCDS(CDS_PiPlus)>0  ) MyHistTools::fillTH("KP_MM_pip", mm);
    if( anaInfo->nCDS(CDS_Kaon)>0    ) MyHistTools::fillTH("KP_MM_km", mm);
    if( anaInfo->nCDS(CDS_Proton)>0  ) MyHistTools::fillTH("KP_MM_p", mm);

    if( anaInfo->nCDS(CDS_Kaon)==1 ){
      CDSInfo *km=anaInfo->CDS(CDS_Kaon, 0);
      double im=(km->lmom()+fcInfo->lmom()).M();
      double mm_km=(beam->lmom()+D_LMOM-fcInfo->lmom()-km->lmom()).M();

      MyHistTools::fillTH("KPkm_MM_Pkm_IM", mm_km, im);

      int seg=-1;
      if( fc_hit->cid()==CID_PC ) seg=34+fc_hit->seg();
      else seg=fc_hit->seg();
      MyHistTools::fillTH(Form("KPkm_MM_%d",seg), mm_km);

      if( 0.9<mm_km && mm_km<0.98 ){
	int seg=-1;
	if( fc_hit->cid()==CID_PC ) seg=34+fc_hit->seg();
	else seg=fc_hit->seg();
	MyHistTools::fillTH("hitpatFC_mmN", seg);
	double offset=calc_offset(nMass, beam->lmom()+D_LMOM-km->lmom(), fcInfo->lmom(), anaInfo);
	if( fc_hit->cid()==CID_PC ) MyHistTools::fillTH(Form("PC%d_offset_mmN", fc_hit->seg()), offset);
	else MyHistTools::fillTH(Form("CVC%d_offset_mmN", fc_hit->seg()), offset);
	MyHistTools::fillTH("FC_offset_mmN", offset);
      }
    }

    for( int i=0; i<anaInfo->nCDS(CDS_PiMinus); i++ ){
      CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, i);
      double im=(fcInfo->lmom()+pim->lmom()).M();
      double mm_pim=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim->lmom()).M();
      MyHistTools::fillTH("Ppim_IM", im);
      MyHistTools::fillTH("KPpim_MM_Ppim_IM", mm_pim, im);
    }
       // cout<<"FP"<<endl;
       // anaInfo->dump();

    if( anaInfo->nCDS(CDS_PiMinus)==2 ){
      //           cout<<" p pi- pi- event"<<endl;
      CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
      CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
      if( pim0->dca()>pim1->dca() ) swap(pim0, pim1);

      double mm_pim0=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim0->lmom()).M();
      double mm_pim1=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim1->lmom()).M();
      bool mmL_flag=false;
      bool mmS0_flag=false;
      if( FC_mmL_MIN<mm_pim0 && mm_pim0<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmL_MIN<mm_pim1 && mm_pim1<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmS0_MIN<mm_pim0 && mm_pim0<FC_mmS0_MAX ) mmS0_flag=true;
      if( FC_mmS0_MIN<mm_pim1 && mm_pim1<FC_mmS0_MAX ) mmS0_flag=true;

      CDS2Info *pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiMinus, 0);
      double mm_pimpim=(beam->lmom()+D_LMOM-fcInfo->lmom()-pipi->lmom()).M();

      MyHistTools::fillTH("KPpim_KPpim_MM", mm_pim0, mm_pim1);
      MyHistTools::fillTH("KP_MM_2pim", mm);
      MyHistTools::fillTH("KPpimpim_MM", mm_pimpim);
      if( mmL_flag ) MyHistTools::fillTH("KPpimpim_MM_mmL", mm_pimpim);
      if( mmS0_flag ) MyHistTools::fillTH("KPpimpim_MM_mmS0", mm_pimpim);

      if( FC_mmP_MIN<mm_pimpim && mm_pimpim<FC_mmP_MAX ){
	int seg=-1;
	if( fc_hit->cid()==CID_PC ) seg=34+fc_hit->seg();
	else seg=fc_hit->seg();
	MyHistTools::fillTH("hitpatFC_mmP", seg);
	MyHistTools::fillTH("KPpim_KPpim_MM_mmP", mm_pim0, mm_pim1);
	MyHistTools::fillTH("KPpim_MM_mmP", mm_pim0);
	MyHistTools::fillTH("KPpim_MM_mmP", mm_pim1);
	MyHistTools::fillTH("KP_MM_2pim_mmP", mm);
	MyHistTools::fillTH("KP_MM_2pim_mmP_wL", mm);

	double offset=calc_offset(pMass, beam->lmom()+D_LMOM-pipi->lmom(), fcInfo->lmom(), anaInfo);
	MyHistTools::fillTH("FC_offset_mmP", offset);
	if( fc_hit->cid()==CID_PC ) MyHistTools::fillTH(Form("PC%d_offset_mmP", fc_hit->seg()), offset);
	else MyHistTools::fillTH(Form("CVC%d_offset_mmP", fc_hit->seg()), offset);

	double diff0=fabs(mm_pim0-lMass);
	double diff1=fabs(mm_pim1-lMass);
	double offset2=-999;
	if( diff0<diff1 && diff0<0.035 ) offset2=calc_offset(lMass, beam->lmom()+D_LMOM-pim0->lmom(), fcInfo->lmom(), anaInfo);
	if( diff1<diff1 && diff1<0.035 ) offset2=calc_offset(lMass, beam->lmom()+D_LMOM-pim1->lmom(), fcInfo->lmom(), anaInfo);
	if( offset2!=-999 ) MyHistTools::fillTH("FC_offset_mmL", offset2);
      }
      if( 0.99<mm_pimpim && mm_pimpim<1.06 ){
	MyHistTools::fillTH("KPpim_KPpim_MM_mmPgamma", mm_pim0, mm_pim1);
	MyHistTools::fillTH("KPpim_MM_mmPgamma", mm_pim0);
	MyHistTools::fillTH("KPpim_MM_mmPgamma", mm_pim1);
	MyHistTools::fillTH("KP_MM_2pim_mmPgamma", mm);
	MyHistTools::fillTH("KP_MM_2pim_mmPgamma_wS0", mm);

	double diff0=fabs(mm_pim0-s0Mass);
	double diff1=fabs(mm_pim1-s0Mass);
	double offset2=-999;
	if( diff0<diff1 && diff0<0.035 ) offset2=calc_offset(s0Mass, beam->lmom()+D_LMOM-pim0->lmom(), fcInfo->lmom(), anaInfo);
	if( diff1<diff1 && diff1<0.035 ) offset2=calc_offset(s0Mass, beam->lmom()+D_LMOM-pim1->lmom(), fcInfo->lmom(), anaInfo);
	if( offset2!=-999 ) MyHistTools::fillTH("FC_offset_mmS0", offset2);
      }
    }
    if( anaInfo->nCDS(CDS_Proton)==1 && anaInfo->nCDS(CDS_PiMinus)>0 ) MyHistTools::fillTH("KP_MM_ppim", mm);
  }
}

void fillFC_gamma(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo)
{
  //  if( !header->IsTrig(Trig_Charged) ) return;
  if( header && !header->IsTrig(Trig_Kaon) ) return;
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag () ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
  //  if( anaInfo->nFCharge()!=1 ) return;
  //  ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
  BeamInfo *beam=anaInfo->beam(0);
  CDSInfo *vtx=anaInfo->minDCA();

  vector<HodoscopeLikeHit*> BVC_hits=MyTools::getHodo(blMan, CID_BVC);
  vector<HodoscopeLikeHit*> CVC_hits=MyTools::getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> PC_hits=MyTools::getHodo(blMan, CID_PC);
  
  if( BVC_hits.size()==0 && beam->nFDC1()==0 ){
    HodoscopeLikeHit *fc_hit=0;
    if( CVC_hits.size()==1 && PC_hits.size()==0 ) fc_hit=CVC_hits[0];
    if( CVC_hits.size()==0 && PC_hits.size()==1 ) fc_hit=PC_hits[0];
    if( fc_hit ){
      int seg=-1;
      if( fc_hit->cid()==CID_PC ) seg=34+fc_hit->seg();
      else seg=fc_hit->seg();

      TVector3 pos;
      conf->GetGeomMapManager()-> GetGPos(fc_hit->cid(), fc_hit->seg(), pos);

      double fl=(pos-vtx->vertexBeam()).Mag();

      double mom_out, tmp_tof;
      ELossTools::CalcElossBeamTGeo(beam->T0pos(), vtx->vertexBeam(), beam->D5mom(), parMass[beam->pid()], mom_out, tmp_tof);
      double tof=fc_hit->ctmean()-beam->T0time()-tmp_tof;
      double calc_tof=fl/(100.*Const);

      double eu=fc_hit->eu(), ed=fc_hit->ed();
      double tu=fc_hit->tu()-beam->T0time()-tmp_tof-calc_tof;
      double td=fc_hit->td()-beam->T0time()-tmp_tof-calc_tof;
      double ctu=fc_hit->ctu()-beam->T0time()-tmp_tof-calc_tof;
      double ctd=fc_hit->ctd()-beam->T0time()-tmp_tof-calc_tof;
      double offset=fc_hit->ctmean()-beam->T0time()-tmp_tof-calc_tof;

      if( offset<-5 || 5<offset ) return;
      MyHistTools::fillTH("hitpatFC_gamma", seg);

      if( fc_hit->cid()==CID_PC ){
	MyHistTools::fillTH(Form("PC%d_offset_gamma", fc_hit->seg()), tof-calc_tof);
	MyHistTools::fillTH(Form("PC%d_eu_offset_gamma", fc_hit->seg()), fc_hit->eu(), tof-calc_tof);
	MyHistTools::fillTH(Form("PC%d_ed_offset_gamma", fc_hit->seg()), fc_hit->ed(), tof-calc_tof);

	TNtuple *tup=(TNtuple*)gFile->Get(Form("PC%d_slewing_info_gamma", fc_hit->seg()));
	tup->Fill(eu, ed, tu, td, ctu, ctd, offset);
      }
      else{
	MyHistTools::fillTH(Form("CVC%d_offset_gamma", fc_hit->seg()), tof-calc_tof);
	MyHistTools::fillTH(Form("CVC%d_eu_offset_gamma", fc_hit->seg()), fc_hit->eu(), tof-calc_tof);
	MyHistTools::fillTH(Form("CVC%d_ed_offset_gamma", fc_hit->seg()), fc_hit->ed(), tof-calc_tof);

	TNtuple *tup=(TNtuple*)gFile->Get(Form("CVC%d_slewing_info_gamma", fc_hit->seg()));
	tup->Fill(eu, ed, tu, td, ctu, ctd, offset);
      }
      //      cout<<" FC seg"<<seg<<"  tof="<<tof<<"  calc tof"<<calc_tof<<"  offset="<<tof-calc_tof<<endl;
    }
  }
}

double calc_offset(const double &mm, const TLorentzVector tot_lmom, const TLorentzVector lmom, AnaInfo *anaInfo)
{
  TLorentzVector fit_lmom=lmom;
  double fit_mm=(tot_lmom-lmom).M();
  double delta=fit_mm-mm;
  //  cout<<"Particle mass : "<<mm<<"  fit mass : "<<fit_mm<<"  delta : "<<delta<<endl;

  int count=0;
  while( true ){
    TVector3 mom=fit_lmom.Vect();
    mom.SetMag(mom.Mag()+delta);
    fit_lmom.SetVectM(mom, lmom.M());
    fit_mm=(tot_lmom-fit_lmom).M();
    delta=fit_mm-mm;

    // cout<<"Particle mass : "<<mm<<"  fit mass : "<<fit_mm<<"  delta : "<<delta<<endl;
    // string input;
    // cin>>input;
    // if( input=="q" ) exit(0);
    count++;
    if( count>100 ) return -99999;
    if( fabs(delta)<0.0001 ) break;
  }
  //  cout<<"mom : "<<lmom.Vect().Mag()<<"  fit mom : "<<fit_lmom.Vect().Mag()<<endl;

  TVector3 vertex=anaInfo->minDCA()->vertexBeam();
  TVector3 FDC1pos=anaInfo->forwardCharge(0)->FDC1().pos();
  double fl=anaInfo->forwardCharge(0)-> fl();
  double tmp_mom, tof1, tof2;
  ELossTools::CalcElossForwardTGeo(vertex, FDC1pos, fl, lmom.Vect().Mag(), lmom.M(), tmp_mom, tof1);
  ELossTools::CalcElossForwardTGeo(vertex, FDC1pos, fl, fit_lmom.Vect().Mag(), fit_lmom.M(), tmp_mom, tof2);

  double mom_out, tmp_tof;
  //  ELossTools::CalcElossBeamTGeo(anaInfo->beam(0)->T0pos(), vertex, anaInfo->beam(0)->D5mom(), parMass[anaInfo->beam(0)->pid()], mom_out, tmp_tof);
  //  double tof=anaInfo->forwardCharge(0)->time()-anaInfo->beam(0)->T0time()-tmp_tof;
  //  cout<<" TOF="<<tof<<"  tof1="<<tof1<<"  tof2="<<tof2<<endl;

  return tof1-tof2;
}

bool kd_pLpim(AnaInfo *anaInfo){
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag () ) return false;
  if( anaInfo->minDCA() ){
    cout<<GeomTools::GetID(anaInfo->minDCA()->vertexBeam())<<endl;
  }
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return false;
  if( anaInfo->nFCharge()!=1 ) return false;
  ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
  BeamInfo *beam=anaInfo->beam(0);
  if( fcInfo->step()<2 ) return false;

  if( FC_P_MIN<fcInfo->mass2byRK() && fcInfo->mass2byRK()<FC_P_MAX ){
    //    cout<<"FP event"<<endl;
    // TVector3 pos=beam->T0pos();
    // anaInfo-> dump();

    double mm=(beam->lmom()+D_LMOM-fcInfo->lmom()).M();
    if( anaInfo->nCDS(CDS_PiMinus)==2 ){
      cout<<"    CDS 2pi-"<<endl;
      CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
      CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
      if( pim0->dca()>pim1->dca() ) swap(pim0, pim1);

      double mm_pim0=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim0->lmom()).M();
      double mm_pim1=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim1->lmom()).M();
      bool mmL_flag=false;
      bool mmS0_flag=false;
      if( FC_mmL_MIN<mm_pim0 && mm_pim0<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmL_MIN<mm_pim1 && mm_pim1<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmS0_MIN<mm_pim0 && mm_pim0<FC_mmS0_MAX ) mmS0_flag=true;
      if( FC_mmS0_MIN<mm_pim1 && mm_pim1<FC_mmS0_MAX ) mmS0_flag=true;

      CDS2Info *pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiMinus, 0);
      double mm_pimpim=(beam->lmom()+D_LMOM-fcInfo->lmom()-pipi->lmom()).M();

      if( FC_mmP_MIN<mm_pimpim && mm_pimpim<FC_mmP_MAX && mmL_flag ) return true;
    }
  }
  return false;
}


bool kd_pS0pim(AnaInfo *anaInfo){
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag () ) return false;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return false;
  if( anaInfo->nFCharge()!=1 ) return false;
  ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
  BeamInfo *beam=anaInfo->beam(0);
  if( fcInfo->step()<2 ) return false;

  if( FC_P_MIN<fcInfo->mass2byRK() && fcInfo->mass2byRK()<FC_P_MAX ){
    double mm=(beam->lmom()+D_LMOM-fcInfo->lmom()).M();
    if( anaInfo->nCDS(CDS_PiMinus)==2 ){
      CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
      CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
      if( pim0->dca()>pim1->dca() ) swap(pim0, pim1);

      double mm_pim0=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim0->lmom()).M();
      double mm_pim1=(beam->lmom()+D_LMOM-fcInfo->lmom()-pim1->lmom()).M();
      bool mmL_flag=false;
      bool mmS0_flag=false;
      if( FC_mmL_MIN<mm_pim0 && mm_pim0<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmL_MIN<mm_pim1 && mm_pim1<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmS0_MIN<mm_pim0 && mm_pim0<FC_mmS0_MAX ) mmS0_flag=true;
      if( FC_mmS0_MIN<mm_pim1 && mm_pim1<FC_mmS0_MAX ) mmS0_flag=true;

      CDS2Info *pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiMinus, 0);
      double mm_pimpim=(beam->lmom()+D_LMOM-fcInfo->lmom()-pipi->lmom()).M();

      if( 0.99<mm_pimpim && mm_pimpim<1.06 && mmS0_flag ) return true;
    }
  }
  return false;
}

