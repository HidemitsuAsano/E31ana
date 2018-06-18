#include "MyHistReadKNpipi_H2.h"

using namespace std;

void initHistReadKNpipi_H2()
{
  //  cout<<" Initialize Histogram for H2 target"<<endl;
  new TH1F("Npipi_MM", "p(K^{-}, #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("Npipi_MM_true", "p(K^{-}, #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);

  new TH1F("CDS_IM_pipi_mmN", "CDS #pi^{-} #pi^{+} IM", 1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_mmN_true", "CDS #pi^{-} #pi^{+} IM", 1000, 0.0, 1.0);

  new TH1F("CDS_IM_pipi_mmN_woS", "CDS #pi^{-} #pi^{+} IM", 1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_mmN_woS_true", "CDS #pi^{-} #pi^{+} IM", 1000, 0.0, 1.0);

  new TH1F("Npim_MM_mmN", "p(K^{-}, #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("Npip_MM_mmN", "p(K^{-}, #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("Npim_MM_mmN_true", "p(K^{-}, #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("Npip_MM_mmN_true", "p(K^{-}, #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH1F("Npim_MM_mmN_woK0", "p(K^{-}, #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("Npip_MM_mmN_woK0", "p(K^{-}, #pi^{+})\"X\"", 2000, 0.0, 2.0);

  new TH2F("Npim_Npip_MM_mmN", "", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH2F("NC_calc_hitpos",      "", 10000, -500, 500, 10000, -500, 500);
  new TH2F("NC_calc_hitpos_whit", "", 10000, -500, 500, 10000, -500, 500);

  new TH1F("nNC_mmN", "", 113, -0.5, 112.5); 

  new TH1F("effNC" ,"", 10, 0, 10);
  new TH1F("trigNC" ,"", 10, 0, 10);
}

void fillHistReadKNpipi_H2(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;
  if( cdstrackMan->nGoodTrack()!=2 ) return;

  if( anaInfo->nCDS(CDS_PiMinus)!=1 || anaInfo->nCDS(CDS_PiPlus)!=1 ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;

  if( !anaInfo->CDS(CDS_PiMinus, 0)->flag() || !anaInfo->CDS(CDS_PiPlus)->flag() ) return;
  if( !anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->flag() ){ cout<<"  !!! CDS2 pi+ pi- flag is false !!!"<<endl; return; }

  if( header->trigmode2(Mode_KCDH2, conf) ){
    CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, 0);
    CDSInfo *pip=anaInfo->CDS(CDS_PiPlus, 0);
    CDS2Info *pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0);
    std::vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
    std::vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
    if( BVChits.size()!=0 ) return;
    if( CVChits.size()!=0 ) return;

    TLorentzVector target_lmom=MyAnaTools::target_lmom();
    TLorentzVector pipi_lmom=pipi->lmom();
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TLorentzVector pim_lmom=pim->lmom();
    TLorentzVector pip_lmom=pip->lmom();

    double kpipi_mm=(beam_lmom+target_lmom-pipi_lmom).M();
    double kpim_mm=(beam_lmom+target_lmom-pim_lmom).M();
    double kpip_mm=(beam_lmom+target_lmom-pip_lmom).M();

    MyHistTools::fillTH("Npipi_MM", kpipi_mm);

    if( 0.9<kpipi_mm && kpipi_mm<0.98 ){
      bool K0_flag=false;
      bool Sm_flag=false;
      bool Sp_flag=false;
      if( H2_K0_MIN<pipi->im() && pipi->im()<H2_K0_MAX ) K0_flag=true;
      if( H2_Sm_MIN<kpip_mm && kpip_mm<H2_Sm_MAX ) Sm_flag=true;
      if( H2_Sp_MIN<kpim_mm && kpim_mm<H2_Sp_MAX ) Sp_flag=true;

      MyHistTools::fillTH("Npipi_MM_true", kpipi_mm);

      MyHistTools::fillTH("CDS_IM_pipi_mmN", pipi->im());
      if( K0_flag )  MyHistTools::fillTH("CDS_IM_pipi_mmN_true", pipi->im());

      if( !Sm_flag && !Sp_flag ){
	MyHistTools::fillTH("CDS_IM_pipi_mmN_woS", pipi->im());
	if( K0_flag )  MyHistTools::fillTH("CDS_IM_pipi_mmN_woS_true", pipi->im());
      }

      MyHistTools::fillTH("Npim_MM_mmN", kpim_mm);
      MyHistTools::fillTH("Npip_MM_mmN", kpip_mm);
      if( Sp_flag ) MyHistTools::fillTH("Npim_MM_mmN_true", kpim_mm);
      if( Sm_flag ) MyHistTools::fillTH("Npip_MM_mmN_true", kpip_mm);
      MyHistTools::fillTH("Npim_Npip_MM_mmN", kpim_mm, kpip_mm);

      if( !K0_flag ){
	MyHistTools::fillTH("Npim_MM_mmN_woK0", kpim_mm);
	MyHistTools::fillTH("Npip_MM_mmN_woK0", kpip_mm);
      }

      if( !Sm_flag && !Sp_flag && K0_flag ){
	TLorentzVector n_lmom=beam_lmom+target_lmom-pipi_lmom;
	TVector3 n_mom=n_lmom.Vect();

	TVector3 nc_pos;
	conf->GetGeomMapManager()->GetGPos(CID_NC, 1, nc_pos);
	TVector3 vertex=pipi->vertexBeam();
	TVector3 n_dir=n_mom.Unit();
	TVector3 calc_nc_pos=vertex+((nc_pos.Z()-vertex.Z())/n_dir.Z())*n_dir;

	vector<HodoscopeLikeHit*> NC_hits=MyTools::getHodo(blMan, CID_NC);
	vector<vector<HodoscopeLikeHit*> > NC_layer_hit=MyTools::getNChits(blMan);

	MyHistTools::fillTH("NC_calc_hitpos", calc_nc_pos.X(), calc_nc_pos.Y());
	MyHistTools::fillTH("nNC_mmN", NC_hits.size());


	//	if( calc_nc_pos.X()


	if( anaInfo->nFNeutral()==1 ){
	  MyHistTools::fillTH("NC_calc_hitpos_whit", calc_nc_pos.X(), calc_nc_pos.Y());
	}
      }
    }
  }
}
