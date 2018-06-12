#include "MyHistReadKNpipi.h"

using namespace std;

const int binWidth=5;

void initHistReadKNpipi()
{
  new TH1F("KN_MM_pipi_woAll", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_KNpipi_MM_woAll", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi", "d(K^{-}, n)\"X\"",          2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN", "d(K^{-}, n)\"X\"",       2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_wK0", "d(K^{-}, n)\"X\"",   2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_wSm", "d(K^{-}, n)\"X\"",   2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_wSp", "d(K^{-}, n)\"X\"",   2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woAll", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);

  new TH1F("KNpipi_MM",      "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wK0",  "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wSm",  "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wSp",  "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_true", "d(K^{-}, n #pi^{+} #pi^{-}", 2000, 0.0, 2.0);

  new TH1F("CDS_IM_pipi_wN", "CDS IM #pi^{+} #pi^{-}", 1000, 0.0, 1.0);
  new TH1F("Npim_IM_pip", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim", "n #pi^{+} #pi^{-}", 1000, 1.0, 2.0);

  new TH1F("CDS_IM_pipi_wN_true", "CDS IM #pi^{+} #pi^{-}", 1000, 0.0, 1.0);
  new TH1F("Npim_IM_pip_true", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_true", "n #pi^{+} #pi^{-}", 1000, 1.0, 2.0);

  new TH1F("CDS_IM_pipi_wN_wN", "CDS IM #pi^{+} #pi^{-}", 1000, 0.0, 1.0);
  new TH2F("CDS_IM_pipi_wN_wN_pipi_mom", "CDS IM #pi^{+} #pi^{-} vs pi^{+} #pi^{-} mom", 1000, 0.0, 1.0, 1000, 0.0, 1.0);
  new TH1F("Npim_IM_pip_wN", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN", "n #pi^{+} #pi^{-}", 1000, 1.0, 2.0);

  new TH1F("CDS_IM_pipi_wN_wN_true", "CDS IM #pi^{+} #pi^{-}", 1000, 0.0, 1.0);
  new TH1F("Npim_IM_pip_wN_true", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN_true", "n #pi^{+} #pi^{-}", 1000, 1.0, 2.0);

  new TH2F("KNpim_KNpip_MM",       "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpim_KNpip_MM_wN",       "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpim_KNpip_MM_wN_wK0",   "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpim_KNpip_MM_wN_wSm",   "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpim_KNpip_MM_wN_wSp",   "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpim_KNpip_MM_wN_woAll", "d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  for( int mm=1350; mm<1900; mm+=binWidth ){
    new TH2F(Form("KNpim_KNpip_MM_wN_woAll_%d_%d", mm, mm+binWidth),  
	     Form("d(K^{-}, n #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\" : d(K^{-}, n)\"X\" %d~%d [GeV/c^{2}]", mm, mm+binWidth),
	     1000, 1.0, 2.0, 1000, 1.0, 2.0);
  }

  new TH1F("Npipi_IM_wN",     "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSm", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSp", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);

  //*** for fitting K- d-> n pi- pi+ "n" ***//
  new TH2F("KN_MM_pipi_wN_woAll_KNpim_MM", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_wN_woAll_KNpip_MM", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_wN_CDS_IM_pipi", "d(K^{-}, n)\"X\" vs CDS IM #pi^{+} #pi^{-}", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_wN_Npim_IM",     "d(K^{-}, n)\"X\" vs n #pi^{-} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_pipi_wN_Npip_IM",     "d(K^{-}, n)\"X\" vs n #pi^{+} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);

  new TH1F("mmN_mom_wK0", "d(K^{-}, n #pi^{-} #pi^{+})\"n\" mom", 1000, 0.0, 1.0);
  new TH1F("mmN_mom_wSm", "d(K^{-}, n #pi^{-} #pi^{+})\"n\" mom", 1000, 0.0, 1.0);
  new TH1F("mmN_mom_wSp", "d(K^{-}, n #pi^{-} #pi^{+})\"n\" mom", 1000, 0.0, 1.0);

  new TH2F("KN_MM_pipi_wN_wK0_mmN_mom",   "d(K^{-}, n)\"X\" vs \"n\" mom",   2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_wN_wSm_mmN_mom",   "d(K^{-}, n)\"X\" vs \"n\" mom",   2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_wN_wSp_mmN_mom",   "d(K^{-}, n)\"X\" vs \"n\" mom",   2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_wN_woAll_mmN_mom", "d(K^{-}, n)\"X\" vs \"n\" mom",   2000, 0.0, 2.0, 1000, 0.0, 1.0);

  new TH2F("KN_MM_pipi_wN_wK0_Npipi_IM",   "d(K^{-}, n)\"X\" vs n #pi^{+} #pi^{-} IM",   2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_pipi_wN_wSm_Npipi_IM",   "d(K^{-}, n)\"X\" vs n #pi^{+} #pi^{-} IM",   2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_pipi_wN_wSp_Npipi_IM",   "d(K^{-}, n)\"X\" vs n #pi^{+} #pi^{-} IM",   2000, 0.0, 2.0, 1000, 1.0, 2.0);

  new TH1F("npipin_K0_mom", "K^{0} mom", 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_wN_wK0_K0_mom",   "d(K^{-}, n)\"X\" vs K^{0} mom",   2000, 0.0, 2.0, 1000, 0.0, 1.0);
}

void fillHistReadKNpipi(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  if( anaInfo->nFNeutral()!=1 ) return;
  if( header ){
    if( !header->trigmode2(Mode_KCDH1N) ) return;
  }
  vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  if( CVChits.size()!=0 || BVChits.size()!=0 ) return;
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;

  ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
  if( fnInfo->pid()!=F_Neutron ) return;

  if( anaInfo->minDCA() && GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
  if( cdstrackMan->nGoodTrack()!=2 ) return;
  if( anaInfo->nCDS(CDS_PiMinus)!=1 || anaInfo->nCDS(CDS_PiPlus)!=1 ) return;
  if( !anaInfo->CDS(CDS_PiMinus, 0)->flag() || !anaInfo->CDS(CDS_PiPlus)->flag() ) return;
  if( !anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->flag() ){ cout<<"  !!! CDS2 pi+ pi- flag is false !!!"<<endl; return; }

  TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
  TLorentzVector target_lmom=MyAnaTools::target_lmom();
  TLorentzVector fn_lmom=fnInfo->lmom();
  TLorentzVector kn_lmom=beam_lmom+target_lmom-fn_lmom;
  CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, 0); TLorentzVector pim_lmom=pim->lmom();
  CDSInfo *pip=anaInfo->CDS(CDS_PiPlus, 0); TLorentzVector pip_lmom=pip->lmom();
  TLorentzVector knpim_lmom=kn_lmom-pim_lmom;
  TLorentzVector knpip_lmom=kn_lmom-pip_lmom;
  double npipi_im=(fn_lmom+pim_lmom+pip_lmom).M();
  TVector3 mmN_mom=(kn_lmom-pim_lmom-pip_lmom).Vect();

  if( pim->CDHseg()==pip->CDHseg() ) return;

  CDS2Info *pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0);

  double kn_pipi_mm=(kn_lmom-pim_lmom-pip_lmom).M();
  MyHistTools::fillTH("KN_MM_pipi", kn_lmom.M());
  MyHistTools::fillTH("KNpipi_MM", kn_pipi_mm);
  TTree *tree= (TTree*)gFile->Get("npipi_event");
  tree-> Fill();

  MyHistTools::fillTH("CDS_IM_pipi_wN", pipi->im());
  MyHistTools::fillTH("Npim_IM_pip", (fn_lmom+pim_lmom).M());
  MyHistTools::fillTH("Npip_IM_pim", (fn_lmom+pip_lmom).M());
  bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
  if( Npipi_K0_MIN<pipi->im() && pipi->im()<Npipi_K0_MAX ){
    MyHistTools::fillTH("KNpipi_MM_wK0", kn_pipi_mm);
    MyHistTools::fillTH("CDS_IM_pipi_wN_true", pipi->im());
    K0_flag=true;
  }
  if( Npipi_Sm_MIN<(fn_lmom+pim_lmom).M() && (fn_lmom+pim_lmom).M()<Npipi_Sm_MAX ){
    MyHistTools::fillTH("KNpipi_MM_wSm", kn_pipi_mm);
    MyHistTools::fillTH("Npim_IM_pip_true", (fn_lmom+pim_lmom).M());
    Sm_flag=true;
  }
  if( Npipi_Sp_MIN<(fn_lmom+pip_lmom).M() && (fn_lmom+pip_lmom).M()<Npipi_Sp_MAX ){
    MyHistTools::fillTH("KNpipi_MM_wSp", kn_pipi_mm);
    MyHistTools::fillTH("Npip_IM_pim_true", (fn_lmom+pip_lmom).M());
    Sp_flag=true;
  }
  MyHistTools::fillTH("KNpim_KNpip_MM", knpim_lmom.M(), knpip_lmom.M());

  if( !K0_flag && !Sm_flag && !Sp_flag ){
    MyHistTools::fillTH("KN_MM_pipi_woAll", kn_lmom.M());
    MyHistTools::fillTH("KN_MM_pipi_KNpipi_MM_woAll", kn_lmom.M(), kn_pipi_mm);
  }

  if( Npipi_N_MIN<kn_pipi_mm && kn_pipi_mm<Npipi_N_MAX ){
    MyHistTools::fillTH("Npipi_IM_wN", npipi_im);

    MyHistTools::fillTH("KNpipi_MM_true", kn_pipi_mm);
    MyHistTools::fillTH("KN_MM_pipi_wN", kn_lmom.M());

    MyHistTools::fillTH("CDS_IM_pipi_wN_wN", pipi->im());
    MyHistTools::fillTH("CDS_IM_pipi_wN_wN_pipi_mom", pipi->im(), pipi->lmom().Vect().Mag());
    MyHistTools::fillTH("Npim_IM_pip_wN", (fn_lmom+pim_lmom).M());
    MyHistTools::fillTH("Npip_IM_pim_wN", (fn_lmom+pip_lmom).M());
    MyHistTools::fillTH("KNpim_KNpip_MM_wN", knpim_lmom.M(), knpip_lmom.M());

    MyHistTools::fillTH("KN_MM_pipi_wN_CDS_IM_pipi", kn_lmom.M(), pipi->im());
    MyHistTools::fillTH("KN_MM_pipi_wN_Npim_IM", kn_lmom.M(), (fn_lmom+pim_lmom).M());
    MyHistTools::fillTH("KN_MM_pipi_wN_Npip_IM", kn_lmom.M(), (fn_lmom+pip_lmom).M());

    if( Npipi_K0_MIN<pipi->im() && pipi->im()<Npipi_K0_MAX ){
      MyHistTools::fillTH("CDS_IM_pipi_wN_wN_true", pipi->im());
      MyHistTools::fillTH("KN_MM_pipi_wN_wK0", kn_lmom.M());
      MyHistTools::fillTH("KNpim_KNpip_MM_wN_wK0", knpim_lmom.M(), knpip_lmom.M());

      MyHistTools::fillTH("Npipi_IM_wN_wK0", npipi_im);

      MyHistTools::fillTH("mmN_mom_wK0", mmN_mom.Mag());
      MyHistTools::fillTH("KN_MM_pipi_wN_wK0_mmN_mom", kn_lmom.M(), mmN_mom.Mag());

      MyHistTools::fillTH("npipin_K0_mom", pipi->lmom.Vect().Mag());
      MyHistTools::fillTH("KN_MM_pipi_wN_wK0_K0_mom", kn_lmom.M(), pipi->lmom().Vect().Mag());

      MyHistTools::fillTH("KN_MM_pipi_wN_wK0_Npipi_IM", kn_lmom.M(), npipi_im);
    }
    if( Npipi_Sm_MIN<(fn_lmom+pim_lmom).M() && (fn_lmom+pim_lmom).M()<Npipi_Sm_MAX ){
      MyHistTools::fillTH("Npim_IM_pip_wN_true", (fn_lmom+pim_lmom).M());
      MyHistTools::fillTH("KN_MM_pipi_wN_wSm", kn_lmom.M());
      MyHistTools::fillTH("KNpim_KNpip_MM_wN_wSm", knpim_lmom.M(), knpip_lmom.M());

      MyHistTools::fillTH("Npipi_IM_wN_wSm", npipi_im);

      MyHistTools::fillTH("mmN_mom_wSm", mmN_mom.Mag());
      MyHistTools::fillTH("KN_MM_pipi_wN_wSm_mmN_mom", kn_lmom.M(), mmN_mom.Mag());

      MyHistTools::fillTH("KN_MM_pipi_wN_wSm_Npipi_IM", kn_lmom.M(), npipi_im);
    }
    if( Npipi_Sp_MIN<(fn_lmom+pip_lmom).M() && (fn_lmom+pip_lmom).M()<Npipi_Sp_MAX ){
      MyHistTools::fillTH("Npip_IM_pim_wN_true", (fn_lmom+pip_lmom).M());
      MyHistTools::fillTH("KN_MM_pipi_wN_wSp", kn_lmom.M());
      MyHistTools::fillTH("KNpim_KNpip_MM_wN_wSp", knpim_lmom.M(), knpip_lmom.M());

      MyHistTools::fillTH("Npipi_IM_wN_wSp", npipi_im);

      MyHistTools::fillTH("mmN_mom_wSp", mmN_mom.Mag());
      MyHistTools::fillTH("KN_MM_pipi_wN_wSp_mmN_mom", kn_lmom.M(), mmN_mom.Mag());

      MyHistTools::fillTH("KN_MM_pipi_wN_wSp_Npipi_IM", kn_lmom.M(), npipi_im);
    }

    if( !K0_flag && !Sm_flag && !Sp_flag ){
      for( int mm=1350; mm<1900; mm+=binWidth){
	if( mm<1000*kn_lmom.M() && 1000*kn_lmom.M()<mm+binWidth ){
	  MyHistTools::fillTH(Form("KNpim_KNpip_MM_wN_woAll_%d_%d", mm, mm+binWidth), kn_lmom.M());
	  break;
	}
      }

      MyHistTools::fillTH("KN_MM_pipi_wN_woAll", kn_lmom.M());
      MyHistTools::fillTH("KNpim_KNpip_MM_wN_woAll", knpim_lmom.M(), knpip_lmom.M());

      //*** for fit ***//
      MyHistTools::fillTH("KN_MM_pipi_wN_woAll_KNpim_MM", kn_lmom.M(), knpim_lmom.M());
      MyHistTools::fillTH("KN_MM_pipi_wN_woAll_KNpip_MM", kn_lmom.M(), knpip_lmom.M());

      MyHistTools::fillTH("KN_MM_pipi_wN_woAll_mmN_mom", kn_lmom.M(), mmN_mom.Mag());
    }
  }
}
