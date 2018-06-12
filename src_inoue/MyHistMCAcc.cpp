#include "MyHistMCAcc.h"

using namespace std;

void initHistMCAcc()
{
  new TH1F("trigY_mass_fn",               "Y^{*} trigger event",   1000, 1.0, 2.0);
  new TH1F("effY_mass_fn_pipi_MC",        "Y^{*} effective event", 1000, 1.0, 2.0);
  new TH1F("effY_mass_fn_pipi_data",      "Y^{*} effective event", 1000, 1.0, 2.0);
  new TH1F("effY_mass_fn_pipi_mmN",       "Y^{*} effective event", 1000, 1.0, 2.0);
  new TH1F("effY_mass_fn_pipi_mmN_woAll", "Y^{*} effective event", 1000, 1.0, 2.0);

  new TNtuple("tup_trig_mass_fn",             "", "m");
  new TNtuple("tup_eff_mass_fn_pipi_MC",      "", "m");
  new TNtuple("tup_eff_mass_fn_pipi_data",    "", "m");
  new TNtuple("tup_eff_mass_fn_pipi_mmN",     "", "m");
  new TNtuple("tup_eff_mass_fn_pipi_mmN_woK0",  "", "m");
  new TNtuple("tup_eff_mass_fn_pipi_mmN_woAll", "", "m");

  new TH2F("fp_gen", "fp_gen", 1000, -0.5, 0.5, 1000, -0.1, 0.1);
  new TH2F("fp_acc", "fp_acc", 1000, -0.5, 0.5, 1000, -0.1, 0.1);

  new TNtuple("tup_acc_fp", "", "m:SA_x:SA_y:flag");
  new TNtuple("tup_trig_mass_fp",             "", "m");
  new TNtuple("tup_eff_mass_fp_pipi_MC",      "", "m");
  new TNtuple("tup_eff_mass_fp_pipi_data",    "", "m");
  new TNtuple("tup_eff_mass_fp_pipi_mm1",     "", "m");
  new TNtuple("tup_eff_mass_fp_pipi_mm2",     "", "m");

  new TH2F("KN_MM_diff_MC", "", 1000, 1.0, 2.0, 1000, -0.1, 0.1);
}

void fillHistMCAcc_fp(DetectorData *detData, MCData* mcData, ReactionData *reacData, AnaInfo *anaInfo,
		      BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan)
{
  TNtuple *tup;
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;

  if( reacData->ReactionID()==3098 ){
    Track *init_p=MyMCTools::initTrack(mcData, 2212);
    if( !init_p ) return;
    TVector3 mom=init_p->momentum();
    TVector3 dir=mom.Unit();
    MyHistTools::fillTH("fp_gen", asin(dir.X()), asin(dir.Y()));
    int flag=0;
    if( MyMCTools::effFC(detData, mcData, blMan) ) flag=1;
    tup=(TNtuple*)gFile->Get("tup_acc_fp"), tup->Fill(MyMCTools::Ystar_mass(reacData), asin(dir.X()), asin(dir.Y()), flag);

    if( MyMCTools::effFC(detData, mcData, blMan) ){
      MyHistTools::fillTH("fp_acc", asin(dir.X()), asin(dir.Y()));
      int mode=-1; // 1 : L pi-  2 : S0 pi-
      if( MyMCTools::trackByPDG(mcData, 3122) ) mode=1;
      if( MyMCTools::trackByPDG(mcData, 3212) ) mode=2;

      if( mode!=1 && mode!=2 ) cout<<"  !!! K-d -> S(1385)- p decay mode not found !!!"<<endl;
      //      cout<<"mode : "<<mode<<endl;
      //      cout<<MyMCTools::Ystar_mass(reacData)<<endl;
      tup=(TNtuple*)gFile->Get("tup_trig_mass_fp"), tup->Fill(MyMCTools::Ystar_mass(reacData));

      list<int> pim_track_id;
      for( int i=0; i<detData->detectorHitSize(); i++ ){
	DetectorHit *hit=detData->detectorHit(i);
	if( hit->pdg()==-211 && hit->detectorID()==CID_CDH ) pim_track_id.push_back(hit->trackID());
      }
      pim_track_id.sort(); pim_track_id.unique();

      if( pim_track_id.size()>2 ){
	//	cout<<"!!! n hit pi- : "<<pim_track_id.size()<<endl;
	for( int i=0; i<detData->detectorHitSize(); i++ ){
	  DetectorHit *hit=detData->detectorHit(i);
	  if( hit->pdg()==-211 && hit->detectorID()==CID_CDH ){
	    //	    cout<<"trackID : "<<hit->trackID()<<endl;
	  }
	}
      }
      if( pim_track_id.size()>1 ) tup=(TNtuple*)gFile->Get("tup_eff_mass_fp_pipi_MC"), tup->Fill(MyMCTools::Ystar_mass(reacData));

      if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
      if( anaInfo->nFCharge()!=1 ) return;
      ForwardChargeInfo *fcInfo = anaInfo->forwardCharge(0);

      if( fcInfo->mass2byRK()<FC_P_MIN || FC_P_MAX<fcInfo->mass2byRK() ) return;

      if( anaInfo->nCDS(CDS_PiMinus)!=2 ) return;
      tup=(TNtuple*)gFile->Get("tup_eff_mass_fp_pipi_data"), tup->Fill(MyMCTools::Ystar_mass(reacData));

      BeamInfo *beam=anaInfo->beam(0);
      CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
      CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
      if( pim0->dca()>pim1->dca() ) swap(pim0, pim1);
      TLorentzVector beam_lmom=beam->lmom();

      TLorentzVector pim_lmom0=pim0->lmom();
      TLorentzVector pim_lmom1=pim1->lmom();

      TLorentzVector fp_lmom=fcInfo->lmom();
      TLorentzVector tgt_lmom=MyAnaTools::target_lmom();

      double mm=(beam_lmom+tgt_lmom-fp_lmom).M();
      double kppimpim_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom0-pim_lmom1).M();
      double kppim0_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom0).M();
      double kppim1_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom1).M();
      bool mmL_flag=false; bool mmS0_flag=false; bool mmP_flag=false; bool mmPgamma_flag=false;
      if( FC_mmL_MIN<kppim0_mm && kppim0_mm<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmL_MIN<kppim1_mm && kppim1_mm<FC_mmL_MAX ) mmL_flag=true;
      if( FC_mmS0_MIN<kppim0_mm && kppim0_mm<FC_mmS0_MAX ) mmS0_flag=true;
      if( FC_mmS0_MIN<kppim1_mm && kppim1_mm<FC_mmS0_MAX ) mmS0_flag=true;
      if( FC_mmP_MIN<kppimpim_mm && kppimpim_mm<FC_mmP_MAX ) mmP_flag=true;
      if( 0.99<kppimpim_mm && kppimpim_mm<1.06 ) mmPgamma_flag=true;

      if( mode==1 ){
	if( mmP_flag ){
	  tup=(TNtuple*)gFile->Get("tup_eff_mass_fp_pipi_mm1"), tup->Fill(MyMCTools::Ystar_mass(reacData));
	  if( mmL_flag ){
	    tup=(TNtuple*)gFile->Get("tup_eff_mass_fp_pipi_mm2"), tup->Fill(MyMCTools::Ystar_mass(reacData));
	  }
	}
      }
      else if( mode==2 ){
	if( mmPgamma_flag ){
	  tup=(TNtuple*)gFile->Get("tup_eff_mass_fp_pipi_mm1"), tup->Fill(MyMCTools::Ystar_mass(reacData));
	  if( mmS0_flag ){
	    tup=(TNtuple*)gFile->Get("tup_eff_mass_fp_pipi_mm2"), tup->Fill(MyMCTools::Ystar_mass(reacData));
	  }
	}
      }      
    }
  }
}

void fillHistMCAcc_fn(DetectorData *detData, MCData* mcData, ReactionData *reacData, AnaInfo *anaInfo,
			BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan)
{
  TNtuple *tup;
  if( MyMCTools::effNC(detData) ){
    MyHistTools::fillTH("trigY_mass_fn", MyMCTools::Ystar_mass(reacData));
    tup=(TNtuple*)gFile->Get("tup_trig_mass_fn"), tup->Fill(MyMCTools::Ystar_mass(reacData));
    bool pim_hit=false, pip_hit=false;
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit=detData->detectorHit(i);
      if( hit->detectorID()==CID_CDH && hit->pdg()==-211 ) pim_hit=true;
      if( hit->detectorID()==CID_CDH && hit->pdg()== 211 ) pip_hit=true;
    }

    if( pim_hit && pip_hit ){
      MyHistTools::fillTH("effY_mass_fn_pipi_MC", MyMCTools::Ystar_mass(reacData));
      tup=(TNtuple*)gFile->Get("tup_eff_mass_fn_pipi_MC"), tup->Fill(MyMCTools::Ystar_mass(reacData));
    }
    vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
    vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
    if( CVChits.size()!=0 || BVChits.size()!=0 ) return;
    if( anaInfo->nFNeutral()!=1 ) return;
    ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
    if( fnInfo->pid()!=F_Neutron ) return;
    if( anaInfo->minDCA() && GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
    if( cdstrackMan->nGoodTrack()!=2 ) return;
    if( anaInfo->nCDS(CDS_PiMinus)!=1 || anaInfo->nCDS(CDS_PiPlus)!=1 ) return;
    if( !anaInfo->CDS(CDS_PiMinus, 0)->flag() || !anaInfo->CDS(CDS_PiPlus, 0)->flag() ) return;
    if( !anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->flag() ){ cout<<"  !!! CDS2 pi+ pi- flag is false !!!"<<endl; return; }

    MyHistTools::fillTH("effY_mass_fn_pipi_data", MyMCTools::Ystar_mass(reacData));
    tup=(TNtuple*)gFile->Get("tup_eff_mass_fn_pipi_data"), tup->Fill(MyMCTools::Ystar_mass(reacData));

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

    bool K0_flag=false; bool Sm_flag=false; bool Sp_flag=false;
    if( Npipi_K0_MIN<pipi->im() && pipi->im()<Npipi_K0_MAX ) K0_flag=true;
    if( Npipi_Sm_MIN<(fn_lmom+pim_lmom).M() && (fn_lmom+pim_lmom).M()<Npipi_Sm_MAX ) Sm_flag=true;
    if( Npipi_Sp_MIN<(fn_lmom+pip_lmom).M() && (fn_lmom+pip_lmom).M()<Npipi_Sp_MAX ) Sp_flag=true;

    if( Npipi_N_MIN<kn_pipi_mm && kn_pipi_mm<Npipi_N_MAX ){
      MyHistTools::fillTH("effY_mass_fn_pipi_mmN", MyMCTools::Ystar_mass(reacData));
      tup=(TNtuple*)gFile->Get("tup_eff_mass_fn_pipi_mmN"), tup->Fill(MyMCTools::Ystar_mass(reacData));
      if( !K0_flag ){
	tup=(TNtuple*)gFile->Get("tup_eff_mass_fn_pipi_mmN_woK0"), tup->Fill(MyMCTools::Ystar_mass(reacData));
      }
      if( !K0_flag && !Sm_flag && !Sp_flag ){
	MyHistTools::fillTH("KN_MM_diff_MC", MyMCTools::Ystar_mass(reacData), MyMCTools::Ystar_mass(reacData)-kn_lmom.M());
	tup=(TNtuple*)gFile->Get("tup_eff_mass_fn_pipi_mmN_woAll"), tup->Fill(MyMCTools::Ystar_mass(reacData));
	MyHistTools::fillTH("effY_mass_fn_pipi_mmN_woAll", MyMCTools::Ystar_mass(reacData));
      }
    }
  }
}
