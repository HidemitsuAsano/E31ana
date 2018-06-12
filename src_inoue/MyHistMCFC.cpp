#include "MyHistMCFC.h"

using namespace std;

void initHistMCFC()
{
  new TH2F("FL_diff_Arc_MC", "", 1000, -10, 10, 1000, 1300, 1500);
  new TH2F("FL_diff_RK_MC", "", 1000, -10, 10, 1000, 1300, 1500);

  new TH2F("FL_diff_diffCou_MC", "", 1000, -10, 10, 1000, -25, 25);
  new TH2F("FL_diff_FDC1_MC", "", 1000, -10, 10, 1000, -25, 25);
  new TH2F("FL_diff_Vtx_MC", "", 1000, -10, 10, 1000, -25, 25);

  new TH2F("mom_Ang_MC", "", 1500, 0.0, 1.5, 1500, 0.0, 1.5);
  new TH2F("mom_RK_MC", "",  1500, 0.0, 1.5, 1500, 0.0, 1.5);
  new TH2F("mom_TOF_MC", "", 1500, 0.0, 1.5, 1500, 0.0, 1.5);

  new TH2F("mom_diff_Ang_MC", "", 1000, -0.1, 0.1,   1500, 0.0, 1.5);
  new TH2F("mom_diff_RK_MC", "",  1000, -0.1, 0.1,   1500, 0.0, 1.5);
  new TH2F("mom_diff_TOF_MC", "",   1000, -0.1, 0.1, 1500, 0.0, 1.5);

  new TH2F("FDC1_pos_diff_MC", "", 1000, -1.0, 1.0, 1000, -1.0, 1.0);
  // new TH2F("fitCouPos_XZ", "", 10000, -500, 500, 10000, 1250, 1750);
  // for( int seg=1; seg<=34; seg++ ){
  //   new TH2F(Form("fitCouPos_XZ_CVC%d", seg), "", 10000, -500, 500, 10000, 1250, 1750);
  // }
  // for( int seg=1; seg<=27; seg++ ){
  //   new TH2F(Form("fitCouPos_XZ_PC%d", seg), "", 10000, -500, 500, 10000, 1250, 1750);
  // }

  new TH2F("hitpos_diff_XZ_RK_MC", "", 1000, -100, 100, 1000, -100, 100);
  new TH1F("hitpos_diff_XZ", "", 1000, -100, 100);
}

void fillHistMCFC(ConfMan *conf, AnaInfo *anaInfo, CDSHitMan *cdsMan, BeamLineHitMan *blMan,
                  CDSTrackingMan *cdstrackMan, BeamLineTrackMan *bltrackMan,
                  ReactionData *reacData, MCData *mcData, DetectorData *detData)
{
  if( anaInfo->nFCharge()==1 && anaInfo->forwardCharge(0)->pid()==F_Proton ){
    ForwardChargeInfo *fcInfo = anaInfo->forwardCharge(0);
    BeamInfo *beam = anaInfo->beam(0);
    //    fcInfo-> dump();
    HodoscopeLikeHit *fc_hit=fcInfo->hodo(blMan);
    DetectorHit *fc_hit_mc=MyMCTools::getHit(detData, fc_hit);

    Track *p_track = MyMCTools::track(mcData, fc_hit_mc->trackID());
    if( !p_track ){
      return;
    }
    if( p_track->pdgID()!=2212 ){
      return;
    }
    double fl_true=0.1*p_track->FlightLength();
    double mom_true=0.001*p_track->momentum().Mag();
   
    TVector3 momRK=fcInfo->momentum();
    momRK.SetMag(fcInfo->momByRK());
    TVector3 momMC=0.001*p_track->momentum();
    // cout<<"mom by RK : "<<momRK.X()<<", "<<momRK.Y()<<", "<<momRK.Z()<<"   "<<momRK.Mag()<<endl;
    // cout<<"mom by MC : "<<momMC.X()<<", "<<momMC.Y()<<", "<<momMC.Z()<<"   "<<momMC.Mag()<<endl;

    MyHistTools::fillTH("FL_diff_Arc_MC", fcInfo->flByArc()-fl_true, fl_true);
    MyHistTools::fillTH("FL_diff_RK_MC",  fcInfo->flByRK()-fl_true, fl_true);
    MyHistTools::fillTH("mom_diff_Ang_MC", fcInfo->momByAng()-mom_true, mom_true);
    MyHistTools::fillTH("mom_diff_RK_MC",  fcInfo->momByRK()-mom_true, mom_true);
    MyHistTools::fillTH("mom_diff_TOF_MC",  fcInfo->momByTOF()-mom_true, mom_true);

    MyHistTools::fillTH("mom_Ang_MC", fcInfo->momByAng(), mom_true);
    MyHistTools::fillTH("mom_RK_MC",  fcInfo->momByRK(),  mom_true);
    MyHistTools::fillTH("mom_TOF_MC", fcInfo->momByTOF(), mom_true);

    double diff_fl_RK=fcInfo->flByRK()-fl_true;
    TVector3 fit_vertex=fcInfo->fitVertex();
    TVector3 vertex=fcInfo->vertex();
    TVector3 FDC1pos=fcInfo->FDC1().pos();
    DetectorHit *hitFDC1[6];
    HodoscopeLikeHit *hit =fcInfo->hodo(blMan);
    DetectorHit *hitMC=MyMCTools::getHit(detData, hit);
    int nFDC1[6]={0, 0, 0, 0, 0, 0};
    for( int i=0; i<detData->detectorHitSize(); i++ ){
      DetectorHit *hit=detData->detectorHit(i);
      if( hit->detectorID()==CID_FDC1 ){
        nFDC1[hit->layerID()]++;
	hitFDC1[hit->layerID()]=hit;
      }
    }
    // for( int i=0; i<6; i++ ){
    //   cout<<"FDC1  layer "<<i<<"  n hit : "<<nFDC1[i]<<endl;
    //   TVector3 pos=0.1*hitFDC1[i]->pos();
    //   cout<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<endl;
    // }
    TVector3 FDC1pos_MC(-9999, -9999, -9999);
    if( hitFDC1[2] && hitFDC1[3] ){
      TVector3 pos1=0.1*hitFDC1[2]->pos();
      TVector3 pos2=0.1*hitFDC1[3]->pos();
      TVector3 dir=(pos2-pos1).Unit();
      FDC1pos_MC=pos1+((FDC1pos.Z()-pos1.Z())/dir.Z())*dir;
    }
    MyHistTools::fillTH("FDC1_pos_diff_MC", (FDC1pos-FDC1pos_MC).X(), (FDC1pos-FDC1pos_MC).Y());

    // cout<<"FDC1 pos by data : "<<FDC1pos.X()<<", "<<FDC1pos.Y()<<", "<<FDC1pos.Z()<<endl;
    // cout<<"FDC1 pos by MC   : "<<FDC1pos_MC.X()<<", "<<FDC1pos_MC.Y()<<", "<<FDC1pos_MC.Z()<<endl;
    TVector3 hitposMC=hitMC->pos();
    TVector3 gpos, grot;
    TVector3 cpos, crot;
    double param[20];
    conf->GetGeomMapManager()-> GetPos(hit->cid(), 0, gpos);
    conf-> GetGeomMapManager()-> GetRot(hit->cid(), 0, grot);
    conf-> GetGeomMapManager()-> GetParam(hit->cid(), hit->seg(), param);
    conf->GetGeomMapManager()-> GetPos(hit->cid(), hit->seg(), cpos);
    conf-> GetGeomMapManager()-> GetRot(hit->cid(), hit->seg(), crot);

    TVector3 fitFDC1=fcInfo->fitFDC1();
    MyHistTools::fillTH("FL_diff_Vtx_MC", diff_fl_RK, (fit_vertex-vertex).Mag());
    MyHistTools::fillTH("FL_diff_FDC1_MC", diff_fl_RK, (FDC1pos-fitFDC1).Mag());
    MyHistTools::fillTH("FL_diff_diffCou_MC", diff_fl_RK, fcInfo->diffCounter());

    TVector3 FChitpos=fcInfo->hitpos();
    TVector3 FChitpos_MC=0.1*fc_hit_mc->pos();
    TVector3 FChitpos_diff=FChitpos-FChitpos_MC;
    double FCdiff=sqrt(FChitpos_diff.X()*FChitpos_diff.X()+FChitpos_diff.Z()*FChitpos_diff.Z());
    MyHistTools::fillTH("hitpos_diff_XZ_RK_MC", FChitpos_diff.X(), FChitpos_diff.Z());
    MyHistTools::fillTH("hitpos_diff_XZ", FCdiff);
    //    cout<<FChitpos_MC.X()<<", "<<FChitpos_MC.Y()<<", "<<FChitpos_MC.Z()<<endl;

    //    MyHistTools::fillTH("fitCouPos_XZ", FChitpos.X(), FChitpos.Z());
    // if( hit->cid()==CID_CVC ){
    //   MyHistTools::fillTH(Form("fitCouPos_XZ_CVC%d", hit->seg()), FChitpos.X(), FChitpos.Z());
    // }
    // else{
    //   MyHistTools::fillTH(Form("fitCouPos_XZ_PC%d", hit->seg()), FChitpos.X(), FChitpos.Z());
    // }

    // cout<<"Track PDG : "<<p_track->pdgID()<<"  parent: "<<p_track->parentTrackID()<<endl;
    // cout<<"FL true : "<<0.1*p_track->FlightLength()<<"  Arc : "<<fcInfo->flByArc()<<"  RK : "<<fcInfo->flByRK()<<endl;
    // cout<<"mom true : "<<0.001*p_track->momentum().Mag()<<"  Arc : "<<fcInfo->momByAng()<<"  RK : "<<fcInfo->momByRK()<<"  TOF : "<<fcInfo->momByTOF()<<endl;
  }
}
