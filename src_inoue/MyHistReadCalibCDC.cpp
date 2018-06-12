#include "MyHistReadCalibCDC.h"

using namespace std;

void initHistReadCalibCDC(ConfMan *conf)
{
  new TH1F("CDC_chi2", "CDC chi2", 1000, 0, 100);

  new TH2F("CDC_dt_dl",   "CDC dt vs dl",  200, -20, 380, 200,  -0.1,   1.1);
  new TH2F("CDC_dt_res",  "CDC dt vs res", 200, -20, 380, 200,  -0.25, 0.25);
  new TH2F("CDC_ob2_res", "CDC dt vs res", 200, 0,   100, 200,  -25, 25);
  new TH2F("CDC_ob2_res2", "CDC dt vs res", 200, 0,   100, 200, -25, 25);

  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    new TH1F(Form("CDC_res_%d", lay), Form("CDC_res_%d", lay), 1000, -1.0, 1.0);
    new TH1F(Form("CDC_res2_%d", lay), Form("CDC_res_%d", lay), 1000, -1.0, 1.0);
    new TH2F(Form("CDC_dt_dl_%d", lay), Form("CDC lay%d dt vs dl ", lay),  200, -20, 380, 200,  -0.1,   1.1);
    new TH2F(Form("CDC_dt_res_%d", lay), Form("CDC lay%d dt vs res", lay), 200, -20, 380, 200,  -0.25, 0.25);
    new TH2F(Form("CDC_ob2_res_%d", lay), Form("CDC lay%d 1/#beta^{2} vs res", lay),  200, 0, 100, 200,  -25, 25);
    new TH2F(Form("CDC_ob2_res2_%d", lay), Form("CDC lay%d 1/#beta^{2} vs res", lay), 200, 0, 100, 200,  -25, 25);

    int nwire=conf->GetCDCWireMapManager()->nw(lay);
    for( int wire=1; wire<=nwire; wire++ ){
      new TH1F(Form("CDC_res_%d_%d", lay, wire), Form("CDC_res_%d_%d", lay, wire), 1000, -1.0, 1.0);
      new TH1F(Form("CDC_res2_%d_%d", lay, wire), Form("CDC_res_%d_%d", lay, wire), 1000, -1.0, 1.0);
      new TH2F(Form("CDC_dt_dl_%d_%d", lay, wire), Form("CDC lay%d wire%d dt vs dl ", lay, wire),  200, -20, 380, 200,  -0.1,   1.1);
      new TH2F(Form("CDC_dt_res_%d_%d", lay, wire), Form("CDC lay%d wire%d dt vs res", lay, wire), 200, -20, 380, 200,  -0.25, 0.25);
      new TH2F(Form("CDC_ob2_res_%d_%d", lay, wire), Form("CDC lay%d wire%d 1/#beta^{2} vs res", lay, wire),  200, 0, 100, 200,  -25, 25);
      new TH2F(Form("CDC_ob2_res2_%d_%d", lay, wire), Form("CDC lay%d wire%d 1/#beta^{2} vs res", lay, wire), 200, 0, 100, 200,  -25, 25);
    }
  }

  new TH1F("CDH_dphi", "CDH_ang", 1000, -30, 30);
  new TH2F("CDC_Z_CDH_tsub", "", 1000, -75, 75, 1000, -5, 5);
  new TH2F("CDC_Z_CDH_hitpos", "", 1000, -75, 75, 1000, -75, 75);

  new TH2F("CDC_Z_CDH_tsub2", "", 1000, -75, 75, 1000, -5, 5);
  new TH2F("CDC_Z_CDH_hitpos2", "", 1000, -75, 75, 1000, -75, 75);
  for( int seg=1; seg<=36; seg++ ){
    new TH1F(Form("CDH%d_dphi",seg), "CDH_ang", 1000, -30, 30);
    new TH2F(Form("CDC_Z_CDH%d_tsub", seg), "", 1000, -75, 75, 1000, -5, 5);
    new TH2F(Form("CDC_Z_CDH%d_hitpos", seg), "", 1000, -75, 75, 1000, -75, 75);

    new TH2F(Form("CDC_Z_CDH%d_tsub2", seg), "", 1000, -75, 75, 1000, -5, 5);
    new TH2F(Form("CDC_Z_CDH%d_hitpos2", seg), "", 1000, -75, 75, 1000, -75, 75);
  }

  new TH1F("IH_dphi", "IH_ang", 1000, -30, 30);
  for( int seg=1; seg<=24; seg++ ){
    new TH1F(Form("IH%d_dphi",seg), "IH_ang", 1000, -30, 30);
  }
}

void fillHistReadCalibCDC(ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;
  BeamInfo *beam=anaInfo->beam(0);
  if( !beam->flag() ) return;

  std::vector<HodoscopeLikeHit*> CDHhits=MyTools::getCDH(cdsMan);
  if( CDHhits.size()==1 ){
    HodoscopeLikeHit *hit = CDHhits[0];
    double min_dphi=DBL_MAX;
    if( cdstrackMan->nGoodTrack()==1 ){
      CDSTrack *track =cdstrackMan->GoodTrack(0);
      double tmpparam[5];
      TVector3 tmp;
      track->GetParameters(CID_CDC,tmpparam,tmp);
      TVector3 CDCpos =MathTools::CalcHelixPosatR(tmpparam, 54.4);
      
      TVector3 pos;
      conf->GetGeomMapManager()->GetGPos(CID_CDH, hit->seg(), pos);
      double tmp_dphi=CDCpos.Phi()-pos.Phi();
      if(tmp_dphi>TMath::Pi()) tmp_dphi=TMath::Pi()*2-tmp_dphi;
      tmp_dphi*=TMath::RadToDeg();
      
      if( fabs(tmp_dphi)<fabs(min_dphi) ) min_dphi=tmp_dphi;

      if( fabs(tmp_dphi)<7 ){
      }
    }
    MyHistTools::fillTH("CDH_dphi", min_dphi);
    MyHistTools::fillTH(Form("CDH%d_dphi",hit->seg()), min_dphi);
  }

  std::vector<HodoscopeLikeHit*> IHhits=MyTools::getIH(cdsMan);
  if( IHhits.size()==1 ){
    HodoscopeLikeHit *hit = IHhits[0];
    double min_dphi=DBL_MAX;
    if( cdstrackMan->nGoodTrack()==1 ){
      CDSTrack *track =cdstrackMan->GoodTrack(0);
      double tmpparam[5];
      TVector3 tmp;
      track->GetParameters(CID_CDC,tmpparam,tmp);
      TVector3 CDCpos =MathTools::CalcHelixPosatR(tmpparam, 14.1);
      
      TVector3 pos;
      conf->GetGeomMapManager()->GetGPos(CID_IH, hit->seg(), pos);
      double tmp_dphi=CDCpos.Phi()-pos.Phi();
      if(tmp_dphi>TMath::Pi()) tmp_dphi=TMath::Pi()*2-tmp_dphi;
      tmp_dphi*=TMath::RadToDeg();

      if( fabs(tmp_dphi)<fabs(min_dphi) ) min_dphi=tmp_dphi;  
    }
    MyHistTools::fillTH("IH_dphi", min_dphi);
    MyHistTools::fillTH(Form("IH%d_dphi",hit->seg()), min_dphi);
  }

  for( int i=0; i<anaInfo->nCDS(); i++ ){
    CDSInfo *cds=anaInfo->CDS(i);
    if( GeomTools::GetID(cds->vertexBeam())!=CID_Fiducial ) continue;
    CDSTrack *track=cds->track(cdstrackMan);
    track-> SetHitPos(cdsMan);
    //    cout<<"Fitting Level : "<<track->FittingLevel()<<endl;

    MyHistTools::fillTH("CDC_chi2", track->Chi());
    if( track->Chi()>30 ) continue;

    bool single=true;
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ) if( track->nTrackHit(layer)!=1 ) single=false;
    if( !single ) continue;

    for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
      CDCHit *cdc = track->TrackHit(cdsMan, lay, 0);
      int wire=cdc->wire();
      double res=cdc->resl();
      double dt=cdc->dt();
      double dl=cdc->dl();
      double dlr=dl-res;
      double dxdt=conf->GetXTMapManager()->CalcDxDt(CID_CDC, lay, wire, dt);
      MyHistTools::fillTH(Form("CDC_res_%d", lay), res);
      MyHistTools::fillTH(Form("CDC_res_%d_%d", lay, wire), res);

      MyHistTools::fillTH("CDC_dt_dl", dt, dlr);
      MyHistTools::fillTH("CDC_dt_res", dt, res);
      MyHistTools::fillTH("CDC_ob2_res", 1./(cds->beta()*cds->beta()), res/dxdt);

      MyHistTools::fillTH(Form("CDC_dt_dl_%d", lay), dt, dlr);
      MyHistTools::fillTH(Form("CDC_dt_res_%d", lay), dt, res);
      MyHistTools::fillTH(Form("CDC_ob2_res_%d", lay), 1./(cds->beta()*cds->beta()) ,res/dxdt);

      MyHistTools::fillTH(Form("CDC_dt_dl_%d_%d", lay, wire), dt, dlr);
      MyHistTools::fillTH(Form("CDC_dt_res_%d_%d", lay, wire), dt, res);
      MyHistTools::fillTH(Form("CDC_ob2_res_%d_%d", lay, wire), 1./(cds->beta()*cds->beta()) ,res/dxdt);
    }

    if( cds->beta()==DEFAULTD ) continue;

    HodoscopeLikeHit *hit=cds->CDH(cdsMan);
    TVector3 CDCpos=track->CDHVertex();
    MyHistTools::fillTH("CDC_Z_CDH_tsub", CDCpos.Z(), hit->ctsub());
    MyHistTools::fillTH(Form("CDC_Z_CDH%d_tsub", hit->seg()), CDCpos.Z(), hit->ctsub());
    MyHistTools::fillTH("CDC_Z_CDH_hitpos", CDCpos.Z(), hit->hitpos());
    MyHistTools::fillTH(Form("CDC_Z_CDH%d_hitpos", hit->seg()), CDCpos.Z(), hit->hitpos());

    track-> Retiming(cdsMan, conf, cds->beta(), true);
    track-> SetHitPos(cdsMan);

    CDCpos=track->CDHVertex();
    MyHistTools::fillTH("CDC_Z_CDH_tsub2", CDCpos.Z(), hit->ctsub());
    MyHistTools::fillTH(Form("CDC_Z_CDH%d_tsub2", hit->seg()), CDCpos.Z(), hit->ctsub());
    MyHistTools::fillTH("CDC_Z_CDH_hitpos2", CDCpos.Z(), hit->hitpos());
    MyHistTools::fillTH(Form("CDC_Z_CDH%d_hitpos2", hit->seg()), CDCpos.Z(), hit->hitpos());
    //    cout<<"Fitting Level2 : "<<track->FittingLevel()<<endl;
    for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
      CDCHit *cdc = track->TrackHit(cdsMan, lay, 0);
      int wire=cdc->wire();
      double res=cdc->resl();
      double dt=cdc->dt();
      double dxdt=conf->GetXTMapManager()->CalcDxDt(CID_CDC, lay, wire, dt);
      MyHistTools::fillTH(Form("CDC_res2_%d", lay), res);
      MyHistTools::fillTH(Form("CDC_res2_%d_%d", lay, wire), res);
      MyHistTools::fillTH("CDC_ob2_res2", 1./(cds->beta()*cds->beta()), res/dxdt);
      MyHistTools::fillTH(Form("CDC_ob2_res2_%d", lay), 1./(cds->beta()*cds->beta()) ,res/dxdt);
      MyHistTools::fillTH(Form("CDC_ob2_res2_%d_%d", lay, wire), 1./(cds->beta()*cds->beta()) ,res/dxdt);
    }
  }
}
