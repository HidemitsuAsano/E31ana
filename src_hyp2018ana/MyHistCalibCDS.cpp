#include "MyHistCalibCDS.h"

using namespace std;

void initHistCalibCDS()
{
  new TH2F("CDS_mass2_mom_Z40", "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH1F("CDS_mass2_plus_Z40", "CDS mass2 plus", 10000, -1.0, 9.0);
  new TH1F("CDS_mass2_minus_Z40", "CDS mass2 minus", 10000, -1.0, 9.0);

  new TH2F("CDS_mass2_mom",     "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_pim", "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_pip", "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_p",   "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_d",   "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);
  new TH2F("CDS_mass2_mom_km",  "CDS PID", 1000, -1.0, 9.0, 3000, -1.5, 1.5);

  new TH1F("CDS_mass2_plus", "CDS mass2 plus", 10000, -1.0, 9.0);
  new TH1F("CDS_mass2_minus", "CDS mass2 minus", 10000, -1.0, 9.0);

  new TH1F("CDC_chi2",     "CDS chi-sqaure", 5000, 0, 500);
  new TH2F("CDC_chi2_DCA", "CDS chi-square vs DCA", 1000, 0, 500, 1000, 0, 25);

  new TH1F("CDS_DCA",     "CDS DCA", 1000, 0, 25);
  new TH1F("CDS_DCA_pim", "CDS DCA", 1000, 0, 25);
  new TH1F("CDS_DCA_pip", "CDS DCA", 1000, 0, 25);
  new TH1F("CDS_DCA_p",   "CDS DCA", 1000, 0, 25);
  new TH1F("CDS_DCA_d",   "CDS DCA", 1000, 0, 25);
  new TH1F("CDS_DCA_km",  "CDS DCA", 1000, 0, 25);

  for( int i=0; i<20; i++ ){
    int min=50*i; int max=50*(i+1);
    new TH1F(Form("CDS_mass2_p%d_%d", min, max), Form("CDS mass2  %d ~ %d", min, max), 10000, -1.0, 9.0);
    new TH1F(Form("CDS_mass2_m%d_%d", min, max), Form("CDS mass2 -%d ~-%d", min, max), 10000, -1.0, 9.0);
  }

  for( int seg=1; seg<=24; seg++ ){
    new TH1F(Form("IH_ADC_%d", seg),     Form("IH seg%d ADC", seg), 4000, 0, 4000);
    new TH1F(Form("IH_ADC_%d_MIP", seg), Form("IH seg%d ADC", seg), 4000, 0, 4000);
  }

  for( int seg=1; seg<=36; seg++ ){
    new TH1F(Form("CDH_%d_offset", seg), Form("CDH seg%d time offset", seg), 10000, -50, 50);
    new TH1F(Form("CDH_%d_offset_pim", seg), Form("CDH seg%d time offset", seg), 10000, -50, 50);
    new TH1F(Form("CDH_%d_offset_pip", seg), Form("CDH seg%d time offset", seg), 10000, -50, 50);
    new TH1F(Form("CDH_%d_offset_km", seg),  Form("CDH seg%d time offset", seg), 10000, -50, 50);
    new TH1F(Form("CDH_%d_offset_p", seg),   Form("CDH seg%d time offset", seg), 10000, -50, 50);
  }
}

void fillHistCalibCDS(EventHeader *header, ConfMan *confMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  for( int i=0; i<cdsMan->nIH(); i++ ){
    HodoscopeLikeHit *hit=cdsMan->IH(i);
    MyHistTools::fillTH(Form("IH_ADC_%d", hit->seg()), hit->adcu());
  }

  for( int i=0; i<anaInfo->nCDS(); i++ ){
    CDSInfo *cds=anaInfo->CDS(i);
    CDSTrack *track=anaInfo->CDS(i)->track(cdstrackMan);
    //    cout<<track->FittingLevel()<<endl;

    MyHistTools::fillTH("CDC_chi2", track->Chi());
    MyHistTools::fillTH("CDC_chi2_DCA", track->Chi(), cds->dca());
    if( !cds->flag() || GeomTools::GetID(cds->vertexBeam())!=CID_Fiducial ) continue;

    TVector3 vtxCDH=track->CDHVertex();
    if( fabs(vtxCDH.Z())>40 ){
      MyHistTools::fillTH("CDS_mass2_mom_Z40", cds->mass2(), cds->mom());
      if( cds->mom()>0 ) MyHistTools::fillTH("CDS_mass2_plus_Z40", cds->mass2());
      if( cds->mom()<0 ) MyHistTools::fillTH("CDS_mass2_minus_Z40", cds->mass2());

      anaInfo->CDS(i)->SetFlag(false);
      continue;
    }

    if( fabs(cds->mom())>0.3 && (cds->pid()==CDS_PiMinus || cds->pid()==CDS_PiPlus) ){
      if( MyTools::searchIHHit(cds->track(cdstrackMan), cdsMan) ){
	HodoscopeLikeHit *ih=cds->track(cdstrackMan)->IHHit(cdsMan, 0);
	MyHistTools::fillTH(Form("IH_ADC_%d_MIP", ih->seg()), ih->adcu());
      }
    }
    HodoscopeLikeHit *CDH=cds->CDH(cdsMan);
    double CDHoffset=CDH->ctmean()-anaInfo->beam(0)->T0time()-cds->offset();
    MyHistTools::fillTH("CDS_mass2_mom", cds->mass2(), cds->mom());
    MyHistTools::fillTH("CDS_DCA", cds->dca());
    int pid=confMan->GetCDSFittingParamManager()->PID(cds->mom(), cds->mass2());
    if( pid==CDS_PiMinus  ){
      MyHistTools::fillTH("CDS_mass2_mom_pim", cds->mass2(), cds->mom());
      MyHistTools::fillTH("CDS_DCA_pim", cds->dca());
      if( cds->dca()<3.0 ){
	MyHistTools::fillTH(Form("CDH_%d_offset", cds->CDHseg()), CDHoffset);
	MyHistTools::fillTH(Form("CDH_%d_offset_pim", cds->CDHseg()), CDHoffset);
      }
    }
    if( pid==CDS_Kaon     ){
      MyHistTools::fillTH("CDS_mass2_mom_km", cds->mass2(), cds->mom());
      MyHistTools::fillTH("CDS_DCA_km", cds->dca());
      if( cds->dca()<3.0 ){
	MyHistTools::fillTH(Form("CDH_%d_offset", cds->CDHseg()), CDHoffset);
	MyHistTools::fillTH(Form("CDH_%d_offset_km", cds->CDHseg()), CDHoffset);
      }
    }
    if( pid==CDS_PiPlus   ){
      MyHistTools::fillTH("CDS_mass2_mom_pip", cds->mass2(), cds->mom());
      MyHistTools::fillTH("CDS_DCA_pip", cds->dca());
      if( cds->dca()<3.0 ){
	MyHistTools::fillTH(Form("CDH_%d_offset", cds->CDHseg()), CDHoffset);
	MyHistTools::fillTH(Form("CDH_%d_offset_pip", cds->CDHseg()), CDHoffset);
      }
    }
    if( pid==CDS_Proton   ){
      MyHistTools::fillTH("CDS_mass2_mom_p", cds->mass2(), cds->mom());
      MyHistTools::fillTH("CDS_DCA_p", cds->dca());
      if( cds->dca()<3.0 ){
	MyHistTools::fillTH(Form("CDH_%d_offset", cds->CDHseg()), CDHoffset);
	MyHistTools::fillTH(Form("CDH_%d_offset_p", cds->CDHseg()), CDHoffset);
      }
    }
    if( pid==CDS_Deuteron ){
      MyHistTools::fillTH("CDS_mass2_mom_d", cds->mass2(), cds->mom());
      MyHistTools::fillTH("CDS_DCA_d", cds->dca());
    }

    if( cds->mom()>0 ){
      MyHistTools::fillTH("CDS_mass2_plus", cds->mass2());
      for( int i=1; i<20; i++ ){
	if( fabs(cds->mom())<0.05*i ){
	  MyHistTools::fillTH(Form("CDS_mass2_p%d_%d",50*(i-1), 50*i), cds->mass2());
	  break;
	}
      }
    }
    else{
      MyHistTools::fillTH("CDS_mass2_minus", cds->mass2());
      for( int i=1; i<20; i++ ){
	if( fabs(cds->mom())<0.05*i ){
	  MyHistTools::fillTH(Form("CDS_mass2_m%d_%d",50*(i-1), 50*i), cds->mass2());
	  break;
	}
      }
    }
  }
}
