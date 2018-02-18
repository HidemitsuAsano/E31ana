// MyAnalysisBLDCTest.cpp

#include "MyAnalysisBLDCTest.h"

MyAnalysisBLDCTest::MyAnalysisBLDCTest(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	Clear();
}

MyAnalysisBLDCTest::~MyAnalysisBLDCTest()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisBLDCTest::Clear()
{
	return;
}

bool MyAnalysisBLDCTest::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
	rtFile->cd();
	DetectorList *dlist=DetectorList::GetInstance();

	const int ndc=6;
	int dccid[ndc]={CID_BLC1a,CID_BLC1b,
		CID_BLC2a,CID_BLC2b,
		CID_BPC,
		CID_FDC1
	};
  char LayConfig[ndc][16] = {
    {"UUVVUUVV"},
    {"UUVVUUVV"},
    {"UUVVUUVV"},
    {"VVUUVVUU"},
    {"XXYYXXYY"},
    {"UUXXVV"},
  };

	int Mult[ndc][8];
	for(int idc=0;idc<ndc; idc++){
		for(int layer=0;layer<8;layer++){
			Mult[idc][layer] = -1;
		}
	}
	for(int idc=0;idc<ndc;idc++){
		const int cid= dccid[idc];
		const int nlays= dlist->GetNlayers(cid);
		for( int layer=1; layer<=nlays; layer++ ){
			Mult[idc][layer-1]=blMan->nBLDC(cid,layer);
		}
	}

  /* Basic data filling */
	for(int idc=0;idc<ndc;idc++){
		const int cid= dccid[idc];
		const char* name= dlist->GetName(cid).data();
		const int nlays= dlist->GetNlayers(cid);
		/* Selection */
		if(Mult[idc][0]<1 || Mult[idc][nlays-1]<1){
      continue;
    }
    for(int layer=1;layer<=nlays;layer++){
      int mult = blMan->nBLDC(cid,layer);
      FillHist(Form("%sl%d_Multiplicity",name,layer),mult);
      for(int i=0;i<mult;i++){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        int wire = hit->wire();
        int tdc = hit->tdc();
        double dt = hit->dt();
        double dl = hit->dl();
        FillHist(Form("%sl%d_HitPattern",name,layer),wire,1);
        FillHist(Form("%sl%d_TDC",name,layer),tdc,1);
        FillHist(Form("%sl%d_Time",name,layer),dt,1);
        FillHist(Form("%sl%d_Length",name,layer),dl,1);
        FillHist(Form("%sl%d_TimevsLength",name,layer),dt,dl,1);
        FillHist(Form("%sl%dw%02d_TDC",name,layer,wire),tdc,1);
        FillHist(Form("%sl%dw%02d_Time",name,layer,wire),dt,1);
        FillHist(Form("%sl%dw%02d_Length",name,layer,wire),dl,1);
        FillHist(Form("%sl%dw%02d_TimevsLength",name,layer,wire),dt,dl,1);
        FillHist(Form("%s_HitPattern%c",name,LayConfig[idc][layer-1]),layer,wire,1);
      }
    }
  }

  /* Hit efficiency evalutation */
  for(int idc=0; idc<ndc; idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    /* 1 - nlays hit */
    if(blMan->nBLDC(cid,1)==1&&blMan->nBLDC(cid,nlays)==1){
      FillHist(Form("%s_HitEfficiency1%d",name,nlays),0,1);
      for(int layer=1; layer<=nlays; layer++){
        if(blMan->nBLDC(cid,layer)>=1){
          FillHist(Form("%s_HitEfficiency1%d",name,nlays),layer,1);
        }
      }
    }
    /* 1 - (nlays-1) hit */
    if(blMan->nBLDC(cid,1)==1&&blMan->nBLDC(cid,nlays-1)==1){
      FillHist(Form("%s_HitEfficiency1%d",name,nlays-1),0,1);
      for(int layer=1; layer<=nlays; layer++){
        if(blMan->nBLDC(cid,layer)>=1){
          FillHist(Form("%s_HitEfficiency1%d",name,nlays-1),layer,1);
        }
      }
    }
    /* 2 - nlays hit */
    if(blMan->nBLDC(cid,2)==1&&blMan->nBLDC(cid,nlays)==1){
      FillHist(Form("%s_HitEfficiency2%d",name,nlays),0,1);
      for(int layer=1; layer<=nlays; layer++){
        if(blMan->nBLDC(cid,layer)>=1){
          FillHist(Form("%s_HitEfficiency2%d",name,nlays),layer,1);
        }
      }
    }
  }


  /* Tracking */
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    FillHist(Form("%s_TrackEfficiency",name),0,1);
    bltrackMan->LocalTracking(blMan,conf,cid,"mwpc");
    FillHist(Form("%s_nTrack",name),bltrackMan->ntrackBLDC(cid),1);
    if(bltrackMan->ntrackBLDC(cid)>=1){
      FillHist(Form("%s_TrackEfficiency",name),bltrackMan->ntrackBLDC(cid),1);
    }
    for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
      LocalTrack* track = bltrackMan->trackBLDC(cid,i);
      double chi2 = track -> chi2all();
      FillHist(Form("%s_ChiSquare",name),chi2,1);

      double xpos=-9999.9, ypos=-9999.9;
      track -> XYLocalPosatZ(0,xpos,ypos);
      double dx = track->dx();
      double dy = track->dy();
      FillHist(Form("%s_TrackPositionXY",name),xpos,ypos,1);
      FillHist(Form("%s_TrackDirectionXY",name),dx,dy,1);
      FillHist(Form("%s_TrackPositionvsDirectionX",name),xpos,dx,1);
      FillHist(Form("%s_TrackPositionvsDirectionY",name),ypos,dy,1);

      int mult[nlays]; for(int layer=1;layer<=nlays;layer++) { mult[layer-1]=0; }
      for(int ihit=0; ihit<track->nhit(); ihit++){
        ChamberLikeHit* hit = track->hit(ihit);
        int layer = hit->layer();
        int wire = hit->wire();
        double dt = hit->dt();
        double resl = hit->resl();
        mult[layer-1]++;
        FillHist(Form("%sl%d_HitPattern_inTrack",name,layer),wire,1);
        FillHist(Form("%sl%d_Time_inTrack",name,layer),dt,1);
        FillHist(Form("%sl%d_Residual_inTrack",name,layer),resl,1);
        FillHist(Form("%sl%d_TimevsResidual_inTrack",name,layer),dt,resl,1);
      }
      for(int layer=1; layer<=nlays; layer++){
        FillHist(Form("%sl%d_Multiplicity_inTrack",name,layer),mult[layer-1],1);
      }
    }
  }


  return true;

}

bool MyAnalysisBLDCTest::FillHist(TString name, double val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){ 
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLDCTest::FillHist(TString name, TString val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLDCTest::FillHist(TString name, double val1, double val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2&&weight>0){
    for(int i=0; i<weight; i++){
      h2 -> Fill(val1,val2,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLDCTest::FillHist(TString name, TString val1, TString val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2&&weight>0){
    for(int i=0; i<weight; i++){
      h2 -> Fill(val1,val2,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLDCTest::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisBLDCTest::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaBLDCTest");
  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  const int ndc=6;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,
    CID_BLC2a,CID_BLC2b,
    CID_BPC,
    CID_FDC1
  };
  int LayConfigNum[ndc] = {2,2,2,2,2,3};
  char LayConfig[ndc][6] = {
    {"UV"},
    {"UV"},
    {"UV"},
    {"VU"},
    {"XY"},
    {"UXV"},
  };

  /* For hit */
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    int maxwire = 0;
    for( int layer=1; layer<=nlays; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire(cid,layer);
      if(maxwire<nwire) maxwire = nwire;
      new TH1F( Form("%sl%d_Multiplicity_inTrack",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_HitPattern_inTrack",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_Time_inTrack",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),400,0,400);
      new TH1F( Form("%sl%d_Residual_inTrack",name,layer),Form("Residual %s-L%d;Residual (um);Counts",name,layer),200,0,200);
      new TH2F( Form("%sl%d_TimevsResidual_inTrack",name,layer),Form("Time vs. Residual %s-L%d;Time (ns);Residual (um)",name,layer),400,0,400,200,0,200);

      new TH1F( Form("%sl%d_Multiplicity",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_HitPattern",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_TDC",name,layer),Form("TDC %s-L%d;TDC ch.;Counts",name,layer),2000,0,2000);
      new TH1F( Form("%sl%d_Time",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),400,0,400);
      new TH1F( Form("%sl%d_Length",name,layer),Form("Length %s-L%d;Length (um);Counts",name,layer),200,0,200);
      new TH2F( Form("%sl%d_TimevsLength",name,layer),Form("Time vs. Length %s-L%d;Time (ns);Length (um)",name,layer),400,0,400,200,0,200);
      new TH1F( Form("%sl%d_Residual",name,layer),Form("Residual %s-L%d;Residual (um);Counts",name,layer),200,0,200);
      new TH2F( Form("%sl%d_TimevsResidual",name,layer),Form("Time vs. Residual %s-L%d;Time (ns);Residual (um)",name,layer),400,0,400,200,0,200);
      for(int wire=1;wire<=nwire;wire++){
        new TH1F( Form("%sl%dw%02d_TDC",name,layer,wire),Form("TDC %s-L%dW%02d;TDC ch.;Counts",name,layer,wire),2000,0,2000);
        new TH1F( Form("%sl%dw%02d_Time",name,layer,wire),Form("Time %s-L%dW%02d;Time (ns);Counts",name,layer,wire),400,0,400);
        new TH1F( Form("%sl%dw%02d_Length",name,layer,wire),Form("Length %s-L%dW%02d;Length (um);Counts",name,layer,wire),200,0,200);
        new TH2F( Form("%sl%dw%02d_TimevsLength",name,layer,wire),Form("Time vs. Length %s-L%dW%02d;Time (ns);Length (um)",name,layer,wire),400,0,400,200,0,200);
      }
    }
    for(int iconf=1;iconf<=LayConfigNum[idc];iconf++){
      new TH2F( Form("%s_HitPattern%c",name,LayConfig[idc][iconf-1]),Form("%s Hit pattern (%c%c' plane);Layer;Wire",name,LayConfig[idc][iconf-1],LayConfig[idc][iconf-1]),nlays+1,-0.5,nlays+0.5,maxwire+1,-0.5,maxwire+0.5);
    }
    new TH1F( Form("%s_HitEfficiency1%d",name,nlays),Form("Tracking efficiency %s (1-%d hit);Layer;Event",name,nlays),nlays+1,-0.5,nlays+0.5);
    new TH1F( Form("%s_HitEfficiency1%d",name,nlays-1),Form("Tracking efficiency %s (2-%d hit);Layer;Event",name,nlays-1),nlays+1,-0.5,nlays+0.5);
    new TH1F( Form("%s_HitEfficiency2%d",name,nlays),Form("Tracking efficiency %s (2-%d hit);Layer;Event",name,nlays),nlays+1,-0.5,nlays+0.5);

    /* For track */
    new TH1F( Form("%s_nTrack",name),Form("Number of tracks %s;Number of tracks;Counts",name),11,-0.5,10.5);
    new TH1F( Form("%s_ChiSquare",name),Form("#Chi^{2}/NDF %s;#Chi^{2}/NDF;Counts",name),100,0,100);
    new TH2F( Form("%s_TrackPositionXY",name),Form("Track position X vs. Y %s;X track position (cm);Y track position (cm)",name),20,-20,20,20,-20,20);
    new TH2F( Form("%s_TrackDirectionXY",name),Form("Track direction X vs. Y %s;X direction;Y direction",name),100,-0.50,0.50,100,-0.50,0.50);
    new TH2F( Form("%s_TrackPositionvsDirectionX",name),Form("Track position vs. direction X %s;X track position (cm);X direction",name),20,-20,20,100,-0.50,0.50);
    new TH2F( Form("%s_TrackPositionvsDirectionY",name),Form("Track position vs. direction Y %s;Y track position (cm);Y direction",name),20,-20,20,100,-0.50,0.50);
    TH1F* h1 = new TH1F( Form("%s_TrackEfficiency",name),Form("Tracking efficiency %s;;Event",name),6,-0.5,5.5);
    h1->GetXaxis()->SetBinLabel(1,"All");
    h1->GetXaxis()->SetBinLabel(2,"nTrack = 1");
    h1->GetXaxis()->SetBinLabel(3,"nTrack = 2");
    h1->GetXaxis()->SetBinLabel(4,"nTrack = 3");
    h1->GetXaxis()->SetBinLabel(5,"nTrack = 4");
    h1->GetXaxis()->SetBinLabel(6,"nTrack = 5");
  }

  return true;
}
