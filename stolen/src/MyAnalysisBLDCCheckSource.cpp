// MyAnalysisBLDCCheckSource.cpp

#include "MyAnalysisBLDCCheckSource.h"

MyAnalysisBLDCCheckSource::MyAnalysisBLDCCheckSource(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	Clear();
}

MyAnalysisBLDCCheckSource::~MyAnalysisBLDCCheckSource()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisBLDCCheckSource::Clear()
{
	return;
}

bool MyAnalysisBLDCCheckSource::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
	rtFile->cd();
	DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  //if(!header->IsTrig(Trig_Beam)){
  //  return true;
  //}
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* ===================================== */
  /* ========== Event selection ========== */
  //double T0Timing = -9999.; 
  //int T0Segment = -1;
  //int MulT0 = 0;
  //for(int i=0; i<blMan->nT0(); i++){
  //  if(blMan->T0(i)->CheckRange()){
  //    MulT0++;
  //    T0Timing = blMan->T0(i)->ctmean();
  //    T0Segment = blMan->T0(i)->seg();
  //  }
  //}
  //if(MulT0!=1){
  //  return true;
  //}
  //if(T0Timing<-1.0||1.0<T0Timing){
  //  return true;
  //}
  //double BHDTiming = -9999.; 
  //int BHDSegment = -1;
  //int MulBHD = 0;
  //double tof = -9999;
  //for(int i=0; i<blMan->nBHD(); i++){
  //  if(blMan->BHD(i)->CheckRange()){
  //    MulBHD++;
  //    BHDTiming = blMan->BHD(i)->ctmean();
  //    BHDSegment = blMan->BHD(i)->seg();
  //    tof = BHDTiming - T0Timing;
  //  }
  //  if(25<tof&&tof<27) break;
  //}
  //if(tof<25||27<tof){
  //  return true;
  //}
  bltrackMan->DoTracking(blMan,conf,true,true);
  //if(bltrackMan->ntrackBLC1()!=1||bltrackMan->ntrackBLC2()!=1){
  //  return true;
  //}
  /* ========== Event selection ========== */
  /* ===================================== */

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


  /* Basic data filling */
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int layer=1;layer<=nlays;layer++){
      int mult = blMan->nBLDC(cid,layer);
      if(mult>0){
        FillHist(Form("%sl%d_Multiplicity",name,layer),mult);
      }
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

  /* Track data filling */
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    FillHist(Form("%s_nTrack",name),bltrackMan->ntrackBLDC(cid),1);
    if(bltrackMan->ntrackBLDC(cid)!=1) continue;
    for(int itra=0; itra<bltrackMan->ntrackBLDC(cid); itra++){
      LocalTrack* track = bltrackMan->trackBLDC(cid,itra);
      double tracktime = track->GetTrackTime();
      if(tracktime<-10||10<tracktime) continue;
      double chi2 = track -> chi2all();
      FillHist(Form("%s_TrackTime",name),tracktime,1);
      FillHist(Form("%s_ChiSquare",name),chi2,1);
      for(int ihit=0; ihit<track->nhit(); ihit++){
        ChamberLikeHit* hit = track->hit(ihit);
        int layer = hit->layer();
        int wire = hit->wire();
        int tdc = hit->tdc();
        double dt = hit->dt();
        double resl = hit->resl();
        FillHist(Form("%sl%d_HitPattern_inTrack",name,layer),wire,1);
        FillHist(Form("%sl%d_TDC_inTrack",name,layer),tdc,1);
        FillHist(Form("%sl%d_Time_inTrack",name,layer),dt,1);
        FillHist(Form("%sl%d_Residual_inTrack",name,layer),resl,1);
        FillHist(Form("%sl%d_TimevsResidual_inTrack",name,layer),dt,resl,1);
        FillHist(Form("%sl%dw%02d_TDC_inTrack",name,layer,wire),tdc,1);
        FillHist(Form("%sl%dw%02d_Time_inTrack",name,layer,wire),dt,1);
        FillHist(Form("%sl%dw%02d_Residual_inTrack",name,layer,wire),resl,1);
        FillHist(Form("%sl%dw%02d_TimevsResidual_inTrack",name,layer,wire),dt,resl,1);
        FillHist(Form("%s_HitPattern%c_inTrack",name,LayConfig[idc][layer-1]),layer,wire,1);
      }
      for(int ixy=0; ixy<2; ixy++){
        for(int iclust=0; iclust<track->ncluster(ixy); iclust++){
          BLDCCluster* cluster = track->cluster(ixy,iclust);
          int layer = 0;
          if(cluster->nhit()!=2) continue;
          ChamberLikeHit* hit = cluster->hit(0);
          layer = hit->layer();
          double timediff = cluster->GetTimeSub();
          double timemean = cluster->GetTimeMean();
          double ctimemean = cluster->GetCTime();
          FillHist(Form("%sl%d_TimeDifferencevsTimeMean_inTrack",name,layer),timediff,timemean,1);
          FillHist(Form("%sl%d_TimeDifferencevsCTimeMean_inTrack",name,layer),timediff,ctimemean,1);
        }
      }
    }
  }


  return true;

}

bool MyAnalysisBLDCCheckSource::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisBLDCCheckSource::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisBLDCCheckSource::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisBLDCCheckSource::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisBLDCCheckSource::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisBLDCCheckSource::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaBLDCCheckSource");
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
      new TH1F( Form("%sl%d_TDC_inTrack",name,layer),Form("TDC %s-L%d;TDC ch.;Counts",name,layer),1000,0,2000);
      new TH1F( Form("%sl%d_Time_inTrack",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),500,-100,400);
      new TH1F( Form("%sl%d_Residual_inTrack",name,layer),Form("Residual %s-L%d;Residual (cm);Counts",name,layer),500,-0.5,0.5);
      new TH2F( Form("%sl%d_TimevsResidual_inTrack",name,layer),Form("Time vs. Residual %s-L%d;Time (ns);Residual (cm)",name,layer),500,-100,400,500,-0.5,0.5);
      if(layer%2==1){
        new TH2F( Form("%sl%d_TimeDifferencevsTimeMean_inTrack",name,layer),Form("Time Difference vs. Time Mean %s-L%d and L%d;Time Difference (ns);Time Mean (ns)",name,layer,layer+1),500,-250,250,500,-100,400);
        new TH2F( Form("%sl%d_TimeDifferencevsCTimeMean_inTrack",name,layer),Form("Time Difference vs. Corrected Time Mean %s-L%d and L%d;Time Difference (ns);Time Mean (ns)",name,layer,layer+1),500,-250,250,500,-100,400);
      }
      for(int wire=1;wire<=nwire;wire++){
        //new TH1F( Form("%sl%dw%02d_TDC_inTrack",name,layer,wire),Form("TDC %s-L%dW%02d;TDC ch.;Counts",name,layer,wire),1000,0,2000);
        new TH1F( Form("%sl%dw%02d_Time_inTrack",name,layer,wire),Form("Time %s-L%dW%02d;Time (ns);Counts",name,layer,wire),500,-100,400);
        new TH1F( Form("%sl%dw%02d_Residual_inTrack",name,layer,wire),Form("Residual %s-L%dW%02d;Residual (um);Counts",name,layer,wire),200,-0.2,0.2);
        new TH2F( Form("%sl%dw%02d_TimevsResidual_inTrack",name,layer,wire),Form("Time vs. Residual %s-L%dW%02d;Time (ns);Residual (um)",name,layer,wire),500,-100,400,200,-0.2,0.2);
      }

      new TH1F( Form("%sl%d_Multiplicity",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_HitPattern",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_TDC",name,layer),Form("TDC %s-L%d;TDC ch.;Counts",name,layer),1000,0,2000);
      new TH1F( Form("%sl%d_Time",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),500,-100,400);
      new TH1F( Form("%sl%d_Length",name,layer),Form("Length %s-L%d;Length (cm);Counts",name,layer),250,0.0,0.5);
      new TH2F( Form("%sl%d_TimevsLength",name,layer),Form("Time vs. Length %s-L%d;Time (ns);Length (cm)",name,layer),500,-100,400,250,0.0,0.5);
      for(int wire=1;wire<=nwire;wire++){
        //new TH1F( Form("%sl%dw%02d_TDC",name,layer,wire),Form("TDC %s-L%dW%02d;TDC ch.;Counts",name,layer,wire),1000,0,2000);
        new TH1F( Form("%sl%dw%02d_Time",name,layer,wire),Form("Time %s-L%dW%02d;Time (ns);Counts",name,layer,wire),500,-100,400);
        //new TH1F( Form("%sl%dw%02d_Length",name,layer,wire),Form("Length %s-L%dW%02d;Length (cm);Counts",name,layer,wire),250,0.0,0.5);
        //new TH2F( Form("%sl%dw%02d_TimevsLength",name,layer,wire),Form("Time vs. Length %s-L%dW%02d;Time (ns);Length (cm)",name,layer,wire),500,-100,400,250,0.0,0.5);
      }
    }

    /* For track */
    new TH1F( Form("%s_nTrack",name),Form("Number of tracks %s;Number of tracks;Counts",name),11,-0.5,10.5);
    new TH1F( Form("%s_ChiSquare",name),Form("#Chi^{2}/NDF %s;#Chi^{2}/NDF;Counts",name),1000,0,100);
    new TH1F( Form("%s_TrackTime",name),Form("Track time %s;Track time (ns);Counts",name),1000,-100,400);
  }

  return true;
}
