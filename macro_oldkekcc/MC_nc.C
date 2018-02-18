const double ADC_TH_CVC = 0.5; // MeV
const double ADC_TH_TOFstop = 1.0; // MeV
const double ADC_TH_NC = 2.0; // MeV

class MCData;
class CDSHitMan;
class CDSTrackingMan;
class BeamLineHitMan;
class BeamLineTrackMan;
class ConfMan;

#define DEBUG 0

void MC_nc()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(0);
  gROOT->cd();

  /*** load library ***/
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("./lib/libAll.so");
  //#######Need knucl lib file!! Set your dir!! ###
  gSystem->Load("~/work/ana/geant/knucl3/libKnuclRootData.so");
  //#######################################

  /*** assign input file & call tree ***/
  TFile *f = new TFile( "tmp.root" );
  TTree *evtree = (TTree*)f->Get("EventTree");

  /*** conf file for new parameters ***/
  ConfMan *conf = new ConfMan("conf/Run43/analyzer.conf");
  conf->Initialize();

  /*** declaration of classes ***/
  MCData* mcData = 0;
  evtree->SetBranchAddress( "MCData", &mcData );
  CDSHitMan* cdsMan = 0;
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  CDSTrackingMan* cdsTrackMan = new CDSTrackingMan();
  evtree->SetBranchAddress( "CDSTrackingMan", &cdsTrackMan );
  BeamLineHitMan* blMan = 0;
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );
  BeamLineTrackMan* blTrackMan = new BeamLineTrackMan();
  evtree->SetBranchAddress( "BeamLineTrackMan", &blTrackMan );

  /*** define histrgrams ***/
  TH1F* adc[2][2][3]; // [time,layer][w/o CVC, w/ CVC][all, gamma, neutron]
  TH1F* tdc[2][2][3]; // [time,layer][w/o CVC, w/ CVC][all, gamma, neutron]
  TH2F* hit[2][2][3]; // [time,layer][w/o CVC, w/ CVC][all, gamma, neutron]
  char com[128];
  for(int i = 0; i<2; i++ ){
    for(int j = 0; j<2; j++ ){
      for(int k = 0; k<3; k++ ){
	sprintf(com, "adc%d_%d_%d", i, j, k);
	adc[i][j][k] = new TH1F(com, com, 100, 0, 100);
	sprintf(com, "tdc%d_%d_%d", i, j, k);
	tdc[i][j][k] = new TH1F(com, com, 2000, 0, 200);
	sprintf(com, "hit%d_%d_%d", i, j, k);
	hit[i][j][k] = new TH2F(com, com, 16, 1, 17, 7, 1, 8);
      }
    }
  }

  /*                */
  /* event analysis */
  /*                */
  int nev = evtree->GetEntries();
  int fev = 0;
  std::cerr<<"# of events = "<<nev<<std::endl;
  for(int iev = 0; iev<nev; iev++ ){
    evtree->GetEvent(iev);
    if( iev%1000==0 )
      std::cout << " Event : " << iev << std::endl;
    //###########################//
    //### check detector hits ###//
    //###########################//

#if DEBUG
    std::cerr<<"###### event num  = "<<iev<<std::endl;
    std::cerr<<"# of T0 hits      = "<<blMan->nT0()<<std::endl;
    std::cerr<<"# of CDH hits     = "<<cdsMan->nCDH()<<std::endl;
    std::cerr<<"# of BLC1 hits    = "<<blMan->nBLC1()<<std::endl;
    std::cerr<<"# of BLC2 hits    = "<<blMan->nBLC2()<<std::endl;
    std::cerr<<"# of CDC hits     = "<<cdsMan->nCDC()<<std::endl;
    std::cerr<<"# of NC hits      = "<<blMan->nNC()<<std::endl;
    std::cerr<<"# of TOFstop hits = "<<blMan->nTOF()<<std::endl;
    std::cerr<<"# of PC hits      = "<<blMan->nPC()<<std::endl;
    std::cerr<<"# of CVC hits     = "<<blMan->nCV()<<std::endl;
#endif

    double th_time  = 100000000000;
    double th_layer = 10;
    int tag_time  = -1;
    int tag_layer = -1;
    //if( !blMan->nCV() && !blMan->nTOF() && blMan->nNC() ){ // forward nuetral paticle event
    if( !blMan->nTOF() && blMan->nNC() ){ // forward nuetral paticle event
      int seg, layer, column;
      for(int i=0; i<blMan->nNC(); i++){
	HodoscopeLikeHit *hod = blMan->NC(i);
	seg    = hod->seg()-1; // 0-111
	layer  = seg/16+1; // 1-7
	column = seg%16+1; // 1-16
#if 0
	std::cerr<<"evnum="<<iev<<", NC : "<<hod->hid()<<" ("<<layer<<","<<column<<") "
		 <<hod->ctmean()<<" "<<hod->emean()<<std::endl;
#endif
	if ( hod->ctmean()<th_time ){
	  th_time = hod->ctmean();
	  tag_time = i;
	}
	if ( layer<th_layer) {
	  th_time = hod->ctmean();
	  th_layer = layer;	  
	  tag_layer = i;
	}
	else if ( layer==th_layer && hod->ctmean()<th_time ) {
	  th_time = hod->ctmean();
	  tag_layer = i;
	}
      }
      HodoscopeLikeHit *hod2[2];
      hod2[0] = blMan->NC(tag_time);
      hod2[1] = blMan->NC(tag_layer);
      int seg2[2], layer2[2], column2[2];

      // [time,layer][w/o CVC, w/ CVC][all, gamma, neutron]
      for(int i=0; i<2; i++){
	seg2[i]    = hod2[i]->seg()-1; // 0-111
	layer2[i]  = seg2[i]/16+1; // 1-7
	column2[i] = seg2[i]%16+1; // 1-16
	adc[i][0][0]->Fill(hod2[i]->emean());
	tdc[i][0][0]->Fill(hod2[i]->ctmean());
	hit[i][0][0]->Fill(column2[i], layer2[i]);
	if(hod2[i]->ctmean()<60){
	  adc[i][0][1]->Fill(hod2[i]->emean());
	  tdc[i][0][1]->Fill(hod2[i]->ctmean());
	  hit[i][0][1]->Fill(column2[i], layer2[i]);
	}else{
	  adc[i][0][2]->Fill(hod2[i]->emean());
	  tdc[i][0][2]->Fill(hod2[i]->ctmean());
	  hit[i][0][2]->Fill(column2[i], layer2[i]);
	}
      }
      if(!blMan->nCV()){
	for(int i=0; i<2; i++){
	  adc[i][1][0]->Fill(hod2[i]->emean());
	  tdc[i][1][0]->Fill(hod2[i]->ctmean());
	  hit[i][1][0]->Fill(column2[i], layer2[i]);
	  if(hod2[i]->ctmean()<60){
	    adc[i][1][1]->Fill(hod2[i]->emean());
	    tdc[i][1][1]->Fill(hod2[i]->ctmean());
	    hit[i][1][1]->Fill(column2[i], layer2[i]);
	  }else{
	    adc[i][1][2]->Fill(hod2[i]->emean());
	    tdc[i][1][2]->Fill(hod2[i]->ctmean());
	    hit[i][1][2]->Fill(column2[i], layer2[i]);
	  }
	}
      }

#if 0
      layer  = hod_time->seg()/16+1; // 1-7
      column = hod_time->seg()%16;   // 1-16
      std::cerr<<"fast timing hit: "<<hod_time->hid()<<" ("<<layer<<","<<column<<") "
	       <<hod_time->ctmean()<<" "<<hod_time->emean()<<std::endl;
      layer  = hod_layer->seg()/16+1; // 1-7
      column = hod_layer->seg()%16;   // 1-16
      std::cerr<<"first layer hit: "<<hod_layer->hid()<<" ("<<layer<<","<<column<<") "
	       <<hod_layer->ctmean()<<" "<<hod_layer->emean()<<std::endl;
      std::cerr<<std::endl;
#endif
      
    fev++;
    }

  } // for(int iev = 0; iev<nev; iev++ ){

  std::cerr<<"filled events = "<<fev<<std::endl;



  //%%%%%%%%%%%%%%%%%%%%%%//
  //%%% plot histogram %%%//
  //%%%%%%%%%%%%%%%%%%%%%%//

  gROOT->cd();

  // [time,layer][w/o CVC, w/ CVC][all, gamma, neutron]
  TCanvas *c1 = new TCanvas("c1", "");
  c1->Divide(2,2);
  c1->cd(1); adc[0][0][0]->Draw();
  adc[0][0][1]->SetLineColor(2); adc[0][0][1]->Draw("same");
  adc[0][0][2]->SetLineColor(4); adc[0][0][2]->Draw("same");
  c1->cd(2); tdc[0][0][0]->Draw();
  tdc[0][0][1]->SetLineColor(2); tdc[0][0][1]->Draw("same");
  tdc[0][0][2]->SetLineColor(4); tdc[0][0][2]->Draw("same");
  c1->cd(3); adc[1][0][0]->Draw();
  adc[1][0][1]->SetLineColor(2); adc[1][0][1]->Draw("same");
  adc[1][0][2]->SetLineColor(4); adc[1][0][2]->Draw("same");
  c1->cd(4); tdc[1][0][0]->Draw();
  tdc[1][0][1]->SetLineColor(2); tdc[1][0][1]->Draw("same");
  tdc[1][0][2]->SetLineColor(4); tdc[1][0][2]->Draw("same");
  c1->Print("tmp.pdf");

  TCanvas *c2 = new TCanvas("c2", "");
  c2->Divide(2,2);
  c2->cd(1); adc[0][1][0]->Draw();
  adc[0][1][1]->SetLineColor(2); adc[0][1][1]->Draw("same");
  adc[0][1][2]->SetLineColor(4); adc[0][1][2]->Draw("same");
  c2->cd(2); tdc[0][1][0]->Draw();
  tdc[0][1][1]->SetLineColor(2); tdc[0][1][1]->Draw("same");
  tdc[0][1][2]->SetLineColor(4); tdc[0][1][2]->Draw("same");
  c2->cd(3); adc[1][1][0]->Draw();
  adc[1][1][1]->SetLineColor(2); adc[1][1][1]->Draw("same");
  adc[1][1][2]->SetLineColor(4); adc[1][1][2]->Draw("same");
  c2->cd(4); tdc[1][1][0]->Draw();
  tdc[1][1][1]->SetLineColor(2); tdc[1][1][1]->Draw("same");
  tdc[1][1][2]->SetLineColor(4); tdc[1][1][2]->Draw("same");
  c2->Print("tmp2.pdf");

  gStyle->SetOptStat(0);

  TCanvas *c3 = new TCanvas("c3", "");
  c3->Divide(2,2);
  c3->cd(1); hit[0][1][1]->Draw("colz");
  c3->cd(2); hit[0][1][2]->Draw("colz");
  c3->cd(3); hit[1][1][1]->Draw("colz");
  c3->cd(4); hit[1][1][2]->Draw("colz");
  c3->Print("tmp3.pdf");

}

