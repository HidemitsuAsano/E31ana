//////////////////////////////////////////
///           onana.c                
///           2010/11/06               
///    2010/11/08 modified to function as 
///               online analyzer
///               changed name from makeroot.c
///               H. Shi              
//////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TGClient.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>

#include "ana_common.h"

#define debug 0

void produce_hist(TString, TString);

int main( int argc, char *argv[])
{

  TString root_file(argv[1]);
  TString hist_file(argv[2]);

  Int_t utime;
  Int_t evid_daq;
  Int_t spill_gpio;
  Int_t evid_gpio;
  Int_t rpv100[8];
  Int_t v1290_stored_event[NV1290];
  Int_t v1290_gheader[NV1290];
  Int_t v1290_gtrailer[NV1290];
  Int_t v1290_theader[NV1290];
  Int_t v1290_ttrailer[NV1290];
  Int_t mtdc[NV1290][32][MAXHITS];

  printf("set branch \n");
  TFile *f_root  = new TFile( root_file, "READ");
  TTree *tree_mtdc = (TTree*)f_root->Get("tree_mtdc");

  tree_mtdc->SetBranchAddress("utime",&utime);
  tree_mtdc->SetBranchAddress("evid_daq",&evid_daq);
  tree_mtdc->SetBranchAddress("spill_gpio",&spill_gpio);
  tree_mtdc->SetBranchAddress("evid_gpio",&evid_gpio);
  tree_mtdc->SetBranchAddress("rpv100",rpv100);
  tree_mtdc->SetBranchAddress("v1290_stored_event",v1290_stored_event);
  tree_mtdc->SetBranchAddress("v1290_gheader",v1290_gheader);
  tree_mtdc->SetBranchAddress("v1290_gtrailer",v1290_gtrailer);
  tree_mtdc->SetBranchAddress("v1290_theader",v1290_theader);
  tree_mtdc->SetBranchAddress("v1290_ttrailer",v1290_ttrailer);
  tree_mtdc->SetBranchAddress("mtdc",mtdc);

  printf("create histo \n");
  TFile *f_hist  = new TFile( hist_file, "RECREATE");

  Int_t i, ii, jj, kk;

  TH1F *h_utime = new TH1F("h_utime","h_utime",1000, 0, 1000);
  TH1F *h_evid_daq = new TH1F("h_evid_daq","h_evid_daq",8000, 0, 8000);
  TH1F *h_spill_gpio = new TH1F("h_spill_gpio","h_spill_gpio",1000, 0, 1000);
  TH1F *h_evid_gpio = new TH1F("h_evid_gpio","h_evid_gpio",10000, 0, 10000);
  TH1F *h_error_flag = new TH1F("h_error_flag","h_error_flag",10, 0, 10);
  TH2F *h_gpio_rpv100 = new TH2F("h_gpio_rpv100","h_gpio_rpv100",5000, 0, 5000, 40000, 0, 40000);

  TH1F *h_rpv100[8];  
  for(ii=0; ii<8; ii++){
    h_rpv100[ii] = new TH1F(Form("h_rpv100_%d", ii),
			    Form("h_rpv100_%d", ii),
			    8000, 0, 8000);
  }

  TH1F *h_mtdc[NV1290][32];
  TH1F *h_mtdc_sub[NV1290][32];
  TH1F *h_nhits[NV1290][32];
  TH1F *h_v1290_stored_event[NV1290];
  TH1F *h_v1290_gheader[NV1290];
  TH1F *h_v1290_gtrailer[NV1290];
  TH1F *h_v1290_theader[NV1290];
  TH1F *h_v1290_ttrailer[NV1290];

  for(ii=0; ii<NV1290; ii++){
    h_v1290_stored_event[ii] = new TH1F(Form("h_v1290_stored_event_%d", ii),
				      Form("h_v1290_stored_event_%d", ii),
				      10, -1, 9);
    h_v1290_gheader[ii] = new TH1F(Form("h_v1290_gheader_%d", ii),
				   Form("h_v1290_gheader_%d", ii),
				   10000, 0, 10000);
    h_v1290_gtrailer[ii] = new TH1F(Form("h_v1290_gtrailer_%d", ii),
				    Form("h_v1290_gtrailer_%d", ii),
				    10000, 0, 10000);
    h_v1290_theader[ii] = new TH1F(Form("h_v1290_theader_%d", ii),
				   Form("h_v1290_theader_%d", ii),
				   10000, 0, 10000);
    h_v1290_ttrailer[ii] = new TH1F(Form("h_v1290_ttrailer_%d", ii),
				    Form("h_v1290_ttrailer_%d", ii),
				    10000, 0, 10000);
  
    for(jj=0; jj<32; jj++){
      h_mtdc[ii][jj] = new TH1F(Form("h_mtdc_%d_%d", ii, jj),
				Form("h_mtdc_%d_%d", ii, jj),
				40000, 0, 40000);
      h_mtdc_sub[ii][jj] = new TH1F(Form("h_mtdc_sub_%d_%d", ii, jj),
				    Form("h_mtdc_sub_%d_%d", ii, jj),
				40000, 0, 40000);
      h_nhits[ii][jj] = new TH1F(Form("h_nhits_%d_%d", ii, jj),
				 Form("h_nhits_%d_%d", ii, jj),
				 MAXHITS, 0, MAXHITS);
    }
    
  }

  TH1F *h_daq_eff = new TH1F("h_daq_eff","h_daq_eff",10000, 0, 10000);

  printf("start reading data \n");
  int nentries = tree_mtdc->GetEntries();
  printf("nentries = %d \n", nentries);

  int error_flag=0;
  int error_counter=0;
  for (i=0;i<nentries;i++) {

    tree_mtdc->GetEntry(i); 

    if (i%10000==0)
      printf("iev = %d \n", i);

    if( rpv100[1]!=(error_counter+1) ){
      error_flag = 1;
      error_counter = rpv100[1];
    }else{
      error_flag = 0;
      error_counter++;
    }
    h_error_flag->Fill(error_flag);

    h_utime->Fill(utime);
    h_evid_daq->Fill(evid_daq);
    h_spill_gpio->Fill(spill_gpio);
    h_evid_gpio->Fill(evid_gpio);
    h_gpio_rpv100->Fill(evid_gpio, rpv100[0]);

    for(ii=0; ii<8; ii++){
      h_rpv100[ii]->Fill(rpv100[ii]);
    }

    for(ii=0; ii<NV1290; ii++){
      h_v1290_stored_event[ii]->Fill(v1290_stored_event[ii]);
      h_v1290_gheader[ii]->Fill(v1290_gheader[ii]);
      h_v1290_gtrailer[ii]->Fill(v1290_gtrailer[ii]);
      h_v1290_theader[ii]->Fill(v1290_theader[ii]);
      h_v1290_ttrailer[ii]->Fill(v1290_ttrailer[ii]);
      for(jj=0; jj<32; jj++){
	int nhits=0;
	for(kk=0; kk<MAXHITS; kk++){
	  if( mtdc[ii][jj][kk] > 0){
	    nhits++;
	    int ref = 0;
	    if(ii == 0){
	      ref = mtdc[ii][15][0];
	    }else {
	      ref = mtdc[ii][31][0];
	    }
	    int tdc = mtdc[ii][jj][kk];
	    h_mtdc[ii][jj]->Fill(tdc);
	    h_mtdc_sub[ii][jj]->Fill(ref-tdc);
	  }
	}
	h_nhits[ii][jj]->Fill(nhits);
      }
    }
  }

  f_hist->Write();
  f_hist->Close();
  f_root->Close();

  return 0;
}
