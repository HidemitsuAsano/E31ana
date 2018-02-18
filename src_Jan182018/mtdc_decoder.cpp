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

static HeaderStruct      header_s; 
static Rpv100Struct      rpv100_s;
static GpioStruct        gpioRM_s;
static V1290Struct       v1290_s[NV1290];
static EventStruct       savs;

void init_tree( TString, TTree*, EventStruct* ev);
void produce_root(FILE *, TString);
void print_hist(TString);

int main( int argc, char *argv[])
{
    FILE *fp=0;
    TString fname;

    // Raw file open ---
    if( argc==2 || argc==3 || argc==4){
        fp = fopen( argv[1], "r");
        if( fp == NULL ){
            printf("### error: cannot open the file %s...\n", argv[1]);
            exit(1);
        }
    }else{
      printf("### error: usage: ./onana file_name");
      exit(1);
    }

    TString rootfile(argv[2]);
    produce_root(fp, rootfile);
    
    fclose(fp);
    return 0;
}

void produce_root(FILE *fp, TString rootfile)
{
  TFile *f  = new TFile( rootfile, "RECREATE" );
  TTree *tree_mtdc= new TTree("tree_mtdc", "VME TKO daq evid sync");
  init_tree( rootfile, tree_mtdc, &savs );
  printf("Initialized ... \n");
  
  int result;
  int i, j, n, m, k;
  
  while(1){

    int gheader[NV1290]={0};
    int gtrailer[NV1290]={0};
    int theader[NV1290]={0};
    int ttrailer[NV1290]={0};
    int stored_event[NV1290]={0};
    int ch[NV1290]={0};
    int nhit[NV1290][32]={{0}};
    int tdc[NV1290][32][MAXHITS]={{{-999}}};
    int scaler[8]={-9};

    result = fread(&header_s, sizeof(HeaderStruct), 1, fp);

    if( result != 1 ){
      if(feof(fp)){
	printf("Reached end of file. Polling...\n");
      }break;
    }
    
    // read rpv100 data
    if(fread(&rpv100_s, sizeof(Rpv100Struct), 1, fp) != 1)break;
    for (i=0; i<8; i++ ){
      char tmp[10];
      int itmp;
      sprintf(tmp,"%02x%02x%02x%02x ",
	      rpv100_s.count[i]>>8&0xff,
	      rpv100_s.count[i]&0xff,
	      rpv100_s.count[i]>>24&0xff,
	      rpv100_s.count[i]>>16&0xff);
      sscanf(tmp,"%d",&itmp);
      //      if(i==0)
      //	printf("itmp = %d \n", itmp);
      scaler[i]=itmp;
    }
    
    // read v1290 data
    for( i=0; i<NV1290; i++){
      if(fread(&v1290_s[i].len1,         sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].tag1,         sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].num,          sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].blank1,       sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].len2,         sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].stored_event, sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].i,            sizeof(int), 1, fp) != 1)break;
      if(fread(&v1290_s[i].tag2,         sizeof(int), 1, fp) != 1)break;

      int buf_length;
      buf_length = v1290_s[i].len1-8;
      if(buf_length == 0) continue;
      if(fread(v1290_s[i].data_buf,   sizeof(int)*buf_length, 1, fp) !=1)break;
      
      for(j=0;j<v1290_s[i].len1-8;j++){ 
	if( (0x1f&(v1290_s[i].data_buf[j]>>27)) == 0x8 ){
	  gheader[i]++;
	}
	if( (0x1f&(v1290_s[i].data_buf[j]>>27)) == 0x1 ){
	  theader[i]++;
	}
	if( (0x1f&(v1290_s[i].data_buf[j]>>27)) == 0x0 ){
	  ch[i] = 0x1f&(v1290_s[i].data_buf[j]>>21);
	  tdc[i][ch[i]][nhit[i][ch[i]]] = 0x1fffff&(v1290_s[i].data_buf[j]);
	  nhit[i][ch[i]]++;
	  if(nhit[i][ch[i]]>10){
	    printf("in event %d  module %d, ch %d too many hits! #nhits = %d\n", header_s.iev, i, ch[i], nhit[i][ch[i]]);
	    //	    printf("force nhits to be 10!!\n");
	    //	    nhit[i][ch[i]]=10;
	    //getchar();
	  }
	}
	if( (0x1f&(v1290_s[i].data_buf[j]>>27)) == 0x3 ){
	  ttrailer[i]++;
	}
	if( (0x1f&(v1290_s[i].data_buf[j]>>27)) == 0x10 ){
	  gtrailer[i]++;

	}

	if( (0x1f&(v1290_s[i].data_buf[j]>>27)) == 0x4 ) printf("TDC error \n");

      }

    }

    for( i=0; i<NV1290; i++){
      if( (gheader[i]==0)|(gtrailer[i]==0) ){
	printf("event %d from module %d has no global header/trailer! \n", header_s.iev, i);
	getchar();
      }
    }

    // read gpioRM data
    if(fread(&gpioRM_s, sizeof(GpioStruct), 1, fp) != 1)break;

    // pass data to EventStruct
    savs.utime  = header_s.tepoch;
    savs.evid_daq  = header_s.iev;
    savs.spill_gpio  = gpioRM_s.spill & 0xff;
    savs.evid_gpio   = gpioRM_s.evid & 0xfff;
    for(n=0; n<8; n++){
      savs.rpv100[n]=scaler[n];
    }

    for(i=0; i<2; i++){
      savs.v1290_stored_event[i]=stored_event[i];
      savs.v1290_gheader[i]=gheader[i];
      savs.v1290_gtrailer[i]=gtrailer[i];
      savs.v1290_theader[i]=theader[i];
      savs.v1290_ttrailer[i]=ttrailer[i];
      
      for(m=0; m<32; m++){
	for(k=0; k<MAXHITS; k++){
	  if(tdc[i][m][k]>0){
	    savs.mtdc[i][m][k]=tdc[i][m][k];
	  }	
	}
      }
      
      savs.v1290_stored_event[i]=stored_event[i];
      savs.v1290_gheader[i]=gheader[i];
      savs.v1290_gtrailer[i]=gtrailer[i];
      savs.v1290_theader[i]=theader[i];
      savs.v1290_ttrailer[i]=ttrailer[i];
    }

    tree_mtdc->Fill();

    //    printf("%d\t%d\t%d\t%d\n",gpioRM_s.spill & 0xff,gpioRM_s.evid & 0xfff, tdc_tmp_cha-tdc_tmp_chb,tdc_tmp_cha-tdc_tmp_chc);

#if debug
    printf("################################################################## \n");
    printf("output of event info for debug \n");
    /*
    printf("sizeof(int) = %d \n", sizeof(int));
    printf("sizeof(unsigned int) = %d \n", sizeof(unsigned int));
    printf("sizeof(long) = %d \n", sizeof(long));
    printf("sizeof(HeaderStruct) = %d \n", sizeof(HeaderStruct));
    printf("size from header = %d \n", header_s.size);
    printf("tepoch from header = %d \n", header_s.tepoch);
    printf("iev from header = %d \n", header_s.iev);
    */
    printf("iev from header = %d \n", header_s.iev);
    printf("gpio event = %d \n", gpioRM_s.evid & 0xfff);
    printf("gpio spill = %d \n", gpioRM_s.spill & 0xff);

    for( i=0; i<NV1290; i++){
      printf("------v1290 output from module %d------ \n", i);
      printf("len1 = %d \n", v1290_s[i].len1);
      printf("tag2 = 0x%04x \n", v1290_s[i].tag2);
      printf("stored event = %d \n", v1290_s[i].stored_event);
    
      for( j=0; j<v1290_s[i].len1-8; j++ ){
	printf("data_buf = 0x%04x \n", v1290_s[i].data_buf[j]);
      }

      printf("global header: %d \n", gheader[i]);
      printf("TDC header: %d \n", theader[i]);

      printf("TDC measurement \n");
      for(m=0; m<32; m++){
	for(k=0; k<MAXHITS; k++){
	  if(tdc[i][m][k]>0){
	    printf("channel = %d, #hit = %d, tdc = %d \n", m, k, tdc[i][m][k]);
	  }	
	}
      }
      printf("TDC trailer: %d \n", ttrailer[i]);
      printf("global trailer: %d \n", gtrailer[i]);
    }

#endif 

    //clear data
    header_s.tepoch = 0;
    header_s.iev = 0;
    gpioRM_s.spill = 0;
    gpioRM_s.evid = 0;
    for(n=0; n<8; n++){
      scaler[n] = 0;
    }
    for(i=0; i<2; i++){
      int jj, kk;
      gheader[i]=0;
      gtrailer[i]=0;
      theader[i]=0;
      ttrailer[i]=0;
      stored_event[i]=0;	  
      ch[i] = 0;
      for(jj=0; jj<32; jj++){
	nhit[i][jj] = 0;
	for(kk=0; kk<MAXHITS; kk++){
	  tdc[i][jj][kk] = -999;
	  savs.mtdc[i][jj][kk]=-999;
	}
      }
    }

  }
  
  printf("write & close file \n");
  f->Write();
  f->Close();
  
}

void init_tree( TString root_file,  TTree *tree_mtdc,  EventStruct *savs )
{
   tree_mtdc->Branch("utime",  &savs->utime, "utime/I");
   //   printf("iev from header = %d \n", header_s.iev);
   tree_mtdc->Branch("evid_daq",   &header_s.iev,  "evid_daq/I");
   //gpio data 
   tree_mtdc->Branch("spill_gpio",  &savs->spill_gpio, "spill_gpio/I");
   tree_mtdc->Branch("evid_gpio",   &savs->evid_gpio,  "evid_gpio/I");
   //rpv100 data 
   tree_mtdc->Branch("rpv100",    savs->rpv100,  "rpv100[8]/I");
   //mtdc data 
   tree_mtdc->Branch("v1290_stored_event", savs->v1290_stored_event, 
	      Form("v1290_stored_event[%d]/I", NV1290));
   tree_mtdc->Branch("v1290_gheader", savs->v1290_gheader,
	      Form("v1290_gheader[%d]/I", NV1290));
   tree_mtdc->Branch("v1290_gtrailer", savs->v1290_gtrailer,
	      Form("v1290_gtrailer[%d]/I", NV1290));
   tree_mtdc->Branch("v1290_theader", savs->v1290_theader,
	      Form("v1290_theader[%d]/I", NV1290));
   tree_mtdc->Branch("v1290_ttrailer", savs->v1290_gtrailer,
	      Form("v1290_ttrailer[%d]/I", NV1290));
   tree_mtdc->Branch("mtdc",    savs->mtdc,  
		     Form("mtdc[%d][32][%d]/I", NV1290, MAXHITS));

}

