#include <iostream>
#include <fstream>
#include <iomanip.h>
#include "TMinuit.h"

#define DISPLAY 0
#define DISPLAY2 0

class GlobalVariables;
class ConfMan;

#define CDH 0
#define T0 0
#define NC 1
#define DRAW 1
#define DEBUG 1
void CalcOffsetNC(std::string InFileName);
void CalcOffsetNC(){
  /*
  CalcOffsetNC("SlewNC34.param");
  CalcOffsetNC("SlewNC61.param");
  CalcOffsetNC("SlewNC67.param");
  CalcOffsetNC("SlewNC73.param");
  CalcOffsetNC("SlewNC82.param");
  CalcOffsetNC("SlewNC90.param");
  CalcOffsetNC("SlewNC103.param");
  CalcOffsetNC("SlewNC104.param");
  CalcOffsetNC("SlewNC105.param");
  CalcOffsetNC("SlewNC118.param");
  CalcOffsetNC("SlewNC124.param");
  CalcOffsetNC("SlewNC136.param");
  CalcOffsetNC("SlewNC142.param");
  CalcOffsetNC("SlewNC148.param");
  CalcOffsetNC("SlewNC165.param");
  CalcOffsetNC("SlewNC176.param");
  */
  CalcOffsetNC("SlewPC170.param");
  CalcOffsetNC("SlewPC176.param");
  CalcOffsetNC("SlewPC182.param");
  CalcOffsetNC("SlewPC188.param");
}

void TDCOffset(TH1F *h1,ConfMan* conf, int cid,int seg,int iud){
  const int ncounter=2;
  TString name[ncounter]={"NC","PC"};
  double tof[ncounter]={0,0};
  TString ud[2]={"U","D"};

  int cr,sl,ch;  
  TH1F* h1;
  //  h1->Rebin(10);
  double tmpm,mean;
  tmpm=h1->GetBinCenter(h1->GetMaximumBin());
  h1->Fit("gaus","LI0q","",tmpm-5,tmpm+5);
  mean=h1->GetFunction("gaus")->GetParameter(1);
  h1->Draw();
  h1->GetFunction("gaus")->Draw("same");
  
  double gain,ped,tmpped;
  conf->GetCounterMapManager()->GetCNA(cid,seg,1,iud,cr,sl,ch);
  conf->GetGainMapManager()->GetParam( cr,sl,ch, gain, tmpped );
  ped = mean*gain;
  conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, ped );
#if DEBUG
  std::cout<<ud[iud]<<seg<<"\t"<<ped<<std::endl;
#endif  
  int ymax=  h1->GetMaximum();
  TLine line;
  line.SetLineColor(2);  
  line.DrawLine( mean, 1, mean, ymax );
}

void ADCOffset(TH1F *h1,ConfMan* conf, int cid,int seg,int iud){
  const int ncounter=2;
  TString name[ncounter]={"NC","PC"};
  double tof[ncounter]={0,0};
  TString ud[2]={"U","D"};

  int cr,sl,ch;  
  TH1F* h1;
  double tmpm,mean;
  tmpm=h1->GetBinCenter(h1->GetMaximumBin());
  h1->Fit("gaus","LI0q","",tmpm-5,tmpm+5);
  mean=h1->GetFunction("gaus")->GetParameter(1);
  h1->GetXaxis()->SetRangeUser(mean-10,mean+10);
#if DRAW
  h1->Draw();
  h1->GetFunction("gaus")->Draw("same");
#endif
  double gain,ped,tmpped;
  conf->GetCounterMapManager()->GetCNA(cid,seg,0,iud,cr,sl,ch);
  conf->GetGainMapManager()->GetParam( cr,sl,ch, gain, tmpped );
  conf->GetGainMapManager()->SetParam( cr,sl,ch, gain, mean );
#if DEBUG
  std::cout<<ud[iud]<<seg<<"\t"<<mean<<std::endl;
#endif
#if DRAW  
  int ymax=  h1->GetMaximum();
  TLine line;
  line.SetLineColor(2);  
  line.DrawLine( mean, 1, mean, ymax );
#endif
}


void CalcOffsetNC(std::string InFileName)
{

  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");
  //########parameter file

  gStyle->SetOptFit();
#define MAXCHAR 256   
  ifstream fp;
  //  string InFileName="SlewNC.param";
  char str[MAXCHAR],str1[MAXCHAR];


  int runnum = 83;
  int Layer[2]={6,7};
  int sseg[2][4]={9,10,11,12,9,10,11,12};
  char* CName[2]={"NC","NC"};
  int CID[2]={CID_NC,CID_NC};
  int ita_ini=15;
  int ita_fin=16;
  int nset=2;
  string conffile;
  string slewfile;
  int writeparam=0;
  std::cout << InFileName <<endl;
  if( fp.open(InFileName.c_str() )==0 )
    {
      std::cerr << " File open fail. [" << InFileName << "]" << std::endl;
      return;
    }

  while( fp.getline(str,MAXCHAR) )
    {
      if( str[0]=='#' )continue;
      if( sscanf(str,"Run: %d", &runnum)==1 ) 
	std::cout<<"#Run= "<<runnum<<std::endl;
      else if( sscanf(str,"NumberOfSet: %d", &nset)==1 ) 
      	std::cout<<"nset= "<<nset<<std::endl;
      else if( sscanf(str,"Layer: %d %d", &Layer[0],&Layer[1])==2 ) 
      	std::cout<<"Layer[0]= "<<Layer[0]<<" Layer[1]= "<<Layer[1]<<std::endl;
      else if( sscanf(str,"segment[0]: %d %d %d %d", &sseg[0][0], &sseg[0][1], &sseg[0][2], &sseg[0][3])==4 )
      	std::cout<<"sseg[0]= "<<sseg[0][0]<<" "<<sseg[0][1]<<" "<<sseg[0][2]<<" "<<sseg[0][3]<<std::endl;
      else if( sscanf(str,"segment[1]: %d %d %d %d", &sseg[1][0], &sseg[1][1], &sseg[1][2], &sseg[1][3])==4 )
      	std::cout<<"sseg[1]= "<<sseg[1][0]<<" "<<sseg[1][1]<<" "<<sseg[1][2]<<" "<<sseg[1][3]<<std::endl;
      else if( sscanf(str,"itanum: %d %d", &ita_ini, &ita_fin)==2 )
	std::cout<<"itanum= "<<ita_ini<<" "<<ita_fin<<std::endl;
      else if( sscanf(str,"Conffile: %s", str1)==1 )
	{
	  conffile=str1;	
	  std::cout<<"conffile= "<<conffile<<std::endl;
	}
      else if( sscanf(str,"Slewfile: %s", str1)==1 )
	{
	  slewfile=str1;	
	  std::cout<<"slewfile= "<<slewfile<<std::endl;
	}
      else if( sscanf(str,"WriteParam: %d", &writeparam)==1 )
	std::cout<<"writeparam= "<<writeparam<<std::endl;
    }

  for(int i=0;i<2;i++)
    if(Layer[i]<0){
      CID[i]=CID_PC; 
      CName[i]="PC";
    }
  
  int rseg[2][4];
  for(int iset=0;iset<nset;iset++){
    if(CID[iset]==CID_NC)
      for(int iseg=0;iseg<4;iseg++){
	rseg[iset][iseg]=(Layer[iset]-1)*16+sseg[iset][iseg];
	std::cout<<"rseg "<<iset <<" "<<iseg<<" "<<CName[iset]<<" "<<rseg[iset][iseg]<<std::endl;;
      }
    else
      for(int iseg=0;iseg<4;iseg++){
	rseg[iset][iseg]=sseg[iset][iseg];
	std::cout<<"rseg "<<iset <<" "<<iseg<<" "<<CName[iset]<<" "<<rseg[iset][iseg]<<std::endl;;
      }
  }

  ConfMan* conf=new ConfMan(conffile.data());
  conf->Initialize();
  GainMapMan* gainmap=conf->GetGainMapManager();
  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TFile *f;
  f = new TFile( Form( "./root/NCcalib/NCr000%d.root",runnum) );
  TH1F *h1;
  TH2F *h2;
  TProfile *prof;
  TCanvas *c1;
  c1=new TCanvas();
  c1->Divide(4,4);
  for(int iset=0;iset<nset;iset++){
    for(int iseg1=0;iseg1<4;iseg1++){
      int seg1=rseg[iset][iseg1];
      for(int ud=0;ud<2;ud++){      
	c1->cd(4*iset+2*iseg1+ud+1);
	if(ud==0)
	  h1 = (TH1F*)f->Get( Form("T%sU%d",CName[iset],seg1) );
	else 
	  h1 = (TH1F*)f->Get( Form("T%sD%d",CName[iset],seg1) );
	TDCOffset(h1,conf,CID[iset],seg1,ud);
      }
    }
  }
  c1->Print(Form("NC_offset_run%d.ps(",runnum));
  c1=new TCanvas();
  c1->Divide(4,4);
  for(int iset=0;iset<nset;iset++){
    for(int iseg1=0;iseg1<4;iseg1++){
      int seg1=rseg[iset][iseg1];
      for(int ud=0;ud<2;ud++){      
	c1->cd(4*iset+2*iseg1+ud+1);
	if(ud==0)
	  h1 = (TH1F*)f->Get( Form("AwoT%sU%d",CName[iset],seg1) );
	else 
	  h1 = (TH1F*)f->Get( Form("AwoT%sD%d",CName[iset],seg1) );
	ADCOffset(h1,conf,CID[iset],seg1,ud);
      }
    }
  }
  c1->Print(Form("NC_offset_run%d.ps)",runnum));

  ofstream ofs("NCGain.param");
  conf->GetGainMapManager()->PrintMapBL(ofs); // print map for BeamLine
  ofs.close();
}
