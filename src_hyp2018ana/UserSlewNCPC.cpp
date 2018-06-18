#include <iostream>
#include <fstream>
#include <TApplication.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TRint.h>
#include <TLorentzVector.h>
#include  <TGraphErrors.h> 

#include "ConfMan.h"
#include "TKO.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "ScalerMan.h"
#include "HodoscopeLikeHit.h"
#include "CherenkovLikeHit.h"
#include "ChamberLikeHit.h"
#include "EventHeader.h"
#include "Display3D.h"

#define DEBUG 0
#define SlewType 1
#define SLEWING 1

int main( int argc, char **argv )
{
  if(argc!=2){
    std::cout<<"Please set \"SlewNC.param\""<<std::endl;
    return 1;
  }  
#define MAXCHAR 256   
  //  std::string InFileName="SlewNC.param";
  std::string InFileName=argv[1];
  std::cout << InFileName <<std::endl;
  ifstream fp(InFileName.c_str(),std::ios::in);
  if(!fp)
    {
      std::cerr << " File open fail. [" << InFileName << "]" << std::endl;
      return 2;
    }

  double slewfac=0.15;
  char str[MAXCHAR],str1[MAXCHAR];
  int CID[2]={CID_NC,CID_NC};
  char* CName[2]={"NC","NC"};
  int runnum,Layer[2],sseg[2][4];
  int nset=2;
  int ita_ini,ita_fin,writeparam=0;;
  std::string conffile,slewfile;
  
  while( fp.getline(str,MAXCHAR) )
    {
      if( str[0]=='#' ) continue;
      if( sscanf(str,"Run: %d", &runnum)==1 ) 
	std::cout<<"#Run= "<<runnum<<std::endl;
      else if( sscanf(str,"NumberOfSet: %d",&nset)==1 ) 
      	std::cout<<"NumberOfSet= "<<nset<<std::endl;
      else if( sscanf(str,"CorrectionFactor: %lf",&slewfac)==1 ) 
      	std::cout<<"CorrectionFactor= "<<slewfac<<std::endl;
      else if( sscanf(str,"Layer: %d %d", &Layer[0],&Layer[1])==2 ) 
      	std::cout<<"Layer[0]= "<<Layer[0]<<" Layer[1]= "<<Layer[1]<<std::endl;
      else if( sscanf(str,"segment[0]: %d %d %d %d", &sseg[0][0], &sseg[0][1], &sseg[0][2], &sseg[0][3])==4 )
      	std::cout<<"sseg[0]= "<<sseg[0][0]<<" "<<sseg[0][1]<<" "<<sseg[0][2]<<" "<<sseg[0][3]<<std::endl;
      else if( sscanf(str,"segment[1]: %d %d %d %d", &sseg[1][0], &sseg[1][1], &sseg[1][2], &sseg[1][3])==4 )
      	std::cout<<"sseg[1]= "<<sseg[1][0]<<" "<<sseg[1][1]<<" "<<sseg[1][2]<<" "<<sseg[1][3]<<std::endl;
      else if( sscanf(str,"itanum: %d %d", &ita_ini, &ita_fin)==2 )
	std::cout<<"itanum= "<<ita_ini<<" "<<ita_fin<<std::endl;
      else if( sscanf(str,"Conffile: %s", str1)==1 ){
	conffile=str1;
	std::cout<<"conffile= "<<conffile<<std::endl;
      }
      else if( sscanf(str,"Slewfile: %s", str1)==1 ){
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

  TH1F *h1;
  TH2F *h2;  
  double submean[2][4];
  double subsigma[2][4];
  double TOFsigma[2][6]={0};
  double TOFsigma_err[2][6]={0};
  double Sigma[2][4][3]={0};
  double Sigma_err[2][4][3]={0};

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */

  TFile *f;
  f = new TFile( Form( "./root/NCcalib/NCr000%d.root",runnum) );
  TTree *evtree = (TTree*)f->Get("EventTree");
  int nev = evtree->GetEntries();
  /*** declaration of classes ***/
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "EventHeader", &head );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );

  /*                */
  /* event analysis */
  /*                */
  TFile *fout;
  
  for(int itanum=ita_ini;itanum<=ita_fin;itanum++)
    {
      fout = new TFile( Form( "./root/NCcalib/NCr000%d_ita%d.root",
			      runnum,itanum),"recreate" );

      //======== define histograms--------------
  for(int iset=0;iset<nset;iset++)
    for(int iseg1=0;iseg1<4;iseg1++){
      int seg1=rseg[iset][iseg1];
      new TH1F(Form("TSub%s%d",CName[iset], seg1 ),Form("TSub%s%d",CName[iset],seg1 ),
	       400,-10,10);
      new TH2F( Form("%sU%dTDCEvNum",CName[iset],seg1), Form("%sU%dTDCEvNum Run%d",CName[iset],seg1,runnum),
		100,    0, nev,  400,   -0.5, 1999.5 );
      new TH2F( Form("%sD%dTDCEvNum",CName[iset],seg1), Form("%sD%dTDCEvNum Run%d",CName[iset],seg1,runnum),
		100,    0, nev,  400,   -0.5, 1999.5 );
      new TH2F( Form("%sU%dADCEvNum",CName[iset],seg1), Form("%sU%dADCEvNum Run%d",CName[iset],seg1,runnum),
		100,    0, nev,  400,   -0.5, 1999.5 );
      new TH2F( Form("%sD%dADCEvNum",CName[iset],seg1), Form("%sD%dADCEvNum Run%d",CName[iset],seg1,runnum),
		100,    0, nev,  400,   -0.5, 1999.5 );
      for(int iseg2=iseg1+1;iseg2<4;iseg2++)
	{
	  int seg2=rseg[iset][iseg2];
	  new TH2F( Form("%sU%d_E_TOF%d_%d",CName[iset],seg1,seg1,seg2),   
		    Form("Eu of %s%d TOF %d-%d corr. ",CName[iset],seg1,seg1,seg2),  
		    550,    -100, 1000,  400,   -5, 5 );
	  new TH2F( Form("%sD%d_E_TOF%d_%d",CName[iset],seg1,seg1,seg2),
		    Form("Ed of %s%d TOF %d-%d corr. ",CName[iset],seg1,seg1,seg2),
		    550,    -100, 1000,  400,   -5, 5 );
	  new TH2F( Form("%sU%d_E_TOF%d_%d",CName[iset],seg2,seg1,seg2),
		    Form("Eu of %s%d TOF %d-%d corr. ",CName[iset],seg2,seg1,seg2),
		    550,    -100, 1000,  400,   -5, 5 );
	  new TH2F( Form("%sD%d_E_TOF%d_%d",CName[iset],seg2,seg1,seg2),
		    Form("Ed of %s%d TOF %d-%d corr. ",CName[iset],seg2,seg1,seg2),
		    550,    -100, 1000,  400,   -5, 5 );
	  new TH1F( Form("%s_TOF%d_%d",CName[iset],seg1,seg2), 
		    Form("TOF %s%d-%d Run%d",CName[iset],seg1,seg2,runnum),  400,   -5, 5 );
	}
    }
  ConfMan *conf = new ConfMan(conffile);
  conf->Initialize();
  
  std::cout << " AllEvent : " << nev << std::endl;
  std::cout << "|0%                  |50%                |100% "     
	    << nev <<" run "<<runnum<<" ita "<<itanum<<  std::endl;

  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);
    if( iev==0 )
      std::cout << "|";
    if( iev%(nev/40)==0 )
      std::cout << "*"<<std::flush;
    if( iev==nev-1 )
      std::cout << "| fin!" << std::endl;

    blMan->Calc(conf);
    for(int cset=0;cset<nset;cset++)
      for( int i=0; i<4; i++ ){
	HodoscopeLikeHit *hit = blMan->Hodo(CID[cset],rseg[cset][i]);
	if(!hit) continue;
	int seg1 = hit->seg();
	if(seg1!=rseg[cset][i]) continue;
	h2=(TH2F*)gDirectory->Get(Form("%sU%dTDCEvNum",CName[cset],seg1));
	h2->Fill(iev,hit->tdcu());      
	h2=(TH2F*)gDirectory->Get(Form("%sD%dTDCEvNum",CName[cset],seg1));
	h2->Fill(iev,hit->tdcd());      
	h2=(TH2F*)gDirectory->Get(Form("%sU%dADCEvNum",CName[cset],seg1));
	h2->Fill(iev,hit->adcu());      
	h2=(TH2F*)gDirectory->Get(Form("%sD%dADCEvNum",CName[cset],seg1));
	h2->Fill(iev,hit->adcd());      
      }
    
    bool IsGood[2]={true,true};
    for(int cset=0;cset<nset;cset++)
      for( int i=0; i<4; i++ ){
	HodoscopeLikeHit *hit = blMan->Hodo(CID[cset],rseg[cset][i]);
	if(!hit){
	  IsGood[cset]=false;
	  continue;
	}
	int seg = hit->seg();
	if(seg!=rseg[cset][i]){
	  IsGood[cset]=false;
	  continue;
	}
	double ctsub = hit->ctsub();
	double eu = hit->eu(), ed = hit->ed();
	if( !(hit->CheckRange() && eu >20 && ed>20) ){
	  IsGood[cset]=false;
	  continue;
	}
	h1=(TH1F*)gDirectory->Get(Form("TSub%s%d",CName[cset],seg ));
	h1->Fill(ctsub);
	//	if((itanum!=ita_ini)&&fabs(ctsub-submean[cset][i])>1.*subsigma[cset][i] ) 
	//	  IsGood[cset]=false;;
      }
  
    //determine current set
    int cset;
    if(IsGood[0])         cset=0;
    else if(IsGood[1])    cset=1;
    else continue;
    //---------start filling histogram----------------------
    for( int i=0; i<4; i++ ){
      HodoscopeLikeHit *hit = blMan->Hodo(CID[cset],rseg[cset][i]);
      int seg1 = hit->seg();
      double eu1 = hit->eu(), ed1 = hit->ed();
      double ctm1 = hit->ctmean();
      for( int ii=i+1; ii<4; ii++ ){
	HodoscopeLikeHit *hit2 = blMan->Hodo(CID[cset],rseg[cset][ii]);
	int seg2 = hit2->seg();
	double eu2 = hit2->eu(), ed2 = hit2->ed();
	double ctm2 = hit2->ctmean();
	h2 = (TH2F*)gDirectory->Get( Form("%sU%d_E_TOF%d_%d",CName[cset],seg1,seg1,seg2) );  h2->Fill( eu1, ctm1-ctm2 );	
	h2 = (TH2F*)gDirectory->Get( Form("%sD%d_E_TOF%d_%d",CName[cset],seg1,seg1,seg2) );  h2->Fill( ed1, ctm1-ctm2 );
	h2 = (TH2F*)gDirectory->Get( Form("%sU%d_E_TOF%d_%d",CName[cset],seg2,seg1,seg2) );  h2->Fill( eu2, ctm1-ctm2 );
	h2 = (TH2F*)gDirectory->Get( Form("%sD%d_E_TOF%d_%d",CName[cset],seg2,seg1,seg2) );  h2->Fill( ed2, ctm1-ctm2 );	
	h1 = (TH1F*)gDirectory->Get( Form("%s_TOF%d_%d",CName[cset],seg1,seg2) ); h1->Fill(ctm1-ctm2 );		
      }
    }    
  }
  ///// histogram fill end
#if SLEWING
  std::cout<<"slewing correction start"<<std::endl;
  /// slewing correction start
#if SlewType
  TF1 *slewing=new TF1("slewing","[0]/sqrt(x)+[1]+[2]*x");
#else
  TF1 *slewing=new TF1("slewing","[0]/sqrt(x)+[1]");
#endif
  slewing->SetLineColor(2);
  for(int iset=0;iset<nset;iset++){
    int iseg12=0;
    for(int iseg1=0;iseg1<4;iseg1++){
      for(int iseg2=iseg1+1;iseg2<4;iseg2++){      
	iseg12++;
	int seg1=rseg[iset][iseg1];
	int seg2=rseg[iset][iseg2];
	for(int ud=0;ud<2;ud++){      	  
	  if(ud==0)
	    h2 = (TH2F*)gDirectory->Get( Form("%sU%d_E_TOF%d_%d",CName[iset],seg1,seg1,seg2) );
	  else 
	    h2 = (TH2F*)gDirectory->Get( Form("%sD%d_E_TOF%d_%d",CName[iset],seg1,seg1,seg2) );
	  if(h2->GetEntries()<100) continue;
	  TProfile* prof=h2->ProfileX();
	  slewing->SetParameters(0,0);	  
	  prof->Fit("slewing","W0q","",20,100);		  
	  prof->Fit("slewing","I0q","",20,100);
	  double par[3];
	  std::vector<double> apar;
	  int nth,type;
	  double chi2,dof;
	  prof->GetFunction("slewing")->GetParameters(par);
	  dof=slewing->GetNDF();
	  chi2=slewing->GetChisquare();
	  
	  conf->GetSlewingMapManager()->GetParam( CID[iset],seg1,ud,0,type,nth,apar ); // get a new parameter into map
	  int itmp;
	  for(itmp=0;itmp<nth;itmp++)
	    par[itmp]=apar[itmp]-slewfac*par[itmp];
#if SlewType 
	  if(nth==2)
	    par[2]=-slewfac*par[itmp];
	  conf->GetSlewingMapManager()->SetParam( CID[iset],seg1,ud,0,2,3,par ); // set a new parameter into map
#else	      
	  conf->GetSlewingMapManager()->SetParam( CID[iset],seg1,ud,0,1,2,par ); // set a new parameter into map
#endif
	  if(ud==0)
	    h2 = (TH2F*)gFile->Get( Form("%sU%d_E_TOF%d_%d",CName[iset],seg2,seg1,seg2) );
	  else 
	    h2 = (TH2F*)gFile->Get( Form("%sD%d_E_TOF%d_%d",CName[iset],seg2,seg1,seg2) );
	  if(h2->GetEntries()<100) continue;
	  prof=h2->ProfileX();
	  slewing->SetParameters(0,0);
	  prof->Fit("slewing","W0q","",20,100);		  
	  prof->Fit("slewing","I0q","",20,100);
	  prof->GetFunction("slewing")->GetParameters(par);
	  dof=slewing->GetNDF();
	  chi2=slewing->GetChisquare();

	  conf->GetSlewingMapManager()->GetParam( CID[iset],seg2,ud,0,type,nth,apar ); // get a new parameter into map
	  for(itmp=0;itmp<nth;itmp++)
	    par[itmp]=apar[itmp]+slewfac*par[itmp];
#if SlewType 
	  if(nth==2)
	    par[2]=+slewfac*par[itmp];
	  conf->GetSlewingMapManager()->SetParam( CID[iset],seg2,ud,0,2,3,par ); // set a new parameter into map
#else
	  conf->GetSlewingMapManager()->SetParam( CID[iset],seg2,ud,0,1,2,par ); // set a new parameter into map
#endif
	  h2->GetXaxis()->SetRangeUser(-50,300);
	  h2->GetYaxis()->SetRangeUser(-2,2);
	}
	h1=(TH1F*)gDirectory->Get( Form("%s_TOF%d_%d",CName[iset],seg1,seg2) );
	if(h1->GetEntries()<100) continue;
	h1->GetXaxis()->SetRangeUser(-2,2);
	h1->Fit("gaus","0q","",-10.2,10.2);
	double tmp=h1->GetFunction("gaus")->GetParameter(1);
	h1->Fit("gaus","0q","",tmp-0.5,tmp+0.5);
	TOFsigma[iset][iseg12-1]=h1->GetFunction("gaus")->GetParameter(2);
	TOFsigma_err[iset][iseg12-1]=h1->GetFunction("gaus")->GetParError(2);
	std::cout<<"tofsigma"<<seg1<<"-"<<seg2<<" = "<<TOFsigma[iset][iseg12-1]<<" / "<<h1->GetEntries()<<" events"<<std::endl;
      }
    }
    
  //solve rSigma
    for(int i=0;i<6;i++)  
      {
	TOFsigma_err[iset][i]=TOFsigma[iset][i]*TOFsigma_err[iset][i]*sqrt(2)*1e6;
	TOFsigma[iset][i]*=TOFsigma[iset][i]*1e6;
      }
    Sigma[iset][1][0]=sqrt( ( (TOFsigma[iset][0]-TOFsigma[iset][1])+TOFsigma[iset][3] )/2. ) ;
    Sigma[iset][2][0]=sqrt( ( -(TOFsigma[iset][0]-TOFsigma[iset][1])+TOFsigma[iset][3] )/2. ) ;
    Sigma[iset][0][0]=sqrt( TOFsigma[iset][0]-Sigma[iset][1][0]*Sigma[iset][1][0] ) ;
    Sigma[iset][3][0]=sqrt( TOFsigma[iset][5]-Sigma[iset][2][0]*Sigma[iset][2][0] ) ;
    
    Sigma_err[iset][1][0]=sqrt(0.5)/sqrt(Sigma[iset][1][0])*sqrt( TOFsigma_err[iset][0]+TOFsigma_err[iset][1]+TOFsigma_err[iset][3])/2.;
    Sigma_err[iset][2][0]=sqrt(0.5)/sqrt(Sigma[iset][2][0])*sqrt( TOFsigma_err[iset][0]+TOFsigma_err[iset][1]+TOFsigma_err[iset][3] )/2.;
    Sigma_err[iset][0][0]=sqrt(0.5)/sqrt(Sigma[iset][0][0])*sqrt( TOFsigma_err[iset][0]+Sigma[iset][1][0]*Sigma_err[iset][1][0]*sqrt(2) ) ;
    Sigma_err[iset][3][0]=sqrt(0.5)/sqrt(Sigma[iset][3][0])*sqrt( TOFsigma_err[iset][5]+Sigma[iset][2][0]*Sigma_err[iset][2][0]*sqrt(2) ) ;
    
    
    Sigma[iset][2][1]=sqrt( ( (TOFsigma[iset][3]-TOFsigma[iset][4])+TOFsigma[iset][5] )/2. ) ;
    Sigma[iset][3][1]=sqrt( ( -(TOFsigma[iset][3]-TOFsigma[iset][4])+TOFsigma[iset][5] )/2. ) ;
    Sigma[iset][1][1]=sqrt( TOFsigma[iset][3]-Sigma[iset][2][1]*Sigma[iset][2][1] );
    Sigma[iset][0][1]=sqrt( TOFsigma[iset][0]-Sigma[iset][1][1]*Sigma[iset][1][1] ) ;
    
    Sigma_err[iset][2][1]=sqrt(0.5)/sqrt(Sigma[iset][2][1])*sqrt( TOFsigma_err[iset][3]+TOFsigma_err[iset][4]+TOFsigma_err[iset][5])/2.;
    Sigma_err[iset][3][1]=sqrt(0.5)/sqrt(Sigma[iset][3][1])*sqrt( TOFsigma_err[iset][3]+TOFsigma_err[iset][4]+TOFsigma_err[iset][5] )/2.;
    Sigma_err[iset][1][1]=sqrt(0.5)/sqrt(Sigma[iset][1][1])*sqrt( TOFsigma_err[iset][3]+Sigma[iset][2][1]*Sigma_err[iset][2][1]*sqrt(2) ) ;
    Sigma_err[iset][0][1]=sqrt(0.5)/sqrt(Sigma[iset][0][1])*sqrt( TOFsigma_err[iset][0]+Sigma[iset][1][1]*Sigma_err[iset][1][1]*sqrt(2) ) ;
  }
  

  int isetseg=0;
  for(int iset=0;iset<nset;iset++)
    {
      std::cout<<"#"<<CName[iset]<<" Layer"<<Layer[iset]<<"  Sigma "<<std::endl;
      for(int iseg=0;iseg<4;iseg++)
	{
	  isetseg++;
	  std::cout<<"#seg "<<sseg[iset][iseg]<<" ";
	  for(int n=0;n<2;n++) 	std::cout<<Sigma[iset][iseg][n]<<" +- "<<Sigma_err[iset][iseg][n]<<" ";
	  std::cout<<std::endl;
	  double xbin[3]={1,2,3};
	  double xerr[3]={0,0,0};	
	  TGraphErrors *g1;
	  g1=new TGraphErrors(2,xbin,Sigma[iset][iseg],xerr,Sigma_err[iset][iseg]);
	  g1->SetName(Form("%slayer%dSeg%d",CName[iset],Layer[iset],sseg[iset][iseg]));
	  g1->SetTitle(Form("%slayer%dSeg%dRun%d",CName[iset],Layer[iset],sseg[iset][iseg],runnum));
	  g1->GetYaxis()->SetRangeUser(30,170);
	  g1->Fit("pol0","wq0");
	  g1->Write();
	}
    }
    

  for( int iset=0;iset<nset; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      h1 = (TH1F*)gDirectory->Get( Form("TSub%s%d",CName[iset],rseg[iset][iseg]) );
      h1->Fit("gaus","0q");
      subsigma[iset][iseg]=h1->GetFunction("gaus")->GetParameter(2);
      submean[iset][iseg]=h1->GetFunction("gaus")->GetParameter(1);
      h1->Fit("gaus","0q","",submean[iset][iseg]-5*subsigma[iset][iseg],submean[iset][iseg]+5*subsigma[iset][iseg]);
      std::cout<<"sub "<<iset<<" " << iseg<<" " <<h1->GetFunction("gaus")->GetParameter(1)<<"\t"<<h1->GetFunction("gaus")->GetParameter(2)<<std::endl;
      subsigma[iset][iseg]=h1->GetFunction("gaus")->GetParameter(2);
      submean[iset][iseg]=h1->GetFunction("gaus")->GetParameter(1);
    }
  }
  
  if(writeparam==1)
    {
      ofstream ofs(slewfile.c_str());
      conf->GetSlewingMapManager()->PrintMapBL(ofs); // print header
      ofs.close();
      ofs.open(Form("SlewingmapBL_NCtmp_run%d_ita%d",runnum,itanum));
      conf->GetSlewingMapManager()->PrintMapBL(ofs); // print header
      ofs.close();
    }
#endif  
  fout->Write();
  fout->Close();
  
    }//end ita loop
  return 0;
}

