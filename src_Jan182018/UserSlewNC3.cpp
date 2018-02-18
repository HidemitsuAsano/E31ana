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

  double slewfac=0.2;
  char str[MAXCHAR],str1[MAXCHAR];
  int CID=CID_NC;
  //  int CName[2]={"NC","NC"};
  int runnum = 83;
  int Layer[2];
  int sseg[2][4];
  int nset=2;
  int ita_ini=15;
  int ita_fin=16;
  std::string conffile;
  std::string slewfile;
  int writeparam=0;
  
  while( fp.getline(str,MAXCHAR) )
    {
      if( str[0]=='#' )
	{
	  //	  ofp<<str<<std::endl;
	} 
      if( sscanf(str,"Run: %d", &runnum)==1 ) 
	std::cout<<"#Run= "<<runnum<<std::endl;
      else if( sscanf(str,"NumberOfSet: %d",&nset)==1 ) 
      	std::cout<<"NumberOfSet= "<<nset<<std::endl;
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
  int rseg[2][4];
  for(int iset=0;iset<nset;iset++){
    for(int iseg=0;iseg<4;iseg++){
      rseg[iset][iseg]=(Layer[iset]-1)*16+sseg[iset][iseg];
      std::cout<<"rseg "<<iset <<" "<<iseg<<" "<<rseg[iset][iseg]<<std::endl;;
    }
  }

  
  /*** conf file for new parameters ***/

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
  //  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  //  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader", &head );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );


  //-------------calc sub mean and sigma for event cut--------
  TH1F *h1;
  TH2F *h2;
  
  double submean[2][4];
  double subsigma[2][4];
  /*
  for( int iset=0;iset<nset; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      h1 = (TH1F*)f->Get( Form("TsubNC%d",rseg[iset][iseg]) );
      h1->Fit("gaus","0q");
      subsigma[iset][iseg]=h1->GetFunction("gaus")->GetParameter(2);
      submean[iset][iseg]=h1->GetFunction("gaus")->GetParameter(1);
      h1->Fit("gaus","0q","",submean[iset][iseg]-2*submean[iset][iseg],submean[iset][iseg]+2*subsigma[iset][iseg]);
      std::cout<<"sub "<<iset<<" " << iseg<<" " <<h1->GetFunction("gaus")->GetParameter(1)<<"\t"<<h1->GetFunction("gaus")->GetParameter(2)<<std::endl;
      subsigma[iset][iseg]=h1->GetFunction("gaus")->GetParameter(2);
      submean[iset][iseg]=h1->GetFunction("gaus")->GetParameter(1);
    }
  }
  */
  //----------------------------------------------------------

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
      int seg1=sseg[iset][iseg1]+(Layer[iset]-1)*16;
      new TH1F(Form("TSubNC%d",seg1 ),Form("TSubNC%d",seg1 ),
	       400,-10,10);
      new TH2F( Form("NCU%dTDCEvNum",seg1), Form("NCU%dTDCEvNum",seg1),
		100,    0, nev,  400,   -0.5, 1999.5 );
      new TH2F( Form("NCD%dTDCEvNum",seg1), Form("NCD%dTDCEvNum",seg1),
		100,    0, nev,  400,   -0.5, 1999.5 );
      new TH2F( Form("NCU%dADCEvNum",seg1), Form("NCU%dADCEvNum",seg1),
		100,    0, nev,  400,   -0.5, 1999.5 );
      new TH2F( Form("NCD%dADCEvNum",seg1), Form("NCD%dADCEvNum",seg1),
		100,    0, nev,  400,   -0.5, 1999.5 );
      for(int iseg2=iseg1+1;iseg2<4;iseg2++)
	{
	  int seg2=sseg[iset][iseg2]+(Layer[iset]-1)*16;
	  //	      std::cout << "Define Histgram " <<seg1<<" "<<seg2<< std::endl;
	  new TH2F( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2),   Form("Eu of NC%d TOF NC%d-%d corr. ",seg1,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	  new TH2F( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2),   Form("Ed of NC%d TOF NC%d-%d corr. ",seg1,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	  new TH2F( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2),   Form("Eu of NC%d TOF NC%d-%d corr. ",seg2,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	  new TH2F( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2),   Form("Ed of NC%d TOF NC%d-%d corr. ",seg2,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	  new TH1F( Form("NC_TOF%d_%d",seg1,seg2),   Form("TOF NC%d-%d",seg1,seg2),  400,   -5, 5 );
	}
    }
  ConfMan *conf = new ConfMan(conffile);
  conf->Initialize();
  
  std::cout << " AllEvent : " << nev << std::endl;
  std::cout << "|0%                  |50%                |100% "     
	    << nev <<" run "<<runnum<<" ita "<<itanum<<  std::endl;

  for( int iev=0; iev<nev; iev++ ){
    evtree->GetEvent(iev);
    //    if( iev>2000 ) break;
    //    if( iev%10000==0 )
    //      std::cout << " Event : " << iev << std::endl;
    if( iev==0 )
      std::cout << "|";
    if( iev%(nev/40)==0 )
      std::cout << "*"<<std::flush;
    if( iev==nev-1 )
      std::cout << "| fin!" << std::endl;

    blMan->Calc(conf);
    ///-----------check trigger hit ----------
    /****
    int ntrig=0;
    for( int i=0; i<blMan->nPC(); i++ ){
      HodoscopeLikeHit *hit = blMan->PC(i);
      int seg = hit->seg();
      int au = hit->adcu();
      int tu = hit->tdcu();
      if( tu>0 )
	{
	  //      std::cout<<"PC seg "<<seg<<std::endl;
	  if(seg==1 && (au>80 && au<200) ) ntrig++;
	  if(seg==2 && (au>80 && au<220)) ntrig++;
	  if(seg==3 &&(au>150 && au<340)) ntrig++;
	  if(seg==4 &&(au>120 && au<250)) ntrig++;
	}
    }
    if(ntrig!=2) continue;    
    ****/

    ///-----------check trigger hit end -------------
    for(int cset=0;cset<nset;cset++)
      for( int i=0; i<4; i++ ){
	//!!!! it should be seg1<seg2 
	HodoscopeLikeHit *hit = blMan->Hodo(CID,rseg[cset][i]);
	int seg1 = hit->seg();
	//    std::cout<<"NC seg "<<seg<<std::endl;
	h2=(TH2F*)gDirectory->Get(Form("NCU%dTDCEvNum",seg1));
	h2->Fill(iev,hit->tdcu());      
	h2=(TH2F*)gDirectory->Get(Form("NCD%dTDCEvNum",seg1));
	h2->Fill(iev,hit->tdcd());      
	h2=(TH2F*)gDirectory->Get(Form("NCU%dADCEvNum",seg1));
	h2->Fill(iev,hit->adcu());      
	h2=(TH2F*)gDirectory->Get(Form("NCD%dADCEvNum",seg1));
	h2->Fill(iev,hit->adcd());      
      }
    
    ///-------require all for counters in a set-------------
    int ncnum[2]={0,0};
    bool NCSUB[2]={true,true};
    for( int i=0; i<blMan->nNC(); i++ ){
      HodoscopeLikeHit *hit = blMan->NC(i);
      int seg = hit->seg();
      //    std::cout<<"NC seg "<<seg<<std::endl;
      double ctsub = hit->ctsub();
      double eu = hit->eu(), ed = hit->ed();
      if( !(hit->CheckRange() && eu >30 && ed>30) ) continue;
      int isegtmp=0,isettmp=0;
      for(int iset=0;iset<nset;iset++)	
	for(int iseg=0;iseg<4;iseg++)	    
	  if(  rseg[iset][iseg]==seg  )	
	    {
	      ncnum[iset]++;isettmp=iset;isegtmp=iseg;
	      h1=(TH1F*)gDirectory->Get(Form("TSubNC%d",seg ));
	      h1->Fill(ctsub);
	    }      
      if((itanum!=ita_ini)&&fabs(ctsub-submean[isettmp][isegtmp])>1.*subsigma[isettmp][isegtmp] ) 
	//      if(fabs(ctsub-submean[isettmp][isegtmp])>0.3 ) 
	{
	  NCSUB[isettmp]=false;
	}    	    
    }
    //determine current set
    int cset;
    if((ncnum[0]==4 && ncnum[1]==0))
      cset=0;
    else if(ncnum[1]==4 && ncnum[0]==0)
      cset=1;
    else
      continue;
    //------------------------------------------------------
    
    //---------start filling histogram----------------------

   
    for( int i=0; i<4; i++ ){
      //!!!! it should be seg1<seg2 
      HodoscopeLikeHit *hit = blMan->Hodo(CID,rseg[cset][i]);
      int seg1 = hit->seg();
      //    std::cout<<"NC seg "<<seg<<std::endl;
      double eu1 = hit->eu(), ed1 = hit->ed();
      double ctm1 = hit->ctmean();
      for( int ii=i+1; ii<4; ii++ ){
	HodoscopeLikeHit *hit2 = blMan->Hodo(CID,rseg[cset][ii]);
	int seg2 = hit2->seg();
	double eu2 = hit2->eu(), ed2 = hit2->ed();
	double ctm2 = hit2->ctmean();
	//	std::cout<<"NC seg1,2 "<<seg1<<" "<<seg2<<std::endl;
	//	    std::cout <<"seg1 seg2  "<<seg1<<" "<<seg2  << std::endl;
	if(NCSUB[cset]){
	h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2) );  h2->Fill( eu1, ctm1-ctm2 );	
	h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2) );  h2->Fill( ed1, ctm1-ctm2 );
	h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2) );  h2->Fill( eu2, ctm1-ctm2 );
	h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2) );  h2->Fill( ed2, ctm1-ctm2 );	

	  h1 = (TH1F*)gDirectory->Get( Form("NC_TOF%d_%d",seg1,seg2) );
	  h1->Fill(ctm1-ctm2 );		
	}
      }
    }
    
  }
  ///// histogram fill end

#if SlewType
  TF1 *slewing=new TF1("slewing","[0]/sqrt(x)+[1]+[2]*x");
#else
  TF1 *slewing=new TF1("slewing","[0]/sqrt(x)+[1]");
#endif
  //  TF1 *pol=new TF1("pol","[0]+[1]*x");

  slewing->SetLineColor(2);

  /// slewing correction
  for(int iset=0;iset<nset;iset++){
    int iseg12=0;
    for(int iseg1=0;iseg1<4;iseg1++){
      for(int iseg2=iseg1+1;iseg2<4;iseg2++){      
	iseg12++;
	int seg1=sseg[iset][iseg1]+(Layer[iset]-1)*16;
	int seg2=sseg[iset][iseg2]+(Layer[iset]-1)*16;	
	for(int ud=0;ud<2;ud++){      	  
	  if(ud==0)
	    h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2) );
	  else 
	    h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2) );
	  //	  std::cout<<"seg1 "<<seg1<<" seg2 "<<seg2<<"h2 en "<<h2->GetEntries()<<std::endl;
	  TProfile* prof=h2->ProfileX();
	  slewing->SetParameters(0,0);
	  
	  prof->Fit("slewing","W0q","",50,150);		  
	  prof->Fit("slewing","I0q","",50,150);
	  double par[3];
	  std::vector<double> apar;
	  //	  double apar[10];
	  int nth,type;
	  double chi2,dof;
	  prof->GetFunction("slewing")->GetParameters(par);
	  dof=slewing->GetNDF();
	  chi2=slewing->GetChisquare();
	  
	  conf->GetSlewingMapManager()->GetParam( CID,seg1,ud,0,type,nth,apar ); // get a new parameter into map
	  int itmp;
	  for(itmp=0;itmp<nth;itmp++)
	    {
	      par[itmp]=apar[itmp]-slewfac*par[itmp];
	    }
#if SlewType 
	  if(nth==2)
	    par[2]=-slewfac*par[itmp];
	  //	  if(iseg1==0 && iseg2==1)
	  conf->GetSlewingMapManager()->SetParam( CID,seg1,ud,0,2,3,par ); // set a new parameter into map
#else	      
	  //	  if(iseg1==0 && iseg2==1)
	  conf->GetSlewingMapManager()->SetParam( CID,seg1,ud,0,1,2,par ); // set a new parameter into map
#endif
	  
	  //	  std::cout<<"tof"<<seg1<<"-"<<seg2<<" NC seg :ud "<<seg1<<" "<<ud <<" par "<<par[0]<<" "<<par[1]<<std::endl;
	  
	  if(ud==0)
	    h2 = (TH2F*)gFile->Get( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2) );
	  else 
	    h2 = (TH2F*)gFile->Get( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2) );
	  prof=h2->ProfileX();
	  slewing->SetParameters(0,0);
	  prof->Fit("slewing","W0q","",50,150);		  
	  prof->Fit("slewing","I0q","",50,150);
	  prof->GetFunction("slewing")->GetParameters(par);
	  dof=slewing->GetNDF();
	  chi2=slewing->GetChisquare();

	  conf->GetSlewingMapManager()->GetParam( CID,seg2,ud,0,type,nth,apar ); // get a new parameter into map
	  for(itmp=0;itmp<nth;itmp++)
	    {
	      par[itmp]=apar[itmp]+slewfac*par[itmp];
	    }
#if SlewType 
	  if(nth==2)
	    par[2]=+slewfac*par[itmp];
	  //	  if(iseg1==0)
	  conf->GetSlewingMapManager()->SetParam( CID,seg2,ud,0,2,3,par ); // set a new parameter into map
#else
	  //	  if(iseg1==0)
	  conf->GetSlewingMapManager()->SetParam( CID,seg2,ud,0,1,2,par ); // set a new parameter into map
#endif
	  
	  //	  std::cout<<"tof"<<seg1<<"-"<<seg2<<" NC seg :ud "<<seg2<<" "<<ud <<" par "<<par[0]<<" "<<par[1]<<std::endl;
	  h2->GetXaxis()->SetRangeUser(-50,300);
	  h2->GetYaxis()->SetRangeUser(-2,2);
	}
	h1=(TH1F*)gDirectory->Get( Form("NC_TOF%d_%d",seg1,seg2) );
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
      //	std::cout<<"TOFsigma "<<TOFsigma[iset][i]<<"  +- "<<TOFsigma_err[iset][i]<<std::endl;;
      TOFsigma_err[iset][i]=TOFsigma[iset][i]*TOFsigma_err[iset][i]*sqrt(2)*1e6;
	TOFsigma[iset][i]*=TOFsigma[iset][i]*1e6;
	//	std::cout<<"TOFsigma "<<TOFsigma[iset][i]<<"  +- "<<TOFsigma_err[iset][i]<<std::endl;
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
    std::cout<<"#Layer"<<Layer[iset]<<"  Sigma "<<std::endl;;
    for(int iseg=0;iseg<4;iseg++)
      {
	isetseg++;
	std::cout<<"#seg "<<sseg[iset][iseg]<<" ";
	for(int n=0;n<3;n++) 	std::cout<<Sigma[iset][iseg][n]<<" +- "<<Sigma_err[iset][iseg][n]<<" ";
	std::cout<<std::endl;
	double xbin[3]={1,2,3};
	double xerr[3]={0,0,0};	
	TGraphErrors *g1;
	g1=new TGraphErrors(2,xbin,Sigma[iset][iseg],xerr,Sigma_err[iset][iseg]);
	g1->SetName(Form("layer%dSeg%d",Layer[iset],sseg[iset][iseg]));
	g1->SetTitle(Form("layer%dSeg%d",Layer[iset],sseg[iset][iseg]));
	g1->GetYaxis()->SetRangeUser(30,170);
	g1->Fit("pol0","wq0");
	g1->Write();
      }
    }

  for( int iset=0;iset<nset; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      h1 = (TH1F*)gDirectory->Get( Form("TSubNC%d",rseg[iset][iseg]) );
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
  
  fout->Write();
  fout->Close();
  
    }//end ita loop
  return 0;
}

