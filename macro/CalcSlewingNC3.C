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

#define SlewType 1

void CalcSlewingNC3()
{

  //#####################


  //#####################

  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");






  //########parameter file

#define MAXCHAR 256   
  ifstream fp;
  string InFileName="SlewNC.param";
  char str[MAXCHAR],str1[MAXCHAR];

  int runnum = 83;
  int Layer[2]={6,7};
  int sseg[2][4]={9,10,11,12,9,10,11,12};
  int ita_ini=15;
  int ita_fin=16;
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

      if( str[0]=='#' )
	{
	  //	  ofp<<str<<std::endl;
	} 

      if( sscanf(str,"Run: %d", &runnum)==1 ) 
	std::cout<<"#Run= "<<runnum<<std::endl;
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
  for(int iset=0;iset<2;iset++){
    for(int iseg=0;iseg<4;iseg++){
      rseg[iset][iseg]=(Layer[iset]-1)*16+sseg[iset][iseg];
      cout<<"rseg "<<iset <<" "<<iseg<<rseg[iset][iseg]<<endl;;
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

  //########  Fill Hist ###############//
  
  /*** declaration of classes ***/
  CDSHitMan *cdsMan = 0;
  BeamLineHitMan *blMan = 0;
  EventHeader *head = 0;
  evtree->SetBranchAddress( "CDSHitMan", &cdsMan );
  evtree->SetBranchAddress( "EventHeader", &head );
  evtree->SetBranchAddress( "BeamLineHitMan", &blMan );

  double submean[2][4];
  double subsigma[2][4];
  TH1F *h1;
  for( int iset=0;iset<2; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      h1 = (TH1F*)f->Get( Form("TsubNC%d",rseg[iset][iseg]) );
      h1->Fit("gaus","0q");
      subsigma[iset][iseg]=h1->GetFunction("gaus")->GetParameter(2);
      submean[iset][iseg]=h1->GetFunction("gaus")->GetParameter(1);
      h1->Fit("gaus","0q","",submean[iset][iseg]-2*submean[iset][iseg],submean[iset][iseg]+2*subsigma[iset][iseg]);
      subsigma[iset][iseg]=h1->GetFunction("gaus")->GetParameter(2);
      submean[iset][iseg]=h1->GetFunction("gaus")->GetParameter(1);

    }
  }




  /*                */
  /* event analysis */
  /*                */
  for(int itanum=ita_ini;itanum<ita_fin;itanum++)
    {


  for(int iset=0;iset<2;iset++)
    {
      for(int iseg1=0;iseg1<4;iseg1++)
	{
	  for(int iseg2=iseg1+1;iseg2<4;iseg2++)
	    {
	      int seg1=sseg[iset][iseg1]+(Layer[iset]-1)*16;
	      int seg2=sseg[iset][iseg2]+(Layer[iset]-1)*16;
	      if(seg1>seg2)
		{
		  int segtmp=seg1;
		  seg1=seg2;
		  seg2=segtmp;
		}
	      //	      std::cout << "Define Histgram " <<seg1<<" "<<seg2<< std::endl;
	      new TH2F( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2),   Form("Eu of NC%d TOF NC%d-%d corr. ",seg1,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	      new TH2F( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2),   Form("Ed of NC%d TOF NC%d-%d corr. ",seg1,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	      new TH2F( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2),   Form("Eu of NC%d TOF NC%d-%d corr. ",seg2,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	      new TH2F( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2),   Form("Ed of NC%d TOF NC%d-%d corr. ",seg2,seg1,seg2),     550,    -100, 1000,  400,   -5, 5 );
	    }
	}
    }

  ConfMan *conf = new ConfMan(conffile);
			      //"analyzerNCcalib2tmp2.conf");
  conf->Initialize();

  int nev = evtree->GetEntries();
  int ntrack=0;
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
      std::cout << "*"<<flush;
    if( iev==nev-1 )
      std::cout << "| fin!" << std::endl;



    int ntrig=0;
    for( int i=0; i<blMan->nPC(); i++ ){
      HodoscopeLikeHit *hit = blMan->PC(i);
      int seg = hit->seg();

      int au = hit->adcu(), ad = hit->adcd();
      double eu = hit->eu(), ed = hit->ed();
      int tu = hit->tdcu(), td = hit->tdcd();
      double ctm = hit->ctmean();		
      if( tu>0 )
	{
	  //      std::cout<<"PC seg "<<seg<<std::endl;
		if(seg==1 && (au>80 && au<200) ) ntrig++;
		if(seg==2 && (au>80 && au<220)) ntrig++;
		if(seg==3 &&(au>150 && au<340)) ntrig++;
		if(seg==4 &&(au>120 && au<250)) ntrig++;
	}
    }

    int ncnum[2]={0,0};
    for( int i=0; i<blMan->nNC(); i++ ){
      if(ntrig!=2) continue;
      HodoscopeLikeHit *hit = blMan->NC(i);
      int seg = hit->seg();
      //    std::cout<<"NC seg "<<seg<<std::endl;
      int au = hit->adcu(), ad = hit->adcd();
      double eu = hit->eu(), ed = hit->ed();
      int tu = hit->tdcu(), td = hit->tdcd();
      double ctm = hit->ctmean();
      double tsub = hit->tsub();
      if( !(eu >30 && ed>30) ) continue;
      int isegtmp=0,isettmp=0;
      for(int n=0;n<2;n++)
	{
	  for(int ii=0;ii<4;ii++)
	    {
	      if(  rseg[n][ii]==seg  )	{ncnum[n]++;isettmp=n;isegtmp=ii;}
	    }
	}

      if(fabs(tsub-submean[isettmp][isegtmp])>subsigma[isettmp][isegtmp] ) 
	{
	  if(isettmp==0) ncnum[0]--;
	  else if(isettmp==1) ncnum[1]--;
	}    	    
    }

    if( !(ncnum[0]==4 && ncnum[1]==0) && !(ncnum[1]==4 && ncnum[0]==0) )
      continue;

    blMan->Calc(conf);
    for( int i=0; i<blMan->nNC(); i++ ){
      HodoscopeLikeHit *hit = blMan->NC(i);

      int seg = hit->seg();
      //    std::cout<<"NC seg "<<seg<<std::endl;
      int au = hit->adcu(), ad = hit->adcd();
      double eu = hit->eu(), ed = hit->ed();
      int tu = hit->tdcu(), td = hit->tdcd();
      double ctm = hit->ctmean();
      if( hit->CheckRange() ){
	for( int ii=i+1; ii<blMan->nNC(); ii++ ){
	  int seg1 = seg;
	  int au1 = au, ad1 = ad;
	  double eu1 = eu, ed1 = ed;
	  int tu1 = tu, td1 = td;
	  double ctm1 = ctm;
	  HodoscopeLikeHit *hit2 = blMan->NC(ii);
	  int seg2 = hit2->seg();
	  int au2 = hit2->adcu(), ad2 = hit2->adcd();
	  double eu2 = hit2->eu(), ed2 = hit2->ed();
	  int tu2 = hit2->tdcu(), td2 = hit2->tdcd();
	  double ctm2 = hit2->ctmean();

	  if( hit2->CheckRange() ){
	    //	    std::cout<<"NC seg1,2 "<<seg1<<" "<<seg2<<std::endl;
	    if(seg2<seg1)
	      {
		int segt = seg2;
		int aut = au2, adt = ad2;
		double eut = eu2, edt = ed2;
		int tut = tu2, tdt = td2;
		double ctmt = ctm2;
		seg2 = seg;	 au2 = au; ad2 =ad;
		eu2 = eu; ed2 = ed;
		tu2 = tu, td2 = td;
		ctm2 = ctm;	
		seg1 = segt; au1 = aut; ad1 =adt;
		eu1 = eut; ed1 = edt;
		tu1 = tut, td1 = tdt;
		ctm1 = ctmt;	        
	      }
	    bool pairflag=false;

	    if( ( rseg[0][0]==seg1 || rseg[0][1]==seg1 ||
		 rseg[0][2]==seg1 || rseg[0][3]==seg1 ) &&
		( rseg[0][0]==seg2 || rseg[0][1]==seg2 ||
		  rseg[0][2]==seg2 || rseg[0][3]==seg2 ) )
	      pairflag=true;
	    if( ( rseg[1][0]==seg1 || rseg[1][1]==seg1 ||
		 rseg[1][2]==seg1 || rseg[1][3]==seg1 ) &&
		( rseg[1][0]==seg2 || rseg[1][1]==seg2 ||
		  rseg[1][2]==seg2 || rseg[1][3]==seg2 ) )
	      pairflag=true;

	    if(!pairflag) continue;		 


	    //	    std::cout <<"seg1 seg2  "<<seg1<<" "<<seg2  << std::endl;
	    h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2) );  h2->Fill( eu1, ctm1-ctm2 );

	    h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2) );  h2->Fill( ed1, ctm1-ctm2 );
	    h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2) );  h2->Fill( eu2, ctm1-ctm2 );
	    h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2) );  h2->Fill( ed2, ctm1-ctm2 );

	  }
	}

      }
      
    }

  }

  
  TF1 *slewing=new TF1("slewing","[0]/sqrt(x)+[1]");
  TF1 *slewing2=new TF1("slewing2","[0]/sqrt(x)+[1]+[2]*x");
  TF1 *pol=new TF1("pol","[0]+[1]*x");

  slewing->SetLineColor(2);

  TH1F *h1;
  TH2F *h2;
  TCanvas *c1=new TCanvas();
  TCanvas *c2=new TCanvas();
  TCanvas *c3=new TCanvas();
  c1->Divide(3,2);
  c2->Divide(3,2);
  c3->Divide(3,2);
  TCanvas *c4=new TCanvas();
  TCanvas *c5=new TCanvas();
  TCanvas *c6=new TCanvas();
  c4->Divide(3,2);
  c5->Divide(3,2);
  c6->Divide(3,2);
  TCanvas *c7=new TCanvas();
  TCanvas *c8=new TCanvas();
  TCanvas *c9=new TCanvas();
  c7->Divide(3,2);
  c8->Divide(3,2);
  c9->Divide(3,2);
  TCanvas *c10=new TCanvas();
  TCanvas *c11=new TCanvas();
  TCanvas *c12=new TCanvas();
  c10->Divide(3,2);
  c11->Divide(3,2);
  c12->Divide(3,2);

  for(int iset=0;iset<2;iset++){
    int iseg12=0;
    for(int iseg1=0;iseg1<4;iseg1++){
      for(int iseg2=iseg1+1;iseg2<4;iseg2++){      
	iseg12++;
	for(int ud=0;ud<2;ud++){      
	  int seg1=sseg[iset][iseg1]+(Layer[iset]-1)*16;
	  int seg2=sseg[iset][iseg2]+(Layer[iset]-1)*16;	
	  if(seg1>seg2)
	    {
	      int segtmp=seg1;
	      seg1=seg2;
	      seg2=segtmp;
	    }
	  if(iset==0)
	    {
	      if(iseg1==0)
		{
		  if(iseg2==1) c1->cd(1+ud);
		  else if(iseg2==2) c2->cd(1+ud);
		  else if(iseg2==3) c3->cd(1+ud);
		}
	      else if(iseg1==1)
		{
		  if(iseg2==2) c4->cd(1+ud);
		  else if(iseg2==3) c5->cd(1+ud);
		}
	      if(iseg1==2)
		{
		  if(iseg2==3) c6->cd(1+ud);
		}
	    }
	  else 
	    {
	      if(iseg1==0)
		{
		  if(iseg2==1) c7->cd(1+ud);
		  else if(iseg2==2) c8->cd(1+ud);
		  else if(iseg2==3) c9->cd(1+ud);
		}
	      else if(iseg1==1)
		{
		  if(iseg2==2) c10->cd(1+ud);
		  else if(iseg2==3) c11->cd(1+ud);
		}
	      if(iseg1==2)
		{
		  if(iseg2==3) c12->cd(1+ud);
		}
	    }

	  if(ud==0)
	    h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg1,seg1,seg2) );
	  else 
	    h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg1,seg1,seg2) );
	  cout<<"seg1 "<<seg1<<" seg2 "<<seg2<<"h2 en "<<h2->GetEntries()<<endl;
	  TProfile *prof=h2->ProfileX();
	  slewing->SetParameters(0,0);
	  slewing2->SetParameters(0,0);
#if SlewType
	  prof->Fit("slewing2","Wq","",50,200);		  
	  prof->Fit("slewing2","Iq","",50,200);
#else
	  prof->Fit("slewing","Wq","",50,200);		  
	  prof->Fit("slewing","Iq","",50,200);
#endif
	  double par[3];
	  std::vector <double> apar;
	  int nth;
	  double chi2,dof;
#if SlewType
	  prof->GetFunction("slewing2")->GetParameters(par);
	  dof=slewing2->GetNDF();
	  chi2=slewing2->GetChisquare();
#else
	  prof->GetFunction("slewing")->GetParameters(par);
	  dof=slewing->GetNDF();
	  chi2=slewing->GetChisquare();
#endif

	  if(iseg1==0)
	    {
	      conf->GetSlewingMapManager()->GetParam( CID_NC,seg1,ud,0,1,nth,apar ); // get a new parameter into map
	      for(int itmp=0;itmp<nth;itmp++)
		{
		  par[itmp]=apar[itmp]-0.8*par[itmp];
		}
#if SlewType 
	      if(nth==2)
 		  par[2]=-0.8*par[itmp];
	      if(iseg1==0 && iseg2==1)
		conf->GetSlewingMapManager()->SetParam( CID_NC,seg1,ud,0,2,3,par ); // set a new parameter into map
#else
	      
	      if(iseg1==0 && iseg2==1)
		conf->GetSlewingMapManager()->SetParam( CID_NC,seg1,ud,0,1,2,par ); // set a new parameter into map
#endif
	    }
	  std::cout<<"tof"<<seg1<<"-"<<seg2<<" NC seg :ud "<<seg1<<" "<<ud <<" par "<<par[0]<<" "<<par[1]<<std::endl;
	  h2->GetXaxis()->SetRangeUser(-50,300);
	  h2->GetYaxis()->SetRangeUser(-2,2);
	  h2->Draw("colz");
#if SlewType 
	  prof->GetFunction("slewing2")->Draw("same");
#else
	  prof->GetFunction("slewing")->Draw("same");
#endif

	  if(ud==0)
	    h2 = (TH2F*)gDirectory->Get( Form("NCU%d_E_TOF%d_%d",seg2,seg1,seg2) );
	  else 
	    h2 = (TH2F*)gDirectory->Get( Form("NCD%d_E_TOF%d_%d",seg2,seg1,seg2) );
	  if(iset==0)
	    {
	      if(iseg1==0)
		{
		  if(iseg2==1) c1->cd(4+ud);
		  else if(iseg2==2) c2->cd(4+ud);
		  else if(iseg2==3) c3->cd(4+ud);
		}
	      else if(iseg1==1)
		{
		  if(iseg2==2) c4->cd(4+ud);
		  else if(iseg2==3) c5->cd(4+ud);
		}
	      if(iseg1==2)
		{
		  if(iseg2==3) c6->cd(4+ud);
		}
	    }
	  else 
	    {
	      if(iseg1==0)
		{
		  if(iseg2==1) c7->cd(4+ud);
		  else if(iseg2==2) c8->cd(4+ud);
		  else if(iseg2==3) c9->cd(4+ud);
		}
	      else if(iseg1==1)
		{
		  if(iseg2==2) c10->cd(4+ud);
		  else if(iseg2==3) c11->cd(4+ud);
		}
	      if(iseg1==2)
		{
		  if(iseg2==3) c12->cd(4+ud);
		}
	    }
	  TProfile *prof=h2->ProfileX();
	  slewing->SetParameters(0,0);
	  slewing2->SetParameters(0,0);
#if SlewType
	  prof->Fit("slewing2","Wq","",50,200);		  
	  prof->Fit("slewing2","Iq","",50,200);
#else
	  prof->Fit("slewing","Wq","",50,200);		  
	  prof->Fit("slewing","Iq","",50,200);
#endif
	  double par[3];
	  std::vector <double> apar;
	  int nth;
	  double chi2,dof;
#if SlewType
	  prof->GetFunction("slewing2")->GetParameters(par);
	  dof=slewing2->GetNDF();
	  chi2=slewing2->GetChisquare();
#else
	  prof->GetFunction("slewing")->GetParameters(par);
	  dof=slewing->GetNDF();
	  chi2=slewing->GetChisquare();
#endif
	  conf->GetSlewingMapManager()->GetParam( CID_NC,seg2,ud,0,1,nth,apar ); // get a new parameter into map
	  for(int itmp=0;itmp<nth;itmp++)
	    {
	      par[itmp]=apar[itmp]+0.8*par[itmp];
	    }
#if SlewType 
	      if(nth==2)
 		  par[2]=+0.8*par[itmp];
	      if(iseg1==0)
		conf->GetSlewingMapManager()->SetParam( CID_NC,seg2,ud,0,2,3,par ); // set a new parameter into map
#else
	      if(iseg1==0)
		conf->GetSlewingMapManager()->SetParam( CID_NC,seg2,ud,0,1,2,par ); // set a new parameter into map
#endif

	  std::cout<<"tof"<<seg1<<"-"<<seg2<<" NC seg :ud "<<seg2<<" "<<ud <<" par "<<par[0]<<" "<<par[1]<<std::endl;
	  h2->GetXaxis()->SetRangeUser(-50,300);
	  h2->GetYaxis()->SetRangeUser(-2,2);
 	  h2->Draw("colz");
#if SlewType 
	  prof->GetFunction("slewing2")->Draw("same");
#else
	  prof->GetFunction("slewing")->Draw("same");
#endif


	  if(iset==0)
	    {
	      if(iseg1==0)
		{
		  if(iseg2==1) c1->cd(3);
		  else if(iseg2==2) c2->cd(3);
		  else if(iseg2==3) c3->cd(3);
		}
	      else if(iseg1==1)
		{
		  if(iseg2==2) c4->cd(3);
		  else if(iseg2==3) c5->cd(3);
		}
	      if(iseg1==2)
		{
		  if(iseg2==3) c6->cd(3);
		}
	    }
	  else 
	    {
	      if(iseg1==0)
		{
		  if(iseg2==1) c7->cd(3);
		  else if(iseg2==2) c8->cd(3);
		  else if(iseg2==3) c9->cd(3);
		}
	      else if(iseg1==1)
		{
		  if(iseg2==2) c10->cd(3);
		  else if(iseg2==3) c11->cd(3);
		}
	      if(iseg1==2)
		{
		  if(iseg2==3) c12->cd(3);
		}
	    }

	  h2->ProjectionY()->GetXaxis()->SetRangeUser(-2,2);
	  h2->ProjectionY()->Fit("gaus","q","",-10.2,10.2);
	  double tmp=gaus->GetParameter(1);
	  h2->ProjectionY()->Fit("gaus","q","",tmp-0.5,tmp+0.5);
	  TOFsigma[iset][iseg12-1]=gaus->GetParameter(2);
	  TOFsigma_err[iset][iseg12-1]=gaus->GetParError(2);

	}
      }
    }
  //solve rSigma
    for(int i=0;i<6;i++)  
      {
	cout<<"TOFsigma "<<TOFsigma[iset][i]<<"  +- "<<TOFsigma_err[iset][i]<<endl;;
	TOFsigma_err[iset][i]=TOFsigma[iset][i]*TOFsigma_err[iset][i]*sqrt(2)*1e6;
	TOFsigma[iset][i]*=TOFsigma[iset][i]*1e6;
	cout<<"TOFsigma "<<TOFsigma[iset][i]<<"  +- "<<TOFsigma_err[iset][i]<<endl;;
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


//     Sigma[iset][2][2]=sqrt( ( (TOFsigma[iset][1]-TOFsigma[iset][2])+TOFsigma[iset][5] )/2. ) ;
//     Sigma[iset][3][2]=sqrt( ( -(TOFsigma[iset][1]-TOFsigma[iset][2])+TOFsigma[iset][5] )/2. ) ;
//     Sigma[iset][0][2]=sqrt( TOFsigma[iset][1]-Sigma[iset][2][2]*Sigma[iset][2][2] ) ;
//     Sigma[iset][1][2]=sqrt( TOFsigma[iset][0]-Sigma[iset][0][2]*Sigma[iset][0][2] ) ;

//     Sigma_err[iset][2][2]=sqrt(0.5)/sqrt(Sigma[iset][2][2])*sqrt( TOFsigma_err[iset][1]+TOFsigma_err[iset][2]+TOFsigma_err[iset][5])/2.;
//     Sigma_err[iset][3][2]=sqrt(0.5)/sqrt(Sigma[iset][3][2])*sqrt( TOFsigma_err[iset][1]+TOFsigma_err[iset][2]+TOFsigma_err[iset][5] )/2.;
//     Sigma_err[iset][0][2]=sqrt(0.5)/sqrt(Sigma[iset][0][2])*sqrt( TOFsigma_err[iset][1]+Sigma[iset][2][2]*Sigma_err[iset][2][2]*sqrt(2) ) ;
//     Sigma_err[iset][1][2]=sqrt(0.5)/sqrt(Sigma[iset][1][2])*sqrt( TOFsigma_err[iset][0]+Sigma[iset][0][2]*Sigma_err[iset][0][2]*sqrt(2) ) ;

  }


  TCanvas *c13=new TCanvas();
  c13->Divide(4,2);

  int isetseg=0;
  for(int iset=0;iset<2;iset++)
    {
    cout<<"#Layer"<<Layer[iset]<<"  Sigma "<<endl;;
    for(int iseg=0;iseg<4;iseg++)
      {
	isetseg++;
	cout<<"#seg "<<sseg[iset][iseg]<<" ";
	for(int n=0;n<3;n++) 	cout<<Sigma[iset][iseg][n]<<" +- "<<Sigma_err[iset][iseg][n]<<" ";
	cout<<endl;
	double xbin[3]={1,2,3};
	double xerr[3]={0,0,0};

	TGraphErrors *g1;
	g1=new TGraphErrors(2,xbin,Sigma[iset][iseg],xerr,Sigma_err[iset][iseg]);
	c13->cd(isetseg);
	g1->GetYaxis()->SetRangeUser(30,170);
	g1->Draw("AP");
	g1->Fit("pol0","wq");


	double fitsigma=g1->GetFunction("pol0")->GetParameter(0);
	double err_sigma=g1->GetFunction("pol0")->GetParError(0);
	TLatex latex;
	latex.SetTextColor(2);
	latex.SetTextSize(0.07);
	latex.DrawLatex( 1, 60, Form("Layer%d Seg%d ", Layer[iset],sseg[iset][iseg]) );
	latex.DrawLatex( 1, 50, Form("#sigma = %2.2lf +- %1.2lf", fitsigma,err_sigma) ); 
	  //;	latex.DrawLatex( ihv[ilay][iseg][ud][2]-60, 80, Form("HV@65[pC] =%4.2lf ", rhv[ilay][iseg][ud]) );

      }
    }


  c1->SaveAs(Form("NC_slew_run%d_ita%d.ps(",runnum,itanum));
  c2->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c3->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c4->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c5->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c6->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c7->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c8->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c9->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c10->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c11->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c12->SaveAs(Form("NC_slew_run%d_ita%d.ps",runnum,itanum));
  c13->SaveAs(Form("NC_slew_run%d_ita%d.ps)",runnum,itanum));

  if(writeparam==1)
    {
      ofstream ofs(slewfile.c_str());
      conf->GetSlewingMapManager()->PrintMapBL(ofs); // print header
      ofs.close();
    }
  
    }//end ita loop
  return;
}
