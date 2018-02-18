//asano memo
//macro to calculate x-t curve of BLDCs ? 
//
//Input
//Output

#include <iostream>
#include <fstream>
#include <iomanip.h>
#include "TMinuit.h"

#define DISPLAY 0
#define DISPLAY2 0

class GlobalVariables;
class ConfMan;
class Display;
class CDSHitMan;
class CDCHit;
class CDSTrackingMan;
class CDSTrack;
class CircleFit;
class HelixFit;

//######## Calc XT #### ///

void CalcXTofBLC()
{

  //#####################

  int itanum = 3;  
  int runnum = 3036;

  //#####################


  /*** load library ***/
  gSystem->Load("libPhysics.so");
  gSystem->Load("lib/libAll.so");


  /*** conf file for new parameters ***/
  //  ConfMan *conf = new ConfMan("conf/Oct2010/analyzer.conf");
  ConfMan *conf = new ConfMan("conf/Oct2010/analyzerBLCtmp.conf");
  conf->Initialize();

  /* --------------- */
  /* ---- start ---- */
  /* --------------- */



  TFile *f;

  
  f = new TFile( Form( "./root/xtout_blc%d_%d.root", itanum,runnum ) );

  //########  Calc XT ###############//
  

  
  TF1 *resl_ita=new TF1("resl_ita","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5");
  //  TF1 *resl_ita=new TF1("resl_ita","[1]*(x-[0])+[2]*(x-[0])*(x-[0])+[3]*(x-[0])*(x-[0])*(x-[0])+[4]*(x-[0])**4+[5]*(x-[0])**5");
  //   TF1 *resl_ita=new TF1("resl_ita","pol4",-100,800.);
  resl_ita->SetLineColor(2);
  //  resl_ita->FixParameter(0,0);
  //  resl_ita->FixParameter(1,0);  
  //  resl_ita->FixParameter(2,0);
  // resl_ita->FixParameter(3,0);
  //  resl_ita->FixParameter(4,0);
  //  resl_ita->FixParameter(5,0);
  
  double range=100;

#define MAXCHAR 256   
  ifstream fp;
  ofstream ofp;
  string InFileName=conf->GetXTMapFileNameBL();
  string OutFileName=Form("./XTforBLC/XTCurveBL_BLCtmp%d.param",itanum+1);
  char str[MAXCHAR];
  int nd,cid,layer,wire,npar;
  double p0,p1,p2,p3,p4,p5;
  double  fitp0[2][8][32];

  std::cout << InFileName <<endl;
  //  cout<<setprecision(6);
  if( fp.open(InFileName.c_str() )==0 )
    {
      std::cerr << " File open fail. [" << InFileName << "]" << std::endl;
      return;
    }

  if( ofp.open(OutFileName.c_str())==0 )
    {
      std::cerr << " File open fail. [" << OutFileName << "]" << std::endl;
      return;
    }

  while( fp.getline(str,MAXCHAR) )
    {
      std::cout << "testa1" <<endl;      
      if( str[0]=='#' )
	{
	  ofp<<str<<std::endl;
	} 
      else if( (nd=sscanf(str,"%d %d %d %d %lf %lf %lf %lf %lf %lf", &cid, &layer, &wire, &npar, &p0, &p1,  &p2, &p3, &p4, &p5 ) ) >=6 ) 
	{
	  
#if 1
	  std::cout << cid << "  " << layer << "  " << wire << "  "
		    << p0 << "  " << p1 << std::endl;
#endif

	  if( cid!=17 && cid!=18) 
	    {
	      ofp<<str<<std::endl;
	      continue;
	    }
	  int blcnum;
	  if(cid==17) blcnum=1;
	  else if(cid==18) blcnum=2;

	  if(wire==0)
	    {

	      ofp << setprecision(7); 
	      ofp << setiosflags(ios::scientific);      
#if 0
	      TH2F *h2d=(TH2F*)f->Get(Form("BLC%ddt_resid%d",blcnum,layer) );	      
	      if(h2d->GetEntries()<100 ) 
		{
		  ofp<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
		     <<p0<<"  "<<p1<<"  "<<p2<<"  "<<p3<<"  "<<p4<<"  "<<p5<<std::endl;
		  //	  fitp0[blcnum-1][layer-1][wire-1]=-1;
		  continue;
		}

	      TProfile *prof=h2d->ProfileX();
	      double range2=range+blcnum*10;
	      //     if(blcnum==1 && (layer==2 && layer==3 &&layer4)) range2=40;
	      prof->Fit("resl_ita","W","",5,range2);		  
	      prof->Fit("resl_ita","I","",5,range2);
	      
	      double par[6];
	      resl_ita->GetParameters(par);
	      double chi2,dof;
	      dof=resl_ita->GetNDF();
	      chi2=resl_ita->GetChisquare();
	      //		  double fac=0.8;
	      double fac=0.8;
	      cout<<"cid : l : w= "<<cid<<" : "<<layer<<" : "<<wire<<"  chi2 = "<<chi2/(double)dof<<endl;	 

	      //    if(par[0]<-1e-4) par[0]=1e-4;
	      if(chi2/(double)dof<5000)
		{
		  //		  fitp0[blcnum-1][layer-1][wire-1]=par[0];
		  ofp<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
		     <<p0-par[0]*fac<<"  "<<p1-par[1]*fac<<"  "<<p2-par[2]*fac<<"  "<<p3-par[3]*fac<<"  "<<p4-par[4]*fac<<"  "<<p5-par[5]*fac<<std::endl;
		}
	      else
		{
		  //		  fitp0[blcnum-1][layer-1][wire-1]=-1;
		  ofp<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
		     <<p0<<"  "<<p1<<"  "<<p2<<"  "<<p3<<"  "<<p4<<"  "<<p5<<std::endl;
		  
		}
#endif
#if 1	      
	      ofp << setprecision(7); 
	      ofp << setiosflags(ios::scientific);      
		for(int iw=1;iw<=32;iw++)
		  {
		    
		    TH2F *h2d=(TH2F*)f->Get(Form("BLC%ddt_resid%d_%d",blcnum,layer,iw) );
		    if(h2d->GetEntries()<100 ) 
		      {
			ofp<<cid<<"  "<<layer<<"  "<<iw<<"  "<<npar<<"  "
			   <<p0<<"  "<<p1<<"  "<<p2<<"  "<<p3<<"  "<<p4<<"  "<<p5<<std::endl;
			fitp0[blcnum-1][layer-1][wire-1]=-1;
			continue;
		      }
		    
		    
		    TProfile *prof=h2d->ProfileX();
		    double range2=range+blcnum*10;
		  // if(blcnum==1 && (layer==2 && layer==3 &&layer4)) range2=40;
		    prof->Fit("resl_ita","W","",5,range2);		  
		    prof->Fit("resl_ita","I","",5,range2);
		    
		    double par[6];
		    resl_ita->GetParameters(par);
		    double chi2,dof;
		    dof=resl_ita->GetNDF();
		    chi2=resl_ita->GetChisquare();
		    //		  double fac=0.8;
		    double fac=0.8;
		    cout<<"cid : l : w= "<<cid<<" : "<<layer<<" : "<<wire<<"  chi2 = "<<chi2/(double)dof<<endl;	 
		    
		    //  if(par[0]<-1e-4) par[0]=1e-4;
		    if(chi2/(double)dof<500)
		      {
			//		fitp0[blcnum-1][layer-1][wire-1]=par[0];
			ofp<<cid<<"  "<<layer<<"  "<<iw<<"  "<<npar<<"  "
			   <<p0-par[0]*fac<<"  "<<p1-par[1]*fac<<"  "<<p2-par[2]*fac<<"  "<<p3-par[3]*fac<<"  "<<p4-par[4]*fac<<"  "<<p5-par[5]*fac<<std::endl;
		      }
		    else
		      {
			//			fitp0[blcnum-1][layer-1][wire-1]=-1;
			ofp<<cid<<"  "<<layer<<"  "<<iw<<"  "<<npar<<"  "
			   <<p0<<"  "<<p1<<"  "<<p2<<"  "<<p3<<"  "<<p4<<"  "<<p5<<std::endl;
			
		      }
		  }
#endif	      
	    }
	  else
	    {
	      TH2F *h2d=(TH2F*)f->Get(Form("BLC%ddt_resid%d_%d",blcnum,layer,wire) );
	      if(h2d->GetEntries()==0 )
		{
 		  fitp0[blcnum-1][layer-1][wire-1]=-1;
 		  ofp<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
 			 <<p0<<"  "<<p1<<"  "<<p2<<"  "<<p3<<"  "<<p4<<"  "<<p5<<std::endl;
		  continue;
		}
 
		else if(h2d->GetEntries()<100 ) 
		{
		  h2d=(TH2F*)f->Get(Form("BLC%ddt_resid%d_16",blcnum,layer) );

		}
	      TProfile *prof=h2d->ProfileX();
	      double range2=range+blcnum*10;
	      if(blcnum==1 && (layer==2 && layer==3 &&layer4)) range2=40;
	      prof->Fit("resl_ita","W","",5,range2);		  
	      prof->Fit("resl_ita","I","",5,range2);
	      
	      double par[6];
	      resl_ita->GetParameters(par);
	      double chi2,dof;
	      dof=resl_ita->GetNDF();
	      chi2=resl_ita->GetChisquare();
	      double fac=0.8;
	      cout<<"cid : l : w= "<<cid<<" : "<<layer<<" : "<<wire<<"  chi2 = "<<chi2/(double)dof<<endl;

	      //	      if(par[0]<-1e-4) par[0]=1e-4;
	      if(chi2/(double)dof<50)
		{
		  fitp0[blcnum-1][layer-1][wire-1]=par[0];
		  cout<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
		     <<p0-par[0]*fac<<"  "<<p1-par[1]*fac<<"  "<<p2-par[2]*fac<<"  "<<p3-par[3]*fac<<"  "<<p4-par[4]*fac<<"  "<<p5-par[5]*fac<<std::endl;
		  ofp<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
		     <<p0-par[0]*fac<<"  "<<p1-par[1]*fac<<"  "<<p2-par[2]*fac<<"  "<<p3-par[3]*fac<<"  "<<p4-par[4]*fac<<"  "<<p5-par[5]*fac<<std::endl;
		}
	      else
		{
		  //		  fitp0[blcnum-1][layer-1][wire-1]=-1;
		  ofp<<cid<<"  "<<layer<<"  "<<wire<<"  "<<npar<<"  "
		     <<p0<<"  "<<p1<<"  "<<p2<<"  "<<p3<<"  "<<p4<<"  "<<p5<<std::endl;
		      
		}
	      cout<<"test hit fin"<<endl;
	    }
	}
      cout<<"test file fin"<<endl;
    }

  ofp.close();  
  fp.close();  

  cout<<"test fin"<<endl;
  return;
}
