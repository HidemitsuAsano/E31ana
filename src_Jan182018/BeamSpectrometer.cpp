// BeamSpectrometer.cpp
#include <iostream>
#include "BeamSpectrometer.h"
#include "TMath.h"

#define DEBUG 0
#define DEBUG2 0
ClassImp(BeamSpectrometer);


// parameters for TMinuit
static const Double_t  FitStep[5] = { 1e-10, 1e-10, 1e-10,1e-10,1e-5 };
static const Double_t LowLimit[5] = { -20, -100, -20,-100, -10 };
static const Double_t  UpLimit[5] = { 20,  100, 20 ,100, 10 };

// global variables for TMinuit
static Int_t gBLC1NumOfHits;
static TVector3 gBLC1HitPos[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC1Weight[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC1Rotation[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC1GX,gBLC1GY,gBLC1GZ,gBLC1GdX,gBLC1GdY;

static Int_t gBLC2NumOfHits;
static TVector3 gBLC2HitPos[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC2Weight[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC2Rotation[MAX_NUM_OF_HITS_BEAM];
static Double_t gBLC2GX,gBLC2GY,gBLC2GZ,gBLC2GdX,gBLC2GdY;

static TMatrixD gD5Matrix1st;
static TMatrixD gD5Matrix2nd[5];

// --------------------------------------------------------------//
// functions for TMinuit
static void BLC1GtoL( Double_t *in, Double_t *out ){
  out[0]=in[0]-gBLC1GX-in[1]*gBLC1GZ;
  out[1]=in[1]-gBLC1GdX;
  out[2]=in[2]-gBLC1GY-in[3]*gBLC1GZ;
  out[3]=in[3]-gBLC1GdY;
  out[4]=0.;
  out[5]=in[4];
  //  std::cout<<"blc1gtol:\t"<<in[4]<<"\t"<<out[5]<<std::endl;
}
static void ParBLC1toBLC2( Double_t *parblc1, Double_t *parblc2){
  double parin[6]={parblc1[0],parblc1[1],parblc1[2],parblc1[3],0.,parblc1[4]};
  TMatrixD in;
  in.Use(6,1,parin);
  //  std::cout<<"blc1toblc:\t"<<parblc1[4]<<"\t"<<in[5][0]<<std::endl;
  TMatrixD out1st;
  out1st.ResizeTo(6,1);
  //  gD5Matrix1st.Print();
  //  in.Print();
  out1st.Mult(gD5Matrix1st,in);

  for(int i=0;i<6;i++){
    double out2nd=0;
    if(i<5)
      for(int j=0;j<6;j++)
	for(int k=0;k<6;k++)
	  out2nd+=gD5Matrix2nd[i][j][k]*in[j][0]*in[k][0];
    parblc2[i]=out1st[i][0]+out2nd;
  }
  parblc2[0]=parblc2[0]-gBLC2GX-parblc2[1]*gBLC2GZ;
  parblc2[1]=parblc2[1]-gBLC2GdX;
  parblc2[2]=parblc2[2]-gBLC2GY-parblc2[3]*gBLC2GZ;
  parblc2[3]=parblc2[3]-gBLC2GdY;
}
static void fcn( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  Double_t chisq=0.;  Int_t dof = 0;
  double parblc1[6];
  BLC1GtoL(par,parblc1);

  for( Int_t i=0; i<gBLC1NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = gBLC1HitPos[i].x();
    fitx = parblc1[0]+TMath::Tan(parblc1[1]/1000.)*gBLC1HitPos[i].z();
    fity = parblc1[2]+TMath::Tan(parblc1[3]/1000.)*gBLC1HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-gBLC1Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/gBLC1Weight[i]/gBLC1Weight[i];
    dof++;
  }

  Double_t parblc2[6];
  ParBLC1toBLC2(par,parblc2);

  for( Int_t i=0; i<gBLC2NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = gBLC2HitPos[i].x();
    fitx = parblc2[0]+TMath::Tan(parblc2[1]/1000.)*gBLC2HitPos[i].z();
    fity = parblc2[2]+TMath::Tan(parblc2[3]/1000.)*gBLC2HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-gBLC2Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/gBLC2Weight[i]/gBLC2Weight[i];
    dof++;
  }
  f = chisq/(dof-5);
}

// --------------------------------------------------------------//
BeamSpectrometer::BeamSpectrometer()
{
  BLC1Track=NULL;
  BLC2Track=NULL;
  for( int i=0; i<5; i++ ){
    Par[i] = Err[i] = -999.;
  }

  BLC1NumOfHits = BLC2NumOfHits =-999;
  FitDof = FitStat = -999;
  FitChi2 = -999.;
  
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    BLC1HitPos[i].SetXYZ(-999,-999,-999);
    BLC1Weight[i] = -999.;
    BLC1Rotation[i]=-999.;
    BLC2HitPos[i].SetXYZ(-999,-999,-999);
    BLC2Weight[i] = -999.;
    BLC2Rotation[i]=-999.;
  }
  
  minuit = new TMinuit(5);
  TROOT minexam("BeamSpectrometer","helix fit using TMinuit");
}

BeamSpectrometer::BeamSpectrometer(ConfMan* conf)
{
  BLC1Track=NULL;
  BLC2Track=NULL;
  for( int i=0; i<5; i++ ) {
    Par[i] = -999.;
    Err[i] = -999.;
  }
  FitDof = FitStat = -999;
  FitChi2 = -999.;
  
  minuit = new TMinuit(5);
  TROOT minexam("BeamSpectrometer","Beam fit using TMinuit");
  
  MomCenter=conf->GetTransferMatrixManager()->GetCentralMomentum();

  D5Matrix1st.ResizeTo(6,6,-1);
  const double *ele = conf->GetTransferMatrixManager()->GetD5Matrix();
  D5Matrix1st.Use(6,6,ele);
  for(int i=0; i<5 ;i++){
    D5Matrix2nd[i].ResizeTo(6,6,-1);
    ele = conf->GetTransferMatrixManager()->GetD5Matrix2nd(i);
    D5Matrix2nd[i].Use(6,6,ele);
  }
}

BeamSpectrometer::~BeamSpectrometer()
{
  if(BLC1Track) delete BLC1Track;
  if(BLC2Track) delete BLC2Track;
  delete minuit;
}

void BeamSpectrometer::SetLocalTrack(LinearTrack *blc1,LinearTrack* blc2)
{
  BLC1Track=new LinearTrack(*blc1);
  BLC2Track=new LinearTrack(*blc2);
  CalcInitPar();
}

void BeamSpectrometer::Clear()
{
  if(BLC1Track) delete BLC1Track;
  BLC1Track = NULL;
  if(BLC2Track) delete BLC2Track;
  BLC2Track = NULL;
}

void BeamSpectrometer::GetBLC2fromBLC1(LocalTrack *blc1,double *parblc,const double &momentum)
{
}
bool BeamSpectrometer::CalcInitPar()
{
  double a,b,c,d;
  blc1track()->gabcd(a,b,c,d);
  double a2,b2,c2,d2;
  blc2track()->gabcd(a2,b2,c2,d2);

  double x,y;
  blc1track()->XYPosatZ(0,x,y);
  double x2,y2;
  blc2track()->XYPosatZ(-130,x2,y2);

  Par[0]=x; //cm
  Par[1]=TMath::ATan(b)*1000; //mrad
  Par[2]=y; //cm
  Par[3]=TMath::ATan(d)*1000; //mrad
  double cang=-1./12.5;
  double cplus=3.75/12.5;
  double cminus=-0.11/12.5;  
  double angle=-TMath::ATan(b)*1000.+TMath::ATan(b2)*1000.;
  double cmom=cang*angle+cplus*(x+x2)+cminus*(x2-x);
  cmom *= -1;
  //  std::cout<<"cmom:\t"<<cmom<<std::endl;
  Par[4]=cmom;
  return true;
}
void BeamSpectrometer::SetParameters( const double *param )
{
  for( int i=0; i<5; i++ ) Par[i] = (Double_t)param[i];
}

void BeamSpectrometer::GetParameters( double *param )
{
  for( int i=0; i<5; i++ ) param[i] = (double)Par[i];
}

void BeamSpectrometer::SetGlobalVariables()
{
  gBLC1NumOfHits = BLC1NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    gBLC1HitPos[i] = BLC1HitPos[i];
    gBLC1Weight[i] = BLC1Weight[i];
    gBLC1Rotation[i]= BLC1Rotation[i];
  }
  gBLC2NumOfHits = BLC2NumOfHits;
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    gBLC2HitPos[i] = BLC2HitPos[i];
    gBLC2Weight[i] = BLC2Weight[i];
    gBLC2Rotation[i]= BLC2Rotation[i];
  }

  gBLC1GX=BLC1GX;
  gBLC1GY=BLC1GY;
  gBLC1GZ=BLC1GZ;
  gBLC1GdX=BLC1GdX;
  gBLC1GdY=BLC1GdY;
  
  gBLC2GX=BLC2GX;
  gBLC2GY=BLC2GY;
  gBLC2GZ=BLC2GZ;
  gBLC2GdX=BLC2GdX;
  gBLC2GdY=BLC2GdY;
  
  gD5Matrix1st.Use(6,6,D5Matrix1st.GetMatrixArray());
  for(int i=0;i<5;i++)
    gD5Matrix2nd[i].Use(6,6,D5Matrix2nd[i].GetMatrixArray());
}

void BeamSpectrometer::TMinuitFit(LinearTrack* blc1,LinearTrack* blc2, ConfMan* conf)
{
#if DEBUG2
  std::cout<<"chi2all: "<<blc1->chi2all()<<"\t"<<blc2->chi2all()<<std::endl;
#endif
  SetLocalTrack(blc1,blc2);
  fit(conf);
}
void BeamSpectrometer::fit(ConfMan* conf)
{
#if 0
  std::cout << "!!! BeamSpectrometer::fit() !!!" << std::endl;
#endif

  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ) {
    if(i<blc1track()->nhit()) {
      ChamberLikeHit* hit=blc1track()->hit(i);
      double dl = hit->dl();
      double temptheta,tempx,tempz;
      hit->gwpos(tempx,temptheta,tempz,1,0);           
      unsigned int lr =hit->leftright();
      if( lr==0 ) tempx += dl;
      else tempx-= dl;
      BLC1HitPos[i].SetXYZ( tempx,0,tempz ); //e in cm
      BLC1Weight[i] = 0.02;
      //conf->GetReslMapManager()->GetResolution(CID_BLC,Weight)
      BLC1Rotation[i]= temptheta/180.*TMath::Pi();
      BLC1GX=hit->gx();
      BLC1GY=hit->gy();
      BLC1GZ=hit->gz();
      BLC1GdX=TMath::ATan(hit->dgx())*1000.;
      BLC1GdY=TMath::ATan(hit->dgy())*1000.;
    }
    else {
      BLC1HitPos[i].SetXYZ(-999,-999,-999);
      BLC1Weight[i] = -999.;
      BLC1Rotation[i]=-999.;
    }
  }
  
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ) {
    if(i<blc2track()->nhit()) {
      ChamberLikeHit* hit=blc2track()->hit(i);
      double dl = hit->dl();
      double temptheta,tempx,tempz;
      hit->gwpos(tempx,temptheta,tempz,1,0);           
      unsigned int lr =hit->leftright();
      if( lr==0 ) tempx += dl;
      else tempx-= dl;

      BLC2HitPos[i].SetXYZ( tempx,0,tempz ); //e in cm
      BLC2Weight[i] = 0.02;
      //conf->GetReslMapManager()->GetResolution(CID_BLC2,Weight)
      BLC2Rotation[i]= temptheta/180.*TMath::Pi();
      BLC2GX=hit->gx();
      BLC2GY=hit->gy();
      BLC2GZ=hit->gz()+130.;
      BLC2GdX=TMath::ATan(hit->dgx())*1000.;
      BLC2GdY=TMath::ATan(hit->dgy())*1000.;
    }
    else {
      BLC2HitPos[i].SetXYZ(-999,-999,-999);
      BLC2Weight[i] = -999.;
      BLC2Rotation[i]=-999.;
    }
  }
  
  BLC1NumOfHits=blc1track()->nhit();
  BLC2NumOfHits=blc2track()->nhit();  
  
  Int_t plevel=-1;
  //  Int_t plevel=1;
  SetGlobalVariables();
  minuit->SetPrintLevel( plevel );
  minuit->SetFCN( fcn );

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  //  minuit->mnexcm("SET ERR", arglist,1,ierflg);
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  // Set starting values and step sizes for parameters
  TString name[5] = {"x_blc1", "theta_blc1", "y_blc1","phi_blc1","dp"};
  for( Int_t i=0; i<5; i++){
    minuit->mnparm(i, name[i],Par[i],FitStep[i],LowLimit[i],UpLimit[i],ierflg);
  }
  //  for( Int_t i=0; i<4; i++)
  //    minuit->FixParameter(i);

  minuit->Command("SET STRategy 0");
  // Now ready for minimization step
  arglist[0] = 1000;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  //  for( Int_t i=0; i<4; i++)
  //    minuit->Release(i);
  //  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit->mnstat(amin, edm,  errdef, nvpar, nparx, icstat);
#if DEBUG
  minuit->mnprin(5,amin);
#endif
  
  Int_t err;
  Double_t bnd1, bnd2;
  for( Int_t i=0; i<5; i++ ){
    minuit->mnpout(i, name[i], Par[i], Err[i], bnd1, bnd2, err);
  }  
  FitStat = icstat;
  CalcChi2();
#if 0
  double a,b,c,d;
  blctrack()->gabcd(a,b,c,d);
  double a1,b1,c1,d1;
  double parblc[6];
  CalcParBLC1toBLC2(Par,parblc);
  blc1track()->gabcd(a1,b1,c1,d1);
  std::cout<<a1<<"\t"<<b1<<"\t"<<c1<<"\t"<<d1<<
    TMath::Tan(Par[1]/1000.)<<"\t"<<Par[0]<<std::endl;
  std::cout<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<
    TMath::Tan(parblc[1]/1000.)<<"\t"<<parblc[0]<<std::endl;
#endif
}

void BeamSpectrometer::CalcParBLC1toBLC2( double *parblc1, double *parblc2, const bool &GTOL){
  TMatrixD in;
  in.Use(6,1,parblc1);
  TMatrixD out1st;
  out1st.ResizeTo(6,1);
  out1st.Mult(D5Matrix1st,in);
  for(int i=0;i<6;i++){
    double out2nd=0;
    if(i<5)
      for(int j=0;j<6;j++)
	for(int k=0;k<6;k++)
	  out2nd+=D5Matrix2nd[i][j][k]*in[j][0]*in[k][0];
    parblc2[i]=out1st[i][0]+out2nd;
  }
  if(GTOL){
    parblc2[0]-=BLC2GX;
    parblc2[1]-=BLC2GdX;
    parblc2[2]-=BLC2GY;
    parblc2[3]-=BLC2GdY;
  }
}
void BeamSpectrometer::CalcBLC1GtoL( Double_t *in, Double_t *out ){
  out[0]=in[0]-BLC1GX;
  out[1]=in[1]-BLC1GdX;
  out[2]=in[2]-BLC1GY;
  out[3]=in[3]-BLC1GdY;
  out[4]=in[4];
  out[5]=in[5];
  //  std::cout<<"blc1gtol:\t"<<in[4]<<"\t"<<out[5]<<std::endl;
}
void BeamSpectrometer::CalcChi2()
{  
  double chisq=0.;
  int dof = 0;
  blc1track()->SetGABCD(Par[0],TMath::Tan(Par[1]/1000.),Par[2],TMath::Tan(Par[3]/1000.));
  
  double parorg[6]={Par[0],Par[1],Par[2],Par[3],0.,Par[4]};
  double parblc2[6];
  double parblc1[6];
  CalcBLC1GtoL(parorg,parblc1);
  CalcParBLC1toBLC2(parorg,parblc2,0);
  blc2track()->SetGABCD(parblc2[0],TMath::Tan(parblc2[1]/1000.),parblc2[2],TMath::Tan(parblc2[3]/1000.));
  CalcParBLC1toBLC2(parorg,parblc2,1);
#if DEBUG2
  std::cout<<"fit:  "<<parblc2[0]<<"\t"<<TMath::Tan(parblc2[1]/1000.)<<"\t"
	   <<parblc2[2]<<"\t"<<TMath::Tan(parblc2[3]/1000.)<<std::endl;
  double a,b,c,d;
  blc2track()->gabcd(a,b,c,d);
  std::cout<<"afterglobal: "<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<std::endl;
  blc2track()->abcd(a,b,c,d);
  std::cout<<"afterlocal: "<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<std::endl;
#endif

  for( Int_t i=0; i<BLC1NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = BLC1HitPos[i].x();
    fitx = parblc1[0]+TMath::Tan(parblc1[1]/1000.)*BLC1HitPos[i].z();
    fity = parblc1[2]+TMath::Tan(parblc1[3]/1000.)*BLC1HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-BLC1Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/BLC1Weight[i]/BLC1Weight[i];
    dof++;
  }
  for( Int_t i=0; i<BLC2NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = BLC2HitPos[i].x();
    fitx = parblc2[0]+TMath::Tan(parblc2[1]/1000.)*BLC2HitPos[i].z();
    fity = parblc2[2]+TMath::Tan(parblc2[3]/1000.)*BLC2HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-BLC2Rotation[i]);
    chisq +=( (hitx-vec2.X())*(hitx-vec2.X()) )/BLC2Weight[i]/BLC2Weight[i];
    dof++;
  }

  FitDof = dof - 5;
  FitChi2 = chisq / FitDof;

  blc1track()->CalcHitPosition();
  blc1track()->CalcResidual();
  blc2track()->CalcHitPosition();
  blc2track()->CalcResidual();

#if 0
  double blc1chi2=blc1track()->GetCalcChisquare();
  double blc1dof= blc1track()->dof();
  double blc2chi2=blc2track()->GetCalcChisquare();
  double blc2dof= blc2track()->dof();
  std::cout<<"mom: "<<mom()<<std::endl;
  std::cout<<"total: "<< FitChi2 <<"\t"<<FitDof<<"\t"<<chisq<<std::endl;
  std::cout<<"blc1: "<<  blc1chi2 <<"\t"<< blc1dof <<"\t"<< blc1chi2*blc1dof <<std::endl;
  std::cout<<"blc2: "<<  blc2chi2 <<"\t"<< blc2dof <<"\t"<< blc2chi2*blc2dof <<std::endl;
#endif
}

