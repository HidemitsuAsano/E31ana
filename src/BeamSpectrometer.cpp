// BeamSpectrometer.cpp
#include <iostream>
#include "BeamSpectrometer.h"

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
  double tmp1=TMath::Tan(in[1]/1000.);
  double tmp3=TMath::Tan(in[3]/1000.);

  TVector3 pos(in[0]-tmp1*gBLC1GZ-gBLC1GX,in[2]-tmp3*gBLC1GZ-gBLC1GY,0);
  pos.RotateY(-gBLC1GdY);
  pos.RotateX(-gBLC1GdX);
  TVector3 dir(tmp1,tmp3,1);
  dir.RotateY(-gBLC1GdY);
  dir.RotateX(-gBLC1GdX);

  out[1]=dir.X()/dir.Z();              //dx
  out[0]=pos.X()-out[1]*pos.Z(); //x
  out[3]=dir.Y()/dir.Z();              //dy
  out[2]=pos.Y()-out[3]*pos.Z(); //y
  out[4]=0.;
  out[5]=in[4];

  out[1]=TMath::ATan(out[1])*1000.;
  out[3]=TMath::ATan(out[3])*1000.;
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
  parblc2[1]=TMath::Tan(parblc2[1]/1000.);
  parblc2[3]=TMath::Tan(parblc2[3]/1000.);

  TVector3 pos(parblc2[0]-parblc2[1]*gBLC2GZ-gBLC2GX,parblc2[2]-parblc2[3]*gBLC2GZ-gBLC2GY,0);
  pos.RotateY(-gBLC2GdY);
  pos.RotateX(-gBLC2GdX);
  TVector3 dir(parblc2[1],parblc2[3],1);
  dir.RotateY(-gBLC2GdY);
  dir.RotateX(-gBLC2GdX);

  parblc2[1]=dir.X()/dir.Z();              //dx
  parblc2[0]=pos.X()-parblc2[1]*pos.Z(); //x
  parblc2[3]=dir.Y()/dir.Z();              //dy
  parblc2[2]=pos.Y()-parblc2[3]*pos.Z(); //y

  parblc2[1]=TMath::ATan(parblc2[1])*1000;
  parblc2[3]=TMath::ATan(parblc2[3])*1000;
}

static void fcn( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag )
{
  // TObject *obj=(TObject*)gMinuit->GetObjectFit();
  // obj->Print();
  // std::cout<<"fcn"<<std::endl;
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
  //  std::cout<<f<<std::endl;
}

// --------------------------------------------------------------//
BeamSpectrometer::BeamSpectrometer()
{
  Clear();
  minuit = new TMinuit(5);
  TROOT minexam("BeamSpectrometer","helix fit using TMinuit");
}

BeamSpectrometer::BeamSpectrometer(ConfMan* conf)
{
  Clear();  
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
  delete minuit;
}

void BeamSpectrometer::SetLocalTrack(const LocalTrack &blc1,const LocalTrack &blc2)
{
  BLC1Track=LocalTrack(blc1);
  BLC2Track=LocalTrack(blc2);
  CalcInitPar();
}

void BeamSpectrometer::Clear()
{
  BLC1Track.Clear();
  BLC2Track.Clear();
  for( int i=0; i<5; i++ ){
    Par[i] = Err[i] = -999.;
  }

  BLC1NumOfHits = BLC2NumOfHits =-999;
  FitDof = FitStat = -999;
  FitChi2 = -999.;
  //  MomCenter=-999.;

  BLC1GX=BLC1GY=BLC1GZ=BLC1GdX=BLC1GdY=-999.;  
  BLC2GX=BLC2GY=BLC2GZ=BLC2GdX=BLC2GdY=-999.;  
  for( int i=0; i<MAX_NUM_OF_HITS_BEAM; i++ ){
    BLC1HitPos[i].SetXYZ(-999,-999,-999);
    BLC1Weight[i] = -999.;
    BLC1Rotation[i]=-999.;
    BLC2HitPos[i].SetXYZ(-999,-999,-999);
    BLC2Weight[i] = -999.;
    BLC2Rotation[i]=-999.;
  }

}

void BeamSpectrometer::GetBLC2fromBLC1(LocalTrack *blc1,double *parblc,const double &momentum)
{
}
bool BeamSpectrometer::CalcInitPar()
{
  double a,b,c,d,tmp;
  blc1track()->gabc(tmp,b,a);
  blc1track()->gdef(tmp,d,c);
  double a2,b2,c2,d2;
  blc2track()->gabc(tmp,b2,a2);
  blc2track()->gdef(tmp,d2,c2);

  double x,y;
  blc1track()->XYPosatZ(0,x,y);
  double x2,y2;
  blc2track()->XYPosatZ(-130,x2,y2);

  Par[0]=x; //cm
  Par[1]=TMath::ATan(-b)*1000; //mrad
  Par[2]=y; //cm
  Par[3]=TMath::ATan(-d)*1000; //mrad
  double cang=-1./12.5;
  double cplus=3.75/12.5;
  double cminus=-0.11/12.5;  
  double angle=-TMath::ATan(-b)*1000.+TMath::ATan(-b2)*1000.;
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

void BeamSpectrometer::TMinuitFit(LocalTrack *blc1,LocalTrack *blc2, ConfMan* conf)
{
  SetLocalTrack(*blc1,*blc2);
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
      double temptheta,tempx,tempz,tmp;
      hit->gwpos(tempx,temptheta,tempz,1,0);           
      unsigned int lr =hit->leftright();
      if( lr==0 ) tempx += dl;
      else tempx -= dl;
      BLC1HitPos[i].SetXYZ( tempx,0,tempz ); //e in cm
      BLC1Weight[i] = 0.02;
      conf->GetReslMapManager()->GetParam(hit->cid(),hit->layer(),0,BLC1Weight[i],tmp);
      //conf->GetReslMapManager()->GetResolution(CID_BLC,Weight)
      BLC1Rotation[i]= temptheta*Deg2Rad;
      BLC1GX=hit->gx();
      BLC1GY=hit->gy();
      BLC1GZ=hit->gz();
      BLC1GdX=hit->dgx()*Deg2Rad; // deg -> rad
      BLC1GdY=hit->dgy()*Deg2Rad;
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
      double temptheta,tempx,tempz,tmp;
      hit->gwpos(tempx,temptheta,tempz,1,0);           
      unsigned int lr =hit->leftright();
      if( lr==0 ) tempx += dl;
      else tempx-= dl;

      BLC2HitPos[i].SetXYZ( tempx,0,tempz ); //e in cm
      BLC2Weight[i] = 0.02;
      conf->GetReslMapManager()->GetParam(hit->cid(),hit->layer(),0,BLC2Weight[i],tmp);
      //conf->GetReslMapManager()->GetResolution(CID_BLC2,Weight)
      BLC2Rotation[i]= temptheta*Deg2Rad;
      BLC2GX=hit->gx();
      BLC2GY=hit->gy();
      BLC2GZ=hit->gz()+130.;
      BLC2GdX=hit->dgx()*Deg2Rad; // deg -> rad // rotation around x axis
      BLC2GdY=hit->dgy()*Deg2Rad;// rotation around y axis
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
  //minuit->mnexcm("SET ERR", arglist,1,ierflg);
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
  //  std::cout<<minuit->fAmin<<"  "<<FitChi2<<std::endl;
#if 0
  double a,b,c,d,tmp;
  blctrack()->gabc(tmp,b,a);
  blctrack()->gdef(tmp,d,c);
  double a1,b1,c1,d1;
  double parblc[6];
  CalcParBLC1toBLC2(Par,parblc);
  blc1track()->gabc(tmp,b1,a1);
  blc1track()->gdef(tmp,d1,c1);
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
    parblc2[1]=TMath::Tan(parblc2[1]/1000.);
    parblc2[3]=TMath::Tan(parblc2[3]/1000.);

    TVector3 pos(parblc2[0]-parblc2[1]*BLC2GZ-BLC2GX,parblc2[2]-parblc2[3]*BLC2GZ-BLC2GY,0);
    pos.RotateY(-BLC2GdY);
    pos.RotateX(-BLC2GdX);
    TVector3 dir(parblc2[1],parblc2[3],1);
    dir.RotateY(-BLC2GdY);
    dir.RotateX(-BLC2GdX);
    
    parblc2[1]=dir.X()/dir.Z();              //dx
    parblc2[0]=pos.X()-parblc2[1]*pos.Z(); //x
    parblc2[3]=dir.Y()/dir.Z();              //dy
    parblc2[2]=pos.Y()-parblc2[3]*pos.Z(); //y

    parblc2[1]=TMath::ATan(parblc2[1])*1000;
    parblc2[3]=TMath::ATan(parblc2[3])*1000;

  }
}
void BeamSpectrometer::CalcBLC1GtoL( Double_t *in, Double_t *out ){
  in[1]=TMath::Tan(in[1]/1000.);
  in[3]=TMath::Tan(in[3]/1000.);

  TVector3 pos(in[0]-in[1]*BLC1GZ-BLC1GX,in[2]-in[3]*BLC1GZ-BLC1GY,0);
  pos.RotateY(-BLC1GdY);
  pos.RotateX(-BLC1GdX);
  TVector3 dir(in[1],in[3],1);
  dir.RotateY(-BLC1GdY);
  dir.RotateX(-BLC1GdX);

  out[1]=dir.X()/dir.Z();              //dx
  out[0]=pos.X()-out[1]*pos.Z(); //x
  out[3]=dir.Y()/dir.Z();              //dy
  out[2]=pos.Y()-out[3]*pos.Z(); //y

  out[4]=in[4];
  out[5]=in[5];
  //  std::cout<<"blc1gtol:\t"<<in[4]<<"\t"<<out[5]<<std::endl;

  in[1]=TMath::ATan(in[1])*1000.;
  in[3]=TMath::ATan(in[3])*1000.;
  out[1]=TMath::ATan(out[1])*1000.;
  out[3]=TMath::ATan(out[3])*1000.;
  //  for(int i=0;i<6;i++) std::cout<<in[i]<<"  "; std::cout<<std::endl;
  //  for(int i=0;i<6;i++) std::cout<<out[i]<<"  "; std::cout<<std::endl;
}

void BeamSpectrometer::CalcChi2()
{  
  double chisq=0.;
  int dof = 0;
#if 0
  std::cout<<"before : "<< blc1track()->GetCalcChisquare()<< "  "<< blc2track()->GetCalcChisquare()<<std::endl;
  double a,b,c,d,tmp;
  blc1track()->gabc(tmp,a,b);
  blc1track()->gdef(tmp,c,d);
  std::cout<<"beforeglobal1: "<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<std::endl;
  blc1track()->abc(tmp,a,b);
  blc1track()->def(tmp,c,d);
  std::cout<<"beforelocal1: "<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<std::endl;
#endif
  blc1track()->SetGParam(Par[0],TMath::Tan(Par[1]/1000.),Par[2],TMath::Tan(Par[3]/1000.));
  
  double parorg[6]={Par[0],Par[1],Par[2],Par[3],0.,Par[4]};
  double parblc2[6];
  double parblc1[6];
  CalcBLC1GtoL(parorg,parblc1);
  CalcParBLC1toBLC2(parorg,parblc2,0);
  blc2track()->SetGParam(parblc2[0],TMath::Tan(parblc2[1]/1000.),parblc2[2],TMath::Tan(parblc2[3]/1000.));
  CalcParBLC1toBLC2(parorg,parblc2,1);

#if 0
  std::cout<<"fit1:  "<<parblc1[0]<<"\t"<<TMath::Tan(parblc1[1]/1000.)<<"\t"
	   <<parblc1[2]<<"\t"<<TMath::Tan(parblc1[3]/1000.)<<std::endl;
  std::cout<<"fit2:  "<<parblc2[0]<<"\t"<<TMath::Tan(parblc2[1]/1000.)<<"\t"
	   <<parblc2[2]<<"\t"<<TMath::Tan(parblc2[3]/1000.)<<std::endl;
  //  double a,b,c,d,tmp;
  blc1track()->gabc(tmp,a,b);
  blc1track()->gdef(tmp,c,d);
  std::cout<<"afterglobal: "<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<std::endl;
  blc2track()->abc(tmp,a,b);
  blc2track()->def(tmp,c,d);
  std::cout<<"afterlocal: "<<a<<"\t"<<b<<"\t"<<c<<"\t"<<d<<std::endl;
#endif

  for( Int_t i=0; i<BLC1NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = BLC1HitPos[i].x();
    fitx = parblc1[0]+TMath::Tan(parblc1[1]/1000.)*BLC1HitPos[i].z();
    fity = parblc1[2]+TMath::Tan(parblc1[3]/1000.)*BLC1HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-BLC1Rotation[i]);
    double tmp=( (hitx-vec2.X())*(hitx-vec2.X()) )/BLC1Weight[i]/BLC1Weight[i];
    chisq +=tmp;
    //    std::cout<<"BLC1 "<<i<<"  "<<hitx<<"  "<<vec2.X()<<"  "<<tmp<<std::endl;
    dof++;
  }
  for( Int_t i=0; i<BLC2NumOfHits; i++ ){
    Double_t hitx, fitx, fity;
    hitx = BLC2HitPos[i].x();
    fitx = parblc2[0]+TMath::Tan(parblc2[1]/1000.)*BLC2HitPos[i].z();
    fity = parblc2[2]+TMath::Tan(parblc2[3]/1000.)*BLC2HitPos[i].z();
    TVector2 vec(fitx,fity);
    TVector2 vec2=vec.Rotate(-BLC2Rotation[i]);
    double tmp=( (hitx-vec2.X())*(hitx-vec2.X()) )/BLC2Weight[i]/BLC2Weight[i];
    chisq +=tmp;
    //    std::cout<<"BLC2 "<<i<<"  "<<hitx<<"  "<<vec2.X()<<"  "<<tmp<<std::endl;
    dof++;
  }

  FitDof = dof - 5;
  FitChi2 = chisq / FitDof;

  blc1track()->CalcHitPosition(true);
  blc1track()->CalcResidual(true);
  blc2track()->CalcHitPosition(true);
  blc2track()->CalcResidual(true);
#if 0
  std::cout<<"after : "<< blc1track()->GetCalcChisquare()<< "  "<< blc2track()->GetCalcChisquare()<<std::endl;
  std::cout<<"before: "<< blc1track()->chi2all()<< "  "<< blc2track()->chi2all()<<std::endl;
  // double blc1chi2=blc1track()->GetCalcChisquare();
  // double blc1dof= blc1track()->dof();
  // double blc2chi2=blc2track()->GetCalcChisquare();
  // double blc2dof= blc2track()->dof();
  std::cout<<"mom: "<<mom()<<std::endl;
  std::cout<<"total: "<< FitChi2 <<"\t"<<FitDof<<"\t"<<chisq<<std::endl;
  //  std::cout<<"blc1: "<<  blc1chi2 <<"\t"<< blc1dof <<"\t"<< blc1chi2*blc1dof <<std::endl;
  //  std::cout<<"blc2: "<<  blc2chi2 <<"\t"<< blc2dof <<"\t"<< blc2chi2*blc2dof <<std::endl;
#endif
}

