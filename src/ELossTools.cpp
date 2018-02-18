// ELossTools.cpp
#include "ELossTools.h"
#define DEBUG 0
bool ELossTools::CalcHelixElossToVertex(const double param[5],const TVector3 &vertex, const double &momin, const double &mass, double &momout,double &tof){
  double rmin,rmax;
  double tmpmom;
  double tmptof;
  TString mat;
  momout=momin;
  tof=0;
  if(mass<=0) return false;
  //CDC Gas 
  tmpmom=momout;  rmax=15.1+(53.-15.1)/3.;  rmin=15.1;  mat="CDCGas";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //CDC CFRP 
  tmpmom=momout;  rmax=15.1;  rmin=15.0;  mat="CFRP";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //Air
  tmpmom=momout;  rmax=15.0;  rmin=14.05;  mat="Air";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //IH Scinti
  tmpmom=momout;  rmax=14.05;  rmin=13.75;  mat="Plastic";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //Air
  tmpmom=momout;  rmax=13.75;  rmin=7.5;  mat="Air";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //Target CFRP
  tmpmom=momout;  rmax=7.5;  rmin=7.4;  mat="CFRP";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //Target Heat Shield(Al)
  tmpmom=momout;  rmax=6.63;  rmin=6.6;  mat="Aluminum";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //Target Cell (Be)
  tmpmom=momout;  rmax=3.43;  rmin=3.4;  mat="Beryllium";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  //Target 3He
  tmpmom=momout;  rmax=3.4;  rmin=0.;  mat="LHelium-3";
  if(!ELossTools::CalcdE(param,vertex,rmax,rmin,tmpmom,mass,mat,momout,tmptof)) return false;
  tof+=tmptof;
  return true;
}

bool ELossTools::CalcHelixElossToVertexTGeo(const double param[5],const TVector3 &vertex, const double &momin, const double &mass, double &momout,double &tof,const double &rstart){
  //  int sign=1; // add energy loss
  TVector3 pos1=MathTools::CalcHelixPosatR(param,rstart);
  return ELossTools::CalcHelixElossPointToPointTGeo(param,pos1,vertex,momin,mass,momout,tof);
}

bool ELossTools::CalcHelixElossToNextBoundary(const double param[5],TVector3 &pos1,const double &momin, const double &mass,
					      double &momout,  double &tof, double &length,int &id, int sign)
{
  if(pos1.Mag()>70) return false;
  TString mat;
  TVector3 pos2;
  if(!GeomTools::HelixStepToNextVolume(param,pos1,pos2,length,mat,id)) return false;
  GeomTools::GetIDMat(pos1,mat);
  pos1=pos2;
  return ELossTools::CalcdE(momin,mass,length,mat,momout,sign,tof);
}

bool ELossTools::CalcElossBeamTGeo(const TVector3 &in, const TVector3 &out, const double &momin, const double &mass, double &momout,double &tof){
#if DEBUG
  std::cout<<"-- calcelossbeamtgeo --"<<std::endl;
#endif
  int sign=-1;
  double tmpmom=momin;
  TString mat;
  TVector3 tmppos=in;
  GeomTools::GetIDMat(in,mat);
  TString newmat;
  double length;
  double sumlength=0;
  tof=0;
  double tmptof=-999;
  int count=0;
  while(GeomTools::StepToNextVolume(tmppos,out,length,newmat)){
    if(!ELossTools::CalcdE(tmpmom,mass,length,mat,momout,sign,tmptof,false))      return false;
    tof+=tmptof;
#if DEBUG
    std::cout<<"  "<<mat<<"  "<<length<<"  "<<(momout-tmpmom)*1000<<"  "<<tmptof<<"  "<<tof<<std::endl;
#endif
    tmpmom=momout;
    mat=newmat;
    sumlength+=length;
    tmppos+=(out-tmppos).Unit()*length;
    count++;
    if(count>10000) break;
  }
  length=(out-tmppos).Mag();//-sumlength;
  if(!ELossTools::CalcdE(tmpmom,mass,length,mat,momout,sign,tmptof,false)) return false;
  tof+=tmptof;
#if DEBUG
  std::cout<<"  "<<mat<<"  "<<length<<"  "<<(momout-tmpmom)*1000<<"  "<<tmptof<<"  "<<tof<<std::endl;
#endif
  //  exit(-1);
  return true;
}

bool ELossTools::CalcElossForwardTGeo(const TVector3 &in, const TVector3 &fdcpos,const double &flength, const double &momin, const double &mass, double &momout,double &tof){
#if DEBUG
  std::cout<<"-- calcelossforwardtgeo --"<<std::endl;
#endif
  double tmpmom=momin;
  double tmptof;
  int sign=-1;
  if(!CalcElossBeamTGeo(in,fdcpos,tmpmom,mass,momout,tmptof))return false;
  tmpmom=momout;
  tof=tmptof;
  TString mat= "Air";
  double length=flength-(fdcpos-in).Mag();
  if(!ELossTools::CalcdE(tmpmom,mass,length,mat,momout,sign,tmptof,false))
    return false;
  tof+=tmptof;
#if DEBUG
  std::cout<<"  "<<mat<<"  "<<length<<"  "<<(momout-tmpmom)*1000<<"  "<<tmptof<<"  "<<tof<<std::endl;
#endif
  //  exit(-1);
  return true;
}

bool ELossTools::CalcHelixElossPointToPointTGeo(const double param[5],const TVector3 &posin, const TVector3 &posout, const double &momin, const double &mass, double &momout,double &tof,double offs){
#if DEBUG
  std::cout<<"-----------CalcHelixElossPointToPointTGeo------------"<<std::endl;
  posin.Print();
  posout.Print();
  std::cout<<momin<<"  "<<mass<<std::endl;
  std::cout<<"-----------"<<std::endl;
#endif
  momout=momin;  
  if(mass<=0){
#if DEBUG
    std::cout<<"calchelixelosspointtopointtgeo :: negative mass is not allowed !!!" << mass <<std::endl;
#endif
    return false;
  }
  double tmpmom=momin;
  TString mat,tmpmat;
  double defaultstep=1*cm;
  double step=defaultstep;
  double tmpstep=0;
  int sign=1; // add energy loss
  double phiin=MathTools::CalcHelixPhi(posin,param);
  double phiout=MathTools::CalcHelixPhi(posout,param);
  if(phiin*param[2]>phiout*param[2]) sign=-1;
  //  if(phiin<phiout) sign=-1;
#if DEBUG
  std::cout<<"phiin,phiout,param[2],sign  "<<phiin<<"  "<<phiout<<"  "<<param[2]<<"  "<<sign<<std::endl;
#endif 
  TVector3 pos1=posin;
  TVector3 pos2=MathTools::CalcHelixStep(param,pos1,sign*step);
  double tmptof,length=0;
  double phi1=MathTools::CalcHelixPhi(pos1,param);
  double phi2=MathTools::CalcHelixPhi(pos2,param);
  int count=0;
  int id1=GeomTools::GetID(pos1);
  int id2=GeomTools::GetID(pos2);
  int tmp;
  tof=0;
  while(1){
#if DEBUG
    std::cout<<"INIT: id1,id2,phi1,phi2,pos1  "<<id1<<"  "<<id2
	     <<"  "<<MathTools::CalcHelixPhi(pos1,param)
	     <<"  "<<MathTools::CalcHelixPhi(pos2,param)<<"  ";
    pos1.Print();
    pos2.Print();
#endif
    if(pos1.Mag()>70){

#if DEBUG
      std::cout<<"ELossTools:: to large rho "<<std::endl;
      pos1.Print();
      posin.Print();
      std::cout<<"  ->  "; posout.Print();
      std::cout<<phiin<<"  ->  "<<phiout<<std::endl;
#endif
      return false;
    }
    length=tmpstep;
    while(step>99*um&&(phi1-phiout)*(phi2-phiout)>0){
      if(GeomTools::IsSameVolume(pos1,pos2)){
	length+=step;
	pos1=pos2;
	pos2=MathTools::CalcHelixStep(param,pos1,sign*step);
	phi1=MathTools::CalcHelixPhi(pos1,param);
	phi2=MathTools::CalcHelixPhi(pos2,param);
      }else{
	step/=10.;
	pos2=MathTools::CalcHelixStep(param,pos1,sign*step);
	phi1=MathTools::CalcHelixPhi(pos1,param);
	phi2=MathTools::CalcHelixPhi(pos2,param);
      }    
      if(length>100*cm)return false;
    }

    if((phi1-phiout)*(phi2-phiout)<0){
#if DEBUG
      std::cout<<"FINISH: id1,id2,phi1,phi2,pos1  "<<id1<<"  "<<id2
	       <<"  "<<MathTools::CalcHelixPhi(pos1,param)
	       <<"  "<<MathTools::CalcHelixPhi(pos2,param)<<"  ";
      pos1.Print();
      pos2.Print();
#endif
      length+=MathTools::CalcHelixArc(param,pos1,posout);
      length-=offs;
      if(length<0) length=0.01;
      GeomTools::GetIDMat(pos1,mat);
      ELossTools::CalcdE(tmpmom,mass,length,mat,momout,sign,tmptof);
      //      std::cout<<mat<<"  "<<mass<<"  "<<length<<"  "<<tmpmom<<"  -> "<<momout<<std::endl;
      tof+=tmptof;
      //      std::cout<<"   "<<tmptof<<"  "<<tof<<std::endl;
      return true;
    }else{
#if DEBUG
      std::cout<<"NEXT VOL: id1,id2,phi1,phi2,pos1  "<<id1<<"  "<<id2
	       <<"  "<<MathTools::CalcHelixPhi(pos1,param)
	       <<"  "<<MathTools::CalcHelixPhi(pos2,param)<<"  ";
      pos1.Print();
      pos2.Print();
#endif
      step*=10;
      GeomTools::GetIDMat(pos1,mat);
      //      std::cout<<pos1.Perp()<<"   "<<pos1.Z()<<"  "<<mat<<std::endl;
      pos2=MathTools::CalcHelixStep(param,pos1,sign*step);
      step*=sign;
      GeomTools::CrossBoundary(pos1,pos2,step,tmpmat,tmpstep,tmp);
      length+=step;
      //      std::cout<<"length "<<length<<std::endl;
      ELossTools::CalcdE(tmpmom,mass,length,mat,momout,sign,tmptof);
      //      std::cout<<mat<<"  "<<mass<<"  "<<length<<"  "<<tmpmom<<"  -> "<<momout<<std::endl;
      tof+=tmptof;
      //      std::cout<<"   "<<tmptof<<"  "<<tof<<std::endl;
      tmpmom=momout;
      //    pos2.Print();
      step=defaultstep;
      pos1=pos2;
      pos2=MathTools::CalcHelixStep(param,pos1,sign*step);
      id1=GeomTools::GetID(pos1);
      id2=GeomTools::GetID(pos2);
      phi1=MathTools::CalcHelixPhi(pos1,param);
      phi2=MathTools::CalcHelixPhi(pos2,param);
      count++;
      if(count>50)return false;
    }
  }
  return true;
}


bool ELossTools::CalcdE(const double param[5], const TVector3 &vertex, const double &rmax, const double &rmin,
		       const double &momin, const double &mass , const TString &material, double &momout,double &tof)
{
  int sign=1;
  double rvtx=vertex.Perp();
  if(rvtx>rmax) return false;
  double length=0;
  if(rvtx>rmin){
    if(1){
      double anglein=MathTools::CalcHelixPhi(vertex,param);
      // TVector3 temp=MathTools::GetPosition(anglein,param);
      // vertex.Print();
      // temp.Print();
      // std::cout<<param[2]<<std::endl;
      if(anglein>TMath::Pi()) anglein-=TMath::Pi();
      double angleout=MathTools::CalcHelixPhiatR(param,rmax);
      length=TMath::Abs(1./param[2] *(angleout-anglein))*sqrt(1+param[4]*param[4]);
      //     std::cout<<"anglein,out"<<anglein<<"  "<<angleout<<std::endl;
    }
    //    length=MathTools::CalcHelixArc(param,rmax,rvtx)*sqrt(1+param[4]*param[4]);
  }else{
    length=MathTools::CalcHelixArc(param,rmax,rmin);
  }
  return ELossTools::CalcdE(momin,mass,length,material,momout,sign,tof);
}

bool ELossTools::CalcdE(const double &momin,const double &mass, const double &length,
			const TString &matname,double &momout,const int &sign,double &tof,bool CORR)
{
  // length in cm
  momout=momin;
  if(mass<=0) return false;
  double Ein=sqrt(mass*mass+momin*momin);
  double pm=1.;
  if(momin<0) pm=-1.;
  double beta=TMath::Abs(momin)/Ein;
  double step=0.05; //cm
  if(mass>0.5&&momin<0.4) step=0.01; 
  if(mass<0.5&&momin<0.2) step=0.01; 
  if(matname.Contains("Air")||matname.Contains("Gas")||matname.Contains("Vacuum")) step*=20;
  //  step=0.001;
  double Eout=0.;
  double eloss=0.;
  double toteloss=0.;
  double totlength=0.;
  tof=0;
  int i=0;
  double betamid=0.;
  for(i=1;i<=length/step;i++){
    eloss=ELossTools::CalcdEdX(beta,matname,CORR)*step;
    //    if(eloss) return false;
    Eout=Ein+sign*eloss;
    if(Eout<mass||beta<=0.0){
#if DEBUG
      std::cout<<"particle stops in " << matname <<" at "<<i*step<< " cm"<<std::endl;
#endif
      momout*=pm;
      return false;
    }    
    momout=sqrt(Eout*Eout-mass*mass);

    betamid=(beta+momout/Eout)*0.5;
    tof+=step/(Const*betamid*100);

    beta=momout/Eout;
    Ein=Eout;
    toteloss+=eloss;
    totlength+=step;
  }
  eloss=ELossTools::CalcdEdX(beta,matname,CORR);
  //  if(eloss) return false;
  eloss*=(length-totlength);
  //  totlength+=(length-step*i);
  Eout=Ein+sign*eloss;
  if(Eout<mass||beta<=0.0){
#if DEBUG
    std::cout<<"particle stops in " << matname <<" at "<<length<< " cm"<<std::endl;
#endif
    momout*=pm;
    return false;
  }
  momout=sqrt(Eout*Eout-mass*mass);

  betamid=(beta+momout/Eout)*0.5;
  tof+=(length-totlength)/(Const*betamid*100);

  toteloss+=eloss;
  momout*=pm;

  
#if DEBUG
  //#if 1
  std::cout<<matname<<"  mass "<<mass<<"  totlength/length "<<totlength<<" / "<<length<<"   "<<momin<<" -> "<<momout<<" eloss = "<<toteloss/MeV<<" time "<<tof<<std::endl;
#endif
  return true;  
}

double ELossTools::CalcdEdX(const double &beta,const TString &matname, bool CORR){
  double rho,I,Z_A;
  double Z,C0,a,M,X1,X0;
  // rho, A,Z from knucl3_ag: KnuclMaterials
  // I,C0,a,M,X1,X0 from Atomic Data and Nuclear Data Tables 30, 261-271 (1984)
  if( !matname.CompareTo("Aluminum") ){
    rho=2.7;//g/cm3
    I=166.;// eV
    Z_A=0.48181;
    C0=-4.2395;
    a=0.08024;
    M=3.6345;
    X1=3.0127;
    X0=0.1708;
    Z=13.;
  }
  else if(!matname.CompareTo("CDCGas") ){
    // Ar-Ethane 50-50, simply take avarage
    rho=(1.356*0.5+1.782*0.5)/1000.;
    I=(188.+45.4)/2.;
    Z_A=(0.45059+0.59861)/2.;  
    C0=-(11.9480+9.1043)/2.;
    a=(0.19714+0.09627)/2.;
    M=(2.9618+3.6095)/2.;
    X1=(4.4855+3.8743)/2.;
    X0=(1.7635+1.5107)/2.;
    Z=(18.+18/8.)/2.;
  }
  else if(!matname.CompareTo("CFRP") ){
    rho=1.700;
    I=78;
    Z_A=0.49954;  
    C0=-3.1550;
    a=0.20762;
    M=2.9532;
    X1=2.5387;
    X0=0.0480;
    Z=6.;
  }
  else if(!matname.CompareTo("Scinti") || !matname.CompareTo("Plastic")){
    rho=1.032;
    I=64.7;
    Z_A=0.54141;  
    C0=-3.1997;
    a=0.16101;
    M=3.2393;
    X1=2.4855;
    X0=0.1464;
    Z=3.5;
  } 
  else if(!matname.CompareTo("Aerogel") ){
    //SiO2, take medium of Si & O
    rho=0.2;
    I=(173+95)/2.;
    Z_A=0.5;  
    C0=-(10.7+4.435)/2.;
    a=(0.15+0.118)/2.;
    M=(3.255+3.29)/2.;
    X1=(2.8715+4.3213)/2.;
    X0=(0.2014+1.7541)/2.;
    Z=11;
  } 
  else if(!matname.CompareTo("Mylar") ){
    rho=1.39;
    I=78.7;
    Z_A=0.52037;
    C0=-3.3262;
    a=0.12679;
    M=3.3076;
    X1=2.6507;
    X0=0.1562;
    Z=50./96.*11.;
  }
  else if(!matname.CompareTo("Beryllium") ){
    rho=1.848;
    I=63.7;
    Z_A=0.44384;  
    C0=-2.7847;
    a=0.80392;
    M=2.4339;
    X1=1.6922;
    X0=0.0592;
    Z=4.;
  }
  else if(!matname.CompareTo("AlBeMet") ){
    //Al 38% Be 62%
    rho=2.071;
    I=63.7*0.62+166.*0.38;
    Z_A=0.44384*0.62+0.48181*0.38;  
    C0=-2.7847*0.62-4.2395*0.38;
    a=0.80392*0.62+0.08024*0.38;
    M=2.4339*0.62+3.6345*0.38;
    X1=1.6922*0.62+3.0127*0.38;
    X0=0.0592*0.62+0.1708*0.38;
    Z=4.*0.62+13*0.38;
  }
  else if(!matname.CompareTo("Iron") ){
    rho=7.874;
    I=286.;
    Z_A=0.46556;  
    C0=-4.2911;
    a=0.14680;
    M=2.9632;
    X1=3.1531;
    X0=-0.0012;
    Z=26.;
  }
  else if(!matname.CompareTo("LHelium-3") ){
    rho=0.08;
    I=41.8;
    Z_A=2./3.;  
    C0=-11.1393;
    a=0.13443;
    M=5.8347;
    X1=3.6122;
    X0=2.2017;
    Z=2.;
  }
  else if(!matname.CompareTo("LHydrogen") ){
    rho=0.071;
    I=21.8;
    Z_A=1./1.;  
    C0=-3.2632;
    a=0.13483;
    M=5.6249;
    X1=1.9215;
    X0=0.4759;
    Z=1.;
  }
  else if(!matname.CompareTo("LDeuterium") ){
    rho=0.168;
    I=21.8;
    Z_A=1./2.;  
    C0=-3.2632;
    a=0.13483;
    M=5.6249;
    X1=1.9215;
    X0=0.4759;
    Z=1.;
  }
  else if(!matname.CompareTo("Air") || !matname.CompareTo("BLDCGas") ){
    rho=1.2048/1000.;
    I=85.7;
    Z_A=0.49919;  
    C0=-10.5961;
    a=0.10914;
    M=3.3994;
    X1=4.2759;
    X0=1.7418;
    Z=6.8;
  }
  else if(!matname.CompareTo("Vacuum") ){
    rho=1.2048/100000.;
    I=85.7;
    Z_A=0.49919;  
    C0=-10.5961;
    a=0.10914;
    M=3.3994;
    X1=4.2759;
    X0=1.7418;
    Z=6.8;
  }
  else{
    std::cout<<"Material name "<<matname<<" is not defined"<<std::endl;
    return 0.;
  }
  //  Double=Mathtools(Beta,I,Z_A);
  double val2=ELossTools::CalcdEdX(beta,rho,I,Z_A,Z,C0,a,M,X0,X1,CORR);
  //std::cout<<matname<<"  beta "<<beta<<"  rho "<<rho<<"  I "<<I<<" eloss = "<<val1<<"  "<<val2<<std::endl;
  return val2;
}

double ELossTools::CalcdEdX(const double &beta,const double &rho, const double &I, const double &Z_A, const double &Z,
			    const double &C0, const double &a, const double &M, const double &X0, const double &X1, bool CORR){
  const double C=0.1535; /*MeVcm^2/g*/
  const double m_e=0.511; // MeV/c^2

  double gamma_2=1/(1-pow(beta,2.0));
  double gamma=sqrt(gamma_2);
  double W_max=2.0*m_e*pow(beta,2.0)*gamma_2;
  //  double C0=-3.2632;
  //  double X0=0.4759;
  //  double X1=1.9215;
  //  double a=0.13483;
  //  double M=5.6249;
  
  double X=0.0;
  double C_shell=0.0;
  double delta=0.0;
  double corr=0.0;
  if(CORR){
    double eta=beta*gamma;
    C_shell=(0.422377/pow(eta,2)+0.0304043*pow(eta,-4)-0.00038106*pow(eta,-6))*1e-6*pow(I,2)
      +(3.850190*pow(eta,-2)-0.1667989*pow(eta,-4)+0.00157955*pow(eta,-6))*1e-9*pow(I,3);
    X = log10(eta);
    if (X<=X0)
      delta=0.0;
    else if (X0<X && X<X1) 
      delta=4.6052*X+C0+a*pow((X1-X),M);
    else if (X>=X1)  
      delta=4.6052*X+C0;
    corr=-delta-2*C_shell/Z;
  }
  
  //  std::cout<<"delta,C_shell   "<<delta<<"  "<<C_shell<<std::endl;
  int z=1;
  double logterm=log(2.0*m_e*gamma_2*pow(beta,2.0)*W_max*pow(10.0,12.0)/pow(I,2.0))-2.0*pow(beta,2.0)+corr;
  double val=C*rho*Z_A*pow((double)z,2.0)*logterm/pow(beta,2.0);
  
  return val/1000.; // MeV->GeV
}
