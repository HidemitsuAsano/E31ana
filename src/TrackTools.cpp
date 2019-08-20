// TrackTools.cpp
#include "TrackTools.h"
#include "TF1.h"

bool TrackTools::CalcBetaMass2TOF(CDSTrack *track, LocalTrack* bpc,double tof,double mass,double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof)
{
  double param[5];
  TVector3 tmp;
  track->GetParameters(CID_CDC,param,tmp);   
  TVector3 cdhvtx=track->CDHVertex();
  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.1);
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,pos1);
  TVector3 vtxb1,vtxb2;
  double dis;
  TrackTools::CalcLineHelixVertex(bpc,track,vtxb1,vtxb2,dis);
  double mom,time_beam;
  TVector3 t0pos=bpc->GetPosatZ(-110-0.5);
  ELossTools::CalcElossBeamTGeo(t0pos,vtxb1,beammom,parMass[pid_beam],mom,time_beam);
  mom = track->Momentum();   
  double tmpl;
  ELossTools::CalcHelixElossToVertexTGeo(param,vtxb2,mom,mass,tmpl,tmptof,15.1);    
  tmptof+=time_beam;
  beta_calc=(cdc_dis)/(tof-tmptof)/(Const*100);
  mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
  return true;
}

bool TrackTools::FindMass2(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof)
{
  double param[5];
  TVector3 tmp;
  track->GetParameters(CID_CDC,param,tmp);   
  TVector3 cdhvtx=track->CDHVertex();
  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.1);
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,pos1);
  TVector3 vtxb1,vtxb2;
  double dis;
  TrackTools::CalcLineHelixVertex(bpc,track,vtxb1,vtxb2,dis);
  double mom,time_beam;
  TVector3 t0pos=bpc->GetPosatZ(-110-0.5); 

  ELossTools::CalcElossBeamTGeo(t0pos,vtxb1,beammom,parMass[pid_beam],mom,time_beam);
  mom = track->Momentum();
  int count =0;
  double tmpl,tmpmass2;
  do{
    tmpmass2=mass2;
    double tmpmass=sqrt(tmpmass2);
    if(tmpmass2<0.0005) tmpmass=0.0005;
    ELossTools::CalcHelixElossToVertexTGeo(param,vtxb2,mom,tmpmass,tmpl,tmptof,15.1);    
    //    std::cout<<tmptof<<"  "<<tmpl<<std::endl;
    tmptof+=time_beam;
    beta_calc=(cdc_dis)/(tof-tmptof)/(Const*100);
    mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
    //    std::cout<<"tmptof,beta,mass2  "<<tmptof<<"  "<<beta_calc<<"  "<<mass2<<std::endl;
    count++;
    if(count>50)   return false;
  } while( TMath::Abs(tmpmass2-mass2)>0.0001);
  return true;
}

bool TrackTools::FindMass2C(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof)
{
  double param[5];
  TVector3 tmp;
  track->GetParameters(CID_CDC,param,tmp);   
  TVector3 cdhvtx=track->CDHVertex();
  //  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.1);
  TVector3 pos1=tmp;
  double cdc_dis=MathTools::CalcHelixArc(param,pos1,cdhvtx);
  //  TVector3 vtxb1,vtxb2;
  //  double dis;
  //  TrackTools::CalcLineHelixVertex(bpc,track,vtxb1,vtxb2,dis);
  double mom,tmpmom,time_beam;
  TVector3 t0pos=bpc->GetPosatZ(-110-0.5); 
  //TVector3 t0pos=bpc->GetPosatZ(-110); 

  int count =0;
  double tmpl,tmpmass2;
  TVector3 lpos,hpos;
  mom = track->Momentum();   
  //  std::cout<<"---------"<<std::endl;
  do{
    if(mass2<0) mass2=0;
    tmpmass2=mass2;
    double tmpmass=sqrt(tmpmass2);
    if(tmpmass2<0.0005) tmpmass=0.0005;
    TVector3 tmpff=bpc->GetPosatZ(0);

    if(!track->CalcVertexTimeLength(tmpff,bpc->GetMomDir(),tmpmass,lpos,hpos,tmptof,tmpl)) {  return false;}
    //    std::cout<<"C  "<<tmptof<<"  "<<tmpl<<std::endl;
    ELossTools::CalcElossBeamTGeo(t0pos,lpos,beammom,parMass[pid_beam],tmpmom,time_beam);
    tmptof+=time_beam;
    beta_calc=(cdc_dis)/(tof-tmptof)/(Const*100);
    mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
    //    std::cout<<"tmptof,beta,mass2  "<<tmptof<<"  "<<beta_calc<<"  "<<mass2<<std::endl;
    mass2=(mass2+tmpmass2)/2.;
    //    std::cout<<count<<" tmpmass2 "<<tmpmass2<<" mass2 "<<mass2<<std::endl;

    count++;
    if(mass2<=0) return true;
    if(count>80)   {std::cout<<"f2c difmass "<<TMath::Abs(tmpmass2-mass2)<<std::endl;
      return false;}
  } while( TMath::Abs(tmpmass2-mass2)>0.0005);
  //  std::cout<<"mass "<<sqrt(mass2)<<" beam1 "<<time_beam<<std::endl;
  return true;
}


bool TrackTools::Calc2HelixVertex(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2)
{
  double par1[5];
  double par2[5];
  TVector3 tmppos1,tmppos2;
  TVector3 tmp1,tmp2;
  int id,id2;
  double margin=0.5;
  double dis, tmpdis=999;
  //  std::cout<<"---------------------------"<<std::endl;
  // std::cout<<cds1->PID()<<"  "<<cds1->CDHFlag()<<"  "<<cds1->Mass()<<std::endl;
  // cds1->Print();
  // std::cout<<cds2->PID()<<"  "<<cds2->CDHFlag()<<"  "<<cds2->Mass()<<std::endl;
  // cds2->Print();
  for(int i=0;i<cds1->nParamSets();i++){
    cds1->GetNthParameters(i,id,par1,tmppos1);
    for(int i2=0;i2<cds2->nParamSets();i2++){
      cds2->GetNthParameters(i2,id2,par2,tmppos2);
      MathTools::HelixToHelix(par1,par2,tmp1,tmp2,dis);
      //      if((GeomTools::GetID(tmp1)==id||GeomTools::GetID(tmp2)==id2)&&dis<tmpdis){
      if((GeomTools::GetID(tmp1)==id||GeomTools::IsSameVolumeHelix(par1,tmppos1,tmp1,margin))){	
	// std::cout<<dis<<"  "<<tmpdis<<std::endl;
	// std::cout<<GeomTools::GetID(tmp1)<<"  "<<id<<"   "<<tmp1.Perp()<<"  "<<tmp1.Z()<<std::endl;
	// std::cout<<GeomTools::GetID(tmp2)<<"  "<<id2<<"   "<<tmp2.Perp()<<"  "<<tmp2.Z()<<std::endl;
	if((GeomTools::GetID(tmp2)==id2||GeomTools::IsSameVolumeHelix(par2,tmppos2,tmp2,margin))
	   &&dis<tmpdis){
	  vtx1=tmp1;
	  vtx2=tmp2;
	  tmpdis=dis;
	  //	  std::cout<<"!!!"<<std::endl;
	}
      }
    }
  }
  //  std::cout<<"-----------!!!"<<std::endl;
  if(tmpdis<100) return true;
  else{
    //    std::cout<<tmpdis<<std::endl;
    return false;
  }
}


bool TrackTools::Calc2HelixVertexWresl(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2)
{
  double par1[5];
  double par2[5];
  TVector3 tmppos1,tmppos2;
  TVector3 tmp1,tmp2;
  int id,id2;
  double margin=0.5;
  double dis, tmpdis=999;
  for(int i=0;i<cds1->nParamSets();i++){
    cds1->GetNthParameters(i,id,par1,tmppos1);
    for(int i2=0;i2<cds2->nParamSets();i2++){
      cds2->GetNthParameters(i2,id2,par2,tmppos2);
      MathTools::HelixToHelixWresl(par1,par2,tmp1,tmp2,dis);
      if((GeomTools::GetID(tmp1)==id||GeomTools::IsSameVolumeHelix(par1,tmppos1,tmp1,margin))){	
	if((GeomTools::GetID(tmp2)==id2||GeomTools::IsSameVolumeHelix(par2,tmppos2,tmp2,margin))
	   &&dis<tmpdis){
	  vtx1=tmp1;
	  vtx2=tmp2;
	  tmpdis=dis;
	  //	  std::cout<<"!!!"<<std::endl;
	}
      }
    }
  }
  //  std::cout<<"-----------!!!"<<std::endl;
  if(tmpdis<100) return true;
  else{
    //    std::cout<<tmpdis<<std::endl;
    return false;
  }
}


bool TrackTools::Calc2HelixVertex2(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2)
{
  double par1[5];
  double par2[5];
  TVector3 tmpvtx;
  double dis;
  cds1->GetParameters(CID_CDC,par1,tmpvtx);
  cds2->GetParameters(CID_CDC,par2,tmpvtx);
  MathTools::HelixToHelix(par1,par2,vtx1,vtx2,dis);
  if(dis<10) return true;
  else return false;
}

bool TrackTools::CalcLineHelixVertex(LocalTrack *track,CDSTrack *cds, TVector3 &vtx1, TVector3 &vtx2, double &dis)
{
  double par[5];
  TVector3 tmpvtx;
  cds->GetParameters(CID_CDC,par,tmpvtx);
  TVector3 lpos=track->GetPosatZ(0);
  TVector3 ldir=track->GetMomDir();
  if(!MathTools::LineToHelix(lpos,ldir,par,vtx1,vtx2,dis)) return false;
  if(dis<10) return true;
  else return false;
}

bool TrackTools::CalcLineHelixVertex(pBeam *beam,CDSTrack *cds,TVector3 &vtx1,TVector3 &vtx2,double &dis, bool ELOSS)
{
  if(ELOSS){
    cds->GetVertex(beam->bpcpos(),beam->bpcdir(),vtx2,vtx1);
    dis=(vtx1-vtx2).Mag();
  }else{
    double par[5];
    TVector3 tmpvtx;
    cds->GetParameters(CID_CDC,par,tmpvtx);
    TVector3 lpos=beam->bpcpos();
    TVector3 ldir=beam->bpcdir();
    if(!MathTools::LineToHelix(lpos,ldir,par,vtx1,vtx2,dis)) return false;
  }
  if(dis<10) return true;
  else return false;
}

void TrackTools::CalcBetaMass(TVector3 vertex,LocalTrack *beam,CDSTrack* cdc, 
			      ConfMan *conf,int beam_pid,double tof, double &beta,double &mass2)
{
  double param[5];
  cdc->GetParameters(param);
  if(param[2]==0 ||param[4]==0   ) return;
  double mom = cdc->Momentum();

  TVector3 gpos;
  conf->GetGeomMapManager()->GetPos(CID_T0,0,gpos);
  double z_t0=gpos.Z(),z_vtxb=vertex.Z();
  TVector3 t0pos=beam->GetPosatZ(z_t0);
  TVector3 vpos=beam->GetPosatZ(z_vtxb);
  double beam_dis=(t0pos-vpos).Mag();

  double beta_b=0.;
  if(beam_pid==Beam_Kaon)  beta_b=1.0/sqrt(1.0*1.0+kpMass*kpMass);
  else if(beam_pid==Beam_Pion)  beta_b=1.0/sqrt(1.0*1.0+piMass*piMass);  
  double time_beam=beam_dis/beta_b/(Const*100);
  
  //#####CDC dis#######		  	
  TVector3 cdhvtx=cdc->CDHVertex();
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,vertex);

  double beta_calc=cdc_dis/(tof-time_beam)/(Const*100);
  double calc_mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
  
  beta = beta_calc;
  mass2 = calc_mass2;
#if 0
  std::cout<<"vertex.x,y,z"<<vertex.X()<<"\t"<<vertex.Y()<<"\t"<<vertex.Z()<<std::endl;
  std::cout<<"cdh x,y,z"<<cdhvtx.X()<<"\t"<<cdhvtx.Y()<<"\t"<<cdhvtx.Z()<<std::endl;
  std::cout<<"tof, beta, mass"<<tof<<"\t"<<beta<<"\t"<<mass<<std::endl;
#endif
}

void TrackTools::CalcBetaMass2(TVector3 vertex,LocalTrack *beam,CDSTrack* cdc, 
			       ConfMan *conf,int beam_pid,double beammom,double tof, double &beta,double &mass2)
{
  double param[5];
  cdc->GetParameters(param);
  if(param[2]==0 ||param[4]==0   ) return;
  double mom = cdc->Momentum();

  TVector3 gpos;
  conf->GetGeomMapManager()->GetPos(CID_T0,0,gpos);
  double z_t0=gpos.Z(),z_vtxb=vertex.Z();
  TVector3 t0pos=beam->GetPosatZ(z_t0);
  TVector3 vpos=beam->GetPosatZ(z_vtxb);
  double beam_dis=(t0pos-vpos).Mag();

  double beta_b=0.;
  if(beam_pid==Beam_Kaon)  beta_b=beammom/sqrt(beammom*beammom+kpMass*kpMass);
  else if(beam_pid==Beam_Pion)  beta_b=beammom/sqrt(beammom*beammom+piMass*piMass);  
  double time_beam=beam_dis/beta_b/(Const*100);
  
  //#####CDC dis#######		  	
  TVector3 cdhvtx=cdc->CDHVertex();
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,vertex);

  double beta_calc=cdc_dis/(tof-time_beam)/(Const*100);
  double calc_mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
  
  beta = beta_calc;
  mass2 = calc_mass2;
#if 0
  std::cout<<"vertex.x,y,z"<<vertex.X()<<"\t"<<vertex.Y()<<"\t"<<vertex.Z()<<std::endl;
  std::cout<<"cdh x,y,z"<<cdhvtx.X()<<"\t"<<cdhvtx.Y()<<"\t"<<cdhvtx.Z()<<std::endl;
  std::cout<<"tof, beta, mass"<<tof<<"\t"<<beta<<"\t"<<mass<<std::endl;
#endif
}


int TrackTools::PID(double mom,double mass2)
{
  int ptype=-1;
  if(mom>0) {
    if(mass2<0.3) ptype=CDS_PiPlus;
    else if(mass2<2.) ptype=CDS_Proton;
    else if(mass2<5.) ptype=CDS_Deuteron;
    else if(mass2<9.) ptype=CDS_Triton;
    else ptype=CDS_Other;
  }
  else{
    if(mass2<0.12) ptype=CDS_PiMinus;
    else if(mass2<0.5) ptype=CDS_Kaon;
    else ptype=CDS_Other;
  }
  return ptype;
}


int TrackTools::PIDcorr(double mom,double mass2)
{
  int ptype=CDS_Other;
  double fmom=fabs(mom);
  TF1 *f_sigma[4];

  double param[4][4]={{0.049154,0.011690,0.090581,piMass},
		      {0.049154,0.011690,0.090581,kpMass},
		      {0.049154,0.011690,0.090581,pMass},
		      {0.049154,0.011690,0.090581,dMass}};
  for(int n=0;n<4;n++)
    {
      f_sigma[n]=new TF1("f_sigma","sqrt( 4*pow([3],4)*[0]*[0]*x*x+4*pow([3],4)*(1+pow([3]/x,2))*[1]*[1]+4*x*x*([3]*[3]+x*x)*[2]*[2])");
      f_sigma[n]->SetParameters(param[n][0],param[n][1],param[n][2],param[n][3],param[n][4]);
    }

  if(mom>0) {
    if(mom<0.1)
      {
	double th1=piMass*piMass-2.5*f_sigma[0]->Eval(fmom);
	double th2=piMass*piMass+2.5*f_sigma[0]->Eval(fmom);
	if(th1<mass2 && mass2<th2) ptype=CDS_PiPlus;
      }
    else if(mom>=0.1 && mom<0.19)
      {
	double th1_pi=piMass*piMass-2.5*f_sigma[0]->Eval(fmom);
	double th2_pi=piMass*piMass+2.5*f_sigma[0]->Eval(fmom);
	double th1_p=pMass*pMass-2.5*f_sigma[2]->Eval(fmom);
	double th2_p=pMass*pMass+2.5*f_sigma[2]->Eval(fmom);

	if(th1_pi<mass2 && mass2<th2_pi) ptype=CDS_PiPlus;
	else if(th1_p<mass2 && mass2<th2_p) ptype=CDS_Proton;
      }
    else if(mom>=0.19)
      {
	double th1_pi=piMass*piMass-2.5*f_sigma[0]->Eval(fmom);
	double th2_pi=piMass*piMass+2.5*f_sigma[0]->Eval(fmom);
	double th1_p=pMass*pMass-2.5*f_sigma[2]->Eval(fmom);
	double th2_p=pMass*pMass+2.5*f_sigma[2]->Eval(fmom);
	double th1_d=dMass*dMass-2.5*f_sigma[3]->Eval(fmom);
	double th2_d=dMass*dMass+2.5*f_sigma[3]->Eval(fmom);

	if(th1_pi<mass2 && mass2<th2_pi && mass2<th1_p) ptype=CDS_PiPlus;
	else if(th1_p<mass2 && mass2<th2_p && mass2>th2_pi) ptype=CDS_Proton;
	else if(th1_d<mass2 && mass2<th2_d) ptype=CDS_Deuteron;

      }
  }
  else if(mom<0) {    
    double th1_pi=piMass*piMass-2.5*f_sigma[0]->Eval(fmom);
    double th2_pi=piMass*piMass+2.5*f_sigma[0]->Eval(fmom);
    double th1_k=kpMass*kpMass-2.5*f_sigma[1]->Eval(fmom);
    double th2_k=kpMass*kpMass+2.5*f_sigma[1]->Eval(fmom);
    if(th1_pi<mass2 && mass2<th2_pi && (mom>-0.05 || mass2<th1_k) ) ptype=CDS_PiMinus;    
    else if(mom<-0.05 && th1_k<mass2 && mass2<th2_k && mass2>th2_pi) ptype=CDS_Kaon;    
  }

  delete f_sigma[0];
  delete f_sigma[1];
  delete f_sigma[2];
  delete f_sigma[3];
  return ptype;
}


int TrackTools::PIDcorr_wide(double mom,double mass2)
{
  const double PID_Param[3] = { 0.00381414, 0.000119896, 0.0113647 };
  const double Kpi_mid_mass2 = 0.1031;
  const double Ppi_mid_mass2 = 0.76247;

  double pi_sigma = sqrt(4.*piMass*piMass*mom*mom*PID_Param[0]+4.*piMass*piMass*piMass*piMass*PID_Param[1]*(1+piMass*piMass/(mom*mom))+4.*mom*mom*(piMass*piMass+mom*mom)*PID_Param[2]);
  double k_sigma = sqrt(4.*kpMass*kpMass*mom*mom*PID_Param[0]+4.*kpMass*kpMass*kpMass*kpMass*PID_Param[1]*(1+kpMass*kpMass/(mom*mom))+4.*mom*mom*(kpMass*kpMass+mom*mom)*PID_Param[2]);
  double p_sigma = sqrt(4.*pMass*pMass*mom*mom*PID_Param[0]+4.*pMass*pMass*pMass*pMass*PID_Param[1]*(1+pMass*pMass/(mom*mom))+4.*mom*mom*(pMass*pMass+mom*mom)*PID_Param[2]);

  bool pim_flag=false, pip_flag=false, km_flag=false, p_flag=false;
  if( mom>0 ){
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pip_flag=true;
    if( pMass*pMass-2.5*p_sigma<mass2 && mass2<pMass*pMass+2.5*p_sigma && mom>0.1 ) p_flag=true;

    if( pip_flag && p_flag ){
      if( mass2<Ppi_mid_mass2 ) return CDS_PiPlus;
      else return CDS_Proton;
    }
    else if( pip_flag ) return CDS_PiPlus;
    else if( p_flag ) return CDS_Proton;
  }
  else{
    if( piMass*piMass-2.5*pi_sigma<mass2 && mass2<piMass*piMass+2.5*pi_sigma ) pim_flag=true;
    if( kpMass*kpMass-2.5*k_sigma<mass2 && mass2<kpMass*kpMass+2.5*k_sigma && mom<-0.03 ) km_flag=true;

    if( pim_flag ) return CDS_PiMinus;
    else if( km_flag ) return CDS_Kaon;
  }
  return CDS_Other;
}


int TrackTools::PIDcorr3(double mom,double mass2)
{
  int ptype=CDS_Other;
  double fmom=fabs(mom);
  //double param[4][4]={{0.00227537,0.000112426,piMass*piMass,0.0130341},
  //                    {0.00227537,0.000112426,kpMass*kpMass,0.0130341},
  //                    {0.00227537,0.000112426,pMass*pMass  ,0.0130341},
  //                    {0.00227537,0.000112426,dMass*dMass  ,0.0130341}};
  //double pi_mass2 = 0.139570*0.139570;
  //double k_mass2 = 0.493600*0.493600;
  //double p_mass2 = 0.938272*0.938272;
  //double d_mass2 = 1.875610*1.875610;
  
  /* 2016.07.05 -----> */
  double pi_mass2 = 0.139570*0.139570;
  double k_mass2 = 0.497234*0.497234;
  double p_mass2 = 0.918312*0.918312;
  double d_mass2 = 1.829990*1.829990;
  /* <----- 2016.07.05 */

  //double param[4][4]={{0.00227537,0.000112426,pi_mass2,0.0130341},
  //                    {0.00227537,0.000112426,k_mass2,0.0130341},
  //                    {0.00227537,0.000112426,p_mass2,0.0130341},
  //                    {0.00227537,0.000112426,d_mass2,0.0130341}};
  //double param[4][4]={{0.00049612,0.000105752,pi_mass2,0.00809271},
  //                    {0.00049612,0.000105752,k_mass2 ,0.00809271},
  //                    {0.00049612,0.000105752,p_mass2 ,0.00809271},
  //                    {0.00049612,0.000105752,d_mass2 ,0.00809271}};
  
  /* 2016.07.05 -----> */
  double param[4][4]={{0.00084734,0.00011146,pi_mass2,0.00674811},
                      {0.00084734,0.00011146,k_mass2 ,0.00674811},
                      {0.00084734,0.00011146,p_mass2 ,0.00674811},
                      {0.00084734,0.00011146,d_mass2 ,0.00674811}};
  /* <----- 2016.07.05 */

  double pi_sigma = sqrt(4.0*param[0][2]*param[0][0]*mom*mom+4.0*param[0][2]*param[0][2]*param[0][1]*(1.0+param[0][2]/(mom*mom))+4.0*param[0][3]*mom*mom*(param[0][2]+mom*mom));
  double k_sigma = sqrt(4.0*param[1][2]*param[1][0]*mom*mom+4.0*param[1][2]*param[1][2]*param[1][1]*(1.0+param[1][2]/(mom*mom))+4.0*param[1][3]*mom*mom*(param[1][2]+mom*mom));
  double p_sigma = sqrt(4.0*param[2][2]*param[2][0]*mom*mom+4.0*param[2][2]*param[2][2]*param[2][1]*(1.0+param[2][2]/(mom*mom))+4.0*param[2][3]*mom*mom*(param[2][2]+mom*mom));
  double d_sigma = sqrt(4.0*param[3][2]*param[3][0]*mom*mom+4.0*param[3][2]*param[3][2]*param[3][1]*(1.0+param[3][2]/(mom*mom))+4.0*param[3][3]*mom*mom*(param[3][2]+mom*mom));
  /*=== positive ===*/
  double n_sigma = 2.5;
  double pi_ll = pi_mass2 - n_sigma*pi_sigma;
  double pi_ul = pi_mass2 + n_sigma*pi_sigma;
  double k_ll = k_mass2 - n_sigma*k_sigma;
  double k_ul = k_mass2 + n_sigma*k_sigma;
  double p_ll = p_mass2 - n_sigma*p_sigma;
  double p_ul = p_mass2 + n_sigma*p_sigma;
  double d_ll = d_mass2 - n_sigma*d_sigma;
  double d_ul = d_mass2 + n_sigma*d_sigma;
  if(mom>0){
    bool pi_flag = false;
    bool p_flag = false;
    bool d_flag = false;
    if(pi_ll<mass2&&mass2<pi_ul){
      pi_flag = true;
    }
    if(p_ll<mass2&&mass2<p_ul){
      if(mom>0.1){
        p_flag = true;
      }
    }
    if(d_ll<mass2&&mass2<d_ul){
      if(mom>0.2){
        d_flag = true;
      }
    }
    /*=== pi+ ===*/
    if(pi_flag){
      if(!p_flag){
        ptype = CDS_PiPlus; return ptype;
      }
    }
    /*=== proton ===*/
    if(p_flag){
      if(!pi_flag&&!d_flag){
        ptype = CDS_Proton; return ptype;
      }
    }
    /*=== deuteron ===*/
    if(d_flag){
      if(!p_flag){
        ptype = CDS_Deuteron; return ptype;
      }
    }
  }
  /*=== negative ===*/
  else {
    bool pi_flag = false;
    bool k_flag = false;
    if(pi_ll<mass2&&mass2<pi_ul){
      pi_flag = true;
    }
    if(k_ll<mass2&&mass2<k_ul){
      if(mom<-0.05){
        k_flag = true;
      }
    }
    /*=== pi- ===*/
    if(pi_flag){
      if(!k_flag){
        ptype = CDS_PiMinus; return ptype;
      }
    }
    /*=== k- ===*/
    if(k_flag){
      if(!pi_flag){
        ptype = CDS_Kaon; return ptype;
      }
    }
  }

  ptype = CDS_Other; return ptype;
}


















































pCDS* TrackTools::CalcSingleAll(pBeam *beam,CDSTrack* cdc,CDSHitMan *cdsMan,bool ELOSS)
{
  double cdhtime,dis;
  int cdhseg;
  if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) return 0;
  TVector3 vtx1,vtx2;
  double par[5];
  cdc->GetParameters(CID_CDC,par,vtx1);
  double mom=cdc->Momentum();
  TrackTools::CalcLineHelixVertex(beam,cdc,vtx1,vtx2,dis);

  pCDS* cdstrack=new pCDS();
  cdstrack->SetParameters(par);
  cdstrack->SetTOF(cdhtime);
  cdstrack->SetCDHSeg(cdhseg);

  double time_beam=(beam->t0pos()-cdstrack->vbeam()).Mag()/beam->beta()/(Const*100);
  TVector3 cdhvtx=cdc->CDHVertex();
  double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,vtx2);
  double beta=cdc_dis/(cdhtime-beam->t0time()-time_beam)/(Const*100);
  double mass2=mom*mom*(1/(beta*beta)-1);
  int pid=TrackTools::PID(mom,mass2);
  // for(int i=0;i<5;i++)  std::cout<<par[i]<<"  "; std::cout<<std::endl;
  // std::cout<<beam->mom()<<"  "<<beam->beta()<<"  "<<time_beam<<"  "<<cdc_dis<<"  "<<mom<<"  "<<beta<<"  "<<mass2<<std::endl;
  cdc->SetPID(pid);
  cdc->CalcELoss();
  if(ELOSS)   TrackTools::CalcLineHelixVertex(beam,cdc,vtx1,vtx2,dis,ELOSS);

  TVector3 Pd;
  if( !cdc->GetMomentum(vtx1, Pd ,ELOSS, ELOSS) ){
    return 0;
    // std::cout<<"error of GetMomentum() track"<<std::endl;;
    // std::cout<<"------------------------------------"<<std::endl;;
  }      

  //  cdstrack->SetTrackID(trackID);
  cdstrack->SetVertexCDC(vtx1);
  cdstrack->SetVertexBeam(vtx2);
  cdstrack->SetVDis(dis);
  cdstrack->SetMomDir(Pd);
  cdstrack->SetPID(pid);
  cdstrack->SetMomentum(mom/TMath::Abs(mom)*Pd.Mag());
  cdstrack->SetBeta(beta);
  cdstrack->SetMass(sqrt(mass2));
  cdstrack->SetMass2(mass2);
  cdstrack->SetAngleLab(beam->bpcdir().Angle(Pd));    
  
  for(int icdh=0;icdh<cdc->nCDHHit();icdh++)
    cdstrack->SetCDHSeg(cdc->CDHHit(cdsMan,icdh)->seg());
  for(int icdh=0;icdh<cdc->nIHHit();icdh++)
    cdstrack->SetIHSeg(cdc->IHHit(cdsMan,icdh)->seg());
  
  return cdstrack;
}

pCDS *TrackTools::Calc2HelixAll(pBeam *beam, pCDS* cds1, pCDS* cds2,CDSTrackingMan *trackMan, bool ELOSS)
{
  CDSTrack *track1=trackMan->Track(cds1->id());	
  CDSTrack *track2=trackMan->Track(cds2->id());	
  //double* paramtmp=cds1->GetParameters();
  TrackVertex vertex;
  TVector3 vtx=DEFVECT,vtx1=DEFVECT,vtx2=DEFVECT;
  
  if(ELOSS){
    if( TrackTools::Calc2HelixVertex(track1,track2,vtx1,vtx2) ){
      vtx=(vtx1+vtx2)*0.5;
    }
  }else if( trackMan->GetVertex(cds1->id(),cds2->id(),vertex) )
    {
      vtx=vertex.GetVertexPos_mean();
      vtx1=vertex.GetVertexPos1();
      vtx2=vertex.GetVertexPos2();
    }
  else return 0;
  TVector3 Pp1,Pp2;
  if( !track1->GetMomentum(vtx1,Pp1,ELOSS,ELOSS) ){
    //    std::cout<<"Track1 !!! "<<std::endl;
    return 0;
  }
  if( !track2->GetMomentum(vtx2,Pp2,ELOSS,ELOSS) ){
    //    std::cout<<"Track2 !!! "<<std::endl;
    return 0;
  }
  TVector3 Pp= Pp1+Pp2;
  TLorentzVector L1; L1.SetVectM( Pp1, cds1->pdgmass() );
  TLorentzVector L2; L2.SetVectM( Pp2, cds2->pdgmass() );
  double im = (L1+L2).M();	
  
  double dist,dltmp=0;
  TVector3 xest,nest;

  MathTools::LineToLine( vtx,Pp.Unit(),beam->bpcpos(), beam->bpcdir(),dltmp,dist,xest,nest );
  pCDS* pro=new pCDS();
    
  pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
  pro->SetMomentum(Pp.Mag());
  pro->SetMass(im);
  pro->SetVertex(vtx);
  pro->SetVertexBeam(xest);
  pro->SetMomDir(Pp.Unit());
  pro->SetVDis((vtx1-vtx2).Mag());
  pro->SetVBDis((xest-vtx).Mag());
  pro->SetPBDCA(dist);
  pro->SetOA(Pp1.Angle(Pp2));
  pro->SetAngleLab(beam->bpcdir().Angle(Pp));
  return pro;
}
