// TrackTools.cpp
#include "TrackTools.h"
#include "TF1.h"
double TrackTools::Mass2Resol(int pid, double mom, double fl)
{
  double tofresol=0.160;
  double lv=29.97;//cm/ns  
  double mass=cdsMass[pid];
  double momres=MomResol(pid,mom);
  double tmp2=4*pow(mass,4)*pow(momres,2)+4*pow(mom,2)*(pow(mom,2)+pow(mass,2))*pow(lv/fl*tofresol,2);
  if(tmp2>0)
    return sqrt(tmp2);
  return 999.;
}

double TrackTools::MomResol(int pid, double mom)
{
  int id=-1;
  if(pid==CDS_PiPlus||pid==CDS_PiMinus) id=0;
  else if(pid==CDS_Proton) id=1;
  else if(pid==CDS_Kaon)  id=2;
  else return 999;
  // for mom resol. [0]+[1]*x+[2]/(x-[3])
#if 0
  // parameter for de corrected
  double param[3][4]={{1.03626e-05, 0.062906, 0.000538243, 0.0421856},
		      {-5.6182e-05, 0.0632374, 0.0002, 0.195},
		      {0.00216906, 0.0596655, 0.000187301, 0.13}
  };
#else  
  //parameter at CDC 
  double param[3][4]={{-0.000154926, 0.0628548, 0.000544216, 0.0203668}, 
		      {-0.0119525, 0.0699056, 0.0058263, 0.00385919}, 
		      {-0.00514371, 0.0653424, 0.00293869, -0.0220217}
  };
#endif
  double tmp=param[id][0]+param[id][1]*mom+param[id][2]/(mom-param[id][3]);
  if(tmp>0)
    return tmp;
  return 999;
}

bool TrackTools::CalcBetaMass2TOF(CDSTrack *track, LocalTrack* bpc,double tof,double mass,double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof)
{
  double param[5];
  TVector3 tmp;
  track->GetParameters(CID_CDC,param,tmp);   
  TVector3 cdhvtx=track->CDHVertex();
  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.4);
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,pos1);
  TVector3 vtxb1,vtxb2;
  double dis;
  TrackTools::CalcLineHelixVertex(bpc,track,vtxb1,vtxb2,dis);
  double mom,time_beam;
  TVector3 t0pos=bpc->GetPosatZ(-110-0.5);
  ELossTools::CalcElossBeamTGeo(t0pos,vtxb1,beammom,parMass[pid_beam],mom,time_beam);
  mom = track->Momentum();   
  double tmpl;
  ELossTools::CalcHelixElossToVertexTGeo(param,vtxb2,mom,mass,tmpl,tmptof,15.4);    
  tmptof+=time_beam;
  beta_calc=(cdc_dis)/(tof-tmptof)/(Const*100);
  mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
  return true;
}

bool TrackTools::FindMass2(CDSTrack *track, LocalTrack* bpc,double tof, double beammom, int pid_beam,double &beta_calc,double &mass2,double &tmptof)
{
  //  std::cout<<"-------------"<<std::endl;
  double param[5];
  TVector3 tmp;
  track->GetParameters(CID_CDC,param,tmp);   
  TVector3 cdhvtx=track->CDHVertex();
  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.4);
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
    if( !ELossTools::CalcHelixElossToVertexTGeo(param,vtxb2,mom,tmpmass,tmpl,tmptof,15.4)){
      //      std::cout<<" !!! error in calchelixelosstovertextgeo !!! "<<std::endl;
      return false;
    }
    //    std::cout<<tmptof<<"  "<<tmpl<<std::endl;
    tmptof+=time_beam;
    beta_calc=(cdc_dis)/(tof-tmptof)/(Const*100);
    mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
    //    std::cout<<"tmptof,beta,cdcdis,mass2  "<<tmptof<<"  "<<beta_calc<<"  "<<cdc_dis<<"  "<<mass2<<std::endl;
    count++;
    if(count>50)   return false;
  } while( TMath::Abs(tmpmass2-mass2)>0.0001);
  return true;
}

bool TrackTools::FindMass2(CDSTrack *track, pBeam *beam,double tof ,double &beta_calc,double &mass2,double &tmptof)
{
  //  std::cout<<"-------------"<<std::endl;
  double param[5];
  TVector3 tmp;
  track->GetParameters(CID_CDC,param,tmp);   
  TVector3 cdhvtx=track->CDHVertex();
  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.4);
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,pos1);
  TVector3 vtxb1,vtxb2;
  double dis;
  TrackTools::CalcLineHelixVertex(beam,track,vtxb1,vtxb2,dis);
  double mom,time_beam;
  TVector3 t0pos=beam->t0pos();
  double beammom=beam->mom();
  int pid_beam=beam->pid();
  ELossTools::CalcElossBeamTGeo(t0pos,vtxb1,beammom,parMass[pid_beam],mom,time_beam);
  mom = track->Momentum();
  int count =0;
  double tmpl,tmpmass2;
  do{
    tmpmass2=mass2;
    double tmpmass=sqrt(tmpmass2);
    if(tmpmass2<0.0005) tmpmass=0.0005;
    if( !ELossTools::CalcHelixElossToVertexTGeo(param,vtxb2,mom,tmpmass,tmpl,tmptof,15.4)){
      //      std::cout<<" !!! error in calchelixelosstovertextgeo !!! "<<std::endl;
      return false;
    }
    //    std::cout<<tmptof<<"  "<<tmpl<<std::endl;
    tmptof+=time_beam;
    beta_calc=(cdc_dis)/(tof-tmptof)/(Const*100);
    mass2=mom*mom*(1/(beta_calc*beta_calc)-1);
    //    std::cout<<"tmptof,beta,cdcdis,mass2  "<<tmptof<<"  "<<beta_calc<<"  "<<cdc_dis<<"  "<<mass2<<std::endl;
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
  TVector3 pos1=MathTools::CalcHelixPosatR(param,15.4);
  double cdc_dis=MathTools::CalcHelixArc(param,cdhvtx,pos1);
  //  TVector3 vtxb1,vtxb2;
  //  double dis;
  //  TrackTools::CalcLineHelixVertex(bpc,track,vtxb1,vtxb2,dis);
  double mom,tmpmom,time_beam;
  TVector3 t0pos=bpc->GetPosatZ(-110-0.5); 

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
    if(!track->CalcVertexTimeLength(t0pos,bpc->GetMomDir(),tmpmass,lpos,hpos,tmptof,tmpl)) return false;
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
    if(count>50)   return false;
  } while( TMath::Abs(tmpmass2-mass2)>0.0001);
  return true;
}

bool TrackTools::Calc2HelixVertex(CDSTrack *cds1, CDSTrack *cds2, TVector3 &vtx1, TVector3 &vtx2)
{
  double par1[5];
  double par2[5];
  TVector3 tmppos1,tmppos2;
  TVector3 tmp1,tmp2;
  int id,id2;
  double margin=1.;
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
      //      std::cout<<GeomTools::GetID(tmp1)<<"  "<<id<<"   "<<tmp1.Perp()<<"  "<<tmp1.Z()<<std::endl;
      if((GeomTools::GetID(tmp1)==id||GeomTools::IsSameVolumeHelix(par1,tmppos1,tmp1,margin))){	
	//	std::cout<<dis<<"  "<<tmpdis<<std::endl;
	//	std::cout<<GeomTools::GetID(tmp2)<<"  "<<id2<<"   "<<tmp2.Perp()<<"  "<<tmp2.Z()<<std::endl;
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
#if DEBUG
    std::cout<<"++++++++++++++++++++++++++++++++"<<tmpdis<<std::endl;
    cds1->Print();
    cds2->Print();
    exit(0);
#endif
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

int TrackTools::PID(double mom,double mass2)
{
  int ptype=-1;
  if(mom>0) {
    if(mass2<0.4) ptype=CDS_PiPlus;
    else if(mass2<2.2) ptype=CDS_Proton;
    else if(mass2<5.) ptype=CDS_Deuteron;
    else if(mass2<9.) ptype=CDS_Triton;
    else ptype=CDS_Other;
  }
  else{
    if(mass2<0.14) ptype=CDS_PiMinus;
    else if(mass2<0.5) ptype=CDS_Kaon;
    else ptype=CDS_Other;
  }
  return ptype;
}


int TrackTools::PIDcorr(double mom,double mass2)
{
  int ptype=CDS_Other;
  double fmom=fabs(mom);
  TF1 *f_sigma=new TF1("f_sigma","sqrt( 4*pow([3],4)*[0]*[0]*x*x+4*pow([3],4)*(1+pow([3]/x,2))*[1]*[1]+4*x*x*([3]*[3]+x*x)*[2]*[2])");

  double param[3][4]={{0.5332,0.09674,0.09533,piMass},
                      {0.000001,0.01164,0.10197,kpMass},
                      {0.1046,0.01037,0.000004,pMass}};
  
  double kpi_th=0;
  int n=0;
  double dis=999;
  double tmpmom[2];
  tmpmom[0]=0.2;
  tmpmom[1]=0.6;
  while(n<20 &&  dis>1e-4)
    {
      n++;
      double distmp[2]={999};
      for(int m=0;m<2;m++)
        {
          f_sigma->SetParameters(param[0][0],param[0][1],param[0][2],param[0][3]);
          double th_pi=piMass*piMass+2.5*f_sigma->Eval(tmpmom[m]);
          f_sigma->SetParameters(param[1][0],param[1][1],param[1][2],param[1][3]);
          double th_k=kpMass*kpMass-2.5*f_sigma->Eval(tmpmom[m]);
          distmp[m]=th_pi-th_k;
        }

      double  midmom=(tmpmom[0]+tmpmom[1])/2.;
      f_sigma->SetParameters(param[0][0],param[0][1],param[0][2],param[0][3]);
      double th_pi=piMass*piMass+2.5*f_sigma->Eval(midmom);
      f_sigma->SetParameters(param[1][0],param[1][1],param[1][2],param[1][3]);
      double th_k=kpMass*kpMass-2.5*f_sigma->Eval(midmom);
      kpi_th=th_k;
      dis=th_pi-th_k;
      if(dis*distmp[1]<0)
        {
          tmpmom[0]=midmom;
        }
      else  tmpmom[1]=midmom;
    }

  //  std::cout<<"th_kpi "<<kpi_th<<std::endl;
  if(mom>0) {
    if(mass2<0.25)
      {
	f_sigma->SetParameters(param[0][0],param[0][1],param[0][2],param[0][3]);
	double th1=piMass*piMass-2.5*f_sigma->Eval(fmom);
        double th2=piMass*piMass+2.5*f_sigma->Eval(fmom);
        if(th1<mass2 && mass2<th2) ptype=CDS_PiPlus;
      }
    else if(mass2<2.2)
      {
        f_sigma->SetParameters(param[2][0],param[2][1],param[2][2],param[2][3]);
        double th1=pMass*pMass-2.5*f_sigma->Eval(fmom);
        double th2=pMass*pMass+2.5*f_sigma->Eval(fmom);
        if(th1<mass2 && mass2<th2) ptype=CDS_Proton;
      }
    else if(mass2<5)
      {
        ptype=CDS_Deuteron;
      }
    else if(mass2<9.) ptype=CDS_Triton;
  }
  else if(mom<0) {
    if(mass2<kpi_th)
      {
        f_sigma->SetParameters(param[0][0],param[0][1],param[0][2],param[0][3]);
        double th1=piMass*piMass-2.5*f_sigma->Eval(fmom);
        double th2=piMass*piMass+2.5*f_sigma->Eval(fmom);
        if(th1<mass2 && mass2<th2) ptype=CDS_PiMinus;
      }
    else
      {
        f_sigma->SetParameters(param[1][0],param[1][1],param[1][2],param[1][3]);
        double th1=kpMass*kpMass-2.5*f_sigma->Eval(fmom);
	double th2=kpMass*kpMass+2.5*f_sigma->Eval(fmom);
        if(th1<mass2 && mass2<th2) ptype=CDS_Kaon;
      }
  }
  delete f_sigma;
  return ptype;
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

  double beta_calc,tofvtxcdc;
  double tof=cdhtime-beam->t0time();
  TrackTools::FindMass2(cdc,beam,tof,beta_calc,mass2,tofvtxcdc);
  int pid=TrackTools::PID(mom,mass2);
  cdc->SetPID(pid);

  double tmpl;
  if(ELOSS) cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],vtx2,vtx1,tofvtxcdc,tmpl,true);
  //  TrackTools::CalcLineHelixVertex(beam,cdc,vtx1,vtx2,dis,ELOSS);

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
  
  TVector3 pos1;
  cdc->GetParameters(CID_CDCCFRP,par,pos1);
  cdc->GetParameters(CID_CDC,par,vtx1);
  cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,pos1);
  double beta1=fabs(mom)/sqrt(cdsMass[pid]*cdsMass[pid]+mom*mom);
  double time_cdc=cdc_dis/beta1/(Const*100);
  double cdcdt=(tof-tofvtxcdc)-time_cdc;

  cdstrack->SetFL(cdc_dis);
  cdstrack->SetDt(cdcdt);
  cdstrack->SetChi2(cdc->Chi());

  return cdstrack;
}

pCDS *TrackTools::Calc2HelixAll(pBeam *beam, pCDS* cds1, pCDS* cds2,CDSTrackingMan *trackMan, bool ELOSS)
{
  CDSTrack *track1=trackMan->Track(cds1->id());	
  CDSTrack *track2=trackMan->Track(cds2->id());	
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
