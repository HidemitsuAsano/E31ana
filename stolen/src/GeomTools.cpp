#include "GeomTools.h"
#include "MathTools.h"

#define DEBUG 0

void GeomTools::MakeGeometry(ConfMan *conf) //run49c
{
  std::cout<<"============================="<<std::endl;
  //  std::cout<<"start Makeing TGeo"<<std::endl;
  new TGeoManager("K1.8BR spectrometer", "geometry for K1.8BR spectrometer");
  ConstructMaterial();
  ConstructHadronHall(conf);
  ConstructHodoscopes(conf);
  ConstructChambers(conf);
  ConstructTarget(conf);
  //--- close the geometry
  gGeoManager->CloseGeometry();
  gGeoManager->CheckOverlaps(0.00001, "d");
  gGeoManager->PrintOverlaps();
  std::cout<<"============================="<<std::endl;
}

void GeomTools::MakeGeometryHeates(ConfMan *conf)
{
  std::cout<<"============================="<<std::endl;
  //  std::cout<<"start Makeing TGeo"<<std::endl;
  if(gGeoManager) delete gGeoManager;
  new TGeoManager("HEATES spectrometer at J-PARC", "geometry for HEATES spectrometer at J-PARC");
  ConstructMaterial();
  ConstructShape(CID_Hall,conf);
  ConstructShape(CID_T0,conf);
  ConstructShape(CID_AC,conf);

  ConstructShape(CID_DegC,conf);
  ConstructShape(CID_DegCu,conf);
  ConstructShape(CID_E0,conf);
  //  ConstructShape(CID_BPC,conf);
  ConstructShape(CID_LC1,conf);

  ConstructShape(CID_TarSys,conf);
  ConstructShape(CID_TarChm,conf);
  ConstructShape(CID_TarCell,conf);
  ConstructShape(CID_CellTube,conf);
  ConstructShape(CID_CellFlange,conf);
  if(0){
    ConstructShape(CID_CellWindow,conf);
  }else{
    ConstructShape(CID_RadS,conf);
  }
  //		 ConstructShape(CID_BeamProf,conf);
  
  ConstructShape(CID_TESSYS,conf);
  ConstructShape(CID_TES,conf);
  ConstructShape(CID_TESCol,conf);
  ConstructShape(CID_TESSi,conf);
  ConstructShape(CID_TESGrid,conf);
  ConstructShape(CID_TESAntico,conf);
  ConstructShape(CID_TESHS,conf);
  ConstructShape(CID_RS60K,conf);
  ConstructShape(CID_RS1K,conf);
  ConstructShape(CID_RS50mK,conf);
  ConstructShape(CID_RSwindow,conf);
  ConstructShape(CID_DEF,conf);
  ConstructShape(CID_BVC,conf);

  //--- close the geometry
  gGeoManager->CloseGeometry();
  gGeoManager->CheckOverlaps(0.00001, "d");
  gGeoManager->PrintOverlaps();
  std::cout<<"============================="<<std::endl;
}

void GeomTools::MakeGeometryPSI(ConfMan *conf)
{
  std::cout<<"============================="<<std::endl;
  //  std::cout<<"start Makeing TGeo"<<std::endl;
  if(gGeoManager) delete gGeoManager;
  new TGeoManager("HEATES spectrometer at PSI", "geometry for HEATES spectrometer at PSI");
  ConstructMaterial();
  ConstructShape(CID_Hall,conf);
  ConstructShape(CID_BC1,conf);
  ConstructShape(CID_BC2,conf);
  ConstructShape(CID_BC3,conf);
  ConstructShape(CID_BC4,conf);
  
  ConstructShape(CID_DegC,conf);
  ConstructShape(CID_DegCu,conf);

  ConstructShape(CID_CellFlange,conf);
  ConstructPSICarbonTarget(CID_Target,conf);
  
  ConstructShape(CID_TESSYS,conf);
  ConstructShape(CID_TESChamber,conf);
  ConstructShape(CID_TESWindow,conf);
  ConstructShape(CID_RS60K,conf);
  ConstructShape(CID_RS1K,conf);
  ConstructShape(CID_RS50mK,conf);
  ConstructShape(CID_TES,conf);
  ConstructShape(CID_TESCol,conf);
  ConstructShape(CID_TESSi,conf);
  ConstructShape(CID_TESGrid,conf); 
  ConstructShape(CID_TESHS,conf);
  ConstructShape(CID_TESAntico,conf);

  //--- close the geometry
  gGeoManager->CloseGeometry();
  gGeoManager->CheckOverlaps(0.00001, "d");
  gGeoManager->PrintOverlaps();
  std::cout<<"============================="<<std::endl;
}


bool GeomTools::IsSameVolume(const TVector3 &pos1, const TVector3 &pos2,double margin)
{
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  double length=gGeoManager->GetStep();
  if( (pos2-pos1).Mag()<(length+margin) ) return true;
  return false;
}

bool GeomTools::IsSameVolumeHelix(const double param[5],const TVector3 &pos1, const TVector3 &pos2,double margin)
{
  TVector3 tmppos;
  double tmpl;
  TString mat;
  int id;
  if(!GeomTools::HelixStepToNextVolume(param,pos1,tmppos,tmpl,mat,id)) return false;
  double tmpl2=MathTools::CalcHelixArc(param,pos1,pos2);
  if(tmpl2<(tmpl+margin)) true;
  return false;
}

bool GeomTools::HelixStepToNextVolume(const double param[5],const TVector3 &in,TVector3 &pos2,double &length, TString &mat, int &id)
{
  double defaultstep=1*cm;
  double step=defaultstep;  
  length=0;
  TVector3 pos1=in;
  pos2=MathTools::CalcHelixStep(param,pos1,step);
  while(step>9*um){
    if(GeomTools::IsSameVolume(pos1,pos2)){
      length+=step;
      pos1=pos2;
      pos2=MathTools::CalcHelixStep(param,pos1,step);
    }else{
      step/=10.;
      pos2=MathTools::CalcHelixStep(param,pos1,step);
    }    
    if(length>100*cm) return false;
  }
  step*=10;
  pos2=MathTools::CalcHelixStep(param,pos1,step);
  double tmpstep;
  GeomTools::CrossBoundary(pos1,pos2,step,mat,tmpstep,id);
  id/=1000;
  length+=step;
  //  std::cout<<length<<std::endl;
  return true;
}
bool GeomTools::StepToNextVolume(const TVector3 &pos1,const TVector3 &pos2,double &length, TString &newmat)
{
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  double step;
  //  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  length=gGeoManager->GetStep()+1*um;
  if(length>(pos2-pos1).Mag()) return false;
#if 0
  const double *pos=gGeoManager->GetCurrentPoint();
  for(int i=0;i<3;i++) std::cout<<pos1[i]<<"  ";
  for(int i=0;i<3;i++) std::cout<<pos2[i]<<"  ";
  std::cout<<length<<std::endl;
  //  std::cout<<std::endl;
#endif
  //  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep();
  if(!node){
#if DEBUG
    std::cout<<"newnode cannot be defined"<<std::endl;
#endif
    return false;
  }
  gGeoManager->CdNext();

  if(node->GetMedium())
    newmat = node->GetMedium()->GetName();
  else
    std::cout<<"no medium in "<<node->GetNumber()<< "  "<<node->GetName()<<std::endl;
  if(newmat=="Concrete"){
    pos1.Print();
    pos2.Print();
    std::cout<<"Concrete !!!"<<std::endl;
    return false;
  }
    
  return true;
}

bool GeomTools::CrossBoundary(const TVector3 &pos1,const TVector3 &pos2,const double &step, TString &mat, double &tmpstep, int &newid){
  TVector3 dir=(pos2-pos1).Unit();
  gGeoManager->InitTrack(pos1.X(),pos1.Y(),pos1.Z(),dir.X(),dir.Y(),dir.Z());
  TGeoNode* node=gGeoManager->FindNextBoundaryAndStep(step,true);
  newid=node->GetNumber();
  mat = node->GetMedium()->GetName();
  tmpstep=gGeoManager->GetStep();
  return true;
}

double GeomTools::CalcLengthinFiducial(const TVector3 &pos,const TVector3 &dir){
  TVector3 dir2=(dir-pos).Unit();
  TVector3 tmppos=pos;
  TString newmat;
  double length;
  double sumlength=0;
  bool FIDUCIAL=false;
  TVector3 in,out;
  while(GeomTools::StepToNextVolume(tmppos,dir,length,newmat)){
    tmppos+=dir2*length;
    int id=GetID(tmppos);
    if(FIDUCIAL&&id!=CID_Fiducial){
      out=tmppos;
      break;
    }
    if(id==CID_Fiducial){
      in=tmppos;
      FIDUCIAL=true;
    }
    if(tmppos.Z()>10) return 0;
  }
  return (out-in).Mag();
}

void GeomTools::ConstructPSICarbonTarget(Int_t CID,ConfMan* conf,TGeoVolume *mother){
  DetectorList* dlist=DetectorList::GetInstance();

  Int_t nsegments=dlist->GetNsegs(CID);
  TString name=dlist->GetName(CID);
  //  std::cout<<CID<<"  "<<name<<std::endl;
  if(name=="NULL"){
    std::cout<<"No counter registerd for CID = "<<CID<<std::endl;
    exit(1);
  }
  //  std::cout<<nsegments<<"  "<<name<<"  "<<dlist->GetMaterial(CID)<<std::endl;
  TGeoMedium *medium1= gGeoManager->GetMedium("Air");
  TGeoMedium *medium2= gGeoManager->GetMedium(dlist->GetMaterial(CID).data());
  if(nsegments==0) medium1=medium2;  

  double param[20];
  int seg=0;
  TGeoVolume* assembly=0;
  if(GeomTools::GetParam(CID,seg,conf,param)){
    double zplane[]={-40,40};
    double rin[]={20,7.5};
    double rout[]={21,8.5};
    assembly = gGeoManager->MakePgon("PGON",medium1, 0.0,360.0,4,2);
    TGeoPgon *pgon = (TGeoPgon*)(assembly->GetShape());
    for(int i=0;i<2;i++)
      pgon->DefineSection(i,zplane[i]/10.,rin[i]/10.,rout[i]/10.);

    if(param[10]==-1){
      gGeoManager->SetTopVolume(assembly);
    }else{
      if(!mother)
	mother=gGeoManager->GetVolume(dlist->GetName((int)param[10]).data());
      if(!mother)
	mother=gGeoManager->GetTopVolume();
      //    std::cout<<"mother= "<<CID<<"  "<<param[10]<<"  "<<mother->GetName()<<std::endl;
      if(mother&&assembly){
	TGeoCombiTrans* assembly_trans=GeomTools::MakeTrans(param);
	assembly->SetLineColor(kBlack);
	assembly->SetTransparency(1);
	mother->AddNode(assembly, CID*1000, assembly_trans);
      }
    }
  }
}

void GeomTools::ConstructHadronHall(ConfMan* conf)
{
  // We need to construct hadronhall first 
  ConstructShape(CID_Hall,conf);
  ConstructShape(CID_Doraemon,conf);
  ConstructShape(CID_USWK,conf);
  ConstructMan();
  ConstructShape(CID_Floor,conf);
  ConstructShape(CID_BeamDump,conf);
  ConstructShape(CID_SideDump,conf);
  ConstructShape(CID_NShield,conf);
  ConstructShape(CID_SideCon,conf);
  ConstructShape(CID_DoorCon,conf);
}
void GeomTools::ConstructTarget(ConfMan *conf)
{
  ConstructShape(CID_CDCCFRP,conf);
  ConstructShape(CID_CDCMylar,conf);
  ConstructShape(CID_TarSys,conf);
  ConstructShape(CID_RadS,conf);
  ConstructShape(CID_TarCFRP,conf);
  ConstructShape(CID_TarCap,conf);
  ConstructShape(CID_TarCell,conf);  
  ConstructShape(CID_TarRing,conf);
  ConstructShape(CID_Target,conf);
  ConstructShape(CID_CellBe,conf);
  ConstructShape(CID_CellAlBe,conf);
  ConstructShape(CID_CellRing,conf);
  ConstructShape(CID_BShield,conf);
  ConstructShape(CID_BFrange,conf);
  ConstructShape(CID_Fiducial,conf);
}
void GeomTools::ConstructHodoscopes(ConfMan* conf){
  //  ConstructShape(CID_BHD,conf);
  ConstructShape(CID_T0,conf);
  ConstructShape(CID_AC,conf);
  ConstructShape(CID_BPD,conf);
  ConstructShape(CID_DEF,conf);
  ConstructShape(CID_BVC,conf);
  ConstructShape(CID_CVC,conf);
  ConstructShape(CID_NC,conf);
  ConstructShape(CID_PC,conf);
  //  ConstructShape(hadron_hall, CID_BD);
  ConstructShape(CID_CDH,conf);
  ConstructShape(CID_IH,conf);
}
void GeomTools::ConstructChambers(ConfMan* conf){
  //  ConstructShape( CID_BLC1a, conf);
  //  ConstructShape( CID_BLC1b);
  ConstructShape( CID_BLC2a, conf);
  ConstructShape( CID_BLC2b, conf);
  ConstructShape( CID_FDC1,  conf);
  ConstructShape( CID_BPC,   conf);  
  ConstructShape( CID_CDC,   conf);
}

bool GeomTools::GetParam(const int &cid,const int &seg, ConfMan* conf, double *param)
{
  bool status=false;
  if(cid==CID_CDC)
    status=conf->GetCDCWireMapManager()->GetGParam(param);
  else if(DetectorList::GetInstance()->IsChamber(cid))
    status=conf->GetBLDCWireMapManager()->GetGParam(cid,0,param);
  else
    status=conf->GetGeomMapManager()->GetParam(cid,seg,param);
  return status;
}

TGeoVolume* GeomTools::ConstructShape(Int_t CID,ConfMan *conf,TGeoVolume *mother){
  //====================================================
  //--- define detector based on CID
  //====================================================

  DetectorList* dlist=DetectorList::GetInstance();

  Int_t nsegments=dlist->GetNsegs(CID);
  TString name=dlist->GetName(CID);
  //  std::cout<<CID<<"  "<<name<<std::endl;
  if(name=="NULL"){
    std::cout<<"No counter registerd for CID = "<<CID<<std::endl;
    return 0;
    exit(1);
  }
  //  std::cout<<nsegments<<"  "<<name<<"  "<<dlist->GetMaterial(CID)<<std::endl;
  TGeoMedium *medium1= gGeoManager->GetMedium("Air");
  TGeoMedium *medium2= gGeoManager->GetMedium(dlist->GetMaterial(CID).data());
  if(nsegments==0) medium1=medium2;  

  double param[20];
  int seg=0;
  TGeoVolume* assembly=0;
  if(GeomTools::GetParam(CID,seg,conf,param)){
    assembly=GeomTools::MakeShape(param,name,medium1);
#if 0
    for(int i=0;i<11;i++)
      std::cout<<"  "<<param[i];
    std::cout<<std::endl;
#endif
    if(param[10]==-1){
      gGeoManager->SetTopVolume(assembly); return assembly;
    }else{
      if(!mother)
	mother=gGeoManager->GetVolume(dlist->GetName((int)param[10]).data());
      if(!mother)
	mother=gGeoManager->GetTopVolume();
      //    std::cout<<"mother= "<<CID<<"  "<<param[10]<<"  "<<mother->GetName()<<std::endl;
      if(mother&&assembly){
	TGeoCombiTrans* assembly_trans=GeomTools::MakeTrans(param);
	assembly->SetLineColor(kBlack);
	assembly->SetTransparency(1);
	mother->AddNode(assembly, CID*1000, assembly_trans);
	if(nsegments==0) return assembly;
      }
    }
  }
  if(!assembly){
    if(!GeomTools::GetParam(CID,1,conf,param)) return 0;
    assembly=gGeoManager->GetVolume(dlist->GetName((int)param[10]).data());
  }
  if(!assembly)
    return 0;
  
  TGeoVolume *segment;
  for (Int_t i=1; i<=nsegments; i++){
    TString segname=Form("%s_segment%d",name.Data(),i);
    //    std::cout<<segname<<std::endl;
    if(!GeomTools::GetParam(CID,i,conf,param)) continue;
    segment = GeomTools::MakeShape(param,segname,medium2); 
    TGeoCombiTrans *segment_trans = GeomTools::MakeTrans(param);
    segment->SetLineColor(kBlue);
    segment->SetTransparency(1);
    assembly->AddNode(segment, CID*1000+i, segment_trans);
  }
  return assembly;
}

TGeoVolume* GeomTools::MakeShape(double *param,const TString &name, TGeoMedium* medium ){
  Double_t rmin = param[6]*cm;
  Double_t rmax = param[7]*cm;
  Double_t phi  = param[8]*deg;
  Double_t z    = param[9]*cm/2.0;
  TGeoVolume *shape=0;
  if(z>0){
    shape= gGeoManager->MakeTubs(name, medium,rmin, rmax, z,0,phi);
  }else{
    Double_t x = param[6]*cm/2.0;
    Double_t y = param[7]*cm/2.0;
    z    = param[8]*cm/2.0;
    if(z>0){
      shape= gGeoManager->MakeBox(name, medium,x,y,z);
    }else{
      return 0;
    }
  }    
  return shape;
}

TGeoCombiTrans* GeomTools::MakeTrans(double *param){
  Double_t pos_x = param[0]*cm;
  Double_t pos_y = param[1]*cm;
  Double_t pos_z = param[2]*cm;
  if(param[9]!=0&&param[8]<90){
    Double_t pos_r   = param[0]*cm;
    Double_t pos_phi = param[1]*Deg2Rad;
    Double_t pos_z   = param[2]*cm;
    TVector3 pos(pos_r,0,pos_z);
    pos.RotateZ(pos_phi);
    pos_x = pos.X();
    pos_y = pos.Y();
    pos_z = pos.Z();
  }
  TGeoRotation *rot   = new TGeoRotation();
  rot->RotateX(param[3]);
  rot->RotateY(param[4]);
  rot->RotateZ(param[5]);
  TGeoCombiTrans *trans = new TGeoCombiTrans(pos_x,pos_y,pos_z,rot);
  return trans;
}

void GeomTools::ConstructMaterial(){
  //====================================================
  //--- define materials
  //====================================================
  TGeoMaterial *mat[20];
  for(int i=0;i<20;i++)
    mat[i]= new TGeoMaterial(Form("mat%d",i),    1.008,  1,  1.e-25);
  
  new TGeoMedium("Vacuum",       0, mat[0]);
  new TGeoMedium("LHelium-3",    1, mat[1]);
  new TGeoMedium("Beryllium",    2, mat[2]);
  new TGeoMedium("Air",          3, mat[3]);
  new TGeoMedium("CFRP",         4, mat[4]);
  new TGeoMedium("Mylar",        5, mat[5]);
  new TGeoMedium("Scinti",       6, mat[6]);
  new TGeoMedium("Plastic",      6, mat[6]);
  new TGeoMedium("Aluminum",     7, mat[7]);
  new TGeoMedium("Concrete",     8, mat[8]);
  new TGeoMedium("Iron",         9, mat[9]);
  new TGeoMedium("Tungsten",     10,mat[10]);
  new TGeoMedium("CDCGas",       11,mat[11]);
  new TGeoMedium("BLDCGas",      12,mat[12]);
  new TGeoMedium("Aerogel",      13,mat[13]);
  new TGeoMedium("AlBeMet",      14,mat[14]);
  new TGeoMedium("G4_STAINLESS-STEEL",      15,mat[15]);
  new TGeoMedium("Graphite",     16,mat[16]);
  new TGeoMedium("G4_Bi",        17,mat[17]);
  new TGeoMedium("LHydrogen",    18,mat[18]);
  new TGeoMedium("LDeuterium",   19,mat[19]);
}

void GeomTools::ConstructMan(){
  TGeoVolume *mother=gGeoManager->GetTopVolume();
  TGeoMedium* Air=gGeoManager->GetMedium("Air");
  //====================================================
  //--- define operator
  //====================================================
  Double_t beam_line_height = 2.0*m;
  Double_t operator_assembly_x = 60*cm/2.0;
  Double_t operator_assembly_y = 180*cm/2.0;
  Double_t operator_assembly_z = 30.0*cm/2.0;
  TGeoVolume *operator_assembly = gGeoManager->MakeBox("operator_assembly", 
						      Air, 
						      operator_assembly_x,
						      operator_assembly_y,
						      operator_assembly_z);
  Double_t operator_assembly_pos_x =  3.0*m;
  Double_t operator_assembly_pos_y = -beam_line_height + operator_assembly_y;
  Double_t operator_assembly_pos_z =  7.0*m;
  Double_t operator_assembly_angle =  90.0*degree; 
  TGeoRotation *operator_assembly_rot = new TGeoRotation();
  operator_assembly_rot->RotateY(operator_assembly_angle);
  TGeoCombiTrans *operator_assembly_trans = new TGeoCombiTrans(operator_assembly_pos_x,
							       operator_assembly_pos_y,
							       operator_assembly_pos_z,
							       operator_assembly_rot);
  operator_assembly->SetLineColor(kBlue);
  mother->AddNode(operator_assembly, 0, operator_assembly_trans);
  Double_t operator_head_rmin     = 0.0*cm;
  Double_t operator_head_rmax     = 15.0*cm;
  Double_t operator_head_thetamin = 0.0*degree;
  Double_t operator_head_thetamax = 180.0*degree;
  Double_t operator_head_phimin   = 0.0*degree;
  Double_t operator_head_phimax   = 360.0*degree;
  TGeoVolume *operator_head = gGeoManager->MakeSphere("operator_head", 
						     Air,
						     operator_head_rmin,
						     operator_head_rmax,
						     operator_head_thetamin,
						     operator_head_thetamax,
						     operator_head_phimin,
						     operator_head_phimax);
  Double_t operator_head_pos_x =  0.0*m;
  Double_t operator_head_pos_y =  0.75*m;
  Double_t operator_head_pos_z =  0.0*m;
  Double_t operator_head_angle =  0.0*degree; 
  TGeoRotation *operator_head_rot = new TGeoRotation();
  operator_head_rot->RotateY(operator_head_angle);
  TGeoCombiTrans *operator_head_trans = new TGeoCombiTrans(operator_head_pos_x,
							   operator_head_pos_y,
							   operator_head_pos_z,
							   operator_head_rot);
  operator_head->SetLineColor(kBlue);
  operator_assembly->AddNode(operator_head, 0, operator_head_trans);
  Double_t operator_body_x = 40.0*cm/2.0;
  Double_t operator_body_y = 50.0*cm/2.0;
  Double_t operator_body_z = 10.0*cm/2.0;
  TGeoVolume *operator_body = gGeoManager->MakeBox("operator_body", 
						  Air,
						  operator_body_x,
						  operator_body_y,
						  operator_body_z);
  Double_t operator_body_pos_x =  0.0*m;
  Double_t operator_body_pos_y =  0.35*m;
  Double_t operator_body_pos_z =  0.0*m;
  Double_t operator_body_angle =  0.0*degree; 
  TGeoRotation *operator_body_rot = new TGeoRotation();
  operator_body_rot->RotateY(operator_body_angle);
  TGeoCombiTrans *operator_body_trans = new TGeoCombiTrans(operator_body_pos_x,
							   operator_body_pos_y,
							   operator_body_pos_z,
							   operator_body_rot);
  operator_body->SetLineColor(kBlue);
  operator_assembly->AddNode(operator_body, 0, operator_body_trans);
  Double_t operator_arm_rmin = 0.0*cm;
  Double_t operator_arm_rmax = 5.0*cm;
  Double_t operator_arm_z    = 70.0*cm/2.0;
  TGeoVolume *operator_arm = gGeoManager->MakeTube("operator_arm", 
						  Air,
						  operator_arm_rmin,
						  operator_arm_rmax,
						  operator_arm_z);
  for(Int_t i=0; i<2; i++){
    Double_t operator_arm_pos_x;
    if( i==0 )
      operator_arm_pos_x =  20.0*cm + operator_arm_rmax;
    if( i==1 )
      operator_arm_pos_x = -20.0*cm - operator_arm_rmax;
    Double_t operator_arm_pos_y =  0.25*m;
    Double_t operator_arm_pos_z =  0.0*m;
    Double_t operator_arm_angle =  90.0*degree; 
    TGeoRotation *operator_arm_rot = new TGeoRotation();
    operator_arm_rot->RotateX(operator_arm_angle);
    TGeoCombiTrans *operator_arm_trans = new TGeoCombiTrans(operator_arm_pos_x,
							    operator_arm_pos_y,
							    operator_arm_pos_z,
							    operator_arm_rot);
    operator_arm->SetLineColor(kBlue);
    operator_assembly->AddNode(operator_arm, i, operator_arm_trans);
  }
  Double_t operator_leg_rmin = 0.0*cm;
  Double_t operator_leg_rmax = 10.0*cm;
  Double_t operator_leg_z    = 100.0*cm/2.0;
  TGeoVolume *operator_leg = gGeoManager->MakeTube("operator_leg", 
						  Air,
						  operator_leg_rmin,
						  operator_leg_rmax,
						  operator_leg_z);
  for(Int_t i=0; i<2; i++){
    Double_t operator_leg_pos_x;
    if( i==0 )
      operator_leg_pos_x =  operator_leg_rmax;
    if( i==1 )
      operator_leg_pos_x = -operator_leg_rmax;
    Double_t operator_leg_pos_y =  -40*cm;
    Double_t operator_leg_pos_z =   0.0*m;
    Double_t operator_leg_angle =   90.0*degree; 
    TGeoRotation *operator_leg_rot = new TGeoRotation();
    operator_leg_rot->RotateX(operator_leg_angle);
    TGeoCombiTrans *operator_leg_trans = new TGeoCombiTrans(operator_leg_pos_x,
							    operator_leg_pos_y,
							    operator_leg_pos_z,
							    operator_leg_rot);
    operator_leg->SetLineColor(kBlue);
    operator_assembly->AddNode(operator_leg, i, operator_leg_trans);
  }
}

