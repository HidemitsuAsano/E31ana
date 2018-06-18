#ifndef PARTICLE_h
#define PARTICLE_h 1

#include <vector>
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "GlobalVariables.h"
#include "HodoscopeLikeHit.h"
#include "ConfMan.h"
#ifndef ROOT_TObject
#include "TObject.h"
#endif


class pBeam : public TObject
{
 public:
  pBeam();
  virtual ~pBeam(){};

 private:
  int PID; //geant particle code
  //  int TrigID[20];
  int T0seg;
  int BHDseg;
  double T0time;
  double BHDX;
  double tofbhdt0;
  TVector3 Vertex;
  TVector3 T0Pos;

  bool BPCTRACK, BLCTRACK, MOM;

  TVector3 BPCPos;//at FF
  TVector3 BPCDir;
  TVector3 BLCPos;//at FF
  TVector3 BLCDir;
  double Momentum;

 public:
  void SetPID(const int &pid) { PID=pid; }
  void SetT0seg( const int &seg) { T0seg = seg; }
  void SetT0Time( const double &time) { T0time = time; }
  void SetT0Pos( const TVector3 &pos) { T0Pos = pos; }
  void SetBHDseg( const int &seg) { BHDseg = seg; }
  void SetBHDT0TOF( const double &tof) { tofbhdt0 = tof; }
  void SetBHDX( const double &x ) { BHDX =x; }

  void SetVertex(TVector3 vec) { Vertex = vec; }
  
  void SetBPCPos(TVector3 vec) { BPCPos = vec; BPCTRACK = true; }
  void SetBPCDir(TVector3 vec) { BPCDir = vec; BPCTRACK = true; }
  void SetBLCPos(TVector3 vec) { BLCPos = vec; BLCTRACK = true; }
  void SetBLCDir(TVector3 &vec) { BLCDir = vec; BLCTRACK = true; }
  
  void SetBPCPos(const double &x, const double &y, const double &z)
  { BPCPos = TVector3(x,y,z); BPCTRACK = true; }
  void SetBPCDir(const double &x, const double &y, const double &z)
  { BPCDir = TVector3(x,y,z); BPCTRACK = true; }
  void SetBLCPos(const double &x, const double &y, const double &z)
  { BLCPos = TVector3(x,y,z); BLCTRACK = true; }
  void SetBLCDir(const double &x, const double &y, const double &z)
  { BLCDir = TVector3(x,y,z); BLCTRACK = true; }
  
  void SetMomentum(const double &mom) { Momentum = mom; MOM = true; } 

  int pid() const { return PID; }
  int t0seg() const { return T0seg; }
  double t0time() const { return T0time; }
  TVector3 t0pos() const { return T0Pos; }
  int bhdseg() const { return BHDseg; }
  double bhdt0tof() const { return tofbhdt0; }
  double bhdx() const { return BHDX; }

  double CalcVertexTime(const TVector3 &vertex);
  double CalcVertexMom( const TVector3 &vertex);

  TLorentzVector GetLorentzVector() 
    { 
    TLorentzVector lv;
    lv.SetVectM(BPCDir.Unit()*Momentum,parMass[PID]);
    return lv;
    }

  TLorentzVector GetLorentzVector(const TVector3 &vertex) 
  { 
    TLorentzVector lv;
    lv.SetVectM(BPCDir.Unit()*CalcVertexMom(vertex),parMass[PID]);
    return lv;
  }
  
  TVector3 bpcpos() { return BPCPos; }
  TVector3 bpcdir() { return BPCDir; }
  TVector3 blcpos() { return BLCPos; }
  TVector3 blcdir() { return BLCDir; }
  TVector3 vertex() { return Vertex; }

  void bpcpos(const double &z,double &x, double &y)
    {
      x = BPCPos.X() + z* BPCDir.X() / BPCDir.Z();
      y = BPCPos.Y() + z* BPCDir.Y() / BPCDir.Z();
    }
  void bpcdir(double &x, double &y)
    {
      x = BPCDir.X() / BPCDir.Z();
      y = BPCDir.Y() / BPCDir.Z();
    }
  void blcpos(const double &z,double &x, double &y)
    {
      x = BLCPos.X() + z* BLCDir.X() / BLCDir.Z();
      y = BLCPos.Y() + z* BLCDir.Y() / BLCDir.Z();
    }
  void blcdir(double &x, double &y)
    {
      x = BLCDir.X() / BLCDir.Z();
      y = BLCDir.Y() / BLCDir.Z();
    }
  double mom() const { return Momentum; }
  double mass() const { return parMass[PID]; }
  double beta() const { return mom()/sqrt(mass()*mass()+mom()*mom()); }
  bool isbpc() const { return BPCTRACK; }
  bool isblc() const { return BLCTRACK; }
  bool ismom() const { return MOM; }
  
  ClassDef( pBeam, 1 );
};

class pNC : public TObject
{
 public:
  pNC();
  pNC(int seg,double time, TVector3 pos);
  virtual ~pNC(){};
  
 protected:
  int Segment;
  int PID;
  TVector3 HitPos;
  TVector3 Mom;
  double Time;
  double Energy;
  double Beta;
  double Mass;
  double TOF;
  double FLength;

 public:
  void SetSegment( const int &seg ) { Segment=seg; }
  void SetHitPosition( const double &x,const double &y, const double &z ) { HitPos.SetXYZ(x,y,z); }
  void SetHitPosition( const TVector3 &pos ) { HitPos=pos; }
  void SetTime(const double &time ) { Time = time; }
  void SetBeta(const double &beta ) { Beta = beta; }
  void SetEnergy(const double &ene ) { Energy = ene; }
  void SetMomentum(TVector3 &mom) { Mom = mom; }

  TLorentzVector GetLorentzVector() 
  { 
    TLorentzVector lv;
    //    TVector3 tmpmom(Mom.X(),Mom.Y(), Mom.Z());
    lv.SetVectM(Mom,mass());
    return lv; 
  }
  
  int seg() const { return Segment; }
  TVector3 hitpos() { return HitPos; }
  double time() const { return Time; }
  double energy() const { return Energy; } 
  TVector3 mom() const { return Mom; }
  double beta() const { return Beta; }
  double mass() const { return Mass; };
  double tof()     const { return TOF; }
  double fl()      const { return FLength; }


  int pid() const { return PID; }
  bool isneutron() const { return PID==F_Neutron; }
  bool isgamma()   const { return PID==F_Gamma; }

  void CalcMom(pBeam* beam, const TVector3 &vertex,double offs=0);

  ClassDef( pNC, 1 );
};

class pPC : public pNC
{
 public:
  pPC();
  pPC(int seg,double time, TVector3 pos);
  virtual ~pPC(){};
  
 private:
  std::vector <HodoscopeLikeHit> BVCHit;
  TVector3 FDC1Pos;
  double Angle;
  double Radius;
  double Mass2;
  /* static const double magUSWK=0.9185; */
  /* static const double offsUSWK=0.042; */
  //  static const double magUSWK=1.0503;
  //  static const double offsUSWK=-0.003;
  static const double magUSWK=0.9782;
  static const double offsUSWK=0.023;
 public:
  void SetFDC1Pos( const TVector3 &pos ) { FDC1Pos=pos; }
  void SetPID(int id,double m) { PID=id; Mass=m; }
  void SetBVCHit(const HodoscopeLikeHit &hit) { BVCHit.push_back(hit); }
  
  TVector3 fdc1pos()     { return FDC1Pos; }
  double r()       const { return Radius; }
  double mass2()   const { return Mass2; }
  double momr()   const { return r()*magUSWK*Const/100.+offsUSWK; }
  double momUSWK(const TVector3 &vertex);
  double angle()   const { return Angle; }
  int nbvchit()    const { return BVCHit.size(); } 
  const HodoscopeLikeHit &bvchit( int i ){ return BVCHit[i]; }

  void CalcMom(pBeam* beam, const TVector3 &vertex,int pid=-1,bool SIM=false);
  void CalcMomBVC(const TVector3 &vertex,int pid=-1,bool SIM=false);
  //  void CalcMom(pBeam* beam);

  ClassDef( pPC, 1 );
};

class pCDS : public TObject
{
 public:
  pCDS();
  virtual ~pCDS(){};

 private:
  int trackID;
  int daughterID1;
  int daughterID2;

  int PID;
  int CombID;

  double Momentum;
  double CMomentum;
  double Beta;
  double Gamma;
  double TOF;
  double Mass;
  double Mass2;
  double PDGMass;
  double VertexDistance;
  double VertexBeamDistance;
  double ProductBeamDCA;
  double Param[5];

  TVector3 Vertex;
  TVector3 VertexCDC;
  TVector3 VertexBeam;
  TVector3 MomDir;

  std::vector<int> CDHseg;
  std::vector<int> IHseg;

  double OpenAngle;
  double AngleLab;
  double FlightLength;
  double Dt;
  double Chi;

 public:
  void SetTrackID(const int &id) { trackID = id; }
  void SetDaughterID1(const int &id) { daughterID1 = id; }
  void SetDaughterID2(const int &id) { daughterID2 = id; }

  void SetCombID(const int &id) { CombID = id; }
  void SetPID(const int &pid) { PID = pid; PDGMass=cdsMass[pid]; }
  void SetMomentum(const double &mom) { Momentum = mom; }
  void SetRawMomentum(const double &mom) { CMomentum = mom; }
  void SetMass(const double &mass) { Mass = mass; }
  void SetMass2(const double &mass) { Mass2 = mass; }
  void SetBeta(const double &beta) { Beta= beta; }
  void SetTOF(const double &tof) { TOF = tof; }
  void SetFL(const double &fl) { FlightLength = fl; }
  void SetDt(const double &dt) { Dt=dt; }
  void SetChi2(const double &chi) { Chi=chi; }

  void SetVertex(     TVector3 vec) { Vertex     =vec; }
  void SetVertexCDC(  TVector3 vec) { VertexCDC  =vec; }
  void SetVertexBeam( TVector3 vec) { VertexBeam =vec; }
  void SetMomDir(     TVector3 vec) { MomDir     =vec; }
  void SetVDis(  const double &vdis ){ VertexDistance     =vdis; }
  void SetVBDis( const double &vdis ){ VertexBeamDistance =vdis; }
  void SetPBDCA( const double &vdis ){ ProductBeamDCA     =vdis; }

  void SetCDHSeg(const int &seg ) { CDHseg.push_back(seg); }
  void SetIHSeg( const int &seg ) { IHseg.push_back(seg); }

  void SetAngleLab( const double &ang ){ AngleLab =ang; }
  void SetOA(       const double &ang ){ OpenAngle=ang; }
  void SetParameters( double *par ) { for(int i=0;i<5;i++) Param[i] =par[i]; }

  int id() const { return trackID; }
  int daughter1() const { return daughterID1; }
  int daughter2() const { return daughterID2; }

  int comb()const { return CombID; }
  int pid() const { return PID; }
  int pid2() const { return daughterID2; }
  double mom() const { return Momentum; }
  double rawmom() const { return CMomentum; }
  double mom(ConfMan *conf);
  double mass() const { return Mass ; }
  double mass2() const { return Mass2 ; }
  double pdgmass() const { return cdsMass[PID]; }
  double beta() const { return Beta; }
  double tof() const { return TOF; };
  double fl() const { return FlightLength; }
  double dt() const { return Dt; }
  double chi() const { return Chi; }

  TVector3 vertex() { return Vertex     ; }
  TVector3 vcdc()   { return VertexCDC  ; }
  TVector3 vbeam()  { return VertexBeam ; }
  TVector3 momdir() { return MomDir.Unit(); }

  double vdis()   const { return VertexDistance; }
  double vbdis()  const { return VertexBeamDistance; }
  double pbdca()  const { return ProductBeamDCA; }

  double oa()        const { return OpenAngle; }
  double* GetParameters() { return Param; }

  TLorentzVector GetLorentzVector() 
  {
    TLorentzVector lv;
    if(PID>=0)
      lv.SetVectM(momdir()*TMath::Abs(mom()),pdgmass());
    else
      lv.SetVectM(momdir()*TMath::Abs(mom()),mass());
    return lv;
  }

  int cdhseg(const int i) const {return CDHseg[i]; }
  bool ncdh() const { return CDHseg.size(); }
  int ihseg(const int i) const {return IHseg[i]; }
  bool nih() const { return IHseg.size(); }

  ClassDef( pCDS, 1 );
};

class pProduct : public TObject
{
 public:
  pProduct();
  virtual ~pProduct(){};

 private:
  int trackID1;
  int trackID2;

  int CombID;
  int PID;
  double Momentum;
  double Mass;
  double CMomentum;
  double CMass;
  double PDGMass;

  double VertexDistance;
  double VertexBeamDistance;
  double ProductBeamDCA;

  double CVertexDistance;
  double CVertexBeamDistance;
  double CProductBeamDCA;

  double Beta;
  double Gamma;
  double TOF;

  double OpenAngle;
  double AngleLab;
  double AngleCM;
  double COpenAngle;
  double CAngleLab;
  double CAngleCM;

  TVector3 Vertex;
  TVector3 MomDir;
  TVector3 VertexBeam;
  TVector3 CVertex;
  TVector3 CMomDir;
  TVector3 CVertexBeam;

 public:
  void SetTrackID1(const int &id) { trackID1 = id; }
  void SetTrackID2(const int &id) { trackID2 = id; }
  void SetCombID(const int &id) { CombID = id; }
  void SetPID(const int &id) { PID = id; }

  void SetMomentum(const double &p) { Momentum = p; }
  void SetMass(const double &m) { Mass = m; }
  void SetCMomentum(const double &p) { CMomentum = p; }
  void SetCMass(const double &m) { CMass = m; }

  void SetPDGMass(const double &m);
  void SetBeta(const double &b) { Beta = b; }
  void SetGamma(const double &g) { Gamma = g; }
  void SetTOF(const double &t) { TOF = t; }

  void SetVertex(     TVector3 vec)     { Vertex      =vec; }
  void SetVertexBeam( TVector3 vec)     { VertexBeam  =vec; }
  void SetMomDir(     TVector3 vec)     { MomDir      =vec; }
  void SetCVertex(    TVector3 vec)     { CVertex     =vec; }
  void SetCVertexBeam(TVector3 vec)     { CVertexBeam =vec; }
  void SetCMomDir(    TVector3 vec)     { CMomDir     =vec; }
  void SetVDis(     const double &vdis ){ VertexDistance     =vdis; }
  void SetVBDis(    const double &vdis ){ VertexBeamDistance =vdis; }
  void SetPBDCA(    const double &vdis ){ ProductBeamDCA     =vdis; }
  void SetOA(       const double &ang  ){ OpenAngle          =ang; }
  void SetAngleLab( const double &ang  ){ AngleLab           =ang; }
  void SetCVDis(    const double &vdis ){ CVertexDistance    =vdis; }
  void SetCVBDis(   const double &vdis ){ CVertexBeamDistance=vdis; }
  void SetCPBDCA(   const double &vdis ){ CProductBeamDCA    =vdis; }
  void SetCOA(      const double &ang  ){ COpenAngle         =ang; }
  void SetCAngleLab(const double &ang  ){ CAngleLab          =ang; }

  int id1() const { return trackID1; }
  int id2() const { return trackID2; }
  int comb()const { return CombID; }
  int pid() const { return PID; }
  double mom()     const { return Momentum; }
  double mass()    const { return Mass ; }
  double cmom()    const { return CMomentum; }
  double cmass()   const { return CMass ; }
  double pdgmass() const { return PDGMass ; }
  double beta()    const { return Beta; }
  double gamma()   const { return Gamma; }
  double tof()     const { return TOF; };

  double oa()        const { return OpenAngle; }
  double anglelab()  const { return AngleLab; }
  double anglecm()   const { return AngleCM; }
  double coa()       const { return COpenAngle; }
  double canglelab() const { return CAngleLab; }
  double canglecm()  const { return CAngleCM; }

  TVector3 vertex()  { return Vertex; }
  TVector3 vertexb() { return VertexBeam; }
  TVector3 momdir()  { return MomDir; }  
  TVector3 cvertex() { return CVertex; }
  TVector3 cvertexb(){ return CVertexBeam; }
  TVector3 cmomdir() { return CMomDir; }  

  double vdis()   const { return VertexDistance; }
  double vbdis()  const { return VertexBeamDistance; }
  double pbdca()  const { return ProductBeamDCA; }
  double cvdis()  const { return CVertexDistance; }
  double cvbdis() const { return CVertexBeamDistance; }
  double cpbdca() const { return CProductBeamDCA; }

  TLorentzVector GetLorentzVector() 
  {  
    TLorentzVector lv;
    lv.SetVectM(MomDir.Unit()*Momentum,Mass);
    return lv;
  }  
  TLorentzVector GetCLorentzVector() 
  {  
    TLorentzVector lv;
    lv.SetVectM(CMomDir.Unit()*CMomentum,CMass);
    return lv;
  }  
  
  void CalcAngleCM( TLorentzVector beam, const double &mass);

  ClassDef( pProduct, 1 );
};

class Particle : public TObject
{
 public:
  Particle();
  virtual ~Particle(){};

 private:
  typedef std::vector<pBeam> pBeamcontainer;
  pBeamcontainer BeamContainer;

  typedef std::vector<pNC> pNCcontainer;
  pNCcontainer NCContainer;

  typedef std::vector<pPC> pPCcontainer;
  pPCcontainer PCContainer;
  pPCcontainer CVCContainer;

  //  typedef std::vector<pProduct> pProductcontainer;
  //  pProductcontainer ProductContainer;
  
  typedef std::vector<pCDS> pCDScontainer;
  pCDScontainer CDSContainer;
  pCDScontainer ProductContainer;
  
  typedef std::vector<int> TrackIDContainer;
  TrackIDContainer PiplusTrackIDContainer;
  TrackIDContainer PiminusTrackIDContainer;
  TrackIDContainer ProtonTrackIDContainer;
  TrackIDContainer KaonTrackIDContainer;
  TrackIDContainer DeuteronTrackIDContainer;
  TrackIDContainer TritonTrackIDContainer;
  TrackIDContainer Helium3TrackIDContainer;
  TrackIDContainer OtherTrackIDContainer;

  std::map<int,int> CDSTrackIDContainer;
  
  typedef std::vector<HodoscopeLikeHit> hodoContainer;
  hodoContainer NCrawContainer;
  hodoContainer PCrawContainer;
  hodoContainer CVCrawContainer;
  hodoContainer BVCrawContainer;
  hodoContainer CDHrawContainer;
  hodoContainer IHrawContainer;

 public:
  int nProduct(){return ProductContainer.size(); }
  
  int nCDS(){return CDSTrackIDContainer.size(); }
  int nKaon(){return KaonTrackIDContainer.size(); }
  int nPiplus(){return PiplusTrackIDContainer.size(); }
  int nPiminus(){return PiminusTrackIDContainer.size(); }
  int nDeuteron(){return DeuteronTrackIDContainer.size(); }
  int nTriton(){return TritonTrackIDContainer.size(); }
  int nHelium3(){return Helium3TrackIDContainer.size(); }
  int nProton(){return ProtonTrackIDContainer.size(); }
  int nOther(){return OtherTrackIDContainer.size(); }

  int nBeam(){return BeamContainer.size(); }
  int nNC(){return NCContainer.size(); }
  int nPC(){return PCContainer.size(); }
  int nCVC(){return CVCContainer.size(); }

  int nrawNC(){return NCrawContainer.size(); }
  int nrawPC(){return PCrawContainer.size(); }
  int nrawIH(){return IHrawContainer.size(); }
  int nrawCDH(){return CDHrawContainer.size(); }
  int nrawCVC(){return CVCrawContainer.size(); }
  int nrawBVC(){return BVCrawContainer.size(); }

  //  pProduct *product(const int &i){ return &ProductContainer[i]; }
  pCDS *product(const int &i){ return &ProductContainer[i]; }
  pBeam *beam(const int &i){ return &BeamContainer[i]; }
  pNC *nc(const int &i){ return &NCContainer[i]; }
  pPC *pc(const int &i){ return &PCContainer[i]; }
  pPC *cvc(const int &i){ return &CVCContainer[i]; }

  HodoscopeLikeHit *ncraw(const int &i){ return &NCrawContainer[i]; }
  HodoscopeLikeHit *ncrawseg(const int &seg);
  HodoscopeLikeHit *pcraw(const int &i){ return &PCrawContainer[i]; }
  HodoscopeLikeHit *ihraw(const int &i){ return &IHrawContainer[i]; }
  HodoscopeLikeHit *cdhraw(const int &i){ return &CDHrawContainer[i]; }
  HodoscopeLikeHit *cvcraw(const int &i){ return &CVCrawContainer[i]; }
  HodoscopeLikeHit *bvcraw(const int &i){ return &BVCrawContainer[i]; }

  void AddPCraw(HodoscopeLikeHit* hit){ PCrawContainer.push_back(*hit); }
  void AddNCraw(HodoscopeLikeHit* hit){ NCrawContainer.push_back(*hit); }
  void AddIHraw(HodoscopeLikeHit* hit){ IHrawContainer.push_back(*hit); }
  void AddCDHraw(HodoscopeLikeHit* hit){ CDHrawContainer.push_back(*hit); }
  void AddCVCraw(HodoscopeLikeHit* hit){ CVCrawContainer.push_back(*hit); }
  void AddBVCraw(HodoscopeLikeHit* hit){ BVCrawContainer.push_back(*hit); }
  pCDS *cdsi(const int &i){ return &CDSContainer[i]; }
  pCDS *cds(const int &id);
  pCDS *kaon(const int &id);
  pCDS *pip(const int &id);
  pCDS *pim(const int &id);
  pCDS *proton(const int &id);
  pCDS *deuteron(const int &id);
  pCDS *helium3(const int &id);
  pCDS *triton(const int &id);
  pCDS *other(const int &id);
  
  //  void AddProduct( const pProduct &pro);
  void AddProduct( const pCDS &pro);
  void AddCDS( const pCDS &track);
  void AddBeam( const pBeam &beam);
  void AddNC( const pNC &nc);
  void AddPC( const pPC &pc);
  void AddCVC( const pPC &pc);
  
  void Clear();
  void CalcAngleCM(const double &targetmass);
  void CalcNCMom(const TVector3 &vertex);
  //  void CalcPCMom(const TVector3 &vertex);
  void CalcNCMom();

  int keycdc(){ return KEY(nPiplus(),nPiminus(),nKaon(),nProton(),nDeuteron(),nTriton()); }

 private:
  static const unsigned int MASK    = 0x0003;      /* K Mask 2 Bits (0-3) */
  static const int          PIPSHIFT   =  0;
  static const int          PIMSHIFT   =  2;
  static const int          KSHIFT   =  4;
  static const int          PSHIFT   =  6;
  static const int          DSHIFT   =  8;
  static const int          TSHIFT   =  10;
  
  inline int KEY(int npip,int npim,int nk,int np,int nd,int nt){
    return ((((npip)&MASK)) | (((npim)&MASK)<<PIMSHIFT) | (((nk)&MASK)<<KSHIFT) | (((np)&MASK)<<PSHIFT) | (((nd)&MASK)<<DSHIFT)| (((nt)&MASK)<<TSHIFT) );
  }

 public:
  ClassDef( Particle, 1 );

};

#endif
