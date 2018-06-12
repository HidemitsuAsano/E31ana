#ifndef SlewingData_h
#define SlewingData_h 1

#include "HodoscopeLikeHit.h"

class SlewingData : public TObject
{
 public:
  SlewingData();
  SlewingData(const int cid1, const int cid2);
  virtual ~SlewingData(){};

 private:
  int BeamPID;

  int CID1;
  int Seg1;
  double Eu1;
  double Ed1;
  double Emean1;
  double Tu1;
  double Td1;
  double Ctu1;
  double Ctd1;
  double Ctmean1;

  int CID2;
  int Seg2;
  double Eu2;
  double Ed2;
  double Emean2;
  double Tu2;
  double Td2;
  double Ctu2;
  double Ctd2;
  double Ctmean2;

 public:
  int beam_pid() const { return BeamPID; };

  int cid1() const { return CID1; };
  int seg1() const { return Seg1; };
  double eu1() const { return Eu1; };
  double ed1() const { return Ed1; };
  double emean1() const { return Emean1; };
  double tu1() const { return Tu1; };
  double td1() const { return Td1; };
  double ctu1() const { return Ctu1; };
  double ctd1() const { return Ctd1; };
  double ctmean1() const { return Ctmean1; };

  int cid2() const { return CID2; };
  int seg2() const { return Seg2; };
  double eu2() const { return Eu2; };
  double ed2() const { return Ed2; };
  double emean2() const { return Emean2; };
  double tu2() const { return Tu2; };
  double td2() const { return Td2; };
  double ctu2() const { return Ctu2; };
  double ctd2() const { return Ctd2; };
  double ctmean2() const { return Ctmean2; };

  void SetBeamPID(const double &pid){ BeamPID=pid; };
  void Set(HodoscopeLikeHit *hit1, HodoscopeLikeHit *hit2);
  void Calc(const double par[4], const double offset[4]);
  void Calc(const double par[4], const double par1[4], const double offset[4]);
  void SetOffset(const int &id, const double &offset);
  void clear();

  ClassDef(SlewingData, 1);
};
#endif
