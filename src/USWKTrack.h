#ifndef USWKTrack_h
#define USWKTrack_h 1

#include "GlobalVariables.h"
#include <TVector3.h>

class USWKTrack : public TObject
{
 public:
  USWKTrack();
  USWKTrack(const USWKTrack &right);
  virtual ~USWKTrack(){};

 private:
  std::string ParticleName;
  double Mass;
  double Charged;
  double InTime;
  TVector3 InPosition;
  TVector3 InMomentum;

  std::vector<double> Time;
  std::vector<TVector3> Position;
  std::vector<TVector3> Momentum;

  double OutTime;
  TVector3 OutPosition;
  TVector3 OutMomentum;

 public:
  // Setter
  void SetParticleName(std::string &name);

  void SetInTime(const double &time){ InTime=time; };
  void SetInPosition(const TVector3 &pos){ InPosition=pos; };
  void SetInMomentum(const TVector3 &mom){ InMomentum=mom; };

  void SetTime(const double &time){ Time.push_back(time); };
  void SetPosition(const TVector3 &pos){ Position.push_back(pos); };
  void SetMomentum(const TVector3 &mom){ Momentum.push_back(mom); };

  void SetOutTime(const double &time){ OutTime=time; };
  void SetOutPosition(const TVector3 &pos){ OutPosition=pos; };
  void SetOutMomentum(const TVector3 &mom){ OutMomentum=mom; };

  // Getter
  std::string particlename() const { return ParticleName; };
  double mass() const { return Mass; };
  double charge() const { return Charged; };

  double intime() const { return InTime; };
  TVector3 inposition() const { return InPosition; };
  TVector3 inmomentum() const { return InMomentum; };

  int NPoint() const { return Time.size(); };
  double time(const int &i) const { return Time[i]; };
  TVector3 position(const int &i) const { return Position[i]; };
  TVector3 momentum(const int &i) const { return Momentum[i]; };

  double time() const { return Time[Time.size()-1]; };
  TVector3 position() const { return Position[Position.size()-1]; };
  TVector3 momentum() const { return Momentum[Momentum.size()-1]; };

  double outtime() const { return OutTime; };
  TVector3 outposition() const { return OutPosition; };
  TVector3 outmomentum() const { return OutMomentum; };

  // Methode
  void Clear();

  ClassDef(USWKTrack, 1);
};
#endif
