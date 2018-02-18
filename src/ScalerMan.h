#ifndef ScalerMan_h
#define ScalerMan_h 1

#include <vector>
#include <iostream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class Scaler : public TObject
{
 public:
  Scaler();
  virtual ~Scaler() {};

 private:
  std::string Name;
  int Value;

 public:
  std::string name() const { return Name; }
  int val()          const { return Value; }

  void SetName( const std::string &name ) { Name = name; }
  void SetValue(const int &i ) { Value = i; }

  void Clear();

  ClassDef(Scaler,1);
};

class ScalerMan : public TObject
{
 public:
  ScalerMan();
  virtual ~ScalerMan() {};

 private:
  int BlockEventNumber;

  typedef std::vector <Scaler> ScalerContainer;
  ScalerContainer ScaContainer;

 public:
  int blev() const { return BlockEventNumber; }
  int nsca() const { return ScaContainer.size(); }
  Scaler *sca( const int &i ) { return &ScaContainer[i]; }

  void SetBlockEventNumber( const int &evnum ) { BlockEventNumber = evnum; }
  void AddHit( const int &val, const std::string &name );

  void Clear();

  ClassDef( ScalerMan, 1 );
};

#endif
