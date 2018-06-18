#ifndef ScalerMapMan_h
#define ScalerMapMan_h 1

#include <vector>
#include <string>
#include <stdlib.h>

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class ScalerMapMan : public TObject
{
 public:
  ScalerMapMan();
  ScalerMapMan( const std::string &name );
  virtual ~ScalerMapMan() {};

  void SetFileName( const std::string &filename );
  bool Initialize();

 private:

  std::string FileName;
  int NumModules;

  typedef std::vector <std::string> ScalerMapContainer;
  ScalerMapContainer ScalerContainer;

 public:
  int GetNumModules() { return NumModules; }
  std::string GetFileName() { return FileName; }
  std::string GetName( const int &i );

  void PrintMap();
  void Clear();
  
  ClassDef( ScalerMapMan, 1 );
};

#endif

    
