#ifndef BLEVENT_H
#define BLEVENT_H 1

#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "TKO.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "EventAlloc.h"
#include "EventTemp.h"

class BLEvent
{
 private:
/*   std::string ConfFileName;// = "conf/analyzer.conf"; */
  ConfMan *confMan;
  std::string RootFileName;// = "tmp.root";

 public:
  BLEvent(){};
  ~BLEvent(){};

 public:
  int BMain( int argc, char **argv, ConfMan *con );
};

#include "analysis.h"

#ifdef __cplusplus
extern "C" {
#endif
  int cmain( int, char** );
  int  uana( int evnum, s_raw *raw );
  void uintr( int isig );
  int  uinit();
  int  ufini();
  int utime( int i );
  int  uopt( int i );
  int usca( int nsca, unsigned int *sca );
  int set_num_modules( int nscach, int nsmp );
  int set_crate_type( int cr, int sl, int type);
  int uscaler( int i );
  int uhead( int i, int j );
  s_run info_run();
  unsigned long uevent_id();

#ifdef __cplusplus
}
#endif

#endif
