#ifndef MKTREE_H
#define MKTREE_H 1

#include "analysis.h"

#ifdef __cplusplus
extern "C" {
#endif
  int  uana( int evnum, s_raw *raw );
  void uintr( int isig );
  int  uinit();
  int  ufini();
  int  uopt( int i );
  int usca( int nsca, unsigned int *sca );

  int uscaler( int i );
  int uhead( int i, int j );
  s_run info_run();
  unsigned long uevent_id();

#ifdef __cplusplus
}
#endif

#endif
