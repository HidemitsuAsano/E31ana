// MyAnalysisAlloc.h
// 2014.12.22
// T. Yamaga

#ifndef MyAnalysisAlloc_h
#define MyhAnalysisAlloc_h

#include <string>

#include "MyAnalysisBase.h"

class MyAnalysisAlloc
{
  public:
    explicit MyAnalysisAlloc();
    ~MyAnalysisAlloc();

    MyAnalysisBase* MyAnalysisAllocator();

};

#endif
