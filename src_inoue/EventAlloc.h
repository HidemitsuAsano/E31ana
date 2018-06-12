// EventAlloc

#ifndef EventAlloc_h
#define EventAlloc_h 1

#include <string>

#include "EventTemp.h"

class EventAlloc
{
  // parivate:
  
  //static EventAlloc *evAlloc;

 public:
  explicit EventAlloc();
  ~EventAlloc();

  EventTemp *EventAllocator();

};

#endif

