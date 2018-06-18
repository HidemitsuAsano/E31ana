#include "EventAnalysisSca.h"

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisSca *event = new EventAnalysisSca();
  return (EventTemp*)event;
}
