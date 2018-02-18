/* ChamberLayer.cpp */
#include "ChamberLayer.h"

ChamberLayer::ChamberLayer()
{
}

ChamberLayer::ChamberLayer( const int &layer, const int &xy ) : Layer(layer), XY(xy)
{
}

void ChamberLayer::Clear()
{
  Hits.clear();
}
