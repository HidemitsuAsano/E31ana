/* ChamberLayer.h */
#ifndef ChamberLayer_h
#define ChamberLayer_h 1

#include "Common.h"

class ChamberLayer : public TObject 
{
 private:
  ChamberLayer();

 public:
  ChamberLayer(const int &layer, const int &xy);
  virtual ~ChamberLayer(){};

 private:
  int Layer;
  int XY;
  std::vector<ChamberLikeHit> Hits;

 public:
  const int GetLayer() const { return Layer; };
  const int GetXY()    const { return XY; };
  const int nHit()     const { return Hits.size(); };
  ChamberLikeHit *GetHit(const int &i){ 
    assert( 0<=i && i<Hits.size() );
    return &Hits[i];
 };

  void AddHit( const ChamberLikeHit &hit ){ Hits.push_back(hit); };
  void DeleteHit( const int &i ){ 
    std::vector<ChamberLikeHit>::iterator it=Hits.begin();
    for( int j=0; j<i; j++ ){
      ++it;
    } 
    Hits.erase(it);
  }
  std::vector<ChamberLikeHit>::iterator DeleteHit( std::vector<ChamberLikeHit>::iterator it ){ 
    it = Hits.erase(it);
    return it;
  }
  void Clear();

  ClassDef( ChamberLayer, 1 )
;};

#endif
