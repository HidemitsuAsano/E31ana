// CDHCluster.h
//

#ifndef CDHCLUSTER_h
#define CDHCLUSTER_h

#include <vector>
#include <map>
#include "TVector3.h"
#include "HodoscopeLikeHit.h"
#include "ConfMan.h"

class CDHCluster
{
	public:
		CDHCluster() {};
		~CDHCluster() {};

	private:
		typedef std::vector<HodoscopeLikeHit> HodoscopeLikeContainer;
		HodoscopeLikeContainer CDHContainer;

	public:
		int nCDH() const { return CDHContainer.size(); }
		HodoscopeLikeHit* CDH( const int &i ) { return &CDHContainer[i]; }
		void Calc();
		void Clear();
		void AddHit(HodoscopeLikeHit* hit){ CDHContainer.push_back(*hit); }

	private:
		double EMean;
		double CTMean;
	public:
		double emean() const { return EMean; }
		double ctmean() const { return CTMean; }
		void SetEMean( const double &in ){ EMean=in; }
		void SetCTMean( const double &in ){ CTMean=in; }

	private:
		TVector3 HitPos;
		double HitPosition;
	public:
		TVector3 pos() const { return HitPos; }
		double x() const { return HitPos.X(); }
		double y() const { return HitPos.Y(); }
		double z() const { return HitPos.Z(); }
		double hitpos() const { return HitPosition; }
		void SetPos( const TVector3 &pos ) { HitPos=pos; }
		void SetHitPosition( const double &in ) { HitPosition=in; }

};

#endif

