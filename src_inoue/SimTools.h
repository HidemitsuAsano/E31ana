#ifndef SIMTOOLS_HH
#define SIMTOOLS_HH 1

#include "KnuclRootData.h"

namespace SimTools
{
  DetectorHit *GetHit(const int &cid, const int &seg, DetectorData *data);
  DetectorHit *GetHit(const int &cid, const int &layer, const int &wire, DetectorData *data);

  Track *GetTrack(const int &id, MCData *data);
  int Generation(const int &id, MCData *data);
  std::vector<Track*> GetDaughterTracks(const int &id, MCData *data);
  std::vector<Track*> GetParentTracks(const int &id, MCData *data);

  std::string ParticleName(const int &pdg);
  std::string ReactionName(const int &rid);

  static TString ElementSymbol[]={ "", "H", "He", "Li", "Be", "B", "C", "N", "F", "Ne",
				   "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
				   "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
				   "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", };

};

#endif
