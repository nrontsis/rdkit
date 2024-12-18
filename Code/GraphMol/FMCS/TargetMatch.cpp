#include "Seed.h"
#include "SubstructMatchCustom.h"
#include "Target.h"
#include "TargetMatch.h"

namespace RDKit {
namespace FMCS {

TargetMatch::TargetMatch() {}

TargetMatch::TargetMatch(const TargetMatch& src) { *this = src; }

TargetMatch& TargetMatch::operator=(const TargetMatch& src) {
  Empty = src.Empty;
  if (!Empty) {
    MatchedAtomSize = src.MatchedAtomSize;
    MatchedBondSize = src.MatchedBondSize;
    TargetAtomIdx = src.TargetAtomIdx;
    TargetBondIdx = src.TargetBondIdx;
    VisitedTargetBonds = src.VisitedTargetBonds;
    VisitedTargetAtoms = src.VisitedTargetAtoms;
  }
  return *this;
}
bool TargetMatch::empty() const { return Empty; }

void TargetMatch::clear() {
  Empty = true;

  TargetAtomIdx.clear();
  TargetBondIdx.clear();
  VisitedTargetBonds.clear();
  VisitedTargetAtoms.clear();
}

void TargetMatch::init(const Seed& seed, const match_V_t& match,
                       const ROMol& query, const Target& target) {
  TargetAtomIdx.clear();
  TargetAtomIdx.resize(query.getNumAtoms(), NotSet);
  TargetBondIdx.clear();
  TargetBondIdx.resize(query.getNumBonds(), NotSet);
  VisitedTargetBonds.resize(target.Molecule->getNumBonds());
  VisitedTargetAtoms.resize(target.Molecule->getNumAtoms());
  VisitedTargetBonds.reset();
  VisitedTargetAtoms.reset();

  MatchedAtomSize = match.size();
  for (const auto& m : match) {
    TargetAtomIdx[seed.MoleculeFragment.Atoms.at(m.first)->getIdx()] =
        m.second;
    VisitedTargetAtoms.set(m.second);
  }

  MatchedBondSize = 0;
  for (const auto bond : seed.MoleculeFragment.Bonds) {
    unsigned int i = bond->getBeginAtomIdx();
    unsigned int j = bond->getEndAtomIdx();
    unsigned int ti = TargetAtomIdx.at(i);
    unsigned int tj = TargetAtomIdx.at(j);
    const auto tb = target.Molecule->getBondBetweenAtoms(ti, tj);
    if (tb) {
      ++MatchedBondSize;
      TargetBondIdx[bond->getIdx()] = tb->getIdx();  // add
      VisitedTargetBonds.set(tb->getIdx());
    }
  }
  Empty = false;
}
}  // namespace FMCS
}  // namespace RDKit