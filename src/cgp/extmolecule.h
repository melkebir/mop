/*
 * extmolecule.h
 *
 *  Created on: 21-sep-2011
 *      Author: M. El-Kebir
 */

#ifndef EXTMOLECULE_H_
#define EXTMOLECULE_H_

#include "molecule.h"
#include "verbose.h"
#include <set>
#include <vector>
#include <lemon/core.h>
#include <iostream>

template<typename GR>
class ExtMolecule : public Molecule<GR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Base class type
  typedef Molecule<GR> Parent;

  using Parent::_g;
  using Parent::_d;
  using Parent::_k;
  using Parent::generate;
  using Parent::isFeasible;
  using Parent::getLabel;
  using Parent::computeError;
  using Parent::printChargeGroup;
  using Parent::getNumberOfAtoms;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public:
  /// Charge group type
  typedef typename Parent::ChargeGroup ChargeGroup;
  /// Charge group iterator
  typedef typename Parent::ChargeGroupIt ChargeGroupIt;
  /// Charge group vector type
  typedef typename Parent::ChargeGroupVector ChargeGroupVector;
  /// Charge group vector iterator
  typedef typename Parent::ChargeGroupVectorIt ChargeGroupVectorIt;
  /// Charge group set type
  typedef typename Parent::ChargeGroupSet ChargeGroupSet;
  /// Charge group set iterator
  typedef typename Parent::ChargeGroupSetIt ChargeGroupSetIt;
  /// Charge group non const set iterator
  typedef typename Parent::ChargeGroupSetNonConstIt ChargeGroupSetNonConstIt;
  /// Vector of charge group sets type
  typedef typename Parent::ChargeGroupSetVector ChargeGroupSetVector;
  /// Vector of charge group sets iterator type
  typedef typename Parent::ChargeGroupSetVectorIt ChargeGroupSetVectorIt;
  /// Vector of charge group sets iterator type
  typedef typename Parent::ChargeGroupSetVectorNonConstIt ChargeGroupSetVectorNonConstIt;

protected:
  ChargeGroupSetVector _chargeGroups;

  void compute(size_t i);

public:
  ExtMolecule()
    : Parent()
    , _chargeGroups()
  {
  }

  void printChargeGroups(std::ostream& out) const;

  size_t getNumberOfChargeGroups() const;

  const ChargeGroupSet& getChargeGroup(size_t cardinality) const
  {
    assert(1 <= cardinality && cardinality <= _chargeGroups.size());
    return _chargeGroups[cardinality-1];
  }

  virtual void initLocalK(double diameter);

  virtual size_t getMaxOverallChargeGroupSize() const
  {
    return _chargeGroups.size();
  }
};

template<typename GR>
inline void ExtMolecule<GR>::initLocalK(double diameter)
{
  assert(diameter >= 0);
  _d = diameter;

  _chargeGroups.clear();
  _chargeGroups.push_back(ChargeGroupSet());

  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    // { v }
    ChargeGroup cgInit;
    cgInit.insert(v);

    // add { { v } } to cgSetVec[0]
    _chargeGroups[0].insert(cgInit);
  }

  // generate all charge groups
  while (generate(_chargeGroups));

  // determine k for each node
  size_t curSize = 1;
  for (ChargeGroupSetVectorNonConstIt cgSetIt = _chargeGroups.begin();
       cgSetIt != _chargeGroups.end(); cgSetIt++, curSize++)
  {
    ChargeGroupSet& cgSet = *cgSetIt;

    for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end();)
    {
      ChargeGroupSetIt nextCgIt = cgIt;
      nextCgIt++;

      const ChargeGroup& cg = *cgIt;
      if (isFeasible(cg))
      {
        for (ChargeGroupIt nodeIt = cg.begin(); nodeIt != cg.end(); nodeIt++)
        {
          assert(cg.size() == curSize);
          _k[*nodeIt] = curSize;
        }
      }
      else
      {
        // remove it, if not feasible
        cgSet.erase(cgIt);
      }

      cgIt = nextCgIt;
    }
  }

  if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
  {
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
      std::cerr << getLabel(v) << "\t" << _k[v] << std::endl;
  }
}

template<typename GR>
size_t ExtMolecule<GR>::getNumberOfChargeGroups() const
{
  size_t res = 0;

  for (size_t i = 0; i < _chargeGroups.size(); i++)
  {
    res += _chargeGroups[i].size();
  }

  return res;
}

template<typename GR>
void ExtMolecule<GR>::printChargeGroups(std::ostream& out) const
{
  for (size_t i = 0; i < _chargeGroups.size(); i++)
  {
    out << "Charge groups of size: " << i+1 << std::endl;

    const ChargeGroupSet& cgSet = _chargeGroups[i];
    for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++)
    {
      printChargeGroup(*cgIt, out);
    }

    out << std::endl;
  }
}

#endif /* EXTMOLECULE_H_ */
