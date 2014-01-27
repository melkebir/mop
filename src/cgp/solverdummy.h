/* 
 * solverdummy.h
 *
 *  Created on: 20-oct-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVERDUMMY_H_
#define SOLVERDUMMY_H_

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <lemon/core.h>
#include "verbose.h"
#include "molecule.h"
#include "solver.h"

namespace nina {
namespace cgp {

template<typename GR>
class SolverDummy : public Solver<GR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Base class type
  typedef Solver<GR> Parent;

  using Parent::_chargeGroupVector;
  using Parent::_color;
  using Parent::_pMolecule;
  using Parent::updateCorrectedPartialCharges;

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

private:
  typedef std::set<ChargeGroup> ChargeGroupSet;
  typedef typename ChargeGroupSet::const_iterator ChargeGroupSetIt;

public:
  SolverDummy(const Molecule<Graph>* pMolecule)
    : Parent(pMolecule)
  {
  }

  bool solve(int resError)
  {
    const Graph& g = _pMolecule->getGraph();

    ChargeGroupSet cgSet;
    for (NodeIt n(g); n != lemon::INVALID; ++n)
    {
      cgSet.insert(_pMolecule->getSameChargeGroup(n));
    }

    _chargeGroupVector.clear();
    int colIdx = 0;
    for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++, colIdx++)
    {
      _chargeGroupVector.push_back(*cgIt);
      for (ChargeGroupIt nodeIt = cgIt->begin(); nodeIt != cgIt->end(); nodeIt++)
      {
        _color[*nodeIt] = colIdx;
      }
    }

    updateCorrectedPartialCharges(resError);
    return true;
  }
};

} // namespace cgp
} // namespace nina

#endif /* SOLVERDUMMY_H_ */
