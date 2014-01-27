/* 
 * solvermatching.h
 *
 *  Created on: 14-sep-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_MATCHING_H_
#define SOLVER_MATCHING_H_

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <vector>
#include <lemon/core.h>
#include <lemon/matching.h>
#include "solver.h"
#include "verbose.h"
#include "molecule.h"

template<typename GR>
class SolverMatching : public Solver<GR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Base class type
  typedef Solver<GR> Parent;

  using Parent::_pMolecule;
  using Parent::_chargeGroupVector;
  using Parent::_color;
  using Parent::updateCorrectedPartialCharges;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  /// Weights on edges
  typedef typename Graph::template EdgeMap<double> WeightEdgeMap;
  /// Node map to Node
  typedef typename Graph::template NodeMap<Node> NodeToNodeMap;
  /// Edge map to Edge
  typedef typename Graph::template EdgeMap<Edge> EdgeToEdgeMap;
  /// Maximum Weight Matching type
  typedef typename lemon::MaxWeightedPerfectMatching<Graph, WeightEdgeMap> MWM;

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
  /// The graph
  Graph _g;
  /// Map to original nodes
  NodeToNodeMap _origNode;
  /// Map from original nodes
  NodeToNodeMap _node;
  /// Map node to its counterpart
  NodeToNodeMap _counterpart;
  /// Map to indicate whether a node is real/fake
  BoolNodeMap _isFake;
  /// The weights
  WeightEdgeMap _weight;
  /// Max weight matching
  MWM _mwm;

  void init();

public:
  SolverMatching(const Molecule<Graph>* pMolecule)
    : Parent(pMolecule)
    , _g()
    , _origNode(_g)
    , _node(_pMolecule->getGraph())
    , _counterpart(_g)
    , _isFake(_g)
    , _weight(_g)
    , _mwm(_g, _weight)
  {
  };

  void solve(int resError);
};

template<typename GR>
inline void SolverMatching<GR>::init()
{
  // make a copy of the molecule
  graphCopy(_pMolecule->getGraph(), _g)
    .nodeRef(_node)
    .nodeCrossRef(_origNode)
    .run();

  if (g_verbosity >= VERBOSE_DEBUG)
    std::cout << "Original edges: " << std::endl;

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    Node orig_u = _origNode[_g.u(e)];
    Node orig_v = _origNode[_g.v(e)];

    // set the edge weights: i.e. for e=(u,v): w(e) = |w(u) + w(v)|
    _weight[e] = -fabs(_pMolecule->getPartialCharge(orig_u) 
        + _pMolecule->getPartialCharge(orig_v));

    if (g_verbosity >= VERBOSE_DEBUG)
      std::cout << "(" << _g.id(_g.u(e)) << "," 
        << _g.id(_g.v(e)) << ")\t" << _weight[e] << std::endl;
  }

  if (g_verbosity >= VERBOSE_DEBUG)
    std::cout << "Fake edges: " << std::endl;

  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    Node fake_n = _g.addNode();
    Edge fake_e = _g.addEdge(n, fake_n);

    _counterpart[fake_n] = n;
    _counterpart[n] = fake_n;
    _isFake[n] = false;
    _isFake[fake_n] = true;

    _weight[fake_e] = -fabs(_pMolecule->getPartialCharge(_origNode[n]));

    if (g_verbosity >= VERBOSE_DEBUG)
      std::cout << "(" << _g.id(n) << "," 
        << _g.id(fake_n) << ")\t" << _weight[fake_e] << std::endl;
  }

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    Node u = _g.u(e);
    Node v = _g.v(e);
    if (!_isFake[u] && !_isFake[v])
    {
      Node fake_u = _counterpart[u];
      Node fake_v = _counterpart[v];
      Edge fake_e = _g.addEdge(fake_u, fake_v);
      _weight[fake_e] = 0;

      if (g_verbosity >= VERBOSE_DEBUG)
        std::cout << "(" << _g.id(fake_u) << "," 
          << _g.id(fake_v) << ")\t" << _weight[fake_e] << std::endl;
    }
  }
}

template<typename GR>
inline void SolverMatching<GR>::solve(int resError)
{
  init();
  _chargeGroupVector.clear();
  _mwm.init();
  _mwm.start();

  if (g_verbosity >= VERBOSE_ESSENTIAL)
    std::cerr << "Matching weight: " << _mwm.matchingWeight() << std::endl;

  // obtain the charge groups from the solution
  int color = 0;
  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    if (!_mwm.matching(e))
      continue;

    Node u = _g.u(e);
    Node v = _g.v(e);

    // if both u and v are fake, simply skip it, e has no meaning then
    if (_isFake[u] && _isFake[v])
      continue;

    ChargeGroup chargeGroup;
    if (!_isFake[u] && !_isFake[v])
    {
      chargeGroup.insert(_origNode[u]);
      chargeGroup.insert(_origNode[v]);
      _color[_origNode[u]] = _color[_origNode[v]] = color++;
    }
    else if (!_isFake[u] && _isFake[v])
    {
      chargeGroup.insert(_origNode[u]);
      _color[_origNode[u]] = color++;
    }
    else if (_isFake[u] && !_isFake[v])
    {
      chargeGroup.insert(_origNode[v]);
      _color[_origNode[v]] = color++;
    }

    _chargeGroupVector.push_back(chargeGroup);
  }

  updateCorrectedPartialCharges(resError);
}

#endif /* SOLVER_MATCHING_H_ */
