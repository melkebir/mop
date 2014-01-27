/*
 * symmetrization.h
 *
 *  Created on: 21-feb-2013
 *      Author: M. El-Kebir
 */

#ifndef SYMMETRIZATION_H
#define SYMMETRIZATION_H

#include <lemon/core.h>
#include "common/verbose.h"
#include "CGP/molecule.h"

namespace nina {
namespace cgp {

template<typename GR>
class Symmetrization
{
public:
  /// Graph type
  typedef GR Graph;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public:
  /// Weights on nodes
  typedef typename Graph::template NodeMap<double> WeightNodeMap;
  /// Molecule type
  typedef Molecule<Graph> MoleculeType;

protected:
  struct NodeVectorOrder
  {
  private:
    const WeightNodeMap& _partialCharge;

  public:
    NodeVectorOrder(const WeightNodeMap& partialCharge)
      : _partialCharge(partialCharge)
    {
    }

    bool operator() (Node v, Node w)
    {
      return _partialCharge[v] < _partialCharge[w];
    }
  };

  /// Node vector type
  typedef typename std::vector<Node> NodeVector;
  /// Node vector iterator type
  typedef typename NodeVector::const_iterator NodeVectorIt;
  /// Node matrix type
  typedef typename std::vector<NodeVector> NodeMatrix;
  /// Node matrix iterator type
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  /// Node matrix non const iterator type
  typedef typename NodeMatrix::iterator NodeMatrixNonConstIt;
  /// Atom type list type
  typedef typename std::list<int> AtomList;
  /// Atom type list iterator type
  typedef typename AtomList::const_iterator AtomListIt;
  /// Atom type list set type
  typedef typename std::multiset<AtomList> AtomListMultiSet;
  /// Atom type list set iterator type
  typedef typename AtomListMultiSet::const_iterator AtomListSetMultiIt;
  /// Atom list set map type
  typedef typename Graph::template NodeMap<AtomListMultiSet> AtomListMultiSetMap;
  /// Node set type
  typedef typename std::set<Node> NodeSet;
  /// Node set iterator type
  typedef typename NodeSet::const_iterator NodeSetIt;

protected:
  void symmetrizeNonHeavyAtoms(const MoleculeType& molecule,
                               WeightNodeMap& partialCharge);

  const int _neighborhoodSize;
  const double _maxAbsDifference;
  NodeMatrix _eqClasses;

private:
  AtomListMultiSetMap* _pAtomListSetMap;

private:
  void dfs(const MoleculeType& molecule,
           Node node,
           NodeSet visited,
           const AtomList& curPath,
           AtomListMultiSet& pathSet);


public:
  Symmetrization(int neighborhoodSize, double maxAbsDifference)
    : _neighborhoodSize(neighborhoodSize)
    , _maxAbsDifference(maxAbsDifference)
    , _eqClasses()
    , _pAtomListSetMap(NULL)
  {
  }

  virtual ~Symmetrization()
  {
    delete _pAtomListSetMap;
  }

  virtual bool equivalent(const MoleculeType& molecule,
                          Node v,
                          Node w);

  virtual void symmetrize(MoleculeType &molecule,
                          WeightNodeMap& partialCharge);

  virtual void identifyDegree1(MoleculeType& molecule);

  void print(const MoleculeType& molecule,
             std::ostream& out) const;
};
} // namespace cgp
} // namespace nina

template<typename GR>
inline void nina::cgp::Symmetrization<GR>::identifyDegree1(MoleculeType& molecule)
{
  const Graph& g = molecule.getGraph();
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    if (lemon::countIncEdges(g, v) > 1)
    {
      NodeSet adjNodes;
      adjNodes.insert(v);

      for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
      {
        Node w = g.oppositeNode(v, e);
        if (lemon::countIncEdges(g, w) == 1)
          adjNodes.insert(w);
      }

      molecule.setSameChargeGroup(adjNodes);
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cerr << "Node " << molecule.getLabel2(v)
                  << " (" << molecule.getLabel(v)
                  << ") has degree 1 neighbors that have to co-occur" << std::endl;
      }
    }
  }
}

template<typename GR>
inline bool nina::cgp::Symmetrization<GR>::equivalent(const MoleculeType& molecule,
                                                      Node v,
                                                      Node w)
{
  assert(v != lemon::INVALID && w != lemon::INVALID);
  const Graph& g = molecule.getGraph();

  if (_pAtomListSetMap == NULL)
  {
    _pAtomListSetMap = new AtomListMultiSetMap(g);
    for (NodeIt u(g); u != lemon::INVALID; ++u)
    {
      dfs(molecule, u, NodeSet(), AtomList(), (*_pAtomListSetMap)[u]);
    }
  }

  if (v == w)
    return true;
  else
    return (*_pAtomListSetMap)[v] == (*_pAtomListSetMap)[w];
}

template<typename GR>
inline void nina::cgp::Symmetrization<GR>::symmetrize(MoleculeType& molecule,
                                                      WeightNodeMap& partialCharge)
{
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Applying symmetrization, atoms v and w are "
              << "considered chemically equivalent if"
              << std::endl
              << "- |charge(v) - charge(w)| <= " << _maxAbsDifference << "; and"
              << std::endl
              << "- neighborhood up to " << _neighborhoodSize << " atoms match"
              << std::endl;
  }

  _eqClasses.clear();

  const Graph& g = molecule.getGraph();

  // preprocess the molecule
  symmetrizeNonHeavyAtoms(molecule, partialCharge);

  // determine equivalence classes
  IntNodeMap eqIdx(g, -1);
  int nEqClasses = 0;
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    for (NodeIt w = v; w != lemon::INVALID; ++w)
    {
      if (w == v) continue;

      //if (fabs(_partialCharge[v] - _partialCharge[w]) <= maxAbsDifference &&
      if (equivalent(molecule, v, w))
      {
        if (eqIdx[v] != -1)
        {
          eqIdx[w] = eqIdx[v];
        }
        else if (eqIdx[w] != -1)
        {
          eqIdx[v] = eqIdx[w];
        }
        else
        {
          eqIdx[v] = eqIdx[w] = nEqClasses++;
        }
      }
    }
  }

  // put everything in eqClasses
  NodeMatrix eqClasses(nEqClasses);
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    if (eqIdx[v] != -1)
    {
      eqClasses[eqIdx[v]].push_back(v);
    }
  }

  // split up equivalence groups taking partial charges into account
  NodeMatrix newEqClasses;
  for (NodeMatrixNonConstIt nodeVectorIt = eqClasses.begin();
       nodeVectorIt != eqClasses.end(); nodeVectorIt++)
  {
    // sort equivalence class according to partial charges
    std::sort(nodeVectorIt->begin(), nodeVectorIt->end(), NodeVectorOrder(partialCharge));

    // determine splitting points
    NodeVectorIt prevIt = nodeVectorIt->begin();
    for (NodeVectorIt nodeIt = nodeVectorIt->begin();
         nodeIt != nodeVectorIt->end(); nodeIt++)
    {
      // do we have a splitting point?
      if (fabs(partialCharge[*prevIt] - partialCharge[*nodeIt]) > _maxAbsDifference)
      {
        newEqClasses.push_back(NodeVector(prevIt, nodeIt));
        prevIt = nodeIt;
      }
    }

    // add last eq class
    if (prevIt != nodeVectorIt->end())
    {
      newEqClasses.push_back(NodeVector(prevIt, (NodeVectorIt)nodeVectorIt->end()));
    }
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Number of equivalence classes before taking charges into account: "
              << eqClasses.size() << std::endl;
    std::cerr << "Number of equivalence classes after taking charges into account: "
              << newEqClasses.size() << std::endl;
  }

  _eqClasses = newEqClasses;

  molecule.setEqClasses(_eqClasses);

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    print(molecule, std::cerr);
  }

  // symmetrize
  for (NodeMatrixIt nodeVectorIt = newEqClasses.begin();
       nodeVectorIt != newEqClasses.end(); nodeVectorIt++)
  {
    // determine the average
    double sum = 0;
    for (NodeVectorIt nodeIt = nodeVectorIt->begin();
         nodeIt != nodeVectorIt->end(); nodeIt++)
    {
      sum += partialCharge[*nodeIt];
    }

    double avg = sum / nodeVectorIt->size();

    // set everything to the average
    for (NodeVectorIt nodeIt = nodeVectorIt->begin();
         nodeIt != nodeVectorIt->end(); nodeIt++)
    {
      partialCharge[*nodeIt] = avg;
    }
  }
}

template<typename GR>
inline void nina::cgp::Symmetrization<GR>::symmetrizeNonHeavyAtoms(const MoleculeType& molecule,
                                                                   WeightNodeMap& partialCharge)
{
  typedef std::set<int> IntSet;
  typedef std::vector<double> DoubleVector;

  const Graph& g = molecule.getGraph();

  // for each atom v, assign its deg-1 neighbors - that have the same atom type -
  // the same partial charge
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    int d = lemon::countIncEdges(g, v);
    if (d >= 2)
    {
      std::map<int, int> consideredAtomTypes;
      NodeMatrix adjNodes;
      DoubleVector sumCharges;

      // determine neighboring nodes with similar atom types
      for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
      {
        Node w = g.oppositeNode(v, e);

        // ensure that the degree of w is 1
        if (lemon::countIncEdges(g, w) != 1)
          continue;

        // check if the atom type is in adjNodes, if not add it
        int t = molecule.getAtomType(w);
        if (consideredAtomTypes.find(t) == consideredAtomTypes.end())
        {
          consideredAtomTypes[t] = static_cast<int>(adjNodes.size());
          adjNodes.push_back(NodeVector());
          sumCharges.push_back(0);
        }

        adjNodes[consideredAtomTypes[t]].push_back(w);
        sumCharges[consideredAtomTypes[t]] += partialCharge[w];
      }

      // assign neighboring nodes with similar atom types the same charge
      for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
      {
        Node w = g.oppositeNode(v, e);

        if (lemon::countIncEdges(g, w) != 1)
          continue;

        int idx = consideredAtomTypes[molecule.getAtomType(w)];
        partialCharge[w] = sumCharges[idx] / adjNodes[idx].size();
      }
    }
  }
}

template<typename GR>
inline void nina::cgp::Symmetrization<GR>::dfs(const MoleculeType& molecule,
                                               Node node,
                                               NodeSet visited,
                                               const AtomList& curPath,
                                               AtomListMultiSet& pathSet)
{
  assert(curPath.size() == 0 || pathSet.find(curPath) != pathSet.end());
  assert(visited.find(node) == visited.end());

  const Graph& g = molecule.getGraph();

  // 0. visit node
  visited.insert(node);

  // 1. add node to curPath and put it in pathSet
  AtomList newPath = curPath;
  newPath.push_back(molecule.getAtomType(node));
  pathSet.insert(newPath);

  // 2a. stop if the maximal path size is reached
  if (static_cast<int>(newPath.size()) == _neighborhoodSize)
    return;

  // 2b. recurse on unvisited neighbors otherwise
  for (IncEdgeIt e(g, node); e != lemon::INVALID; ++e)
  {
    Node childNode = g.oppositeNode(node, e);
    if (visited.find(childNode) == visited.end())
    {
      dfs(molecule, childNode, visited, newPath, pathSet);
    }
  }
}

template<typename GR>
inline void nina::cgp::Symmetrization<GR>::print(const MoleculeType& molecule,
                                                 std::ostream& out) const
{
  for (NodeMatrixIt it = _eqClasses.begin(); it != _eqClasses.end(); it++)
  {
    bool first = true;
    const NodeVector& nodes = *it;

    for (NodeVectorIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        out << ", ";
      }
      out << molecule.getLabel(*nodeIt) << " ("
          << molecule.getLabel2(*nodeIt) << ")";
    }
    out << std::endl;
  }
}

#endif // SYMMETRIZATION_H
