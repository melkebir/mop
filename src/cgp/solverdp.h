/* 
 * solverdp.h
 *
 *  Created on: 27-sep-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_DP_H_
#define SOLVER_DP_H_

#include "solver.h"
#include "molecule.h"
#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>
#include <iostream>
#include <algorithm>

template<typename GR>
class SolverDP : public Solver<GR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Base class type
  typedef Solver<GR> Parent;

  using Parent::_chargeGroupVector;
  using Parent::_color;
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
  /// Charge group matrix type
  typedef typename std::vector<ChargeGroupVector> ChargeGroupMatrix;
  /// Charge group matrix iterator
  typedef typename ChargeGroupMatrix::const_iterator ChargeGroupMatrixIt;
  /// Charge group set type
  typedef typename ExtMolecule<Graph>::ChargeGroupSet ChargeGroupSet;
  /// Charge group set iterator type
  typedef typename ExtMolecule<Graph>::ChargeGroupSetIt ChargeGroupSetIt;

protected:
  /// Orienter type
  typedef typename lemon::Orienter<const Graph, BoolEdgeMap> RootedGraph;
  /// Rooted arc iterator type
  typedef typename RootedGraph::ArcIt RootedArcIt;
  /// Rooted out arc iterator type
  typedef typename RootedGraph::OutArcIt RootedOutArcIt;
  /// Charge group matrix map type
  typedef typename Graph::template NodeMap<ChargeGroupMatrix> ChargeGroupMatrixMap;
  /// Charge group map type
  typedef typename Graph::template NodeMap<ChargeGroup> ChargeGroupMap;
  /// Node vector type
  typedef typename std::vector<Node> NodeVector;
  /// Node vector iterator type
  typedef typename NodeVector::const_iterator NodeVectorIt;
  /// Node vector map type
  typedef typename Graph::template NodeMap<NodeVector> NodeVectorMap;
  /// Node matrix type
  typedef typename std::vector<NodeVector> NodeMatrix;
  /// Node matrix type iterator
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  /// Int vector type
  typedef typename std::vector<int> IntVector;

private:
  /// The molecule
  const ExtMolecule<Graph>* _pExtMolecule;
  /// The unrooted graph
  const Graph& _unrootedG;
  /// Maximum BFS level
  int _maxBfsLevel;
  /// BFS level of a specified node
  IntNodeMap _bfsLevel;
  /// Maximal charge group size at a specified node
  IntNodeMap _maxCgSize;

  /// Holds the optimal (sub)solution
  DoubleNodeMap _opt;
  /// Holds the charge group belonging to the optimal solution
  ChargeGroupMap _optCg; 

  /// Map of children of a node
  NodeVectorMap _children;

  /// Holds all nodes with similar specified BFS level
  NodeMatrix _nodeMatrix;

  /// Holds all charge groups whose representative is the node in question
  ChargeGroupMatrixMap _cgMatrixMap;

  /// Used for rooting the molecule graph
  BoolEdgeMap _orientation;
  /// The rooted graph
  RootedGraph _rootedG;

  /// Root the tree
  /// Initializes _rootedG, _orientation, _bfsLevel, _maxBfsLevel, _children, _nodeMatrix
  void root();
  /// Initializes _cgMatrixMap, _opt and _optCg
  void init();

  /// Computes objective function given a node n,
  /// which is a representative of a charge group cg
  /// \pre: all nodes higher BFS level have been computed
  double computeObj(const ChargeGroup& cg, const Node n) const;

  /// Reconstructs the solution from _optCg
  void reconstruct(const ChargeGroup& cg, const Node n, const int color);

  /// returns false when there is no next
  static bool next(const IntVector& limits, 
                   const int k,
                   IntVector& current);

  /// returns false when there is no next
  static bool next(const IntVector& limits, 
                   IntVector& current);

public:
  SolverDP(const ExtMolecule<Graph>* pExtMolecule)
    : Parent(pExtMolecule)
    , _pExtMolecule(pExtMolecule)
    , _unrootedG(pExtMolecule->getGraph())
    , _maxBfsLevel(-1)
    , _bfsLevel(_unrootedG)
    , _maxCgSize(_unrootedG)
    , _opt(_unrootedG)
    , _optCg(_unrootedG)
    , _children(_unrootedG)
    , _nodeMatrix()
    , _cgMatrixMap(_unrootedG)
    , _orientation(_unrootedG)
    , _rootedG(_unrootedG, _orientation)
  {
  }

  bool solve(int resError);

  void printRootedTree(std::ostream& out) const;
};

template<typename GR>
void SolverDP<GR>::init()
{
  // for every node n there must be k entries in _cgMatrixMap[n]
  for (NodeIt n(_unrootedG); n != lemon::INVALID; ++n)
  {
    const size_t k = _pExtMolecule->getMaxChargeGroupSize(n);

    _cgMatrixMap[n] = ChargeGroupMatrix(k);
    _maxCgSize[n] = 0;

    // initialize DP table
    _opt[n] = std::numeric_limits<double>::max();
    _optCg[n] = ChargeGroup();
  }

  // now go through all charge groups
  const size_t overallK = _pExtMolecule->getMaxOverallChargeGroupSize();
  for (size_t i = 0; i < overallK; i++)
  {
    const ChargeGroupSet& cgSet = _pExtMolecule->getChargeGroup(i+1);
    for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++)
    {
      // determine the representative
      Node repNode = lemon::INVALID;
      int repLevel = _maxBfsLevel + 1;
      for (ChargeGroupIt nodeIt = cgIt->begin(); nodeIt != cgIt->end(); nodeIt++)
      {
        if (_bfsLevel[*nodeIt] < repLevel)
        {
          repNode = *nodeIt;
          repLevel = _bfsLevel[*nodeIt];
        }
      }

      // add charge group to repNode at the correct slot (based on the cardinality)
      _cgMatrixMap[repNode][i].push_back(*cgIt);
      if (static_cast<int>(cgIt->size()) > _maxCgSize[repNode])
        _maxCgSize[repNode] = cgIt->size();
    }
  }
}

template<typename GR>
void SolverDP<GR>::root()
{
  // do a stupid BFS here and store the BFS levels
  lemon::Bfs<Graph> bfs(_unrootedG);
  bfs.distMap(_bfsLevel);
  bfs.run(NodeIt(_unrootedG));

  // determine max BFS level
  _maxBfsLevel = lemon::mapMaxValue(_unrootedG, _bfsLevel);

  // add nodes to _nodeMatrix
  _nodeMatrix = NodeMatrix(_maxBfsLevel+1);
  for (NodeIt n(_unrootedG); n != lemon::INVALID; ++n)
  {
    _nodeMatrix[_bfsLevel[n]].push_back(n);
  }

  // orient the digraph correctly
  for (RootedArcIt a(_rootedG); a != lemon::INVALID; ++a)
  {
    if (_rootedG.source(a) != bfs.predNode(_rootedG.target(a)))
    {
      _rootedG.reverseArc(a);
    }
  }

  // now add the children
  for (NodeIt n(_unrootedG); n != lemon::INVALID; ++n)
  {
    for (RootedOutArcIt a(_rootedG, n); a != lemon::INVALID; ++a)
      _children[n].push_back(_rootedG.target(a));
  }
}

template<typename GR>
bool SolverDP<GR>::solve(int resError)
{
  //assert(lemon::tree(g));
  if (!lemon::tree(_unrootedG))
  {
    std::cerr << "Input graph is not a tree!" << std::endl;
    exit(1);
  }

  root();
  init();

  if (g_verbosity >= VERBOSE_DEBUG)
  {
    //printRootedTree(std::cerr);

    // shows the charge groups stored at every node in the tree
    for (NodeIt n(_unrootedG); n != lemon::INVALID; ++n)
    {
      std::cout << "Node " << _pExtMolecule->getLabel(n) << std::endl;
      const size_t k = _pExtMolecule->getMaxChargeGroupSize(n);
      for (size_t i = 0; i < k; i++)
      {
        for (ChargeGroupVectorIt cgIt = _cgMatrixMap[n][i].begin(); 
            cgIt != _cgMatrixMap[n][i].end(); cgIt++)
        {
          _pExtMolecule->printChargeGroup(*cgIt, std::cout);
        }
      }
      std::cout << std::endl;
    }
  }

  // this is complicated stuff...
  // time to solve the DP, 
  // start from nodes with the highest BFS levels and work towards the root
  for (int i = _maxBfsLevel; i >= 0; i--)
  {
    const NodeVector& nodes = _nodeMatrix[i];
    for (NodeVectorIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
    {
      const Node n = *nodeIt;
      const size_t k = _pExtMolecule->getMaxChargeGroupSize(n);
      const NodeVector& children = _children[n];
      const int nChildren = static_cast<int>(children.size());

      // now we have to consider all combinations starting from card = (0, 0, ..., 0)
      // where every element corresponds to a child of node n
      // these combinations are about cg size, and must respect two things:
      // 1) the sum of the elements is at most k-1
      // 2) every element with index i must be at most _maxCgSize[i-th child of n]
      IntVector cardLimits(nChildren, 0);

      for (int j = 0; j < nChildren; j++)
      {
        // at most k-1 because node n is going to be part of the charge group 
        cardLimits[j] = std::min(static_cast<int>(k-1), _maxCgSize[children[j]]);
      }

      IntVector card(nChildren, 0);
      do
      {
        // here we have to consider all cg combinations 
        // where every cg i has size card[i]
        // we start by computing the limits
        IntVector cgIdxLimits(nChildren, 0);
        IntVector idx(nChildren, 0);
        for (int j = 0; j < nChildren; j++)
        {
          // two cases...
          if (card[j] == 0)
          {
            // we pick nothing from child j, we flag this using -1
            cgIdxLimits[j] = -1;
            idx[j] = -1;
          }
          else
          {
            // we pick exactly card[j] > 0 nodes from child j, 
            // then we simply set the max idx 
            // by looking up how many cg-s there are with that size
            cgIdxLimits[j] = 
              static_cast<int>(_cgMatrixMap[children[j]][card[j]-1].size()-1);
          }
        }

        do
        {
          // build charge group as specified by idx and card
          ChargeGroup cg;
          cg.insert(n);

          for (int j = 0; j < nChildren; j++)
          {
            Node child_j = children[j];
            const int idx_j = idx[j];
            const int card_j = card[j];

            if (idx_j != -1)
            {
              assert(card_j != 0);
              const ChargeGroupMatrix& cgMatrix = _cgMatrixMap[child_j];
              if (cgMatrix[card_j-1].size() != 0)
              {
                const ChargeGroup& cg_j = cgMatrix[card_j-1][idx_j];
                cg.insert(cg_j.begin(), cg_j.end());
              }
              else
              {
                std::cerr << "STRANGE" << std::endl;
              }
            }
          }
          assert(cg.size() <= k);

          // now evaluate whether cg is better than what we currently have
          // we need to do a DFS starting from node n
          double cg_obj = _pExtMolecule->computeError(cg) + computeObj(cg, n);
          if (cg_obj < _opt[n])
          {
            _opt[n] = cg_obj;
            _optCg[n] = cg;
          }
        } while (next(cgIdxLimits, idx));

      } while (next(cardLimits, k-1, card));
    }
  }

  const Node root = _nodeMatrix[0].front();
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Obj value: " << _opt[root] << std::endl;
  }

  _chargeGroupVector.clear();
  _chargeGroupVector.push_back(_optCg[root]);
  reconstruct(_optCg[root], root, 0);
  updateCorrectedPartialCharges(resError);

  //IntVector limits, current(3, 0);
  //limits.push_back(2);
  //limits.push_back(2);
  //limits.push_back(0);

  //do
  //{
  //  for (int i = 0; i < 3; i++)
  //    std::cout << current[i] << ", ";
  //  std::cout << std::endl;
  //} while (next(limits, 4, current));
  return true;
}

template<typename GR>
void SolverDP<GR>::printRootedTree(std::ostream& out) const
{
  out << "digraph G {" << std::endl;
  
  // nodes per level
  for (int i = 0; i <= _maxBfsLevel; i++)
  {
    out << "\t{" << std::endl
      << "\t\trank=same" << std::endl;

    for (NodeIt n(_unrootedG); n != lemon::INVALID; ++n)
    {
      if (_bfsLevel[n] == i)
      {
        out << "\t\t" << _unrootedG.id(n) 
          << " [label=\"" << _pExtMolecule->getLabel(n) << "\"]" << std::endl;
      }
    }

    out << "\t}" << std::endl;
  }

  // arcs
  for (RootedArcIt a(_rootedG); a != lemon::INVALID; ++a)
  {
    out << "\t" << _unrootedG.id(_rootedG.source(a)) << " -> "
      << _unrootedG.id(_rootedG.target(a)) << std::endl;
  }

  out << "}" << std::endl;
}

template<typename GR>
double SolverDP<GR>::computeObj(const ChargeGroup& cg,
                                const Node n) const
{
  if (cg.find(n) == cg.end())
  {
    // we're done
    return _opt[n];
  }
  else
  {
    // recurse over the children of n
    double obj = 0;

    const NodeVector& children = _children[n];
    for (NodeVectorIt childIt = children.begin(); childIt != children.end(); childIt++)
    {
      obj += computeObj(cg, *childIt);
    }

    return obj;
  }
}

template<typename GR>
void SolverDP<GR>::reconstruct(const ChargeGroup& cg, 
                               const Node n,
                               const int color)
{
  const NodeVector& children = _children[n];
  //_color[n] = static_cast<int>(_chargeGroupVector.size() - 1);
  _color[n] = color;

  int i = 1;
  for (NodeVectorIt childIt = children.begin(); childIt != children.end(); childIt++, i++)
  {
    Node child = *childIt;
    if (cg.find(child) == cg.end())
    {
      _chargeGroupVector.push_back(_optCg[child]);
      reconstruct(_optCg[child], child, color + i);
    }
    else
    {
      reconstruct(cg, child, color);
    }
  }
}

template<typename GR>
bool SolverDP<GR>::next(const IntVector& limits,
                        IntVector& current)
{
  assert(current.size() == limits.size());

  const int n = static_cast<int>(current.size());
  IntVector cpyCurrent = current;

  while (true)
  {
    bool update = false;
    for (int i = n - 1; i >= 0; i--)
    {
      if (cpyCurrent[i] == limits[i])
        continue;
      
      update = true;
      cpyCurrent[i]++;
      // and reset everything to the right
      for (int j = i+1; j < n; j++)
      {
        if (limits[j] != -1)
          cpyCurrent[j] = 0;
      }

      break;
    }

    // nothing new could be found
    if (!update)
      return false;

    // we're probably good, 
    // we just have to check whether the limit's been exceeded
    bool all_bigger_or_equal = true;
    bool one_strictly_bigger = false;
    for (int i = 0; i < n; i++)
    {
      all_bigger_or_equal &= (cpyCurrent[i] >= limits[i]);
      one_strictly_bigger |= cpyCurrent[i] > limits[i];
    }

    if (all_bigger_or_equal && one_strictly_bigger)
      return false;
    else
    {
      current = cpyCurrent;
      return true;
    }
  }
}

template<typename GR>
bool SolverDP<GR>::next(const IntVector& limits,
                        const int k,
                        IntVector& current)
{
  assert(current.size() == limits.size());

  const int n = static_cast<int>(current.size());
  IntVector cpyCurrent = current;

  while (true)
  {
    bool update = false;
    for (int i = n - 1; i >= 0; i--)
    {
      if (cpyCurrent[i] == limits[i])
        continue;
      
      update = true;
      cpyCurrent[i]++;
      // and reset everything to the right
      for (int j = i+1; j < n; j++)
        cpyCurrent[j] = 0;

      break;
    }

    // nothing new could be found
    if (!update)
      return false;

    // we're probably good, 
    // we just have to check whether the limit's been exceeded
    // and also the sum of course
    bool all_bigger_or_equal = true;
    bool one_strictly_bigger = false;
    int sum = 0;
    for (int i = 0; i < n; i++)
    {
      all_bigger_or_equal &= (cpyCurrent[i] >= limits[i]);
      one_strictly_bigger |= cpyCurrent[i] > limits[i];
      sum += cpyCurrent[i];
    }

    if (all_bigger_or_equal && one_strictly_bigger)
      return false;
    else if (sum <= k)
    {
      current = cpyCurrent;
      return true;
    }
  }
}

#endif /* SOLVER_DP_H_ */
