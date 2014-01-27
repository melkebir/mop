/*
 * treedecomposition.h
 *
 *  Created on: 13-feb-2013
 *      Author: M. El-Kebir
 */

#ifndef TREEDECOMPOSITION_H
#define TREEDECOMPOSITION_H

#include <ostream>
#include <vector>
#include <set>
#include <list>
#include <limits>
#include <string>
#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/bfs.h>

namespace nina
{

template<typename GR>
class TreeDecomposition
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Tree decomposition criteria
  typedef enum {
    TREE_DECOMPOSITION_GREEDY_DEGREE,
    TREE_DECOMPOSITION_GREEDY_FILL_IN,
  } TreeDecompositionCriterion;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  /// Labels on nodes type
  typedef typename Graph::template NodeMap<std::string> LabelNodeMap;
  /// Orienter type
  typedef typename lemon::Orienter<const Graph, BoolEdgeMap> RootedTree;
  /// Rooted arc iterator type
  typedef typename RootedTree::ArcIt RootedArcIt;
  /// Rooted out arc iterator type
  typedef typename RootedTree::OutArcIt RootedOutArcIt;
  /// Rooted out arc iterator type
  typedef typename RootedTree::InArcIt RootedInArcIt;
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

protected:
  /// Node set type
  typedef std::set<Node> NodeSet;
  /// Node set const iterator type
  typedef typename NodeSet::const_iterator NodeSetIt;
  /// Node list type
  typedef std::list<Node> NodeList;
  /// Node list iterator type
  typedef typename NodeList::const_iterator NodeListIt;
  /// Node set map type
  typedef typename Graph::template NodeMap<NodeSet> NodeSetMap;
  /// Node map type
  typedef typename Graph::template NodeMap<Node> NodeMap;
  /// Arc lookup on graph type
  typedef typename lemon::DynArcLookUp<Graph> ArcLookUp;
  /// Subgraph type
  typedef lemon::FilterNodes<Graph, BoolNodeMap> SubGraph;
  /// Subgraph node iterator type
  typedef typename SubGraph::NodeIt SubNodeIt;
  /// Subgraph incident edge iterator type
  typedef typename SubGraph::IncEdgeIt SubIncEdgeIt;
  /// Arc lookup on subgraph type
  typedef typename lemon::DynArcLookUp<SubGraph> SubArcLookUp;

public:
  /// Constructor
  TreeDecomposition(const Graph& g);

  /// Destructor
  virtual ~TreeDecomposition() {}

  /// Computes a tree decomposition
  virtual void init(TreeDecompositionCriterion criterion);

  /// Print the tree decomposition
  void print(std::ostream& out) const;

  /// Return treewidth
  int getTreewidth() const { return _treewidth; }

  /// Return number of bags
  int getNumberOfBags() const { return _nBags; }

  /// Return tree
  const Graph& getTree() { return _tree; }

  /// Return root bag
  Node getRootBag() const { return _rootBag; }

  /// Return nodes in bag
  const NodeSet& getNodesInBag(Node bag) const
  {
    return _bag[bag];
  }

  /// Return bag size
  int getBagSize(Node bag) const
  {
    assert(bag != lemon::INVALID);
    return static_cast<int>(_bag[bag].size());
  }

  /// Return maximum ancestral bag size
  int getMaxAncestralBagSize(Node bag) const
  {
    assert(bag != lemon::INVALID);
    return _maxAncestralBagSize[bag];
  }

  /// Return number of ancestors
  int getNumberOfAncestors(Node bag) const
  {
    assert(bag != lemon::INVALID);
    return _numberOfAncestors[bag];
  }

  /// Return nodes in the subtree rooted at bag
  const NodeSet& getNodesSubTree(Node bag) const
  {
    assert(bag != lemon::INVALID);
    return _nodesSubTree[bag];
  }

  /// Return number of nodes in the subtree rooted at bag
  int getNodesSubTreeSize(Node bag) const
  {
    assert(bag != lemon::INVALID);
    return static_cast<int>(_nodesSubTree[bag].size());
  }

  /// Return rooted tree decomposition
  const RootedTree& getRootedTree() { return _rootedTree; }

  /// Return max BFS level
  int getMaxBfsLevel() const { return _maxBfsLevel; }

  /// Return bags by BFS level
  const NodeVector& getBagsByBfsLevel(int level) const
  {
    assert(0 <= level && level < static_cast<int>(_bagsByBfsLevel.size()));
    return _bagsByBfsLevel[level];
  }

  /// Root the tree decomposition with rootBag as the root
  /// Initializes _rootedG, _orientation, _bfsLevel, _maxBfsLevel, _bagsByBfsLevel
  void root(Node rootBag = lemon::INVALID);

  /// Print tree decomposition
  void printTreeDecomposition(std::ostream& out) const;

  /// Print rooted tree decomposition
  void printRootedTreeDecomposition(std::ostream& out) const;


private:
  /// Initializes _maxAncestralBagSize, _numberOfAncestors
  /// and _nodesSubTree of the sub tree rooted at bag
  void initRootedTree(Node bag);

protected:
  /// Original graph
  const Graph& _g;
  /// Tree decomposition
  Graph _tree;
  /// Bags of the tree decomposotion
  NodeSetMap _bag;
  /// Treewidth
  int _treewidth;
  /// Number of bags
  int _nBags;

  /// Used for rooting the molecule graph
  BoolEdgeMap _orientation;
  /// The rooted tree decomposition
  RootedTree _rootedTree;
  /// The root bag
  Node _rootBag;

  /// Nodes in the bags in the sub tree rooted at the given bag
  NodeSetMap _nodesSubTree;
  /// Number of ancestors of a specfied bag
  IntNodeMap _numberOfAncestors;
  /// Max ancestral bag size
  IntNodeMap _maxAncestralBagSize;

  /// Maximum BFS level
  int _maxBfsLevel;
  /// BFS level of a specified bag
  IntNodeMap _bfsLevel;
  /// Holds all bags with similar specified BFS level
  NodeMatrix _bagsByBfsLevel;

  /// Computes tree decomposition from an elimination ordering
  void transformPermutation(const NodeVector& perm,
                            const IntNodeMap& permMap);

  /// Makes a copy of the internal graph
  void makeCopy(Graph& h,
                NodeMap& crossRef,
                const IntNodeMap* pPermMapG = NULL,
                IntNodeMap* pPermMapH = NULL) const;

  /// Get depth if the tree were to be rooted at bag
  int getDepth(Node bag, Node parentBag = lemon::INVALID) const;

  double getAvgChildDepth(Node bag) const;

private:
  /// Computes tree decomposition from an elimination ordering
  void transformPermutation(const NodeVector& perm,
                            const IntNodeMap& permMap,
                            const NodeMap& crossRef,
                            Graph& g);
  /// Eliminate node
  int eliminate(Graph& g,
                const ArcLookUp& arcLookUp,
                Node v,
                bool shaduw = false);

  /// Choose the next node to add to the elimination order
  Node choose(Graph& h,
              const ArcLookUp& arcLookUp,
              TreeDecompositionCriterion criterion);
};
} // namespace nina

template<typename GR>
inline double nina::TreeDecomposition<GR>::getAvgChildDepth(Node bag) const
{
  typedef std::vector<int> IntVector;
  typedef typename IntVector::const_iterator IntVectorIt;

  std::vector<int> degVector;
  for (IncEdgeIt e(_tree, bag); e != lemon::INVALID; ++e)
  {
    Node childBag = _tree.oppositeNode(bag, e);
    degVector.push_back(getDepth(childBag));
  }

  if (degVector.size() == 0)
    return 0;
  else if (degVector.size() == 1)
    return degVector.front();

  int n = 0;
  double absDiffSum = 0;
  for (IntVectorIt degIt1 = degVector.begin(); degIt1 != degVector.end(); degIt1++)
  {
    for (IntVectorIt degIt2 = degIt1 + 1; degIt2 != degVector.end(); degIt2++)
    {
      absDiffSum += abs(*degIt1 - *degIt2);
      n++;
    }
  }

  return absDiffSum / n;
}

template<typename GR>
inline int nina::TreeDecomposition<GR>::getDepth(Node bag, Node parentBag) const
{
  int maxDepth = 0;
  for (IncEdgeIt e(_tree, bag); e != lemon::INVALID; ++e)
  {
    Node childBag = _tree.oppositeNode(bag, e);
    if (childBag != parentBag)
      maxDepth = std::max(maxDepth, 1 + getDepth(childBag, bag));
  }

  return maxDepth;
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::initRootedTree(Node bag)
{
  RootedInArcIt a(_rootedTree, bag);
  if (a != lemon::INVALID)
  {
    Node parentBag = _rootedTree.source(a);
    _numberOfAncestors[bag] = 1 + _numberOfAncestors[parentBag];
    _maxAncestralBagSize[bag] = std::max(_maxAncestralBagSize[parentBag],
      static_cast<int>(_bag[parentBag].size()));
  }

  // recurse on the children
  _nodesSubTree[bag] = _bag[bag];
  for (RootedOutArcIt a(_rootedTree, bag); a != lemon::INVALID; ++a)
  {
    Node childBag = _rootedTree.target(a);
    initRootedTree(childBag);
    _nodesSubTree[bag].insert(_nodesSubTree[childBag].begin(), _nodesSubTree[childBag].end());
  }
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::root(Node rootBag)
{
  if (rootBag == lemon::INVALID)
  {
    // determine the root such that the tree is most balanced
    double avgDeg = std::numeric_limits<double>::max();
    for (NodeIt bag(_tree); bag != lemon::INVALID; ++bag)
    {
      double curAvgDeg = getAvgChildDepth(bag);
      if (curAvgDeg < avgDeg)
      {
        avgDeg = curAvgDeg;
        _rootBag = bag;
      }
    }
  }
  else
  {
    _rootBag = rootBag;
  }

  // do a stupid BFS here and store the BFS levels
  lemon::Bfs<Graph> bfs(_tree);
  bfs.distMap(_bfsLevel);
  bfs.run(_rootBag);

  // determine max BFS level
  _maxBfsLevel = lemon::mapMaxValue(_tree, _bfsLevel);

  // add nodes to _bagsByBfsLevel
  _bagsByBfsLevel = NodeMatrix(_maxBfsLevel+1);
  for (NodeIt n(_tree); n != lemon::INVALID; ++n)
  {
    _bagsByBfsLevel[_bfsLevel[n]].push_back(n);
  }

  // orient the digraph correctly
  for (RootedArcIt a(_rootedTree); a != lemon::INVALID; ++a)
  {
    if (_rootedTree.source(a) != bfs.predNode(_rootedTree.target(a)))
    {
      _rootedTree.reverseArc(a);
    }
  }

  // update _maxAncestralBagSize, _numberOfAncestors and _nodesSubTree
  _numberOfAncestors[_rootBag] = 0;
  _maxAncestralBagSize[_rootBag] = 0;
  initRootedTree(_rootBag);
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::init(TreeDecompositionCriterion criterion)
{
  Graph h;
  NodeMap crossRef(h);

  makeCopy(h, crossRef);
  ArcLookUp arcLookUp(h);

  NodeVector perm;
  IntNodeMap permMap(_g);

  int n = lemon::countNodes(_g);
  for (int i = 0; i < n; i++)
  {
    Node v = choose(h, arcLookUp, criterion);

    permMap[crossRef[v]] = static_cast<int>(perm.size());
    perm.push_back(crossRef[v]);

    eliminate(h, arcLookUp, v);
  }

  transformPermutation(perm, permMap);
}

template<typename GR>
inline nina::TreeDecomposition<GR>::TreeDecomposition(const Graph& g)
  : _g(g)
  , _tree()
  , _bag(_tree)
  , _treewidth(-1)
  , _nBags(0)
  , _orientation(_tree)
  , _rootedTree(_tree, _orientation)
  , _rootBag(lemon::INVALID)
  , _nodesSubTree(_tree)
  , _numberOfAncestors(_tree)
  , _maxAncestralBagSize(_tree)
  , _maxBfsLevel(-1)
  , _bfsLevel(_tree)
  , _bagsByBfsLevel()
{
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::makeCopy(Graph& h,
                                                  NodeMap& crossRef,
                                                  const IntNodeMap* pPermMapG,
                                                  IntNodeMap* pPermMapH) const
{
  if (pPermMapG && pPermMapH)
  {
    lemon::graphCopy(_g, h)
      .nodeCrossRef(crossRef)
      .nodeMap(*pPermMapG, *pPermMapH)
      .run();
  }
  else
  {
    lemon::graphCopy(_g, h)
      .nodeCrossRef(crossRef)
      .run();
  }
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::transformPermutation(const NodeVector& perm,
                                                              const IntNodeMap& permMap)
{
  Graph cpyG;

  NodeMap crossRef(cpyG);
  IntNodeMap cpyPermMap(cpyG);

  makeCopy(cpyG, crossRef, &permMap, &cpyPermMap);

  // re-create perm
  NodeVector cpyPerm(perm.size(), lemon::INVALID);
  for (NodeIt v(cpyG); v != lemon::INVALID; ++v)
  {
    Node orgV = crossRef[v];
    cpyPerm[permMap[orgV]] = v;
  }

  transformPermutation(cpyPerm, cpyPermMap, crossRef, cpyG);
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::transformPermutation(const NodeVector& perm,
                                                              const IntNodeMap& permMap,
                                                              const NodeMap& crossRef,
                                                              Graph& g)
{
  ArcLookUp arcLookUp(g);
  NodeMap repBag(_g, lemon::INVALID);

  _treewidth = -1;
  _nBags = 0;

  NodeMap lowestNeighborMap(_tree);

  const int n = static_cast<int>(perm.size());
  for (int idx = 0; idx < n; idx++)
  {
    const Node curNode = perm[idx];
    const Node orgCurNode = crossRef[curNode];

    if (idx == n - 1)
    {
      Node newBag = _tree.addNode();
      _bag[newBag].insert(orgCurNode);
      repBag[orgCurNode] = newBag;

      if (_treewidth < 0) _treewidth = 0;
      _nBags++;
    }
    else if (idx == n - 2)
    {
      Node newBag = _tree.addNode();
      _bag[newBag].insert(orgCurNode);
      _bag[newBag].insert(crossRef[perm[idx+1]]);
      repBag[orgCurNode] = newBag;
      repBag[crossRef[perm[idx+1]]] = newBag;

      if (_treewidth < 1) _treewidth = 1;
      _nBags++;
      break;
    }
    else
    {
      Node newBag = _tree.addNode();
      _bag[newBag].insert(orgCurNode);
      repBag[orgCurNode] = newBag;

      // compute contents of newBag
      // and determine lowest neighbor of curNode
      int lowestNeighborIdx = std::numeric_limits<int>::max();
      lowestNeighborMap[newBag] = lemon::INVALID;
      for (IncEdgeIt e(g, curNode); e != lemon::INVALID; ++e)
      {
        Node neighbor = g.oppositeNode(curNode, e);

        if (permMap[neighbor] < lowestNeighborIdx)
        {
          lowestNeighborMap[newBag] = crossRef[neighbor];
          lowestNeighborIdx = permMap[neighbor];
        }

        _bag[newBag].insert(crossRef[neighbor]);
      }

      // eliminate curNode
      eliminate(g, arcLookUp, curNode);

      if (_treewidth < static_cast<int>(_bag[newBag].size()) - 1)
      {
        _treewidth = static_cast<int>(_bag[newBag].size()) - 1;
      }
      _nBags++;
    }
  }

  // now let's add the edges
  for (int idx = 0; idx < (n >= 2 ? n - 2 : n - 1); idx++)
  {
    Node orgCurNode = crossRef[perm[idx]];
    Node curBag = repBag[orgCurNode];
    _tree.addEdge(curBag, repBag[lowestNeighborMap[curBag]]);
  }
}

template<typename GR>
inline int nina::TreeDecomposition<GR>::eliminate(Graph& g,
                                                  const ArcLookUp& arcLookUp,
                                                  Node v,
                                                  bool shadow)
{
  assert(v != lemon::INVALID);
  int res = 0;

  // determine nodes incident to v
  NodeList incNodes;
  for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
  {
    incNodes.push_back(g.oppositeNode(v, e));
  }

  // add edges between nodes (u,w) adjacent to v that are not adjacent themselves
  for (NodeListIt nodeIt1 = incNodes.begin(); nodeIt1 != incNodes.end(); nodeIt1++)
  {
    for (NodeListIt nodeIt2 = nodeIt1; nodeIt2 != incNodes.end(); nodeIt2++)
    {
      if (nodeIt1 == nodeIt2)
        continue;

      if (arcLookUp(*nodeIt1, *nodeIt2) == lemon::INVALID)
      {
        res++;

        if (!shadow)
          g.addEdge(*nodeIt1, *nodeIt2);
      }
    }
  }

  // finally remove v
  if (!shadow)
    g.erase(v);

  return res;
}

template<typename GR>
inline void nina::TreeDecomposition<GR>::print(std::ostream& out) const
{
  out << "# " << lemon::countNodes(_tree) << std::endl;
  for (NodeIt x(_tree); x != lemon::INVALID; ++x)
  {
    out << _tree.id(x);

    const NodeSet& bag = _bag[x];
    for (NodeSetIt vIt = bag.begin(); vIt != bag.end(); ++vIt)
    {
      out << "\t" << _g.id(*vIt);
    }

    out << std::endl;
  }

  out << "# " << lemon::countEdges(_tree) << std::endl;
  for (EdgeIt e(_tree); e != lemon::INVALID; ++e)
  {
    out << _tree.id(_tree.u(e)) << "\t" << _tree.id(_tree.v(e)) << std::endl;
  }
}

template<typename GR>
inline typename nina::TreeDecomposition<GR>::Node
nina::TreeDecomposition<GR>::choose(Graph& h,
                                    const ArcLookUp& arcLookUp,
                                    TreeDecompositionCriterion criterion)
{
  switch (criterion)
  {
  case TREE_DECOMPOSITION_GREEDY_DEGREE:
    {
      // lowest degree
      Node minNode = lemon::INVALID;

      int minDegree = std::numeric_limits<int>::max();
      for (NodeIt v(h); v != lemon::INVALID; ++v)
      {
        int d = 0;
        for (IncEdgeIt e(h, v); e != lemon::INVALID; ++e) d++;

        if (d < minDegree)
        {
          minNode = v;
          minDegree = d;
        }
      }

      return minNode;
    }
  case TREE_DECOMPOSITION_GREEDY_FILL_IN:
    {
      // lowest number of non-adjacent neigbors
      Node minNode = lemon::INVALID;

      int minFilledEdges = std::numeric_limits<int>::max();
      for (NodeIt v(h); v != lemon::INVALID; ++v)
      {
        int filledEdges = eliminate(h, arcLookUp, v, true);

        if (filledEdges < minFilledEdges)
        {
          minNode = v;
          minFilledEdges = filledEdges;
        }
        else if (filledEdges == 0)
        {
          break;
        }
      }

      return minNode;
    }
  }

  return lemon::INVALID;
}

template<typename GR>
void nina::TreeDecomposition<GR>::printRootedTreeDecomposition(std::ostream& out) const
{
  out << "digraph G {" << std::endl;

  // print the nodes
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    out << "\t" << _tree.id(v) << " [label=\"";

    const NodeSet& bag = _bag[v];
    bool first = true;
    for (NodeSetIt nodeIt = bag.begin(); nodeIt != bag.end(); nodeIt++)
    {
      if (!first)
        out << ", ";
      else
        first = false;

      out << _g.id(*nodeIt);
    }

    out << "\\n" << _maxAncestralBagSize[v]
      << "\\n" << _numberOfAncestors[v]
      << "\\n" << getMaxBlacklistSize(v) << "\\n";

    const NodeSet& subBags = _nodesSubTree[v];
    first = true;
    for (NodeSetIt nodeIt = subBags.begin(); nodeIt != subBags.end(); nodeIt++)
    {
      if (!first)
        out << ", ";
      else
        first = false;

      out << _g.id(*nodeIt);
    }

    out << "\"]" << std::endl;
  }

  // print the edges
  for (RootedArcIt a(_rootedTree); a != lemon::INVALID; ++a)
  {
    Node u = _rootedTree.source(a);
    Node v = _rootedTree.target(a);

    out << "\t" << _tree.id[u]
      << " -> "
      << _tree.id(v) << std::endl;
  }

  out << "}" << std::endl;
}

template<typename GR>
void nina::TreeDecomposition<GR>::printTreeDecomposition(std::ostream& out) const
{
  out << "graph G {" << std::endl;

  // print the nodes
  for (NodeIt v(_tree); v != lemon::INVALID; ++v)
  {
    out << "\t" << _tree.id(v) << " [label=\"";
    const NodeSet& bag = _bag[v];
    bool first = true;
    for (NodeSetIt nodeIt = bag.begin(); nodeIt != bag.end(); nodeIt++)
    {
      if (!first)
        out << ", ";
      else
        first = false;

      out << _g.id(*nodeIt);
    }
    out << "\"]" << std::endl;
  }

  // print the edges
  for (EdgeIt e(_tree); e != lemon::INVALID; ++e)
  {
    Node u = _tree.u(e);
    Node v = _tree.v(e);

    out << "\t" << _tree.id(u)
      << " -- "
      << _tree(v) << std::endl;
  }

  out << "}" << std::endl;
}

#endif // TREEDECOMPOSITION_H
