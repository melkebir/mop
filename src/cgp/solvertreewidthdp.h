/*
 * solvertreewidthdp.h
 *
 *  Created on: 25-oct-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVERTREEWIDTHDP_H_
#define SOLVERTREEWIDTHDP_H_

#include "solver.h"
#include "treedecomposition.h"
#include "molecule.h"
#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <set>

template<typename GR>
class SolverTreewidthDP : public Solver<GR>
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

public:
  /// Charge group type
  typedef typename Parent::ChargeGroup ChargeGroup;
  /// Charge group iterator
  typedef typename Parent::ChargeGroupIt ChargeGroupIt;
  /// Charge group vector type
  typedef typename Parent::ChargeGroupVector ChargeGroupVector;
  /// Charge group vector iterator type
  typedef typename Parent::ChargeGroupVectorIt ChargeGroupVectorIt;
  /// Charge group matrix type
  typedef typename std::vector<ChargeGroupVector> ChargeGroupMatrix;
  /// Charge group matrix iterator type
  typedef typename ChargeGroupMatrix::const_iterator ChargeGroupMatrixIt;
  /// Molecule type
  typedef Molecule<Graph> MoleculeType;
  /// String to node map type
  typedef typename MoleculeType::StringToNodeMap StringToNodeMap;
  /// String to node map iterator
  typedef typename MoleculeType::StringToNodeMapIt StringToNodeMapIt;
  /// Labels on nodes type
  typedef typename MoleculeType::LabelNodeMap LabelNodeMap;
  /// Tree decomposition type
  typedef TreeDecomposition<Graph> TreeDecompositionType;

protected:
  /// Orienter type
  typedef typename TreeDecompositionType::RootedTree RootedTree;
  /// Rooted arc iterator type
  typedef typename TreeDecompositionType::RootedArcIt RootedArcIt;
  /// Rooted out arc iterator type
  typedef typename TreeDecompositionType::RootedOutArcIt RootedOutArcIt;
  /// Rooted in arc iterator type
  typedef typename TreeDecompositionType::RootedInArcIt RootedInArcIt;
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
  /// Node set type
  typedef typename std::set<Node> NodeSet;
  /// Node set iterator type
  typedef typename NodeSet::const_iterator NodeSetIt;
  /// Node set map type
  typedef typename Graph::template NodeMap<NodeSet> NodeSetMap;
  /// Node set set type
  typedef typename std::set<NodeSet> NodeSetSet;
  /// Node set set type iterator
  typedef typename NodeSetSet::const_iterator NodeSetSetIt;
  /// Node set matrix type
  typedef typename std::vector<NodeSetSet> NodeSetSetMatrix;
  /// Node set matrix type iterator
  typedef typename NodeSetSetMatrix::const_iterator NodeSetSetMatrixIt;
  /// Node set vector type
  typedef typename std::vector<NodeSet> NodeSetVector;
  /// Node set vector type iterator
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
  /// Node set matrix type
  typedef typename std::vector<NodeSetVector> NodeSetMatrix;
  /// Node set matrix type iterator
  typedef typename NodeSetMatrix::const_iterator NodeSetMatrixIt;
  /// Node set together with index pair type
  typedef typename std::pair<NodeSet, int> NodeSetIntPair;
  /// Partition type
  typedef typename std::set<NodeSet> Partition;
  /// Partition type iterator
  typedef typename Partition::const_iterator PartitionIt;
  /// Partition vector type
  typedef typename std::vector<Partition> PartitionVector;
  /// Partition vector type iterator
  typedef typename PartitionVector::const_iterator PartitionVectorIt;
  /// Partition matrix map type
  typedef typename Graph::template NodeMap<PartitionVector> PartitionVectorMap;

  struct NodeSetIntPairComp
  {
    bool operator() (const NodeSetIntPair& lhs, const NodeSetIntPair& rhs) const
    {
      return lhs.first < rhs.first;
    }
  };

  /// Set of node sets type
  typedef typename std::set<NodeSetIntPair, NodeSetIntPairComp> NodeSetIntPairSet;
  /// Set of node sets type iterator
  typedef typename NodeSetIntPairSet::const_iterator NodeSetIntPairSetIt;
  /// Matrix of sets of node sets type
  typedef typename std::vector<NodeSetIntPairSet> NodeSetIntPairSetMatrix;
  /// Matrix of sets of node sets type iterator
  typedef typename NodeSetIntPairSetMatrix::const_iterator NodeSetIntPairSetMatrixIt;
  /// Matrix of sets of node sets map type
  typedef typename Graph::template NodeMap<NodeSetIntPairSetMatrix> NodeSetIntPairSetMatrixMap;
  /// Int vector type
  typedef typename std::vector<int> IntVector;
  /// Int vector map type
  typedef typename Graph::template NodeMap<IntVector> IntVectorMap;
  /// Int vector type
  typedef typename std::vector<double> DoubleVector;
  /// Int vector map type
  typedef typename Graph::template NodeMap<DoubleVector> DoubleVectorMap;
  /// Size_t vector type
  typedef typename std::vector<size_t> SizeTVector;

  struct TracebackEntry
  {
    Node _bag;
    int _blacklistSize;
    int _blacklistIdx;

    TracebackEntry()
      : _bag(lemon::INVALID)
      , _blacklistSize(-1)
      , _blacklistIdx(-1)
    {
    }

    TracebackEntry(const Node bag, int blacklistSize, int blacklistIdx)
      : _bag(bag)
      , _blacklistSize(blacklistSize)
      , _blacklistIdx(blacklistIdx)
    {
    }
  };

  /// Traceback vector type
  typedef typename std::vector<TracebackEntry> TracebackVector;
  /// Traceback vector type iterator
  typedef typename TracebackVector::const_iterator TracebackVectorIt;
  /// Traceback matrix type
  typedef typename std::vector<TracebackVector> TracebackMatrix;
  /// Traceback map type
  typedef typename Graph::template NodeMap<TracebackMatrix> TracebackMatrixMap;

private:
  /// The original graph
  const Graph& _g;
  /// Tree decomposition
  TreeDecompositionType _treeDecomposition;

  /// All possible blacklists of a given bag
  NodeSetIntPairSetMatrixMap _blacklistMatrix;

  /// Holds optimal partition given a bag and a blacklist
  PartitionVectorMap _optPartition;
  /// Holds optimal score given a bag and a blacklist
  DoubleVectorMap _optCost;
  /// Traceback map
  TracebackMatrixMap _backTrace;

  /// Initializes tree decomposition
  /// Initializes _cgMatrixMap, _optPartition and _optPartitionCg
  void init();

  void solve(Node bag, const NodeSetIntPair& blacklistPair);

  /// Returns the next part having same cardinality, returns false if there is none
  static bool next(const NodeVector& nodes, SizeTVector& currentPart);

  void makeAndProcessPartition(const Node bag,
                               const NodeSet& blacklist,
                               const int blacklistIdx,
                               NodeVector availableNodes, 
                               const Partition& currentPartition);

  void extendPartition(const Node bag,
                       const NodeSet& blacklist,
                       const int blacklistIdx,
                       NodeSet allowedNodes,
                       Partition partitionToExtend, 
                       Partition extendedPartition);

  double compute(const Node bag,
                 const NodeSet& blacklist,
                 const int blacklistIdx,
                 const Partition& extendedPartition);

  void findNonEmptyBags(const Node bag,
                        const NodeSet& toRemove,
                        BoolNodeMap& emptyFlag,
                        NodeVector& nonEmptyBags);

  void reconstruct(const Node bag,
                   const int blacklistSize,
                   const int blacklistIdx,
                   int& colors);

  double computeCost(const Partition& extendedPartition) const
  {
    double cost = 0;
    for (PartitionIt partIt = extendedPartition.begin(); 
      partIt != extendedPartition.end(); partIt++)
    {
      cost += _pMolecule->computeError(*partIt);
    }
    return cost;
  }

public:
  SolverTreewidthDP(const Molecule<Graph>* pMolecule)
    : Parent(pMolecule)
    , _g(pMolecule->getGraph())
    , _treeDecomposition(_g)
    , _blacklistMatrix(_treeDecomposition.getTree())
    , _optPartition(_treeDecomposition.getTree())
    , _optCost(_treeDecomposition.getTree())
    , _backTrace(_treeDecomposition.getTree())
  {
  }

  void solve(int resError);

  int getMaxBlacklistSize(const Node bag) const
  {
    /*int a = _maxAncestralBagSize[bag] 
        * std::min(_numberOfAncestors[bag], getBagSize(bag) -1) 
        * static_cast<int>(_pMolecule->getMaxChargeGroupSize());
    int b =
      static_cast<int>(_nodesSubTree[bag].size()) - 1;*/

    //if (a < 0 || b < 0)
    // std::cout << "Problem with bag " << getLabel(bag) << std::endl;
    //assert(a >= 0 && b >= 0);
    
    return std::min(
        _treeDecomposition.getMaxAncestralBagSize(bag)
        * std::min(_treeDecomposition.getNumberOfAncestors(bag),
                   _treeDecomposition.getBagSize(bag) - 1)
        * static_cast<int>(_pMolecule->getMaxChargeGroupSize())
      ,
        static_cast<int>(_treeDecomposition.getNodesSubTreeSize(bag) - 1));
  }

  void printBlacklistMatrix(std::ostream& out) const;

  int getBlacklistSize(const Node bag) const
  {
    int res = 0;

    const NodeSetIntPairSetMatrix& blacklist = _blacklistMatrix[bag];
    for (size_t i = 0; i < blacklist.size(); i++)
      res += static_cast<int>(blacklist[i].size());

    return res;
  }
};

template<typename GR>
void SolverTreewidthDP<GR>::init()
{
  // construct tree decomposition
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Computing tree decomposition..." << std::flush;
  }

  _treeDecomposition.init(TreeDecompositionType::TREE_DECOMPOSITION_GREEDY_FILL_IN);
  _treeDecomposition.root();

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << " Done! Treewidth is " << _treeDecomposition.getTreewidth() << std::endl;
  }

  const Graph& tree = _treeDecomposition.getTree();

  // initialize _blacklistMatrix[][0] to {}
  int count = 1;
  for (NodeIt bag(tree); bag != lemon::INVALID; ++bag, count++)
  {
    int idx = 0;
    const NodeSet& nodes = _treeDecomposition.getNodesSubTree(bag);
    const NodeSet& nodesInBag = _treeDecomposition.getNodesInBag(bag);
    const int n = getMaxBlacklistSize(bag);

    NodeSetIntPairSetMatrix& blacklistMatrix = _blacklistMatrix[bag];
    blacklistMatrix.clear();
    for (int i = 0; i <= n; i++)
      blacklistMatrix.push_back(NodeSetIntPairSet());
    blacklistMatrix[0].insert(NodeSetIntPair(NodeSet(), idx++));

    // build black list
    for (int i = 1; i <= n; i++)
    {
      const NodeSetIntPairSet& old_blacklist = blacklistMatrix[i-1];
      NodeSetIntPairSet& new_blacklist = blacklistMatrix[i];

      // pre-condition everything up to size i-1 has been constructed
      for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
      {       
        for (NodeSetIntPairSetIt nodeSetIt = old_blacklist.begin(); 
          nodeSetIt != old_blacklist.end(); nodeSetIt++)
        {
          NodeSetIntPair subset(nodeSetIt->first, idx);
          if (subset.first.find(*nodeIt) == subset.first.end())
          {
            subset.first.insert(*nodeIt);
            if (new_blacklist.find(subset) == new_blacklist.end() &&
                !std::includes(subset.first.begin(), subset.first.end(), 
                               nodesInBag.begin(), nodesInBag.end()))
            {
              if (g_verbosity >= VERBOSE_ESSENTIAL)
                std::cerr << "\rGenerating blacklist for bag "
                  << tree.id(bag) << " "
                  << count << "/" << _treeDecomposition.getBagSize(bag)
                  << ", " << i << "/" << n << ", " << idx << std::flush;
              
              new_blacklist.insert(subset);
              idx++;
            }
          }
        }
      }
    }
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
    std::cerr << std::endl;

  // initialize _optPartition
  for (NodeIt bag(tree); bag != lemon::INVALID; ++bag)
  {
    int n = getBlacklistSize(bag);
    _optPartition[bag] = PartitionVector(n);
    _optCost[bag] = DoubleVector(n, std::numeric_limits<double>::max());
    _backTrace[bag] = TracebackMatrix(n);
  }
}

template<typename GR>
void SolverTreewidthDP<GR>::findNonEmptyBags(const Node bag,
                                             const NodeSet& toRemove,
                                             BoolNodeMap& emptyFlag,
                                             NodeVector& nonEmptyBags)
{
  const NodeSet& nodesInBag = _treeDecomposition.getNodesInBag(bag);
  const RootedTree& rootedTree = _treeDecomposition.getRootedTree();
  
  NodeSet set;
  std::insert_iterator<NodeSet> insertIt(set, set.begin()); 
  std::set_difference(nodesInBag.begin(), nodesInBag.end(), 
                      toRemove.begin(), toRemove.end(),
                      insertIt);

  if (set.size() == 0)
  {
    emptyFlag[bag] = true;
  }
  else
  {
    emptyFlag[bag] = false;
    
    // get the parent
    const Node parentBag = rootedTree.source(RootedInArcIt(rootedTree, bag));
    if (emptyFlag[parentBag])
    {
      nonEmptyBags.push_back(bag);
    }
  }

  // recurse on children
  for (RootedOutArcIt a(rootedTree, bag); a != lemon::INVALID; ++a)
  {
    const Node childBag = rootedTree.target(a);
    findNonEmptyBags(childBag, toRemove, emptyFlag, nonEmptyBags);
  }
}

template<typename GR>
double SolverTreewidthDP<GR>::compute(const Node bag,
                                      const NodeSet& blacklist,
                                      const int blacklistIdx,
                                      const Partition& extendedPartition)
{
  const Graph& tree = _treeDecomposition.getTree();

  // we do a dfs on bag, to identify all non-empty bags j
  // (after removing from X_j, (A \cup P) in the subtree whose parent is empty
  NodeSet nodesInPartition;
  for (PartitionIt partIt = extendedPartition.begin(); 
    partIt != extendedPartition.end(); partIt++)
  {
    nodesInPartition.insert(partIt->begin(), partIt->end());
  }
  
  NodeSet toRemove = blacklist;
  toRemove.insert(nodesInPartition.begin(), nodesInPartition.end());
  
  NodeVector nonEmptyBags;
  BoolNodeMap emptyFlag(tree, false);
  findNonEmptyBags(bag, toRemove, emptyFlag, nonEmptyBags);

  // construct A \cup P
  NodeSet blacklistMask = blacklist;
  blacklistMask.insert(nodesInPartition.begin(), nodesInPartition.end());

  double cost = computeCost(extendedPartition);
  TracebackVector backTrace;
  for (NodeVectorIt bagIt = nonEmptyBags.begin(); bagIt != nonEmptyBags.end(); bagIt++)
  {
    const Node subBag = *bagIt;
    const NodeSet& subBagTree = _treeDecomposition.getNodesSubTree(subBag);

    // make the new blacklist = blacklistMask \cap subBagTree
    NodeSet subBagBlacklist;
    std::insert_iterator<NodeSet> insertIt(subBagBlacklist, subBagBlacklist.begin());
    std::set_intersection(blacklistMask.begin(), blacklistMask.end(), 
                          subBagTree.begin(), subBagTree.end(),
                          insertIt);

    // find the index corresponding to subBagBlacklist
    assert(_blacklistMatrix[subBag][subBagBlacklist.size()].find(std::make_pair(subBagBlacklist, 0)) != _blacklistMatrix[subBag][subBagBlacklist.size()].end());
    int subBagBlacklistIdx = _blacklistMatrix[subBag][subBagBlacklist.size()].find(std::make_pair(subBagBlacklist, 0))->second;
    
    cost += _optCost[*bagIt][subBagBlacklistIdx];
    backTrace.push_back(TracebackEntry(subBag, static_cast<int>(subBagBlacklist.size()), subBagBlacklistIdx));
  }

  if (cost < _optCost[bag][blacklistIdx])
  {
    _optCost[bag][blacklistIdx] = cost;
    _optPartition[bag][blacklistIdx] = extendedPartition;
    _backTrace[bag][blacklistIdx] = backTrace;
  }

  return cost;
}

template<typename GR>
void SolverTreewidthDP<GR>::extendPartition(const Node bag,
                                            const NodeSet& blacklist,
                                            const int blacklistIdx,
                                            NodeSet allowedNodes,
                                            Partition partitionToExtend, 
                                            Partition extendedPartition)
{
  if (partitionToExtend.size() == 0)
  {
    // base case
    if (g_verbosity >= VERBOSE_DEBUG)
    {
      std::cout << "Extended { ";
      bool firstPart = true;
      for (PartitionIt partIt = extendedPartition.begin(); partIt != extendedPartition.end(); partIt++)
      {
        if (firstPart)
          firstPart = false;
        else
          std::cout << ", ";
        
        std::cout << "{ ";
        bool firstNode = true;
        for (NodeSetIt nodeIt = partIt->begin(); nodeIt != partIt->end(); nodeIt++)
        {
          if (firstNode)
            firstNode = false;
          else
            std::cout << ", ";

          std::cout << _pMolecule->getLabel(*nodeIt);
        }
        std::cout << " }";
      }
      std::cout << " }" << std::flush;
    }

    double cost = compute(bag, blacklist, blacklistIdx, extendedPartition);
    
    if (g_verbosity >= VERBOSE_DEBUG)
      std::cout << " Cost: " << cost << std::endl;
  }
  else
  {
    const int k = static_cast<int>(_pMolecule->getMaxChargeGroupSize());
    
    // pick a part to extend and remove it from partitionToExtend
    NodeSet partToExtend = *partitionToExtend.begin();
    partitionToExtend.erase(partitionToExtend.begin());

    // pick a node inside this part to extend
    assert(partToExtend.size() > 0);
    Node v = *partToExtend.begin();

    NodeSetSetMatrix partsMatrix;
    {
      // partsMatrix[0] = { {} }
      NodeSetSet partSet;
      partSet.insert(NodeSet());
      partsMatrix.push_back(partSet);
    }
    {
      // partsMatrix[1] = { {v} }
      NodeSetSet partSet;
      NodeSet part;
      part.insert(v);
      partSet.insert(part);
      partsMatrix.push_back(partSet);
    }
    
    // make all parts containing v and up to size k
    for (NodeSetIt nodeIt = partToExtend.begin(); nodeIt != partToExtend.end(); nodeIt++)
    {
      allowedNodes.insert(*nodeIt);
    }
    
    for (int i = 2; i <= k; i++)
    {
      partsMatrix.push_back(NodeSetSet());
      const NodeSetSet& prevParts = partsMatrix[i-1];
      for (NodeSetSetIt partIt = prevParts.begin(); partIt != prevParts.end(); partIt++)
      {
        const NodeSet& part = *partIt;
        for (NodeSetIt nodeIt = part.begin(); nodeIt != part.end(); nodeIt++)
        {
          const Node node = *nodeIt;
          for (OutArcIt a(_g, node); a != lemon::INVALID; ++a)
          {
            Node newNode = _g.target(a);
            
            // only add newNode to a new part if it's not in current part
            // and it's in allowedNodes
            if (part.find(newNode) == part.end() &&
                allowedNodes.find(newNode) != allowedNodes.end())
            {
              NodeSet newPart = part;
              newPart.insert(newNode);
              partsMatrix[i].insert(newPart);
            }
          }
        }
      }
    }

    // now we have all parts, we simply have to check whether they're good and if so recurse
    for (int i = static_cast<int>(partToExtend.size()); i <= k; i++)
    {
      const NodeSetSet& parts = partsMatrix[i];
      for (NodeSetSetIt partIt = parts.begin(); partIt != parts.end(); partIt++)
      {
        const NodeSet& extendedPart = *partIt;
        if (std::includes(extendedPart.begin(), extendedPart.end(), partToExtend.begin(), partToExtend.end()))
        {
          // yes we got one!
          // let's recurse
          NodeSet newAllowedNodes;
          std::insert_iterator<NodeSet> insertIt(newAllowedNodes, newAllowedNodes.begin());
          std::set_difference(allowedNodes.begin(), allowedNodes.end(), 
                              extendedPart.begin(), extendedPart.end(), 
                              insertIt);

          Partition newExtendedPartition = extendedPartition;
          newExtendedPartition.insert(extendedPart);
          extendPartition(bag, blacklist, blacklistIdx, 
            newAllowedNodes, partitionToExtend, newExtendedPartition);
        }
      }
    }
  }
}

template<typename GR>
void SolverTreewidthDP<GR>::makeAndProcessPartition(const Node bag,
                                                    const NodeSet& blacklist,
                                                    const int blacklistIdx,
                                                    NodeVector availableNodes, 
                                                    const Partition& currentPartition)
{
  if (availableNodes.size() == 0)
  {
    // base case
    if (g_verbosity >= VERBOSE_DEBUG)
    {
      // for now: print the current partition
      std::cout << "{ ";
      bool firstPart = true;
      for (PartitionIt partIt = currentPartition.begin(); partIt != currentPartition.end(); partIt++)
      {
        if (firstPart)
          firstPart = false;
        else
          std::cout << ", ";
        
        std::cout << "{ ";
        bool firstNode = true;
        for (NodeSetIt nodeIt = partIt->begin(); nodeIt != partIt->end(); nodeIt++)
        {
          if (firstNode)
            firstNode = false;
          else
            std::cout << ", ";

          std::cout << _pMolecule->getLabel(*nodeIt);
        }
        std::cout << " }";
      }
      std::cout << " }" << std::endl;
    }

    NodeSet allowedNodes, tmpAllowedNodes;
    std::insert_iterator<NodeSet> insertIt(allowedNodes, allowedNodes.begin());
    std::insert_iterator<NodeSet> insertTmpIt(tmpAllowedNodes, tmpAllowedNodes.begin());

    const NodeSet& nodesInSubTreeBag = _treeDecomposition.getNodesSubTree(bag);
    std::set_difference(nodesInSubTreeBag.begin(), nodesInSubTreeBag.end(),
                        blacklist.begin(), blacklist.end(),
                        insertTmpIt);

    const NodeSet& nodesInBag = _treeDecomposition.getNodesInBag(bag);
    std::set_difference(tmpAllowedNodes.begin(), tmpAllowedNodes.end(),
                        nodesInBag.begin(), nodesInBag.end(),
                        insertIt);
    
    extendPartition(bag, blacklist, blacklistIdx, allowedNodes, currentPartition, Partition());
  }
  else
  {
    int maxPartSize = -1 +
      static_cast<int>(std::min(_pMolecule->getMaxChargeGroupSize(), availableNodes.size()));

    // pick a node
    Node v = availableNodes.front();
    availableNodes.erase(availableNodes.begin());

    // first process v in its own part
    Partition newPartition = currentPartition;
    NodeSet part;
    part.insert(v);
    newPartition.insert(part);
    makeAndProcessPartition(bag, blacklist, blacklistIdx, availableNodes, newPartition);

    for (int partSize = 1; partSize <= maxPartSize; partSize++)
    {
      SizeTVector currentPart;
      for (int j = 0; j < partSize; j++) 
        currentPart.push_back(j);

      do {
        NodeVector newAvailableNodes;

        NodeSet part;
        part.insert(v);

        for (int j = 0, k = 0; k < maxPartSize; k++)
        {
          if (j == partSize || k != static_cast<int>(currentPart[j]))
            newAvailableNodes.push_back(availableNodes[k]);
          else
          {
            part.insert(availableNodes[k]);
            j++;
          }
        }
        
        Partition newPartition = currentPartition;
        newPartition.insert(part);
        makeAndProcessPartition(bag, blacklist, blacklistIdx, newAvailableNodes, newPartition);
      } while (next(availableNodes, currentPart));
    }
  }
}

template<typename GR>
void SolverTreewidthDP<GR>::solve(Node bag, const NodeSetIntPair& blacklistPair)
{
  const Graph& tree = _treeDecomposition.getTree();
  const NodeSet& blacklist = blacklistPair.first;
  const int blacklistIdx = blacklistPair.second;

  const NodeSet& nodesInBag = _treeDecomposition.getNodesInBag(bag);
  NodeVector availableNodes;
  
  std::insert_iterator<NodeVector> insertIt(availableNodes, availableNodes.begin());
  // availableNodes = nodesInBag \ blacklist
  std::set_difference(nodesInBag.begin(), nodesInBag.end(), 
                      blacklist.begin(), blacklist.end(), 
                      insertIt);

  if (g_verbosity >= VERBOSE_DEBUG)
    std::cout << "Solving " << tree.id(bag) << " with blacklist " << blacklistIdx << std::endl;
  
  makeAndProcessPartition(bag, blacklist, blacklistIdx, availableNodes, Partition());
  
  // for each partition P over bag \ blacklist:
  //   for each extension P' of P (extension of subtree[bag] \ blacklist \ P)
}

template<typename GR>
void SolverTreewidthDP<GR>::reconstruct(const Node bag,
                                        const int blacklistSize,
                                        const int blacklistIdx,
                                        int& color)
{
  const Partition& optPartition = _optPartition[bag][blacklistIdx];

  for (PartitionIt partIt = optPartition.begin(); partIt != optPartition.end(); partIt++)
  {
    _chargeGroupVector.push_back(*partIt);
    // now color the nodes
    for (NodeSetIt nodeIt = partIt->begin(); nodeIt != partIt->end(); nodeIt++)
    {
      _color[*nodeIt] = color;
    }
    color++;
  }

  // recurse on the trace back stuff
  const TracebackVector& traceback = _backTrace[bag][blacklistIdx];
  for (TracebackVectorIt entryIt = traceback.begin(); entryIt != traceback.end(); entryIt++)
  {
    reconstruct(entryIt->_bag, entryIt->_blacklistSize, entryIt->_blacklistIdx, color);
  }
}

template<typename GR>
void SolverTreewidthDP<GR>::solve(int resError)
{
  init();

  // time to solve the DP, 
  // start from bags with the highest BFS levels and work towards the root
  const int maxBfsLevel = _treeDecomposition.getMaxBfsLevel();
  const int nBags = _treeDecomposition.getNumberOfBags();

  int count = 1;
  for (int i = maxBfsLevel; i >= 0; i--, count++)
  {
    const NodeVector& bags = _treeDecomposition.getBagsByBfsLevel(i);
    for (NodeVectorIt bagIt = bags.begin(); bagIt != bags.end(); bagIt++)
    {
      const NodeSetIntPairSetMatrix& blacklistMatrix = _blacklistMatrix[*bagIt];
      const size_t n = blacklistMatrix.size();
      for (size_t i = 0; i < n; i++)
      {
        const NodeSetIntPairSet& blacklistSet = blacklistMatrix[i];
        int idx = 1;
        for (NodeSetIntPairSetIt blacklistIt = blacklistSet.begin(); 
          blacklistIt != blacklistSet.end(); blacklistIt++, idx++)
        {
          if (g_verbosity >= VERBOSE_ESSENTIAL)
            std::cerr << "\rSolving bag "
              //<< getLabel(*bagIt) << " "
              << count << "/" << nBags
              << ", " << i << "/" << n << ", " << idx << std::flush;
          solve(*bagIt, *blacklistIt);
        }
      }
    }
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
    std::cerr << std::endl;


  // reconstruct optimal solution
  const Node rootBag = _treeDecomposition.getRootBag();

  _chargeGroupVector.clear();
  int color = 0;
  reconstruct(rootBag, 0, 0, color);
  updateCorrectedPartialCharges(resError);

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cerr << "Obj value: " << _optCost[rootBag][0] << std::endl;
  }
}

template<typename GR>
void SolverTreewidthDP<GR>::printBlacklistMatrix(std::ostream& out) const
{
  const Graph& tree = _treeDecomposition.getTree();
  for (NodeIt bag(tree); bag != lemon::INVALID; ++bag)
  {
    const NodeSetIntPairSetMatrix& blacklistMatrix = _blacklistMatrix[bag];
    const int n = getMaxBlacklistSize(bag);

    out << "Bag: " << getLabel(bag) 
      << ", max_card = " << n 
      << ", size = " << getBlacklistSize(bag)
      << std::endl;
    
    for (int i = 0; i <= n; i++)
    {
      const NodeSetIntPairSet& blacklist = blacklistMatrix[i];
      for (NodeSetIntPairSetIt nodeSetIt = blacklist.begin(); 
        nodeSetIt != blacklist.end(); nodeSetIt++)
      {
        const NodeSet& nodeSet = nodeSetIt->first;
        out << "{ ";
        bool first = true;
        for (NodeSetIt nodeIt = nodeSet.begin(); nodeIt != nodeSet.end(); nodeIt++)
        {
          if (!first)
            out << ", ";
          else
            first = false;
          out << _pMolecule->getLabel(*nodeIt);
        }
        out << " } " << nodeSetIt->second << std::endl;
      }
    }
  }
}

template<typename GR>
bool SolverTreewidthDP<GR>::next(const NodeVector& nodes, SizeTVector& currentPart)
{
  bool done = true;
  for (size_t i = nodes.size() - currentPart.size(); i < nodes.size(); i++)
    done &= (nodes[i] == nodes[currentPart[i - (nodes.size() - currentPart.size())]]);

  if (done)
    return false;

  // which index to increment
  size_t idxToIncr;
  size_t j = 1;
  for (size_t i = currentPart.size() - 1; i >= 0; i--, j++)
  {
    if (currentPart[i] < nodes.size() - j)
    {
      idxToIncr = i;
      break;
    }
  }

  // increment the index
  currentPart[idxToIncr]++;
  for (size_t i = idxToIncr + 1; i < currentPart.size(); i++)
    currentPart[i] = currentPart[i-1] + 1;

  return true;
}

#endif /* SOLVERTREEWIDTHDP_H_ */
