/* 
 * solvertreewidthmemdp.h
 *
 *  Created on: 27-oct-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVERTREEWIDTHMEMDP_H_
#define SOLVERTREEWIDTHMEMDP_H_

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
#include <map>

namespace nina
{
namespace cgp
{

template<typename GR>
class SolverTreewidthMemDP : public Solver<GR>
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
  /// Partition type
  typedef typename std::set<NodeSet> Partition;
  /// Partition type iterator
  typedef typename Partition::const_iterator PartitionIt;
  /// Partition map type
  typedef typename std::map<NodeSet, Partition> PartitionMap;
  /// Partition vector type iterator
  typedef typename PartitionMap::const_iterator PartitionMapIt;
  /// Partition matrix map type
  typedef typename Graph::template NodeMap<PartitionMap> PartitionMapMap;
  /// Int vector type
  typedef typename std::vector<int> IntVector;
  /// Int vector map type
  typedef typename Graph::template NodeMap<IntVector> IntVectorMap;
  /// Int vector type
  typedef typename std::map<NodeSet, double> DoubleMap;
  /// Int vector map type
  typedef typename Graph::template NodeMap<DoubleMap> DoubleMapMap;
  /// Size_t vector type
  typedef typename std::vector<size_t> SizeTVector;

  struct TracebackEntry
  {
    Node _bag;
    NodeSet _blacklist;

    TracebackEntry()
      : _bag(lemon::INVALID)
      , _blacklist()
    {
    }

    TracebackEntry(const Node bag, NodeSet blacklist)
      : _bag(bag)
      , _blacklist(blacklist)
    {
    }
  };

  /// Traceback vector type
  typedef typename std::vector<TracebackEntry> TracebackVector;
  /// Traceback vector type iterator
  typedef typename TracebackVector::const_iterator TracebackVectorIt;
  /// Traceback map vector type
  typedef typename std::map<NodeSet, TracebackVector> TracebackVectorMap;
  /// Traceback map type
  typedef typename Graph::template NodeMap<TracebackVectorMap> TracebackVectorMapMap;

private:
  /// The original graph
  const Graph& _g;
  /// Tree decomposition
  TreeDecompositionType _treeDecomposition;

  /// Holds optimal partition given a bag and a blacklist
  PartitionMapMap _optPartition;
  /// Holds optimal score given a bag and a blacklist
  DoubleMapMap _optCost;
  /// Traceback map
  TracebackVectorMapMap _backTrace;

  /// Initializes tree decomposition
  /// Initializes _cgMatrixMap, _optPartition and _optPartitionCg
  void init();

  void solve(Node bag, const NodeSet& blacklist);

  /// Returns the next part having same cardinality, returns false if there is none
  static bool next(const NodeVector& nodes, SizeTVector& currentPart);

  void makeAndProcessPartition(const Node bag,
                               const NodeSet& blacklist,
                               NodeVector availableNodes, 
                               const Partition& currentPartition);

  void extendPartition(const Node bag,
                       const NodeSet& blacklist,
                       NodeSet allowedNodes,
                       Partition partitionToExtend, 
                       Partition extendedPartition);

  double compute(const Node bag,
                 const NodeSet& blacklist,
                 const Partition& extendedPartition);

  void findNonEmptyBags(const Node bag,
                        const NodeSet& toRemove,
                        BoolNodeMap& emptyFlag,
                        NodeVector& nonEmptyBags);

  void reconstruct(const Node bag,
                   const NodeSet& blacklist,
                   int& color);

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

  bool isComputed(const Node bag, const NodeSet& blacklist) const
  {
    const DoubleMap& map = _optCost[bag];
    return map.find(blacklist) != map.end();
  }

public:
  SolverTreewidthMemDP(const Molecule<Graph>* pMolecule)
    : Parent(pMolecule)
    , _g(pMolecule->getGraph())
    , _treeDecomposition(_g)
    , _optPartition(_treeDecomposition.getTree())
    , _optCost(_treeDecomposition.getTree())
    , _backTrace(_treeDecomposition.getTree())
  {
  }

  bool solve(int resError);
};
} // namespace cgp
} // namespace nina

template<typename GR>
void nina::cgp::SolverTreewidthMemDP<GR>::init()
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
}

template<typename GR>
void nina::cgp::SolverTreewidthMemDP<GR>::findNonEmptyBags(const Node bag,
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

    // recurse on children
    for (RootedOutArcIt a(rootedTree, bag); a != lemon::INVALID; ++a)
    {
      const Node childBag = rootedTree.target(a);
      findNonEmptyBags(childBag, toRemove, emptyFlag, nonEmptyBags);
    }
  }
  else
  {
    emptyFlag[bag] = false;
    
    // get the parent
    const Node parentBag = (bag == _treeDecomposition.getRootBag() ?
        lemon::INVALID :
        rootedTree.source(RootedInArcIt(rootedTree, bag)));

    if (parentBag == lemon::INVALID || emptyFlag[parentBag])
    {
      nonEmptyBags.push_back(bag);
    }
    else
    {
      // recurse on children
      for (RootedOutArcIt a(rootedTree, bag); a != lemon::INVALID; ++a)
      {
        const Node childBag = rootedTree.target(a);
        findNonEmptyBags(childBag, toRemove, emptyFlag, nonEmptyBags);
      }
    }
  }
}

template<typename GR>
double nina::cgp::SolverTreewidthMemDP<GR>::compute(const Node bag,
                                                    const NodeSet& blacklist,
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
    
    solve(subBag, subBagBlacklist);
    cost += _optCost[subBag][subBagBlacklist];
    backTrace.push_back(TracebackEntry(subBag, subBagBlacklist));
  }

  if (!isComputed(bag, blacklist) || cost < _optCost[bag][blacklist])
  {
    _optCost[bag][blacklist] = cost;
    _optPartition[bag][blacklist] = extendedPartition;
    _backTrace[bag][blacklist] = backTrace;
  }

  return cost;
}

template<typename GR>
void nina::cgp::SolverTreewidthMemDP<GR>::extendPartition(const Node bag,
                                                          const NodeSet& blacklist,
                                                          NodeSet allowedNodes,
                                                          Partition partitionToExtend,
                                                          Partition extendedPartition)
{
  if (partitionToExtend.size() == 0)
  {
    const Graph& tree = _treeDecomposition.getTree();
    double cost = compute(bag, blacklist, extendedPartition);

    // base case
    if (g_verbosity >= VERBOSE_DEBUG)
    {
      std::cout << "Bag: " << tree.id(bag) << ", Blacklist: { ";
      bool firstNode = true;
      for (NodeSetIt nodeIt = blacklist.begin(); nodeIt != blacklist.end(); nodeIt++)
      {
        if (firstNode)
          firstNode = false;
        else
          std::cout << ", ";

        std::cout << _pMolecule->getLabel(*nodeIt);
      }
      
      std::cout << " }, Extended: { ";
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
      std::cout << " }" << " Cost: " << cost << std::endl;
    }
  }
  else
  {
    // pick a part to extend and remove it from partitionToExtend
    const NodeSet partToExtend = *partitionToExtend.begin();
    partitionToExtend.erase(partitionToExtend.begin());

    // determine k
    int k = std::numeric_limits<int>::max();
    for (NodeSetIt nodeIt = partToExtend.begin(); nodeIt != partToExtend.end(); nodeIt++)
    {
      if (_pMolecule->getMaxChargeGroupSize(*nodeIt) < k)
      {
        k = _pMolecule->getMaxChargeGroupSize(*nodeIt);
      }
    }

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
    allowedNodes.insert(partToExtend.begin(), partToExtend.end());

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

              if (_pMolecule->isPotFeasible(newPart))
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

        if (std::includes(extendedPart.begin(), extendedPart.end(), partToExtend.begin(), partToExtend.end()) &&
            _pMolecule->isFeasible(extendedPart))
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
          extendPartition(bag, blacklist,
            newAllowedNodes, partitionToExtend, newExtendedPartition);
        }
      }
    }
  }
}

template<typename GR>
void nina::cgp::SolverTreewidthMemDP<GR>::makeAndProcessPartition(const Node bag,
                                                                  const NodeSet& blacklist,
                                                                  NodeVector availableNodes,
                                                                  const Partition& currentPartition)
{
  if (availableNodes.size() == 0)
  {
    // base case
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

    extendPartition(bag, blacklist, allowedNodes, currentPartition, Partition());
  }
  else
  {
    // pick a node
    Node v = availableNodes.front();

    int maxPartSize = -1 +
      static_cast<int>(std::min(_pMolecule->getMaxChargeGroupSize(v), availableNodes.size()));

    availableNodes.erase(availableNodes.begin());

    // first process v in its own part
    Partition newPartition = currentPartition;
    NodeSet part;
    part.insert(v);
    newPartition.insert(part);
    makeAndProcessPartition(bag, blacklist, availableNodes, newPartition);

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
          {
            newAvailableNodes.push_back(availableNodes[k]);
          }
          else
          {
            part.insert(availableNodes[k]);
            j++;
          }
        }
        
        Partition newPartition = currentPartition;
        newPartition.insert(part);
        makeAndProcessPartition(bag, blacklist, newAvailableNodes, newPartition);
      } while (next(availableNodes, currentPart));
    }
  }
}

template<typename GR>
void nina::cgp::SolverTreewidthMemDP<GR>::solve(Node bag, const NodeSet& blacklist)
{
  if (!isComputed(bag, blacklist))
  {
    const NodeSet& nodesInBag = _treeDecomposition.getNodesInBag(bag);
    NodeVector availableNodes;
    
    std::insert_iterator<NodeVector> insertIt(availableNodes, availableNodes.begin());
    // availableNodes = nodesInBag \ blacklist
    std::set_difference(nodesInBag.begin(), nodesInBag.end(), 
                        blacklist.begin(), blacklist.end(), 
                        insertIt);

    makeAndProcessPartition(bag, blacklist, availableNodes, Partition());
    
    // for each partition P over bag \ blacklist:
    //   for each extension P' of P (extension of subtree[bag] \ blacklist \ P)
  }
}

template<typename GR>
void nina::cgp::SolverTreewidthMemDP<GR>::reconstruct(const Node bag,
                                                      const NodeSet& blacklist,
                                                      int& color)
{
  const Partition& optPartition = _optPartition[bag][blacklist];

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
  const TracebackVector& traceback = _backTrace[bag][blacklist];
  for (TracebackVectorIt entryIt = traceback.begin(); entryIt != traceback.end(); entryIt++)
  {
    reconstruct(entryIt->_bag, entryIt->_blacklist, color);
  }
}

template<typename GR>
bool nina::cgp::SolverTreewidthMemDP<GR>::solve(int resError)
{
  init();

  const Node rootBag = _treeDecomposition.getRootBag();
  solve(rootBag, NodeSet());

  // reconstruct optimal solution
  _chargeGroupVector.clear();
  int color = 0;
  reconstruct(rootBag, NodeSet(), color);
  updateCorrectedPartialCharges(resError);

  bool solved = true;
  for (NodeIt v(_pMolecule->getGraph()); v != lemon::INVALID; ++v)
    solved &= _color[v] != -1;

  if (g_verbosity >= VERBOSE_ESSENTIAL && solved)
  {
    std::cerr << "Obj value: " << _optCost[rootBag][NodeSet()] << std::endl;
  }

  return solved;
}

template<typename GR>
bool nina::cgp::SolverTreewidthMemDP<GR>::next(const NodeVector& nodes, SizeTVector& currentPart)
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

#endif /* SOLVERTREEWIDTHMEMDP_H_ */
