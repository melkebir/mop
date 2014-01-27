/*
 * symmetrizationautomorphism.h
 *
 *  Created on: 4-apr-2013
 *      Author: M. El-Kebir
 */

#ifndef SYMMETRIZATIONAUTOMORPHISM_H
#define SYMMETRIZATIONAUTOMORPHISM_H

#include <lemon/smart_graph.h>
#include <lemon/bfs.h>
#include "symmetrization.h"
#include "GNA/io/input/matchinggraph.h"
#include "GNA/io/input/bpidentityparser.h"
#include "GNA/natalie.h"
#include "common/io/input/identityparser.h"

namespace nina {
namespace cgp {

template<typename GR>
class SymmetrizationAutomorphism : public Symmetrization<GR>
{
public:
  /// Graph type
  typedef GR Graph;
  /// Parent type
  typedef Symmetrization<GR> Parent;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Bipartite graph type
  typedef lemon::SmartBpGraph BpGraph;

public:
  typedef typename Parent::WeightNodeMap WeightNodeMap;
  typedef typename Parent::MoleculeType MoleculeType;

protected:
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeMatrix NodeMatrix;
  typedef typename Parent::NodeMatrixIt NodeMatrixIt;
  typedef typename Parent::NodeMatrixNonConstIt NodeMatrixNonConstIt;
  typedef typename Parent::AtomList AtomList;
  typedef typename Parent::AtomListIt AtomListIt;
  typedef typename Parent::AtomListMultiSet AtomListMultiSet;
  typedef typename Parent::AtomListSetMultiIt AtomListSetMultiIt;
  typedef typename Parent::AtomListMultiSetMap AtomListMultiSetMap;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;

  typedef IdentityParser<Graph> IdentityParserType;
  typedef gna::BpIdentityParser<Graph, BpGraph> BpIdentityParserType;
  typedef gna::MatchingGraph<Graph, BpGraph> MatchingGraphType;
  typedef gna::Natalie<Graph, BpGraph> NatalieType;
  typedef gna::GlobalProblemConstrained<Graph, BpGraph> GlobalProblemConstrainedType;
  typedef lemon::Bfs<Graph> BfsType;
  typedef typename BpGraph::Node BpNode;
  typedef typename BpGraph::Edge BpEdge;
  typedef typename BpGraph::NodeIt BpNodeIt;
  typedef typename BpGraph::EdgeIt BpEdgeIt;
  typedef typename BpGraph::IncEdgeIt BpIncEdgeIt;
  typedef typename BpGraph::RedNode BpRedNode;
  typedef typename BpGraph::BlueNode BpBlueNode;
  typedef typename BpGraph::RedIt BpRedIt;
  typedef typename BpGraph::BlueIt BpBlueIt;
  typedef typename BpGraph::template NodeMap<Node> MatchNodeToOrigNodeMap;
  typedef typename Graph::template NodeMap<BpNode> OrigNodeToMatchNodeMap;
  typedef typename Graph::template NodeMap<Node> NodeNodeMap;

  using Parent::_neighborhoodSize;
  using Parent::_maxAbsDifference;
  using Parent::_eqClasses;

private:
  NatalieType* construct(const MoleculeType& molecule,
                         Node v,
                         Node w);

  void constructGraph(const MoleculeType& molecule,
                      Node root,
                      Graph& g_out,
                      IntNodeMap& atomType,
                      NodeNodeMap& map);

  BpEdge constructBpGraph(const Graph& g1,
                          const IntNodeMap& atomTypeG1,
                          const Graph& g2,
                          const IntNodeMap& atomTypeG2,
                          const Node root1,
                          const Node root2,
                          BpGraph& gm,
                          MatchNodeToOrigNodeMap& gmToG12,
                          OrigNodeToMatchNodeMap& g1ToGm,
                          OrigNodeToMatchNodeMap& g2ToGm);

public:
  SymmetrizationAutomorphism(int neighborhoodSize, double maxAbsDifference)
    : Parent(neighborhoodSize, maxAbsDifference)
  {
  }

  ~SymmetrizationAutomorphism()
  {
  }

  virtual bool equivalent(const MoleculeType& molecule,
                          Node v,
                          Node w);
};

template<typename GR>
inline void SymmetrizationAutomorphism<GR>::constructGraph(const MoleculeType& molecule,
                                                           Node root,
                                                           Graph& g_out,
                                                           IntNodeMap& atomType,
                                                           NodeNodeMap& map)
{
  const Graph& g = molecule.getGraph();

  BfsType bfs(g);
  bfs.init();
  bfs.run(root);

  g_out.clear();
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    if (_neighborhoodSize == -1 || bfs.dist(v) <= _neighborhoodSize)
    {
      Node new_v = g_out.addNode();
      atomType[new_v] = molecule.getAtomType(v);
      map[v] = new_v;
    }
    else
    {
      map[v] = lemon::INVALID;
    }
  }

  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = map[g.u(e)];
    Node v = map[g.v(e)];
    if (u != lemon::INVALID && v != lemon::INVALID)
    {
      g_out.addEdge(u, v);
    }
  }
}

template<typename GR>
inline typename SymmetrizationAutomorphism<GR>::BpEdge
SymmetrizationAutomorphism<GR>::constructBpGraph(const Graph& g1,
                                                 const IntNodeMap& atomTypeG1,
                                                 const Graph& g2,
                                                 const IntNodeMap& atomTypeG2,
                                                 const Node root1,
                                                 const Node root2,
                                                 BpGraph& gm,
                                                 MatchNodeToOrigNodeMap& gmToG12,
                                                 OrigNodeToMatchNodeMap& g1ToGm,
                                                 OrigNodeToMatchNodeMap& g2ToGm)
{
  gm.clear();

  for (NodeIt v(g1); v != lemon::INVALID; ++v)
  {
    BpNode new_v = gm.addRedNode();
    gmToG12[new_v] = v;
    g1ToGm[v] = new_v;
  }

  for (NodeIt v(g2); v != lemon::INVALID; ++v)
  {
    BpNode new_v = gm.addBlueNode();
    gmToG12[new_v] = v;
    g2ToGm[v] = new_v;
  }

  BpEdge res = lemon::INVALID;

  for (BpRedIt x(gm); x != lemon::INVALID; ++x)
  {
    for (BpBlueIt y(gm); y != lemon::INVALID; ++y)
    {
      Node org_x = gmToG12[x];
      Node org_y = gmToG12[y];

      if (atomTypeG1[org_x] == atomTypeG2[org_y])
      {
        BpEdge e = gm.addEdge(x, y);
        if (org_x == root1 && org_y == root2)
          res = e;
      }
    }
  }

  return res;
}

template<typename GR>
inline typename SymmetrizationAutomorphism<GR>::NatalieType*
SymmetrizationAutomorphism<GR>::construct(const MoleculeType& molecule,
                                          Node v,
                                          Node w)
{
  const Graph& g = molecule.getGraph();

  Graph g1;
  IntNodeMap atomTypeG1(g1);
  NodeNodeMap mapToG1(g);

  Graph g2;
  IntNodeMap atomTypeG2(g2);
  NodeNodeMap mapToG2(g);

  BpGraph gm;
  MatchNodeToOrigNodeMap mapGmToG12(gm);
  OrigNodeToMatchNodeMap mapG1ToGm(g1);
  OrigNodeToMatchNodeMap mapG2ToGm(g2);

  constructGraph(molecule, v, g1, atomTypeG1, mapToG1);
  constructGraph(molecule, w, g2, atomTypeG2, mapToG2);
  BpEdge toFix = constructBpGraph(g1,
                                  atomTypeG1,
                                  g2,
                                  atomTypeG2,
                                  mapToG1[v],
                                  mapToG2[w],
                                  gm,
                                  mapGmToG12,
                                  mapG1ToGm,
                                  mapG2ToGm);

  // now let's initialize _pMatchingGraph
  typename NatalieType::Options options;
  options._discretizeWeight = true;
  options._integral = true;
  NatalieType* pResult = new NatalieType(options);

  IdentityParserType parser_g1(g1, NULL, NULL, NULL);
  IdentityParserType parser_g2(g2, NULL, NULL, NULL);
  BpIdentityParserType parser_gm(gm, mapGmToG12, NULL, &parser_g1, &parser_g2);

  pResult->initMatchingGraph(&parser_g1, &parser_g2, &parser_gm);

  toFix = parser_gm.map(toFix);
  std::set<BpEdge> fixedBpEdge;
  fixedBpEdge.insert(toFix);

  GlobalProblemConstrainedType* pGlobal =
      new GlobalProblemConstrainedType(*pResult->getMatchingGraph(),
                                       *pResult->getScoreModel(),
                                       fixedBpEdge);

  pResult->setGlobalProblem(pGlobal);
  pResult->init();

  return pResult;
}

template<typename GR>
inline bool SymmetrizationAutomorphism<GR>::equivalent(const MoleculeType& molecule,
                                                       Node v,
                                                       Node w)
{
  if (molecule.getAtomType(v) == molecule.getAtomType(w))
  {
    const Graph& g = molecule.getGraph();
    int deg_v = 0, deg_w = 0;
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e) deg_v++;
    for (IncEdgeIt e(g, w); e != lemon::INVALID; ++e) deg_w++;

    if (deg_v == 1 && deg_w == 1 && g.oppositeNode(v, IncEdgeIt(g, v)) == g.oppositeNode(w, IncEdgeIt(g, w)))
    {
      // degree 1 nodes with same parent and type are symmetric
      return true;
    }
    else if (deg_v != deg_w)
    {
      // if the degree is different, no symmetry
      return false;
    }

    VerbosityLevel verbosityCpy = g_verbosity;
    g_verbosity = VERBOSE_NONE;
    NatalieType* pNatalie = construct(molecule, v, w);
    double nEdgesG1 = pNatalie->getMatchingGraph()->getEdgeCountG1();
    double nEdgesG2 = pNatalie->getMatchingGraph()->getEdgeCountG2();

    bool exists = false;
    if (pNatalie->exists(std::min(nEdgesG1, nEdgesG2), exists))
    {
      std::cerr << "Unable to determine whether " << molecule.getLabel(v)
                << " and " << molecule.getLabel(w) << " belong to the same symmetry group"
                << std::endl;
    }

    //pNatalie->solve();
    //bool exists = std::min(nEdgesG1, nEdgesG2) == pNatalie->getScore();

    g_verbosity = verbosityCpy;

    return exists;
  }
  else
  {
    return false;
  }
}

} // namespace cgp
} // namespace nina

#endif // SYMMETRIZATIONAUTOMORPHISM_H
