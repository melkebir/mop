/*
 * product.h
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 */

// TODO:
// - Remove deg-1 nodes

#ifndef PRODUCT_H
#define PRODUCT_H

#include "molecule.h"
#include <set>
#include <vector>

namespace nina {
namespace cgp {

template<typename GR>
class Product
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  typedef Molecule<Graph> MoleculeType;

private:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Graph::template NodeMap<Node> NodeNodeMap;
  typedef std::multiset<int> IntSet;
  typedef typename Graph::template NodeMap<IntSet> IntSetNodeMap;
  typedef std::vector<Node> NodeVector;
  typedef typename Graph::template NodeMap<NodeVector> NodeVectorMap;

public:
  Product(const MoleculeType& mol1, const MoleculeType& mol2, int shell)
    : _mol1(mol1)
    , _mol2(mol2)
    , _shell(shell)
    , _g()
    , _mol1ToG(mol1.getGraph(), lemon::INVALID)
    , _mol2ToG(mol2.getGraph(), lemon::INVALID)
    , _gToMol1(_g)
    , _gToMol2(_g)
    , _connectivityEdge(_g)
    , _g1ToDeg1Neighbors(mol1.getGraph())
    , _g2ToDeg1Neighbors(mol2.getGraph())
    , _numNodes(0)
    , _numEdges(0)
  {
    generate();
  }

  void printDOT(std::ostream& out) const;

  void printProductNodeJSON(Node uv, std::ostream& out) const
  {
    Node u = _gToMol1[uv];
    Node v = _gToMol2[uv];

    out << "            {" << std::endl
        << "              id1: " << _mol1.getLabel(u) << "," << std::endl
        << "              id2: " << _mol2.getLabel(v) << "," << std::endl
        << "              charge: " << _mol2.getPartialCharge(v) << std::endl
        << "            }";

    const NodeVector& uNeighbors = _g1ToDeg1Neighbors[u];
    const NodeVector& vNeighbors = _g2ToDeg1Neighbors[v];

    assert(uNeighbors.size() == vNeighbors.size());
    for (size_t i = 0; i < uNeighbors.size(); ++i)
    {
      out << "," << std::endl;
      out << "            {" << std::endl
          << "              id1: " << _mol1.getLabel(uNeighbors[i]) << "," << std::endl
          << "              id2: " << _mol2.getLabel(vNeighbors[i]) << "," << std::endl
          << "              charge: " << _mol2.getPartialCharge(vNeighbors[i]) << std::endl
          << "            }";
    }
  }

  void printProductNode(Node uv, std::ostream& out) const
  {
    Node u = _gToMol1[uv];
    Node v = _gToMol2[uv];

    out << "[" << _g.id(uv) << ": " << _mol1.getLabel2(u)
        << " (" << _mol1.getLabel(u) << ") , "
        << _mol2.getLabel2(v)
        << " (" << _mol2.getLabel(v) << ")]";

    const NodeVector& uNeighbors = _g1ToDeg1Neighbors[u];
    const NodeVector& vNeighbors = _g2ToDeg1Neighbors[v];

    assert(uNeighbors.size() == vNeighbors.size());
    for (size_t i = 0; i < uNeighbors.size(); ++i)
    {
      out << ", ";
      out << "[" << _mol1.getLabel2(uNeighbors[i])
          << " (" << _mol1.getLabel(uNeighbors[i]) << ") , "
          << _mol2.getLabel2(vNeighbors[i])
          << " (" << _mol2.getLabel(vNeighbors[i]) << ")]";
    }
  }

  void printProductNodeVector(const std::vector<Node>& nodes,
                              std::ostream& out) const
  {
    out << nodes.size() << ": ";
    bool first = true;
    for (typename std::vector<Node>::const_iterator it = nodes.begin();
         it != nodes.end(); ++it)
    {
      if (first)
      {
        first = false;
      }
      else
      {
        out << ", ";
      }
      printProductNode(*it, out);
    }
    out << std::endl;
  }

  void printProductNodeVectorJSON(const std::vector<Node>& nodes,
                                  std::ostream& out,
                                  const bool& fst) const
  {
    if (fst)
    {
      out << "          pairs: [";
    }
    else
    {
      out << "," << std::endl << "          pairs: [";
    }
    bool first = true;
    for (typename std::vector<Node>::const_iterator it = nodes.begin();
         it != nodes.end(); ++it)
    {
      if (first)
      {
        first = false;
        out << std::endl;
      }
      else
      {
        out << "," << std::endl;
      }
      printProductNodeJSON(*it, out);
    }
    out << std::endl << "          ]";
  }

  int getNumNodes() const { return _numNodes; }
  int getNumEdges() const { return _numEdges; }

  const Graph & getGraph() const { return _g; }

  bool connectivityEdge(Edge e) const { return _connectivityEdge[e]; }

private:
  const MoleculeType& _mol1;
  const MoleculeType& _mol2;
  const int _shell;
  Graph _g;
  NodeNodeMap _mol1ToG;
  NodeNodeMap _mol2ToG;
  NodeNodeMap _gToMol1;
  NodeNodeMap _gToMol2;
  BoolEdgeMap _connectivityEdge;
  NodeVectorMap _g1ToDeg1Neighbors;
  NodeVectorMap _g2ToDeg1Neighbors;

  int _numNodes;
  int _numEdges;

  void generate();

  void generateAtomTypes(const MoleculeType& mol, IntSetNodeMap& intSet);

  void generateDeg1NeighborSet(const Graph& g,
                               const IntNodeMap& deg,
                               NodeVectorMap& deg1NeighborMap);

  void dfs(const Node v, const int depth,
           const MoleculeType& mol,
           BoolNodeMap& visited,
           IntSet& s);

  void determineDegrees(const Graph& g, IntNodeMap& deg);
};

template<typename GR>
inline void Product<GR>::dfs(const Node v, const int depth,
                             const MoleculeType& mol,
                             BoolNodeMap& visited,
                             IntSet& s)
{
  const Graph& g = mol.getGraph();
  visited[v] = true;
  s.insert(mol.getAtomType(v));

  if (depth < _shell)
  {
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      Node w = g.oppositeNode(v, e);
      if (!visited[w])
      {
        dfs(w, depth + 1, mol, visited, s);
      }
    }
  }
}

template<typename GR>
inline void Product<GR>::determineDegrees(const Graph& g, IntNodeMap& deg)
{
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      ++deg[v];
    }
  }
}

template<typename GR>
inline void Product<GR>::generateDeg1NeighborSet(const Graph& g,
                                                 const IntNodeMap& deg,
                                                 NodeVectorMap& deg1NeighborMap)
{
  for (NodeIt u(g); u != lemon::INVALID; ++u)
  {
    if (deg[u] > 1)
    {
      for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e)
      {
        Node v = g.oppositeNode(u, e);
        if (deg[v] == 1) {
          deg1NeighborMap[u].push_back(v);
        }
      }
    }
  }

}

template<typename GR>
inline void Product<GR>::generate()
{
  const Graph& g1 = _mol1.getGraph();
  const Graph& g2 = _mol2.getGraph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  IntSetNodeMap set1(g1);
  IntSetNodeMap set2(g2);

  IntNodeMap deg1(g1, 0);
  determineDegrees(g1, deg1);
  IntNodeMap deg2(g2, 0);
  determineDegrees(g2, deg2);

  // determine degrees
  generateAtomTypes(_mol1, set1);
  generateAtomTypes(_mol2, set2);

  // generate nodes
  for (NodeIt u(g1); u != lemon::INVALID; ++u)
  {
    for (NodeIt v(g2); v != lemon::INVALID; ++v)
    {
      if (_mol1.getAtomType(u) == _mol2.getAtomType(v) && set1[u] == set2[v])
      {
        assert(deg1[u] == deg2[v]);
        // don't add product nodes for a pair of deg-1 nodes unless _shell == 0
        if (_shell == 0 || deg1[u] > 1)
        {
          Node uv = _g.addNode();
          _mol1ToG[u] = uv;
          _mol2ToG[v] = uv;
          _gToMol1[uv] = u;
          _gToMol2[uv] = v;
          ++_numNodes;
        }
      }
    }
  }

  // generate deg1 sets
  generateDeg1NeighborSet(g1, deg1, _g1ToDeg1Neighbors);
  generateDeg1NeighborSet(g2, deg2, _g2ToDeg1Neighbors);

  // generate edges
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1)
  {
    Node u1 = _gToMol1[u1v1];
    Node v1 = _gToMol2[u1v1];
    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2)
    {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _gToMol1[u2v2];
      Node v2 = _gToMol2[u2v2];

      assert(_mol1.getAtomType(u1) == _mol2.getAtomType(v1));

      if (u1 != u2 && v1 != v2)
      {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 == v1v2)
        {
          _connectivityEdge[_g.addEdge(u1v1, u2v2)] = u1u2;
          ++_numEdges;
        }
      }
    }
  }
}

template<typename GR>
inline void Product<GR>::generateAtomTypes(const MoleculeType& mol,
                                           IntSetNodeMap& intSet)
{
  const Graph& g = mol.getGraph();
  BoolNodeMap visited(g, false);
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    lemon::mapFill(g, visited, false);
    IntSet& s = intSet[v];
    dfs(v, 0, mol, visited, s);
  }
}

template<typename GR>
inline void Product<GR>::printDOT(std::ostream& out) const
{
  // header
  out << "graph G {" << std::endl
      << "\toverlap=scale" << std::endl
      << "\tlayout=neato" << std::endl;

  // nodes
  for (NodeIt uv(_g); uv != lemon::INVALID; ++uv)
  {
    Node u = _gToMol1[uv];
    Node v = _gToMol2[uv];
    out << "\t" << _g.id(uv) << " [label=\"[" << _mol1.getLabel2(u) << " (" << _mol1.getLabel(u) << ") , "
        << _mol2.getLabel2(v) << " (" << _mol2.getLabel(v) << ")]" << "\\n" << _g.id(uv) << "\"]" << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
    if (connectivityEdge(e))
    {
      out << " [color=red]";
    }
    out << std::endl;
  }

  out << "}" << std::endl;
}


} // namespace cgp
} // namespace nina

#endif // PRODUCT_H
