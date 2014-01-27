/* moleculegraphviz.h
 *
 *  Created on: 30-Jul-2012
 *      Author: M. El-Kebir
 */

#ifndef MOLECULEGRAPHVIZ_H_
#define MOLECULEGRAPHVIZ_H_

#include "graphviz/gvc.h"
#include "moleculeviz.h"

template<typename GR>
class MoleculeGraphviz : public MoleculeViz<GR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Base class type
  typedef MoleculeViz<GR> Parent;
  /// Base of the base class type
  typedef Molecule<GR> GrandParent;

  using Parent::_g;
  using Parent::_partialCharge;
  using Parent::_label;
  using Parent::_label2;
  using Parent::_maxChargeGroupSize;
  using Parent::_radius;
  using Parent::_coordX;
  using Parent::_coordY;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Mapping to graphviz nodes type
  typedef typename Graph::template NodeMap<Agnode_t*> GraphvizNodeMap;
  /// Mapping to graphviz edges type
  typedef typename Graph::template EdgeMap<Agedge_t*> GraphvizEdgeMap;

  GraphvizNodeMap _graphvizNode;
  GraphvizEdgeMap _graphvizEdge;

public:
  MoleculeGraphviz(size_t k, double radius)
    : Parent(k)
    , _graphvizNode(_g)
    , _graphvizEdge(_g)
  {
    _radius = radius;
  }

  virtual void readLGF(std::istream& in);
  bool performLayout();
};

template<typename GR>
void MoleculeGraphviz<GR>::readLGF(std::istream& in)
{
  GrandParent::readLGF(in);
  performLayout();
}

template<typename GR>
bool MoleculeGraphviz<GR>::performLayout()
{
  bool res = true;

  GVC_t* pGvc;
  pGvc = gvContext();

  Agraph_t* pGraph = agopen("molecule", AGRAPHSTRICT);
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    _graphvizNode[n] = agnode(pGraph, const_cast<char*>(_label[n].c_str()));
    agset(_graphvizNode[n], "label", const_cast<char*>(_label2[n].c_str()));

    char buf[1024];
    sprintf(buf, "%lf", _radius / 72.0);
    agset(_graphvizNode[n], "width", buf);
    agset(_graphvizNode[n], "height", buf);
    agset(_graphvizNode[n], "fixedsize", "true");
  }

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    _graphvizEdge[e] = 
      agedge(pGraph, _graphvizNode[_g.u(e)], _graphvizNode[_g.v(e)]);
  }

  gvLayout(pGvc, pGraph, "neato");
  gvRender(pGvc, pGraph, "xdot", NULL);

  // read the node position information
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    const char* pos = agget(_graphvizNode[n], "pos");
    double x, y;
    if (sscanf(pos, "%lf,%lf", &x, &y) != 2)
    {
      std::cerr << "Unable to parse position information '" 
                << pos << "' of node '" << _label[n] << "'" << std::endl;
      res &= false;
    }
    else
    {
      _coordX[n] = x;
      _coordY[n] = y;
    }
  }

  // clean up the mess
  gvFreeLayout(pGvc, pGraph);
  agclose(pGraph);
  gvFreeContext(pGvc);

  return res;
}

#endif /* MOLECULEGRAPHVIZ_H_ */
