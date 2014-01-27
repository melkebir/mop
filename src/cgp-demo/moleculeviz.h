/* moleculeviz.h
 *
 *  Created on: 30-Jul-2012
 *      Author: M. El-Kebir
 */

#ifndef MOLECULEVIZ_H_
#define MOLECULEVIZ_H_

#include "molecule.h"
#include <ostream>

template<typename GR>
class MoleculeViz : public nina::cgp::Molecule<GR>
{
public:
  /// The graph type of the input graph
  typedef GR Graph;
  /// Base class type
  typedef nina::cgp::Molecule<GR> Parent;

  using Parent::_g;
  using Parent::_partialCharge;
  using Parent::_label;
  using Parent::_label2;
  using Parent::_nAtoms;
  using Parent::_labelToNode;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  size_t _maxChargeGroupSize;
  double _radius;
  DoubleNodeMap _coordX;
  DoubleNodeMap _coordY;

public:
  MoleculeViz(size_t k)
    : Parent()
    , _maxChargeGroupSize(k)
    , _radius()
    , _coordX(_g)
    , _coordY(_g)
  {
  }

  virtual void readLGF(std::istream& in);
  virtual void writeLGF(std::ostream& out);

  double getCoordX(Node n) const { return _coordX[n]; }
  double getCoordY(Node n) const { return _coordY[n]; }
  void setCoordX(Node n, double x) { _coordX[n] = x; }
  void setCoordY(Node n, double y) { _coordY[n] = y; }

  size_t getMaxChargeGroupSize() const { return _maxChargeGroupSize; }
  void setMaxChargeGroupSize(size_t k) { _maxChargeGroupSize = k; }
};

template<typename GR>
void MoleculeViz<GR>::readLGF(std::istream& in)
{
  assert(in.good());
  _g.clear();

  lemon::graphReader(_g, in)
    .nodeMap("partial_charge", _partialCharge)
    .nodeMap("label", _label)
    .nodeMap("label2", _label2)
    .nodeMap("coordX", _coordX)
    .nodeMap("coordY", _coordY)
    .attribute("radius", _radius)
    .run();

  _nAtoms = 0;
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    _labelToNode[_label[n]] = n;
    _nAtoms++;
  }
}

template<typename GR>
void MoleculeViz<GR>::writeLGF(std::ostream& out)
{
  assert(out.good());
  lemon::graphWriter(_g, out)
    .nodeMap("partial_charge", _partialCharge)
    .nodeMap("label", _label)
    .nodeMap("label2", _label2)
    .nodeMap("coordX", _coordX)
    .nodeMap("coordY", _coordY)
    .attribute("radius", _radius)
    .run();
}

#endif /* MOLECULEVIZ_H_ */
