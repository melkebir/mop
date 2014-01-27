/*
 * solver.h
 *
 *  Created on: 21-sep-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <assert.h>
#include <stdlib.h>
#include <limits>
#include <math.h>
#include <set>
#include <lemon/core.h>
#include <iomanip>
#include "verbose.h"
#include "molecule.h"

namespace nina {
namespace cgp {

template<typename GR>
class Solver
{
public:
  /// The graph type of the input graph
  typedef GR Graph;

private:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public:
  /// Charge group type
  typedef typename Molecule<Graph>::ChargeGroup ChargeGroup;
  /// Charge group iterator
  typedef typename Molecule<Graph>::ChargeGroupIt ChargeGroupIt;
  /// Color map type
  typedef typename Molecule<Graph>::ColorMap ColorMap;
  /// Charge group vector type
  typedef typename Molecule<Graph>::ChargeGroupVector ChargeGroupVector;
  /// Charge group vector iterator
  typedef typename Molecule<Graph>::ChargeGroupVectorIt ChargeGroupVectorIt;
  /// Weights on nodes
  typedef typename Molecule<Graph>::WeightNodeMap WeightNodeMap;

protected:
  /// The molecule
  const Molecule<Graph>* _pMolecule;
  /// Charge groups
  ChargeGroupVector _chargeGroupVector;
  /// Color map
  ColorMap _color;
  /// Corrected partial charges map
  WeightNodeMap _correctedPartialCharge;

public:
  Solver(const Molecule<Graph>* pMolecule)
    : _pMolecule(pMolecule)
    , _chargeGroupVector()
    , _color(pMolecule->getGraph(), -1)
    , _correctedPartialCharge(pMolecule->getGraph())
  {
  }

  virtual ~Solver() {}
  
  virtual bool solve(int resError) = 0;

  void updateCorrectedPartialCharges(int resError);

  const ChargeGroupVector& getChargeGroupVector() const { return _chargeGroupVector; }

  const ColorMap& getColorMap() const { return _color; }

  int getColor(Node n) const { return _color[n]; }

  const WeightNodeMap& getCorrectedPartialChargeMap() const
  {
    return _correctedPartialCharge;
  }

  double getCorrectedPartialCharge(Node n) const
  {
    return _correctedPartialCharge[n];
  }

  void printChargeGroups(std::ostream& out) const
  {
    for (ChargeGroupVectorIt it = _chargeGroupVector.begin(); 
        it != _chargeGroupVector.end(); it++)
    {
      _pMolecule->printChargeGroup(*it, out);
    }
  }

  void getTotalCharge(const ChargeGroup& cg,std::ostream& out) const
  {
    _pMolecule->printCGTotalCharge(cg, out);

  }

  void printSolution(std::ostream& out) const
  {
    const Graph& g = _pMolecule->getGraph();

    for (NodeIt n(g); n != lemon::INVALID; ++n)
    {
      out << _pMolecule->getLabel(n) << "\t"
        << std::setiosflags(std::ios::fixed)
        << std::setprecision(6) << getCorrectedPartialCharge(n) << "\t"
        << getColor(n) << std::endl;
    }
  }

  double getError() const
  {
    double err = 0;
    for (ChargeGroupVectorIt cgIt = _chargeGroupVector.begin(); 
        cgIt != _chargeGroupVector.end(); cgIt++)
    {
      err += _pMolecule->computeError(*cgIt, true, false);
    }

    return err;
  }
};
} // namespace cgp
} // namespace nina

template<typename GR>
inline void nina::cgp::Solver<GR>::updateCorrectedPartialCharges(int resError)
{
  const Graph& g = _pMolecule->getGraph();

  for (NodeIt n(g); n != lemon::INVALID; ++n)
    _correctedPartialCharge[n] = _pMolecule->getPartialCharge(n);

  if (resError == 1 || resError == 2 || resError == 3)
  {
    for (ChargeGroupVectorIt cgIt = _chargeGroupVector.begin(); 
        cgIt != _chargeGroupVector.end(); cgIt++)
    {
      double err = _pMolecule->computeError(*cgIt, false);

      if (err == 0)
        continue;

      ChargeGroup nodes;
      if (resError == 3)
      {
        // in case resError == 3, residual error needs to be distributed over all atoms
        nodes = *cgIt;
      }
      if (resError == 2)
      {
        // in this case we distribute the residual error
        // over all atoms whose partial charges have the same sign
        for (ChargeGroupIt nodeIt = cgIt->begin(); 
            nodeIt != cgIt->end(); nodeIt++)
        {
          if (err < 0 && _pMolecule->getPartialCharge(*nodeIt) < 0)
            nodes.insert(*nodeIt);
          else if (err > 0 && _pMolecule->getPartialCharge(*nodeIt) > 0)
            nodes.insert(*nodeIt);
        }
      }
      if (resError == 1 || nodes.size() == 0)
      {
        // in this case we add the residual error 
        // to the atom whose partial charge is maximal/minimal
        double extremum = err > 0 ? 
          -1 * std::numeric_limits<double>::max() : std::numeric_limits<double>::max();
        Node n = lemon::INVALID;

        for (ChargeGroupIt nodeIt = cgIt->begin(); 
            nodeIt != cgIt->end(); nodeIt++)
        {
          if (err < 0 && _pMolecule->getPartialCharge(*nodeIt) < extremum)
          {
            n = *nodeIt;
            extremum = _pMolecule->getPartialCharge(n);
          }
          else if (err > 0 && _pMolecule->getPartialCharge(*nodeIt) > extremum)
          {
            n = *nodeIt;
            extremum = _pMolecule->getPartialCharge(n);
          }
        }

        if (n != lemon::INVALID)
          nodes.insert(n);
      }

      for (ChargeGroupIt nodeIt = nodes.begin(); 
          nodeIt != nodes.end(); nodeIt++)
      {
        //std::cerr << "Node " << _pMolecule->getLabel(*nodeIt) << ", " << _pMolecule->getLabel2(*nodeIt) 
        //  << ": " << _correctedPartialCharge[*nodeIt]
        //  << std::flush;
        _correctedPartialCharge[*nodeIt] += err / nodes.size();
        //std::cerr << ", " << _correctedPartialCharge[*nodeIt] << std::endl;
      }
    }
  }
}

#endif /* SOLVER_H_ */
