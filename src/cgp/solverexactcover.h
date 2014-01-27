/* 
 * solverexactcover.h
 *
 *  Created on: 21-sep-2011
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_EXACT_COVER_H_
#define SOLVER_EXACT_COVER_H_

#include <ilcplex/ilocplex.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include "solver.h"
#include "cgp/extmolecule.h"

template<typename GR>
class SolverExactCover : public Solver<GR>
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
  /// Charge group set type
  typedef typename ExtMolecule<Graph>::ChargeGroupSet ChargeGroupSet;
  /// Charge group set iterator type
  typedef typename ExtMolecule<Graph>::ChargeGroupSetIt ChargeGroupSetIt;

private:
  const ExtMolecule<Graph>* _pExtMolecule;

public:
  SolverExactCover(const ExtMolecule<Graph>* pExtMolecule)
    : Parent(pExtMolecule)
    , _pExtMolecule(pExtMolecule)
  {
  };

  void solve(int resError);
};

template<typename GR>
void SolverExactCover<GR>::solve(int resError)
{
  const size_t nChargeGroups = _pExtMolecule->getNumberOfChargeGroups();
  const size_t maxChargeGroupSize = _pExtMolecule->getMaxChargeGroupSize();
  const Graph& g = _pExtMolecule->getGraph();

  IloEnv env;
  IloModel model(env);
  IloBoolVarArray x(env, nChargeGroups);
  //IloFloatVar y(env);
  IloExpr obj(env);
  IloCplex* pCplex = NULL;

  // generate model here
  IloExpr sum(env);
  for (NodeIt n(g); n != lemon::INVALID; ++n)
  {
    Node node = n;

    size_t j = 0;
    for (size_t i = 0; i < maxChargeGroupSize; i++)
    {
      const ChargeGroupSet& cgSet = _pExtMolecule->getChargeGroup(i+1);
      for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++, j++)
      {
        assert(j < nChargeGroups);
        if (cgIt->find(node) != cgIt->end())
        {
          sum += x[j];
        }
      }
    }

    model.add(sum == 1);
    sum.clear();
  }
  sum.end();

  // objective function
  size_t j = 0;
  for (size_t i = 0; i < maxChargeGroupSize; i++)
  {
    const ChargeGroupSet& cgSet = _pExtMolecule->getChargeGroup(i+1);
    for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++, j++)
    {
      double err = _pExtMolecule->computeError(*cgIt);
      obj += x[j] * err;
    }
  }
  model.add(IloMinimize(env, obj));

  //size_t j = 0;
  //for (size_t i = 0; i < maxChargeGroupSize; i++)
  //{
  //  const ChargeGroupSet& cgSet = _pExtMolecule->getChargeGroup(i+1);
  //  for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++, j++)
  //  {
  //    double err = _pExtMolecule->computeError(*cgIt);
  //    model.add(y >= x[j] * err);
  //    //obj += x[j] * err;
  //  }
  //}
  //model.add(IloMinimize(env, y));

  // UGLY hack, just to suppress the ILOG License Manager spamming
  int bakfd, newfd;
  fflush(stderr);
  bakfd = dup(2);
  newfd = open("/dev/null", O_WRONLY);
  dup2(newfd, 2);
  close(newfd);

  pCplex = new IloCplex(model);
  if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
  {
    pCplex->setOut(std::cerr);
    pCplex->setWarning(std::cerr);
    pCplex->setError(std::cerr);
  }
  else
  {
    pCplex->setOut(env.getNullStream());
    pCplex->setWarning(env.getNullStream());
    pCplex->setError(env.getNullStream());
  }

  fflush(stderr);
  dup2(bakfd, 2);
  close(bakfd);

  if (g_verbosity > VERBOSE_NONE)
  {
    std::cerr << "Solving exact cover ILP..." << std::endl;
  }

  if (pCplex->solve())
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
      std::cerr << "Obj value: " << pCplex->getObjValue() << std::endl;

    // get solution
    _chargeGroupVector.clear();
    j = 0;
    int color = 0;
    for (size_t i = 0; i < maxChargeGroupSize; i++)
    {
      const ChargeGroupSet& cgSet = _pExtMolecule->getChargeGroup(i+1);
      for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++, j++)
      {
        if (g_verbosity > VERBOSE_NONE)
        {
          std::cerr << "\rReconstructing solution: " 
            << j+1 << "/" << nChargeGroups 
            << std::flush;
        }
        if (pCplex->getValue(x[j]))
        {
          _chargeGroupVector.push_back(*cgIt);
          // assigns all nodes in *cgIt the same color
          for (ChargeGroupIt nodeIt = cgIt->begin(); 
              nodeIt != cgIt->end(); 
              ++nodeIt)
          {
            _color[*nodeIt] = color;
          }
          color++;
        }
      }
    }
    if (g_verbosity > VERBOSE_NONE)
    {
      std::cerr << std::endl;
    }

    updateCorrectedPartialCharges(resError);
  }
  else
  {
    std::cerr << "Infeasible!" << std::endl;
  }

  delete pCplex;
  env.end();
}

#endif /* SOLVER_EXACT_COVER_H_ */
