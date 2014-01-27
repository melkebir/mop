/* 
 * readsolution.cpp
 *
 *  Created on: 20-oct-2011
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include <lemon/smart_graph.h> 
#include <lemon/list_graph.h> 
#include <limits>
#include "molecule.h"
#include "cgp/solverdummy.h"

typedef lemon::SmartGraph Graph;

using namespace nina::cgp;
using namespace nina;

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);

  int verbosityLevel = static_cast<int>(VERBOSE_ESSENTIAL);
  int residualError = 0;
  int legend = 1;

  ap.refOption("r", "Specifies how to deal with the residual error per charge group:\n"
                    "     0 - Ignore\n"
                    "     1 - Add to atom with highest error\n"
                    "     2 - Divide equally over all atoms in the charge group\n",
                    residualError, false)
    .refOption("l", "Specifies whether or not the dot file will contain a legend:\n"
                      "     0 - without legend\n"
                      "     1 - with legend\n",
                      legend, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false);
  ap.parse();
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  // initialize molecule and solver...
  Molecule<Graph>* pMolecule = new Molecule<Graph>();

  if (ap.files().size() > 0)
  {
    std::ifstream lgf(ap.files()[0].c_str());
    if (!lgf.good())
    {
      std::cerr << "Failed to open " << ap.files()[0] << std::endl;
      return 1;
    }
    pMolecule->readLGF(lgf);
  }
  else
  {
    pMolecule->readLGF(std::cin);
  }

  Solver<Graph>* pSolver = new SolverDummy<Graph>(pMolecule);
  pSolver->solve(residualError);

  // print solution
  //pSolver->printChargeGroups(std::cerr);
  pSolver->printSolution(std::cerr);
  if (legend)
  {
    pMolecule->printDOT(pSolver, std::cout);
  }
  else
  {
    pMolecule->printDOTnoLegend(pSolver, std::cout);
  }

  delete pSolver;
  delete pMolecule;

  return 0;
}
