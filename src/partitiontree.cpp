/* 
 * partitiontree.cpp
 *
 *  Created on: 26-feb-2013
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include <lemon/smart_graph.h> 
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <limits>
#include "molecule.h"
#include "extmolecule.h"
#include "solverdp.h"

typedef lemon::ListGraph Graph;

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);

  int verbosityLevel = static_cast<int>(VERBOSE_ESSENTIAL);
  int maxPartSize;
  int residualError = 0;
  double diameter = 2;
  std::string dotFile;

  ap.refOption("k", "Maximum number of atoms in a charge group", maxPartSize, false)
    .refOption("d", "Maximum diameter", diameter, false)
    .refOption("dot", "Specifies the dot output file", dotFile, false)
    .refOption("r", "Specifies how to deal with the residual error per charge group:\n"
                    "     0 - Ignore (default)\n"
                    "     1 - Add to atom with highest error\n"
                    "     2 - Divide equally over all atoms in the charge group\n"
                    "         with the same sign\n"
                    "     3 - Divide equally over all atoms in the charge group",
                    residualError, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false);
  ap.parse();
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  // parse molecule
  ExtMolecule<Graph>* pMolecule = new ExtMolecule<Graph>();

  // let's read some lgf
  if (ap.files().size() > 0)
  {
    std::ifstream lgfIn(ap.files()[0].c_str());
    if (!lgfIn.good())
    {
      std::cerr << "Failed to open " << ap.files()[0] << std::endl;
      return 1;
    }
    pMolecule->readLGF(lgfIn);
  }
  else
  {
    pMolecule->readLGF(std::cin);
  }

  // initialize molecule based on diameter or k
  if (!ap.given("k") && !ap.given("d"))
  {
    std::cerr << "Either k (maximum charge group size) or d (diameter) must be specified"
              << std::endl;
    return 1;
  }

  if (ap.given("k"))
  {
    if (maxPartSize < 1)
    {
      std::cerr << "Value of parameter -k needs to be at least 1" << std::endl;
      return 1;
    }
    pMolecule->initGlobalK(maxPartSize);
  }

  if (ap.given("d"))
  {
    if (diameter < 0)
    {
      std::cerr << "Value of parameter -d needs to be nonnegative" << std::endl;
      return 1;
    }
    pMolecule->initLocalK(diameter);
  }

  // initialize solver
  Solver<Graph>* pSolver = new SolverDP<Graph>(pMolecule);

  // solve
  lemon::Timer t;
  bool solved = pSolver->solve(residualError);
  if (solved && g_verbosity >= VERBOSE_ESSENTIAL)
    std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  // print solution
  if (solved)
  {
    if (!dotFile.empty() && dotFile != "-")
    {
      std::ofstream outDot(dotFile.c_str());
      if (outDot.good())
      {
        pSolver->printSolution(std::cout);
        pMolecule->printDOT(pSolver, outDot);
      }
    }
    else if (dotFile == "-")
    {
      pSolver->printSolution(std::cerr);
      pMolecule->printDOT(pSolver, std::cout);
    }
    else
    {
      pSolver->printSolution(std::cout);
    }
  }
  else
  {
    std::cerr << "Unsolved, infeasible problem instance!" << std::endl;
  }

  delete pSolver;
  delete pMolecule;

  return 0;
}
