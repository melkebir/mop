/*
 * partition.cpp
 *
 *  Created on: 14-sep-2011
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include <lemon/list_graph.h>
#include <limits>
#include "molecule.h"
#include "symmetrization.h"
#include "symmetrizationautomorphism.h"

typedef lemon::ListGraph Graph;
using namespace nina::cgp;
using namespace nina;

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);

  int verbosityLevel = static_cast<int>(VERBOSE_ESSENTIAL);
  double symMaxAbsDiff = 0.125;
  int n = 4;
  std::string mtbFile;
  std::string pdbFile;
  bool deg1 = false;
  bool path = false;

  ap
    .refOption("mtb", "Specifies the MTB file to use", mtbFile, false)
    .refOption("pdb", "Specifies the PDB file to use for spatial coordinates",
                    pdbFile, false)
    .refOption("symMaxDiff", "Maximum absolute charge difference allowed for atoms\n"
                    "     to be considered chemically equivalent (default: 0.125)",
                    symMaxAbsDiff, false)
    .refOption("n", "Neighborhood size taken into account for chemical equivalence\n"
                    "     (default: 4)",
                    n, false)
    .refOption("path", "Use path heuristic (default: false)", path, false)
    .refOption("deg1", "Set degree 1 neighbors to co-occur with their adjacent node (default: false)", deg1, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false);
  ap.parse();

  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  // parse molecule
  Molecule<Graph>* pMolecule = new Molecule<Graph>();

  if (!mtbFile.empty())
  {
    std::ifstream mtbIn(mtbFile.c_str());
    if (!mtbIn.good())
    {
      std::cerr << "Failed to open " << mtbFile << std::endl;
      return 1;
    }

    pMolecule->readMTB(mtbIn);

    if (!pdbFile.empty())
    {
      std::ifstream pdbIn(pdbFile.c_str());
      if (!pdbIn.good())
      {
        std::cerr << "Failed to open " << pdbFile << std::endl;
        return 1;
      }
      pMolecule->readPDB(pdbIn);
    }
  }
  else
  {
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
  }

  Symmetrization<Graph>* pSym = NULL;
  if (path)
    pSym = new Symmetrization<Graph>(n, symMaxAbsDiff);
  else
    pSym = new SymmetrizationAutomorphism<Graph>(n, symMaxAbsDiff);

  pSym->symmetrize(*pMolecule, pMolecule->getPartialChargeMap());
  if (deg1) pSym->identifyDegree1(*pMolecule);

  delete pSym;

  pMolecule->writeSymLGF(std::cout);
}
