/*
 * atb2lgf.cpp
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

typedef lemon::ListGraph Graph;
using namespace nina::cgp;
using namespace nina;

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);

  int verbosityLevel = static_cast<int>(VERBOSE_ESSENTIAL);
  std::string mtbFile;
  std::string pdbFile;

  ap
    .refOption("mtb", "Specifies the MTB file to use", mtbFile, false)
    .refOption("pdb", "Specifies the PDB file to use for spatial coordinates",
                    pdbFile, false)
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

  pMolecule->writeLGF(std::cout);
  delete pMolecule;
}
