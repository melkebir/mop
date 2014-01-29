/*
 * fragments.cpp
 *
 *  Created on: 20-jan-2014
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include "molecule.h"
#include "fragments/product.h"
#include "fragments/bronkerboschconnected.h"
#include "verbose.h"

typedef lemon::ListGraph Graph;

using namespace nina::cgp;
using namespace nina;

typedef Molecule<Graph> MoleculeType;
typedef Product<Graph> ProductType;

bool readLGF(const std::string& filename, MoleculeType& mol)
{
  std::ifstream lgfIn(filename.c_str());
  if (lgfIn.good())
  {
    mol.readLGF(lgfIn);
    return true;
  }
  else
  {
    std::cerr << "Failed to open " << filename << std::endl;
    return false;
  }
}

void output(const ProductType& prod,
            const std::vector< std::vector<Graph::Node> >& cliques,
            std::ostream& out)
{
  for (size_t i = 0; i < cliques.size(); ++i)
  {
    prod.printProductNodeVector(cliques[i], out);
  }
}

void outputJSON(const ProductType& prod,
                const std::vector< std::vector<Graph::Node> >& cliques,
                const std::string& filename,
                std::ostream& out)
{
  out << "    {" << std::endl
      << "      atb_id: " << filename << "," << std::endl
      << "      fragments: [" << std::endl
      << "        {" << std::endl;

  bool first = true;
  for (size_t i = 0; i < cliques.size(); ++i)
  {
    prod.printProductNodeVectorJSON(cliques[i], out, first);
    if (first)
    {
      first = false;
    }
  }

  out << std::endl << "        }" << std::endl
      << "      ]" << std::endl
      << "    }";
}

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);

  int verbosityLevel = static_cast<int>(VERBOSE_ESSENTIAL);
  int shell = 1;
  int bkType = 0;
  bool noJSON = false;
  std::string atb_id;

  ap.refOption("atb_id", "Specifies atb id (used for output)", atb_id, false)
    .refOption("a", "Specifies Bron-Kerbosch algorithm type", bkType, false)
    .refOption("s", "Specifies shell size (default: 1)", shell, false)
    .refOption("no-json", "No JSON output", noJSON, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false);
  ap.parse();
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  if (ap.files().size() != 2)
  {
    std::cerr << "Expected two LGF files as input" << std::endl;
    return 1;
  }

  if (!ap.given("atb_id"))
  {
    atb_id = ap.files()[1];
  }

  MoleculeType mol1, mol2;
  if (readLGF(ap.files()[0], mol1) && readLGF(ap.files()[1], mol2))
  {
    ProductType prod(mol1, mol2, shell);
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
      std::cerr << "Molecule 1 has " << mol1.getNumberOfAtoms()
                << " atoms and " << mol1.getNumberOfBonds()
                << " bonds" << std::endl;
      std::cerr << "Molecule 2 has " << mol2.getNumberOfAtoms()
                << " atoms and " << mol2.getNumberOfBonds()
                << " bonds" << std::endl;
      std::cerr << "Product graph has " << prod.getNumNodes()
                << " nodes and " << prod.getNumEdges()
                << " edges" << std::endl;
    }
    if (g_verbosity >= VERBOSE_DEBUG)
    {
      prod.printDOT(std::cout);
    }

    BronKerboschConnected<Graph> bk(prod);
    lemon::Timer t;
    bk.run(static_cast<BronKerbosch<Graph>::SolverType>(bkType));
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
      std::cerr << "Time: " << t.realTime() << "s" << std::endl;
      std::cerr << "#max-cliques: " << bk.getNumberOfMaxCliques() << std::endl;
    }

    const std::vector< std::vector<Graph::Node> >& cliques = bk.getMaxCliques();
    if (noJSON)
    {
      output(prod, cliques, std::cerr);
    }
    else
    {
      outputJSON(prod, cliques, atb_id, std::cout);
    }
  }
  else
  {
    return 1;
  }

  return 0;
}
