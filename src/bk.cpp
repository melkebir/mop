/*
 * bk.cpp
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include "verbose.h"
#include "bronkerbosch.h"

typedef lemon::ListGraph Graph;
typedef nina::BronKerbosch<Graph> BronKerbosch;

void buildGraph(Graph& g)
{
  // build our own graph
  Graph::Node a = g.addNode();
  Graph::Node b = g.addNode();
  Graph::Node c = g.addNode();
  Graph::Node d = g.addNode();
  Graph::Node e = g.addNode();
  Graph::Node f = g.addNode();
  Graph::Node h = g.addNode();

  g.addEdge(a, b);
  g.addEdge(a, c);
  g.addEdge(b, c);
  g.addEdge(b, d);
  g.addEdge(c, d);
  g.addEdge(d, e);
  g.addEdge(e, f);
  g.addEdge(e, h);
}

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);
  Graph g;

  int verbosityLevel = static_cast<int>(nina::VERBOSE_ESSENTIAL);
  int bkType = 0;

  ap.refOption("a", "Specifies Bron-Kerbosch algorithm type", bkType, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false);
  ap.parse();
  nina::g_verbosity = static_cast<nina::VerbosityLevel>(verbosityLevel);

  if (ap.files().size() == 0) {
    buildGraph(g);
  }
  else
  {
    // read in from file
    std::ifstream inFile(ap.files()[0].c_str());
    if (inFile.good())
    {
      lemon::graphReader(g, inFile).run();
    }
    else
    {
      std::cerr << "Failed to open " << argv[1] << std::endl;
      return 1;
    }
  }

  BronKerbosch bk(g);
  bk.run(static_cast<BronKerbosch::SolverType>(bkType));
  bk.print(std::cout);

  return 0;
}
