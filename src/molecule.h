/*
 * molecule.h
 *
 *  Created on: 14-sep-2011
 *      Author: M. El-Kebir
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <list>
#include <set>
#include <string>
#include <map>
#include <limits>
#include <math.h>
#include <assert.h>
#include <lemon/core.h>
#include <lemon/lgf_reader.h>
#include "verbose.h"

namespace nina {
namespace cgp {

// forward class declaration
template<typename GR>
class Solver;

template<typename GR>
class Molecule
{
public:
  /// The graph type of the input graph
  typedef GR Graph;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public: 
  /// Weights on nodes
  typedef typename Graph::template NodeMap<double> WeightNodeMap;
  /// Labels on nodes
  typedef typename Graph::template NodeMap<std::string> LabelNodeMap;
  /// Charge group type
  typedef typename std::set<Node> ChargeGroup;
  /// Charge group iterator
  typedef typename std::set<Node>::const_iterator ChargeGroupIt;
  /// Charge group vector type
  typedef typename std::vector<ChargeGroup> ChargeGroupVector;
  /// Charge group vector iterator
  typedef typename ChargeGroupVector::const_iterator ChargeGroupVectorIt;
  /// Charge group set type
  typedef typename std::set<ChargeGroup> ChargeGroupSet;
  /// Charge group set iterator
  typedef typename ChargeGroupSet::const_iterator ChargeGroupSetIt;
  /// Charge group non const set iterator
  typedef typename ChargeGroupSet::iterator ChargeGroupSetNonConstIt;
  /// Vector of charge group sets type
  typedef typename std::vector<ChargeGroupSet> ChargeGroupSetVector;
  /// Vector of charge group sets iterator type
  typedef typename ChargeGroupSetVector::const_iterator ChargeGroupSetVectorIt;
  /// Vector of charge group sets iterator type
  typedef typename ChargeGroupSetVector::iterator ChargeGroupSetVectorNonConstIt;
  /// Color map type
  typedef typename Graph::template NodeMap<int> ColorMap;
  /// Cardinality map type
  typedef typename Graph::template NodeMap<int> CardMap;
  /// String to node map type
  typedef typename std::map<std::string, Node> StringToNodeMap;
  /// String to node map iterator
  typedef typename StringToNodeMap::const_iterator StringToNodeMapIt;
  /// Node vector type
  typedef typename std::vector<Node> NodeVector;
  /// Node vector iterator type
  typedef typename NodeVector::const_iterator NodeVectorIt;
  /// Node matrix type
  typedef typename std::vector<NodeVector> NodeMatrix;
  /// Node matrix iterator type
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  /// Node matrix non const iterator type
  typedef typename NodeMatrix::iterator NodeMatrixNonConstIt;
  /// Node set type
  typedef typename std::set<Node> NodeSet;
  /// Node set iterator type
  typedef typename NodeSet::const_iterator NodeSetIt;
  /// Node set vector type
  typedef typename std::vector<NodeSet> NodeSetVector;
  /// Node set vector iterator type
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
  /// Node set vector non const iterator type
  typedef typename NodeSetVector::iterator NodeSetVectorNonConstIt;
  /// Node set map type
  typedef typename Graph::template NodeMap<NodeSet> NodeSetMap;
  /// Atom type list type
  typedef typename std::list<int> AtomList;
  /// Atom type list iterator type
  typedef typename AtomList::const_iterator AtomListIt;
  /// Atom type list set type
  typedef typename std::multiset<AtomList> AtomListMultiSet;
  /// Atom type list set iterator type
  typedef typename AtomListMultiSet::const_iterator AtomListSetMultiIt;
  /// Atom list set map type
  typedef typename Graph::template NodeMap<AtomListMultiSet> AtomListMultiSetMap;

  struct Coord
  {
  private:
    double _x;
    double _y;
    double _z;

  public:
    Coord(double x, double y, double z)
      : _x(x)/*
   * fragments.h
   *
   *  Created on: 20-jan-2014
   *      Author: M. El-Kebir
   */
      , _y(y)
      , _z(z)
    {
    }

    Coord()
      : _x(0)
      , _y(0)
      , _z(0)
    {
    }

    double dist(const Coord& b) const
    {
      double dx = _x - b._x;
      double dy = _y - b._y;
      double dz = _z - b._z;

      return sqrt(dx*dx + dy*dy + dz*dz);
    }

    double getX() const { return _x; }
    double getY() const { return _y; }
    double getZ() const { return _z; }
  };

  typedef typename Graph::template NodeMap<Coord> CoordNodeMap;

protected:/*
 * fragments.h
 *
 *  Created on: 20-jan-2014
 *      Author: M. El-Kebir
 */
  /// The molecule
  Graph _g;
  /// Partial charges map
  WeightNodeMap _partialCharge;
  /// Label map
  LabelNodeMap _label;
  /// Label map
  LabelNodeMap _label2;
  /// Atom type map
  IntNodeMap _atomType;
  /// Coordinates
  CoordNodeMap _coord;
  /// Number of atoms
  size_t _nAtoms;
  /// Map from label to node
  StringToNodeMap _labelToNode;
  /// Local maximum charge group size
  IntNodeMap _k;
  /// Diameter
  double _d;
  /// Global maximum charge group size
  size_t _globalK;
  /// Node set map, containing nodes that must occur together in a charge group
  NodeSetMap _sameCg;
  /// Node matrix contining equivalence classes
  NodeMatrix _eqClasses;

  bool generate(ChargeGroupSetVector& chargeGroups) const;

  void initLocalK(ChargeGroupSetVector& cgSetVec);

public:
  explicit Molecule()
    : _g()
    , _partialCharge(_g)
    , _label(_g)
    , _label2(_g)
    , _atomType(_g)
    , _coord(_g)
    , _nAtoms()
    , _labelToNode()
    , _k(_g)
    , _d(std::numeric_limits<double>::max())
    , _globalK(std::numeric_limits<size_t>::max())
    , _sameCg(_g)
  {
  }

  virtual ~Molecule() {}

  void setEqClasses(const NodeMatrix& eqClasses)
  {
    _eqClasses = eqClasses;
  }

  void setSameChargeGroup(const NodeSet& nodeSet)
  {
    for (NodeSetIt nodeIt = nodeSet.begin(); nodeIt != nodeSet.end(); nodeIt++)
    {
      _sameCg[*nodeIt] = nodeSet;
    }
  }

  const NodeSet& getSameChargeGroup(Node node) const
  {
    assert(node != lemon::INVALID);
    return _sameCg[node];
  }

  virtual void initLocalK(double diameter);

  virtual void initGlobalK(size_t k);

  size_t getMaxDegree() const;

  size_t getMaxChargeGroupSize(Node v) const
  {
    return static_cast<size_t>(_k[v]);
  }

  const Node getNodeByLabel(const std::string& label) const
  {
    StringToNodeMapIt it = _labelToNode.find(label);
    if (it != _labelToNode.end())
      return it->second;
    else
      return lemon::INVALID;
  }

  const LabelNodeMap& getLabelMap() const
  {
    return _label;
  }

  const std::string& getLabel(Node n) const
  {
    return _label[n];
  }

  const LabelNodeMap& getLabel2Map() const
  {
    return _label2;
  }

  const std::string& getLabel2(Node n) const
  {
    return _label2[n];
  }

  const Graph& getGraph() const
  {
    return _g;
  }

  int getAtomType(Node n) const
  {
    return _atomType[n];
  }

  WeightNodeMap& getPartialChargeMap()
  {
    return _partialCharge;
  }

  const WeightNodeMap& getPartialChargeMap() const
  {
    return _partialCharge;
  }

  double getPartialCharge(Node n) const
  {
    return _partialCharge[n];
  }

  double getPartialCharge(const ChargeGroup& cg) const
  {
    double res = 0;
    for (ChargeGroupIt nodeIt = cg.begin(); nodeIt != cg.end(); nodeIt++)
    {
      res += _partialCharge[*nodeIt];
    }
    return res;
  }

  size_t getNumberOfAtoms() const
  {
    return _nAtoms;
  }

  size_t getNumberOfBonds() const
  {
    return static_cast<size_t>(lemon::countEdges(_g));
  }

  virtual bool readMTB(std::istream& in);

  virtual void readLGF(std::istream& in);

  virtual bool readPDB(std::istream& in);

  virtual void writeLGF(std::ostream& out) const;

  virtual void writeSymLGF(std::ostream& out) const;

  double computeError(const ChargeGroup& cg,
                      bool abs = true,
                      bool eps = true) const;
  void printChargeGroup(const ChargeGroup& cg,
                        std::ostream& out) const;
  void printCGTotalCharge(const ChargeGroup& cg,
                        std::ostream& out) const;
  void printDOT(const Solver<GR>* pSolver,
                std::ostream& out) const;
  void printDOTnoLegend(const Solver<GR>* pSolver,
                std::ostream& out) const;

  double getDiameter(const ChargeGroup& cg) const;

  bool isPotFeasible(const ChargeGroup& cg) const
  {
    return getDiameter(cg) <= _d && cg.size() <= _globalK;
  }

  bool isFeasible(const ChargeGroup& cg) const
  {
    if (isPotFeasible(cg))
    {
      for (NodeSetIt nodeIt = cg.begin(); nodeIt != cg.end(); nodeIt++)
      {
        const NodeSet& sameCg = _sameCg[*nodeIt];
        if (!std::includes(cg.begin(), cg.end(), sameCg.begin(), sameCg.end()))
          return false;
      }
      return true;
    }
    else
    {
      return false;
    }
  }
};

template<typename GR>
inline void Molecule<GR>::initGlobalK(size_t k)
{
  _globalK = k;

  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    _k[v] = static_cast<int>(k);
  }
}

template<typename GR>
inline void Molecule<GR>::initLocalK(ChargeGroupSetVector& cgSetVec)
{
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    // { v }
    ChargeGroup cgInit;
    cgInit.insert(v);

    // add { { v } } to cgSetVec[0]
    cgSetVec[0].insert(cgInit);
  }

  // generate all charge groups
  while (generate(cgSetVec));

  // determine k for each node
  size_t curSize = 1;
  for (ChargeGroupSetVectorIt cgSetIt = cgSetVec.begin();
       cgSetIt != cgSetVec.end(); cgSetIt++, curSize++)
  {
    const ChargeGroupSet& cgSet = *cgSetIt;

    for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++)
    {
      const ChargeGroup& cg = *cgIt;
      for (ChargeGroupIt nodeIt = cg.begin(); nodeIt != cg.end(); nodeIt++)
      {
        assert(cg.size() == curSize);
        _k[*nodeIt] = curSize;
      }
    }
  }

  if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
  {
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
      std::cerr << _label[v] << "\t" << _k[v] << std::endl;
  }
}

template<typename GR>
inline void Molecule<GR>::initLocalK(double diameter)
{
  assert(diameter >= 0);
  _d = diameter;

  ChargeGroupSetVector cgSetVec;
  cgSetVec.push_back(ChargeGroupSet());

  initLocalK(cgSetVec);
}

template<typename GR>
inline bool Molecule<GR>::generate(ChargeGroupSetVector& chargeGroups) const
{
  assert(chargeGroups.size() > 0);

  bool result = false;

  const size_t curSize = chargeGroups.size();
  if (curSize == _globalK)
    return false;

  chargeGroups.push_back(ChargeGroupSet());
  ChargeGroupSet& newCgSet = chargeGroups.back();
  const ChargeGroupSet& cgSet = chargeGroups[chargeGroups.size() - 2];
  for (ChargeGroupSetIt cgIt = cgSet.begin(); cgIt != cgSet.end(); cgIt++)
  {
    for (ChargeGroupIt nodeIt = cgIt->begin(); nodeIt != cgIt->end(); nodeIt++)
    {
      // takes d |V| time
      const Node u = *nodeIt;

      for (IncEdgeIt e(_g, u); e != lemon::INVALID; ++e)
      {
        Node v = _g.runningNode(e);

        // skip if u is already in cg
        if (cgIt->find(v) != cgIt->end()) continue;

        bool diameterOK = true;
        for (ChargeGroupIt nodeIt2 = cgIt->begin(); nodeIt2 != cgIt->end(); nodeIt2++)
        {
          if (_coord[v].dist(_coord[*nodeIt2]) > _d)
          {
            diameterOK = false;
            break;
          }
        }

        if (diameterOK)
        {
          assert(cgIt->size() < _globalK);

          ChargeGroup cg = *cgIt;
          cg.insert(v);

          // we found one
          newCgSet.insert(cg);
          result = true;
        }
      }
    }
  }

  return result;
}

template<typename GR>
inline double Molecule<GR>::getDiameter(const ChargeGroup& cg) const
{
  double max = 0;

  for (ChargeGroupIt cgIt1 = cg.begin(); cgIt1 != cg.end(); cgIt1++)
  {
    for (ChargeGroupIt cgIt2 = cgIt1; cgIt2 != cg.end(); cgIt2++)
    {
      double d = _coord[*cgIt1].dist(_coord[*cgIt2]);
      if (d > max) max = d;
    }
  }

  return max;
}

template<typename GR>
inline size_t Molecule<GR>::getMaxDegree() const
{
  size_t maxDegree = 0;
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    size_t degree = 0;
    for (IncEdgeIt e(_g, n); e != lemon::INVALID; ++e)
      degree++;

    if (degree > maxDegree)
      maxDegree = degree;
  }

  return maxDegree;
}


template<typename GR>
inline void Molecule<GR>::printCGTotalCharge(const ChargeGroup& cg,
                                                      std::ostream& out) const
{
  double sum = 0;
  for (ChargeGroupIt it = cg.begin(); it != cg.end(); it++)
  {
    sum += _partialCharge[*it];
  }
  out << std::fixed << std::setprecision( 3 ) << sum;
}

template<typename GR>
inline void Molecule<GR>::printChargeGroup(const ChargeGroup& cg,
                                                      std::ostream& out) const
{
  double sum = 0;
  bool first = true;
  for (ChargeGroupIt it = cg.begin(); it != cg.end(); it++)
  {
    sum += _partialCharge[*it];

    if (first)
    {
      first = false;
      out << getLabel(*it);
      if (getLabel2(*it) != "")
        out << " (" << getLabel2(*it) << ")";
    }
    else
    {
      out << ", " << getLabel(*it);
      if (getLabel2(*it) != "")
        out << " (" << getLabel2(*it) << ")";
    }
  }

  double error = computeError(cg, false);

  if (error > 0)
  {
    out << "\t"
        << ceil(sum) << " - " << sum << " = "
        << computeError(cg, false)
        << " (d= " << getDiameter(cg)
        << ",k=" << cg.size() << ")"
        << std::endl;
  }
  else
  {
    out << "\t"
        << floor(sum) << " - " << sum << " = "
        << computeError(cg, false)
        << " (d= " << getDiameter(cg)
        << ",k=" << cg.size() << ")"
        << std::endl;
  }
}

template<typename GR>
inline double Molecule<GR>::computeError(const ChargeGroup& cg,
                                                    bool abs,
                                                    bool eps) const
{
  double partial_charge = 0;
  for (ChargeGroupIt it = cg.begin(); it != cg.end(); it++)
  {
    partial_charge += _partialCharge[*it];
  }

  double error =
    std::min(ceil(partial_charge) - partial_charge, partial_charge - floor(partial_charge));

  if (abs)
  {
    // epsilon to ensure that as few charge groups as possible are used
    const double epsVal = 0.000001;
    return (eps ? epsVal : 0) + error;
  }
  else
  {
    if (ceil(partial_charge) - partial_charge < partial_charge - floor(partial_charge))
    {
      return + error;
    }
    else
    {
      return -error;
    }
  }
}

template<typename GR>
inline bool Molecule<GR>::readPDB(std::istream& in)
{
  assert(in.good());

  std::string line;
  while (std::getline(in, line))
  {
    if (line.size() >= 6 && line.substr(0, 6) == "HETATM")
    {
      std::stringstream lineStream(line);
      std::string dummy, label, label2;
      double x, y, z;

      lineStream >> dummy >> label >> label2 >> dummy >> dummy >> x >> y >> z;

      if (_labelToNode.find(label) == _labelToNode.end())
      {
        return false;
      }

      _coord[_labelToNode[label]] = Coord(x, y, z);
    }
  }

  // compute pairwise distances, todo: construct look-up table
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    for (NodeIt w(_g); w != lemon::INVALID; ++w)
    {
      if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
      {
        std::cerr << _label[v] << "\t"
                  << _label[w] << "\t"
                  << _coord[v].dist(_coord[w])
                  << std::endl;
      }
    }
  }

  return true;
}

template<typename GR>
inline bool Molecule<GR>::readMTB(std::istream& in)
{
  assert(in.good());
  _g.clear();
  _labelToNode.clear();
  _nAtoms = 0;

  bool readAtoms = false;
  int nAtoms = -1;
  bool readBonds = false;
  int nBonds = -1;

  bool readTitle = false;
  std::string title;

  std::string line;
  while (std::getline(in, line))
  {
    // skip comments
    if (line.empty() || line[0] == '#') continue;

    // let's read on until we find 'MTBUILDBLSOLUTE'
    if (line.substr(0, 15) == "MTBUILDBLSOLUTE")
    {
      readTitle = true;
    }
    else if (readTitle)
    {
      title = line;
      readAtoms = true;
      readTitle = false;
    }
    else if (readAtoms && nAtoms == -1)
    {
      std::stringstream lineStream(line);

      // read the number of atoms
      lineStream >> nAtoms;
    }
    else if (readAtoms && nAtoms > 0)
    {
      std::stringstream lineStream(line);

      // read an atom
      std::string label, label2, dummy;
      int atomType;
      double charge;

      lineStream >> label >> label2 >> atomType >> dummy >> charge;

      nAtoms--;
      _nAtoms++;

      Node atom = _g.addNode();
      _partialCharge[atom] = charge;
      _label[atom] = label;
      _label2[atom] = label2;
      _atomType[atom] = atomType;
      _labelToNode[label] = atom;
    }
    else if (!readBonds && nAtoms == 0)
    {
      std::stringstream lineStream(line);
      lineStream >> nBonds;

      readBonds = true;
    }
    else if (readBonds && nBonds > 0)
    {
      std::stringstream lineStream(line);

      std::string atom1, atom2;
      lineStream >> atom1 >> atom2;

      if (_labelToNode.find(atom1) != _labelToNode.end() &&
          _labelToNode.find(atom2) != _labelToNode.end())
      {
        _g.addEdge(_labelToNode[atom1], _labelToNode[atom2]);
        nBonds--;
      }
      else
      {
        return false;
      }
    }
    else if (nAtoms == 0 && nBonds == 0)
    {
      break;
    }
  }

  if (nAtoms == 0 && nBonds == 0)
  {
    for (NodeIt n(_g); n != lemon::INVALID; ++n)
    {
      _sameCg[n].insert(n);
    }

    return true;
  }
  else
  {
    return false;
  }
}

template<typename GR>
inline void Molecule<GR>::readLGF(std::istream& in)
{
  assert(in.good());
  _g.clear();

  DoubleNodeMap coordX(_g);
  DoubleNodeMap coordY(_g);
  DoubleNodeMap coordZ(_g);
  IntNodeMap initColorIdx(_g);

  lemon::graphReader(_g, in)
    .nodeMap("partial_charge", _partialCharge)
    .nodeMap("label", _label)
    .nodeMap("label2", _label2)
    .nodeMap("atomType", _atomType)
    .nodeMap("coordX", coordX)
    .nodeMap("coordY", coordY)
    .nodeMap("coordZ", coordZ)
    .nodeMap("initColor", initColorIdx)
    .run();

  _nAtoms = 0;
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    _labelToNode[_label[n]] = n;
    _nAtoms++;
    _coord[n] = Coord(coordX[n], coordY[n], coordZ[n]);
  }

  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    NodeSet same;

    for (NodeIt w(_g); w != lemon::INVALID; ++w)
    {
      if (initColorIdx[v] == initColorIdx[w])
      {
        same.insert(w);
      }
    }

    setSameChargeGroup(same);
  }
}

template<typename GR>
inline void Molecule<GR>::writeLGF(std::ostream& out) const
{
  assert(out.good());

  DoubleNodeMap coordX(_g);
  DoubleNodeMap coordY(_g);
  DoubleNodeMap coordZ(_g);

  IntNodeMap initColorIdx(_g, -1);
  int sameIdx = 0;

  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    coordX[n] = _coord[n].getX();
    coordY[n] = _coord[n].getY();
    coordZ[n] = _coord[n].getZ();

    if (initColorIdx[n] == -1)
    {
      const NodeSet& same = _sameCg[n];
      for (NodeSetIt nodeIt = same.begin(); nodeIt != same.end(); nodeIt++)
      {
        initColorIdx[*nodeIt] = sameIdx;
      }
      sameIdx++;
    }
  }

  lemon::graphWriter(_g, out)
    .nodeMap("partial_charge", _partialCharge)
    .nodeMap("label", _label)
    .nodeMap("label2", _label2)
    .nodeMap("atomType", _atomType)
    .nodeMap("coordX", coordX)
    .nodeMap("coordY", coordY)
    .nodeMap("coordZ", coordZ)
    .nodeMap("initColor", initColorIdx)
    .run();
}

template<typename GR>
inline void Molecule<GR>::writeSymLGF(std::ostream& out) const
{
  assert(out.good());

  DoubleNodeMap coordX(_g);
  DoubleNodeMap coordY(_g);
  DoubleNodeMap coordZ(_g);

  IntNodeMap initColorIdx(_g, -1);
  int sameIdx = 0;

  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    coordX[n] = _coord[n].getX();
    coordY[n] = _coord[n].getY();
    coordZ[n] = _coord[n].getZ();

    if (initColorIdx[n] == -1)
    {
      const NodeSet& same = _sameCg[n];
      for (NodeSetIt nodeIt = same.begin(); nodeIt != same.end(); nodeIt++)
      {
        initColorIdx[*nodeIt] = sameIdx;
      }
      sameIdx++;
    }
  }

  IntNodeMap symGrps(_g, -1);
  int symId = 0;

  for (NodeMatrixIt it = _eqClasses.begin(); it != _eqClasses.end(); it++)
  {

    const NodeVector& nodes = *it;

    for (NodeVectorIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
    {
      symGrps[*nodeIt] = symId;
    }
    symId++;
  }

  lemon::graphWriter(_g, out)
    .nodeMap("partial_charge", _partialCharge)
    .nodeMap("label", _label)
    .nodeMap("label2", _label2)
    .nodeMap("atomType", _atomType)
    .nodeMap("coordX", coordX)
    .nodeMap("coordY", coordY)
    .nodeMap("coordZ", coordZ)
    .nodeMap("initColor", initColorIdx)
    .nodeMap("eqClasses", symGrps)
    .run();
}

template<typename GR>
inline void Molecule<GR>::printDOT(const Solver<GR>* pSolver,
                                   std::ostream& out) const
{
  static const int nColors = 11;
  static std::string color[nColors] = 
  { 
    "darkolivegreen1", 
    "red", 
    "cyan", 
    "green", 
    "yellow", 
    "crimson", 
    "gold", 
    "violet",
    "aquamarine",
    "crimson",
    "chartreuse"
  };

  out << "graph G {" << std::endl;
  out << "\tlabel=\"Total error: " << std::fixed << pSolver->getError() 
    //<< ", k=" << getMaxChargeGroupSize()
    << "\"" << std::endl
    << "\tlabelloc=t" << std::endl;
  out << "\tnode [style=filled]" << std::endl;

  int cgIdx = 0;
  const ChargeGroupVector& cgVec = pSolver->getChargeGroupVector();
  for (ChargeGroupVectorIt cgIt = cgVec.begin(); cgIt != cgVec.end(); cgIt++, cgIdx++)
  {
    out << "\tsubgraph cluster_" << cgIdx << " {" << std::endl; 
    //out << "\t\tlabel=" << std::fixed << computeError(*cgIt) << std::endl;
    out << "\t\tlabel=\"\"" << std::endl;

    //out << "\t\t" << "err_" << cgIdx 
    //  << "[label=" << computeError(*cgIt) 
    //  << ", shape=plaintext]" << std::endl;

    for (ChargeGroupIt nodeIt = cgIt->begin(); nodeIt != cgIt->end(); nodeIt++)
    {
      out << "\t\t" << _label[*nodeIt] 
        << " [label=\"" << _label2[*nodeIt] << " (" << _label[*nodeIt] << ")\\n"
        << std::setprecision(4) << pSolver->getCorrectedPartialCharge(*nodeIt) << "\""
        << ",color=" << color[pSolver->getColor(*nodeIt) % nColors] << "]"
        //<< ", colorscheme=paired12]"
        << std::endl;
    }
    out << "\t}" << std::endl;
  }

  //for (NodeIt n(_g); n != lemon::INVALID; ++n)
  //{
  //  out << "\t" << _label[n] 
  //    << " [color=" << color[pSolver->getColor(n) % nColors] << "]"
  //    //<< ", colorscheme=paired12]"
  //    << std::endl;
  //}

  // add a node containing detailed information
  std::stringstream ss; 
  pSolver->printChargeGroups(ss);
  std::string str = ss.str();

  // replace newline characters by \n in str
  size_t pos;
  while ((pos = str.find("\n")) != std::string::npos)
    str.replace(pos, 1, "\\n");

  out << "\tkey [labeljust=l,style=solid,shape=plaintext,label=\"" 
    << str << "\"]" << std::endl;

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    out << "\t" << _label[_g.u(e)] << " -- " << _label[_g.v(e)] << std::endl;
  }

  out << "}" << std::endl;
}

template<typename GR>
inline void Molecule<GR>::printDOTnoLegend(const Solver<GR>* pSolver,
                                   std::ostream& out) const
{
  static const int nColors = 11;
  static std::string color[nColors] =
  {
    "darkolivegreen1",
    "red",
    "cyan",
    "green",
    "yellow",
    "crimson",
    "gold",
    "violet",
    "aquamarine",
    "crimson",
    "chartreuse"
  };

  out << "graph G {" << std::endl;
  out << "\tlabelloc=t" << std::endl;
  out << "\tnode [style=filled]" << std::endl;

  int cgIdx = 0;
  const ChargeGroupVector& cgVec = pSolver->getChargeGroupVector();
  for (ChargeGroupVectorIt cgIt = cgVec.begin(); cgIt != cgVec.end(); cgIt++, cgIdx++)
  {
    out << "\tsubgraph cluster_" << cgIdx << " {" << std::endl;

    std::stringstream ss;
    ss.precision(3);
    pSolver->getTotalCharge(*cgIt, ss);
    std::string str = ss.str();
    // replace newline characters by \n in str
    //size_t pos;
    //while ((pos = str.find("\n")) != std::string::npos)
    //  str.replace(pos, 1, "\\n");

    out << "\t\tlabel=\"" << str
                          << "\""
                          << std::endl;
    out << "\t\tlabeljust=\"r\""
        << std::endl;

    for (ChargeGroupIt nodeIt = cgIt->begin(); nodeIt != cgIt->end(); nodeIt++)
    {
      out << "\t\t" << _label[*nodeIt]
        << " [label=\"" << _label2[*nodeIt] << " (" << _label[*nodeIt] << ")\\n"
        << std::setprecision(4) << pSolver->getCorrectedPartialCharge(*nodeIt) << "\""
        << ",color=" << color[pSolver->getColor(*nodeIt) % nColors] << "]"
        //<< ", colorscheme=paired12]"
        << std::endl;
    }
    out << "\t}" << std::endl;
  }

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    out << "\t" << _label[_g.u(e)] << " -- " << _label[_g.v(e)] << std::endl;
  }

  out << "}" << std::endl;
}

} // namespace cgp
} // namespace nina

#endif /* MOLECULE_H_ */
