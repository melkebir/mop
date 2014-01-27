/* graphwidget.h
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H

#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGraphicsScene>
#include <QtWidgets>
#include <math.h>
#include <lemon/list_graph.h>

#include "moleculeviz.h"

// forward class declarations
class NodeItem;
class EdgeItem;

class GraphWidget : public QGraphicsView
{
  Q_OBJECT

public:
  /// Graph type
  typedef lemon::ListGraph Graph;
  /// Molecule type
  typedef MoleculeViz<Graph> MoleculeVizType;
  /// Charge group type
  typedef MoleculeVizType::ChargeGroup ChargeGroup;
  /// Charge group iterator
  typedef MoleculeVizType::ChargeGroupIt ChargeGroupIt;
  /// Charge group vector type
  typedef MoleculeVizType::ChargeGroupVector ChargeGroupVector;
  /// Charge group vector iterator
  typedef MoleculeVizType::ChargeGroupVectorIt ChargeGroupVectorIt;
  /// Weights on nodes
  typedef MoleculeVizType::WeightNodeMap WeightNodeMap;

private:
  /// LEMON graph typedefs
  GRAPH_TYPEDEFS(Graph);
  /// Index set type
  typedef std::set<size_t> IndexSet;
  /// Index set iterator type
  typedef IndexSet::const_iterator IndexSetIt;


protected:
  /// Rescales the whole thing upon a resize event
  void resizeEvent(QResizeEvent* pEvent);
  /// In case the background is clicked, the current selection is reset
  void mouseReleaseEvent(QMouseEvent* pEvent);

public:
  /// Node graphics item map type
  typedef Graph::NodeMap<NodeItem*> NodeItemMap;

  /// Constructor
  ///
  /// \param pParent is the parent of this object
  GraphWidget(QWidget* pParent = 0);
  /// Destructor
  ~GraphWidget();

  /// Prints all the charge groups
  ///
  /// \param out is the output stream
  void printChargeGroupVector(std::ostream& out) const;
  /// Returns the node (graphics) item map
  const NodeItemMap& getNodeItemMap() const;
  /// Returns the molecule
  const MoleculeVizType& getMolecule() const;
  /// Returns the molecule
  MoleculeVizType& getMolecule();
  /// Returns the charge groups
  const ChargeGroupVector& getChargeGroups() const;
  /// Returns the selected (graphics) node item
  NodeItem* getSelectedNodeItem();
  /// Returns the optimal objective value
  double getOptObj() const;
  /// Returns whether the solution is shown
  bool getShowSol() const;
  /// Returns the number of correct charge groups
  size_t getNumberOfCorrectChargeGroups() const;

private:
  /// Arc look up type
  typedef lemon::ArcLookUp<Graph> ArcLookUpType;

  /// Initializes the scene by drawing
  /// the nodes, the edges and computing an initial coloring
  ///
  /// \param cgVec are the charge groups
  /// \param pScene is the scene on which everything is drawn
  /// \param pNodeItem is the node item map,
  ///        mapping nodes to (graphics) node items
  /// \param pSelectedNodeItem is the currently selected node item
  void init(ChargeGroupVector& cgVec,
            QGraphicsScene*& pScene,
            NodeItemMap*& pNodeItem,
            NodeItem*& pSelectedNodeItem);
  /// Draws the nodes on the scene
  ///
  /// \param cgVec are the charge groups
  /// \param pScene is the scene on which everything is drawn
  /// \param pNodeItem is the node item map,
  ///        mapping nodes to (graphics) node items
  void initNodes(ChargeGroupVector& cgVec,
                 QGraphicsScene*& pScene,
                 NodeItemMap*& pNodeItem);
  /// Draws the edges on the scene
  ///
  /// \param pScene is the scene on which everything is drawn
  /// \param pNodeItem is the node item map,
  ///        mapping nodes to (graphics) node items
  void initEdges(QGraphicsScene*& pScene, NodeItemMap*& pNodeItem);
  /// Computes an initial coloring using a BFS
  ///
  /// \param cgVec are the charge groups
  /// \param pNodeItem is the node item map,
  ///        mapping nodes to (graphics) node items
  void initColors(ChargeGroupVector& cgVec, NodeItemMap*& pNodeItem);

  /// Sets the color of specified charge group
  ///
  /// \param cgIdx is the charge group to color
  /// \param colorIdx is the color to use
  void setColor(size_t cgIdx, size_t colorIdx);

  /// Gets available colors for the given node
  ///
  /// \param node is the LEMON node
  IndexSet getAvailableColors(Node node) const;

  /// Gets a (random) available color
  ///
  /// \param node is the LEMON node
  size_t getAvailableColor(Node node) const;

  /// Recolors the nodes by performing a BFS from the specified node,
  /// only nodes that have conflicting colors are recolored (recursively)
  ///
  /// \param visited is a map indicating which nodes have already been visited
  /// \param node is from which the BFS is performed
  void recolor(BoolNodeMap& visited, Node node);
  /// Shorthand function where only a node needs to be specified
  ///
  /// \param node is from which the BFS is performed
  void recolor(Node node);

  /// Finds an empty charge group
  ///
  /// \param cgIdx initially contains the charge group index
  ///        from which the search is initiated,
  ///        upon termination it either corresponds to an empty charge group
  ///        or it is set to |V| if no empty slot was found
  void findEmptyCg(size_t& cgIdx) const;
  /// Breaks up the charge group at the specified node
  ///
  /// \param node is the node to break up
  void breakUp(Node node);
  /// Returns whether the specified node
  /// can be added to the specified charge group
  ///
  /// \param cgIdx is the charge group
  /// \param node is the node
  bool allowed(size_t cgIdx, Node node) const;

  /// Selects the specified charge group
  ///
  /// \param cgIdx is the charge group
  void selectChargeGroup(size_t cgIdx);
  /// Deselects the specified charge group
  ///
  /// \param cgIdx is the charge group
  void deselectChargeGroup(size_t cgIdx);

  /// Sets the currently selected node (graphics) item
  ///
  /// \param pNodeItem is a node item that is selected
  void setSelectedNodeItem(NodeItem* pNodeItem);

  /// Returns the charge groups
  ChargeGroupVector& getChargeGroupsNonConst();

  /// Returns the node (graphics) item map
  NodeItemMap& getNodeItemMapNonConst();

  /// Shows the solution scene
  void showSolutionScene();
  /// Shows the normal scene
  void showNormalScene();

  /// Updates the number of correct charge groups
  void updateNrCorrectCGs();

  /// Molecule that is visualized
  MoleculeVizType* _pMolecule;
  /// Allows for efficiently determining whether two nodes are adjacent
  ArcLookUpType* _pArcLookUp;
  /// Node item map of the normal scene
  NodeItemMap* _pNodeItem;
  /// Node item map of the solution scene
  NodeItemMap* _pNodeItemSol;
  /// Selected node item in the normal scene
  NodeItem* _pSelectedNodeItem;
  /// Selected node item in the solution scene
  NodeItem* _pSelectedNodeItemSol;
  /// Charge groups used in the normal scene
  ChargeGroupVector _cgVec;
  /// Charge groups used in the solution scene
  ChargeGroupVector _cgVecSol;
  /// The normal scene
  QGraphicsScene* _pScene;
  /// The solution scene, nothing can be changed here
  QGraphicsScene* _pSceneSol;
  /// Optimal objective value
  double _optObj;
  /// Indicates whether the solution scene is on
  bool _showSol;
  /// Number of correct charge groups
  size_t _nCorrectCGs;
  /// Indicates whether it is solved
  bool _solved;

public:
  /// Default maximum charge group size
  static const size_t _defaultMaxChargeGroupSize;
  /// Default number of colors available
  static const size_t _nColors;
  /// Available colors
  static const Qt::GlobalColor _colors[];
  /// Available (dark) colors
  static const Qt::GlobalColor _darkColors[];

public slots:
  /// Parses the specified LGF file
  ///
  /// \param lgf is the filename of a LGF file
  /// \return \c true if the file is successfully parsed
  bool parseLGF(const QString& lgf);
  /// Writes current molecule to specified filename
  ///
  /// \param filename to which the molecule is written
  /// \return \c true if the file is successfully created
  bool writeLGF(const QString& filename);
  /// Computes a solution
  ///
  /// \param twFile contains the tree decomposition
  /// \return \c true if a solution is computed successfully
  bool initSolution(const QString& twFile);
  /// Sets the maximum charge group size
  ///
  /// \param k is the maximum charge group size
  void setK(size_t k);
  /// Changes the currently selected node
  ///
  /// \param pNodeItem is the node item that should be selected
  /// \param addRemove indicates whether the charge groups may be changed
  void selectNodeItem(NodeItem* pNodeItem, bool addRemove);
  /// Toggles the scene between solution and normal
  void toggleScene();
  /// Clears all charge groups
  void clear();

signals:
  /// Emitted when the key needs to be updated
  ///
  /// \param molecule is the visualized molecule
  /// \param nodeItem is the node item map, mapping LEMON nodes
  ///        to graphics node items
  /// \param cgVec are the charge groups
  /// \param pNodeItem is the currently selected node graphics item
  /// \param optObj is the optimal objective value
  /// \param solved will contain whether we are done
  void updateKey(const MoleculeVizType& molecule,
                 const NodeItemMap& nodeItem,
                 const ChargeGroupVector& cgVec,
                 const NodeItem* pNodeItem,
                 double optObj,
                 bool& solved);
  /// Emitted when the error label needs to be updated
  void updateErrLabel();
  /// Emitted when the number of correct charge groups label needs to be updated
  void updateCorrectCgLabel();
  /// Emitted when solved
  void solved();
};

inline void GraphWidget::recolor(Node node)
{
  assert(_pMolecule);

  BoolNodeMap visited(_pMolecule->getGraph());
  recolor(visited, node);
}

inline double GraphWidget::getOptObj() const
{
  return _optObj;
}

inline bool GraphWidget::getShowSol() const
{
  return _showSol;
}

inline const GraphWidget::NodeItemMap& GraphWidget::getNodeItemMap() const
{
  return _showSol ? *_pNodeItemSol : *_pNodeItem;
}

inline GraphWidget::NodeItemMap& GraphWidget::getNodeItemMapNonConst()
{
  return _showSol ? *_pNodeItemSol : *_pNodeItem;
}

inline const GraphWidget::MoleculeVizType& GraphWidget::getMolecule() const
{
  assert(_pMolecule);
  return *_pMolecule;
}

inline GraphWidget::MoleculeVizType& GraphWidget::getMolecule()
{
  assert(_pMolecule);
  return *_pMolecule;
}

inline const GraphWidget::ChargeGroupVector& GraphWidget::getChargeGroups() const
{
  return _showSol ? _cgVecSol : _cgVec;
}

inline GraphWidget::ChargeGroupVector& GraphWidget::getChargeGroupsNonConst()
{
  return _showSol ? _cgVecSol : _cgVec;
}

inline NodeItem* GraphWidget::getSelectedNodeItem()
{
  return _showSol ? _pSelectedNodeItemSol : _pSelectedNodeItem;
}

inline void GraphWidget::setSelectedNodeItem(NodeItem* pNodeItem)
{
  if (_showSol)
    _pSelectedNodeItemSol = pNodeItem;
  else
    _pSelectedNodeItem = pNodeItem;
}

inline size_t GraphWidget::getNumberOfCorrectChargeGroups() const
{
  if (_showSol)
    return _cgVecSol.size();
  else
    return _nCorrectCGs;
}

#endif /* GRAPHWIDGET_H */
