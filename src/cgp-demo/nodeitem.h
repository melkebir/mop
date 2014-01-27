/* nodeitem.h
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#ifndef NODEITEM_H
#define NODEITEM_H

#include <QGraphicsItem>
#include <QList>
#include <lemon/list_graph.h>
#include "moleculeviz.h"

// forward class declarations
class EdgeItem;
class GraphWidget;

QT_BEGIN_NAMESPACE
class QGraphicsSceneMouseEvent;
QT_END_NAMESPACE

class NodeItem : public QObject, public QGraphicsItem
{
  Q_OBJECT

public:
  /// Graph type
  typedef lemon::ListGraph Graph;
  /// Molecule type
  typedef MoleculeViz<Graph> MoleculeVizType;

private:
  /// LEMON graph typedefs
  GRAPH_TYPEDEFS(Graph);

public:
  /// Constructor
  ///
  /// \param molecule is the visualized molecule
  /// \param node is an atom of the visualized molecule
  /// \param pGraphWidget is a reference to the widget on which the node is drawn
  /// \param cgIdx is the initial charge group index to which the node belongs
  NodeItem(const MoleculeVizType& molecule,
           Node node,
           GraphWidget* pGraphWidget,
           int cgIdx);

  /// Adds an edge to the node
  ///
  /// \param pEdge is the edge that is added to the node
  void addIncEdge(EdgeItem* pEdge);

  /// Returns all the incident edges
  QList<EdgeItem*> getIncEdges() const { return _edgeList; }

  /// Returns type of the graphics item
  int type() const { return UserType + 1; }

  /// Returns the bounding rectangle of this node
  QRectF boundingRect() const;

  /// Returns bounding shape of the node
  QPainterPath shape() const;

  /// Paints the node
  void paint(QPainter* pPainter,
             const QStyleOptionGraphicsItem* pOption,
             QWidget* pWidget);

  /// Returns the corresponding LEMON node
  Graph::Node getNode() const;

  /// Returns the charge group index to which the node belongs
  int getCgIdx() const;
  /// Sets the charge group index to which the node belongs
  ///
  /// \param cgIdx is the new charge group index
  void setCgIdx(int cgIdx);

  /// Returns the color of the node
  int getColorIdx() const;
  /// Sets the node color
  ///
  /// \param colorIdx is the new node color
  void setColorIdx(int colorIdx);

  /// Selects the node
  void select();
  /// Deselects the node
  void deselect();
  /// Selects the node as belonging to the selected charge group
  void selectCG();
  /// Deselects the node as belonging to the selected charge group
  void deselectCG();

  /// Returns whether the node is selected
  bool isSelected() const;
  /// Returns whether the node belongs to the selected charge group
  bool isSelectedCG() const;

  /// Returns the graph widget to which the node belongs
  GraphWidget* getGraphWidget();

protected:
  /// Handles node movement (by updating all incident edges)
  QVariant itemChange(GraphicsItemChange change,
                      const QVariant& value);
  /// Results in an update as to show the node as being selected
  void mousePressEvent(QGraphicsSceneMouseEvent* pEvent);
  /// Notifies that the node has been clicked
  void mouseReleaseEvent(QGraphicsSceneMouseEvent* pEvent);

  /// Returns the current color
  ///
  /// \param light indicates whether to return a light variant of the current color
  QColor getColor(bool light) const;

protected:
  /// Molecule that is visualized
  const MoleculeVizType& _molecule;
  /// List of incident edges
  QList<EdgeItem*> _edgeList;
  /// LEMON node that this class represents
  Node _node;
  /// Widget on which everything is drawn
  GraphWidget* _pGraph;
  /// Charge group index
  int _cgIdx;
  /// Node color
  int _colorIdx;
  /// Indicates whether the node is selected
  bool _selected;
  /// Indicates whether the charge group to which the node belongs is selected
  bool _selectedCG;

public:
  /// Radius of a node
  static const double _radius;

signals:
  /// Emitted upon clicking
  void selected(NodeItem* pNodeItem, bool addRemove);
};

inline NodeItem::Graph::Node NodeItem::getNode() const
{
  return _node;
}

inline int NodeItem::getCgIdx() const
{
  return _cgIdx;
}

inline void NodeItem::setCgIdx(int cgIdx)
{
  _cgIdx = cgIdx;
  update();
}

inline void NodeItem::select()
{
  if (!_selected)
  {
    _selected = true;
    update();
  }
}

inline void NodeItem::deselect()
{
  if (_selected)
  {
    _selected = false;
    update();
  }
}

inline void NodeItem::selectCG()
{
  if (!_selectedCG)
  {
    _selectedCG = true;
    update();
  }
}

inline void NodeItem::deselectCG()
{
  if (_selectedCG)
  {
    _selectedCG = false;
    update();
  }
}

inline bool NodeItem::isSelected() const
{
  return _selected;
}

inline bool NodeItem::isSelectedCG() const
{
  return _selectedCG;
}


inline int NodeItem::getColorIdx() const
{
  return _colorIdx;
}

inline GraphWidget* NodeItem::getGraphWidget()
{
  return _pGraph;
}


#endif /* NODEITEM_H */
