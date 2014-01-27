/* nodekeyitem.h
 *
 *  Created on: 14-Aug-2012
 *      Author: M. El-Kebir
 */

#ifndef NODEKEYITEM_H
#define NODEKEYITEM_H

#include <QGraphicsItem>
#include <QList>
#include <lemon/list_graph.h>

#include "moleculeviz.h"
#include "nodeitem.h"
#include "graphwidget.h"

class NodeKeyItem : public QObject, public QGraphicsItem
{
  Q_OBJECT

public:
  /// Graph type
  typedef lemon::ListGraph Graph;
  /// Molecule type
  typedef MoleculeViz<Graph> MoleculeVizType;

private:
  GRAPH_TYPEDEFS(Graph);

public:
  /// Constructor
  ///
  /// \param pNodeItem is the visualized node
  NodeKeyItem(NodeItem* pNodeItem);

  /// Returns type of the graphics item
  int type() const { return UserType + 3; }

  /// Returns the bounding rectangle of this node
  QRectF boundingRect() const;
  /// Returns bounding shape of this shape
  QPainterPath shape() const;
  /// Paints the node
  void paint(QPainter* pPainter,
             const QStyleOptionGraphicsItem* pOption,
             QWidget* pWidget);

protected:
  /// Results in an update as to show the node as being selected
  void mousePressEvent(QGraphicsSceneMouseEvent* pEvent);
  /// Notifies that the node has been clicked
  void mouseReleaseEvent(QGraphicsSceneMouseEvent* pEvent);

private:
  /// Node
  NodeItem* _pNodeItem;

signals:
  /// Emitted upon clicking
  void selected(NodeItem* pNodeItem, bool addRemove);
};

#endif // NODEKEYITEM_H
