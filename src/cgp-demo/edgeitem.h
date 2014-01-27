/* edgeitem.h
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#ifndef EDGEITEM_H
#define EDGEITEM_H

#include <QGraphicsItem>
#include <lemon/list_graph.h>

// forward class declarations
class NodeItem;

class EdgeItem : public QGraphicsItem
{
public:
  /// Graph type
  typedef lemon::ListGraph Graph;

private:
  /// LEMON graph typedefs
  GRAPH_TYPEDEFS(Graph);

public:
  /// Constructor
  ///
  /// \param pSourceNode is the source node
  /// \param pTargetNode is the target node
  EdgeItem(NodeItem* pSourceNode, NodeItem* pTargetNode);

  /// Updates the line coordinates
  void adjust();

  /// Returns type of the graphics item
  int type() const { return UserType + 2; }

protected:
  /// Returns the bounding rectangle of this edge
  QRectF boundingRect() const;
  /// Paints the edge
  void paint(QPainter* pPainter,
             const QStyleOptionGraphicsItem* pOptions,
             QWidget* pWidget);

private:
  /// Source node
  NodeItem* _pSourceNode;
  /// Target node
  NodeItem* _pTargetNode;
  /// Source node coordinates
  QPointF _sourcePoint;
  /// Target node coordinates
  QPointF _targetPoint;
};

#endif /* EDGEITEM_H */
