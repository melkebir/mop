/* edgeitem.cpp
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#include <QPainter>
#include <math.h>

#include "edgeitem.h"
#include "nodeitem.h"

EdgeItem::EdgeItem(NodeItem* pSourceNode, NodeItem* pTargetNode)
  : _pSourceNode(pSourceNode)
  , _pTargetNode(pTargetNode)
  , _sourcePoint()
  , _targetPoint()
{
  setAcceptedMouseButtons(0);

  _pSourceNode->addIncEdge(this);
  _pTargetNode->addIncEdge(this);

  adjust();
}

void EdgeItem::adjust()
{
  if (!_pSourceNode || !_pTargetNode)
    return;

  QLineF line(mapFromItem(_pSourceNode, 0, 0), mapFromItem(_pTargetNode, 0, 0));
  qreal length = line.length();

  prepareGeometryChange();

  if (length > qreal(NodeItem::_radius * 2))
  {
    QPointF edgeOffset((line.dx() * NodeItem::_radius) / length,
                       (line.dy() * NodeItem::_radius) / length);
    _sourcePoint = line.p1() + edgeOffset;
    _targetPoint = line.p2() - edgeOffset;
  }
  else
  {
    _sourcePoint = _targetPoint = line.p1();
  }
}

QRectF EdgeItem::boundingRect() const
{
  if (!_pSourceNode || !_pTargetNode)
    return QRectF();

  qreal penWidth = 1;
  qreal extra = penWidth / 2.0;

  return QRectF(_sourcePoint, QSizeF(_targetPoint.x() - _sourcePoint.x(),
                                     _targetPoint.y() - _sourcePoint.y()))
    .normalized()
    .adjusted(-extra, -extra, extra, extra);
}

void EdgeItem::paint(QPainter* pPainter, const QStyleOptionGraphicsItem*, QWidget*)
{
  if (!_pSourceNode || !_pTargetNode)
    return;

  QLineF line(_sourcePoint, _targetPoint);
  if (qFuzzyCompare(line.length(), qreal(0)))
    return;

  // Draw the line itself
  pPainter->setPen(QPen(Qt::black, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  pPainter->drawLine(line);
}
