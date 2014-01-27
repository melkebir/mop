/* nodekeyitem.cpp
 *
 *  Created on: 14-Aug-2012
 *      Author: M. El-Kebir
 */

#include "nodekeyitem.h"

NodeKeyItem::NodeKeyItem(NodeItem* pNodeItem)
  : QObject(pNodeItem)
  , _pNodeItem(pNodeItem)
{
  connect(this, SIGNAL(selected(NodeItem*, bool)),
          pNodeItem->getGraphWidget(), SLOT(selectNodeItem(NodeItem*, bool)));
  setFlag(ItemIsSelectable);
}

QRectF NodeKeyItem::boundingRect() const
{
  return _pNodeItem->boundingRect();
}

QPainterPath NodeKeyItem::shape() const
{
  return _pNodeItem->shape();
}

void NodeKeyItem::paint(QPainter* pPainter,
                        const QStyleOptionGraphicsItem* pOption,
                        QWidget* pWidget)
{
  return _pNodeItem->paint(pPainter, pOption, pWidget);
}

void NodeKeyItem::mousePressEvent(QGraphicsSceneMouseEvent* pEvent)
{
  update();
  QGraphicsItem::mousePressEvent(pEvent);
}

void NodeKeyItem::mouseReleaseEvent(QGraphicsSceneMouseEvent*)
{
  emit selected(_pNodeItem, false);
}
