/* nodeitem.cpp
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include "nodeitem.h"
#include "edgeitem.h"
#include "graphwidget.h"

const double NodeItem::_radius = 23;

NodeItem::NodeItem(const MoleculeVizType& molecule,
                   Node node,
                   GraphWidget* pGraphWidget,
                   int cgIdx)
  : QObject(pGraphWidget)
  , _molecule(molecule)
  , _edgeList()
  , _node(node)
  , _pGraph(pGraphWidget)
  , _cgIdx(cgIdx)
  , _colorIdx(-1)
  , _selected(false)
  , _selectedCG(false)
{
  connect(this, SIGNAL(selected(NodeItem*, bool)),
          pGraphWidget, SLOT(selectNodeItem(NodeItem*, bool)));

#ifdef GRAPHVIZ
  setFlag(ItemIsMovable);
  setFlag(ItemSendsGeometryChanges);
#else
  setFlag(ItemIsSelectable);
#endif

  setCacheMode(DeviceCoordinateCache);
  setZValue(-1);
}

void NodeItem::addIncEdge(EdgeItem* pEdge)
{
  _edgeList << pEdge;
  pEdge->adjust();
}

QRectF NodeItem::boundingRect() const
{
  qreal adjust = 2;
  return QRectF(-_radius - adjust, -_radius - adjust,
                2 * _radius + adjust, 2 * _radius + adjust);
}

QPainterPath NodeItem::shape() const
{
  QPainterPath path;
  path.addEllipse(-_radius, -_radius, _radius * 2, _radius * 2);

  return path;
}

void NodeItem::paint(QPainter* pPainter,
                     const QStyleOptionGraphicsItem* pOption,
                     QWidget*)
{
  QRadialGradient gradient(-3, -3, _radius);
  if (pOption->state & QStyle::State_Sunken || _selected || _selectedCG)
  {
    gradient.setCenter(3, 3);
    gradient.setFocalPoint(3, 3);
    gradient.setColorAt(1, getColor(true).light(150));
    gradient.setColorAt(0, getColor(false).light(150));

    if (_selected)
      pPainter->setPen(QPen(Qt::red, 0));
  }
  else
  {
    gradient.setColorAt(0, getColor(true));
    gradient.setColorAt(1, getColor(false));
    pPainter->setPen(QPen(Qt::black, 0));
  }
  pPainter->setBrush(gradient);

  pPainter->drawEllipse(-_radius, -_radius, _radius * 2, _radius * 2);
  pPainter->setBrush(Qt::NoBrush);

  QFont font = pPainter->font();
  font.setPointSize(10);
  pPainter->setFont(font);
  pPainter->setPen(Qt::black);

  QString label;
  QString partialCharge;

  partialCharge = QString("%1").arg(_molecule.getPartialCharge(_node));
  label = QString("%1").arg(_molecule.getLabel2(_node).c_str());

  QFontMetrics fm = pPainter->fontMetrics();
  pPainter->drawText(-fm.width(label) / 2, -5, label);

  fm = pPainter->fontMetrics();
  pPainter->drawText(-fm.width(partialCharge) / 2, 10, partialCharge);
}

QVariant NodeItem::itemChange(GraphicsItemChange change, const QVariant& value)
{
  switch (change)
  {
  case ItemPositionHasChanged:
    foreach (EdgeItem* pEdge, _edgeList)
      pEdge->adjust();
    break;
  default:
    break;
  };

  return QGraphicsItem::itemChange(change, value);
}

void NodeItem::mousePressEvent(QGraphicsSceneMouseEvent* pEvent)
{
  update();
  QGraphicsItem::mousePressEvent(pEvent);
}

void NodeItem::mouseReleaseEvent(QGraphicsSceneMouseEvent* pEvent)
{
  emit selected(this, true);
  update();
  QGraphicsItem::mouseReleaseEvent(pEvent);
}

void NodeItem::setColorIdx(int colorIdx)
{
  assert(colorIdx < static_cast<int>(GraphWidget::_nColors));
  _colorIdx = colorIdx;
  update();
}

QColor NodeItem::getColor(bool light) const
{
  if (_colorIdx < 0)
    return QColor(Qt::white);

  if (light)
    return QColor(GraphWidget::_colors[_colorIdx]);
  else
    return QColor(GraphWidget::_darkColors[_colorIdx]);
}
