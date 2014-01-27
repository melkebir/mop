/* errorkeywidget.cpp
 *
 *  Created on: 14-Aug-2012
 *      Author: M. El-Kebir
 */

#include <QPainter>
#include <QStyleOption>

#include "nodeitem.h"
#include "errorkeyitem.h"

ErrorKeyItem::ErrorKeyItem(const QString& text, QObject* pParent)
  : QObject(pParent)
  , _text(text)
  , _boundingBox()
{
}

void ErrorKeyItem::paint(QPainter* pPainter,
                         const QStyleOptionGraphicsItem*,
                         QWidget*)
{
  QFont font = pPainter->font();
  font.setPointSize(14);
  pPainter->setFont(font);
  pPainter->setPen(Qt::black);

  pPainter->drawText(0, 0, _text);

  QFontMetrics fm = pPainter->fontMetrics();
  _boundingBox = fm.boundingRect(_text);
}

QRectF ErrorKeyItem::boundingRect() const
{
  return _boundingBox;
}
