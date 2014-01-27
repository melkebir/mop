/* errorkeywidget.h
 *
 *  Created on: 14-Aug-2012
 *      Author: M. El-Kebir
 */

#ifndef ERRORKEYITEM_H
#define ERRORKEYITEM_H

#include <QGraphicsItem>
#include <QList>

class ErrorKeyItem : public QObject, public QGraphicsItem
{
  Q_OBJECT

public:
  /// Constructor
  ///
  /// \param text that is shown
  /// \param pParent parent object
  ErrorKeyItem(const QString& text, QObject* pParent);

  /// Returns type of the graphics item
  int type() const { return Type + 4; }

  /// Paints the label
  void paint(QPainter* pPainter,
             const QStyleOptionGraphicsItem* pOption,
             QWidget* pWidget);
  /// Returns the bounding box
  QRectF boundingRect() const;

private:
  /// Text
  QString _text;
  /// Rectangle
  QRect _boundingBox;
};

#endif // ERRORKEYITEM_H
