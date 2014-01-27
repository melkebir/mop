/* graphkeywidget.cpp
 *
 *  Created on: 14-Aug-2012
 *      Author: M. El-Kebir
 */

#include "graphkeywidget.h"
#include "graphwidget.h"
#include "errorkeyitem.h"
#include <algorithm>

GraphKeyWidget::GraphKeyWidget(QWidget* pParent)
  : QGraphicsView(pParent)
  , _pNodeKeyItemMap(NULL)
  , _curObj(0)
{
  setCacheMode(CacheBackground);
  setViewportUpdateMode(BoundingRectViewportUpdate);
  setRenderHint(QPainter::Antialiasing);
  setTransformationAnchor(AnchorUnderMouse);
  scale(qreal(0.8), qreal(0.8));

  QGraphicsScene* pScene = new QGraphicsScene(this);
  QRectF sceneRect(-NodeItem::_radius,
                   -NodeItem::_radius,
                   NodeItem::_radius * GraphWidget::_defaultMaxChargeGroupSize,
                   NodeItem::_radius);
  pScene->setSceneRect(sceneRect);
  setScene(pScene);
}

void GraphKeyWidget::resizeEvent(QResizeEvent*)
{
  QGraphicsScene* pScene = scene();
  fitInView(pScene->sceneRect(), Qt::KeepAspectRatio);
}

void GraphKeyWidget::update(const MoleculeVizType& molecule,
                            const NodeItemMap& nodeItem,
                            const ChargeGroupVector& cgVec,
                            const NodeItem* pSelectedNodeItem,
                            double optObj,
                            bool& solved)
{
  scene()->clear();

  const Graph& g = molecule.getGraph();
  const size_t k = molecule.getMaxChargeGroupSize();
  static const double horSpacing = 5;
  static const double verSpacing = 10;

  delete _pNodeKeyItemMap;
  _pNodeKeyItemMap = new NodeKeyItemMap(g, NULL);

  CgComp compare(molecule);
  ChargeGroupVector cpyCgVec = cgVec;

  // sort on error, but make sure that selected charge group is on top
  if (pSelectedNodeItem && cpyCgVec.size() > 1)
  {
    std::swap(cpyCgVec[0], cpyCgVec[pSelectedNodeItem->getCgIdx()]);
    std::sort(cpyCgVec.begin() + 1, cpyCgVec.end(), compare);
  }
  else
  {
    std::sort(cpyCgVec.begin(), cpyCgVec.end(), compare);
  }

  // get all atoms
  ChargeGroup atoms;
  for (NodeIt n(g); n != lemon::INVALID; ++n)
    atoms.insert(n);

  // array of size (k + 1) * l where l is the number of charge groups
  _curObj = 0;
  int row = 0;
  int maxNumberCol = 0;

  // add captions
  ErrorKeyItem* pCaption1 = new ErrorKeyItem("Charge group", this);
  pCaption1->setPos(-NodeItem::_radius, row * (NodeItem::_radius * 2 + verSpacing));
  scene()->addItem(pCaption1);

  ErrorKeyItem* pCaption2 = new ErrorKeyItem("Error", this);
  pCaption2->setPos(k * (NodeItem::_radius * 2 + horSpacing),
                    row * (NodeItem::_radius * 2 + verSpacing));
  scene()->addItem(pCaption2);
  row++;

  for (ChargeGroupVectorIt cgIt = cpyCgVec.begin(); cgIt != cpyCgVec.end(); cgIt++)
  {
    int col = 0;
    for (ChargeGroupIt nodeIt = cgIt->begin(); nodeIt != cgIt->end(); nodeIt++, col++)
    {
      atoms.erase(*nodeIt);

      NodeItem* pNodeItem = nodeItem[*nodeIt];
      NodeKeyItem* pNodeKeyItem = new NodeKeyItem(pNodeItem);

      pNodeKeyItem->setPos(col * (NodeItem::_radius * 2 + horSpacing),
                           row * (NodeItem::_radius * 2 + verSpacing));

      scene()->addItem(pNodeKeyItem);
    }

    if (col > 0)
    {
      double localErr = molecule.computeError(*cgIt, false, false);
      _curObj += fabs(localErr);

      ErrorKeyItem* pErrorItem = new ErrorKeyItem(QString::number(localErr, 'f', 2), this);

      pErrorItem->setPos(k * (NodeItem::_radius * 2 + horSpacing),
                         row * (NodeItem::_radius * 2 + verSpacing));

      scene()->addItem(pErrorItem);
    }

    if (maxNumberCol < col) maxNumberCol = col;
    if (cgIt->size() > 0) row++;
  }

  // update total error
  for (ChargeGroupIt nodeIt = atoms.begin(); nodeIt != atoms.end(); ++nodeIt)
  {
    ChargeGroup cg;
    cg.insert(*nodeIt);

    _curObj += molecule.computeError(cg, true, false);
  }

  // add total error
  ErrorKeyItem* pTotalErrorCaption = new ErrorKeyItem("Total error:", this);
  pTotalErrorCaption->setPos(-NodeItem::_radius, row * (NodeItem::_radius * 2 + verSpacing));
  scene()->addItem(pTotalErrorCaption);
  ErrorKeyItem* pTotalError = new ErrorKeyItem(QString("%1").arg(_curObj), this);
  pTotalError->setPos(k * (NodeItem::_radius * 2 + horSpacing),
                      row * (NodeItem::_radius * 2 + verSpacing));
  scene()->addItem(pTotalError);
  row++;

  // add optimal error
  //ErrorKeyItem* pOptErrorCaption = new ErrorKeyItem("Opt error:", this);
  //pOptErrorCaption->setPos(-NodeItem::_radius, row * (NodeItem::_radius * 2 + verSpacing));
  //scene()->addItem(pOptErrorCaption);
  //ErrorKeyItem* pOptError = new ErrorKeyItem(QString("%1").arg(optObj), this);
  //pOptError->setPos(k * (NodeItem::_radius * 2 + horSpacing),
  //                    row * (NodeItem::_radius * 2 + verSpacing));
  //scene()->addItem(pOptError);
  //row++;

  QRectF sceneRect(-NodeItem::_radius - horSpacing,
                   -NodeItem::_radius - verSpacing,
                   (NodeItem::_radius * 2 + horSpacing) * (k + 2),
                   (NodeItem::_radius * 2 + verSpacing) * row);
  scene()->setSceneRect(sceneRect);

  fitInView(scene()->sceneRect(), Qt::KeepAspectRatio);

  solved = optObj == _curObj;
}
