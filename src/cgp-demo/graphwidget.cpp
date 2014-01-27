/* graphwidget.cpp
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

// TODO:
// - available colors kan slimmer: het is nu te strict,
//   alleen kleuren die gebruikt worden in adj. nodes zijn unavailable
// - SmartGraph global typedef

#include "graphwidget.h"
#include "nodeitem.h"
#include "edgeitem.h"
#include "solvertreewidthmemdp.h"

#include <QtGui>
#include <math.h>
#include <limits>
#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>

#ifdef GRAPHVIZ
#include "moleculegraphviz.h"
#endif

const size_t GraphWidget::_defaultMaxChargeGroupSize = 5;

const size_t GraphWidget::_nColors = 6;

const Qt::GlobalColor GraphWidget::_colors[] =
{
  Qt::red,
  Qt::green,
  Qt::gray,
  Qt::yellow,
  Qt::magenta,
  Qt::cyan
};

const Qt::GlobalColor GraphWidget::_darkColors[] =
{
  Qt::darkRed,
  Qt::darkGreen,
  Qt::darkGray,
  Qt::darkYellow,
  Qt::darkMagenta,
  Qt::darkCyan
};

GraphWidget::GraphWidget(QWidget* pParent)
  : QGraphicsView(pParent)
  , _pMolecule(NULL)
  , _pArcLookUp(NULL)
  , _pNodeItem(NULL)
  , _pNodeItemSol(NULL)
  , _pSelectedNodeItem(NULL)
  , _pSelectedNodeItemSol(NULL)
  , _cgVec()
  , _cgVecSol()
  , _pScene(NULL)
  , _pSceneSol(NULL)
  , _optObj(0)
  , _showSol(false)
  , _nCorrectCGs(0)
  , _solved(false)
{
  setMinimumSize(400, 400);
  setWindowTitle(tr("Charge Group Partioning Visualization"));
  setCacheMode(CacheBackground);
  setViewportUpdateMode(BoundingRectViewportUpdate);
  setRenderHint(QPainter::Antialiasing);
  setTransformationAnchor(AnchorUnderMouse);
  scale(qreal(0.8), qreal(0.8));

#ifdef GRAPHVIZ
  _pMolecule = new MoleculeGraphviz<Graph>(_defaultMaxChargeGroupSize, NodeItem::_radius);
#else
  _pMolecule = new MoleculeVizType(_defaultMaxChargeGroupSize);
#endif
}

GraphWidget::~GraphWidget()
{
  delete _pNodeItem;
  delete _pNodeItemSol;
  delete _pArcLookUp;
  delete _pMolecule;
}

void GraphWidget::init(ChargeGroupVector& cgVec,
                       QGraphicsScene*& pScene,
                       NodeItemMap*& pNodeItem,
                       NodeItem*& pSelectedNodeItem)
{
  pSelectedNodeItem = NULL;
  _solved = false;
  _nCorrectCGs = 0;

  delete pScene;
  pScene = new QGraphicsScene(this);

  // let's add the nodes
  initNodes(cgVec, pScene, pNodeItem);

  // and the edges
  initEdges(pScene, pNodeItem);

  initColors(cgVec, pNodeItem);
}

void GraphWidget::initNodes(ChargeGroupVector& cgVec,
                            QGraphicsScene*& pScene,
                            NodeItemMap*& pNodeItem)
{
  assert(_pMolecule);

  QPointF minPos(std::numeric_limits<qreal>::max(),
                 std::numeric_limits<qreal>::max());
  QPointF maxPos(0, 0);

  const Graph& g = _pMolecule->getGraph();
  delete pNodeItem;
  pNodeItem = new NodeItemMap(g);

  // clear all the charge groups
  cgVec.clear();

  lemon::Bfs<Graph> bfs(g);
  bfs.init();
  bfs.addSource(Graph::NodeIt(g));

  //int colorIdx = 0;
  while (!bfs.emptyQueue())
  {
    Graph::Node n = bfs.processNextNode();

    // new charge group only containing the node
    ChargeGroup cg;
    //cg.insert(n);
    cgVec.push_back(cg);

    NodeItem* pNode = new NodeItem(*_pMolecule, n, this, -1);

    pNodeItem->set(n, pNode);

    QPointF pos(_pMolecule->getCoordX(n), _pMolecule->getCoordY(n));
    pNode->setPos(pos);

    if (pos.x() < minPos.x()) minPos.setX(pos.x());
    if (pos.y() < minPos.y()) minPos.setY(pos.y());
    if (pos.x() > maxPos.x()) maxPos.setX(pos.x());
    if (pos.y() > maxPos.y()) maxPos.setY(pos.y());

    pScene->addItem(pNode);
  }

  minPos -= QPointF(NodeItem::_radius, NodeItem::_radius);
  maxPos += QPointF(NodeItem::_radius, NodeItem::_radius);

  QRectF sceneRect(minPos, maxPos);
  _pScene->setSceneRect(sceneRect);
  fitInView(sceneRect, Qt::KeepAspectRatio);
}

void GraphWidget::initEdges(QGraphicsScene*& pScene, NodeItemMap*& pNodeItem)
{
  assert(_pMolecule);
  const Graph& g = _pMolecule->getGraph();

  for (Graph::EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    EdgeItem* pEdge = new EdgeItem((*pNodeItem)[g.u(e)], (*pNodeItem)[g.v(e)]);
    pScene->addItem(pEdge);
  }
}

void GraphWidget::initColors(ChargeGroupVector& cgVec,
                             NodeItemMap*& pNodeItem)
{
  size_t cgIdx = 0;
  for (ChargeGroupVectorIt it = cgVec.begin(); it != cgVec.end(); it++, cgIdx++)
  {
    for (ChargeGroupIt nodeIt = it->begin(); nodeIt != it->end(); nodeIt++)
    {
      (*pNodeItem)[*nodeIt]->setCgIdx(cgIdx);
      (*pNodeItem)[*nodeIt]->setColorIdx(cgIdx % GraphWidget::_nColors);
    }
  }

  update();
}

void GraphWidget::resizeEvent(QResizeEvent*)
{
  if (scene())
    fitInView(scene()->sceneRect(), Qt::KeepAspectRatio);
}

bool GraphWidget::parseLGF(const QString& lgf)
{
  assert(_pMolecule);

  std::ifstream in(lgf.toLatin1().data());
  if (!in.good())
  {
#ifndef ANDROID
    QMessageBox msgBox(this);
    msgBox.setText(QString("Failed to open '%1'").arg(lgf));
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.exec();
#endif

    return false;
  }
  else
  {
    _pMolecule->readLGF(in);

    delete _pArcLookUp;
    _pArcLookUp = new ArcLookUpType(_pMolecule->getGraph());

    init(_cgVec, _pScene, _pNodeItem, _pSelectedNodeItem);

    if (!_showSol)
      showNormalScene();

    return true;
  }
}

void GraphWidget::setK(size_t k)
{
  assert(_pMolecule);
  _pMolecule->setMaxChargeGroupSize(k);
}

void GraphWidget::selectNodeItem(NodeItem* pNodeItem, bool addRemove)
{
  assert(_pMolecule);

  ChargeGroupVector& cgVec = getChargeGroupsNonConst();

  NodeItem* pSelectedNodeItem = getSelectedNodeItem();
  if (pSelectedNodeItem == pNodeItem)
  {
    // deselect
    if (pSelectedNodeItem)
    {
      const size_t cgIdx = pSelectedNodeItem->getCgIdx();
      if (cgVec[cgIdx].size() > 0 && addRemove)
      {
        deselectChargeGroup(cgIdx);
        breakUp(pSelectedNodeItem->getNode());
        selectChargeGroup(cgIdx);

        updateNrCorrectCGs();
      }
    }
  }
  else
  {
    // deselect currently selected node
    if (pSelectedNodeItem)
      pSelectedNodeItem->deselect();

    if (pNodeItem)
    {
      size_t cgIdx = pNodeItem->getCgIdx();

      // if no cgIdx, assign cgIdx
      if (pNodeItem->getCgIdx() < 0)
      {
        cgIdx = 0;
        findEmptyCg(cgIdx);
        pNodeItem->setCgIdx(static_cast<int>(cgIdx));

        pNodeItem->setColorIdx(static_cast<int>(getAvailableColor(pNodeItem->getNode())));

        ChargeGroup cg;
        cg.insert(pNodeItem->getNode());

        _cgVec[cgIdx] = cg;

        updateNrCorrectCGs();
      }

      if (pSelectedNodeItem)
      {
        const size_t selectedCgIdx = pSelectedNodeItem->getCgIdx();
        ChargeGroup& selectedCg = cgVec[selectedCgIdx];
        const Node selectedNode = pNodeItem->getNode();

        // assign to the same charge group as the previously selected one
        if (selectedCgIdx != cgIdx && allowed(selectedCgIdx, selectedNode) && addRemove)
        {
          breakUp(selectedNode);

          cgVec[cgIdx].erase(selectedNode);

          pNodeItem->setCgIdx(selectedCgIdx);
          pNodeItem->setColorIdx(pSelectedNodeItem->getColorIdx());
          selectedCg.insert(selectedNode);

          recolor(selectedNode);

          updateNrCorrectCGs();
          selectChargeGroup(selectedCgIdx);
        }
        else
        {
          // select charge group
          deselectChargeGroup(selectedCgIdx);
          selectChargeGroup(cgIdx);
        }
      }
      else
      {
        // select charge group
        selectChargeGroup(cgIdx);
      }

      pNodeItem->select();
    }

    setSelectedNodeItem(pNodeItem);
  }

  emit updateKey(*_pMolecule, getNodeItemMap(), cgVec, pNodeItem, _optObj, _solved);
  emit updateErrLabel();
  emit updateCorrectCgLabel();
}

void GraphWidget::mouseReleaseEvent(QMouseEvent* pEvent)
{
  assert(_pMolecule);

  // COMMENT OUT BELOW: HACK FOR THE STUPID TOUCH TABLE

  NodeItem* pSelectedNodeItem = getSelectedNodeItem();
  if (!itemAt(pEvent->pos()) && pSelectedNodeItem)
  {
    deselectChargeGroup(pSelectedNodeItem->getCgIdx());
    selectNodeItem(NULL, true);
    emit updateKey(*_pMolecule,
                   getNodeItemMap(),
                   getChargeGroups(),
                   getSelectedNodeItem(),
                   _optObj,
                   _solved);
    update();
  }

  QGraphicsView::mouseReleaseEvent(pEvent);

  if (!_showSol && _solved)
    emit solved();

}

void GraphWidget::printChargeGroupVector(std::ostream& out) const
{
  assert(_pMolecule);

  const ChargeGroupVector& cgVec = getChargeGroups();
  for (ChargeGroupVectorIt it = cgVec.begin(); it != cgVec.end(); it++)
  {
    if ((*it).size())
      _pMolecule->printChargeGroup(*it, out);
  }
}

bool GraphWidget::allowed(size_t cgIdx, Node node) const
{
  assert(_pMolecule);

  if (_showSol)
    return false;

  const size_t k = _pMolecule->getMaxChargeGroupSize();
  assert(cgIdx < _pMolecule->getNumberOfAtoms());

  const ChargeGroup& cg = getChargeGroups()[cgIdx];

  // first we check whether cg can be extended
  if (cg.size() == k)
    return false;

  // now we check whether node is connected to cg
  for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
  {
    if ((*_pArcLookUp)(*cgIt, node) != lemon::INVALID)
      return true;
  }

  return false;
}

void GraphWidget::deselectChargeGroup(size_t cgIdx)
{
  assert(_pMolecule);
  assert(cgIdx < _pMolecule->getNumberOfAtoms());

  const ChargeGroup& cg = getChargeGroups()[cgIdx];
  for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
  {
    NodeItem* pNodeItem = getNodeItemMap()[*cgIt];
    pNodeItem->deselectCG();
  }
}

void GraphWidget::selectChargeGroup(size_t cgIdx)
{
  assert(cgIdx < _pMolecule->getNumberOfAtoms());

  const ChargeGroup& cg = getChargeGroups()[cgIdx];
  for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
  {
    NodeItem* pNodeItem = getNodeItemMap()[*cgIt];
    pNodeItem->selectCG();
  }
}

void GraphWidget::findEmptyCg(size_t& cgIdx) const
{
  const ChargeGroupVector& cgVec = getChargeGroups();
  for (; cgIdx < cgVec.size() && cgVec[cgIdx].size() != 0; cgIdx++);
}

void GraphWidget::breakUp(Node node)
{
  assert(_pMolecule);

  size_t cgIdx = (getNodeItemMap())[node]->getCgIdx();
  const size_t nAtoms = _pMolecule->getNumberOfAtoms();
  assert(cgIdx < nAtoms);

  const Graph& g = _pMolecule->getGraph();
  BoolEdgeMap edgeMap(g, true);
  BoolNodeMap nodeMap(g, false);

  // make a new subgraph induced by the nodes in the selected charge group
  ChargeGroup& cg = getChargeGroupsNonConst()[cgIdx];
  for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
  {
    if (*cgIt != node)
      nodeMap[*cgIt] = true;
  }

  lemon::SubGraph<const Graph> subG(g, nodeMap, edgeMap);
  lemon::SubGraph<const Graph>::NodeMap<int> comp(subG);

  // find all connected components
  int nComp = lemon::connectedComponents(subG, comp);
  std::vector<size_t> compCgIdx(nComp, nAtoms);

  size_t emptyCgIdx = 0;
  for (lemon::SubGraph<const Graph>::NodeIt n(subG); n != lemon::INVALID; ++n)
  {
    size_t newCgIdx = compCgIdx[comp[n]];
    if (newCgIdx == nAtoms)
    {
      findEmptyCg(emptyCgIdx);
      compCgIdx[comp[n]] = emptyCgIdx;
      newCgIdx = emptyCgIdx;
    }

    getNodeItemMap()[n]->setCgIdx(newCgIdx);

    // remove from current cg and put in newCgIdx
    cg.erase(n);
    getChargeGroupsNonConst()[newCgIdx].insert(n);
  }

  recolor(node);
}

void GraphWidget::setColor(size_t cgIdx, size_t colorIdx)
{
  assert(_pMolecule);
  assert(cgIdx < _pMolecule->getNumberOfAtoms());
  assert(colorIdx < GraphWidget::_nColors);

  const ChargeGroup& cg = getChargeGroups()[cgIdx];
  for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
  {
    getNodeItemMap()[*cgIt]->setColorIdx(colorIdx);
  }
}

GraphWidget::IndexSet GraphWidget::getAvailableColors(Node node) const
{
  const Graph& g = _pMolecule->getGraph();

  IndexSet usedColors;
  for (IncEdgeIt e(g, node); e != lemon::INVALID; ++e)
  {
    Node adjNode = g.oppositeNode(node, e);
    NodeItem* pAdjNodeItem = getNodeItemMap()[adjNode];

    if (pAdjNodeItem->getColorIdx() >= 0)
    {
      usedColors.insert(static_cast<size_t>(pAdjNodeItem->getColorIdx()));
    }
  }

  IndexSet result;
  for (size_t color = 0; color < GraphWidget::_nColors; color++)
  {
    if (usedColors.find(color) == usedColors.end())
    {
      result.insert(color);
    }
  }

  return result;
}

size_t GraphWidget::getAvailableColor(Node node) const
{
  IndexSet availableColors = getAvailableColors(node);

  assert(availableColors.size() > 0);

  size_t colorIdx;
  int r = rand() % availableColors.size();
  for (IndexSetIt colorIt = availableColors.begin();
       colorIt != availableColors.end() && r >= 0; colorIt++, r--)
  {
    colorIdx = (*colorIt);
  }

  return colorIdx;
}

void GraphWidget::recolor(BoolNodeMap& visited, Node node)
{
  assert(_pMolecule);

  if (visited[node])
    return;

  const Graph& g = _pMolecule->getGraph();
  const size_t curCgIdx = getNodeItemMap()[node]->getCgIdx();
  const size_t curColorIdx = getNodeItemMap()[node]->getColorIdx();
  const ChargeGroupVector& cgVec = getChargeGroups();
  const ChargeGroup& cg = cgVec[curCgIdx];

  // determine adjacent charge groups and set visited flag to all nodes in curCgIdx
  visited[node] = true;

  std::set<std::pair<size_t, size_t> > adjCg;
  for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
  {
    for (IncEdgeIt e(g, *cgIt); e != lemon::INVALID; ++e)
    {
      Node v = g.oppositeNode(*cgIt, e);

      int cgIdx = getNodeItemMap()[v]->getCgIdx();
      int colorIdx = getNodeItemMap()[v]->getColorIdx();

      if (cgIdx >= 0 && curCgIdx != cgIdx)
        adjCg.insert(std::make_pair(static_cast<size_t>(cgIdx), static_cast<size_t>(colorIdx)));
      else
        visited[v] = true;
    }
  }

  // determine the available colors
  std::set<size_t> availableColor;
  for (size_t color = 0; color < GraphWidget::_nColors; color++)
    availableColor.insert(color);

  for (std::set<std::pair<size_t, size_t> >::const_iterator it = adjCg.begin();
       it != adjCg.end(); it++)
    availableColor.erase(it->second);

  if (availableColor.find(curColorIdx) != availableColor.end())
  {
    // no need to recolor and therefore no need to recurse
    return;
  }
  else
  {
    // recolor
    size_t colorIdx;

    int r = rand() % availableColor.size();
    for (std::set<size_t>::const_iterator colorIt = availableColor.begin();
         colorIt != availableColor.end() && r >= 0; colorIt++, r--)
    {
      colorIdx = (*colorIt);
    }

    for (ChargeGroupIt cgIt = cg.begin(); cgIt != cg.end(); cgIt++)
      getNodeItemMap()[*cgIt]->setColorIdx(colorIdx);
  }

   // recurse on adjacent charge groups
   for (std::set<std::pair<size_t, size_t> >::const_iterator it = adjCg.begin();
       it != adjCg.end(); it++)
   recolor(visited, *cgVec[it->first].begin());
}

bool GraphWidget::initSolution(const QString& twFile)
{
  assert(_pMolecule);

  std::ifstream tw(twFile.toLatin1().data());
  if (!tw.good())
  {
#ifndef ANDROID
    QMessageBox msgBox(this);
    msgBox.setText(QString("Failed to open '%1'").arg(QString(twFile)));
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.exec();
#endif

    return false;
  }

  nina::cgp::SolverTreewidthMemDP<Graph>* pSolver = new nina::cgp::SolverTreewidthMemDP<Graph>(_pMolecule);

  _pMolecule->initGlobalK(_pMolecule->getMaxChargeGroupSize());
  pSolver->solve(0);
  _optObj = pSolver->getError();

  // update _pSolutionScene
  init(_cgVecSol, _pSceneSol, _pNodeItemSol, _pSelectedNodeItemSol);

  bool oldShowSol = _showSol;
  _showSol |= true;
  _cgVecSol = pSolver->getChargeGroupVector();

  for (NodeIt n(_pMolecule->getGraph()); n != lemon::INVALID; ++n)
    ;

  initColors(_cgVecSol, _pNodeItemSol);
  recolor(NodeIt(_pMolecule->getGraph()));

  _showSol = oldShowSol;

  if (_showSol)
    showSolutionScene();
  else
  {
    emit updateKey(*_pMolecule,
                   getNodeItemMap(),
                   getChargeGroups(),
                   getSelectedNodeItem(),
                   _optObj,
                   _solved);
    emit updateErrLabel();
    emit updateCorrectCgLabel();
  }

  return true;
}

void GraphWidget::showNormalScene()
{
  assert(_pMolecule);

  _showSol = false;
  setScene(_pScene);
  emit updateKey(*_pMolecule,
                 getNodeItemMap(),
                 getChargeGroups(),
                 getSelectedNodeItem(),
                 _optObj,
                 _solved);
  emit updateErrLabel();
  emit updateCorrectCgLabel();
}

void GraphWidget::showSolutionScene()
{
  assert(_pMolecule);

  _showSol = true;
  setScene(_pSceneSol);
  emit updateKey(*_pMolecule,
                 getNodeItemMap(),
                 getChargeGroups(),
                 getSelectedNodeItem(),
                 _optObj,
                 _solved);
  emit updateErrLabel();
  emit updateCorrectCgLabel();
}

void GraphWidget::toggleScene()
{
  if (_showSol)
    showNormalScene();
  else
    showSolutionScene();
}

bool GraphWidget::writeLGF(const QString& filename)
{
  assert(_pMolecule);

  std::ofstream out(filename.toLatin1().data());
  if (!out.good())
  {
#ifndef ANDROID
    QMessageBox msgBox(this);
    msgBox.setText(QString("Failed to open '%1' for writing").arg(filename));
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.exec();
#endif

    return false;
  }
  else
  {
    // make sure that positions are correct
    const Graph& g = _pMolecule->getGraph();

    for (NodeIt n(g); n != lemon::INVALID; ++n)
    {
      NodeItemMap& nodeItemMap = getNodeItemMapNonConst();
      _pMolecule->setCoordX(n, nodeItemMap[n]->pos().x());
      _pMolecule->setCoordY(n, nodeItemMap[n]->pos().y());
    }

    _pMolecule->writeLGF(out);

    return true;
  }
}

void GraphWidget::updateNrCorrectCGs()
{
  _nCorrectCGs = 0;

  // let's do a stupid n^2 loop
  for (ChargeGroupVectorIt cgIt1 = _cgVec.begin(); cgIt1 != _cgVec.end(); cgIt1++)
  {
    if (cgIt1->size() == 0)
      continue;

    for (ChargeGroupVectorIt cgIt2 = _cgVecSol.begin(); cgIt2 != _cgVecSol.end(); cgIt2++)
    {
      if (*cgIt1 == *cgIt2)
        _nCorrectCGs++;
    }
  }
}

void GraphWidget::clear()
{
  assert(_pMolecule);
  if (!_showSol)
  {
    _cgVec = ChargeGroupVector(_pMolecule->getNumberOfAtoms(), ChargeGroup());
    init(_cgVec, _pScene, _pNodeItem, _pSelectedNodeItem);
    setScene(_pScene);

    emit updateKey(*_pMolecule,
                   getNodeItemMap(),
                   getChargeGroups(),
                   getSelectedNodeItem(),
                   _optObj,
                   _solved);
    emit updateErrLabel();
    emit updateCorrectCgLabel();
  }
}
