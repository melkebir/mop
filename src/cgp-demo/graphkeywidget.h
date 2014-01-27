/* graphkeywidget.h
 *
 *  Created on: 14-Aug-2012
 *      Author: M. El-Kebir
 */

#ifndef GRAPHKEYWIDGET_H
#define GRAPHKEYWIDGET_H

#include <QtWidgets/QGraphicsView>
#include <QtWidgets>
#include <lemon/list_graph.h>

#include "moleculeviz.h"
#include "nodekeyitem.h"

class GraphKeyWidget : public QGraphicsView
{
  Q_OBJECT

public:
  /// Graph type
  typedef lemon::ListGraph Graph;
  /// Molecule type
  typedef MoleculeViz<Graph> MoleculeVizType;
  /// Charge group type
  typedef MoleculeVizType::ChargeGroup ChargeGroup;
  /// Charge group iterator
  typedef MoleculeVizType::ChargeGroupIt ChargeGroupIt;
  /// Charge group vector type
  typedef MoleculeVizType::ChargeGroupVector ChargeGroupVector;
  /// Charge group vector iterator
  typedef MoleculeVizType::ChargeGroupVectorIt ChargeGroupVectorIt;
  /// Weights on nodes
  typedef MoleculeVizType::WeightNodeMap WeightNodeMap;
  /// NodeItem map type
  typedef Graph::NodeMap<NodeItem*> NodeItemMap;

  /// Constructor
  ///
  /// \param pParent is the parent of this object
  GraphKeyWidget(QWidget* pParent = 0);

  /// Returns the (current) objective value
  double getCurObj() const;

private:
  /// LEMON graph typedefs
  GRAPH_TYPEDEFS(Graph);

protected:
  /// Rescales the whole thing upon a resize event
  void resizeEvent(QResizeEvent* pEvent);

  /// Used for sorting the charge groups based on their absolute error
  struct CgComp
  {
    private:
      const MoleculeVizType& _molecule;

    public:
      CgComp(const MoleculeVizType& molecule)
        : _molecule(molecule)
      {
      }

      bool operator() (const ChargeGroup& a, const ChargeGroup& b)
      {
        return _molecule.computeError(a) > _molecule.computeError(b);
      }
  };

private:
  /// Node graphics item map type
  typedef Graph::NodeMap<NodeKeyItem*> NodeKeyItemMap;
  /// Node item map
  NodeKeyItemMap* _pNodeKeyItemMap;
  /// Objective value
  double _curObj;

public slots:
  /// Updates the key
  ///
  /// \param molecule is the visualized molecule
  /// \param nodeItem is the node item map, mapping LEMON nodes
  ///        to graphics node items
  /// \param cgVec are the charge groups
  /// \param pNodeItem is the currently selected node graphics item
  /// \param optObj is the optimal objective value
  /// \param solved will contain whether we are done
  void update(const MoleculeVizType& molecule,
              const NodeItemMap& nodeItem,
              const ChargeGroupVector& cgVec,
              const NodeItem* pSelectedNodeItem,
              double optObj,
              bool& solved);
};

inline double GraphKeyWidget::getCurObj() const
{
  return _curObj;
}

#endif // GRAPHKEYWIDGET_H
