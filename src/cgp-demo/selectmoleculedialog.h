#ifndef SELECTMOLECULEDIALOG_H
#define SELECTMOLECULEDIALOG_H

#include <QString>
#include <QDialog>
#include <vector>

namespace Ui {
class SelectMoleculeDialog;
}

class SelectMoleculeDialog : public QDialog
{
  Q_OBJECT
  
public:
  /// Input triple containing molecule name,
  /// LGF filename and TW filename
  struct InputTriple
  {
    QString _name;
    QString _lgfFile;
    QString _twFile;
  };

  /// Type of a vector of input triples
  typedef std::vector<InputTriple> InputTripleVector;
  /// Type of const iterator of a vector of input triplets
  typedef InputTripleVector::const_iterator InputTripleVectorIt;

  /// Constructor
  ///
  /// \param list are the input triples
  /// \param selectedIdx is currently selected triple
  /// \param pParent is the parent widget
  SelectMoleculeDialog(const InputTripleVector& list,
                       int selectedIdx,
                       QWidget* pParent = 0);
  ~SelectMoleculeDialog();

  /// Called when the dialog is closed
  void done(int r);

private:
  Ui::SelectMoleculeDialog *ui;
};

#endif // SELECTMOLECULEDIALOG_H
