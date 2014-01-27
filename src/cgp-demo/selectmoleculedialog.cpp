#include "selectmoleculedialog.h"
#include "ui_selectmoleculedialog.h"
#include <iostream>

SelectMoleculeDialog::SelectMoleculeDialog(const InputTripleVector& list,
                                           int selectedIdx,
                                           QWidget *parent)
  : QDialog(parent)
  , ui(new Ui::SelectMoleculeDialog)
{
  ui->setupUi(this);

  // initialize list
  int n = 0;
  for (InputTripleVectorIt it = list.begin(); it != list.end(); it++)
  {
    new QListWidgetItem(it->_name, ui->listWidget);
    n++;
  }

  if (0 <= selectedIdx && selectedIdx < n)
    ui->listWidget->setCurrentRow(selectedIdx);
}

SelectMoleculeDialog::~SelectMoleculeDialog()
{
  delete ui;
}

void SelectMoleculeDialog::done(int r)
{
  // result = 0 : cancel
  // result > 0 : load with index result
  r = r != 0 ? ui->listWidget->currentRow() + 1 : 0;
  QDialog::done(r);
}
