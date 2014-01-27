#include <QDir>
#include <iostream>

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "selectmoleculedialog.h"

MainWindow::MainWindow(QWidget* pParent)
  : QMainWindow(pParent)
  , ui(new Ui::MainWindow)
  , _list()
  , _selectedIdx(0)
  , _pTimer(new QTimer(this))
  , _started(false)
{
  ui->setupUi(this);
  connect(ui->graphWidget, SIGNAL(updateKey(MoleculeVizType,NodeItemMap,ChargeGroupVector,const NodeItem*,double,bool&)),
          ui->keyWidget, SLOT(update(MoleculeVizType,NodeItemMap,ChargeGroupVector,const NodeItem*,double,bool&)));
  connect(ui->graphWidget, SIGNAL(solved()),this,SLOT(solved()));
  connect(ui->graphWidget, SIGNAL(updateErrLabel()), this, SLOT(updateErrLabel()));
  connect(ui->graphWidget, SIGNAL(updateCorrectCgLabel()), this, SLOT(updateNrCorrectGroupsLabel()));
  connect(ui->pushButtonClear, SIGNAL(clicked()), ui->graphWidget, SLOT(clear()));
  connect(ui->pushButtonClear, SIGNAL(clicked()), this, SLOT(nextLevel()));
  connect(ui->pushButtonStop, SIGNAL(clicked()), this, SLOT(toggleStartStop()));
  ui->graphWidget->setFocus();

  // some keyboard shortcuts
  QShortcut* pExport = new QShortcut(QKeySequence(Qt::Key_F12), this);
  pExport->setContext(Qt::ApplicationShortcut);

  QShortcut* pFullscreen = new QShortcut(QKeySequence(Qt::Key_F11), this);
  pFullscreen->setContext(Qt::ApplicationShortcut);

  QShortcut* pQuit = new QShortcut(QKeySequence(Qt::Key_Q), this);
  pQuit->setContext(Qt::ApplicationShortcut);

  QShortcut* pIncK1 = new QShortcut(QKeySequence("Ctrl++"), this);
  QShortcut* pIncK2 = new QShortcut(QKeySequence("Ctrl+="), this);
  pIncK1->setContext(Qt::ApplicationShortcut);
  pIncK2->setContext(Qt::ApplicationShortcut);

  QShortcut* pDecK = new QShortcut(QKeySequence("Ctrl+-"), this);
  pDecK->setContext(Qt::ApplicationShortcut);

  QShortcut* pLoad = new QShortcut(QKeySequence("Ctrl+l"), this);
  pLoad->setContext(Qt::ApplicationShortcut);

  QShortcut* pShowSolution = new QShortcut(QKeySequence(Qt::Key_F10), this);
  pShowSolution->setContext(Qt::ApplicationShortcut);

  connect(pExport, SIGNAL(activated()), this, SLOT(writeLGF()));
  connect(pFullscreen, SIGNAL(activated()), this, SLOT(toggleFullScreen()));
  connect(pQuit, SIGNAL(activated()), this, SLOT(close()));
  connect(pIncK1, SIGNAL(activated()), this, SLOT(incK()));
  connect(pIncK2, SIGNAL(activated()), this, SLOT(incK()));
  connect(pDecK, SIGNAL(activated()), this, SLOT(decK()));
  connect(pLoad, SIGNAL(activated()), this, SLOT(load()));
  connect(pShowSolution, SIGNAL(activated()), this, SLOT(toggleSolution()));
  connect(_pTimer, SIGNAL(timeout()), this, SLOT(updateTime()));

  showRules(true);
}

MainWindow::~MainWindow()
{
  delete ui;
  delete _pTimer;
}

void MainWindow::showRules(bool start)
{
  QString rules("<b>Rules of the game</b>\n"
                "<ul>\n"
                "<li>Each atom has a charge</li>\n"
                "<li>Form groups of atoms with total charges close to 0</li>\n"
                "<li>Groups have size at most 5</li>\n"
                "<li>You can create a new group by tapping a node that is not adjacent to the current group</li>\n"
                "</ul>"
                "To see how good you are, check the number of correct groups in the lower left corner.\n\n"
                "Be quick, there is a timer running!");

  if (start)
  {
    rules += "<br><br>Press start to begin.";
  }
  else
  {
    rules += "<br>";
  }

  ui->textBrowser->setText(rules);
}

void MainWindow::toggleStartStop()
{
  if (_started)
    stop();
  else
    start();
}

void MainWindow::start()
{
  _started = true;

  _pTimer->start(1000);
  ui->lcdNumber->display(0);
  ui->pushButtonStop->setText("Stop");
  ui->pushButtonClear->setEnabled(true);
  ui->graphWidget->setEnabled(true);

  loadLevel(0);
  showRules(false);
}

void MainWindow::stop()
{
  _started = false;

  _pTimer->stop();
  ui->pushButtonStop->setText("Start");
  ui->pushButtonClear->setEnabled(false);
  ui->graphWidget->scene()->clear();
  ui->keyWidget->scene()->clear();
  ui->keyWidget->repaint();

  ui->labelMolecule->setText("Press start to begin");
  ui->labelCorrectGroups->setText("");
  ui->labelError->setText("");
  ui->lcdNumber->display(0);

  showRules(true);
  _selectedIdx = 0;
}

void MainWindow::updateErrLabel()
{
  ui->labelError->setText(QString("%1/%2")
                          .arg(ui->keyWidget->getCurObj(), 0, 'f', 3)
                          .arg(ui->graphWidget->getOptObj(), 0, 'f', 3));

}

void MainWindow::updateNrCorrectGroupsLabel()
{
  ui->labelCorrectGroups->setText(QString("%1")
                                  .arg(ui->graphWidget->getNumberOfCorrectChargeGroups()));
}

void MainWindow::incK()
{
  GraphWidget::MoleculeVizType& molecule = ui->graphWidget->getMolecule();
  size_t k = molecule.getMaxChargeGroupSize();

  if (k == molecule.getNumberOfAtoms())
    return;

  molecule.setMaxChargeGroupSize(++k);

  ui->graphWidget->initSolution(_list[_selectedIdx]._twFile);

  ui->labelMolecule->setText(QString("%1, k = %2")
                             .arg(_list[_selectedIdx]._name)
                             .arg(k));
}

void MainWindow::decK()
{
  // clear charge groups
  GraphWidget::MoleculeVizType& molecule = ui->graphWidget->getMolecule();
  size_t k = molecule.getMaxChargeGroupSize();

  if (k == 2)
    return;

  molecule.setMaxChargeGroupSize(--k);

  bool dummy = false;
  load(_list[_selectedIdx]);
  ui->keyWidget->update(molecule,
                        ui->graphWidget->getNodeItemMap(),
                        ui->graphWidget->getChargeGroups(),
                        ui->graphWidget->getSelectedNodeItem(),
                        ui->graphWidget->getOptObj(),
                        dummy);
  ui->labelMolecule->setText(QString("%1, k = %2")
                             .arg(_list[_selectedIdx]._name)
                             .arg(k));
}

void MainWindow::load()
{
  SelectMoleculeDialog sel(_list, _selectedIdx, this);
  sel.exec();

  if (sel.result() != 0)
  {
    _selectedIdx = sel.result() - 1;
    load(_list[_selectedIdx]);
  }
}

bool MainWindow::parse(const QString& listFilename)
{
  std::ifstream in(listFilename.toLatin1().data());
  if (!in.good())
  {
#ifndef ANDROID
    QMessageBox msgBox(this);
    msgBox.setText(QString("Failed to open '%1'").arg(listFilename));
    msgBox.setIcon(QMessageBox::Critical);
    msgBox.exec();
#endif

    return false;
  }

  _list.clear();
  while (in.good())
  {
    SelectMoleculeDialog::InputTriple triple;
    std::string name, lgfFile, twFile;
    while (name.empty() && in.good())
      std::getline(in, name);
    while (lgfFile.empty() && in.good())
      std::getline(in, lgfFile);
    while (twFile.empty() && in.good())
      std::getline(in, twFile);

    triple._name = name.c_str();
    triple._lgfFile = lgfFile.c_str();
    triple._twFile = twFile.c_str();

    _list.push_back(triple);
  }

  if (_list.size() == 0)
    return false;

  QFileInfo fileInfo(listFilename);
  QDir::setCurrent(fileInfo.absolutePath());

  return true;
}

bool MainWindow::load(const SelectMoleculeDialog::InputTriple& triple)
{
  if (ui->graphWidget->parseLGF(triple._lgfFile))
  {
    const GraphWidget::MoleculeVizType& molecule = ui->graphWidget->getMolecule();
    ui->labelMolecule->setText(QString("%1, k = %2")
                               .arg(triple._name)
                               .arg(molecule.getMaxChargeGroupSize()));

    return ui->graphWidget->initSolution(triple._twFile);
  }

  return false;
}

void MainWindow::toggleSolution()
{
  ui->graphWidget->toggleScene();
}

void MainWindow::writeLGF()
{
  QString name = _list[_selectedIdx]._lgfFile;
  QString path = QDir::tempPath() + QDir::separator() + name;
  ui->graphWidget->writeLGF(path);
}

void MainWindow::toggleFullScreen()
{
  if (this->isFullScreen())
    this->showNormal();
  else
    this->showFullScreen();
}

void MainWindow::updateTime()
{
  ui->lcdNumber->display(ui->lcdNumber->intValue() + 1);
}

bool MainWindow::loadLevel(int level)
{
  assert(0 <= level && level < static_cast<int>(_list.size()));

  // reset counter
  ui->lcdNumber->display(0);
  ui->labelCorrectGroups->setText(0);
  return load(_list[level]);
}

void MainWindow::solved()
{
  _pTimer->stop();

  if (_selectedIdx + 1 != static_cast<int>(_list.size()))
  {
    ui->textBrowser->setText(ui->textBrowser->toHtml() +=
        QString("<font color=\"red\">Nice! Finished level %2 in %1 seconds. Press next level to continue.</font>").arg(ui->lcdNumber->value()).arg(_selectedIdx + 1));
    ui->pushButtonClear->setText("Next level");
    ui->graphWidget->setEnabled(false);
  }
  else
  {
    ui->textBrowser->setText(ui->textBrowser->toHtml() +=
        QString("<font color=\"red\">Nice! Finished level %2 in %1 seconds.</font>").arg(ui->lcdNumber->value()).arg(_selectedIdx + 1));
    ui->graphWidget->setEnabled(false);
  }
}

void MainWindow::nextLevel()
{
  if (ui->pushButtonClear->text() == "Next level")
  {
    ui->graphWidget->setEnabled(true);
    ui->pushButtonClear->setText("Clear");
    _pTimer->start(1000);
    if (_selectedIdx + 1 != static_cast<int>(_list.size()))
    {
      _selectedIdx++;
      loadLevel(_selectedIdx);
    }
  }
}

void MainWindow::showExpanded()
{
#if defined(Q_WS_SIMULATOR)
    showFullScreen();
#elif defined(Q_WS_MAEMO_5)
    showMaximized();
#else
    show();
#endif
}
