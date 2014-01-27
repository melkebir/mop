#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "selectmoleculedialog.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
    enum ScreenOrientation {
        ScreenOrientationLockPortrait,
        ScreenOrientationLockLandscape,
        ScreenOrientationAuto
    };
  /// Constructor
  ///
  /// \param pParent is the parent widget
  explicit MainWindow(QWidget* pParent = 0);
  ~MainWindow();

public:
  /// Parses a file containing input triples
  ///
  /// \param listFilename contains input triples
  /// \return \c true if the parsing is successfully completed
  bool parse(const QString& listFilename);

public slots:
  /// Increments the maximum charge group size
  void incK();
  /// Decrements the maximum charge group size
  void decK();
  /// Shows the load molecule dialog
  void load();
  /// Updates the objective function label
  void updateErrLabel();
  /// Updates the number of correct groups label
  void updateNrCorrectGroupsLabel();
  /// Toggles between solution/normal scene
  void toggleSolution();
  /// Export current molecule in LGF format
  void writeLGF();
  /// Toggle fullscreen
  void toggleFullScreen();

  /// Start/stop the game
  void toggleStartStop();

  /// Level solved
  void solved();
  /// Some android stuff
  void showExpanded();
  /// Next level
  void nextLevel();

  /// Show rules
  void showRules(bool start);

private:
  /// User interface
  Ui::MainWindow *ui;
  /// Input triples
  SelectMoleculeDialog::InputTripleVector _list;
  /// Selected input triple
  int _selectedIdx;
  /// Timer
  QTimer* _pTimer;
  /// Started
  bool _started;

  /// Loads specified triple
  ///
  /// \param triple is the specified input triple
  bool load(const SelectMoleculeDialog::InputTriple& triple);

  /// Load specified level
  ///
  /// \param level is the specified level
  bool loadLevel(int level);

  /// Start the game
  void start();
  /// Stop the game
  void stop();

private slots:
  void updateTime();
};

#endif // MAINWINDOW_H
