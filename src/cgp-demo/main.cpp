/* main.cpp
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>
#include <lemon/list_graph.h>
#include <QtWidgets>

#include "mainwindow.h"
#include "moleculeviz.h"
#include "common/verbose.h"

typedef lemon::ListGraph Graph;

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

#ifdef ANDROID
  QString filename = "/sdcard/CWI/list.dat";
#else
  lemon::ArgParser ap(argc, argv);
  ap.parse();

  QString filename;
  if (ap.files().size() == 0)
  {
    filename = QFileDialog::getOpenFileName(NULL, "Open file", "", "*.dat");
  }
  else
  {
    filename = ap.files()[0].c_str();
  }
#endif

  nina::g_verbosity = nina::VERBOSE_NONE;

  MainWindow mainWindow;
  if (!mainWindow.parse(filename))
    return 1;

#ifdef ANDROID
  mainWindow.setFocus();
  mainWindow.showFullScreen();
#else
  mainWindow.showNormal();
#endif

  return app.exec();
}
