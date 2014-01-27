/* main.cpp
 *
 *  Created on: 31-Jul-2012
 *      Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/list_graph.h>
#include <QtWidgets>

#include "mainwindow.h"
#include "moleculeviz.h"
#include "verbose.h"

typedef lemon::ListGraph Graph;

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

#ifdef ANDROID
  QString filename = "/sdcard/CWI/list.dat";
#else
  QString filename;
  if (argc < 2)
  {
    filename = QFileDialog::getOpenFileName(NULL, "Open file", "", "*.dat");
  }
  else
  {
    filename = argv[1];
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
