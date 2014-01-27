#CONFIG += graphviz
#CONFIG += static

QT += widgets

HEADERS += \
    graphwidget.h \
    nodeitem.h \
    moleculeviz.h \
    edgeitem.h \
    mainwindow.h \
    nodekeyitem.h \
    errorkeyitem.h \
    selectmoleculedialog.h \
    graphkeywidget.h

SOURCES += \
    main.cpp \
    ../../lemon/arg_parser.cc \
    ../../common/verbose.cpp \
    nodekeyitem.cpp \
    nodeitem.cpp \
    errorkeyitem.cpp \
    graphwidget.cpp \
    edgeitem.cpp \
    mainwindow.cpp \
    graphkeywidget.cpp \
    selectmoleculedialog.cpp

INCLUDEPATH += ../../lemon/include \
        ../../ \
        ../

FORMS += \
    mainwindow.ui \
    selectmoleculedialog.ui

# Please do not modify the following two lines. Required for deployment.
include(deployment.pri)
qtcAddDeployment()

graphviz {
    HEADERS += moleculegraphviz.h
    CONFIG += link_pkgconfig
    PKGCONFIG += libgvc
    DEFINES += GRAPHVIZ
}

OTHER_FILES += \
    android/res/values-fa/strings.xml \
    android/res/values-ro/strings.xml \
    android/res/values-nl/strings.xml \
    android/res/values/libs.xml \
    android/res/values/strings.xml \
    android/res/drawable/logo.png \
    android/res/drawable/icon.png \
    android/res/values-pl/strings.xml \
    android/res/values-ja/strings.xml \
    android/res/values-rs/strings.xml \
    android/res/layout/splash.xml \
    android/res/values-it/strings.xml \
    android/res/values-et/strings.xml \
    android/res/drawable-mdpi/icon.png \
    android/res/values-ms/strings.xml \
    android/res/drawable-hdpi/icon.png \
    android/res/values-el/strings.xml \
    android/res/values-zh-rCN/strings.xml \
    android/res/values-pt-rBR/strings.xml \
    android/res/values-fr/strings.xml \
    android/res/drawable-ldpi/icon.png \
    android/res/values-id/strings.xml \
    android/res/values-nb/strings.xml \
    android/res/values-es/strings.xml \
    android/res/values-ru/strings.xml \
    android/res/values-de/strings.xml \
    android/res/values-zh-rTW/strings.xml \
    android/src/org/kde/necessitas/origo/QtActivity.java \
    android/src/org/kde/necessitas/origo/QtApplication.java \
    android/src/org/kde/necessitas/ministro/IMinistro.aidl \
    android/src/org/kde/necessitas/ministro/IMinistroCallback.aidl \
    android/version.xml \
    android/AndroidManifest.xml \ 
    android/src/org/qtproject/qt5/android/bindings/QtActivity.java \
    android/src/org/qtproject/qt5/android/bindings/QtApplication.java \
    android/src/org/kde/necessitas/ministro/IMinistro.aidl \
    android/src/org/kde/necessitas/ministro/IMinistroCallback.aidl \
    android/AndroidManifest.xml \
    android/res/values/libs.xml \
    android/res/values/strings.xml \
    android/res/layout/splash.xml \
    android/res/values-et/strings.xml \
    android/res/values-zh-rCN/strings.xml \
    android/res/values-zh-rTW/strings.xml \
    android/res/values-es/strings.xml \
    android/res/values-nb/strings.xml \
    android/res/values-ms/strings.xml \
    android/res/values-it/strings.xml \
    android/res/values-pl/strings.xml \
    android/res/values-id/strings.xml \
    android/res/values-pt-rBR/strings.xml \
    android/res/values-rs/strings.xml \
    android/res/values-fr/strings.xml \
    android/res/values-ro/strings.xml \
    android/res/values-nl/strings.xml \
    android/res/values-de/strings.xml \
    android/res/values-ja/strings.xml \
    android/res/values-ru/strings.xml \
    android/res/values-fa/strings.xml \
    android/res/values-el/strings.xml \
    android/version.xml \
    android/src/org/qtproject/qt5/android/bindings/QtActivity.java \
    android/src/org/qtproject/qt5/android/bindings/QtApplication.java \
    android/src/org/kde/necessitas/ministro/IMinistro.aidl \
    android/src/org/kde/necessitas/ministro/IMinistroCallback.aidl \
    android/AndroidManifest.xml \
    android/res/values/libs.xml \
    android/res/values/strings.xml \
    android/res/layout/splash.xml \
    android/res/values-et/strings.xml \
    android/res/values-zh-rCN/strings.xml \
    android/res/values-zh-rTW/strings.xml \
    android/res/values-es/strings.xml \
    android/res/values-nb/strings.xml \
    android/res/values-ms/strings.xml \
    android/res/values-it/strings.xml \
    android/res/values-pl/strings.xml \
    android/res/values-id/strings.xml \
    android/res/values-pt-rBR/strings.xml \
    android/res/values-rs/strings.xml \
    android/res/values-fr/strings.xml \
    android/res/values-ro/strings.xml \
    android/res/values-nl/strings.xml \
    android/res/values-de/strings.xml \
    android/res/values-ja/strings.xml \
    android/res/values-ru/strings.xml \
    android/res/values-fa/strings.xml \
    android/res/values-el/strings.xml \
    android/version.xml
