TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    JacobiRot.cpp

HEADERS += \
    Proj2Lib.h

LIBS += -llapack -lblas -larmadillo

INCLUDEPATH += "/home/gullik/Documents/Fys3150Gullik/Project_2"
