TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -llapack -lblas -larmadillo

INCLUDEPATH += "/home/gullik/Test/gnuplot-cpp"
INCLUDEPATH += "/home/gullik/Downloads/compphys/programs/cppLibrary"
