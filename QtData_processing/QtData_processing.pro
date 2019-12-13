TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fPIC
QMAKE_CXXFLAGS += -pthread

QMAKE_LFLAGS += -fPIC
QMAKE_LFLAGS += -pthread
QMAKE_LFLAGS += -pipe

SOURCES += \
        ../source/AllEventsResults.cpp \
        ../source/AnalysisManager.cpp \
        ../source/ExperimentArea.cpp \
        ../source/GlobalParameters.cpp \
        ../source/GraphCollection.cpp \
        ../source/main.cpp \
        ../source/MTAnalysisManager.cpp \
        ../source/PeakInfo.cpp \
        ../source/Polynom2Order.cpp \
        ../source/PolynomialFit.cpp \
        ../source/Savitzky_Golay_filter.cpp \
        ../source/SignalOperations.cpp \
        ../source/SingleEventData.cpp

HEADERS += \
        ../include/AllEventsResults.h \
        ../include/AnalysisManager.h \
        ../include/ExperimentArea.h \
        ../include/GlobalDefinitions.h \
        ../include/GlobalParameters.h \
        ../include/GraphCollection.h \
        ../include/MTAnalysisManager.h \
        ../include/Polynom2Order.h \
        ../include/PolynomialFit.h \
        ../include/Savitzky_Golay_filter.h \
        ../include/SignalOperations.h \
        ../include/SingleEventData.h

unix:!macx: LIBS += -L$(HOME)/Software/root_v6.14.06/build/lib/ \
        -lpthread \
        -lGui \
        -lCore \
        -lRIO \
        -lNet \
        -lHist \
        -lGraf \
        -lGraf3d \
        -lGpad \
        -lTree \
        -lRint \
        -lPostscript \
        -lMatrix \
        -lPhysics \
        -lMathCore \
        -lThread \
        -lMultiProc \
        -lGeom \
        -lm \
        -lSpectrum \
        -lThread

INCLUDEPATH += $(HOME)/Software/root_v6.14.06/build/lib
INCLUDEPATH += $(HOME)/Software/root_v6.14.06/build/include
INCLUDEPATH += $(HOME)/Software/boost_1_67_0
INCLUDEPATH += $(HOME)/Documents/Data_processing/include
DEPENDPATH += $(HOME)/Software/root_v6.14.06/build/lib
DEPENDPATH += $(HOME)/Software/root_v6.14.06/build/include
DEPENDPATH += $(HOME)/Software/boost_1_67_0
