QT += core
QT -= gui

CONFIG += c++11

TARGET = RIMD
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

DEFINES += _USE_MATH_DEFINES

SOURCES += main.cpp \
    util_3drotation_log_exp.cpp \
    rimd_reconstruction.cpp

HEADERS += \
    util_3drotation_log_exp.h \
    rimd_reconstruction.h \
    trimesh_types_hh.h \

INCLUDEPATH += /home/chern/Downloads/OpenMesh-6.3/build/Build/src \

LIBS+= -L/usr/local/lib \

LIBS += -lOpenMeshCore \
        -lOpenMeshTools \

#LIBS += -lOpenMeshCored \
#        -lOpenMeshToolsd \
