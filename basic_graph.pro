QMAKE_CXXFLAGS += -Werror
HEADERS       = window.h \
		functions.h \
		approx2.h

SOURCES       = main.cpp \
                window.cpp \
		functions2.cpp \
		approx2.cpp

QMAKE_CXXFLAGS += -std=c++11
LIBS += -lm
QT += widgets
