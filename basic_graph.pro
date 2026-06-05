QMAKE_CXXFLAGS += -Werror
HEADERS       = window.h \
		functions.h \
		method33.h \
		method43.h \

SOURCES       = main.cpp \
                window.cpp \
		functions2.cpp \
		method33.cpp \
		method43.cpp \

QMAKE_CXXFLAGS += -std=c++11
LIBS += -lm
QT += widgets
