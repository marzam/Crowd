    ################################################################################
#

#
#  Universidade Federal Fluminense - Medialab
#  por: Marcelo Zamith - mzamith@ic.uff.br
#  para compilar com target release use o comando: make TARGET=release
#  -framework Cocoa
#
################################################################################
EXEFILE     = CA
VERSION     = -D_VERSION=\"0.01\"
APPLICATION = -D_APPLICATION=\"$(EXEFILE)\"
CPUCC     = g++
CPPFLAGS  = -g -std=c++11
DEFS      =  -DALIGN=64

INCLUDES  = -I.                           \
            -I/opt/glew/include \
            -I/usr/include/GL

LIBDIR   =  -L/usr/lib                          \
            -L/opt/glew/lib


LIBS     =   -lm  -lGLEW -lGLU -lglut -lGL

LINK     =  $(LIBDIR) $(LIBS)

C_COMPILE = $(CPUCC) $(DEFS) $(INCLUDES) $(CPPFLAGS)


all:   Crowd main
	$(C_COMPILE) Crowd.o main.o    $(LINK) -o $(EXEFILE)

main:
	$(C_COMPILE) -c main.cpp

Crowd:
	$(C_COMPILE) -c Crowd.cpp


clean:
	rm *.o; rm $(EXEFILE)
