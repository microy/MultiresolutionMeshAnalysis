###########################################################################
#                                                                         #
#                                                                         #
#                                                                         #
#                                  Makefile                               #
#                             -------------------                         #
#    begin                : 2009/02/19                                    #
#    copyright            : (C) 2003-2009 by Michaël Roy                  #
#    email                : michael.roy@u-bourgogne.fr                    #
#                                                                         #
#                                                                         #
#                                                                         #
###########################################################################








###########################################################################
#                                                                         #
#                                                                         #
#                                                                         #

COMPILER  = g++

BIN       = mma

LIBS      = -lGL -lGLU -lfltk -lfltk_gl -ljpeg
LIBDIR    = -L/usr/X11R6/lib -L/usr/lib
INCDIR    = -I/usr/X11R6/include -I/usr/include

# Release flags
FLAGS     = -s -Wall -W -ansi -pedantic -DNDEBUG -O3 $(INCDIR)
# Debug flags
#FLAGS     = -g -Wall -W -ansi -pedantic $(INCDIR)

#                                                                         #
#                                                                         #
#                                                                         #
###########################################################################


#
# generating a list of the object files
#

OBJS = $(notdir $(patsubst %.cpp,%.o,$(wildcard *.cpp)))


#
# defining an explicit rule for object file dependencies
#

%.o : %.cpp
	$(COMPILER) $(FLAGS) -o $(notdir $@) -c $<


#
# targets which don't actually refer to files
#

.PHONY : all clean


###########################################################################
#                                                                         #
#                                                                         #
#                                                                         #

all: $(OBJS)
	$(COMPILER) $(OBJS) -o $(BIN) $(FLAGS) $(LIBDIR) $(LIBS)

clean:
	-rm -f *.o *~ $(BIN)

#                                                                         #
#                                                                         #
###########################################################################

