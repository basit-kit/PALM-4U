#--------------------------------------------------------------------------------#
# This file is part of PALM.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2014  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------#
#
# Current revisions:
# -----------------
#
# Former revisions:
# -----------------
# $Id: Makefile 1310 2014-03-14 08:01:56Z raasch $
#
# 1046 2012-11-09 14:38:45Z maronga
# code put under GPL (PALM 3.9)
#
# utility programs are created separately for each block in the configuration
# file,
# compare_palm_logs added
# added PROG3 in 'all:' statement (line 39, bugfix)
#
# Initial revision somewhere in 2003/04
#
# Description:
# ------------
# Makefile for generating the utility programs needed by mrun and palm
#--------------------------------------------------------------------------------#

PROG1 =  ../SCRIPTS/combine_plot_fields$(BLOCK).x
PROG2 =  ../SCRIPTS/compare_palm_logs$(BLOCK).x
PROG3 =  ../SCRIPTS/interpret_config$(BLOCK).x

OBJS1 =  combine_plot_fields.o
OBJS2 =  compare_palm_logs.o
OBJS3 =  interpret_config.o

CC = cc
CFLAGS = -O

F90 =
F90_SER =
COPT =
F90FLAGS =
LDFLAGS =

.SUFFIXES: $(SUFFIXES) .f90


all: $(PROG1) $(PROG2) $(PROG3)

$(PROG1): $(OBJS1)
	$(F90_SER) -o $(PROG1) $(OBJS1) $(LDFLAGS)

$(PROG2): $(OBJS2)
	$(F90_SER) -o $(PROG2) $(OBJS2) $(LDFLAGS)

$(PROG3): $(OBJS3)
	$(F90_SER) -o $(PROG3) $(OBJS3) $(LDFLAGS)

.f90.o:
	$(F90_SER) $(F90FLAGS) $(COPT) -c $<


