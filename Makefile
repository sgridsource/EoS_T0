# EoS_T0/Makefile

NAME := EoS_T0
OBJS := sgrid_$(NAME).o EoS_T0.o PwP.o tab1d_Of_rho0_AtT0.o

include $(TOP)/Makefile.subdirs
