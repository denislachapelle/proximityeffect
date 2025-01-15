#
# DL241118
# this file compile all programs developped while experimenting with mfem.
#

MFEM_DIR ?= /mnt/c/mfem-4.7
MFEM_BUILD_DIR ?= /mnt/c/mfem-4.7

GMSH_DIR ?= /
#MFEM_DIR ?= /home/denislachapelle2003/fem/mfem-4.6
#MFEM_BUILD_DIR ?= /home/denislachapelle2003/fem/mfem-4.6

#COMMON_LIB = -L$(MFEM_BUILD_DIR)/miniapps/common -lmfem-common

all: proximity2dblock proximity2dblockb se2dblock proximity2dblockc proximity2dblockd

proximity2dblock: proximity2dblock.cpp mytools.cpp mytools.hpp
	g++ -g -o proximity2dblock  -std=c++11 -I$(MFEM_DIR) proximity2dblock.cpp mytools.cpp  -L$(MFEM_BUILD_DIR) -lmfem -lrt

proximity2dblockb: proximity2dblockb.cpp mytools.cpp mytools.hpp
	g++ -g -o proximity2dblockb  -std=c++11 -I$(MFEM_DIR) proximity2dblockb.cpp mytools.cpp  -L$(MFEM_BUILD_DIR) -lmfem -lrt

se2dblock: se2dblock.cpp mytools.cpp mytools.hpp
	g++ -g -o se2dblock  -std=c++11 -I$(MFEM_DIR) se2dblock.cpp mytools.cpp  -L$(MFEM_BUILD_DIR) -lmfem -lrt


proximity2dblockd: proximity2dblockd.cpp mytools.cpp mytools.hpp
	g++ -g -o proximity2dblockd  -std=c++11 -I$(MFEM_DIR) proximity2dblockd.cpp mytools.cpp  -I/usr/local/include -L/usr/lib/libgmsh.so -lgmsh -L$(MFEM_BUILD_DIR) -lmfem -lrt 

