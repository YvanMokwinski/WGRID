PLATFORM=Linux
CC=gcc
CPP=g++

CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=1     -DWMESH_MKL_BLAS
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=1 -DWMESH_MKL_BLAS
CFLAGS+=-std=c99
CPPFLAGS+=-std=c++11

BLASLIB=-L/usr/lib/openblas/lib -lopenblas
BLASLIB=-lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core  -lpthread

