PLATFORM=Linux
CC=gcc
CPP=g++

CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=0     -DWMESH_OPEN_BLAS
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=0 -DWMESH_OPEN_BLAS
BLASLIB=-L/usr/lib/openblas/lib -lopenbla

CFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=0     -DWMESH_OPEN_BLAS -fno-stack-protector
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=0 -DWMESH_OPEN_BLAS -fno-stack-protector
BLASLIB=-L/usr/lib/ -llapack -lblas -lgfortran


CFLAGS= -pipe -m64 -Wall -Werror  -funroll-loops -funroll-all-loops -DWMESH_ILP64=1     -DWMESH_MKL_BLAS  -I/home/mokwinski/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/
CPPFLAGS= -pipe -m64 -Wall -Werror -funroll-loops -funroll-all-loops -DWMESH_ILP64=1 -DWMESH_MKL_BLAS  -I/home/mokwinski/intel/compilers_and_libraries_2020.1.217/linux/mkl/include/
BLASLIB= -L/home/mokwinski/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core  -lpthread

CFLAGS+=-std=c99
CPPFLAGS+=-std=c++11


