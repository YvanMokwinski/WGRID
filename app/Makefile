VPATH=$(DISTRIB)/app

INCLUDES= -DWMESH_ILP64=0  -I$(DISTRIB)/include -I$(DISTRIB)/../WCOMMON/include
LIBS=-L$(DISTRIB)/lib -lwmesh-$(CC)$(DEBUG) -lMeshb.7 -L/usr/lib/openblas/lib -lopenblas #  -lmkl_intel_ilp64 -lmkl_intel_ilp64 -lmkl_lapack95_ilp64  -lmkl_sequential -lmkl_core   -lpthread

cpp_source_lib=$(wildcard $(VPATH)/*.cpp)
CPP_OBJ=$(notdir $(cpp_source_lib:%.cpp=%.cpp.o))
CPP_EXE=$(notdir $(cpp_source_lib:%.cpp=%.cpp.exe))

include $(DISTRIB)/config.mk

install: $(CPP_OBJ) $(CPP_EXE)

%.cpp.exe:%.cpp.o libwmesh.a
	$(CPP) $(CPPFLAGS) $< -o $@  $(LIBS)
	ln -fs $(DISTRIB)/$(PLATFORM)-$(CC)$(DEBUG)/$@ $(DISTRIB)/bin/$*$(DEBUG)

%.cpp.o:%.cpp
	$(CPP) $(CPPFLAGS) -c  $< -o $@ $(INCLUDES)

