VPATH=$(DISTRIB)/app

INCLUDES=-I$(DISTRIB)/include -I$(DISTRIB)/src/include
cpp_source_lib=$(wildcard $(VPATH)/*.cpp)
CPP_OBJ=$(notdir $(cpp_source_lib:%.cpp=%.cpp.o))
CPP_EXE=$(notdir $(cpp_source_lib:%.cpp=%.cpp.exe))

include $(DISTRIB)/config.mk

LIBS=-L$(DISTRIB)/lib -lwmesh-$(CC)$(DEBUG) -lMeshb.7 $(BLASLIB) -lm

install: $(CPP_OBJ) $(CPP_EXE)

%.cpp.exe:%.cpp.o libwmesh.a
	$(CPP) $(CPPFLAGS) $(FLAGS) $< -o $@  $(LIBS)
	ln -fs $(DISTRIB)/$(PLATFORM)-$(CC)$(DEBUG)/$@ $(DISTRIB)/bin/$*$(DEBUG)

%.cpp.o:%.cpp
	$(CPP) $(CPPFLAGS) $(FLAGS) -c  $< -o $@ $(INCLUDES)

