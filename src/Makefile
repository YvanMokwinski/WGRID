VPATH=$(DISTRIB)/src

HEADERS= -I$(DISTRIB)/include -I$(DISTRIB)/src/include


source_lib=$(wildcard $(VPATH)/*.cpp)
OBJ=$(notdir $(source_lib:%.cpp=%.o))

include $(DISTRIB)/config.mk

install: $(OBJ) libwmesh.a

%.o:%.c
	$(CC) $(CFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 
 
%.o:%.cpp
	$(CPP) $(CPPFLAGS) $(FLAGS) $(HEADERS) -c  $<   -o $@ 

libwmesh.a:$(OBJ)
	\rm -f $@
	ar crsu $@ $(OBJ)
	ranlib $@
	echo "'$(PLATFORM)_$(CC)'"
	ln -fs $(DISTRIB)/$(PLATFORM)-$(CC)$(DEBUG)/$@  $(DISTRIB)/lib/libwmesh-$(CC)$(DEBUG).a
