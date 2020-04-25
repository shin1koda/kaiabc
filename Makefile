FC = ifort

TARGET = example
OBJS = model.o main.o

FFLAGS += -O3 -qopenmp
LDLIBS += -mkl=parallel

.SUFFIXES: .f90 .mod

vpath %.f90 src

%.o: %.f90
	$(COMPILE.f) $<

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(LINK.f) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

main.o: model.mod
