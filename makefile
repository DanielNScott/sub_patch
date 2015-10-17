FCOMP = gfortran-4.9
FCFLAGS = -g -fbounds-check -O3
FCFLAGS += -I/usr/include

all: main

main.o: ed_state_vars.o
main: ed_state_vars.o

%: %.o ; ${FCOMP} ${FCFLAGS} -o $@ $^ ${LDFLAGS}

%.o: %.f08 ; $(FCOMP) $(FCFLAGS) -c $<

.PHONY: clean

clean: ; rm -f *.o *.mod *.MOD
