CFL=-g
F77=f77 -v
.f.o: ; $(F77) $(CFL) -c $*.f
.c.o: ; ${CC} ${CFLAGS}  -c $*.c
CFLAGS          =  ${INCLUDE} $(LDFLAGS)

CC = gcc -g -Wall $(CFLAGS)


multitap.o:	multitap.c
mult_tap_spec.o:	mult_tap_spec.c
hires.o:	hires.c
adwait.o:	adwait.c
ftest.o:	ftest.c
drivetap.o:	drivetap.c
sigstuff.o:	sigstuff.c
mtm_getline.o:	mtm_getline.c
jrealft.o:	jrealft.c 
jtridib.o:	jtridib.c
jtinvit.o:	jtinvit.c
nrutil.o:	nrutil.c
jfour1.o:	jfour1.c

TAP=multitap.o mult_tap_spec.o hires.o adwait.o ftest.o sigstuff.o mtm_getline.o jrealft.o jtridib.o  jtinvit.o nrutil.o jfour1.o
testmt:	drivetap.o $(TAP)
	$(CC)  -o   testmt  drivetap.o $(TAP) -lm


clean:
	rm -f testmt *.o core

all:	testmt
