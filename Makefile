f90=pgf90
objects=inpred.o mytecio.o subcode0.o subcode1.o \
        subcode2.o subcode3.o subcode4.o subcode5.o \
        subcode6.o subcode7.o subcode8.o supplement01.o \
        supplement02.o dislin.o
modfile=constants.mod  fido.mod  funs.mod  \
        mmytecplt.mod paraset1.mod paraset3.mod  \
        errcontrol.mod files.mod lineread.mod  \
        mtecvar.mod paraset2.mod paraset4.mod heart.mod \
        dislin.mod mplot.mod min4deck.mod paraset5.mod
modules=subcode0.o mytecio.o dislin.o
# compile option for pgf90 or ifort
copt=-fast
# for gfortran
# copt=-O3 --ffree-line-length-0
inpred: $(objects) 
	$(f90) $(copt)  -o penmshxp \
        $(objects) dislin-9.2-amd64.a -L/usr/X11R6/lib64 -lX11  
inpred.o:inpred.f90  $(modules)
	$(f90) $(copt) -c  inpred.f90
mytecio.o: mytecio.f90 
	$(f90) $(copt) -c  mytecio.f90
dislin.o: dislin.f90 
	$(f90) $(copt) -c  dislin.f90
subcode0.o: subcode0.f90
	$(f90) $(copt) -c subcode0.f90
subcode1.o: subcode1.f90 $(modules)
	$(f90) $(copt) -c subcode1.f90
subcode2.o: subcode2.f90 $(modules)
	$(f90) $(copt) -c subcode2.f90
subcode3.o: subcode3.f90 $(modules)
	$(f90) $(copt) -c subcode3.f90
subcode4.o: subcode4.f90 $(modules)	
	$(f90) $(copt) -c subcode4.f90
subcode5.o: subcode5.f90 $(moudles)	
	$(f90) $(copt) -c subcode5.f90
subcode6.o: subcode6.f90 $(modules)	
	$(f90) $(copt) -c subcode6.f90
subcode7.o: subcode7.f90 $(modules)	
	$(f90) $(copt) -c subcode7.f90
subcode8.o: subcode8.f90 $(modules)	
	$(f90) $(copt) -c subcode8.f90
supplement01.o: supplement01.f90 $(modules)	
	$(f90) $(copt) -c supplement01.f90
supplement02.o: supplement02.f90 $(modules)	
	$(f90) $(copt) -c supplement02.f90
clean:
	rm -f $(objects)
	rm -f $(modfile)
	rm -f penmshxp
clear:
	rm -f $(objects)
	rm -f $(modfile)
