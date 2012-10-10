QEMINSTALLDIR =	$(shell	cat ../QEMINSTALLDIR)

ALL: mex

mex: 
	mex -fortran hwm93.F hwm93_main.f

clean:
	-rm *.mex* *~ *.o

install:
	-cp *.mex* $(QEMINSTALLDIR)
	-cp *.m $(QEMINSTALLDIR)

dist: clean
	-rm hwm93???.f
