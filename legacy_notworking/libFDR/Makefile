

all: mincmagnitude inverse3dfft forward3dfft complexmultiply antialiasdata fftshift full3difft full3dfft full2dfft series1dfft full2difft fullmagnitude buildrollofffilter buildbwfilter buildpowerspectrum buildwienerfilter buildinversefdr normalizepsf psfstackfft limitfdr complexrealmultiply rfftshift comparecomplexdata
debug: dofdrfilter-debug mincmagnitude-debug ffttest
clean: 
	rm -f *.so *.o mincmagnitude inverse3dfft forward3dfft complexmultiply antialiasdata fftshift full3difft full3dfft full2dfft series1dfft full2difft fullmagnitude buildrollofffilter buildbwfilter buildpowerspectrum buildwienerfilter buildinversefdr normalizepsf psfstackfft limitfdr complexrealmultiply rfftshift comparecomplexdata


CC = gcc
CFLAGS=-m64 -fPIC -DHAVE_CONFIG_H -I/projects/mice/share/arch/linux64/include -O2
LDFLAGS= -m64 -fPIC -O2
LIBS = -L . -L /projects/mice/share/arch/linux64/lib -lfftw3f -lfftw3 -lm -lminc2 -lhdf5 -lz -lnetcdf


%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

ffttest: ffttest.c
	$(CC) $(LDFLAGS) -o $@ $< $(LIBS)

dofdrfilter: dofdrfilter.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

mincmagnitude: mincmagnitude.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

inverse3dfft: inverse3dfft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

forward3dfft: forward3dfft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

complexmultiply: complexmultiply.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

complexrealmultiply: complexrealmultiply.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

buildcomparefdr: buildcomparefdr.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

buildinversefdr: buildinversefdr.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

buildnoisedampenfilter: buildnoisedampenfilter.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

gaussianmultiply: gaussianmultiply.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

buildwienerfilter: buildwienerfilter.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

buildrollofffilter: buildrollofffilter.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

test3dfft: test3dfft.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

dumpdata: dumpdata.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

antialiasdata: antialiasdata.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

fftshift: fftshift.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

rfftshift: rfftshift.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

full3dfft: full3dfft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

full3difft: full3difft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

full2dfft: full2dfft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

series1dfft: series1dfft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

full2difft: full2difft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

fullmagnitude: fullmagnitude.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

psfstackfft: psfstackfft.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

normalizepsf: normalizepsf.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

limitfdr: limitfdr.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

comparecomplexdata: comparecomplexdata.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

comparerealdata: comparerealdata.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

printrealvalue: printrealvalue.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

printcomplexvalue: printcomplexvalue.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

buildbwfilter: buildbwfilter.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

buildpowerspectrum: buildpowerspectrum.o libfdr2.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr2.so $(LIBS)

mincphase: mincphase.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

mincsubvol: mincsubvol.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

resample_for_comparison: resample_for_comparison.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

do_complete_fdr: do_complete_fdr.o libfdr.so
	$(CC) $(LDFLAGS) -o $@ $< libfdr.so $(LIBS)

libfdr.so: libfdr.o 
	$(CC) $(LDFLAGS) -shared -o $@ $< $(LIBS)

libfdr2.so: libfdr2.o
	$(CC) $(LDFLAGS) -shared -o $@ $< $(LIBS)


