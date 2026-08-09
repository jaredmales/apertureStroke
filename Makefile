APPS = apertureStroke zernikeTheory

.PHONY: all clean

all: $(APPS)

apertureStroke: apertureStroke.cpp basisFitters.hpp strokeUtils.hpp
	$(MAKE) -f $(MXMAKEFILE) t=$@

zernikeTheory: zernikeTheory.cpp strokeUtils.hpp
	$(MAKE) -f $(MXMAKEFILE) t=$@

clean:
	$(RM) $(APPS) *.o
