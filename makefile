EXEC = D32

# CLIB = -I./lib/portaudio/include ./lib/portaudio/lib/.libs/libportaudio.a -lrt -lasound -ljack -pthread 
CLIB =  ./lib/miniaudio/miniaudio.c -I./lib/miniaudio

MAIN = main.cpp

$(EXEC): $(MAIN)
	g++ -std=c++11 -o $@ $^ $(CLIB)

clean:
	rm -f $(EXEC)
.PHONY: clean

run:
	./$(EXEC)
.PHONY: run

cm:
	make clean
	make
.PHONY: cm

mr:
	make
	make run
.PHONY: cmr

cmr:
	make clean
	make
	make run
.PHONY: cmr