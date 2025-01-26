EXEC = D32

CLIB = -I./lib/portaudio/include ./lib/portaudio/lib/.libs/libportaudio.a -lrt -lasound -ljack -pthread 

MAIN = main4.cpp

$(EXEC): $(MAIN)
	g++ -o $@ $^ $(CLIB)

clean:
	rm -f $(EXEC)
.PHONY: clean

run:
	./$(EXEC)
.PHONY: run

cm:
	rm -f $(EXEC)
	g++ -o $(EXEC) $(MAIN) $(CLIB)
.PHONY: cm

mr:
	g++ -o $(EXEC) $(MAIN) $(CLIB)
	./$(EXEC)
.PHONY: cmr

cmr:
	rm -f $(EXEC)
	g++ -o $(EXEC) $(MAIN) $(CLIB)
	./$(EXEC)
.PHONY: cmr