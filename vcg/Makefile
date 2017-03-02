S_OBJECTS=GraphBasic.cxx GraphNode.cxx GraphEdge.cxx PrintGraph.cxx
T_OBJECTS=TestGraph.cxx $(S_OBJECTS)

CC=g++
CFLAGS=-c -fPIC
LDFLAGS=
SOURCES=$(T_OBJECTS)
OBJECTS=$(SOURCES:.cxx=.o)
LIB_OBJECTS=$(S_OBJECTS:.cxx=.o)
EXECUTABLE=TestGraph
LIBRARY=libPrintGraph.so
LINK=ld

all: $(SOURCES) $(EXECUTABLE) $(LIBRARY)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(LIBRARY): $(LIB_OBJECTS)
	$(LINK) $(LDFLAGS) -G $(LIB_OBJECTS) -o $(LIBRARY)

%.o: %.cxx
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf $(OBJECTS)
