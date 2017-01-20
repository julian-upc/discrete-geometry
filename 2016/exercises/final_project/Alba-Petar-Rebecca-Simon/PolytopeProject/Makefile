PROJECT=main
SOURCES=Point.cpp Matrix.cpp Gale_Diagram.cpp intersectionsTest.cpp main.cpp
LIBRARY=nope
INCPATHS=. /opt/igraph/include/igraph
LIBPATHS=/opt/igraph/lib /usr/local/lib
LDFLAGS=-lgmp -lgmpxx -ligraph -lstdc++
CFLAGS=-c -Wall -Wextra
CC=gcc

# ------------ MAGIC BEGINS HERE -------------

# Automatic generation of some important lists
OBJECTS=$(SOURCES:.cpp=.o)
INCFLAGS=$(foreach TMP,$(INCPATHS),-I$(TMP))
LIBFLAGS=$(foreach TMP,$(LIBPATHS),-L$(TMP))

# Set up the output file names for the different output types
ifeq "$(LIBRARY)" "shared"
    BINARY=lib$(PROJECT).so
    LDFLAGS += -shared
else ifeq "$(LIBRARY)" "static"
    BINARY=lib$(PROJECT).a
else
    BINARY=$(PROJECT)
endif

all: $(SOURCES) $(BINARY)

$(BINARY): $(OBJECTS)
    # Link the object files, or archive into a static library
    ifeq "$(LIBRARY)" "static"
	ar rcs $(BINARY) $(OBJECTS)
    else
	$(CC) $(LIBFLAGS) $(OBJECTS) $(LDFLAGS) -o $@
    endif

.cpp.o:
	$(CC) $(INCFLAGS) $(CFLAGS) -fPIC $< -o $@

distclean: clean
	rm -f $(BINARY)

clean:
	rm -f $(OBJECTS)
