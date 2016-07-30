# Makefile

CC = clang

#CFLAGS = -Wall -Wextra -g -ggdb -O0 -pthread -std=gnu11
CFLAGS = -Wall -Wextra -O3 -pthread -std=gnu11

LIBRARY = -lcrypto -lrt -lm
#LIBRARY = -lcrypto -lrt -lm -ljemalloc

MODULES = table coding mempool debug bloom db rwlock stat conc cmap generator

SOURCES = $(patsubst %, %.c, $(MODULES))

HEADERS = $(patsubst %, %.h, $(MODULES))

DEPS = $(SOURCES) $(HEADERS)

BINARYS = table_test bloom_test rwlock_test generator_test mixed_test cmap_test cm_util io_util staged_read seqio_util

table_test: table_test.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

bloom_test: bloom_test.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

rwlock_test: rwlock_test.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

generator_test: generator_test.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

mixed_test: mixed_test.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

cmap_test: cmap_test.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

cm_util: cm_util.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

io_util: io_util.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

staged_read: staged_read.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

seqio_util: seqio_util.c $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

.PHONY : ess all util clean check
ess : table_test mixed_test
util : io_util cm_util seqio_util

all : $(BINARYS)

% : %.c $(DEPS)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) $(LIBRARY)

clean :
	rm -rf $(BINARYS) *.o

check :
	cppcheck -I /usr/include -DDUMMY --enable=all .
