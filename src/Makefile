# Makefile

CC = clang

#CFLAGS = -Wall -Wextra -g -ggdb -O0 -pthread -std=gnu11
CFLAGS = -Wall -Wextra -O3 -pthread -std=gnu11

LIBRARY = -lcrypto -lrt -lm
#LIBRARY = -lcrypto -lrt -lm -ljemalloc

MODULES = table coding mempool debug bloom db rwlock stat conc cmap generator

LIBOBJS = $(patsubst %, %.o, $(MODULES))

HEADERS = $(patsubst %, %.h, $(MODULES))

BINARYS = table_test bloom_test rwlock_test generator_test mixed_test cmap_test cm_util io_util staged_read seqio_util

table_test: table_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

bloom_test: bloom_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

rwlock_test: rwlock_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

generator_test: generator_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

mixed_test: mixed_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

cmap_test: cmap_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

cm_util: cm_util.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

io_util: io_util.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

staged_read: staged_read.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

seqio_util: seqio_util.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LIBRARY)

.PHONY : ess all util clean check
ess : table_test mixed_test
util : io_util cm_util seqio_util

all : $(BINARYS)

%.o: %.c
	$(CC) -c $(CFLAGS) -c $< -o $@

clean :
	rm -rf $(BINARYS) *.o

check :
	cppcheck -I /usr/include -DDUMMY --enable=all .
