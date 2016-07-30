# Makefile

CC = clang

#CFLAGS = -Wall -Wextra -g -ggdb -O0 -pthread -std=gnu11
CFLAGS = -Wall -Wextra -O3 -pthread -std=gnu11

LDFLAGS = -lcrypto -lrt -lm
#LDFLAGS = -lcrypto -lrt -lm -ljemalloc

MODULES = table coding mempool debug bloom db rwlock stat conc cmap generator

LIBOBJS = $(patsubst %, %.o, $(MODULES))

HEADERS = $(patsubst %, %.h, $(MODULES))

BINARYS = table_test bloom_test rwlock_test generator_test mixed_test cmap_test cm_util io_util staged_read seqio_util

table_test: table_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

bloom_test: bloom_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

rwlock_test: rwlock_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

generator_test: generator_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

mixed_test: mixed_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

cmap_test: cmap_test.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

cm_util: cm_util.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

io_util: io_util.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

staged_read: staged_read.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

seqio_util: seqio_util.c $(LIBOBJS)
	$(CC) $(CFLAGS) -o $@ $< $(LIBOBJS) $(LDFLAGS)

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
