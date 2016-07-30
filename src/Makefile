CC = cc

#OPT = -g -ggdb -O0
OPT = -O3

CFLAGS = -Wall -Wextra $(OPT) -pthread -std=gnu11

ENABLE_JEMALLOC =
#ENABLE_JEMALLOC = -ljemalloc

LDFLAGS = -lcrypto -lrt -lm $(ENABLE_JEMALLOC)

MODULES = table coding mempool debug bloom db rwlock stat conc cmap generator

LIBOBJS = $(patsubst %, %.o, $(MODULES))

HEADERS = $(patsubst %, %.h, $(MODULES))

TESTS = table_test bloom_test rwlock_test generator_test mixed_test cmap_test cm_util io_util staged_read seqio_util

default: all
all: $(TESTS)

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

.PHONY: ess all util clean check
ess: table_test mixed_test
util: io_util cm_util seqio_util

%.o: %.c
	$(CC) -c $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(TESTS) *.o

check:
	cppcheck -I /usr/include -DDUMMY --enable=all .
