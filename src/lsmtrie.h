/*
 * Copyright (c) 2014  Wu, Xingbo <wuxb45@gmail.com>
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#pragma once

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/time.h>

struct BloomFilter {
    uint32_t bytes;  // bytes = bits >> 3 (length of filter)
    uint32_t nr_keys;
    uint8_t filter[];
};

// compact bloom_table
// format: encoded bits, bits
#define BLOOMTABLE_INTERVAL ((16u))
struct BloomTable {
    uint8_t *raw_bf;
    uint32_t nr_bf;
    uint32_t nr_bytes;  // size of raw_bf
    uint32_t offsets[];
};

// Container: storing boxes for multiple tables
// Box: multiple related bloom-filter in one box
struct BloomContainer {
    int raw_fd;
    uint64_t off_raw;        // === off_main in MetaFileHeader
    uint32_t nr_barrels;     // === nr_main
    uint32_t nr_bf_per_box;  // 1 -- 8
    uint32_t nr_index;
    uint64_t mtid;
    uint16_t index_last[];  // the LAST barrel_id in each box
};

struct ContainerMap {
    uint64_t nr_units;
    uint64_t nr_used;
    uint64_t total_cap;
    bool discard;
    int raw_fd;
    pthread_mutex_t mutex_cm;  // lock on operating on ContainerMap
    uint8_t bits[];
};

enum GeneratorType {
    GEN_CONSTANT = 0,  // always a constant number
    GEN_COUNTER,       // +1 each fetch
    GEN_DISCRETE,  // gen within a set of values with its weight as probability
    GEN_EXPONENTIAL,  // exponential
    GEN_FILE,         // string from lines in file
    GEN_HISTOGRAM,    // histogram
    GEN_HOTSET,       // hot/cold set
    GEN_ZIPFIAN,      // Zipfian, 0 is the most popular.
    GEN_XZIPFIAN,  // ScrambledZipfian. scatters the "popular" items across the
                   // itemspace.
    GEN_LATEST,  // Generate a popularity distribution of items, skewed to favor
                 // recent items.
    GEN_UNIFORM,  // Uniformly distributed in an interval [a,b]
};

struct GenInfo_Constant {
    uint64_t constant;
};

struct GenInfo_Counter {
    uint64_t counter;
};

struct Pair64 {
    uint64_t a;
    uint64_t b;
};

struct GenInfo_Discrete {
    uint64_t nr_values;
    struct Pair64 *pairs;
};

struct GenInfo_Exponential {
    double gamma;
};

struct GenInfo_File {
    FILE *fin;
};

struct GenInfo_Histogram {
    uint64_t block_size;
    uint64_t area;
    uint64_t weighted_area;
    double main_size;
    uint64_t *buckets;
};

struct GenInfo_HotSet {
    uint64_t lower_bound;
    uint64_t upper_bound;
    uint64_t hot_interval;
    uint64_t cold_interval;
    double hot_set_fraction;
    double hot_op_fraction;
};

#define ZIPFIAN_CONSTANT ((0.99))
struct GenInfo_Zipfian {
    uint64_t nr_items;
    uint64_t base;
    double zipfian_constant;
    double theta;
    double zeta2theta;
    double alpha;
    double zetan;
    double eta;
    uint64_t countforzeta;
    uint64_t min;
    uint64_t max;
};

struct GenInfo_Latest {
    struct GenInfo_Zipfian zipfion;
    uint64_t max;
};

struct GenInfo_Uniform {
    uint64_t min;
    uint64_t max;
    double interval;
};

struct GenInfo {
    uint64_t (*next)(struct GenInfo *const);
    enum GeneratorType type;
    union {
        struct GenInfo_Constant constant;
        struct GenInfo_Counter counter;
        struct GenInfo_Discrete discrete;
        struct GenInfo_Exponential exponential;
        struct GenInfo_File file;
        struct GenInfo_Histogram histogram;
        struct GenInfo_HotSet hotset;
        struct GenInfo_Zipfian zipfian;
        struct GenInfo_Latest latest;
        struct GenInfo_Uniform uniform;
    } gen;
};

struct ReaderLock {
    uint64_t nr_readers;
    bool open;
    pthread_cond_t cond_reader;
    uint64_t pad1[8];
};

struct RWLock {
    pthread_mutex_t mutex_any;
    uint64_t next_ticket;  // sell to writers
    uint64_t reader_ticket;
    uint64_t writer_ticket;
    uint64_t pad1[8];
    pthread_cond_t cond_writer;
    uint64_t pad2[8];
    struct ReaderLock rl[2];  // mod 2
};

struct Stat {
    uint64_t nr_get;
    uint64_t nr_get_miss;
    uint64_t nr_get_at_hit[2];
    uint64_t nr_get_vc_hit[64];

    uint64_t nr_fetch_barrel;
    uint64_t nr_fetch_bc;

    uint64_t nr_true_negative;
    uint64_t nr_false_positive;
    uint64_t nr_true_positive;

    uint64_t nr_set;
    uint64_t nr_set_retry;

    uint64_t nr_compaction;
    uint64_t nr_active_dumped;

    uint64_t nr_write[64];
    uint64_t nr_write_bc;
};

struct KeyValue {
    uint16_t klen;
    uint16_t vlen;
    uint8_t *pk;
    uint8_t *pv;
    uint8_t kv[];  // don't access it
};

#define HASHBYTES ((20))

#define TABLE_MAX_BARRELS ((UINT64_C(8192)))
// a Prime number
#define TABLE_NR_BARRELS ((UINT64_C(8191)))
// 4KB
#define BARREL_ALIGN ((UINT64_C(4096)))
// 32MB
#define TABLE_ALIGN ((BARREL_ALIGN * TABLE_MAX_BARRELS))
// 8MB
#define TABLE_NR_IO ((UINT64_C(2048)))

#define TABLE_ILOCKS_NR ((UINT64_C(64)))

struct Table {
    uint64_t volume;
    uint64_t capacity;
    struct Mempool *mempool;  // store items
    struct Barrel *barrels;
    uint8_t *io_buffer;
    uint64_t nr_mi;
    struct MetaIndex *mis;
    struct BloomTable *bt;
    pthread_mutex_t
        ilocks[TABLE_ILOCKS_NR];  // used for parallel compaction feed
};

struct MetaFileHeader {
    uint64_t off;
    uint64_t volume;
    uint64_t nr_mi;
} __attribute__((packed));

struct MetaTable {
    struct MetaFileHeader mfh;
    int raw_fd;
    uint64_t mtid;
    struct MetaIndex *mis;
    struct BloomTable *bt;
    struct Stat *stat;
};

// bloom.c
struct BloomFilter *bloom_create(uint32_t nr_keys,
                                 struct Mempool *const mempool);
void bloom_update(struct BloomFilter *const bf, uint64_t hv);
bool bloom_match(const struct BloomFilter *const bf, uint64_t hv);
struct BloomTable *bloomtable_build(struct BloomFilter *const *const bfs,
                                    uint64_t nr_bf);
bool bloomtable_dump(struct BloomTable *const bt, FILE *fo);
struct BloomTable *bloomtable_load(FILE *const fi);
bool bloomtable_match(struct BloomTable *const bt, uint32_t index,
                      uint64_t hv);
void bloomtable_free(struct BloomTable *const bt);
struct BloomContainer *bloomcontainer_build(struct BloomTable *const bt,
                                            int raw_fd,
                                            uint64_t off_raw,
                                            struct Stat *const stat);
struct BloomContainer *bloomcontainer_update(struct BloomContainer *const bc,
                                             struct BloomTable *const bt,
                                             int new_raw_fd,
                                             uint64_t new_off_raw,
                                             struct Stat *const stat);
bool bloomcontainer_dump_meta(struct BloomContainer *const bc, FILE *const fo);
struct BloomContainer *bloomcontainer_load_meta(FILE *const fi,
                                                int raw_fd);
uint64_t bloomcontainer_match(struct BloomContainer *const bc,
                              uint32_t index, uint64_t hv);
void bloomcontainer_free(struct BloomContainer *const bc);

// cmap.c
struct ContainerMap *containermap_create(const char *const raw_fn,
                                         uint64_t cap_hint);
struct ContainerMap *containermap_load(const char *const meta_fn,
                                       const char *const raw_fn);
void containermap_dump(struct ContainerMap *const cm,
                       const char *const meta_fn);
void containermap_show(struct ContainerMap *const cm);
uint64_t containermap_alloc(struct ContainerMap *const cm);
bool containermap_release(struct ContainerMap *const cm, uint64_t offset);
void containermap_destroy(struct ContainerMap *const cm);
uint64_t containermap_unused(const struct ContainerMap *const cm);

// coding.c
uint8_t *encode_uint64(uint8_t *const dst, uint64_t v);
const uint8_t *decode_uint64(const uint8_t *const src, uint64_t *const value);

inline uint8_t *encode_uint16(uint8_t *const dst, uint16_t v) {
    return encode_uint64(dst, (uint64_t)v);
}

inline uint8_t *encode_uint32(uint8_t *const dst, uint32_t v) {
    return encode_uint64(dst, (uint64_t)v);
}

inline const uint8_t *decode_uint16(const uint8_t *const src,
                                    uint16_t *const value) {
    uint64_t v = 0;
    const uint8_t *const p = decode_uint64(src, &v);
    // assert(v < UINT16_MAX)
    *value = (typeof(*value))v;
    return p;
}

inline const uint8_t *decode_uint32(const uint8_t *const src,
                                    uint32_t *const value) {
    uint64_t v = 0;
    const uint8_t *const p = decode_uint64(src, &v);
    *value = (typeof(*value))v;
    return p;
}

// conc.c
void conc_set_affinity_0(void);
void conc_set_affinity_n(uint64_t cpu);
void conc_fork_reduce(uint64_t nr, void *(*func)(void *),
                      void *const arg);

// db.c
struct DB;
struct DB *db_touch(const char *const meta_dir, const char *const cm_conf_fn);
void db_close(struct DB *const db);
bool db_insert(struct DB *const db, struct KeyValue *const kv);
bool db_multi_insert(struct DB *const db, uint64_t nr_items,
                     const struct KeyValue *const kvs);
struct KeyValue *db_lookup(struct DB *const db, uint16_t klen,
                           const uint8_t *const key);
void db_force_dump_meta(struct DB *const db);
void db_stat_show(struct DB *const db, FILE *const fo);
void db_stat_clean(struct DB *const db);
bool db_doing_compaction(struct DB *const db);

// debug.c
uint64_t debug_time_usec(void);
double debug_time_sec(void);
uint64_t debug_diff_usec(uint64_t last);
double debug_diff_sec(double last);
uint64_t debug_tv_diff(const struct timeval *const t0,
                       const struct timeval *const t1);
void debug_print_tv_diff(char *tag, const struct timeval t0,
                         const struct timeval t1);
void debug_trace(void);

// generator.c

uint64_t random_uint64(void);
struct GenInfo *generator_new_constant(uint64_t constant);
struct GenInfo *generator_new_counter(uint64_t start);
struct GenInfo *generator_new_exponential(double percentile,
                                          double range);
struct GenInfo *generator_new_zipfian(uint64_t min, uint64_t max);
struct GenInfo *generator_new_xzipfian(uint64_t min, uint64_t max);
struct GenInfo *generator_new_uniform(uint64_t min, uint64_t max);
void generator_destroy(struct GenInfo *const geninfo);

// mempool.c
struct Mempool;
void *huge_alloc(uint64_t cap);
void huge_free(void *const ptr, uint64_t cap);
struct Mempool *mempool_new(size_t cap);
uint8_t *mempool_alloc(struct Mempool *const p, size_t cap);
void mempool_free(struct Mempool *const p);
void mempool_show(struct Mempool *const p);

// rwlock.c
void rwlock_show(struct RWLock *bo);
void rwlock_initial(struct RWLock *bo);
uint64_t rwlock_reader_lock(struct RWLock *bo);
void rwlock_reader_unlock(struct RWLock *bo, uint64_t ticket);
uint64_t rwlock_writer_lock(struct RWLock *bo);
void rwlock_writer_unlock(struct RWLock *bo, uint64_t ticket);

// stat.c
void stat_inc(uint64_t *const p);
void stat_inc_n(uint64_t *const p, uint64_t n);
void stat_show(struct Stat *const stat, FILE *const out);
uint32_t *latency_initial(void);
void latency_record(uint64_t usec, uint32_t *const counters);
void latency_show(const char *const tag, uint32_t *const counters,
                  FILE *const out);
void latency_95_99_999(uint32_t *const counters, FILE *const out);

// table.c
uint16_t table_select_barrel(const uint8_t *const hash);
bool table_retain(struct Table *const table);
struct Table *table_alloc_new(double cap_percent,
                              double mempool_factor);
struct Table *table_alloc_default(double mempool_factor);
bool table_insert_kv_safe(struct Table *const table,
                          const struct KeyValue *const kv);
bool table_full(const struct Table *const table);
struct KeyValue *table_lookup(struct Table *const table, uint16_t klen,
                              const uint8_t *const key,
                              const uint8_t *const hash);
bool table_build_bloomtable(struct Table *const table);
bool table_dump_meta(struct Table *const table, const char *const metafn,
                     uint64_t off);
uint64_t table_dump_barrels(struct Table *const table, int fd,
                            uint64_t off);
void table_free(struct Table *const table);
void table_analysis_verbose(struct Table *const table, FILE *const out);
void table_analysis_short(struct Table *const table, char *const buffer);
void table_show(struct Table *const table, FILE *const fo);
struct MetaTable *metatable_load(const char *const metafn, int raw_fd,
                                 bool load_bf, struct Stat *const stat);
struct KeyValue *metatable_lookup(struct MetaTable *const mt,
                                  uint16_t klen, const uint8_t *const key,
                                  const uint8_t *const hash);
void metatable_free(struct MetaTable *const mt);
bool metatable_feed_barrels_to_tables(
    struct MetaTable *const mt, uint16_t start, uint16_t nr,
    uint8_t *const arena, struct Table *const *const tables,
    uint64_t (*select_table)(const uint8_t *const, uint64_t),
    uint64_t arg2);
