/*
 * Copyright (c) 2014  Wu, Xingbo <wuxb45@gmail.com>
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */

#define _GNU_SOURCE
#define _LARGEFILE64_SOURCE

#include <inttypes.h>
#include "lsmtrie.h"

uint8_t* encode_uint64(uint8_t* dst, uint64_t v) {
    static const uint64_t B = 0x80;
    static const uint64_t M = 0x7f;
    uint8_t* ptr = dst;
    uint64_t t = v;
    while (t >= B) {
        *ptr = (uint8_t)((t & M) | B);
        ++ptr;
        t >>= 7;
    }
    *ptr = (uint8_t)t;
    ++ptr;
    return ptr;
}

const uint8_t* decode_uint64(const uint8_t* src, uint64_t* value) {
    uint64_t result = 0;
    static const uint64_t B = 0x80;
    static const uint64_t M = 0x7f;
    const uint8_t* p = src;

    for (uint32_t shift = 0; shift <= 63; shift += 7) {
        const uint64_t byte = (uint64_t)(*p);
        ++p;
        if (byte & B) {
            // More bytes are present
            result |= ((byte & M) << shift);
        } else {
            result |= (byte << shift);
            *value = result;
            return p;
        }
    }
    *value = 0;
    return src;
}
