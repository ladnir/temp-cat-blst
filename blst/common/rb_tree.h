/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef CATBLST_COMMON_RB_TREE_H
#define CATBLST_COMMON_RB_TREE_H

#include <stddef.h>

namespace blst
{

enum BLST_UNIQ_RESULT {
    BLST_UNIQ_INSERTED = 0,
    BLST_UNIQ_DUPLICATE = 1,
    BLST_UNIQ_CAPACITY_EXCEEDED = 2,
    BLST_UNIQ_DEPTH_EXCEEDED = 3,
    BLST_UNIQ_INVALID_INPUT = 4
};

size_t blst_uniq_sizeof(size_t capacity);
size_t blst_uniq_alignof(void);
bool blst_uniq_init(void *storage, size_t capacity);
BLST_UNIQ_RESULT blst_uniq_insert(void *storage, const void *data, size_t len);

} // namespace blst

#endif
