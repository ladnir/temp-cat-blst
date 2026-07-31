/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stddef.h>
#include "blst/common/rb_tree.h"
/*
 * Red-black tree tailored for uniqueness test. Amount of messages to be
 * checked is known prior context initialization, implementation is
 * insert-only, and message storage remains owned by the caller for the
 * lifetime of the tree.
 */

namespace blst
{

struct node {
    struct node *leafs[2];
    const void *data;
    size_t len_n_colour;    /* len<<1 | colour */
};

struct rb_tree {
    struct node *root;
    size_t n_nodes;
    size_t capacity;
    struct node nodes[1];
};

static int bytes_compare(const unsigned char *ptr0, size_t len0,
                         const unsigned char *ptr1, size_t len1)
{
    size_t i, len = len0<len1 ? len0 : len1;

    for (i=0; i<len; i++) {
        if (ptr0[i] != ptr1[i])
            return ptr0[i] < ptr1[i] ? -1 : 1;
    }

    return (len0 > len1) - (len0 < len1);
}

#define PAINT_BLACK(p)  ((p)->len_n_colour &= ~(size_t)1)
#define PAINT_RED(p)    ((p)->len_n_colour |= 1)
#define IS_RED(p)       ((p)->len_n_colour & 1)

static BLST_UNIQ_RESULT rb_tree_insert(struct rb_tree *tree,
                                      const void *data, size_t len)
{
    /*
     * A red-black tree with n nodes has height at most 2*log2(n+1).
     * Sixty-four entries cover every representable 32-bit allocation. On
     * 64-bit targets, exceeding this path would require more than 2^32 nodes
     * in the worst case and is reported explicitly before the array write.
     * Keeping the original 64-bit stack footprint avoids penalizing normal
     * insertion for an impractically large caller-owned tree.
     */
    static const size_t max_depth = 64;
    struct node *nodes[64];              /* visited nodes    */
    unsigned char dirs[64];              /* taken directions */
    size_t k = 0;                           /* walked distance  */
    struct node *p, *y, *z;

    if (tree == NULL || (data == NULL && len != 0) || len > ((size_t)-1 >> 1))
        return BLST_UNIQ_INVALID_INPUT;

    for (p = tree->root; p != NULL; k++) {
        int cmp = bytes_compare((const unsigned char*)data, len,
                                (const unsigned char*)p->data,
                                p->len_n_colour>>1);
        if (cmp == 0)
            return BLST_UNIQ_DUPLICATE;

        /* record the step */
        if (k == max_depth)
            return BLST_UNIQ_DEPTH_EXCEEDED;
        nodes[k] = p;
        p = p->leafs[(dirs[k] = cmp>0)];
    }

    /* allocate new node */
    if (tree->n_nodes >= tree->capacity)
        return BLST_UNIQ_CAPACITY_EXCEEDED;
    z = &tree->nodes[tree->n_nodes++];
    z->leafs[0] = z->leafs[1] = NULL;
    z->data = data;
    z->len_n_colour = len<<1;
    PAINT_RED(z);

    /* graft |z| */
    if (k > 0)
        nodes[k-1]->leafs[dirs[k-1]] = z;
    else
        tree->root = z;

    /* re-balance |tree| */
    while (k >= 2 && IS_RED(y = nodes[k-1])) {
        size_t ydir = dirs[k-2];
        struct node *x = nodes[k-2],        /* |z|'s grandparent    */
                    *s = x->leafs[ydir^1];  /* |z|'s uncle          */

        if (s != NULL && IS_RED(s)) {
            PAINT_RED(x);
            PAINT_BLACK(y);
            PAINT_BLACK(s);
            k -= 2;
        } else {
            if (dirs[k-1] != ydir) {
                /*    |        |
                 *    x        x
                 *   / \        \
                 *  y   s -> z   s
                 *   \      /
                 *    z    y
                 *   /      \
                 *  ?        ?
                 */
                struct node *t = y;
                y = y->leafs[ydir^1];
                t->leafs[ydir^1] = y->leafs[ydir];
                y->leafs[ydir] = t;
            }

            /*      |        |
             *      x        y
             *       \      / \
             *    y   s -> z   x
             *   / \          / \
             *  z   ?        ?   s
             */
            x->leafs[ydir] = y->leafs[ydir^1];
            y->leafs[ydir^1] = x;

            PAINT_RED(x);
            PAINT_BLACK(y);

            if (k > 2)
                nodes[k-3]->leafs[dirs[k-3]] = y;
            else
                tree->root = y;

            break;
        }
    }

    PAINT_BLACK(tree->root);

    return BLST_UNIQ_INSERTED;
}

#undef IS_RED
#undef PAINT_RED
#undef PAINT_BLACK

size_t blst_uniq_sizeof(size_t capacity)
{
    const size_t header = offsetof(struct rb_tree, nodes);

    if (capacity == 0 ||
        capacity > ((size_t)-1 - header) / sizeof(struct node))
        return 0;
    return header + sizeof(struct node)*capacity;
}

size_t blst_uniq_alignof(void)
{
    return alignof(struct rb_tree);
}

bool blst_uniq_init(void *storage, size_t capacity)
{
    if (storage == NULL ||
        reinterpret_cast<size_t>(storage) % alignof(struct rb_tree) != 0 ||
        blst_uniq_sizeof(capacity) == 0)
        return false;

    struct rb_tree *tree = static_cast<struct rb_tree*>(storage);
    tree->root = NULL;
    tree->n_nodes = 0;
    tree->capacity = capacity;
    return true;
}

BLST_UNIQ_RESULT blst_uniq_insert(void *storage, const void *data, size_t len)
{
    if (storage == NULL ||
        reinterpret_cast<size_t>(storage) % alignof(struct rb_tree) != 0)
        return BLST_UNIQ_INVALID_INPUT;
    return rb_tree_insert(static_cast<struct rb_tree*>(storage), data, len);
}

} // namespace blst
