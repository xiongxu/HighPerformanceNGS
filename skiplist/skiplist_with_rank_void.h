/*
 * Copyright (C) 2015, Leo Ma <begeekmyfriend@gmail.com>
 */

#ifndef _SKIPLIST_H
#define _SKIPLIST_H
#include <stdlib.h>

struct sk_link {
        struct sk_link *next, *prev;
        unsigned int span;
};

static inline void list_init(struct sk_link *link) {
        link->prev = link;
        link->next = link;
}

static inline void __list_add(struct sk_link *link, struct sk_link *prev, struct sk_link *next) {
        link->next = next;
        link->prev = prev;
        next->prev = link;
        prev->next = link;
}

static inline void __list_del(struct sk_link *prev, struct sk_link *next) {
        prev->next = next;
        next->prev = prev; 
}

static inline void list_add(struct sk_link *link, struct sk_link *next) {
        __list_add(link, next->prev, next);
}

static inline void list_del(struct sk_link *link) {
        __list_del(link->prev, link->next);
        list_init(link);
}

static inline int list_empty(struct sk_link *link) {
        return link->next == link;
}

#define list_entry(ptr, type, member) \
        ((type *)((char *)(ptr) - (size_t)(&((type *)0)->member)))

#define skiplist_foreach_forward(pos, end) \
        for (; pos != end; pos = pos->next)

#define skiplist_foreach_forward_safe(pos, n, end) \
        for (n = pos->next; pos != end; pos = n, n = pos->next)

#define skiplist_foreach_backward(pos, end) \
        for (; pos != end; pos = pos->prev)

#define skiplist_foreach_backward_safe(pos, n, end) \
        for (n = pos->prev; pos != end; pos = n, n = pos->prev)

#define MAX_LEVEL 32  /* Should be enough for 2^32 elements */

struct skiplist {
        int level;
        int count;
        struct sk_link head[MAX_LEVEL];
};

// void *(*keyDup)(const void *key);
// void *(*valDup)(const void *obj);
// int (*keyCompare)(const void *key1, const void *key2);
// void (*keyDestructor)(void *key);
// void (*valDestructor)(void *obj);

#define sklFreeEntryKey(keyDestructor, node) \
    if (keyDestructor) \
        keyDestructor((node)->key)

#define sklSetKey(keyDup, node, _key_) do { \
    if (keyDup) \
        node->key = keyDup(_key_); \
    else \
        node->key = (_key_); \
} while(0)

#define sklFreeEntryVal(valDestructor, node) \
    if (valDestructor) \
        valDestructor((node)->value)

#define sklSetVal(valDup, node, _val_) do { \
    if (valDup) \
        node->value = valDup(_val_); \
    else \
        node->value = (_val_); \
} while(0)

#define sklCompareKeys(keyCompare, key1, key2) \
    ((keyCompare) ? \
        keyCompare(key1, key2) : \
        (key1) == (key2))

struct skipnode {
        void *key;
        void *value;
        struct sk_link link[0];
};

static struct skipnode *skipnode_new(int level,void *key,void *value,void *(*keyDup)(const void *key),void *(*valDup)(const void *obj)) {
        struct skipnode *node = malloc(sizeof(*node) + level * sizeof(struct sk_link));
        if (node != NULL) {
                sklSetKey(keyDup,node,key);
                sklSetVal(valDup,node,value);
        }
        return node;
}

static void skipnode_delete(struct skipnode *node,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj)) {
        sklFreeEntryKey(keyDestructor,node);
        sklFreeEntryVal(valDestructor,node);
        free(node);
}

static struct skiplist *skiplist_new(void) {
        int i;
        struct skiplist *list = malloc(sizeof(*list));
        if (list != NULL) {
                list->level = 1;
                list->count = 0;
                for (i = 0; i < sizeof(list->head) / sizeof(list->head[0]); i++) {
                        list_init(&list->head[i]);
                        list->head[i].span = 0;
                }
        }
        return list;
}

static void skiplist_delete(struct skiplist *list,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj)) {
        struct sk_link *n;
        struct sk_link *pos = list->head[0].next;
        skiplist_foreach_forward_safe(pos, n, &list->head[0]) {
                struct skipnode *node = list_entry(pos, struct skipnode, link[0]);
                skipnode_delete(node,keyDestructor,valDestructor);
        }
        free(list);
}

static int random_level(void) {
        const double p = 0.25;
        int level = 1;
        while ((random() & 0xffff) < (long)(0xffff * p)) {
                level++;
        }
        return level > MAX_LEVEL ? MAX_LEVEL : level;
}

static struct skipnode *skiplist_insert(struct skiplist *list, void *key,void *value, void *(*keyDup)(const void *key),void *(*valDup)(const void *obj), int (*keyCompare)(const void *key1, const void *key2)) {
        int rank[MAX_LEVEL];
        struct sk_link *update[MAX_LEVEL];
        int level = random_level();
        if (level > list->level) {
                list->level = level;
        }
        struct skipnode *node = skipnode_new(level, key, value,keyDup,valDup);
        if (node != NULL) {
                int i = list->level - 1;
                struct sk_link *pos = &list->head[i];
                struct sk_link *end = &list->head[i];
                for (; i >= 0; i--) {
                        rank[i] = i == list->level - 1 ? 0 : rank[i + 1];
                        pos = pos->next;
                        skiplist_foreach_forward(pos, end) {
                                struct skipnode *nd = list_entry(pos, struct skipnode, link[i]);
                                if (keyCompare(nd->key,key)>=0) {
                                        end = &nd->link[i];
                                        break;
                                }
                                rank[i] += nd->link[i].span;
                        }
                        update[i] = end;
                        pos = end->prev;
                        pos--;
                        end--;
                }
                for (i = 0; i < list->level; i++) {
                        if (i < level) {
                                list_add(&node->link[i], update[i]);
                                node->link[i].span = rank[0] - rank[i] + 1;
                                update[i]->span -= node->link[i].span - 1;
                        } else {
                                update[i]->span++;
                        }
                }
                list->count++;
        }
        return node;
}

static void __remove(struct skiplist *list, struct skipnode *node, int level, struct sk_link **update,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj) ) {
        int i;
        int remain_level = list->level;
        for (i = 0; i < list->level; i++) {
                if (i < level) {
                        list_del(&node->link[i]);
                        update[i] = node->link[i].next;
                        update[i]->span += node->link[i].span - 1;
                } else {
                        update[i]->span--;
                }
                if (list_empty(&list->head[i])) {
                        if (remain_level == list->level) {
                                remain_level = i + 1;
                        }
                }
        }
        skipnode_delete(node,keyDestructor,valDestructor);
        list->count--;
        list->level = remain_level;
}

static void skiplist_remove(struct skiplist *list,const void *key,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj), int (*keyCompare)(const void *key1, const void *key2) ) {
        int i;
        struct skipnode *node;
        struct sk_link *pos, *n, *end, *update[MAX_LEVEL];
        i = list->level - 1;
        pos = &list->head[i];
        end = &list->head[i];
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward_safe(pos, n, end) {
                        node = list_entry(pos, struct skipnode, link[i]);
                        if (keyCompare(node->key,key) >0) {
                                end = &node->link[i];
                                break;
                        } else if (keyCompare(node->key,key)==0) {
                                __remove(list, node, i + 1, update,keyDestructor,valDestructor);
                                return;
                        }
                }
                update[i] = end;
                pos = end->prev;
                pos--;
                end--;
        }
}

struct range_spec {
        void *min, *max;
        int minex, maxex;
};

static int key_gte_min(const void *key, struct range_spec *range, int (*keyCompare)(const void *key1, const void *key2) ) {
        return range->maxex ? keyCompare(key,range->max)>0 : keyCompare(key,range->max)>=0;
}

static int key_lte_max(const void *key, struct range_spec *range,  int (*keyCompare)(const void *key1, const void *key2)) {
        return range->minex ? keyCompare(key,range->min)<0 : keyCompare(key,range->min) <= 0;
}

/* Returns if there is node key in range */
static int key_in_range(struct skiplist *list, struct range_spec *range,int (*keyCompare)(const void *key1, const void *key2) ) {
        struct skipnode *node;
        struct sk_link *link;
        if (keyCompare(range->min,range->max)>0 || (keyCompare(range->min,range->max)==0 && (range->minex || range->maxex) ) ) 
                return 0;
        if (list_empty(&list->head[0]))
                return 0;
        link = list->head[0].next;
        node = list_entry(link, struct skipnode, link[0]);
        if (!key_lte_max(node->key, range,keyCompare))
                return 0;
        link = list->head[0].prev;
        node = list_entry(link, struct skipnode, link[0]);
        if (!key_gte_min(node->key, range,keyCompare))
                return 0;
        return 1;
}

/* search the first node key that is contained in the specified range
 * where min and max are inclusive. */
static struct skipnode *first_in_range(struct skiplist *list, struct range_spec *range,  int (*keyCompare)(const void *key1, const void *key2)) {
        if (!key_in_range(list, range,keyCompare))
                return NULL;
        int i = list->level - 1;
        struct sk_link *pos = &list->head[i];
        struct sk_link *end = &list->head[i];
        struct skipnode *node = NULL;
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward(pos, end) {
                        node = list_entry(pos, struct skipnode, link[i]);
                        if (key_gte_min(node->key, range,keyCompare)) {
                                pos = node->link[i].prev;
                                end = node->link[i].next;
                                goto CONTINUE;
                        }
                }
                pos = end->prev;
CONTINUE:
                pos--;
                end--;
        }
        return node;
}

/* search the last node key that is contained in the specified range
 * where min and max are inclusive. */
static struct skipnode *last_in_range(struct skiplist *list, struct range_spec *range,  int (*keyCompare)(const void *key1, const void *key2)) {
        if (!key_in_range(list, range,keyCompare))
                return NULL;
        int i = list->level - 1;
        struct sk_link *pos = &list->head[i];
        struct sk_link *end = &list->head[i];
        struct skipnode *node = NULL;
        for (; i >= 0; i--) {
                pos = pos->prev;
                skiplist_foreach_backward(pos, end) {
                        node = list_entry(pos, struct skipnode, link[i]);
                        if (key_lte_max(node->key, range,keyCompare)) {
                                pos = node->link[i].next;
                                end = node->link[i].prev;
                                goto CONTINUE;
                        }
                }
                pos = end->next;
CONTINUE:
                pos--;
                end--;
        }
        return node;
}

/* remove all the nodes with key in range
 * where min and max are inclusive. */
static unsigned int remove_in_range(struct skiplist *list, struct range_spec *range,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj),  int (*keyCompare)(const void *key1, const void *key2) ) {
        int i;
        unsigned int removed = 0;
        struct sk_link *pos, *n, *end, *update[MAX_LEVEL];
        if (key_in_range(list, range,keyCompare))
                return 0;
        i = list->level - 1;
        pos = &list->head[i];
        end = &list->head[i];
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward_safe(pos, n, end) {
                        struct skipnode *node = list_entry(pos, struct skipnode, link[i]);
                        if (!key_lte_max(node->key, range,keyCompare)) {
                                end = &node->link[i];
                                break;
                        } else if (key_gte_min(node->key, range,keyCompare)) {
                                /* No return because we allow nodes with same key. */
                                __remove(list, node, i + 1, update,keyDestructor,valDestructor);
                                removed++;
                        }
                }
                update[i] = end;
                pos = end->prev;
                pos--;
                end--;
        }
        return removed;
}

/* remove all the nodes with key rank in range
 * where start and stop are inclusive. */
static unsigned int remove_in_rank(struct skiplist *list, unsigned int start, unsigned int stop,void (*keyDestructor)(const void *key), void (*valDestructor)(const void *obj)) {
        int i;
        unsigned int removed = 0, traversed = 0;
        struct sk_link *pos, *n, *end, *update[MAX_LEVEL];
        if (start <= 0 || stop <= 0 || start > list->count)
                return 0;
        i = list->level - 1;
        pos = &list->head[i];
        end = &list->head[i];
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward_safe(pos, n, end) {
                        struct skipnode *node = list_entry(pos, struct skipnode, link[i]);
                        if (traversed + node->link[i].span > stop) {
                                end = &node->link[i];
                                break;
                        } else if (traversed + node->link[i].span >= start) {
                                /* No return because we allow nodes with same key. */
                                __remove(list, node, i + 1, update,keyDestructor,valDestructor);
                                removed++;
                                continue;
                        }
                        traversed += node->link[i].span;
                }
                update[i] = end;
                pos = end->prev;
                pos--;
                end--;
        }
        fprintf(stderr,"removed %d nodes\n",removed);
        return removed;
}

/* Get the node key rank */
static unsigned int skiplist_key_rank(struct skiplist *list, const void *key, int (*keyCompare)(const void *key1, const void *key2)) {
        unsigned int rank = 0;
        int i = list->level - 1;
        struct sk_link *pos = &list->head[i];
        struct sk_link *end = &list->head[i];
        struct skipnode *node;
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward(pos, end) {
                        node = list_entry(pos, struct skipnode, link[i]);
                        if (keyCompare(node->key , key)>=0) {
                                end = &node->link[i];
                                break;
                        }
                        rank += node->link[i].span;
                }
                if (keyCompare(node->key,key)==0) {
                        return rank + node->link[i].span;
                }
                pos = end->prev;
                pos--;
                end--;
        }
        return 0;
}

/* search the node with specified key. */
static struct skipnode *skiplist_search_by_key(struct skiplist *list,const void *key, int (*keyCompare)(const void *key1, const void *key2)) {
        int i = list->level - 1;
        struct sk_link *pos = &list->head[i];
        struct sk_link *end = &list->head[i];
        struct skipnode *node;
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward(pos, end) {
                        node = list_entry(pos, struct skipnode, link[i]);
                        if (keyCompare(node->key,key)>=0) {
                                end = &node->link[i];
                                break;
                        }
                }
                if (keyCompare(node->key,key) == 0) {
                        return node;
                }
                pos = end->prev;
                pos--;
                end--;
        }
        return NULL;
}

/* search the node with specified key rank. */
static struct skipnode *skiplist_search_by_rank(struct skiplist *list, unsigned int rank) {
        if (rank == 0 || rank > list->count) {
                return NULL;
        }
        int i = list->level - 1;
        unsigned int traversed = 0;
        struct sk_link *pos = &list->head[i];
        struct sk_link *end = &list->head[i];
        struct skipnode *node;
        for (; i >= 0; i--) {
                pos = pos->next;
                skiplist_foreach_forward(pos, end) {
                        node = list_entry(pos, struct skipnode, link[i]);
                        if (traversed + node->link[i].span >= rank) {
                                end = &node->link[i];
                                break;
                        }
                        traversed += node->link[i].span;
                }
                if (rank == traversed + node->link[i].span) {
                        fprintf(stderr,"%d %d %d %d\n",i,rank,traversed ,node->link[i].span);
                        return node;
                }
                pos = end->prev;
                pos--;
                end--;
        }
        return NULL;
}

static void skiplist_dump(struct skiplist *list) {
        int i = list->level - 1;
        unsigned int traversed = 0;
        struct sk_link *pos = &list->head[i];
        struct sk_link *end = &list->head[i];
        printf("\nTotal %d nodes: \n", list->count);
        for (; i >= 0; i--) {
                traversed = 0;
                pos = pos->next;
                skiplist_foreach_forward(pos, end) {
                        struct skipnode *node = list_entry(pos, struct skipnode, link[i]);
                        traversed += node->link[i].span;
                        printf("level:%d key:0x%08x value:0x%08x rank:%u\n",
                                i + 1, node->key, node->value, traversed);
                }
                pos = &list->head[i];
                pos--;
                end--;
        }
}

#endif  /* _SKIPLIST_H */
