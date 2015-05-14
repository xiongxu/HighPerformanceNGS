//  gcc -g -O3 -Wall -fPIC -c  list.c && ar -rcs liblist.a list.o && rm list.o
//  list.c
//  adjacencylist
//
//  Created by Sam Goldman on 6/21/11.
//  Copyright 2011 Sam Goldman. All rights reserved.
//

#include "list.h"
#include <stdlib.h>

Node *node_create(void *data);

List *list_create(node_data_free_callback_t node_data_free_callback,int (*cmp)(const void *a, const void *b)) {
	List *list = malloc(sizeof(List));
	list->head = NULL;
	list->count = 0;
	list->node_data_free_callback = node_data_free_callback;
	list->cmp = cmp;
	return list;
}

Node *node_create(void *data) {
	Node *node = malloc(sizeof(Node));
	node->data = data;
	node->next = NULL;
	return node;
}

void list_add_data(List *list, void *data) {
	Node *node = node_create(data);
	node->next = list->head;
	list->head = node;
	list->count++;
}

void list_add_data_sorted(List *list, void *data) {
	Node *node = node_create(data);
	Node *n = list->head;
	Node *prev_n = NULL;
	while (n && list->cmp(n->data, data) < 0) {
		prev_n = n;
		n = n->next;
	}
	node->next = n;
	if (!prev_n) {
		list->head = node;
	}
	else {
		prev_n->next = node;
	}
	list->count++;
}

void list_remove_data(List *list, void *data) {
	Node *n = list->head;
	Node *prev_n = NULL;
	while (n) {
		if (!list->cmp(n->data, data)) {
			if (!prev_n) {
				list->head = n->next;
			}
			else {
				prev_n->next = n->next;
			}
			list->node_data_free_callback(n->data);
			list->count--;
			free(n);
			break;
		}
		prev_n = n;
		n = n->next;
	}
}

int list_pop_front(List *list){
	if (!list || !list->count) return -1;
	Node *n = list->head;
	list->head = n->next ? n->next : NULL;
	list->node_data_free_callback(n->data);
	list->count--;
	return list->count;
}

void *shift(List *list){
	if (!list || !list->count) return NULL;
	Node *n = list->head;
	list->head = n->next ? n->next : NULL;
	list->count--;
	return n->data;
}

int list_empty(List *list) {
	if (!list) return -1;
	return !list->count;
}

void list_free(List *list) {
	Node *n = list->head;
	while (n) {
		Node *next_n = n->next;
		list->node_data_free_callback(n->data);
		free(n);
		n = next_n;
	}
	free (list);
}
/*
ListNode *mergeSortedList(ListNode *L1, ListNode *L2)
{
	ListNode dummy(-1), *p1 = &dummy, *p2 = L2;  //L1的辅助头结点dummy，因为可能在头部插入
	dummy.next = L1;
	while(p1->next != NULL && p2 != NULL)  //停止条件，也包括了判断两个链表是否为空
	{
		if(p1->next->value >= p2->value)
		{
			L2 = p2->next;
			p2->next = p1->next;
			p1->next = p2;
			p1 = p2;
			p2 = L2;
		}
		else
		{
			p1 = p1->next;
		}
	}
	if(p1->next == NULL)    //L2可能还有未处理的结点，直接加在L1尾部即可
	{
		p1->next = p2;
	}

	return dummy.next;
}

List *listMergeSort(List *head){
	if(head == NULL || head->next == NULL)   //链表为空，或者只有一个结点，直接返回
		return head;
	
	List *slow = head, *fast = head;
	while(fast->next != NULL && fast->next->next != NULL)
	{
		fast = fast->next->next;
		slow = slow->next;
	}

	List *leftHead = head, *rightHead = slow->next;
	slow->next = NULL;      //需要把左半链表的尾结点的next赋空值，所以用一个变量来记录右半链表的头

	leftHead  = listMergeSort(leftHead);
	rightHead = listMergeSort(rightHead);

	return mergeSortedList(leftHead, rightHead);
}
*/
void list_sort(List *list) {
	Node *p, *q, *e, *head, *tail;
	int insize, nmerges, psize, qsize, i;
	
	if (!list || !list->head)
		return;

	head = list->head;
	insize = 1;
	while (1) {
		p = head;
		head = NULL;
		tail = NULL;
		nmerges = 0;
		while (p) {
			nmerges++;
			q = p;
			psize = 0;
			for (i = 0; i < insize; i++) {
				psize++;
				q = q->next;
				if (!q) break;
			}
			qsize = insize;
			while (psize > 0 || (qsize > 0 && q)) {
				if (psize == 0) {
					e = q; q = q->next; qsize--;
				}
				else if (qsize == 0 || !q) {
					e = p; p = p->next; psize--;
				}
				else if (list->cmp(p->data, q->data) <= 0) {
					e = p; p = p->next; psize--;
				}
				else {
					e = q; q = q->next; qsize--;
				}
				if (tail) {
					tail->next = e;
				}
				else {
					head = e;
				}
				tail = e;
			}
			p = q;
		}
		tail->next = NULL;
		if (nmerges <= 1) {
			list->head = head;
			return;
		}
		insize *= 2;
	}
}
