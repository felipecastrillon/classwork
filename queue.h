#include <stdio.h>
#include "graph.h"
#ifndef _QUEUE_H
#define _QUEUE_H

#define QUEUE_EMPTY    -1 /* keys cannot equal HASH_EMPTY */
#define QUEUE_NOTFOUND -1 /* data items cannot equal HASH_NOTFOUND */

/*------------------------------------------------------------------------------
-------------Define structure of a Queue abstract data type
-------------def_vert stores array of vertices that have been defined by algorithm
-------------dist stores array of the distance from source vertex to other vertices
------------------------------------------------------------------------------*/

typedef struct
{
    int *undef_vert;
    int *def_vert;
    float *dist;
}
Queue;

Queue* QueueCreate(int nvert);
void QueueDestroy(Queue *q);
int QueueGetMin(Queue *p, Graph *g,int nvert, int *parent,int *child);
void Update_Queue(Queue *q, Graph *g,int nvert, int *parent_p,int *child_p);
void Print_Queue(Queue *q, int nvert);
void Start_Queue(Queue *q,int nvert,int source);
#endif /* _QUEUE_H */
