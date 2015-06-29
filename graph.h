#include <stdio.h>

#ifndef _GRAPH_H
#define _GRAPH_H

#define GRAPH_EMPTY    -1 /* keys cannot equal HASH_EMPTY */
#define GRAPH_NOTFOUND -1 /* data items cannot equal HASH_NOTFOUND */


/*------------------------------------------------------------------------------
-------------Define structure of a Graph abstract data type
-------------"links" array saves matrix array with each cell representing a
-------------connection from vertex i (row) to vertex j (col). Array is saved as 1-D
------------------------------------------------------------------------------*/

typedef struct
{
    float *links;  /* list of keys, used for rehashing */
}
Graph;

/*------------------------------------------------------------------------------
-------------Declare a Graph abstract data type and initialize to default
------------------------------------------------------------------------------*/

Graph *GraphCreate(int nvert);
void GraphGetLinks(int nvert, int parent, int child, float cost, Graph *h);
void GraphDestroy(Graph *h);
void PrintGraph(Graph *h, int nvert);

#endif /* _GRAPH_H */
