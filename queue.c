#include <stdio.h>
#include "queue.h"

/*------------------------------------------------------------------------------
-------------Declare a Queue data type
------------------------------------------------------------------------------*/
Queue* QueueCreate(int nvert){
    Queue *q = (Queue*) malloc(sizeof(Queue));
    float *p;
    int *o,*n,i;
    
    q->undef_vert = (int*) malloc(nvert * sizeof(int));
    q->def_vert = (int*) malloc(nvert * sizeof(int));  /* vertices with undefined dist */
    q->dist = (float*) malloc(nvert * sizeof(float));   /* distance between source and vertex */
    /* Initialize undefind vertices */
    /* Initialize the table to empty */
    p = q->dist;
    o = q->def_vert;
    n = q->undef_vert;
    for (i=0; i<nvert; i++){
        *p++ = 0;
        *n++ = i+1;
        *o++ = QUEUE_EMPTY;
    }
        
    return q;
}

/*------------------------------------------------------------------------------
-------------Frees memory of Queue abstract data type
------------------------------------------------------------------------------*/

void QueueDestroy(Queue *q){
    free(q->def_vert);
    free(q->dist);
    free(q);
}

/*------------------------------------------------------------------------------
-------------Gets minmum distance from defined vertice set to undefined vertices
------------------------------------------------------------------------------*/

int QueueGetMin(Queue *p, Graph *g,int nvert, int *parent,int *child){
    //printf("start Queue Get Min\n");
    int i,j, loc, undefined, defined;
    float cost, dist;
    float min_dist =-1;
    for (j=1; j<=nvert; j++){//for every defined vertice
        defined = p->def_vert[j-1];
        //printf("defined j, %d, %d\n",j,p->def_vert[j-1]);
        if (defined == -1)//next if distance not defined
            break;
            
        for (i=1; i<=nvert; i++){//for every undefined vertice
            //printf("undefined i, %d, %d\n",i,p->undef_vert[i-1]);
            undefined = p->undef_vert[i-1];
            loc = (defined-1)*nvert+(undefined-1);
            cost = g->links[loc];
            
            if (undefined == -1)//break if end of list
                break;            
            if (cost > -1 && defined != undefined){//if undefined and if it is neighbor to j  
                //printf("test,undef %d, loc %d, cost %f\n",undefined,loc,cost);                                                      
                dist = p->dist[defined-1] + cost;
                //printf("test,undef %d, loc %d, cost %f, dist %f\n",undefined,loc,cost,dist); 
                if ((dist < min_dist) || (min_dist ==-1)){
                    min_dist = dist;
                    *child = undefined;
                    *parent = defined;
                    //printf ("min dist %f, child %d, parent %d\n", min_dist, *child,*parent);
                }            
            }
        }
    }
    if (min_dist == -1)
        return 0;
    else 
        return 1;
}

/*------------------------------------------------------------------------------
-------------Starts the queue
------------------------------------------------------------------------------*/

void Start_Queue(Queue *q, int nvert,int source){
    q->def_vert[0] = source;
    int i;
    //move source from undefined to defined list and update lists
     for (i=0; i<nvert; i++){
        if (q-> undef_vert[i] == source){
            int j;
            for (j=i; j<=nvert;j++){
                if(j==nvert-1){
                    q->undef_vert[j] = -1;  //last is undefined
                    break; 
                } 
                q->undef_vert[j] = q->undef_vert[j+1];
            }
        }
    }
}

/*------------------------------------------------------------------------------
-------------Updates new distances and defined/undefined vertex
------------------------------------------------------------------------------*/

void Update_Queue(Queue *q, Graph *g,int nvert, int *parent_p,int *child_p){
    int parent = *parent_p;
    int child = *child_p;
    int skip = 0;
    int i;
    //update distances
    q->dist[child-1] = q->dist[parent-1] +  g->links[(parent-1)*nvert+(child-1)];
     //move child node from undefined to defined list and update lists
    for (i=0; i<nvert; i++){//remove child node and update list
        if (q-> undef_vert[i] == child){
            int j;
            for (j=i; j<=nvert;j++){
                if (q->undef_vert[j] == -1) //if end of list is reached
                    break;
                if(j==nvert){
                    q->undef_vert[j] = -1; //if last element is reached 
                    break; 
                } 
                q->undef_vert[j] = q->undef_vert[j+1];//push lists to beginning
            }
        }
        //move child node to defined list
        if (skip == 0){
            if ((q->def_vert[i] == -1)){
                q->def_vert[i] = child;
                skip = 1;
            }
        }
    }    
}

/*------------------------------------------------------------------------------
-------------Print Queue
------------------------------------------------------------------------------*/

void Print_Queue(Queue *q, int nvert){
    int i;
    printf ("printing queue....\n\n\n");
    
    printf ("def_ver:\n");
    for (i=0; i<nvert; i++){
        printf("%d, ", q->def_vert[i]);
    } 
    printf ("\n\n");
    printf ("undef_ver:\n");
    for (i=0; i<nvert; i++){
        printf("%d, ", q->undef_vert[i]);
    } 
    printf ("done printing queue....\n\n\n");
    
    
}



