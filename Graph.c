#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "graph.h"
#include "queue.h"


int main (int argc, char *argv[]){
    
   /*--------------------------------------------------------------------------
    -------------read and save arguments
    --------------------------------------------------------------------------*/ 
    
   if (argc != 4){
        printf("incorrect number of arguments!\n"); getchar();
        exit(EXIT_FAILURE);
   }
   int source = atoi(argv[2]);
   int sink = atoi(argv[3]);
   
    /*--------------------------------------------------------------------------
    -------------open graph file and read initial parameters
    --------------------------------------------------------------------------*/
    
    //open input file 
    FILE *fp;
    if ((fp = fopen(argv[1],"r")) == NULL){
       printf("Error:cannot read the file %s\n", argv[1]); 
       return 0;     
    }else{
       printf("success opening file %s\n", argv[1]); 
    }
    
    //read number of vertices, and links
    int nvert;
    int nlinks;
    fscanf(fp,"%d",&nvert);
    fscanf(fp,"%d",&nlinks);
    //printf ("nvert %d nlinks %d\n", nvert, nlinks);
    
    /*--------------------------------------------------------------------------
    -------------read input file, save to graph structure, and print graph
    --------------------------------------------------------------------------*/
    
    int parent;
    int child;
    int j;
    float cost;
    
    Graph *g = GraphCreate(nvert);
    for (j=0; j < nlinks; j++){//read file for each column
          fscanf(fp,"%d",&parent);              
          fscanf(fp,"%d",&child); 
          fscanf(fp,"%f",&cost);
          //printf ("line %d, %d, %f\n", parent,child,cost); 
          GraphGetLinks(nvert, parent, child, cost, g);
    }

    PrintGraph(g,nvert); 
    
    /*--------------------------------------------------------------------------
    -------------create priority queue
    --------------------------------------------------------------------------*/
    
    int *parent_p = (int*)malloc (sizeof(int));//parent pointer
    int *child_p = (int*)malloc (sizeof(int));//child pointer
     
    Queue *q = QueueCreate(nvert);
    printf ("queue structure created\n");
    
    /*--------------------------------------------------------------------------
    -------------Run Djikstra's algorithm
    --------------------------------------------------------------------------*/
    
    printf("begin Djikstra's algorithm from source %d to sink %d\n",source,sink);
    q->def_vert[source-1] = 1; 
    do{
        QueueGetMin(q, g, nvert,parent_p,child_p);
        Update_Queue(q, g,nvert,parent_p,child_p);
        if (*child_p == sink)
            break;
    }while (1);
    
    /*--------------------------------------------------------------------------
    -------------Print results
    --------------------------------------------------------------------------*/
    
    printf("\n\nResults: the distance from source %d to sink %d is %f\n", source,sink,q->dist[sink-1]); 
    getchar();
    
    GraphDestroy(g);
    QueueDestroy(q);
    
    printf("End of program, press enter to exit\n"); getchar();
    return 0;
}



