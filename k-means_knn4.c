//created by Felipe Castrillon

#include<stdio.h>
#include<stdlib.h>
#include <unistd.h>
#define BUFSIZE 1024 //max line size
#include "list.h"
#include <math.h>

/*=============================================================================
---------------declare linked list Data
==============================================================================*/

struct Data {
        float *att;
		int kclus;
		struct Data *next;
};
typedef struct Data Node;

/*=============================================================================
---------------addData- adds data to nodes
==============================================================================*/


Node *addData(int cluster, float *temp){
       Node *tp = (Node *)malloc(sizeof(Node));  
       tp->kclus = cluster;
       tp->att = temp;
       tp->next = NULL;
       return tp;
}


/*=============================================================================
---------------append node to end of list
==============================================================================*/


void append(Node *head, Node *node){
     
     Node *tmp1 = head; 
     //Node *tmp2 = node;
     //printf("node %lf, %lf, %d\n", head->col_one,head->col_two,head->kclus);
     
     if(head->next == NULL){
        //printf("head is null\n");
        head->next = node; 
     }else{
          //printf("head is not null\n"); 
          while(tmp1->next){
              tmp1 = tmp1->next;
              
          }
          tmp1->next = node;
     }     
     //printf("get out of here\n");

}


/*=============================================================================
---------------k-means algorithm
==============================================================================*/


int kmeans(Node *head, float *means, int kclus, float *last_rmse, int attributes,
                float *att_sum, float *sqerr_sum, float *rmse, float *rmsediff, int *cnt){
    
    /*--------------------------------------------------------------------------
    ------------- set starting values
    --------------------------------------------------------------------------*/
    
    int i,j;//counters
    for (i=0;i<kclus;i++){//set starting values
        sqerr_sum[i] = 0;
        rmse[i] = 0;
        cnt [i] = 0;  
        rmsediff[i]=0;
        for(j=0;j<attributes;j++){
            att_sum[attributes*i+j] = 0;
        }
    }
    
    /*--------------------------------------------------------------------------
    -------------get distance and closest cluster for each index
    --------------------------------------------------------------------------*/

    float dist,min_dist;
    Node *cnode = head;
    do{
        min_dist = 10000;
        cnode = cnode->next;   
        int clust = 3;
        for (i=0;i<kclus;i++){
            dist=0;
            //printf("cluster is %d:", i); 
            for (j=0;j<attributes;j++){
               dist += pow(cnode->att[j] - means[(attributes)*(i)+j],2);
               //printf("node att %f, means %f\n", cnode->att[j],means[(attributes)*(i)+j]);    
            }    
            dist = pow(dist,.5);   
            //printf("dist %f\n", dist);           
            if (dist < min_dist){
               min_dist = dist;
               clust = i;  
               cnode->kclus = clust ;  
                 
            }
                        
        }
        //printf ("chosen clust %d\n",clust); 
        
        cnt[clust] +=1; //count elements on each cluster
        sqerr_sum[clust] += min_dist*min_dist;//sq err sum for each cluster
        

        for (j=0;j<attributes;j++){
            att_sum[clust*attributes + j] += cnode->att[j];
            //printf("\natt %d, clus %d,attsum %f, cnt %d, sqerr %f\n",j,clust,att_sum[clust*attributes + j],cnt[clust],sqerr_sum[clust]);
        }
        
        
       //printf("chosen cluster is %d, sum %f %f, cnt %d\n",clust, x_sum[clust], y_sum[clust], cnt[clust]); getchar();
    
    }while(cnode->next != NULL);
    
    /*--------------------------------------------------------------------------
    -------------get new centroids and calculate rmse
    --------------------------------------------------------------------------*/
    
    int diff_test = 1;
    for (i=0;i<kclus;i++){    
        if (cnt[i] != 0){   
           rmse[i] = pow(sqerr_sum[i]/cnt[i],.5);
        }else{
              rmse[i] = 0;
        }
        //printf("clus %i, err %f, cnt %d\n", i,sqerr_sum[i],cnt[i]);
        if (last_rmse != NULL){
           rmsediff[i] = last_rmse[i] - rmse[i]; //compare last to current rmse
           if (fabs(rmsediff[i]) > 0.0001){//test to see if rmse is sufficiently small
              diff_test = 0;
           }
        }
        else {
             diff_test = 0;
        }
        last_rmse[i] = rmse[i];
    }   
    
    
    
    for (i=0;i<kclus;i++){
        
        printf("new center for cluster %d:\n",i);
        //printf("iteration is %d, means is %f %f\n", i, means[i]); 
        for (j=0;j<attributes;j++){
            //printf("attsum %f, cnt %f\n", att_sum[(attributes)*(i)+j], cnt[i]);
            if(cnt[i] != 0){
                      means[(attributes)*(i)+j] = att_sum[(attributes)*(i)+j]/(float)cnt[i];
            }else{
                  means[(attributes)*(i)+j] = 0;   
            }
            printf("attribute %d: %f\n", j ,means[(attributes)*(i)+j]);
        }
        printf("cluster %d, counts %d, rmse %f rmse diff from last iter %f\n\n", i ,cnt[i],rmse[i],rmsediff[i]);
    }
    
   
    if (diff_test == 0){
       return 0;
    } else{
           return 1;
    }  
    
    
}



/*=============================================================================
---------------k-nearest neighbor algorithm (only works for 2 attributes)
==============================================================================*/



int knn(Node *head, int nclus){
    
    /*--------------------------------------------------------------------------
    -------------declare and define starting variables
    --------------------------------------------------------------------------*/
    
    printf ("starting k-nearest neighbor algorithm\n");
    float val1, val2; //values to be tested
    int knum = 5; //number of k neighbors
    float tempdist; //temporary distance for each node
    float *dist = (float *)malloc(knum*sizeof(float)); //distances of array
    int *clussum = (int *)malloc(knum*sizeof(int)); //number of nearest neighbors for each item
    int *cluster = (int *)malloc(knum*sizeof(int)); //cluster for each neighbor
    float *att1 = (float *)malloc(knum*sizeof(float)); //values for neighbor
    float *att2 = (float *)malloc(knum*sizeof(float)); 
    int i=0, l=0; //default counters
    Node *cnode = head;
    
    //set starting values
    for (i=1;i<=knum;i++){
        dist[i] = -1;
        clussum[i] = 0;
        cluster[i] = 0;  
    }
        
    
    /*--------------------------------------------------------------------------
    -------------get target data items
    --------------------------------------------------------------------------*/
    
    printf("val 1?\n"); 
     if (scanf("%f", &val1) != 1){ //read in value and make sure it is valid number
               printf("error: number invalid!\n");
               exit(EXIT_FAILURE); 
    }
    printf("value 1 is %f\n", val1);
    printf("val 2?"); 
     if (scanf("%f", &val2) != 1){ //read in value and make sure it is valid number
               printf("error: number invalid!\n");
               exit(EXIT_FAILURE); 
    }
    printf("value 2 is %f\n", val2);
    
    
    /*--------------------------------------------------------------------------
    -------------start algorithm
    --------------------------------------------------------------------------*/
    
    int itemp;
    float ftemp;
    int testbreak;
    do{ //loop through all data
         cnode = cnode->next; //loop        
         tempdist = pow(pow(cnode->att[0] - val1,2) + pow(cnode->att[1] - val2,2),0.5); // calc distance    
         
         testbreak = 0;//test if break is needed
                  
         for (i=1;i<=knum;i++){ //get list of smallest distances filled
                                //and in order 
             if (dist[i] == -1 ){//if list not filled, then fill list
                dist[i] = tempdist;
                cluster[i] = cnode->kclus;
                att1[i] = cnode->att[0];
                att2[i] =  cnode->att[1];
                testbreak = 1;
                
                break;
             }
         }
         if (testbreak == 1) 
            continue;
         
         //for (i=1;i<=knum;i++){
           //printf("i %d, cola %f, colb %f dist %f\n", i,cola[i],colb[i],dist[i]);
           //} //getchar();
                     
         testbreak = 0;
         //for (i=1;i<=knum;i++){
         if(tempdist < dist[knum]){//test if new tempdist is smaller than values in dist
                    dist[knum] = tempdist; //replace values for i = largest value in dist
                    cluster[knum] = cnode->kclus;
                    att1[knum] = cnode->att[0];
                    att2[knum] =  cnode->att[1];
         }                 
              
         for (l=knum;l>1;l--){//order list by distances 
                //printf("iter %d\n", l);
                if (dist[l] == -1){ //if distance has not been set
                   //printf("continue\n");
                   continue;
                }
                if (dist[l] < dist[l-1]){ //switch to get in right order
                   //printf("switch %d and %d\n", l, l-1);
                    
                    ftemp = dist[l-1];
                    dist[l-1] = dist[l];
                    dist[l] = ftemp;
                    
                    itemp = cluster[l-1];
                    cluster[l-1] = cluster[l];
                    cluster[l] = itemp;
                    
                    ftemp = att1[l-1];
                    att1[l-1] = att1[l];
                    att1[l] = ftemp;
                    
                    ftemp = att2[l-1];
                    att2[l-1] = att2[l];
                    att2[l] = ftemp;

                 }
          }    
              
                   
              
    }while(cnode->next != NULL);
    
    for (i=1;i<=knum;i++){
        clussum[cluster[i]] += 1;    
    }
    int max = 0;
    int knnc = 0; //knn-cluster
    for (i=1;i<=knum;i++){
        printf("neighbor %d, values %f %f, dist %f, cluster %d\n", i, att1[i], att2[i],dist[i],cluster[i]);
        
        if (clussum[i] == 0)
           continue;
        if (max == 0 || clussum[i] < max){
           knnc = cluster[i];
           max = clussum[i];
        }      
    }  
    printf("knn is %d", knnc);  
    
    free(dist); 
    free(clussum); 
    free(cluster); 
    free(att1); 
    free(att2);
    return 1;
}  


/*=============================================================================
---------------main
==============================================================================*/


int main (){
    char file[50];//name of file
    int nclus;//number of clusters
    //float *Cone, *Ctwo, *temp, *temp2;
    
    
    /*--------------------------------------------------------------------------
    -------------get input parameters
    --------------------------------------------------------------------------*/
    
    
    //get name of file to open
    printf("what is the name of your file?");    
    if (scanf("%s", &file[0]) != 1){ //read in value and make sure it is valid name
               printf("error: filemane invalid!\n");
               exit(EXIT_FAILURE); 
    }
    printf("filename is %s\n", file);
    
    //get number of clusters
    printf("how many clusters do you want to find?"); 
     if (scanf("%d", &nclus) != 1){ //read in value and make sure it is valid number
               printf("error: number invalid!\n");
               exit(EXIT_FAILURE); 
    }
    printf("number of k clusters is %d\n", nclus);
        
        
    /*--------------------------------------------------------------------------
    -------------open and read input file
    --------------------------------------------------------------------------*/
        
    //open input file 
    FILE *fp;
    if ((fp = fopen(file,"r")) == NULL){
       printf("Error:cannot read the file %s\n", file);        
    }else{
       printf("success opening file %s\n", file); 
    }
    
    //read number of items and attributes
    int items, attributes;
    fscanf(fp,"%d",&items);
    fscanf(fp,"%d",&attributes);
    printf ("items %d attributes %d\n", items, attributes);
    getchar();
    
    //Create nodes    
    Node *head = (Node *) malloc (sizeof(Node)); 
    Node *node = (Node *) malloc (sizeof(Node));
    
    float *temp  = (float *)malloc(items*attributes*sizeof(float));
    //read file
    printf("reading file\n");
    int i, j; //counters
    for (i=0; i < items; i++){//for each item
        printf("reading item %d\n",i);      
          for (j=0; j < attributes; j++){
              fscanf(fp,"%f",&temp[(attributes)*i+j]);              
              //printf("read attribute %d, %f\n",j,temp[(attributes)*i+j]);
              
          }     

          node = addData(nclus+1,&temp[(attributes)*i]);  
          append(head,node);
    }
       
       
    /*--------------------------------------------------------------------------
    -------------get initial means
    --------------------------------------------------------------------------*/
        
    float *means = (float *)malloc(nclus*attributes*sizeof(float));
    //float *means_b = (float *)malloc(nclus*attributes*sizeof(float)); 
    
    int index;
    node = head;
    //printf("node...%f %f\n", node->att[0], node->att[1]); 
    for (i=0;i<nclus;i++){
        node = node->next;
        //printf("node...%f %f\n", node->att[0], node->att[1]); 
        printf("\n\ninitial center for cluster %d:\n",i);
        for (j=0;j<attributes;j++){ 
            index =  (attributes)*(i)+j;    
            means[index] = node->att[j];       
            printf("attribute %d: %f\n", j, node->att[j]);    
        }
    }
    printf("press enter\n");
    getchar();
    
   /*--------------------------------------------------------------------------
    -------------run k-means algorithm
    --------------------------------------------------------------------------*/
   
   //declare variables for k-means algo
   float *att_sum = (float *)malloc(nclus*attributes*sizeof(float)); //sum of each column for each cluster
   float *sqerr_sum = (float *)malloc(nclus*sizeof(float)); //keep track of rmse
   float *rmse = (float *)malloc(nclus*sizeof(float)); //keep track of rmse
   float *rmsediff = (float *)malloc(nclus*sizeof(float)); //differences in rmse
   int *cnt = (int *)malloc(nclus*sizeof(int)); //number of elements in each cluster
   
   
   float *last_rmse = (float *)malloc(nclus*sizeof(float));  
   int test = 0;
   int count = 0;
   do{//repeat until test = 1 (RMSE converges)
      test = kmeans(head,&means[0],nclus,&last_rmse[0], attributes,
           &att_sum[0], &sqerr_sum[0], &rmse[0], &rmsediff[0], &cnt[0]);
      printf("press enter for next iteration\n");
      getchar();
      //printf("test %d\n", test); 
      count ++;
      if (count > 100) 
         test = 1; 
   }while(test == 0);
          
    
   /*--------------------------------------------------------------------------
    -------------print results in output file
    --------------------------------------------------------------------------*/
    FILE *out;
    if ((out = fopen("output_gbm.csv","w")) == NULL){
       printf("Error:cannot read the file output.txt\n");        
    }else{
       printf("success opening output file output.txt for output\n"); 
    }
    node = head;
    while(node->next != NULL){
        node = node->next;      

         for (j=0;j<attributes;j++) { 
             fprintf(out, "%f,", node->att[j]);
         }
         fprintf(out,"%d\n",node->kclus);
    }
    fclose(out);
    fclose (fp) ; 
    
    
    /*--------------------------------------------------------------------------
    -------------run k-nearest neighbor algorithm
    --------------------------------------------------------------------------*/
    int response;
    printf("if you want to run k-nearest neighbor type the number 1 and press enter, otherwise press any other key: "); 
     if (scanf("%d", &response) != 1){ //read in value and make sure it is valid number
               printf("error: number invalid!\n");
               exit(EXIT_FAILURE); 
    }
    printf("response is %d\n", response);
    //run k-nn algorithm
    
    if (response == 1){
         knn(head, nclus); //knn call function
    }
    
    
    free(last_rmse);
    free(means);
    free(att_sum); 
    free(sqerr_sum); 
    free(rmse); 
    free(rmsediff); 
    free(cnt);
    

    return 1; 
}




