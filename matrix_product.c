#include<stdio.h>
#include<stdlib.h>

int matrix_multiply(int n1, int n2, int n3, double *a, double *b, double *c){
    
    //printf ("arguments passed %d %d %d\n", n1, n2, n3);
    //printf ("value of pointers %f %f %f\n", *a, *b, *c);
    double *at; //temporary pointers
    double *bt;
    double *ct;
    int c_size = n1*n3;
    
    int sum = 0;
    int i = 0, j = 0, k = 0;
    for (i=0; i < n1; i++){ //num a
        for (j=0; j < n2; j++){ //num 
            for (k=0; k < n3; k++){ //num b cols
                at = a + n2*i + j;
                bt = b + n3*j + k;
                ct = c + n3*i + k;            
                *ct = *ct + *at * *bt;
            }
        }
    }    
    
    printf ("\n final C matrix: ");
    double *tt;
    tt = c;
    
    matrix_print(n1,n3,&ct[0]);
    
 return 0; 
    
}


int matrix_fill_random(int n1, int n2, double *a){
       int as = n1*n2;
       int i;       
       for (i=0;i < as;i++){
           a[i] = 2*rand()/(double)RAND_MAX; 
           printf("i %d: %f \n",i,a[i]);
       }
       return 0;
}

int matrix_print(int n1, int n2, double *a){
       printf("\n matrix form: \n");
       int as = n1*n2;
       int i;       
       for (i=0;i < as;i++){
           printf("%f  ",a[i]);
           if ((i+1) % n1 == 0){
              printf ("\n");
           }
           
       }
       return 0;
}



int main(void){
       
       int n1, n2, n3; /*declare matrix size parameters*/ 
       int as, bs, cs; // declare size of array
       double *am, *bm, *cm; //declare and define pointers
       
       printf("choose n1\n");
       scanf("%d", &n1);
       getchar();
       printf("choose n2\n");
       scanf("%d", &n2);
       getchar();
       printf("choose n3\n");
       scanf("%d", &n3);
       getchar();

       as = n1 * n2; //define size of arrays
       bs = n2 * n3; 
       cs = n1 * n3;
       
       printf("array size of matrix a, b, c: %d %d %d\n", as, bs, cs); 
       
       am = (double *) malloc(as*sizeof(double));
       bm = (double *) malloc(bs*sizeof(double));
       cm = (double *) malloc(cs*sizeof(double));
       
       /*set value of a matrix*/
       printf("\nget values of matrix A and press enter after every value: \n");
       matrix_fill_random(n1, n2, &am[0]);
       matrix_print(n1,n2,&am[0]);
       
       
       /*set value of b matrix*/
       printf("\nget values of matrix B and press enter after every value: \n");
       int j;       
       for (j=0;j < bs;j++){
           printf("j %d: ",j);
           scanf("%lf", &bm[j]);
           getchar();
       }
       matrix_print(n2,n3,&bm[0]);
       
       int k;       
       for (k=0;k < cs;k++){
           cm[k] = 0;
           //printf("set %d %d to zero\n", k, &cm[k]);
       }
       
           
       /*send values to matrix multiplication*/
       
       matrix_multiply(n1,n2,n3,&am[0],&bm[0],&cm[0]);
       
       getchar();
}
