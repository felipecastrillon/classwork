//created by Felipe Castrillon on December 12th, 2012

#include <stdio.h>
#include<stdlib.h>

/*----------------------------------------------------------------------
print the tic tac toe grid   
-----------------------------------------------------------------------*/
void print (int *square){
     printf ("printing graph...\n\n\n\n");
     char symbol;
     int i, j; //counters
     int loc;  //temp location index
     for (i = 0; i< 3; i++){ //rows
         for (j= 0; j< 3; j++) { //columns
             loc = (i)*3 + j;
             if (square[loc] == 0)
                symbol = ' '; 
             if (square[loc] == 1)
                symbol = 'O'; 
             if (square[loc] == 2)
                symbol = 'X'; 
             printf (" %c ",symbol);
             if (j == 0 || j == 1)
                printf ("|");// col separator    
         }
         if (i == 0 || i == 1){
            printf ("\n---|---|---\n");//row separator
         }else{
            printf("\n");   
         }
     }
     
     printf("\n\n\n");
}

/*----------------------------------------------------------------------
computer chooses a random location in the tic tac toe grid   
-----------------------------------------------------------------------*/
void comp_move (int *square){
    int move = 9 * rand()/(double)RAND_MAX;
    //printf ("move in %d\n", move);
    while (square[move] != 0){
       move ++;
       if (move == 9){ //loop around
          move = 0;
       }                 
    }
    printf ("computer move in %d ... \n", move+1);
    square[move] = 1;
}

/*----------------------------------------------------------------------
user is asked for their move  
-----------------------------------------------------------------------*/
void person_move (int *square){    
    int loc;
    do{
        printf("choose a grid from 1 to 9 to make your move (remember: 1 -top left, 9- bottom right):\n");
        scanf("%d", &loc);//get user prompt
        getchar();
        if (loc > 9 || loc < 1){ //make sure location is not out of range
           printf("location out of range\n");
           continue;
        }
        loc = loc-1;
        if (square[loc] == 1 || square[loc] == 2){ //make sure location has 
                                                   //not been assigned
           printf("that grid has already been assigned\n");
        }
        
        else{break;}
    } while(1);
    square[loc] = 2;
}


/*----------------------------------------------------------------------
check if there is a winner and return winner
-----------------------------------------------------------------------*/
int win(int *square){
    int winner = 0;
    
    //8 possible ways of winning
    if (square[0] == square[1] && square[1] == square[2] && square[0] != 0) 
       winner = square[0];
    if (square[3] == square[4] && square[4] == square[5] && square[3] != 0) 
       winner = square[3];
    if (square[6] == square[7] && square[7] == square[8] && square[6] != 0) 
       winner = square[6];
    if (square[0] == square[3] && square[3] == square[6] && square[0] != 0) 
       winner = square[0];
    if (square[1] == square[4] && square[4] == square[7] && square[1] != 0) 
       winner = square[1];
    if (square[2] == square[5] && square[5] == square[8] && square[2] != 0) 
       winner = square[2];
    if (square[0] == square[4] && square[4] == square[8] && square[0] != 0) 
       winner = square[0];
    if (square[2] == square[4] && square[4] == square[6] && square[2] != 0) 
       winner = square[2];
       
    return winner;
}


int main (){
       printf ("welcome to tic tac toe!\n");
       int move;
       printf ("choose who starts:\n(1)computer \n(2)you\n\n");
       scanf("%d", &move);
       if (move > 2 || move < 1){
          printf ("Not an option!\n"); getchar();
          exit(EXIT_FAILURE);         
       }
       /*----------------------------------------------------------------------
       create grid and initialize to zero    
       -----------------------------------------------------------------------*/
       
       int *tictac = (int *)malloc(9*sizeof(int));
       int i; //counter
       int winner = 0;
       for (i=0;i<9;i++){
           tictac[i] = 0;    
       }
       printf("initial grid:\n");
       print (&tictac[0]);
       
       int cnt = 0;
       
       /*----------------------------------------------------------------------
       start playing   
       -----------------------------------------------------------------------*/
       do{
                                                                                                                                                
           //computer move  
           if(move == 1){
               comp_move (&tictac[0]);
               cnt ++;
               winner = win(&tictac[0]);
               print (&tictac[0]);
               if (winner)
                  continue;
               move =2;
           } 
           //person move
           if (move == 2){
               person_move(&tictac[0]);
               cnt++;
               winner = win(&tictac[0]);
               print (&tictac[0]);
               if (winner)
                  continue;
               move = 1;
           }
           
       }while(!winner && cnt < 9); //loop while there is not a winner and while  
                                   //all cells have not been filled
       
       
       /*----------------------------------------------------------------------
       print winner   
       -----------------------------------------------------------------------*/
       if (winner == 0){
              printf ("\n\n\n\n\nTie, no winner!!!!!!\n");
       } else if (winner == 1){
              printf ("\n\n\n\n\nWinner is COMPUTER, you need practice!!!!!!!\n");
       }else if (winner == 2){
              printf ("\n\n\n\n\nWinner is YOU!!!!!!!\n");
       }
       printf("press enter to exit\n"); getchar();  
       
       return 0;
}
