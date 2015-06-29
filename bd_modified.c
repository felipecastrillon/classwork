//created by Felipe Castrillon

#include<stdio.h>
#include<stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

struct box
{
    int head;
};

// it is possible to use smaller boxes and more complex neighbor patterns
#define NUM_BOX_NEIGHBORS 13
int box_neighbors[NUM_BOX_NEIGHBORS][3] =
{
    {-1,-1,-1},
    {-1,-1, 0},
    {-1,-1,+1},
    {-1, 0,-1},
    {-1, 0, 0},
    {-1, 0,+1},
    {-1,+1,-1},
    {-1,+1, 0},
    {-1,+1,+1},
    { 0,-1,-1},
    { 0,-1, 0},
    { 0,-1,+1},
    { 0, 0,-1}
};


/*=============================================================================
---------------interactions
==============================================================================*/
//interactions(npos,L,4,4,distances2,pairs,1000,&numpairs,xpos,ypos,zpos);
int interactions(int npos, double L, int boxdim,double cutoff2, double *distances2, 
                 int *pairs, int maxnumpairs, int *numpairs_p,double *xpos, double *ypos, double *zpos)
    {
    
    double x_unitv,y_unitv,z_unitv; // coordinates of unit vector
    double mag; //magnitude of vector
    
    if (boxdim < 4 || cutoff2 > (L/boxdim)*(L/boxdim))
    {
        printf("interactions: bad input parameters\n");
        return 1;
    }

    struct box b[boxdim][boxdim][boxdim];
    struct box *bp;
    struct box *neigh_bp;

    // box indices
    int idx, idy, idz;
    int neigh_idx, neigh_idy, neigh_idz;

    // allocate memory for particles in each box
    for (idx=0; idx<boxdim; idx++)
    for (idy=0; idy<boxdim; idy++)
    for (idz=0; idz<boxdim; idz++)
        b[idx][idy][idz].head = -1;

    // allocate implied linked list
    int *next = (int *) malloc(npos*sizeof(int));
    if (next == NULL)
    {
        printf("interactions: could not malloc array for %d particles\n", npos);
        return 1;
    }
    
    
    // traverse all particles and assign to boxes
    int i;
    for (i=0; i<npos; i++)
    {
        // initialize entry of implied linked list
        next[i] = -1;

        // which box does the particle belong to?
        // assumes particles have positions within [0,L]^3
        idx = (int)(xpos[i]/L*boxdim);
        idy = (int)(ypos[i]/L*boxdim);
        idz = (int)(zpos[i]/L*boxdim);

        // add to beginning of implied linked list
        bp = &b[idx][idy][idz];
        next[i] = bp->head;
        bp->head = i;
    }
    
    int numpairs = 0;
    int p1, p2;
    double d2, dx, dy, dz;
    double dist;
    double disp[npos*3];//array of displacements
    
    //set displacements equal to zero
    for (i=0;i<npos*3;i++){
        disp[i]=0;   
    }
    
    for (idx=0; idx<boxdim; idx++)
    {
        for (idy=0; idy<boxdim; idy++)
        {
            for (idz=0; idz<boxdim; idz++)
            {
                bp = &b[idx][idy][idz];

                // within box interactions
                p1 = bp->head;
                while (p1 != -1)
                {
                    p2 = next[p1];
                    //printf("p1 %d, p2 %d\n", p1,p2);
                    while (p2 != -1)
                    {
                        if (numpairs >= maxnumpairs)
                            return -1;

                        // do not need minimum image since we are in same box
                        dx = xpos[p1] - xpos[p2];
                        dy = ypos[p1] - ypos[p2];
                        dz = zpos[p1] - zpos[p2];

                        if ((d2 = dx*dx+dy*dy+dz*dz) < cutoff2)
                        {
                            dist = pow(d2,0.5);
                            //printf("dist %lf\n", dist); getchar();
                            pairs[2*numpairs]   = p1 + 1;
                            pairs[2*numpairs+1] = p2 + 1;
                            distances2[numpairs] = dist;
                            numpairs++;
                            
                            x_unitv = (xpos[p1]-(xpos[p2]))/dist;
                            y_unitv = (ypos[p1]-(ypos[p2]))/dist;
                            z_unitv = (zpos[p1]-(zpos[p2]))/dist;
                            mag = 125*(2-dist);
                        
                            disp[3*p1+0] += (x_unitv * mag)*0.002;//displacement after forces are applied
                            disp[3*p1+1]  += (y_unitv * mag)*0.002;
                            disp[3*p1+2]  += (z_unitv * mag)*0.002;
                            
                            //printf("i %d, index %d,movements %f,%f,%f disp %lf,%lf,%lf\n", 
                            //p1, 3*p1,(xv * mag)*0.002,(yv * mag)*0.002,(zv * mag)*0.002,disp[3*p1+0],disp[3*p1+1],disp[3*p1+2]); getchar();
                            
                        }

                        p2 = next[p2];
                    }                    
                    p1 = next[p1];
                }
                
                // interactions with other boxes
                int j;
                for (j=0; j<NUM_BOX_NEIGHBORS; j++)
                {
                    neigh_idx = (idx + box_neighbors[j][0] + boxdim) % boxdim;
                    neigh_idy = (idy + box_neighbors[j][1] + boxdim) % boxdim;
                    neigh_idz = (idz + box_neighbors[j][2] + boxdim) % boxdim;

                    neigh_bp = &b[neigh_idx][neigh_idy][neigh_idz];

                    // when using boxes, the minimum image computation is 
                    // known beforehand, thus we can  compute position offsets 
                    // to compensate for wraparound when computing distances
                    double xoffset = 0.;
                    double yoffset = 0.;
                    double zoffset = 0.;
                    if (idx + box_neighbors[j][0] == -1)     xoffset = -L;
                    if (idy + box_neighbors[j][1] == -1)     yoffset = -L;
                    if (idz + box_neighbors[j][2] == -1)     zoffset = -L;
                    if (idx + box_neighbors[j][0] == boxdim) xoffset =  L;
                    if (idy + box_neighbors[j][1] == boxdim) yoffset =  L;
                    if (idz + box_neighbors[j][2] == boxdim) zoffset =  L;

                    p1 = neigh_bp->head;
                    while (p1 != -1)
                    {
                        //xf =0; yf =0; zf=0;
                        p2 = bp->head;
                        while (p2 != -1)
                        {
                            if (numpairs >= maxnumpairs)
                                return -1;
                            
                            // compute distance vector
                            dx = xpos[p1] - xpos[p2] + xoffset;
                            dy = ypos[p1] - ypos[p2] + yoffset;
                            dz = zpos[p1] - zpos[p2] + zoffset;

                            if ((d2 = dx*dx+dy*dy+dz*dz) < cutoff2)
                            {
                                
                                dist = pow(d2,0.5);
                                pairs[2*numpairs]   = p1 + 1;
                                pairs[2*numpairs+1] = p2 + 1;
                                distances2[numpairs] = d2;
                                numpairs++;
                                
                                if (dist == 0 ) printf ("dist %f\n",dist);
                                                                
                                x_unitv = (xpos[p1]-xpos[p2]+xoffset)/dist;
                                y_unitv = (ypos[p1]-ypos[p2]+yoffset)/dist;
                                z_unitv = (zpos[p1]-zpos[p2]+zoffset)/dist;
                                mag = 125*(2-dist);
                                
                                //printf("p1 %d, index %d,disp %lf,%lf, %lf\n", p1, 3*p1,disp[3*p1+0],disp[3*p1+1],disp[3*p1+2]);
                                
                                disp[3*p1+0] += (x_unitv * mag)*0.002;//displacement after forces are applied
                                disp[3*p1+1] += (y_unitv * mag)*0.002;
                                disp[3*p1+2] += (z_unitv * mag)*0.002;
                                
                                //printf("p1 %d, index %d,disp %lf,%lf, %lf\n", p1, 3*p1,disp[3*p1+0],disp[3*p1+1],disp[3*p1+2]); 
                                
                            }
                            p2 = next[p2];
                        }
                        
                        p1 = next[p1];
                    }
                }
                
            }
        }
    }
    
    /*for (i=0;i<npos*3;i++){
        printf("i %d, disp %lf\n",i,disp[i]);   
    }
    getchar();*/
    
    
    //update positions
    for(i=0;i<npos;i++){
            //if (i == 0) printf("posa %lf %lf %lf,\n",xpos[i],ypos[i],zpos[i]);
            xpos[i] += disp[3*i + 0];
            ypos[i] += disp[3*i + 1];
            zpos[i] += disp[3*i + 2];  

            //update positions if they are outside of the box
            if (xpos[i] >= L)
                xpos[i] -= L;
            if (xpos[i] < 0)
                xpos[i] += L;
            if (ypos[i] >= L)
                ypos[i] -= L;
            if (ypos[i] < 0)
                ypos[i] += L;
            if (zpos[i] >= L)
                zpos[i] -= L;
            if (zpos[i] < 0)
                zpos[i] += L;
    }


    free(next);

    *numpairs_p = numpairs;

    return 0;
}


/*=============================================================================
---------------updates particles positions for one time-step
==============================================================================*/

int bd_step(double *xpos, double *ypos, double *zpos,int npos,double L) {
    
    //double dist; //calc dist between two points
    int i,k; //counters
    int pairs[1000*2];
    double distances2[1000];
    int numpairs;
    int error;
    
    /*----------------------------------------------------------------------
    -------------displacement by repulsive forces
    ----------------------------------------------------------------------*/
    error = interactions(npos,L,4,4.,distances2,pairs,1000,&numpairs,xpos,ypos,zpos);
    /*int interactions(int npos, const double *pos, double L, int boxdim, 
                 double cutoff2, double *distances2, 
                 int *pairs, int maxnumpairs, int *numpairs_p,
                 double *xpos, double *ypos, double *zpos)*/
    
    if (error == 1){//error 1
          printf ("\nerror: matrix dimensions are negative!\n");    
          exit(EXIT_FAILURE);      
       }   
        
    /*----------------------------------------------------------------------
    -------------Brownian displacement
    ----------------------------------------------------------------------*/    
        
        
    
    for (i =0; i<npos; i++){//for each particle
        double rsumx = 0,rsumy = 0, rsumz = 0;
        //printf("particle %i, pos %lf %lf %lf,\n",i,xpos[i],ypos[i],zpos[i]);
        
        
        for (k=0;k<12;k++){
                rsumx += rand()/(double)RAND_MAX;
                rsumy += rand()/(double)RAND_MAX;
                rsumz += rand()/(double)RAND_MAX;
        }
            
        //printf("rsumx %lf, rsumy %lf, rsumz %lf\n", rsumx,rsumy,rsumz); 
        //getchar();
        xpos[i] += pow(2*0.002,0.5)*(rsumx - 6);
        ypos[i] += pow(2*0.002,0.5)*(rsumy - 6);
        zpos[i] += pow(2*0.002,0.5)*(rsumz - 6); 
        
        //update positions if they are outside of the box
        if (xpos[i] >= L)
            xpos[i] -= L;
        if (xpos[i] < 0)
            xpos[i] += L;
        if (ypos[i] >= L)
            ypos[i] -= L;
        if (ypos[i] < 0)
            ypos[i] += L;
        if (zpos[i] >= L)
            zpos[i] -= L;
        if (zpos[i] < 0)
            zpos[i] += L;
    }   
    
    return 0; 
}





/******************************************************************************
===============================================================================
---------------main
================================================================================
*******************************************************************************/


int main (int argc, char *argv[]){

    int npos, nsteps;
    double L;
    double tm = 0;
    time_t time1;
    time_t time2;
    double timediff;

    /*--------------------------------------------------------------------------
    -------------get input parameters
    --------------------------------------------------------------------------*/
    
    //get number of particles
   int p;
   if (argc != 4){
        printf("error: too many arguments!\n");
        exit(EXIT_FAILURE);
   }else{
        for(p = 0 ; p<argc ; p++){
            printf("\nArgument %d: %s\n", p, argv[p]);
        }
        npos = atoi(argv[1]);
        L = strtod(argv[2],NULL);
        nsteps = atoi(argv[3]);
   }
    printf("number of particles is %d, length of box is %lf, number of time steps is %d\n", npos,L,nsteps);
    
    
    /*--------------------------------------------------------------------------
    -------------get initial positions
    --------------------------------------------------------------------------*/
    
    //get an array of positions
    double *xpos = (double *) malloc (npos*sizeof(double)); 
    double *ypos = (double *) malloc (npos* sizeof(double));
    double *zpos = (double *) malloc (npos* sizeof(double));
    
    int i; //counter
    for (i =0; i<npos; i++){
        xpos[i] = rand()*(L)/(double)RAND_MAX;  //get initial positions   
        ypos[i] = rand()*(L)/(double)RAND_MAX;  
        zpos[i] = rand()*(L)/(double)RAND_MAX;  
    }    
    
    
    /*--------------------------------------------------------------------------
    -------------open outfile
    --------------------------------------------------------------------------*/
    
    FILE *out;
    if ((out = fopen("traj_mod","w")) == NULL){
       printf("Error:cannot read the file output.csv\n");        
    }else{
       printf("success opening output file output.csv for output\n"); 
    }
    
    
    
    /*--------------------------------------------------------------------------
    -------------update positions and print
    --------------------------------------------------------------------------*/ 
    
    //print initial positions
    fprintf(out,"%d\n",npos);
    fprintf(out,"%lf\n",tm);
        for (i =0; i<npos; i++){
            fprintf(out,"X,%lf,%lf,%lf\n",xpos[i],ypos[i],zpos[i]); 
    } 
    time1 = time(NULL);
    int k;//counter
    for (k = 0; k <= nsteps; k++){
        
        tm += 0.002;
        
        printf("updating time step %d\n", k);
        bd_step(&xpos[0],&ypos[0],&zpos[0],npos,L); //update positions of particles
        
        fprintf(out,"%d\n",npos);
        fprintf(out,"%lf\n",tm);
        for (i =0; i<npos; i++){
            fprintf(out,"X,%lf,%lf,%lf\n",xpos[i],ypos[i],zpos[i]); //print current positions
        } 
    
    }    
    time2 = time(NULL);
    
    timediff = difftime(time2,time1); //get timediff
    printf ("timediff %f time 1 %ld time 2 %ld\n", timediff, time1, time2);
    
    free(xpos); 
    free(ypos);
    free(zpos);
    
    fclose(out);
    
    printf("end of main\n");
    return 0; 
}
