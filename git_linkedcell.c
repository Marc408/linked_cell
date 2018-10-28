

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//  Dimension
#define DIM 2
//  define parameters of the Lennard-Jones-Potential
#define sigma 1
#define epsilon 5
//  define cutradius for particles
#define r_cut (2.5 * epsilon)
//  define sqr(x)
#define sqr(x) ((x) * (x))


//  macro to determine grid index by cells position
#if 1 == DIM
#define index(ic, nc) ((ic)[0])
#elif 2 == DIM
#define index(ic, nc) ((ic)[0] + (nc)[0] * (ic)[1])
#elif 3 == DIM
#define index(ic, nc) ((ic)[0] + (nc)[0] * ((ic)[1] + (nc)[1] * (ic)[2]))
#endif

//  datatype real
typedef double real;

typedef struct
{

    //  Mass
    real m;
    //  Position
    real x[DIM];
    //  Velocity
    real v[DIM];
    //  Force
    real F[DIM];
    //  Force_old
    real F_old[DIM];
    //  Bodytype
    char body;

} Particle;

//  Datastructure of linked list
typedef struct ParticleList
{

    Particle p;
    struct ParticleList *next;
} ParticleList;

//  A cell is defined as the anchor of a linked list
typedef ParticleList *Cell;

//  insert into linked list as first element
void insertList(ParticleList **root_list, ParticleList *i)
{

    i->next = *root_list;
    *root_list = i;
}

//  delete from linked list
void deleteList(ParticleList **q)
{

    *q = (*q)->next; // (*q)->next points to the deleted element
}

//  print all list elements recursive
void printList(ParticleList *pList)
{
    
    //  Check if pList is NULL
    if(pList == NULL){

        printf("ERROR: List does not exist!");
    }

    if (NULL != pList->next){

        printf("Particle of body %c at %lf %lf pointing to Particle %c\n", pList->p.body, pList->p.x[0], pList->p.x[1], pList->next->p.body);
        printList(pList->next);
    }
    else
    {

        printf("\n End of cell \n");
    }
    
}

void freeLists_LC(Cell* grid, int pnc){

    for(int i = 0; i < pnc; i++){

        free(grid[i]);
    }
}

//  Update position for a single particle
void updateX(Particle *p, real delta_t)
{

    real a = delta_t * .5 / p->m;
    for (int d = 0; d < DIM; d++)
    {

        p->x[d] += delta_t * (p->v[d] + a * p->F[d]);
        p->F_old[d] = p->F[d];
    }
}

//  Update velocity for a single particle
void updateV(Particle *p, real delta_t)
{

    real a = delta_t * .5 / p->m;
    for (int d = 0; d < DIM; d++)
    {

        p->v[d] += a * (p->F[d] + p->F_old[d]);
    }
}


//  computes the force of particle i on particle j 
void force2(Particle *i, Particle *j)
{

    real r = 0;
    for (int d = 0; d < DIM; d++)
    {

        r += sqr(j->x[d] - i->x[d]);
    }

    //  if distance of particle i and particle j is > r_cut the force is always 0
    if (r <= r_cut)
    {
        real s = sqr(sigma) / r;
        s = sqr(s) * s;
        real f = 24 * epsilon * s / r * (1 - 2 * s);
        for (int d = 0; d < DIM; d++)
        {
            i->F[d] += f * (j->x[d] - i->x[d]);
            j->F[d] -= f * (j->x[d] - i->x[d]);
        }
    }
}

// Sorts all particles into the correct list, depending on their positions
void moveParticles_LC(Cell *grid, int *nc, real *l)
{

    int ic[DIM], kc[DIM];

    //  iterate over every cell (/list)
    for (ic[0] = 0; ic[0] <= nc[0]; ic[0]++)

        for (ic[1] = 0; ic[1] <= nc[1]; ic[1]++)

#if 3 == DIM
            for (ic[2] = 0; ic[2] <= nc[2]; ic[2]++)
#endif

            {

                ParticleList **q = &grid[index(ic, nc)];
                ParticleList *i = *q;
                while (NULL != i->next)
                {
                    printf("Index %d: Particle at %f | %f pointing to %c\n", index(ic,nc), i->p.x[0], i->p.x[1], i->next->p.body);

                    //  if a particle leaves the simulation area, update its position to the opposite side
                    for(int d = 0; d < DIM; d++){

                        if(i->p.x[d] < 0)
                            i->p.x[d] = l[d] -1;
                        if(i->p.x[d] >= l[d])
                            i->p.x[d] = 0;
                    }

                    for (int d = 0; d < DIM; d++)
                    {

                        kc[d] = (int)floor(i->p.x[d] * nc[d] / l[d]);
                    }

                    if ((ic[0] != kc[0]) || (ic[1] != kc[1])
#if 3 == DIM
                        || (ic[2] != kc[2])
#endif
                    ){

                        deleteList(q);
                        insertList(&grid[index(kc, nc)], i);
                    }
                    else{

                        q = &i->next;
                    }
                    i = *q;
                }
            }
}

// iterates through every cell and updates each particles position X. Then calls moveParticles_LC() to sort into correct cell, based on the
// new position
void compX_LC(Cell *grid, int *nc, real *l, real delta_t)
{

    int ic[DIM];

    for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
        for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
#if 3 == DIM
            for (ic[2] = 0; ic[2] < nc[2]; ic[2]++)
#endif
                for (ParticleList *i = grid[index(ic, nc)]; NULL != i; i = i->next)
                    updateX(&i->p, delta_t);

    moveParticles_LC(grid, nc, l);
}

//  iterate through every cell to update each particles velocity v 
void compV_LC(Cell *grid, int *nc, real *l, real delta_t)
{

    int ic[DIM];

    for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
        for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
#if 3 == DIM
            for (ic[2] = 0; ic[2] < nc[2]; ic[2]++)
#endif
                for (ParticleList *i = grid[index(ic, nc)]; NULL != i; i = i->next)
                    updateV(&i->p, delta_t);
}


//  iterates through every cell to update each particles force f
void compF_LC(Cell *grid, int *nc)
{

    //  ic is the current cell, kc is a cell NEIGHBORING the current cell or the current cell (every cell has 8 neighbours + itself)
    int ic[DIM], kc[DIM];

    for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
        for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
#if 3 == DIM
            for (ic[2] = 0; ic[2] < nc[2]; ic[2]++)
#endif
                {
                    printf("\nIn compF_LC: Index is %d\n", index(ic,nc));
                    for (ParticleList *i = grid[index(ic, nc)]; NULL != i; i = i->next)
                    {
                        printf("Inside this for particle %c at %.1f | %.1f pointing to a Particle %c\n", i->p.body, i->p.x[0], i->p.x[1], i->next->p.body);
                        for (int d = 0; d < DIM; d++)
                            i->p.F[d] = 0;

                        for (kc[0] = ic[0] - 1; kc[0] <= ic[0] + 1; kc[0]++)
                            for (kc[1] = ic[1] - 1; kc[1] <= ic[1] + 1; kc[1]++)
    #if 3 == DIM
                                for (kc[2] = ic[2] - 1; kc[2] <= ic[2] + 1; kc[2]++)
    #endif
                                {
                                    // if the neighbor cell kc is out of range of the grid
                                    // take the cell on the opposite side of it
                                    for (int d = 0; d < DIM; d++)
                                    {
                                        if (kc[d] < 0)
                                            kc[d] = nc[d] - 1;
                                        if (kc[d] >= nc[d])
                                            kc[d] = 0;
                                    }
                                    //  TODO check distance beforehand (rcut around particle)
                                    
                                    //  for every particle j in cell kc, compute the force on particle i
                                    for (ParticleList *j = grid[index(kc, nc)]; NULL != j; j = j->next)
                                        if (i != j)
                                        {

                                            real r = 0;
                                            for (int d = 0; d < DIM; d++)
                                                r += sqr(j->p.x[d] - i->p.x[d]);
                                            if (r <= sqr(r_cut))
                                                force2(&i->p, &j->p);
                                        }
                                }
                        printf("Done with this one\n");
                    }
                }
}



//  get parameters from user input
void inputParameters_LC(real *delta_t, real *t_end, int *N_A, int *N_B, int *nc, real *l)
{

    real delta_t_temp, t_end_temp;
    int N_A_temp, N_B_temp;
    int nc_temp[DIM];
    real l_temp[DIM];

    printf("Insert delta_t:\n");
    scanf("%lf", &delta_t_temp);

    printf("Insert t_end:\n");
    scanf("%lf", &t_end_temp);

    printf("Insert Particlecount of body A: \n");
    scanf("%d", &N_A_temp);

    printf("Insert Particlecount of body B: \n");
    scanf("%d", &N_B_temp);

    printf("Insert Length of Simulationarea in the x-dimension: \n");
    scanf("%lf", &l_temp[0]);

    printf("Insert Length of Simulationarea in the y-dimension: \n");
    scanf("%lf", &l_temp[1]);

    printf("test l %lf", l_temp[1]);

    printf("\n\t\t\t----------Processing parameters----------\n\n");

    for (int d = 0; d < DIM; d++)
        nc_temp[d] = floor((l_temp[d] / r_cut));
    *delta_t = delta_t_temp;
    *t_end = t_end_temp;
    *N_A = N_A_temp;
    *N_B = N_B_temp;
    for (int d = 0; d < DIM; d++)
    {
        l[d] = l_temp[d];
        nc[d] = nc_temp[d];
    }

    printf("\n\t\t\t---------Parameters initialized--------\n\n");
}

//  generate all particles in their starting position (2 squares, B is above A, B moving towards A)
void initData_LC(Cell *grid, int N_A, int N_B, int *nc, real *l)
{

    
    int kc[DIM];

    for(int i = 0; i < sqrt(N_A); i++){
        for(int j = 0; j < sqrt(N_A); j++){

            ParticleList *temp = malloc(sizeof(ParticleList));

            temp->p.body = 'A';
            temp->p.m = 1;
            temp->p.v[0] = 0.1;
            temp->p.v[1] = 0.1;
            temp->p.x[0] = i;
            temp->p.x[1] = j;
            
            temp->next = malloc(sizeof(ParticleList));
            temp->next = NULL;

            for(int d = 0; d < DIM; d++)
                kc[d] = (int)floor(temp->p.x[d] * nc[d] / l[d]);

            insertList(&grid[index(kc, nc)], temp);
        }
    }

    for (int i = (sqrt(N_A) / 2) - (sqrt(N_B) / 2); i < (sqrt(N_A) / 2) - (sqrt(N_B) / 2) + sqrt(N_B); i++){

        for (int j = sqrt(N_A) + 10; j < (sqrt(N_A) + 10) + sqrt(N_B); j++){

            ParticleList *temp = malloc(sizeof(ParticleList));

            temp->p.body = 'B';
            temp->p.m = 1;
            temp->p.v[0] = 0.1;
            temp->p.v[1] = -15;
            temp->p.x[0] = i;
            temp->p.x[1] = j;

            temp->next = malloc(sizeof(ParticleList));
            temp->next = NULL;

            printf("\nEIN B GEMACHT\n");

            for(int d = 0; d < DIM; d++)
                kc[d] = (int)floor(temp->p.x[d] * nc[d] / l[d]);

            insertList(&grid[index(kc, nc)], temp);        }
    }

}


//  compute the simulation for each time step delta_t until t_end is reached
void timeIntegration_LC(real t, real delta_t, real t_end, Cell *grid, int *nc, real *l)
{
    // i is used for outputfile (TODO)
    int i = 0;

    printf("\nStarting the integration\n");

    compF_LC(grid, nc);

    printf("Did compf once!");

    while (t < t_end)
    {

        t += delta_t;
        compX_LC(grid, nc, l, delta_t);
        compF_LC(grid, nc);
        compV_LC(grid, nc, l, delta_t);
    

        i++;
        if (t < t_end)
            printf("\n==================================================================================\n\n");
        else
            printf("\n\t\t\t----------Integration finished----------\n\n");
    }
}


// MAIN
int main()
{

    int nc[DIM];
    int N_A, N_B, N, pnc;

    real l[DIM];
    real delta_t, t_end;

    inputParameters_LC(&delta_t, &t_end, &N_A, &N_B, nc, l);

    N = N_A + N_B;
    pnc = 1;

    for (int d = 0; d < DIM; d++)
        pnc *= nc[d];

    Cell *grid = (Cell *)malloc(pnc * sizeof(*grid));
    
    if (*grid == NULL) {
        
        printf("ERROR: Memory fault");
        return 1;
    }
    
    for(int i = 0; i < pnc; i++){

        grid[i] = malloc(sizeof(ParticleList));
        grid[i]->next = NULL;
        //grid[i] = NULL;
    }

    //DEBUG
    printf("nc1 = %d \nnc2 = %d \nN_A = %d \nN_B = %d \npnc = %d \n\n", nc[0], nc[1], N_A, N_B, pnc);

    //DEBUG END

    initData_LC(grid, N_A, N_B, nc, l);    

    //timeIntegration_LC(0, delta_t, t_end, grid, nc, l);
    
/*
    for(ParticleList *z = grid[0]; NULL != z; z = z->next){

        z->p.x[0] += 20;
        z->p.x[1] += 30;
    }

        printf("\nMOVED PARTICLES\n");
*/
    //moveParticles_LC(grid, nc, l);

    
    for(int i = 0; i < pnc; i++){

        printList(grid[i]);
    }    
    compF_LC(grid,nc);
    //printf("anchor is %c at %d %d \n", grid[0]->p.body,grid[0]->p.x[0],grid[0]->p.x[1]);

    freeLists_LC(grid, pnc);
    free(grid);
    
    //debug
    printf("\nSUCCESS\n");
    return 0;
}