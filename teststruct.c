#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cl/CL.h>
#include <assert.h>

//real is of type cl_double
typedef cl_double real;
typedef cl_double2 real2;

typedef struct
{

    //  Mass
    real m;
    //  Position
    real2 x;
    //  Velocity
    real2 v;
    //  Force
    real2 F;
    //  Force_old
    real2 F_old;
    //  Bodytype
    cl_char body;

} Particle;

//  Datastructure of linked list
typedef struct ParticleList
{

    Particle p;
    struct ParticleList *next;
} ParticleList;


int main(){

    printf("%u", sizeof(ParticleList));

    return 0;
    
}