
//real is of type cl_double
typedef double real;
typedef double2 real2;

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
    char body;

} Particle;

//  Datastructure of linked list
typedef struct ParticleList
{

    Particle p;
    __global struct ParticleList *next;
} ParticleList;

//  A cell is defined as the anchor of a linked list
/*typedef ParticleList *Cell;

cl_int index(cl_int2 ic, cl_int2 nc){

    return (ic.x + nc.y * ic.y);
}*/

__kernel void test(__global Particle *pList){


   int gid = (int)get_global_id(0);

   pList[gid].body = 'Z';
}