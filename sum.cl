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


__kernel void sum(
    __global ParticleList *a,
    __global ParticleList *b,
    __global ParticleList *answer){

        int gid = get_global_id(0);
        printf("%p FOR %d and %p for %d \n", &a[gid], a[gid], &a[gid+1], a[gid]);

        answer[gid].p.x.s[0] = a[gid].p.x.s[0] + b.p.x.s[0];
    }