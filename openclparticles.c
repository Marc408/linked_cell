#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <cl/CL.h>
#include <assert.h>

//  Dimension
#define DIM 2
//  define parameters of the Lennard-Jones-Potential
#define sigma 1
#define epsilon 5
//  define cutradius for particles
#define r_cut (2.5 * epsilon)
//  define sqr(x)
#define sqr(x) ((x) * (x))

//real is of type cl_double
typedef cl_double real;
typedef cl_double2 real2;
#if 3 == DIM
    typedef cl_double3 real3;
#endif



#define CL_CHECK(_expr)                                                         \
   do {                                                                         \
     cl_int _err = _expr;                                                       \
     if (_err == CL_SUCCESS)                                                    \
       break;                                                                   \
     fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);   \
     abort();                                                                   \
   } while (0)

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

//  A cell is defined as the anchor of a linked list
typedef ParticleList *Cell;

cl_int index(cl_int2 ic, cl_int2 nc){

    return (ic.s[0] + nc.s[0] * ic.s[1]);
}

static char* Read_Source_File(const char *filename)
{
    long int
        size = 0,
        res  = 0;

    char *src = NULL;

    FILE *file = fopen(filename, "rb");

    if (!file)  return NULL;

    if (fseek(file, 0, SEEK_END))
    {
        fclose(file);
        return NULL;
    }

    size = ftell(file);
    if (size == 0)
    {
        fclose(file);
        return NULL;
    }

    rewind(file);

    src = (char *)calloc(size + 1, sizeof(char));
    if (!src)
    {
        src = NULL;
        fclose(file);
        return src;
    }

    res = fread(src, 1, sizeof(char) * size, file);
    if (res != sizeof(char) * size)
    {
        fclose(file);
        free(src);

        return src;
    }

    src[size] = '\0'; /* NULL terminated */
    fclose(file);

    return src;
}

//  print all list elements recursive
void printList(ParticleList *pList)
{
    //  Check if pList is NULL
    if (pList == NULL)
    {

        printf("This cell is empty!\n");
    }
    else if (pList != NULL && pList->next == NULL)
    { // check if only 1 element in list

        printf("Particle of body %c at %.1f | %.1f pointing to NULL\n", pList->p.body, pList->p.x.s[0], pList->p.x.s[1]);
        printf("-> End of cell \n");
    }
    else
    {

        while (pList->next != NULL)
        {

            printf("Particle of body %c at %.1f | %.1f pointing to Particle %c at %.1f | %.1f\n", pList->p.body, pList->p.x.s[0], pList->p.x.s[1], pList->next->p.body, pList->next->p.x.s[0], pList->next->p.x.s[1]);

            pList = pList->next;
            if (pList->next == NULL)
            {

                printf("Particle of body %c at %.1f | %.1f pointing to NULL\n", pList->p.body, pList->p.x.s[0], pList->p.x.s[1]);
                printf("-> End of cell \n");
            }
        }
    }
}

// insert into linked list as first element
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


void inputParameters(real *delta_t, real *t_end, cl_int *N_A, cl_int *N_B, cl_int *nc1, cl_int *nc2, real *l1, real*l2)
{

    real delta_t_temp, t_end_temp;
    cl_int N_A_temp, N_B_temp;
    cl_int nc_temp;
    real l_temp;

    printf("  Insert delta_t:\n");
    scanf("%lf", delta_t);

    printf("  Insert t_end:\n");
    scanf("%lf", t_end);

    printf("  Insert particlecount of body A: \n");
    scanf("%d", N_A);

    printf("  Insert particlecount of body B: \n");
    scanf("%d", N_B);

    printf("  Insert length of simulationdomain in the x-dimension: \n");
    scanf("%lf", l1);
    while(*l1 < (sqrt(*N_A) + sqrt(*N_B)) || (cl_int)*l1 % 2 != 0){

        printf("  ERROR: Please enter an even length > N \n");
        printf("  Insert length of simulationdomain in the x-dimension: \n");
        scanf("%lf", l1);
    }

    printf("  Insert length of simulationdomain in the y-dimension: \n");
    scanf("%lf", l2);
    while(*l2 < (sqrt(*N_A) + sqrt(*N_B)) || (cl_int)*l2 % 2 != 0){

        printf("  ERROR: Please enter an even length > N \n");
        printf("  Insert length of simulationdomain in the y-dimension: \n");
        scanf("%lf", l2);
    }

    printf("\n\t\t\t----------Processing parameters----------\n\n");

    //set number of cells in each DIM
    *nc1 = floor((*l1 / r_cut));   
    *nc2 = floor((*l2 / r_cut));    

    printf("\n\t\t\t----------Parameters initialized---------\n\n");
}

void initData(cl_context *context, Cell *grid, cl_int N_A, cl_int N_B, cl_int2 *nc, real2 *l){

    cl_int2 kc;

    for(cl_int i = 0; i < sqrt(N_A); i++){

        for(cl_int j = 0; j < sqrt(N_A); j++){

            ParticleList *temp = (ParticleList *)clSVMAlloc(*context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER, sizeof(ParticleList), 0);
            if(temp == NULL){

                printf("  ERROR: Memory fault - allocating memory for particle in A failed\n");
            }
            temp->p.body = 'A';
            temp->p.m = 1;
            temp->p.v.s[0] = 0.1;
            temp->p.v.s[1] = 0.1;
            temp->p.x.s[0] = i;
            temp->p.x.s[1] = j;

            // init next as NULL
            temp->next = NULL;

            for(cl_int d = 0; d < DIM; d++){

                kc.s[d] = (cl_int)floor(temp->p.x.s[d] * nc->s[d] / l->s[d]);
            }

            insertList(&grid[index(kc, *nc)], temp);
        }
    }

    // build b body
    for(cl_int i = (sqrt(N_A) / 2) - (sqrt(N_B) / 2); i < (sqrt(N_A) / 2) - (sqrt(N_B) / 2) + sqrt(N_B); i++){

        for(cl_int j = sqrt(N_A) + 10; j < sqrt(N_A) + 10 + sqrt(N_B); j++){

            ParticleList *temp = (ParticleList *)clSVMAlloc(*context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER, sizeof(ParticleList), 0);
            if(temp == NULL){

                printf("  ERROR: Memory fault - allocating memory for particle in B failed\n");
            }
            temp->p.body = 'B';
            temp->p.m = 1;
            temp->p.v.s[0] = 0.1;
            temp->p.v.s[1] = -15;
            temp->p.x.s[0] = i;
            temp->p.x.s[1] = j;

            // init next as NULL
            temp->next = NULL;

            for(cl_int d = 0; d < DIM; d++){

                kc.s[d] = (cl_int)floor(temp->p.x.s[d] * nc->s[d] / l->s[d]);
            }

            insertList(&grid[index(kc, *nc)], temp);
        }
    }

}

void timeIntegration(cl_context *context, cl_command_queue *cmd_queue,
                    /*cl_kernel kernel,*/ real t, real delta_t, real t_end,
                    Cell *grid, cl_int2 *nc, real2 *l){

    while(t < t_end){

        t += delta_t;
        printf("do compx\ndo compf\ndo comv\n\n");
    }

}


Particle *prepList(ParticleList *pList, cl_context *context, size_t *array_size){

    *array_size = 0;
    ParticleList *head = pList;

    assert(pList != NULL);
    while(pList != NULL){

        pList = pList->next;
        *array_size = *array_size + 1;
    }
    
    Particle *list_array = (Particle *)clSVMAlloc(*context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER,
                                                 *array_size, 0);

    // helper
    int k = 0;

    // loop through list and copy the Particle p
    for(ParticleList *i = head; NULL != i; i = i->next){

        if(k >= *array_size){

            printf("fail\n");
        }
        list_array[k] = i->p;
        k++;
    }

    return list_array;   
}


Particle **prepGrid(ParticleList **grid, size_t pnc, size_t *sizes, cl_context *context){

        Particle **grid_array = (Particle **)clSVMAlloc(*context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER,
                                                        pnc, 0);

        printf("  vor for loop und pnc ist %d\n", pnc);
        assert(grid_array != NULL);


        for(int i = 0; i < pnc; i++){
            
            sizes[i] = 0;
            ParticleList *head = grid[i];
            printf("  grid[%d] getestet\n", i);
            ParticleList *pList = grid[i];

            while(pList != NULL){

                sizes[i] = sizes[i] + 1;
                pList = pList->next;
            }

            grid_array[i] = (Particle *)clSVMAlloc(*context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER,
                                                 sizes[i], 0);

            Particle current_array[sizes[i]];
            // helper
            int k = 0;

            // loop through list and copy the Particle p
            for(ParticleList *j = head; NULL != j; j = j->next){

                if(k >= sizes[i]){

                    printf("fail\n");
                }
                current_array[k] = j->p;

                printf(" copying %c\n", current_array[k].body);
                k++;
            }

            grid_array[i] = current_array;
            printf("%d array size\n", sizes[i]);
        }

    Particle *testarr = grid_array[0];
    printf("------%c", testarr[0].body);
    return grid_array;
}

int main(){


//---------------------------------------------------------------
//---------------------- OpenCL Setup----------------------------
//---------------------------------------------------------------
    cl_program program[1];
    cl_kernel kernel[2];

    cl_command_queue cmd_queue;
    cl_context context;

    cl_device_id cpu = NULL, device = NULL;
    cl_uint devices_n = 0;

    cl_int err = 0;
    size_t returned_size = 0;

    
	cl_platform_id platforms[100];
	cl_uint platforms_n = 0;
    int pref_platform;

	CL_CHECK(clGetPlatformIDs(100, platforms, &platforms_n));

	printf("\n=== %d OpenCL platform(s) found: ===\n", platforms_n);
	for (int i=0; i<platforms_n; i++)
	{
		char buffer[10240];
		printf("  -- %d --\n", i);
		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, 10240, buffer, NULL));
		printf("  PROFILE = %s\n", buffer);
		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, 10240, buffer, NULL));
		printf("  VERSION = %s\n", buffer);
		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 10240, buffer, NULL));
		printf("  NAME = %s\n", buffer);
		CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 10240, buffer, NULL));
		printf("  VENDOR = %s\n", buffer);
	}

    // choose platform
    printf("\n  Please enter prefered OpenCL platform:\n");
    scanf("%d", &pref_platform);
    
    CL_CHECK(clGetDeviceIDs(platforms[pref_platform], CL_DEVICE_TYPE_CPU, 100, &cpu, &devices_n));


    // Find GPU CL device
    // if there is no GPU, use the CPU
    err = clGetDeviceIDs(platforms[pref_platform], CL_DEVICE_TYPE_GPU, 100, &device, &devices_n);
    if(err != CL_SUCCESS)
        device = cpu;
    assert(device); 

    // Get information about the returned device
    cl_char vendor_name[1024] = {0};
    cl_char device_name[1024] = {0};
    cl_device_svm_capabilities svm_caps;
    cl_ulong device_max_mem;

    CL_CHECK(clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(vendor_name), 
                        vendor_name, &returned_size));
    CL_CHECK(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name),
                        device_name, &returned_size));
    err = clGetDeviceInfo(device, CL_DEVICE_SVM_CAPABILITIES, sizeof(svm_caps), &svm_caps, 0);

    clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(device_max_mem),
                        &device_max_mem, &returned_size);
    
    

    if(err == CL_SUCCESS){

        printf("  Selected device supports shared virtual memory\n");
        if(svm_caps & CL_DEVICE_SVM_FINE_GRAIN_BUFFER){

            printf("  -> Fine-grained buffer SVM available\n");
        }

        if(svm_caps & CL_DEVICE_SVM_FINE_GRAIN_SYSTEM){

            printf("  -> Fine-grained system SVM available\n");
        }

        if(svm_caps & CL_DEVICE_SVM_ATOMICS){

            printf("  -> SVM atomics supported\n");
        }
    }else{

        printf("  WARNING: Selected device does not support shared virtual memory");
    }

    printf("  Connecting to %s %s...\n", vendor_name,device_name);
    printf("  Maximum mem size is %ld\n", device_max_mem);

    // setup context and command queue
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    assert(err == CL_SUCCESS);

    cmd_queue = clCreateCommandQueueWithProperties(context, device, 0, &err);
    assert(err == CL_SUCCESS);

    // load program source from disk
    const char *filename = "force.cl";
    char *program_source = Read_Source_File(filename);
    program[0] = clCreateProgramWithSource(context, 1,(const char**)&program_source,
                                            NULL,&err);
    assert(err == CL_SUCCESS);

    // build the program
    err = clBuildProgram(program[0], 0, NULL, NULL, NULL, NULL);
    
    // get error log if clBuildProgram fails (jit compiler)
    if(err == CL_BUILD_PROGRAM_FAILURE) {
    // Determine the size of the log
    size_t log_size;
    clGetProgramBuildInfo(program[0], device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    // Allocate memory for the log
    char *log = (char *) malloc(log_size);

    // Get the log
    clGetProgramBuildInfo(program[0], device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    // Print the log
    printf("%s\n", log);
}

    kernel[0] =  clCreateKernel(program[0], "test", &err);
    assert(err == CL_SUCCESS);

//---------------------------------------------------------------
//---------------------- OpenCL Setup ENDE-----------------------
//---------------------------------------------------------------

#if 2 == DIM
    cl_int2 nc;
    real2 l;
#elif 3 == DIM
    cl_int3 nc;
    real3 l;
#endif

    cl_int N_A, N_B, N, nc1, nc2;
    size_t pnc = 1;
    real l1, l2;

    real delta_t, t_end;
    real timeused;
    size_t buffer_size;

    inputParameters(&delta_t, &t_end, &N_A, &N_B, &nc1, &nc2, &l1, &l2);

    nc.s[0] = nc1;
    nc.s[1] = nc2;

    l.s[0] = l1;
    l.s[1] = l2;

    for (int d = 0; d < DIM; d++)
            pnc *= nc.s[d];

    buffer_size = pnc * sizeof(Cell);

    printf("nc1 = %d \nnc2 = %d \nN_A = %d \nN_B = %d \npnc = %d \n\n", nc.s[0], nc.s[1], N_A, N_B, pnc);


    // allocate memory in SVM for the grid array
    Cell *grid = (Cell *)clSVMAlloc(context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER, buffer_size, 0);
    // clSVMAlloc returns a non-NULL svm adress if allocation is successful
    if (grid == NULL)
    {

        printf("  ERROR: Memory fault - Allocating memory for grid-array failed");
        return 1;
    }

    // initialize the array values
    for (int i = 0; i < pnc; i++)
    {

        grid[i] = NULL;
    }

    initData(&context, grid, N_A, N_B, &nc, &l);

    

// prep
    size_t sizes[pnc];
    Particle **argarray = prepGrid(grid, pnc, sizes, &context);

    printf("%d size\n", sizes[0]);
    for(int i = 0; i < pnc; i++){

        printf("is size %d\n", sizes[i]);
    }
    

    printf("  ended prep\n");
// prep end

    for(int i = 0; i < pnc; i++){

        Particle *ptr = argarray[i];
        printf("now arg, size is %d\n", sizes[i]);
        CL_CHECK(clSetKernelArgSVMPointer(kernel[0], 0, ptr));
        printf("did arg\n");
        CL_CHECK(clEnqueueNDRangeKernel(cmd_queue, kernel[0], 1, NULL, &sizes[i],
                                    NULL, 0, NULL, NULL));
    }
    clFinish(cmd_queue);
    printf("cmd queue finished\n");

    /*
    for(int i = 0; i < pnc; i++){
        
        ParticleList *temp = grid[i];
        Particle temp_arr[sizes[i]];

        temp_arr = argarray[i];

        for(int j = 0; j < sizes[i]; j++){
            printf("j is %d\n", j);
            temp->p = temp_arr[j];
            temp = temp->next;
    }
*/
    timeIntegration(&context, &cmd_queue, 0, delta_t, t_end, grid, &nc, &l);

    for(int i = 0; i < pnc; i++){

        printf("\n\n  ---PRINTING CELL No. %d:---\n", i);
        printList(grid[i]);
    }



    // free memory
    clSVMFree(context, argarray);    
    clReleaseCommandQueue(cmd_queue);
    clReleaseContext(context);



    printf("\n  Success!");
    return 0;
}