#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <cl/CL.h>
#include <assert.h>

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

//  A cell is defined as the anchor of a linked list
typedef ParticleList *Cell;

cl_int index(cl_int2 ic, cl_int2 nc){

    return (ic.s[0] + nc.s[0] * ic.s[1]);
}


#define CL_CHECK(_expr)                                                         \
   do {                                                                         \
     cl_int _err = _expr;                                                       \
     if (_err == CL_SUCCESS)                                                    \
       break;                                                                   \
     fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);   \
     abort();                                                                   \
   } while (0)

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

int runCL(int *a, int *b, int *result, int N){

    cl_program program[1];
    cl_kernel kernel[2];

    cl_command_queue cmd_queue;
    cl_context context;

    cl_device_id cpu = NULL, device = NULL;
    cl_uint devices_n = 0;

    cl_int err = 0;
    size_t returned_size = 0;
    size_t buffer_size;

    cl_mem a_mem, b_mem, res_mem;

    
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
		//CL_CHECK(clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, 10240, buffer, NULL));
		//printf("  EXTENSIONS = %s\n", buffer);
	}

    // choose platform
    printf("\n  Please enter prefered OpenCL platform:\n");
    scanf("%d", &pref_platform);
    
    CL_CHECK(clGetDeviceIDs(platforms[pref_platform], CL_DEVICE_TYPE_CPU, 100, &cpu, &devices_n));
    //err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 1, &cpu, NULL);
    //assert(err == CL_SUCCESS);

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

    CL_CHECK(clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(vendor_name), 
                        vendor_name, &returned_size));
    CL_CHECK(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name),
                        device_name, &returned_size));
    err = clGetDeviceInfo(device, CL_DEVICE_SVM_CAPABILITIES, sizeof(svm_caps), &svm_caps, 0);

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


    // setup context and command queue
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    assert(err == CL_SUCCESS);

    cmd_queue = clCreateCommandQueueWithProperties(context, device, 0, &err);
    assert(err == CL_SUCCESS);

    // load program source from disk
    const char *filename = "sum.cl";
    char *program_source = Read_Source_File(filename);
    program[0] = clCreateProgramWithSource(context, 1,(const char**)&program_source,
                                            NULL,&err);
    assert(err == CL_SUCCESS);

    // build the program
    err = clBuildProgram(program[0], 0, NULL, NULL, NULL, NULL);
    
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

    kernel[0] =  clCreateKernel(program[0], "sum", &err);
    assert(err == CL_SUCCESS);

    //allocate memory on the device to hold and store data
    buffer_size = sizeof(int) * N;

    ParticleList *asvm = (ParticleList *)clSVMAlloc(context, CL_MEM_READ_ONLY | CL_MEM_SVM_FINE_GRAIN_BUFFER, buffer_size, 0);
    ParticleList *bsvm = (ParticleList *)clSVMAlloc(context, CL_MEM_READ_ONLY | CL_MEM_SVM_FINE_GRAIN_BUFFER, buffer_size, 0);
    ParticleList *resultsvm = (ParticleList *)clSVMAlloc(context, CL_MEM_READ_WRITE | CL_MEM_SVM_FINE_GRAIN_BUFFER, buffer_size, 0);

    for(int i = 0; i < N; i++){

        asvm[i].p.x.s[0] = a[i];
        bsvm[i].p.x.s[0] = b[i];
        //printf("%d\n", asvm[i]);
    }

    printf("%p and %p \n", (void*)&asvm[0], (void*)&asvm[1]);


    

    //input array
    /*
    a_mem = clCreateBuffer(context, CL_MEM_READ_ONLY, buffer_size, NULL, NULL);
    CL_CHECK(clEnqueueWriteBuffer(cmd_queue, a_mem, CL_TRUE, 0, buffer_size,
                                (void*)a, 0, NULL, NULL));
    
    b_mem = clCreateBuffer(context,CL_MEM_READ_ONLY, buffer_size, NULL, NULL);
    CL_CHECK(clEnqueueWriteBuffer(cmd_queue,b_mem,CL_TRUE, 0, buffer_size,
                                (void*)b, 0, NULL, NULL));
    
    res_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);

    //wait for commands in queue to complete
    clFinish(cmd_queue);
    */
    // set kernel arguments
    CL_CHECK(clSetKernelArgSVMPointer(kernel[0], 0, asvm));
    CL_CHECK(clSetKernelArgSVMPointer(kernel[0], 1, bsvm));
    CL_CHECK(clSetKernelArgSVMPointer(kernel[0], 2, resultsvm));

    /*
    CL_CHECK(clSetKernelArg(kernel[0], 0, sizeof(cl_mem), &a_mem));
    CL_CHECK(clSetKernelArg(kernel[0], 1, sizeof(cl_mem), &b_mem));
    CL_CHECK(clSetKernelArg(kernel[0], 2, sizeof(cl_mem), &res_mem));
    */
    // execute the kernel
    size_t global_work_size = N;
    CL_CHECK(clEnqueueNDRangeKernel(cmd_queue,kernel[0],1,NULL,
                                    &global_work_size, NULL, 0, NULL, NULL));
                                
    clFinish(cmd_queue);

    // read back the results
    //CL_CHECK(clEnqueueReadBuffer(cmd_queue, resultsvm, CL_TRUE, 0, buffer_size,
    //                            result, 0, NULL, NULL));

    //clFinish(cmd_queue);
    for(int i = 1; i < N; i *= 2){

        printf("%d\n", resultsvm[i].p.x.s[0]);
    }

    // free memory
    clReleaseCommandQueue(cmd_queue);
    clReleaseContext(context);
    clSVMFree(context, asvm);
    clSVMFree(context, bsvm);
    clSVMFree(context, resultsvm);

    

    return CL_SUCCESS;
}

int main(){

    int N = 100;

    int *a = (int *)malloc(N * sizeof(int));
    int *b = (int *)malloc(N * sizeof(int));
    int *result = (int *)malloc(N * sizeof(int));
    
    for(int i = 0; i < N; i++){

        a[i] = i;
        b[i] = i;
    }

    runCL(a, b, result, N);

    

    free(a);
    free(b);
    free(result);

    return 0;

}