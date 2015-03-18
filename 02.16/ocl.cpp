#include <iostream>
#include <CL/cl.h>
#include <cstdlib>
#define READ(...) # __VA_ARGS__
#define MEM_SIZE 1024
using namespace std;
int main(){
cl_platform_id pid = 0;
cl_device_id did = 0;
cl_context context = 0;
cl_command_queue cq = 0;
cl_mem memobj = 0;
cl_program program = 0;
cl_kernel kernel = 0;
cl_uint rnd;
cl_uint rnp;
cl_int ret;
double mem[MEM_SIZE];
char dummy[]=READ(
double roz(__global double* a){
return 1.f;
}
__kernel void main2(__global double* a,__global double* par)
{

int gid = get_global_id(0);
double nn=a[gid]-a[gid-1]+1;
barrier(CLK_LOCAL_MEM_FENCE);
a[gid] = gid+1;
});


size_t source_size=sizeof(dummy);
char *source_str=dummy;

ret = clGetPlatformIDs(1, &pid, &rnp);
ret = clGetDeviceIDs(pid, CL_DEVICE_TYPE_DEFAULT, 1, &did, &rnd);
context = clCreateContext( 0, 1, &did, 0, 0, &ret);
cq = clCreateCommandQueue(context, did, 0, &ret);


memobj = clCreateBuffer(context, CL_MEM_READ_WRITE, MEM_SIZE * sizeof(cl_double), 0, &ret);
ret = clEnqueueWriteBuffer(cq, memobj, CL_TRUE, 0, MEM_SIZE * sizeof(cl_double), mem, 0, 0, 0);


program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);

ret = clBuildProgram(program, 1, &did, 0, 0, 0);
kernel = clCreateKernel(program, "main2", &ret);



ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobj);
//tutaj koleś przyjmuje argumenty. Zmieńmy to.



size_t global_work_size[3] = {MEM_SIZE, 0, 0};
size_t local_work_size[3]  = {MEM_SIZE, 0, 0};
ret = clEnqueueNDRangeKernel(cq, kernel, 1, 0, global_work_size, local_work_size, 0, 0, 0);
ret = clEnqueueReadBuffer(cq, memobj, CL_TRUE, 0, MEM_SIZE * sizeof(cl_double), mem, 0, 0, 0);//on tu zwraca wartości
for(int i=0; i<MEM_SIZE; i++)cout<<"mem["<<i<<"] = "<<mem[i]<<endl;
clFlush(cq);
clFinish(cq);
clReleaseKernel(kernel);
clReleaseProgram(program);
clReleaseMemObject(memobj);
clReleaseCommandQueue(cq);
clReleaseContext(context);


}
