#include "clbessel.h"
#include <OpenCL/opencl.h>
#include <stdio.h>

#include "arrayfire.h"
#include "af/util.h"

#include "helper.h"

cl_context context = 0;
cl_program program = 0;
cl_command_queue queue = 0;

cl_kernel kbessj0 = 0;
cl_kernel kbessj1 = 0;
cl_kernel kbessy0 = 0;

static void launch_kernel(af::array &in, af::array &out, int num, cl_kernel kernel)
{
    //setup af memory
    out = af::array(in.dims());

    float *_in = in.device<float>();
    float *_out = out.device<float>();

    //sync any unfinished commands
    af::sync();

    if (!context)
        context = get_context((cl_mem)_in);

    if (!queue)
        queue = create_queue(context);

    //setup args
    cl_int err = CL_SUCCESS;
    int arg = 0;

    err |= clSetKernelArg(kernel, arg++, sizeof(cl_mem), &_in);
    err |= clSetKernelArg(kernel, arg++, sizeof(cl_mem), &_out);
    err |= clSetKernelArg(kernel, arg++, sizeof(int), &num);

    //launch kernel
    size_t local = 32;
    size_t global = local * (num / local + ((num % local) ? 1 : 0));

    clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
    clFinish(queue);

    //release arrays back to AF
    in.unlock();
    out.unlock();
}

void bessj0(af::array &in, af::array &out)
{
    launch_kernel(in, out, in.elements(), kbessj0);
}

void bessj1(af::array &in, af::array &out)
{
    launch_kernel(in, out, in.elements(), kbessj1);
}

void initKernels(af::array dummy)
{
    std::string kernel = "besselj0.cl";
    //read in kernel
    std::string kernel_str = get_kernel_string(kernel.c_str());

    float *_dum = dummy.device<float>();

    //compile kernel
    cl_context context = get_context((cl_mem)_dum);

    program = build_program(context, kernel_str);

    kbessj0 = create_kernel(program, "besslj0");
    kbessj1 = create_kernel(program, "besslj1");
    kbessy0 = create_kernel(program, "bessly0");

    dummy.unlock();
}

void bessy0(af::array &in, af::array &out)
{
    launch_kernel(in, out, in.elements(), kbessy0);
}
