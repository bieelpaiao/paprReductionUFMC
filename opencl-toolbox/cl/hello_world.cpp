#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS

#include <CL/cl.hpp>
#include <fstream>
#include <iostream>

int main()
{
    // Get Platform and Device
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    auto platform = platforms.front();
    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    auto device = devices.front();

    //Create context
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform(), 0 };
    cl::Context context(device, properties);

    std::ifstream helloWorldFile("hello_world.cl");
    std::string src(std::istreambuf_iterator<char>(helloWorldFile), (std::istreambuf_iterator<char>()));

    //Build Program
    cl::Program program(context, src);
    program.build("-cl-std=CL2.1");

    char buf[14];
    cl::Buffer outputBuffer = cl::Buffer(context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, sizeof(buf));

    cl::Kernel kernel(program, "HelloWorld");
    kernel.setArg(0, outputBuffer);

    cl::CommandQueue commandQueue(context, device);
    auto result = commandQueue.enqueueNDRangeKernel(kernel, 0, 1, 1);
    commandQueue.enqueueReadBuffer(outputBuffer, CL_TRUE, 0, sizeof(buf), buf);

    std::cout << buf;
    std::cin.get();

    return EXIT_SUCCESS;

}