clear all ; clc

ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('cl/hello_world.cl');
ocl.build();


outputBuffer = clbuffer('rw', 'uint8', 13);

helloKernel = clkernel('helloWorld', [1, 0, 0], [1, 0, 0]);

helloKernel(outputBuffer);

outputData = char(outputBuffer.get());

disp(outputData);

clear outputBuffer;

