ocl = opencl();
ocl.initialize(1,1);
ocl.addfile('cl/simple_add.cl');
ocl.build();

% Test buffer:
n = int32(10);

% Validate clkernel
addkernel = clkernel('add', [n,0,0], [n,0,0]);

% Test clobject
x = clobject(single(1:10));
y = clobject(single(11:20));
z = clobject(zeros(1,10, 'single'));

addkernel(x,y,z, uint32(10));
z1 = z.get();
x1 = x.get();
y1 = y.get();


