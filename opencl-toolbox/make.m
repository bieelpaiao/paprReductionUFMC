% %%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opencl_include_dir = '/usr/include';
%opencl_lib_dir = '/usr/lib';
opencl_include_dir = 'D:\OCL_SDK_Light\include';
opencl_lib_dir = 'D:\OCL_SDK_Light\lib\x86_64';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mex('src/openclcmd.cpp', '-Iinclude', ['-I' opencl_include_dir], ...
    ['-L' opencl_lib_dir], '-lOpenCL' );
