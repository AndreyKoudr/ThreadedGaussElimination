# ThreadedGaussElimination

This VS 2019 project demonstrates the use of Intel intrinsics and multithreading for speedup of Gauss elimination code for the solution of linear systems. Can be easily (?) converted into Linux as the code is only standard C++ 11/14, no Windows events etc. as in my orginal version. 

Multithreading
--------------
Simple std::thread::join(), no syncronisation, no Windows events, no WaitForMultipleObjects(). The main thread waits for all created threads to complete. Despite multiple thread recreation, the main calculation routines ThreadProc4()/ThreadProc8() take 95% of CPU time (thanks for easy profiling on VS 2019). So, it's OK.

SIMD
----
The code uses Intel intrinsics with the compiler optimiser which together produce very good result. Only old SSE2 is used.

Allocator
---------
We need a big chunk of memory in one extent to place the system matrix in it.

Class Allocator allocates a big buffer to hold the whole matrix in memory and save/restore its contents in a file. The file is created in current directory named as <I>VirtSquareMatrix.bin</I>. This directory must be writable - no checks made. You may change it to a path in a temporary directory. The file writing mechanism is made with fopen/fwrite functions - not in favor now and VS tolerate it with only <I>#define _CRT_SECURE_NO_WARNINGS</I> set in defines.h. It stores/reads file by big 10M memory chunks what maybe (I do not know) faster than to use std::fstream.

Allocation is 16-byte aligned by a normal malloc() in 64-bit code. This enables use of 16-byte aligned XMM intructions which, in my experiments give some poor improvement in speed (5-10%). Not very much, despite non-aligned instructions would require two memory reads instead of one in loading/writing XMM register contents from/to memory.

More complicated allocator (actually not used in the this console project at all) is SharedMem class derived from Allocator. It allocates buffer in shared memory which theoretically can be bigger than that created by malloc(). It works only in Windows as it creates shared memory by Windows-specific functions. To use it, initialise VirtSquareMatrix() with a SharedMem object instead of Allocator object. 

Files
-----
VirtMatrixThreads files contain thread code for floats and doubles. VirtMatrix files keeps contents of matrix and solves the system by solveSystemSimple() (regular C code with no threads and SIMD, just for speed comparison, do not use) and solveSystem() which is much faster due to use of threads and SIMD.

Test
----
The project itself is a Windows console which makes two tests, one 5x5 and one 10000x10000 (for speed).

Times for solution of system 10000 x 10000 on Aspire A515-44 with CPU AMD Ryzen 5 4500U 2.38GHz, 8Gb of RAM.

The code is compiled in Release on VS 2019, optimised.

  4-byte floats
  -------------

  solveSystemSimple() (regular C, no SIMD, no multithreading)     2655 sec

  solveSystem() (SIMD, multithreading)
  
  num threads 1                                                    172 sec
  
  num threads 2                                                    137
  
  num threads 3                                                    145
  
  num threads 4                                                    147
  
  num threads 5                                                    152
  
  num threads 6                                                    156
  
  num threads 7                                                    158
  
  num threads 8                                                    151
  
  num threads 9                                                    151
  
  num threads 10                                                   152

  Max speedup 19.4


  8-byte floats
  -------------

  solveSystemSimple() (regular C, no SIMD, no multithreading)     3624 sec

  solveSystem() (SIMD, multithreading)
  
  num threads 1                                                    330 sec
  
  num threads 2                                                    268
  
  num threads 3                                                    283
  
  num threads 4                                                    286
  
  num threads 5                                                    312
  
  num threads 6                                                    311
  
  num threads 7                                                    306
  
  num threads 8                                                    300
  
  num threads 9                                                    299
  
  num threads 10                                                   300

  Max speedup 13.5

95% of CPU time is spent on ThreadProc4()/ThreadProc8(). Simple C++ 11/14 thread construction/destruction does not take much time.

Increase in number of threads does not demonstrate sufficient speedup due to current CPU - 6 logical cores. My previuos tests on a desktop with 2 Xeons (32 cores and multiple real sets of XMM registers) made speedup of up to 140 with multiple threads.

