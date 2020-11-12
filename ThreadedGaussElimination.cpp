#include <vector>
#include <chrono>
#include "Allocator.h"
#include "VirtMatrix.h"

/** Times for solution of system 10000 x 10000 on Aspire A515-44 with 
  CPU AMD Ryzen 5 4500U 2.38GHz, 8Gb of RAM.

  The code is compiled in Release on VS 2019, optimised.

  4-byte floats
  -------------

  solveSystemSimple() (regular C, no SIMD, no multithreading)     2655 sec

  solveSystem() (SIMD, multithreading)
  num threads 1                                                    172
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
  num threads 1                                                    330
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


    95% of CPU time spent on ThreadProc4()/ThreadProc8(). Simple C++ 11/14 thread 
  construction/destruction does not take much time.

    Increase of threads does not demonstrate sufficient speedup due to current 
  CPU - 6 logical cores. My tests on a desktop with 2 Xeons (32 cores and multiple 
  real sets of XMM registers) made speedup of up to 140 with multiple threads.
*/

int main()
{
  using namespace std::chrono;
  high_resolution_clock::time_point t1,t2;

                              // report progress
  progressprint = true;

  {
    /** 
      This is a test of solution for 5x5 system with known solution vector.
      solveSystemSimple() - regular C, no SIMD, no multithreading
      solveSystem() - SIMD, multithreading with 1..10 threads
    */
    #define T float
    #define ORDER 5
    #define NUM_THREADS 10

    T tolerance = static_cast<T>(0.00000001);

    T a[ORDER * ORDER] = {
       4,  1,  2,  -3,  5,
      -3,  3, -1,   4, -2,
      -1,  2,  5,   1,  3,
       5,  4,  3,  -1,  2,
       1, -2,  3,  -4,  5 };
    T b[ORDER] = {
     -16, 20, -4, -10,  3 };
    T solution[ORDER] = {
     -15.354, 15.813, -1.770, -22.148,  -6.660 };

    Allocator alloc;

    VirtSquareMatrix<T> matrix(&alloc,ORDER);

    for (int i = 0; i < ORDER * ORDER; i++)
      *matrix[i] = a[i];

    matrix.storeMatrix();

    for (int t = 0; t <= NUM_THREADS; t++)
    {
      std::vector<T> B((ORDER / 4 + 1) * 4);
      for (int i = 0; i < ORDER; i++)
        B[i] = b[i];

      bool res = false;

      t1 = high_resolution_clock::now();
      if (t == 0)
      {
        res = matrix.solveSystemSimple(&B[0], tolerance);
      } else
      {
        res = matrix.solveSystem(&B[0], tolerance, t, PIVOTING_NONE);
      }
      t2 = high_resolution_clock::now();

      assert(res);

      // compare solutions
      for (int i = 0; i < ORDER; i++)
        assert(std::abs(B[i] - solution[i]) < 0.001);

      // residual
      matrix.restoreMatrix(false);
      std::vector<T> R((ORDER / 4 + 1) * 4);
      matrix.multiply(&B[0], &R[0], t);
      VectorsSubtract(&R[0], b, (ORDER / 4 + 1) * 4);

      // max residual
      T maxr = 0;
      for (int i = 0; i < ORDER; i++)
        if (std::abs(R[i]) > maxr)
          maxr = std::abs(R[i]);

      duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

      printf("Order %d num threads %d max residual %20.15f time %f\n", ORDER, t, maxr, time_span.count());
    };

    printf("\n");

    #undef NUM_THREADS
    #undef ORDER
    #undef T
  }

  {
    /** 
      This is a speed test of solution for 10000x10000 system.
      solveSystemSimple() - regular C, no SIMD, no multithreading
      solveSystem() - SIMD, multithreading with 1..10 threads
    */

    #define T float
    #define ORDER 10000
    #define NUM_THREADS 10

    T tolerance = static_cast<T>(0.00000001);

    Allocator alloc;

    VirtSquareMatrix<T> matrix(&alloc,ORDER);

    std::vector<T> b(ORDER);
    for (int i = 0; i < ORDER; i++)
    {
      for (int j = 0; j < ORDER; j++)
      {
        *matrix.getElement(i,j) = static_cast<T>(1.0) / static_cast<T>(std::abs(i - j) + 1);
      }
      b[i] = static_cast<T>(i) / static_cast<T>(ORDER);
    }

    matrix.storeMatrix();

    for (int t = 0; t <= NUM_THREADS; t++) 
    {
      std::vector<T> B((ORDER / 4 + 1) * 4);
      for (int i = 0; i < ORDER; i++)
        B[i] = b[i];

      bool res = false;

      t1 = high_resolution_clock::now();
      if (t == 0)
      {
        res = matrix.solveSystemSimple(&B[0], tolerance);
      } else
      {
        res = matrix.solveSystem(&B[0], tolerance, t, PIVOTING_NONE);
      }
      t2 = high_resolution_clock::now();

      assert(res);

      // residual
      matrix.restoreMatrix(false);
      std::vector<T> R((ORDER / 4 + 1) * 4);
      matrix.multiply(&B[0], &R[0], t);
      VectorsSubtract(&R[0], &b[0], (ORDER / 4 + 1) * 4);

      // max residual
      T maxr = 0;
      for (int i = 0; i < ORDER; i++)
        if (std::abs(R[i]) > maxr)
          maxr = std::abs(R[i]);

      duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

      printf("Order %d num threads %d max residual %20.15f time %f\n", ORDER, t, maxr, time_span.count());
    };

    printf("\n");
    printf("Press [ENTER]\n");

    getchar();

    #undef NUM_THREADS
    #undef ORDER
    #undef T
  }

}

