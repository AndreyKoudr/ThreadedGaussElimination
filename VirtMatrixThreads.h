#pragma once

#include <assert.h>
#include "Allocator.h"

using namespace std;


/**
  Thread code for VirtSquareMatrix
*/

                                  // data for thread
typedef struct ThreadData {
  size_t numthreads;              // number of threads
  void *upperrow;                 // pointer to lead row
  void *b;                        // pointer to RHS for leading row
  size_t stepsright;              // active row length in XMM registers
  size_t shiftdown;               // offset to next row
  size_t numstepsdown;            // number of rows below lead row
  size_t rowoffset;               // offset down in rows from leading row

  void init(const size_t numThreads)
  {
    numthreads = numThreads;
    upperrow = b = NULL;
    stepsright = shiftdown = numstepsdown = 0;
    rowoffset = 1;
  }
} ThreadData;

                              // thread function
void ThreadProc4(ThreadData *data);
                              // thread function
void ThreadProc8(ThreadData *data);
                              // thread function
void ThreadProcFloat(ThreadData *data);
                              // thread function
void ThreadProcDouble(ThreadData *data);

