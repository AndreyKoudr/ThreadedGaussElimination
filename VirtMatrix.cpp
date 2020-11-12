#include "defines.h"
#include "Types.h"

                              // class header
#include "VirtMatrix.h"
#include "VirtMatrixThreads.h"
                              // intrinsics
#include <xmmintrin.h>
#include <emmintrin.h>
#include <assert.h>

bool progressprint = false;

// Generates ranges array of numranges + 1 elements
template <class T>
void getRanges(T max, int numranges, std::vector<T>& ranges)
{
  T step = max / numranges;
  LIMIT_MIN(step, 1);

  ranges.resize(numranges + 1);
  T value = T(0);

  for (size_t i = 0; i <= numranges; i++)
  {
    if (value >= max)
      value = max;

    ranges[i] = value;

    value += step;
  }

  ranges[numranges] = max;
}

SquareMatrix<float>::SquareMatrix(Allocator* pmem, size_t n)
{
  mem = pmem;
                              // compute real matrix band width which should be
                              // at least 4 greater than height
  height = n;
  width = ((height / 4) + 2) * 4;
  assert((width % 4) == 0 && width >= height);
                              // allocate floats
  size_t size = width * height * sizeof(float);
                              // get matrix size in bytes
                              //!!! this will create anonymous shared memory
  OK = mem->init("",size,true);
}

SquareMatrix<double>::SquareMatrix(Allocator* pmem, size_t n)
{
  mem = pmem;
                              // compute real matrix band width which should be
                              // at least 2 greater than height
  height = n;
  width = ((height / 2) + 2) * 2;
  assert((width % 2) == 0 && width >= height);
                              // allocate floats
  size_t size = width * height * sizeof(double);
                              // get matrix size in bytes
                              //!!! this will create anonymous shared memory
  OK = mem->init("", size, true);
}

float *SquareMatrix<float>::getElement(size_t i, size_t j)
{
  assert(OK);
  assert(i >= 0 && i < height);
  assert(j >= 0 && j < height);

	return reinterpret_cast<float *>(mem->buffer) + i * width + j;
}

double *SquareMatrix<double>::getElement(size_t i, size_t j)
{
  assert(OK);
  assert(i >= 0 && i < height);
  assert(j >= 0 && j < height);

	return reinterpret_cast<double *>(mem->buffer) + i * width + j;
}

void SquareMatrix<float>::setElement(size_t i, size_t j, float value)
{
  *getElement(i,j) = value;
}

void SquareMatrix<double>::setElement(size_t i, size_t j, double value)
{
  *getElement(i,j) = value;
}

void SquareMatrix<float>::degenerateEquation(size_t index)
{
  assert(index >= 0 && index < height);
  unsigned char *e = mem->buffer + width * index * sizeof(float);
  memset(e,0,width * sizeof(float));
  setElement(index,index,float(1.0));
}

void SquareMatrix<double>::degenerateEquation(size_t index)
{
  assert(index >= 0 && index < height);
  unsigned char *e = mem->buffer + width * index * sizeof(double);
  memset(e,0,width * sizeof(double));
  setElement(index,index,double(1.0));
}

void SquareMatrix<float>::addToElement(size_t i, size_t j, float value)
{
  *getElement(i,j) += value;
}

void SquareMatrix<double>::addToElement(size_t i, size_t j, double value)
{
  *getElement(i,j) += value;
}

float *SquareMatrix<float>::operator[](size_t k)
{
  size_t i = k / height;
  size_t j = k % height;
  return getElement(i,j);
}

double *SquareMatrix<double>::operator[](size_t k)
{
  size_t i = k / height;
  size_t j = k % height;
  return getElement(i,j);
}

bool VirtSquareMatrix<float>::solveSystem(float *B, float zero, LONG numthreads, int pivoting)
{
  if (!OK)
    return false;

  assert(numthreads > 0);
                              // data to pass to threads
  std::vector<ThreadData> data(numthreads);
  data[0].init(numthreads);
                              // create threads
  if (numthreads > 1)
  {
                              // create threads
    for (int i = 0; i < numthreads; i++)
    {
                              // arguments to threads
      data[i].init(numthreads);
    }
  }
                              // for pivoting
  std::vector<float> temp(width);

															// get first RHS vector
	float *b = B;
															// pointer to diagonal element of leading row
	float *upperrow = getElement(0,0);
															// go down
	size_t height1 = height - 1;
                              // offset to the next row in floats
	size_t shiftdown = width;

	for (size_t k = 0; k < height; k++)
	{
                              // pivoting
    if (pivoting > PIVOTING_NONE)
    {
      pivotReshuffle(k,(pivoting == PIVOTING_ALL) ? (height - 1) : (k + 1),B,temp);
    }

		if (std::abs(*upperrow) < zero)
		{
      if (numthreads > 1)
      {
        for (int i = 0; i < numthreads; i++)
        {
          data[i].upperrow = nullptr;
        }
      }

			return false;
		}
															// make triangle
		size_t numstepsdown = height - k;
		size_t stepsright = numstepsdown / 4 + 1;
    numstepsdown--;
                              // divide leading row by pivot
                              // inverted diagonal element for leading row
	  float diag = float(1.0) / *upperrow;
                                // divide leading row by its pivot
    register __m128 xdiag = _mm_set1_ps(diag);
    float *x = upperrow;
    for (size_t j = 0; j < stepsright; j++)
    {
      _mm_storeu_ps(x,_mm_mul_ps(_mm_loadu_ps(x),xdiag));
      x += 4; 
    }
                              // divide RHS by pivot
    *b *= diag;
                              // here we have :
                              // upperrow     - pointer to lead row
                              // b            - pointer to RHS for leading row
                              // stepsright   - active row length in XMM registers
                              // shiftdown    - offset to next row
                              // numstepsdown - number of rows below lead row
    if (numthreads <= 1)
    {
      data[0].upperrow = upperrow;
      data[0].b = b;
      data[0].stepsright = stepsright;
      data[0].shiftdown = shiftdown;
      data[0].numstepsdown = numstepsdown;
      data[0].rowoffset = 1;

      ThreadProc4(&data[0]);
    } else
    {

                              // divide numstepsdown rows between threads
      std::vector<int> ranges;
      getRanges<int>(static_cast<int>(numstepsdown),numthreads,ranges);
                              // prepare data for each thread and signal to start
      for (int i = 0; i < numthreads; i++)
      {
        data[i].upperrow = upperrow;
        data[i].b = b;
        data[i].stepsright = stepsright;
        data[i].shiftdown = shiftdown;
        data[i].numstepsdown = ranges[i + 1] - ranges[i];
        data[i].rowoffset = ranges[i] + 1;
      }

      // threads
      std::vector<std::thread> threads(numthreads);

      for (int i = 0; i < numthreads; i++)
      {
        // create threads
        threads[i] = std::thread(&ThreadProcFloat, &data[i]);
      }

      for (auto& t : threads)
      {
        t.join();
      }
    }
															// get next diagonal element
		upperrow += (width + 1);
															// get next RHS vector
		b++;

    if (progressprint && k % 100 == 0)
      printf("%d\r",static_cast<int>(k));
	}

  if (progressprint)
    printf("\n");
                              // stop all threads
  if (numthreads > 0)
  {
                              // set data to stop all threads
    for (int i = 0; i < numthreads; i++)
    {
      data[i].upperrow = nullptr;
    }
  }
															// back substitution, not yet parallel
	b--;
	float *bj;
	float *aj;
	for (LINT k = (height - 2); k >= 0; k--)
	{
		b--;
		aj = getElement(k,k + 1);
		bj = b + 1;
		size_t numsteps = height1 - k;
		LIMIT_MAX(numsteps,height1);
		for (size_t j = 0; j < numsteps; j++) 
		{
			*b -= *bj * (*aj);
			aj++;
			bj++;
		}
	}

	return true;
};

bool VirtSquareMatrix<double>::solveSystem(double *B, double zero, LONG numthreads, int pivoting)
{
  if (!OK)
    return false;

  assert(numthreads > 0);
                              // data to pass to threads
  std::vector<ThreadData> data(numthreads);
  data[0].init(numthreads);
                              // create threads
  if (numthreads > 1)
  {
                              // create threads
    for (int i = 0; i < numthreads; i++)
    {
                              // arguments to threads
      data[i].init(numthreads);
    }
  }

                              // for pivoting
  std::vector<double> temp(width);

															// get first RHS vector
	double *b = B;
															// pointer to diagonal element of leading row
	double *upperrow = getElement(0,0);
															// go down
	size_t height1 = height - 1;
                              // offset to the next row in floats
	size_t shiftdown = width;

	for (size_t k = 0; k < height; k++)
	{
                              // pivoting
    if (pivoting > PIVOTING_NONE)
    {
      pivotReshuffle(k,(pivoting == PIVOTING_ALL) ? (height - 1) : (k + 1),B,temp);
    }

    if (std::abs(*upperrow) < zero)
    {
      if (numthreads > 1)
      {
        for (int i = 0; i < numthreads; i++)
        {
          data[i].upperrow = nullptr;
        }
      }

      return false;
    }
                              // make triangle
		size_t numstepsdown = height - k;
		size_t stepsright = numstepsdown / 2 + 1;
    numstepsdown--;
                              // divide leading row by pivot
                              // inverted diagonal element for leading row
	  double diag = double(1.0) / *upperrow;
                                // divide leading row by its pivot
    register __m128d xdiag = _mm_set1_pd(diag);
    double *x = upperrow;
    for (size_t j = 0; j < stepsright; j++)
    {
      _mm_storeu_pd(x,_mm_mul_pd(_mm_loadu_pd(x),xdiag));
      x += 2; 
    }
                              // divide RHS by pivot
    *b *= diag;
                              // here we have :
                              // upperrow     - pointer to lead row
                              // b            - pointer to RHS for leading row
                              // stepsright   - active row length in XMM registers
                              // shiftdown    - offset to next row
                              // numstepsdown - number of rows below lead row
    if (numthreads <= 1)
    {
      data[0].upperrow = upperrow;
      data[0].b = b;
      data[0].stepsright = stepsright;
      data[0].shiftdown = shiftdown;
      data[0].numstepsdown = numstepsdown;
      data[0].rowoffset = 1;

      ThreadProc8(&data[0]);
    } else
    {
                              // divide numstepsdown rows between threads
      std::vector<int> ranges;
      getRanges<int>(static_cast<int>(numstepsdown),numthreads,ranges);
                              // prepare data for each thread and signal to start
      for (int i = 0; i < numthreads; i++)
      {
        data[i].upperrow = upperrow;
        data[i].b = b;
        data[i].stepsright = stepsright;
        data[i].shiftdown = shiftdown;
        data[i].numstepsdown = ranges[i + 1] - ranges[i];
        data[i].rowoffset = ranges[i] + 1;
      }

      // threads
      std::vector<std::thread> threads(numthreads);

      for (int i = 0; i < numthreads; i++)
      {
        // create threads
        threads[i] = std::thread(&ThreadProcDouble, &data[i]);
      }

      for (auto& t : threads)
      {
        t.join();
      }
    }
															// get next diagonal element
		upperrow += (width + 1);
															// get next RHS vector
		b++;

    if (progressprint && k % 100 == 0)
      printf("%zd\r",k);
	}

  if (progressprint)
    printf("\n");
                              // stop all threads
  if (numthreads > 0)
  {
                              // set data to stop all threads
    for (int i = 0; i < numthreads; i++)
    {
      data[i].upperrow = nullptr;
    }
  }
															// back substitution, not yet parallel
	b--;
	double *bj;
	double *aj;
	for (LINT k = (height - 2); k >= 0; k--)
	{
		b--;
		aj = getElement(k,k + 1);
		bj = b + 1;
		size_t numsteps = height1 - k;
		LIMIT_MAX(numsteps,height1);
		for (size_t j = 0; j < numsteps; j++) 
		{
			*b -= *bj * (*aj);
			aj++;
			bj++;
		}
	}

	return true;
};

                              // data for thread operation over matrix and vector
typedef struct ThreadVectorData {
  size_t numthreads;          // total number of threads
  void *matrix;               // pointer to VirtSquareMatrix<float> 
  void *vector;               // vector to multiply by (float *)
  void *result;               // result vector (float *)
  LINT row;                   // starting row number
  LINT numrows;               // num rows to multiply

  void init(size_t numThreads)
  {
    numthreads = numThreads;
    matrix = vector = result = NULL;
    row = numrows = 0;
  }
} ThreadVectorData;

															// multiply numrows of matrix by vector starting
                              // from row
void ThreadMultFloat(ThreadVectorData* data)
{
  if (data->numrows == 0)
    return;
                              // pointer to matrix
  VirtSquareMatrix<float> *matrix = (VirtSquareMatrix<float> *) data->matrix;
  float *vector = (float *) data->vector;
  float *result = (float *) data->result;

	LINT j1,j2;
	LINT height1 = static_cast<LINT>(matrix->height) - 1;
                              // starting result pointer
  float *rs = result + data->row;
                              // loop along rows
	for (LINT n = 0; n < data->numrows; n++)
	{
                              // actual row number
    LINT i = data->row + n;

		j1 = 0;
		j2 = height1;
		LINT numsteps = (j2 - j1 + 1) / 4;
		if ((j2 - j1 + 1) % 4) numsteps++;

		float *as = matrix->getElement(i,j1);
		float *vs = vector + j1;
		*rs = 0;

    register __m128 sum = _mm_set1_ps(0.0);
    float *a = as;
    float *v = vs;
    for (int j = 0; j < numsteps; j++)
    {
      register __m128 m = _mm_mul_ps(_mm_loadu_ps(a),_mm_loadu_ps(v));
                              // dot product
      register __m128 m1 = m;
      register __m128 m2 = m;
      register __m128 m3 = m;
      m1 = _mm_shuffle_ps(m1,m1,1);
      m2 = _mm_shuffle_ps(m2,m2,2);
      m3 = _mm_shuffle_ps(m3,m3,3);
                              // increment sum
      sum = _mm_add_ss(_mm_add_ss(_mm_add_ss(_mm_add_ss(sum,m),m1),m2),m3);

      a += 4; 
      v += 4; 
    }
                              // store
    _mm_store_ss(rs,sum);
                              // next XMM register
		rs++;
	}

  return;
}
															// multiply numrows of matrix by vector starting
                              // from row
void ThreadMultDouble(ThreadVectorData* data)
{
  if (data->numrows == 0)
    return;
                              // pointer to matrix
  VirtSquareMatrix<double> *matrix = (VirtSquareMatrix<double> *) data->matrix;
  double *vector = (double *) data->vector;
  double *result = (double *) data->result;

	LINT j1,j2;
	LINT height1 = static_cast<LINT>(matrix->height) - 1;
                              // starting result pointer
  double *rs = result + data->row;
                              // loop along rows
	for (LINT n = 0; n < data->numrows; n++)
	{
                              // actual row number
    LINT i = data->row + n;

		j1 = 0;
		j2 = height1;
		LINT numsteps = (j2 - j1 + 1) / 2;
		if ((j2 - j1 + 1) % 2) numsteps++;

		double *as = matrix->getElement(i,j1);
		double *vs = vector + j1;
		*rs = 0;

    register __m128d sum = _mm_set1_pd(0.0);
    double *a = as;
    double *v = vs;
    for (int j = 0; j < numsteps; j++)
    {
      register __m128d m = _mm_mul_pd(_mm_loadu_pd(a),_mm_loadu_pd(v));
                              // dot product
      register __m128d m1 = m;
      m1 = _mm_shuffle_pd(m1,m1,1);
                              // increment sum
      sum = _mm_add_sd(_mm_add_sd(sum,m),m1);

      a += 2;
      v += 2;
    }
                              // store
    _mm_store_sd(rs,sum);
                              // next XMM register
		rs++;
	}

  return;
}
															// multiply by vector V with the result in R;
															// vectors V and R MUST BE 4 FLOATS LONGER than the height
															// of the matrix, with the rest filled by zeroes
bool VirtSquareMatrix<float>::multiply(float *V, float *R, LONG numthreads)
{
  if (!OK)
    return false;

  if (numthreads <= 1)
  {
    ThreadVectorData data;
    data.init(numthreads);
    data.matrix = this;
    data.vector = V;
    data.result = R;
    data.row = 0;
    data.numrows = height;

    ThreadMultFloat(&data);
  } else
  {
                              // threads
    std::vector<thread> threads(numthreads);
                              // data to pass to threads
    std::vector<ThreadVectorData> data(numthreads);
                              // divide numstepsdown rows between threads
    std::vector<int> ranges;
    getRanges<int>(static_cast<int>(height),numthreads,ranges);
                              // prepare data for each thread and signal to start
    for (int i = 0; i < numthreads; i++)
    {
      data[i].matrix = this;
      data[i].vector = V;
      data[i].result = R;
      data[i].row = ranges[i];
      data[i].numrows = ranges[i + 1] - ranges[i];
      threads[i] = std::thread(&ThreadMultFloat, &data[i]);
    }

    // start thread
    for (auto& t : threads)
    {
      t.join();
    }
  }

  return true;
}
															// multiply by vector V with the result in R;
															// vectors V and R MUST BE 2 DOUBLES LONGER than the height
															// of the matrix, with the rest filled by zeroes
bool VirtSquareMatrix<double>::multiply(double *V, double *R, LONG numthreads)
{
  if (!OK)
    return false;

  if (numthreads <= 1)
  {
    ThreadVectorData data;
    data.init(numthreads);
    data.matrix = this;
    data.vector = V;
    data.result = R;
    data.row = 0;
    data.numrows = height;

    ThreadMultDouble(&data);
  } else
  {
                              // threads
    std::vector<thread> threads(numthreads);
                              // data to pass to threads
    std::vector<ThreadVectorData> data(numthreads);
                              // divide numstepsdown rows between threads
    std::vector<int> ranges;
    getRanges<int>(static_cast<int>(height),numthreads,ranges);
                              // prepare data for each thread and signal to start
    for (int i = 0; i < numthreads; i++)
    {
      data[i].matrix = this;
      data[i].vector = V;
      data[i].result = R;
      data[i].row = ranges[i];
      data[i].numrows = ranges[i + 1] - ranges[i];
                              // start worker threads
      threads[i] = std::thread(&ThreadMultDouble, &data[i]);
    }
                              // start threads
    for (auto& t : threads)
    {
      t.join();
    }
  }

  return true;
}

template<> void VectorsSubtract(float* v0, float* v1, size_t length)
{
  assert(length % 4 == 0);

  size_t numsteps = length / 4;
  for (size_t i = 0; i < numsteps; i++)
  {
    _mm_storeu_ps(v0, _mm_sub_ps(_mm_loadu_ps(v0), _mm_loadu_ps(v1)));
    v0 += 4;
    v1 += 4;
  }
}

template<> void VectorsSubtract(double* v0, double* v1, size_t length)
{
  assert(length % 2 == 0);

  size_t numsteps = length / 2;
  for (size_t i = 0; i < numsteps; i++)
  {
    _mm_storeu_pd(v0, _mm_sub_pd(_mm_loadu_pd(v0), _mm_loadu_pd(v1)));
    v0 += 2;
    v1 += 2;
  }
}
