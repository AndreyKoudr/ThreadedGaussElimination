#pragma once

#include <assert.h>
#include <string>
#include <vector>
#include <thread>

#include "Types.h"
#include "Allocator.h"

#ifdef _WIN32
  #include "SharedMem.h"
#endif

/** Printout progress */
extern bool progressprint;

/**
    Big square matrix based on shared memory. 
  Matrix is stored by rows from top to bottom. Each row contains a number
  of float fours or pairs of doubles defined by width (width is not less
  than height).
    Maximum height of matrix is defined by shared memory available, e.g. this
  is 2GB for a 2GB computer with 1GB memory free. 
  Solutions to systems up to height 10000(simple)-20000 are real.
*/

template <class T> class SquareMatrix 
{
public:
  /** Memory to hold big matrix */
  Allocator *mem = nullptr;

  /** Memory allocated? */
  bool OK = false;

  /** Matrix height */
  size_t height = 0;

  /** Width is row length in floats (maybe greater than height) */
  size_t width = 0;

  /** Default constructor */
  SquareMatrix() = default;

  /** Constructor, n - matrix height */
  SquareMatrix(Allocator* pmem, size_t n);

  /** Destructor */
  ~SquareMatrix() = default;

	/** Get element address; i,j are positions in GLOBAL
    square matrix, each being in the range 
    0..(height - 1) */
	T *getElement(size_t i, size_t j);

  /** Set element (no checks on i,j, be careful) */
	void setElement(size_t i, size_t j, T value);

  /** Add to element value (no checks on i,j, be careful) */
	void addToElement(size_t i, size_t j, T value);

  /** k = i * N + j (i,j are for square matrix) */
  T *operator[](size_t k);

  /** Change equation to all zeroes and 1.0 at diagonal */
  void degenerateEquation(size_t index);
};

  // types of pivoting
enum {
  PIVOTING_NONE,
  PIVOTING_PARTIAL,
  PIVOTING_ALL
};

template <class T> class VirtSquareMatrix : public SquareMatrix<T> 
{
public:

  /** Constructor */
	VirtSquareMatrix() : SquareMatrix<T>() { }

  /** Constructor; n is number of unknowns */
  VirtSquareMatrix(Allocator* pmem, size_t n) : SquareMatrix<T>(pmem,n) { }

	/** Destructor */
	~VirtSquareMatrix() { }

	/** Solve system by Gauss elimination;
		system matrix IS SPOILT; use StoreMatrix()/RestoreMatrix()
		to avoid the problem; returns false if close to zero 
    pivot is encountered */
	bool solveSystem(T *B, T zero, LONG numthreads, int pivoting);

  /** MATRIX MULTIPLIER
    use it to e.g. check residuals after solution multiply matrix by vector V 
    with the result in R; vectors V and R MUST BE 4 FLOATS (2 DOUBLES) 
    LONGER than the height of the matrix, with the rest filled by zeroes */
	bool multiply(T *V, T *R, LONG numthreads);

	/** Store matrix in file VirtSquareMatrix.bin, current directory must BE WRITABLE */
	bool storeMatrix();

  /** Restore matrix */
	bool restoreMatrix(bool freestored);


  /** Solve system by regular C without threads and SIMD */
  bool solveSystemSimple(T B[], T tolerance);

private:

  /** Reshuffle matrices A and B below row; set rowEnd to height - 1 to reshuffle ALL rows
    below row; set rowEnd to row + 1 to swap only row equation with one with biggest diag
    element */
  void pivotReshuffle(size_t row, size_t rowEnd, T B[], std::vector<T> &temp);
};

template <class T> void VirtSquareMatrix<T>::pivotReshuffle(size_t row, size_t rowEnd, T B[], 
  std::vector<T>& temp)
{
  // Limit
  if (rowEnd > SquareMatrix<T>::height - 1)
    rowEnd = SquareMatrix<T>::height - 1;

  for (size_t k = row; k < rowEnd; k++)
  {
    // Find equation with the largest diagonal element
    size_t index = k;
    T dmax = std::abs(*SquareMatrix<T>::getElement(k, k));
    for (size_t i = k + 1; i < SquareMatrix<T>::height; i++)
    {
      T element = std::abs(*SquareMatrix<T>::getElement(i, k));
      if (element > dmax)
      {
        index = i;
        dmax = element;
      }
    }

    // Make permutation
    if (index != k)
    {
      // Swap matrix rows
      size_t size = sizeof(T) * SquareMatrix<T>::width;
      assert(temp.size() >= SquareMatrix<T>::width);

      memmove(&temp[0], SquareMatrix<T>::getElement(k, 0), size);
      memmove(SquareMatrix<T>::getElement(k, 0), SquareMatrix<T>::getElement(index, 0), size);
      memmove(SquareMatrix<T>::getElement(index, 0), &temp[0], size);

      // Swap right-hand side values
      T btemp = B[k];
      B[k] = B[index];
      B[index] = btemp;
    }
  }
}

template <class T> bool VirtSquareMatrix<T>::solveSystemSimple(T B[], T tolerance)
{
  LINT K,K1,J,I;
  T AKK;

  LINT N = SquareMatrix<T>::height;

  for (K = 0; K < N; K++)
  {
    K1 = K + 1;
    AKK = *SquareMatrix<T>::getElement(K,K);

    if (std::abs(AKK) < tolerance)
	  {	
	    return false;
    }

    B[K] /= AKK;
    if (K == (N-1)) 
      break;

    for (J = K1; J < N; J++)
    {
      *SquareMatrix<T>::getElement(K,J) /= AKK;
      for (I = K1; I < N; I++) 
      {
        *SquareMatrix<T>::getElement(I,J) -= *SquareMatrix<T>::getElement(I,K) * 
          (*SquareMatrix<T>::getElement(K,J));
      }

      B[J] -= (*SquareMatrix<T>::getElement(J,K)  * B[K]);
    }

    if (progressprint && K % 100 == 0)
      printf("%zd\r",K);
	}

  if (progressprint)
    printf("\n");

Back:
  K1 = K;
  K -= 1;
  if (K < 0) goto Success;
  for (J = K1; J < N; J++) B[K] -= (*SquareMatrix<T>::getElement(K,J)) * B[J];
  goto Back;

Success :
  return true;
}

template <class T> bool VirtSquareMatrix<T>::storeMatrix()
{
  return SquareMatrix<T>::mem->storeCopy("VirtSquareMatrix.bin");
}

template <class T> bool VirtSquareMatrix<T>::restoreMatrix(bool freestored)
{
  return SquareMatrix<T>::mem->restoreCopy(freestored);
}

/** Fast SIMD subtraction of floats arbitrarily aligned,
  v0 = v0 - v1, size must be % 4 = 0 for floats
  and % 2 = 0 for doubles */
template <class T> void VectorsSubtract(T* v0, T* v1, size_t length);

