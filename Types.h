#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <assert.h>

#pragma pack(push,1)

namespace MyProjects1 { namespace Common {

//===== INT : QWORD or DWORD ===================================================

// LINT is always 64-bit
#define LINT long long int

//!!!
//#ifdef x64
//  #define LINT long long int
//#else
//  #define LINT int
//#endif

//===== PI - associated ========================================================

#ifndef M_PI
  #define M_PI    3.14159265358979323846
#endif

#ifndef PI05
  #define PI05 (M_PI * 0.5)
#endif

#ifndef PI10
  #define PI10 M_PI
#endif

#ifndef PI20
  #define PI20 (M_PI * 2.0)
#endif

#ifndef PI40
  #define PI40 (M_PI * 4.0)
#endif

#ifndef PCI
  #define PCI (180.0 / M_PI)
#endif

#ifndef CPI
  #define CPI (M_PI / 180.0)
#endif

//===== very small value to use in templates ===================================
#ifndef NEARLY_ZERO
  #define NEARLY_ZERO ((std::numeric_limits<T>::min)( ) * T(10.0))
#endif

//===== very large =============================================================
#ifndef VERY_LARGE
  #define VERY_LARGE ((std::numeric_limits<T>::max)( ) / T(10.0))
#endif

//===== tolerance ==============================================================
#define TOLERANCE std::numeric_limits<T>::epsilon() * static_cast<T>(10.0)

//===== 0.0 ====================================================================
#define ZERO T(0.0)

//===== 1.0 ====================================================================
#define ONE T(1.0)

//===== 0.5 ====================================================================
#define HALF T(0.5)

//===== constants ==============================================================
#define T0 T(0.0)
#define T1 T(1.0)
#define T2 T(2.0)
#define T3 T(3.0)
#define T6 T(6.0)
#define T9 T(9.0)
#define T12 T(12.0)

#define T05 T(0.5)
#define T025 T(0.25)
#define T13 T(0.33333333333333)
#define T23 T(0.66666666666667)

//===== useful macros ==========================================================                                

#define NOT_DEFINED -1

#define DELETE_CLASS(class_inst) if (class_inst != NULL) delete class_inst
#define DELETE_CLASS_SETNULL(class_inst) if (class_inst != NULL) delete class_inst; class_inst = NULL
#define FREE(buffer) if (buffer != NULL) free(buffer)
#define FREE_SETNULL(buffer) if (buffer != NULL) { free(buffer); buffer = NULL; }

#define ROUND(x)  (int) (floor(x + .5))
#define SQR(x)    (pow(x,2))

#define RGBR(c) (c & 0x000000FF)
#define RGBG(c) ((c & 0x0000FF00) >> 8)
#define RGBB(c) ((c & 0x00FF0000) >> 16)

#define LIMIT(x,xmin,xmax) if (x < xmin) x = xmin; if (x > xmax) x = xmax
#define LIMIT_MIN(x,xmin) if (x < xmin) x = xmin
#define LIMIT_MAX(x,xmax) if (x > xmax) x = xmax

#define MAX(x1,x2) ((x1 > x2) ? x1 : x2)
#define MIN(x1,x2) ((x1 < x2) ? x1 : x2)

#define MAX3(a,b,c) ((a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c))
#define MIN3(a,b,c) ((a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c))

#define BETWEEN_EXC(x,x1,x2) ((x > x1) && (x < x2))
#define BETWEEN_INC(x,x1,x2) ((x >= x1) && (x <= x2))

#define SWAP(T,x1,x2) { T temp = x1; x1 = x2; x2 = temp; }

#define INSIDE_TOLERANT(u,umin,umax,tol) (u > (umin - tol) && u < (umax + tol))
#define OUTSIDERIGHT_TOLERANT(u,umax,tol) (u > (umax - tol))
#define OUTSIDELEFT_TOLERANT(u,umin,tol) (u < (umin + tol))

//===== Endianness =============================================================

                              // swap 4 bytes to change endianness
void SwapBytes4(void *data);
                              // swap 2 bytes to change endianness
void SwapBytes2(void *data);

}}

#pragma pack(pop)
