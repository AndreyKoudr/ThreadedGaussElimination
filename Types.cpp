#include "defines.h"
                              // class header
#include "Types.h"
                              // memory leaks
#include "memoryleaktrace.h"

namespace MyProjects1 { namespace Common {

                              // swap 4 bytes to change endianness
void SwapBytes4(void *data) 
{
  char *pdata = (char *) data;
  char idata[4];
  idata[0] = pdata[3];
  idata[1] = pdata[2];
  idata[2] = pdata[1];
  idata[3] = pdata[0];
  memmove(data,idata,4);
};
                              // swap 2 bytes to change endianness
void SwapBytes2(void *data) 
{
  char *pdata = (char *) data;
  char idata[2];
  idata[0] = pdata[1];
  idata[1] = pdata[0];
  memmove(data,idata,2);
};
 
}}