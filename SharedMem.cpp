#include "defines.h"
#include <assert.h>
                              // class header
#include "SharedMem.h"

bool SharedMem::init(const char* pname, const size_t psize, const bool clearMemory)
{
  initMapping(pname, psize, clearMemory);

  return (buffer != nullptr);
}

SharedMem::~SharedMem()
{
  freeStoredCopy();
  closeMapping();
}

bool SharedMem::initMapping(const char* pname, const size_t psize, const bool clearMemory)
{
  DWORD high = psize >> 32;
  DWORD low = psize & 0x00000000FFFFFFFF;
  mapping_ = CreateFileMappingA(INVALID_HANDLE_VALUE,NULL,PAGE_READWRITE,high,low,pname);
  if (mapping_ == NULL) 
    return false;

  buffer = (unsigned char *) MapViewOfFile(mapping_,FILE_MAP_ALL_ACCESS,0,0,psize);
  if (buffer == nullptr)
  {
    CloseHandle(mapping_);
    mapping_ = NULL;
    return false;
  }
                              // fix size
  size = psize;
                              // clear memory
  if (clearMemory)
    memset(buffer,0,size);

  return true;
}

void SharedMem::closeMapping()
{
  if (mapping_ != NULL) 
  {
    if (buffer != nullptr)
    {
      UnmapViewOfFile(buffer);
      buffer = nullptr;
    };
    CloseHandle(mapping_);
    mapping_ = NULL;
    size = 0;
  }
}



