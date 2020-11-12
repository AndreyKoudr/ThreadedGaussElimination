#pragma once
                              // server is defined for char * only
#undef UNICODE

static_assert(_WIN32,"This class compiles under Windows only (due to shared memory specifics)");

#include <windows.h>
#include <string>

#include "Allocator.h"

/**
    Allocator for a big chunk of memory in shared memory. It is expected to be able to 
  allocate more than simple malloc() as in Allocator.
    This class is available only in Windows!
*/

class SharedMem : public Allocator {
private:

  /** Mapping */
  HANDLE mapping_ = NULL;

  /** Init mapping, buffer is nullptr in case of failure  */
  bool initMapping(const char* pname, const size_t psize, const bool clearmemory);

  /** Close mapping */
  void closeMapping();

public:

  /** Default constructor */
  SharedMem() = default;

  /** Initialise */
  virtual bool init(const char* pname, const size_t psize, const bool clearMemory);

  /** Destructor */
  ~SharedMem();
};
