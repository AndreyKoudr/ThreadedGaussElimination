#include "defines.h"
#include <assert.h>
                              // class header
#include "Allocator.h"

using namespace std;

                              // report progress?
extern bool progressprint;

bool Allocator::init(const char* pname, const size_t psize, const bool clearMemory)
{
  assert(buffer == nullptr);
                              // on x64, buffer is 16-byte aligned by default
  buffer = (unsigned char*) malloc(psize);
                              // just in case check
  assert((reinterpret_cast<size_t>(buffer) & 0x000000000000000F) == 0);

  if (buffer != nullptr)
  {
    size = psize;

    if (clearMemory)
      memset(buffer, 0, size);

    return true;
  } else
  {
    return false;
  }
}

Allocator::Allocator(Allocator&& other)
{
  buffer = other.buffer;
  size = other.size;
  storedFileName_ = other.storedFileName_;

  other.buffer = nullptr;
  other.size = 0;
  other.storedFileName_.clear();
}

Allocator& Allocator::operator=(Allocator&& other)
{
  if (this != &other)
  {
    buffer = other.buffer;
    size = other.size;
    storedFileName_ = other.storedFileName_;

    other.buffer = nullptr;
    other.size = 0;
    other.storedFileName_.clear();
  }

  return *this;
}

Allocator::~Allocator()
{
  freeStoredCopy();

  if (buffer != nullptr)
  {
    free(buffer);
    buffer = nullptr;
    size = 0;
  }
}

/** File exists? */
static bool FileExists(const string& path)
{
  FILE* file = fopen(path.c_str(), "rb");
  if (file != NULL)
  {
    fclose(file);
    return true;
  }
  else
  {
    return false;
  }
}

void Allocator::freeStoredCopy()
{
  if (storedFileName_.length() && FileExists(storedFileName_))
  {
    if (progressprint)
      printf("Allocator stored file deleted : %s\n",storedFileName_.c_str());

    remove(storedFileName_.c_str());
  }
}

bool Allocator::storeCopy(const std::string &tempFileName)
{
  assert(buffer != NULL);
  if (buffer == NULL)
    return false;
														  // free stored matrix, only one copy can be stored				
  freeStoredCopy();
														  // allocate
  storedFileName_ = tempFileName;

  if (progressprint)
    printf("Allocator stored file : %s, size %zd\n",storedFileName_.c_str(), size);
														  // save to file
  FILE *fp = fopen(storedFileName_.c_str(),"wb");
                              // save by portions of 16M
  if (fp != NULL)
  {
    size_t portion = 16777216;
    size_t left = size;

    while (left > 0)
    {
      size_t writesize = (left > portion) ? portion : left;
      size_t written = fwrite(buffer + size - left,1,writesize,fp);
      if (written != writesize)
      {
        fclose(fp);
        return false;
      }
      left -= written;
    }

    fclose(fp);
  }

  return true;
}

/** Get file size */
static size_t fileSize(string &path)
{
  struct stat st;
  stat(path.c_str(), &st);
  size_t size = st.st_size;
  return size;
}

bool Allocator::restoreCopy(bool freestored)
{
  if (storedFileName_.length() && FileExists(storedFileName_))
  {
                            // get stored file size
    size_t filesize = fileSize(storedFileName_);

    if (filesize == size)
    {
                            // load from file
      FILE *fp = fopen(storedFileName_.c_str(),"rb");

      if (fp != NULL)
      {
                            // read by portions of 16M
        size_t portion = 16777216;
        size_t left = size;

        while (left > 0)
        {
          size_t readsize = (left > portion) ? portion : left;
          size_t read = fread(buffer + size - left,1,readsize,fp);
          if (read != readsize)
          {
            fclose(fp);
            return false;
          }   
          left -= read;
        }

        fclose(fp);
      }
													  // free stored matrix				
	    if (freestored) 
        freeStoredCopy();

      return true;
    } else
    {
      return false;
    }
  } else
  {
	  return false;
  }
}

