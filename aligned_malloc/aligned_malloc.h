#ifndef ALIGNED_MALLOC_H
#define ALIGNED_MALLOC_H

#include <stdlib.h>

void* aligned_malloc(size_t size, size_t alignment);
void aligned_free(void* ptr);

#endif // MM_MALLOC_H

