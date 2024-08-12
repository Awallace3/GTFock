#include "aligned_malloc.h"
#include <stdio.h>
#include <stdlib.h>

void* aligned_malloc(size_t size, size_t alignment) {
    void* ptr = NULL;
    int result;

    result = posix_memalign(&ptr, alignment, size);
    if (result != 0) {
        // posix_memalign returns nonzero on error
        ptr = NULL; // Ensure ptr is NULL on failure, for consistency with malloc and aligned_malloc
    }
    return ptr;
}

void aligned_free(void* ptr) {
    // if (!ptr) {
    //     // Handle NULL pointer
    //     return;
    // }
    //
    // // Retrieve the original pointer which was stored just before the aligned address
    // void* original = *((void**)((size_t)ptr - sizeof(void*)));
    // free(original);
    free(ptr);
}

