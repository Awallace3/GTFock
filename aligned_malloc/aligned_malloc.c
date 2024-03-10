#include "aligned_malloc.h"
#include <stdio.h>
#include <stdlib.h>

void* aligned_malloc(size_t size, size_t alignment) {
    // if (alignment & (alignment - 1)) {
    //     // The alignment must be a power of two
    //     fprintf(stderr, "Alignment must be a power of two.\n");
    //     return NULL;
    // }
    //
    // if (alignment < sizeof(void*)) {
    //     // The alignment must be at least the size of a pointer
    //     alignment = sizeof(void*);
    // }
    //
    // void* original = NULL;
    // void* aligned = NULL;
    //
    // // Allocate extra memory to ensure we can align the memory and store the offset
    // original = malloc(size + alignment + sizeof(void*));
    // if (original) {
    //     // Calculate the aligned memory address
    //     aligned = (void*)(((size_t)original + sizeof(void*) + alignment - 1) & ~(alignment - 1));
    //     // Store the original pointer just before the aligned address
    //     *((void**)((size_t)aligned - sizeof(void*))) = original;
    // }
    //
    // return aligned;
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

