#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#include "GeneralFunctions.h"

uint32_t excess64(uint32_t x){
    if ((x % 64) != 0) 
        x = (x - (x%64)) + 64;
    return x;
}
uint32_t trunc4(uint32_t x){
    return x & 0xFFFFFFFC;
}

void errorAndExit(char *data, char * text) {
    fprintf(stderr, "%s: %s\n", data, text);
    exit(1);
}


void view512iAsChar(__m512i vector){
    char data[64] __attribute__((aligned(64)));
    _mm512_store_si512 ((__m512i*) data, vector);
    for(int i=0; i<64; i++)
        printf("%c", data[i]);
    printf("\n");
}
void view512iAsBytes(__m512i vector){
    char data[16]__attribute__((aligned(64)));
    _mm512_store_si512 ((__m512i*) data, vector);
    for(int i=0; i<64; i++)
        printf("%d", data[i]);
    printf("\n");
}
