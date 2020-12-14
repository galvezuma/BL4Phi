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
uint32_t trunc3(uint32_t x){
	return x / 3 * 3;
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


/*
 * COLLAPSE LESS ABUNDANT FIRST
 */

char shrinkLetter15(char letter){
	if (letter == 'L' || letter == 'I') return 'L';
	else if (letter == 'S' || letter == 'T') return 'S';
	else if (letter == 'E' || letter == 'D') return 'E';
	else if (letter == 'K' || letter == 'R' || letter == 'H') return 'K';
	else if (letter == 'F' || letter == 'Y') return 'F';
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter16(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return 'S';
	else if (letter == 'E' || letter == 'D') return 'E';
	else if (letter == 'K' || letter == 'R' || letter == 'H') return 'K';
	else if (letter == 'F' || letter == 'Y') return 'F';
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter17(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return 'E';
	else if (letter == 'K' || letter == 'R' || letter == 'H') return 'K';
	else if (letter == 'F' || letter == 'Y') return 'F';
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter18(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'K' || letter == 'R' || letter == 'H') return 'K';
	else if (letter == 'F' || letter == 'Y') return 'F';
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter19(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'R' || letter == 'H') return 'K';
	else if (letter == 'F' || letter == 'Y') return 'F';
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'K' || letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter20(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'R' || letter == 'H') return 'K';
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'K' || letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter21(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'R' || letter == 'H') return letter;
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'K' || letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z') return letter;
	else if (letter == 0) return letter;
	return '*';
}


char (* selectedShrinkletter)(char);

void loadCluster(int c){
	if (c>21 || c<15) errorAndExit("", "Invalid cluster number");
	printf("Using cluster %d\n", c);
	char (*functions[])(char) = {shrinkLetter15, shrinkLetter16, shrinkLetter17, shrinkLetter18, shrinkLetter19, shrinkLetter20, shrinkLetter21};
	selectedShrinkletter = functions[c-15];
}

char shrinkLetter(char letter){
	return selectedShrinkletter(letter);
}

