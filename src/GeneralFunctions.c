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
 * COLLAPSE MORE ABUNDANT FIRST
 */
/*
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
	if (letter == 'L' || letter == 'I') return 'L';
	else if (letter == 'S' || letter == 'T') return 'S';
	else if (letter == 'E' || letter == 'D') return 'E';
	else if (letter == 'K' || letter == 'R') return 'K';
	else if (letter == 'F' || letter == 'Y') return 'F';
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z' || letter == 'H') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter17(char letter){
	if (letter == 'L' || letter == 'I') return 'L';
	else if (letter == 'S' || letter == 'T') return 'S';
	else if (letter == 'E' || letter == 'D') return 'E';
	else if (letter == 'K' || letter == 'R') return 'K';
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z' || letter == 'H') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter18(char letter){
	if (letter == 'L' || letter == 'I') return 'L';
	else if (letter == 'S' || letter == 'T') return 'S';
	else if (letter == 'E' || letter == 'D') return 'E';
	else if (letter == 'K' || letter == 'R') return letter;
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z' || letter == 'H') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter19(char letter){
	if (letter == 'L' || letter == 'I') return 'L';
	else if (letter == 'S' || letter == 'T') return 'S';
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'K' || letter == 'R') return letter;
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z' || letter == 'H') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter20(char letter){
	if (letter == 'L' || letter == 'I') return 'L';
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'K' || letter == 'R') return letter;
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z' || letter == 'H') return letter;
	else if (letter == 0) return letter;
	return '*';
}

char shrinkLetter21(char letter){
	if (letter == 'L' || letter == 'I') return letter;
	else if (letter == 'S' || letter == 'T') return letter;
	else if (letter == 'E' || letter == 'D') return letter;
	else if (letter == 'K' || letter == 'R') return letter;
	else if (letter == 'F' || letter == 'Y') return letter;
	else if (letter == 'C'  || letter == 'P' || letter == 'W') return letter;
	else if (letter == 'G' || letter == 'M' || letter == 'V' || letter == 'B' || letter == 'X') return letter;
	else if (letter == 'A' || letter == 'N' || letter == 'Q' || letter == 'Z' || letter == 'H') return letter;
	else if (letter == 0) return letter;
	return '*';
}
*/

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

/*
 * CLusters made from kmean
 */
/*

char cluster[128];

void loadCluster14();
void loadCluster15();
void loadCluster16();
void loadCluster17();
void loadCluster18();
void loadCluster19();
void loadCluster20();
void loadCluster21();

void loadCluster(int c){
	if (c>21 || c<14) errorAndExit("", "Invalid cluster number");
	printf("Using cluster %d\n", c);
	void (*functions[])() = {loadCluster14, loadCluster15, loadCluster16, loadCluster17, loadCluster18, loadCluster19, loadCluster20, loadCluster21};
	functions[c-14]();
	cluster[0] = 0;
	for(int i=1; i<128; i++) cluster[i] = cluster['*'];
	functions[c-14]();
}

char shrinkLetter(char letter){
	return cluster[letter];
}

void loadCluster21() {
	cluster['A']=	1 ;
	cluster['R']=	2;
	cluster['N']=	3;
	cluster['D']=	4 ;
	cluster['C']=	5;
	cluster['Q']=	6 ;
	cluster['E']=	7 ;
	cluster['G']=	8 ;
	cluster['H']=	9;
	cluster['I']=	10;
	cluster['L']=	11;
	cluster['K']=	12;
	cluster['M']=	13;
	cluster['F']=	14;
	cluster['P']=	15;
	cluster['S']=	16;
	cluster['T']=	17;
	cluster['W']=	18;
	cluster['Y']=	19;
	cluster['V']=	20;
	cluster['B']=	21;
	cluster['J']=	22;
	cluster['Z']=	23;
	cluster['X']=	24;
	cluster['*']=	25;
}

void loadCluster20() {
	cluster['A']=	7 ;
	cluster['R']=	10;
	cluster['N']=	20;
	cluster['D']=	3 ;
	cluster['C']=	15;
	cluster['Q']=	5 ;
	cluster['E']=	1 ;
	cluster['G']=	4 ;
	cluster['H']=	14;
	cluster['I']=	11;
	cluster['L']=	11;
	cluster['K']=	17;
	cluster['M']=	13;
	cluster['F']=	6 ;
	cluster['P']=	18;
	cluster['S']=	8 ;
	cluster['T']=	9 ;
	cluster['W']=	16;
	cluster['Y']=	19;
	cluster['V']=	11;
	cluster['B']=	3 ;
	cluster['J']=	11;
	cluster['Z']=	1 ;
	cluster['X']=	2 ;
	cluster['*']=	12;
}

void loadCluster19(){
	cluster['A']=	4 ;
	cluster['R']=	11;
	cluster['N']=	17;
	cluster['D']=	19;
	cluster['C']=	7 ;
	cluster['Q']=	18;
	cluster['E']=	9 ;
	cluster['G']=	12;
	cluster['H']=	5 ;
	cluster['I']=	2 ;
	cluster['L']=	6 ;
	cluster['K']=	11;
	cluster['M']=	3 ;
	cluster['F']=	14;
	cluster['P']=	13;
	cluster['S']=	1 ;
	cluster['T']=	16;
	cluster['W']=	15;
	cluster['Y']=	10;
	cluster['V']=	2 ;
	cluster['B']=	19;
	cluster['J']=	6 ;
	cluster['Z']=	9 ;
	cluster['X']=	4 ;
	cluster['*']=	8 ;
}

void loadCluster18(){
	cluster['A']=	2 ;
	cluster['R']=	8 ;
	cluster['N']=	17;
	cluster['D']=	4 ;
	cluster['C']=	16;
	cluster['Q']=	6 ;
	cluster['E']=	15;
	cluster['G']=	14;
	cluster['H']=	3 ;
	cluster['I']=	13;
	cluster['L']=	13;
	cluster['K']=	8 ;
	cluster['M']=	12;
	cluster['F']=	1 ;
	cluster['P']=	9 ;
	cluster['S']=	18;
	cluster['T']=	7 ;
	cluster['W']=	11;
	cluster['Y']=	5 ;
	cluster['V']=	13;
	cluster['B']=	4 ;
	cluster['J']=	13;
	cluster['Z']=	15;
	cluster['X']=	2 ;
	cluster['*']=	10;
}

void loadCluster17(){
	cluster['A']=	11;
	cluster['R']=	6 ;
	cluster['N']=	2 ;
	cluster['D']=	2 ;
	cluster['C']=	1 ;
	cluster['Q']=	3 ;
	cluster['E']=	3 ;
	cluster['G']=	16;
	cluster['H']=	5 ;
	cluster['I']=	4 ;
	cluster['L']=	4 ;
	cluster['K']=	6 ;
	cluster['M']=	7 ;
	cluster['F']=	10;
	cluster['P']=	8 ;
	cluster['S']=	13;
	cluster['T']=	15;
	cluster['W']=	9 ;
	cluster['Y']=	17;
	cluster['V']=	4 ;
	cluster['B']=	2 ;
	cluster['J']=	4 ;
	cluster['Z']=	3 ;
	cluster['X']=	12;
	cluster['*']=	14;
}

void loadCluster16(){
	cluster['A']=	6 ;
	cluster['R']=	10;
	cluster['N']=	2 ;
	cluster['D']=	4 ;
	cluster['C']=	3 ;
	cluster['Q']=	5 ;
	cluster['E']=	5 ;
	cluster['G']=	8 ;
	cluster['H']=	1 ;
	cluster['I']=	9 ;
	cluster['L']=	9 ;
	cluster['K']=	10;
	cluster['M']=	7 ;
	cluster['F']=	15;
	cluster['P']=	12;
	cluster['S']=	13;
	cluster['T']=	13;
	cluster['W']=	16;
	cluster['Y']=	11;
	cluster['V']=	9 ;
	cluster['B']=	4 ;
	cluster['J']=	9 ;
	cluster['Z']=	5 ;
	cluster['X']=	6 ;
	cluster['*']=	14;
}

void loadCluster15(){
	cluster['A']=	7 ;
	cluster['R']=	4 ;
	cluster['N']=	1 ;
	cluster['D']=	8 ;
	cluster['C']=	6 ;
	cluster['Q']=	9 ;
	cluster['E']=	9 ;
	cluster['G']=	12;
	cluster['H']=	11;
	cluster['I']=	3 ;
	cluster['L']=	3 ;
	cluster['K']=	4 ;
	cluster['M']=	10;
	cluster['F']=	15;
	cluster['P']=	2 ;
	cluster['S']=	13;
	cluster['T']=	13;
	cluster['W']=	5 ;
	cluster['Y']=	15;
	cluster['V']=	3 ;
	cluster['B']=	8 ;
	cluster['J']=	3 ;
	cluster['Z']=	9 ;
	cluster['X']=	7 ;
	cluster['*']=	14;
}

void loadCluster14(){
	cluster['A']=	3 ;
	cluster['R']=	5 ;
	cluster['N']=	1 ;
	cluster['D']=	1 ;
	cluster['C']=	2 ;
	cluster['Q']=	8 ;
	cluster['E']=	8 ;
	cluster['G']=	10;
	cluster['H']=	6 ;
	cluster['I']=	14;
	cluster['L']=	14;
	cluster['K']=	5 ;
	cluster['M']=	9 ;
	cluster['F']=	4 ;
	cluster['P']=	7 ;
	cluster['S']=	12;
	cluster['T']=	12;
	cluster['W']=	13;
	cluster['Y']=	4 ;
	cluster['V']=	14;
	cluster['B']=	1 ;
	cluster['J']=	14;
	cluster['Z']=	8 ;
	cluster['X']=	3 ;
	cluster['*']=	11;
}
*/
