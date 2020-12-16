/* 
 * File:   GeneralFunctions.h
 * Author: galvez
 *
 * Created on 5 de septiembre de 2020, 19:15
 */

#ifndef GENERALFUNCTIONS_H
#define	GENERALFUNCTIONS_H

#include <inttypes.h>
#include <immintrin.h>

uint32_t excess64(uint32_t x);
uint32_t trunc4(uint32_t x);
uint32_t trunc3(uint32_t x);
void errorAndExit(char *data, char * text);

void view512iAsChar(__m512i vector);
void view512iAsBytes(__m512i vector);
void loadCluster(int c);
char shrinkLetter(char letter);

#define MIN(x,y) ((x)<(y))?(x):(y)

#endif	/* GENERALFUNCTIONS_H */

