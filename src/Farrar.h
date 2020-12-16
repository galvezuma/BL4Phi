/*
 * Farrar.h
 *
 *  Created on: 12 de sept. de 2020
 *      Author: galvez
 */

#ifndef FARRAR_H_
#define FARRAR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

typedef struct {
	int32_t * columnaActual_Max;
	//int32_t * columna_Up;
	int32_t * columna_Left;
	int32_t * columnaPrevia_Max;
} FarrarObject;

int8_t * load_matrix(char* filename);
void transformLettersIntoPositions();
void calculateFarrarProfile();
void prepareForNextQuery();
void prepareFarrarObject(FarrarObject * o);
void freeFarrarObject(FarrarObject * o);
int16_t smith_waterman_farrar(FarrarObject* o, char *sec_database, int16_t long_sec_database);

#endif /* FARRAR_H_ */
