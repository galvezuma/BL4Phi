/*
 * Farrar.c
 *
 *  Created on: 12 de sept. de 2016
 *      Author: galvez
 */

#include "Types.h"
#include "Farrar.h"
#include <immintrin.h>


/* Esta función utiliza en la memoria compartida:
 * un uint8_t con el número de letras del alfabeto
 * 128 uint_8 con la posición de cada letra en la matriz
 * tantos int8_t como sea (num_letras+1)*(num_letras+1)
*/
int8_t * load_matrix(char* filename) {
	int rows;
	char buffer_text[MAX_LINE_LENGTH];

	FILE * myFile = fopen(filename, "r");
	if (!myFile) {
		fprintf(stderr, "Cannot open %s.\n", filename);
		return NULL ;
	}
	for (int i = 0; i < 128; i++)
		Context.letters[i] = 255;
	rows = 0;
	while (fgets(buffer_text, MAX_LINE_LENGTH, myFile)) {
		if (buffer_text[0] == '#')
			continue; // Las líneas que empiezan por # son comentarios.
		if (buffer_text[0] == ' ') { //Es el alfabeto. Se asume que es el mismo en horizontal y en vertical
			int longitud = strlen(buffer_text);
			for (int i = 0; i < longitud; i++) {
				if (buffer_text[i] != ' ' && buffer_text[i] != '\n'
						&& buffer_text[i] != '\r') { //Suponemos que es una letra
					Context.letters[(unsigned)(buffer_text[i])] = Context.num_letters; // Se guarda la posición de la letra
					Context.num_letters++;
				}
			}
			// La última letra se asume que es * (cualquier otra cosa).
			int ultima_pos = Context.num_letters - 1;
			for (int i = 0; i < 128; i++)
				if (Context.letters[i] == 255)
					Context.letters[i] = ultima_pos;
			Context.matrix = (int8_t *)malloc(Context.num_letters*Context.num_letters*sizeof(int8_t));
			continue;
		}
		// Si se llega aquí es porque se están leyendo las filas de la matriz
		// Los datos se graban en memoria compartida conforme se van leyendo

		//Descartamos el primer elemento (una letra en columna) y el resto vamos añadiendo
		char * line = buffer_text;
		char * newLine = line + 1;
		int cols = 0;
		//Mientras siga habiendo enteros (linea!=nuevalinea), añadimos a la matriz
		short end_iteration = 0;
		while (!end_iteration) {
			line = newLine;
			int8_t value = (int8_t) strtol(line, &newLine, 0);
			if (line != newLine) {
				Context.matrix[rows*Context.num_letters+cols] = value;
				cols++;
			} else {
				end_iteration = 1;
			}
		}
		if (cols != Context.num_letters) {
			fprintf(stderr,
					"Row number %d from %s has a number of values (%d) that does not match with that of letters (%d).\n",
					rows, filename, cols, Context.num_letters);
			fclose(myFile);
			free(Context.matrix);
			return Context.matrix=NULL;
		}
		rows++;
	}
	if (rows != Context.num_letters) {
		fprintf(stderr,
				"The number of rows of %s is %d and it does not match with the number of letters (%d).\n",
				filename, rows, Context.num_letters);
		fclose(myFile);
		free(Context.matrix);
		return Context.matrix=NULL;
	}
	fclose(myFile);
	return Context.matrix;
}

void transformLettersIntoPositions() {
	Context.Farrar.letras_sec_ref = (char *)malloc(Context.sec_ref->dataLength*sizeof(char));
	// Transformation from the letters of the query sequence to positions in the matrix of letters
	for(int j=0; j<Context.sec_ref->dataLength; j++) {
		  Context.Farrar.letras_sec_ref[j] = Context.letters[(int)(Context.sec_ref->data[j])];
	}
}

void calculateFarrarProfile() {
	/* Calculus of the Farrar's Profile */
	// We fill the profile NO using Farrar (but using a longer profile multiple of 64).
	int segLen = (64 / sizeof(int32_t));
	Context.Farrar.long_profile = ((Context.sec_ref->dataLength + (63 / sizeof(int32_t))) / segLen) * segLen;
	Context.Farrar.profile = (int32_t *) _mm_malloc(Context.num_letters*Context.Farrar.long_profile*sizeof(int32_t), 64);
	printf("Tam profile %d\n", Context.num_letters*Context.Farrar.long_profile);
	int offset=0;
	for(int x=0; x<Context.num_letters;x++){
		offset = x * Context.Farrar.long_profile;
		int y;
		for(y=0; y<Context.sec_ref->dataLength; y++){
			int8_t pos_letra_ref = Context.Farrar.letras_sec_ref[y];
			Context.Farrar.profile[offset] = Context.matrix[pos_letra_ref*Context.num_letters+x];
			offset++;
		}
		// Rellenar a 0 el resto
		for( ; y < Context.Farrar.long_profile; y++){
			Context.Farrar.profile[offset] = 0;
			offset++;
		}
	}

	// Barajar el perfil segun Farrar (Fig. 1 de su artículo).
	int32_t columnaAux[Context.Farrar.long_profile] __attribute__((aligned(64)));
	for(int x=0; x<Context.num_letters; x++){
		offset = x * Context.Farrar.long_profile;
		for(int y=1; y<=Context.Farrar.long_profile; y++){
			int segLen = (64 / sizeof(int32_t));
			int numSeg = Context.Farrar.long_profile / segLen;
			int base = 1+(y-1)/segLen; // base y tt son la clave del barajar.
			int tt = (y-1) % segLen;
			columnaAux[y-1] = Context.Farrar.profile[offset + base + tt * numSeg - 1];
			//printf("orij[%d]=%d farrar[%d]=%d\n", y-1, profile[offset + y-1], base + tt * numSeg - 1, columnaAux[y-1]);
		}
		memcpy(Context.Farrar.profile + offset, columnaAux, Context.Farrar.long_profile * sizeof(int32_t));
	}
}

void prepareForNextQuery() {
	_mm_free(Context.Farrar.profile); Context.Farrar.profile = NULL;
	free(Context.Farrar.letras_sec_ref); Context.Farrar.letras_sec_ref = NULL;
}


void prepareFarrarObject(FarrarObject * o) {
	o->columnaActual_Max =  (int32_t *)_mm_malloc((Context.Farrar.long_profile)*sizeof(int32_t), 64);
	//estado_interno.columna_Up =   (int32_t *)_mm_malloc((long_profile+1)*sizeof(int32_t), 64);
	o->columna_Left = (int32_t *)_mm_malloc((Context.Farrar.long_profile)*sizeof(int32_t), 64);
	o->columnaPrevia_Max =  (int32_t *)_mm_malloc((Context.Farrar.long_profile)*sizeof(int32_t), 64);
}

void freeFarrarObject(FarrarObject * o) {
    _mm_free(o->columnaActual_Max); o->columnaActual_Max = NULL;
    _mm_free(o->columnaPrevia_Max); o->columnaPrevia_Max = NULL;
    _mm_free(o->columna_Left); o->columna_Left = NULL;
}

inline __m512i shiftRight32(__m512i a, __m512i * vZero){
	int x;
	int32_t rbuffer[64/sizeof(int32_t)] __attribute__((aligned(64)));
	// Stores into memory
	_mm512_store_epi32(rbuffer, a);
	a = _mm512_loadunpackhi_epi32(*vZero, rbuffer+15);
	return a;
}

int16_t smith_waterman_farrar(FarrarObject* o, char *sec_database, int16_t long_sec_database){
	int32_t * aux_Max;
	int16_t ret_max = 0;

	__m512i vGapOpen, vGapExtend, vZero;
	__m512i vF, vH, vMax, vE_j, vAux0;
	int segLen = (64 / sizeof(int32_t));
	int numSeg = (Context.sec_ref->dataLength + 63 / sizeof(int32_t)) / (64 / sizeof(int32_t));

		int32_t cog[segLen] __attribute__((aligned(64)));
		int32_t ceg[segLen] __attribute__((aligned(64)));
		//
		for(int x=0;x<segLen;x++) {
			cog[x] = Context.open_gap_cost;
			ceg[x] = Context.extend_gap_cost;
		}
		vGapOpen = _mm512_load_epi32(cog);
		vGapExtend = _mm512_load_epi32(ceg);
		vZero = _mm512_xor_epi32(vZero, vZero);

	vMax = _mm512_xor_epi32(vMax, vMax); // vMax = <0, 0, ..., 0>

	for(int j=0; j<Context.Farrar.long_profile; j++){
		o->columnaPrevia_Max[j] = 0;
		//o->columna_Up[j] = 0;
		o->columna_Left[j] = 0;
	}
	for(int x=0; x<long_sec_database; x++){
		// vF = <0, 0, ..., 0>
		vF = _mm512_xor_epi32(vF, vF);

		// vH = vHStore[numSeg - 1] << 1
		vH = _mm512_load_epi32(o->columnaPrevia_Max + (numSeg - 1) * segLen);
		vH = shiftRight32(vH, &vZero);

		//
		int8_t pos_letra_database = Context.letters[(int)(sec_database[x])];
		//printf("Letra %d %c %d\n", x, sec_database[x], pos_letra_database);
		int32_t offset = pos_letra_database * Context.Farrar.long_profile;
		int j;
		for(j=0; j<numSeg; j++){
			// vH = vH + vProfile[letra][j]
			int32_t * valor_match = Context.Farrar.profile + offset;
			offset += segLen;
			vAux0 = _mm512_load_epi32(valor_match);
			vH = _mm512_add_epi32(vH, vAux0);

			// vMax = max(vMax, vH);
			vMax = _mm512_max_epi32(vMax, vH);

			// vE[j] = max(vH, vE[j])
			// vH = max(vH, vF)
			vE_j = _mm512_load_epi32(o->columna_Left + j*segLen);
			vH = _mm512_max_epi32(vH, vE_j);
			vH = _mm512_max_epi32(vH, vF);

			// vHStore[j] = vH
			_mm512_store_epi32(o->columnaActual_Max + j*segLen, vH);

			// vAux = vH - vGapOpen
			vAux0 = _mm512_sub_epi32(vH, vGapOpen);
			vAux0 = _mm512_max_epi32(vAux0, vZero);
			// vE[j] = vE[j] - vGapExtend
			vE_j = _mm512_sub_epi32(vE_j, vGapExtend);
			vE_j = _mm512_max_epi32(vE_j, vZero);
			// vE[j] = max(vE[j], vAux)
			vE_j = _mm512_max_epi32(vE_j, vAux0);
			_mm512_store_epi32(o->columna_Left + j*segLen, vE_j);
			// vF = vF - vGapExtend
			vF = _mm512_sub_epi32(vF, vGapExtend);
			vF = _mm512_max_epi32(vF, vZero);
			// vF = max(vF, vAux)
			vF = _mm512_max_epi32(vF, vAux0);

			// vH = vHLoad[j]
			vH = _mm512_load_epi32(o->columnaPrevia_Max + j*segLen);
		}
		// Optimización de SWAT
		/*
		for(int x=0; x<long_profile; x++){
			printf("vMax[%d]=%d\n", x, o->columnaActual_Max[x]);
		}
		printf("Numseg: %d\n", numSeg);
		displayV("F", vF);
		*/
		// vF = vF << 1
		vF = shiftRight32(vF, &vZero);

		j = 0;
		do { // while(AnyElement(vF > vHStore[j] - vGapOpen
			vH = _mm512_load_epi32(o->columnaActual_Max + j*segLen);
			vAux0 = _mm512_sub_epi32(vH, vGapOpen);
			vAux0 = _mm512_max_epi32(vAux0, vZero);
			__mmask16 mascara = _mm512_cmpgt_epi32_mask (vF, vAux0);
			if (mascara == 0) break;
			// vHStore[j] = max(vHStore[j], vF)
			vH = _mm512_max_epi32(vH, vF);
			_mm512_store_epi32(o->columnaActual_Max + j*segLen, vH);

			// vF = vF - vGapExtend
			vF = _mm512_sub_epi32(vF, vGapExtend);
			vF = _mm512_max_epi32(vF, vZero);
			if (++j >= numSeg) {
				// vF = vF << 1
				vF = shiftRight32(vF, &vZero);
				j = 0;
			}

		} while(1);

		//
		aux_Max = o->columnaActual_Max;
		o->columnaActual_Max = o->columnaPrevia_Max;
		o->columnaPrevia_Max = aux_Max;
		//
	}

	int32_t max[segLen] __attribute__((aligned(64)));
	_mm512_store_epi32(max, vMax);
	for(int x=1;x<segLen;x++) {
		if(max[0] < max[x]) max[0] = max[x];
	}
	if (max[0] > 32767) max[0] = 32767;
	ret_max = max[0];
	return ret_max;
}

