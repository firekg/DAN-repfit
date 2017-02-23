/*
	WaveAccess.h -- equates for WaveAccess XOP
*/

/* WaveAccess custom error codes */

#define OLD_IGOR 1 + FIRST_XOP_ERR
#define NON_EXISTENT_WAVE 2 + FIRST_XOP_ERR
#define NEEDS_1D_WAVE 3 + FIRST_XOP_ERR
#define NEEDS_2D_WAVE 4 + FIRST_XOP_ERR
#define NEEDS_3D_WAVE 5 + FIRST_XOP_ERR

/* Prototypes */
HOST_IMPORT void main(IORecHandle ioRecHandle);


