#include "XOPStandardHeaders.r"

resource 'vers' (1) {						/* XOP version info */
	0x01, 0x00, final, 0x00, 0,				/* version bytes and country integer */
	"1.00",
	"1.00, � 1993 WaveMetrics, Inc., all rights reserved."
};

resource 'vers' (2) {						/* Igor version info */
	0x05, 0x00, release, 0x00, 0,			/* version bytes and country integer */
	"5.00",
	"(for Igor 5.00 or later)"
};

resource 'STR#' (1100) {					/* custom error messages */
	{
		/* [1] */
		"MicroArrayFitFunc requires Igor Pro 5.0 or later.",
		/* [2] */
		"Wave does not exist.",
		/* [3] */
		"This function requires a 1D wave.",
		/* [4] */
		"This function requires a 2D wave.",
		/* [5] */
		"This function requires a 3D wave.",
	}
};

/* no menu item */

resource 'XOPI' (1100) {
	XOP_VERSION,							// XOP protocol version.
	DEV_SYS_CODE,							// Development system information.
	0,										// Obsolete - set to zero.
	0,										// Obsolete - set to zero.
	XOP_TOOLKIT_VERSION,					// XOP Toolkit version.
};

resource 'XOPF' (1100) {
	{
		"MicroArrayFitFunc",				/* function name */
		F_EXP | F_EXTERNAL,					/* function category == exponential */
		NT_FP64,							/* return value type */			
		{
			NT_FP64 + WAVE_TYPE,			/* double precision wave (coefficient wave) */
			NT_FP64 + WAVE_TYPE,			/* double precision wave (coefficient wave) */
			NT_FP64 + WAVE_TYPE,			/* double precision wave (coefficient wave) */
			NT_FP64 + WAVE_TYPE,			/* double precision wave (coefficient wave) */
			NT_FP64 + WAVE_TYPE,			/* double precision wave (coefficient wave) */
			NT_FP64 + WAVE_TYPE,			/* double precision wave (coefficient wave) */
			NT_FP64,						/* double precision x */
			NT_FP64,						/* double precision x */
			NT_FP64,						/* double precision x */
		},
		
//		"WAFill3DWaveDirectMethod",			/* function name */
//		F_WAVE | F_EXTERNAL,				/* function category */
//		NT_FP64,							/* return value type */			
//		{
//			WAVE_TYPE,						/* parameter types */
//		},

//		"WAFill3DWavePointMethod",			/* function name */
//		F_WAVE | F_EXTERNAL,				/* function category */
//		NT_FP64,							/* return value type */			
//		{
//			WAVE_TYPE,						/* parameter types */
//		},

//		"WAFill3DWaveStorageMethod",		/* function name */
//		F_WAVE | F_EXTERNAL,				/* function category */
//		NT_FP64,							/* return value type */			
//		{
//			WAVE_TYPE,						/* parameter types */
//		},

//		"WAModifyTextWave",					/* function name */
//		F_WAVE | F_EXTERNAL,				/* function category */
//		NT_FP64,							/* return value type */			
//		{
//			WAVE_TYPE,						/* parameter types */
//			HSTRING_TYPE,
//			HSTRING_TYPE,
//		},
	}
};

