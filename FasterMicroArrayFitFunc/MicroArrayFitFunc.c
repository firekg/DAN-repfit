/*	MicroArratFitFunc.c

calling sequence in Igor: MicroArrayFitFunc(ori_mat, glo_w, it_intgr, Lngt, LnGxt_dense, n_ori, n_probe_start, n_pt_probe)

*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include "MicroArrayFitFunc.h"

// Global Variables
//static int gCallSpinProcess = 1;		// Set to 1 to all user abort (cmd dot) and background processing.


/*static int
AddCStringToHandle(						// Concatenates C string to handle.
	char *theStr,
	Handle theHand)
{
	return PtrAndHand(theStr, theHand, strlen(theStr));
}*/

/*#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned
struct WAGetWaveInfoParams {
	waveHndl w;
	Handle strH;
};
typedef struct WAGetWaveInfoParams WAGetWaveInfoParams;
#include "XOPStructureAlignmentReset.h"
*/



#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned
struct MicroArrayFitFuncParams {
	double n_pt_probe;
	double n_probe_start;
	double n_ori;
	waveHndl miarray_theo_yw_hndl;
	waveHndl LnGxt_dense_hndl;
	waveHndl Lngt_hndl;
	waveHndl it_intgr_hndl;
	waveHndl glo_w_hndl;
	waveHndl ori_mat_hndl;
	double result;
};
typedef struct MicroArrayFitFuncParams MicroArrayFitFuncParams;
#include "XOPStructureAlignmentReset.h"



 
static int
MicroArrayFitFunc(MicroArrayFitFuncParams* p)
{
	waveHndl wavH1, wavH2, wavH3, wavH4, wavH5, wavH6;					// wavehandles for ori_mat, glo_w, t_intgr, Lngt and LnGxt 
	int waveType1, waveType2, waveType3, waveType4, waveType5, waveType6;
	long numDimensions1, numDimensions2, numDimensions3, numDimensions4, numDimensions5, numDimensions6;
	long dimensionSizes1[MAX_DIMENSIONS+1], dimensionSizes2[MAX_DIMENSIONS+1], dimensionSizes3[MAX_DIMENSIONS+1], dimensionSizes4[MAX_DIMENSIONS+1], dimensionSizes5[MAX_DIMENSIONS+1], dimensionSizes6[MAX_DIMENSIONS+1];
	char* dataStartPtr1; 
	char* dataStartPtr2; 
	char* dataStartPtr3; 
	char* dataStartPtr4; 
	char* dataStartPtr5; 
	char* dataStartPtr6; 
	long dataOffset1, dataOffset2, dataOffset3, dataOffset4, dataOffset5, dataOffset6;

	double *ori_mat_c0;							// pointer to the start of ori_mat, column 0
	double *ori_mat_c1, *ori_mat_c2;			// pointers to the start of columns 1 and 2
	double *glo_w;								// pointer to the start of glo_w
	double tp, bg, nor_fac, v, glo_tfac, glo_rate;
	double n_nonzero;
	double *it_intgr;							// pointer to the start of it_intgr
	double *Lngt;								// pointer to the start of Lngt
	double *LnGxt_dense;						// pointer to the start of LnGxt_dense
	double *miarray_theo_yw;
	int i, j, k, ori_x;
	double t_fac, rate;
	double local_n;
	double stop, start;
	double h, addto, array_ele;

	int hState1, hState2, hState3, hState4, hState5, hState6;
	int result;
	
	char buffer[256];

	p->result = 0;				// The Igor function result is always zero.
	
	wavH1 = p->ori_mat_hndl;
	wavH2 = p->glo_w_hndl;
	wavH3 = p->it_intgr_hndl;
	wavH4 = p->Lngt_hndl;
	wavH5 = p->LnGxt_dense_hndl;
	wavH6 = p->miarray_theo_yw_hndl;

	if (wavH1 == NIL)
		return NOWAV;
	if (wavH2 == NIL)
		return NOWAV;
	if (wavH3 == NIL)
		return NOWAV;
	if (wavH4 == NIL)
		return NOWAV;
	if (wavH5 == NIL)
		return NOWAV;
	if (wavH6 == NIL)
		return NOWAV;

	waveType6 = WaveType(wavH6);
	waveType1 = WaveType(wavH1);
	waveType2 = WaveType(wavH2);
	waveType3 = WaveType(wavH3);
	waveType4 = WaveType(wavH4);
	waveType5 = WaveType(wavH5);

	if (waveType1 & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType2 & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType3 & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType4 & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType5 & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType6 & NT_CMPLX)
		return NO_COMPLEX_WAVE;

	if (waveType1==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	if (waveType2==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	if (waveType3==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	if (waveType4==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	if (waveType5==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	if (waveType6==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;

	if (result = MDGetWaveDimensions(wavH1, &numDimensions1, dimensionSizes1))
		return result;
	if (result = MDGetWaveDimensions(wavH2, &numDimensions2, dimensionSizes2))
		return result;
	if (result = MDGetWaveDimensions(wavH3, &numDimensions3, dimensionSizes3))
		return result;
	if (result = MDGetWaveDimensions(wavH4, &numDimensions4, dimensionSizes4))
		return result;
	if (result = MDGetWaveDimensions(wavH5, &numDimensions5, dimensionSizes5))
		return result;
	if (result = MDGetWaveDimensions(wavH6, &numDimensions6, dimensionSizes6))
		return result;

	if (numDimensions1 != 2)
		return NEEDS_2D_WAVE;
	if (numDimensions2 != 1)
		return NEEDS_1D_WAVE;
	if (numDimensions3 != 2)
		return NEEDS_2D_WAVE;
	if (numDimensions4 != 2)
		return NEEDS_2D_WAVE;
	if (numDimensions5 != 2)
		return NEEDS_2D_WAVE;
	if (numDimensions6 != 1)
		return NEEDS_1D_WAVE;

	if (result = MDAccessNumericWaveData(wavH1, kMDWaveAccessMode0, &dataOffset1))
		return result;
	if (result = MDAccessNumericWaveData(wavH2, kMDWaveAccessMode0, &dataOffset2))
		return result;
	if (result = MDAccessNumericWaveData(wavH3, kMDWaveAccessMode0, &dataOffset3))
		return result;
	if (result = MDAccessNumericWaveData(wavH4, kMDWaveAccessMode0, &dataOffset4))
		return result;
	if (result = MDAccessNumericWaveData(wavH5, kMDWaveAccessMode0, &dataOffset5))
		return result;
	if (result = MDAccessNumericWaveData(wavH6, kMDWaveAccessMode0, &dataOffset6))
		return result;
			
	hState1 = MoveLockHandle(wavH1);		// So wave data can't move. Remember to call HSetState when done.
	hState2 = MoveLockHandle(wavH2);		// So wave data can't move. Remember to call HSetState when done.
	hState3 = MoveLockHandle(wavH3);		// So wave data can't move. Remember to call HSetState when done.
	hState4 = MoveLockHandle(wavH4);		// So wave data can't move. Remember to call HSetState when done.
	hState5 = MoveLockHandle(wavH5);		// So wave data can't move. Remember to call HSetState when done.
	hState6 = MoveLockHandle(wavH6);		// So wave data can't move. Remember to call HSetState when done.


	dataStartPtr1 = (char*)(*wavH1) + dataOffset1;
	dataStartPtr2 = (char*)(*wavH2) + dataOffset2;
	dataStartPtr3 = (char*)(*wavH3) + dataOffset3;
	dataStartPtr4 = (char*)(*wavH4) + dataOffset4;
	dataStartPtr5 = (char*)(*wavH5) + dataOffset5;
	dataStartPtr6 = (char*)(*wavH6) + dataOffset6;
	
	result = 0;
	
	/* FUNCTION STARTS HERE */
//	sprintf(buffer, "n_ori: %f"CR_STR, p->n_ori);
//	XOPNotice(buffer);
//	sprintf(buffer, "n_probe_start: %f"CR_STR, p->n_probe_start);
//	XOPNotice(buffer);
//	sprintf(buffer, "n_pt_probe: %f"CR_STR, p->n_pt_probe);
//	XOPNotice(buffer);
	
	ori_mat_c0=(double*)dataStartPtr1;
	ori_mat_c1=ori_mat_c0+dimensionSizes1[0];
	ori_mat_c2=ori_mat_c1+dimensionSizes1[0];

	glo_w=(double*)dataStartPtr2;
	tp=glo_w[0];
	bg=glo_w[1];
	nor_fac=glo_w[2];
	v=glo_w[3];
//	glo_tfac=glo_w[4];
//	glo_rate=glo_w[5];

	n_nonzero=v*tp;
	
	it_intgr=(double*)dataStartPtr3;
	Lngt=(double*)dataStartPtr4;
	LnGxt_dense=(double*)dataStartPtr5;
	miarray_theo_yw=(double*)dataStartPtr6;

	for(i=0; i<p->n_ori; i++) {
		
		//sigmoid model (SM)
		t_fac = ori_mat_c1[i];
		rate = ori_mat_c2[i];
		for(j=0; j<dimensionSizes3[1]; j++){
			h=j*tp/1000;
			it_intgr[j]=1/(1+pow(t_fac/h, rate));
		}
		
		//multiple-initiator model (MIM)
/*		local_n = ori_mat_c1[i];
		for(j=0; j<dimensionSizes3[1]; j++){
			h=j*tp/1000;
			it_intgr[j]= 1 - pow( (1 - 1/(1+pow(glo_tfac/h, glo_rate)) ), local_n);
		}
*/
		// xw relatively well spaced 1kb, not the most efficient code for sparse xw
		stop = n_nonzero;
		if(n_nonzero>199){
			stop = 200;
			n_nonzero = 199;
		}
//			sprintf(buffer, "nonzero: %f"CR_STR, stop);
//			XOPNotice(buffer);
//		for(j=0; j<n_nonzero; j++){
		for(j=0; j<stop; j++){
			k=floor(999-999*j/(v*tp));
			if(it_intgr[k]>0.999){
				Lngt[j] = -100000;				
			}
			else{
				Lngt[j] = log(1-it_intgr[k]);					
			}
		}
		
		ori_x = (int)(ori_mat_c0[i]-p->n_probe_start);
//			sprintf(buffer, "ori_x: %d"CR_STR, ori_x);
//			XOPNotice(buffer);

		stop = ori_x + n_nonzero;
		if(ori_x+n_nonzero > dimensionSizes5[1]){
			stop=dimensionSizes5[1];
		}
//		for(j=ori_x; j<=ori_x+n_nonzero-1; j++){
		for(j=ori_x; j<stop; j++){
//			if(j<=dimensionSizes5[1]){
				LnGxt_dense[i+j*dimensionSizes5[0]]= Lngt[(j-ori_x)];
//			}
		}
		start = ori_x-n_nonzero;
		if(ori_x-n_nonzero < 0){
			start=0;
		}
//		for(j=ori_x-n_nonzero-1; j<=ori_x; j++){
		for(j=start; j<=ori_x; j++){
//			if(j>=0) {
				LnGxt_dense[i+j*dimensionSizes5[0]]= Lngt[(ori_x-j)];
//			}
		}
	}
	
	for(j=0; j<p->n_pt_probe; j++){
		array_ele=0;
			for(i=0; i<p->n_ori; i++){
				addto=LnGxt_dense[i+j*dimensionSizes5[0]];
				array_ele=array_ele + addto;
			}
		miarray_theo_yw[j]=((1-exp(array_ele) )*100 + bg)*nor_fac;
	}
	
	

	
	/* FUNCTION ENDS HERE */
	
	HSetState((Handle)wavH1, hState1);
	HSetState((Handle)wavH2, hState2);
	HSetState((Handle)wavH3, hState3);
	HSetState((Handle)wavH4, hState4);
	HSetState((Handle)wavH5, hState5);
	HSetState((Handle)wavH6, hState6);

	WaveHandleModified(wavH1);			// Inform Igor that we have changed the wave.
	WaveHandleModified(wavH2);			// Inform Igor that we have changed the wave.
	WaveHandleModified(wavH3);			// Inform Igor that we have changed the wave.
	WaveHandleModified(wavH4);			// Inform Igor that we have changed the wave.
	WaveHandleModified(wavH5);			// Inform Igor that we have changed the wave.
	WaveHandleModified(wavH6);			// Inform Igor that we have changed the wave.

	return result;
}

			
			
			
/*	WAFill3DWaveDirectMethod()
 
 This example shows how to access the data in a multi-dimensional wave
 using the direct method.
 
 See the top of the file for instructions on how to invoke this function
 from Igor Pro 3.0 or later.
 */

/*#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned
struct WAFill3DWaveDirectMethodParams {
	waveHndl w;
	double result;
};
typedef struct WAFill3DWaveDirectMethodParams WAFill3DWaveDirectMethodParams;
#include "XOPStructureAlignmentReset.h"

*/




/*
static int
WAFill3DWaveDirectMethod(WAFill3DWaveDirectMethodParams* p)
{
	waveHndl wavH;
	int waveType;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	char* dataStartPtr;
	/*	Pointer terminology
		dp0 points to start of double data.
		dlp points to start of double data for a layer.
		dcp points to start of double data for a column.
		dp points to a particular point of double data.
	*/
/*	double *dp0, *dlp, *dcp, *dp;				// Pointers used for double data.
	float *fp0, *flp, *fcp, *fp;				// Pointers used for float data.
	long *lp0, *llp, *lcp, *lp;					// Pointers used for long data.
	short *sp0, *slp, *scp, *sp;				// Pointers used for short data.
	char *cp0, *clp, *ccp, *cp;					// Pointers used for char data.
	unsigned long *ulp0, *ullp, *ulcp, *ulp;	// Pointers used for unsigned long data.
	unsigned short *usp0, *uslp, *uscp, *usp;	// Pointers used for unsigned short data.
	unsigned char *ucp0, *uclp, *uccp, *ucp;	// Pointers used for unsigned char data.
	long dataOffset;
	long numRows, numColumns, numLayers;
	long row, column, layer;
	long pointsPerColumn, pointsPerLayer;
	int hState;
	int result;
	
	p->result = 0;				// The Igor function result is always zero.
	
	wavH = p->w;
	if (wavH == NIL)
		return NOWAV;

	waveType = WaveType(wavH);
	if (waveType & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	
	if (result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes))
		return result;
	
	if (numDimensions != 3)
		return NEEDS_3D_WAVE;

	numRows = dimensionSizes[0];
	numColumns = dimensionSizes[1];
	numLayers = dimensionSizes[2];
	pointsPerColumn = numRows;
	pointsPerLayer = pointsPerColumn*numColumns;
	
	if (result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset))
		return result;
	
	hState = MoveLockHandle(wavH);		// So wave data can't move. Remember to call HSetState when done.
	dataStartPtr = (char*)(*wavH) + dataOffset;
	
	result = 0;
	switch (waveType) {
		case NT_FP64:
			dp0 = (double*)dataStartPtr;							// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				dlp = dp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					dcp = dlp + column*pointsPerColumn;				// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						dp = dcp + row;
						*dp = row + 1100*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_FP32:
			fp0 = (float*)dataStartPtr;								// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				flp = fp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					fcp = flp + column*pointsPerColumn;				// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						fp = fcp + row;
						*fp = row + 1220*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_I32:
			lp0 = (long*)dataStartPtr;								// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				llp = lp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					lcp = llp + column*pointsPerColumn;				// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						lp = lcp + row;
						*lp = row + 1300*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_I16:
			sp0 = (short*)dataStartPtr;								// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				slp = sp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					scp = slp + column*pointsPerColumn;				// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						sp = scp + row;
						*sp = row + 1400*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_I8:
			cp0 = (char*)dataStartPtr;								// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				clp = cp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					ccp = clp + column*pointsPerColumn;				// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						cp = ccp + row;
						*cp = row + 1500*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_I32 | NT_UNSIGNED:
			ulp0 = (unsigned long*)dataStartPtr;					// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				ullp = ulp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					ulcp = ullp + column*pointsPerColumn;			// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						ulp = ulcp + row;
						*ulp = row + 1600*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_I16 | NT_UNSIGNED:
			usp0 = (unsigned short*)dataStartPtr;					// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				uslp = usp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					uscp = uslp + column*pointsPerColumn;			// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						usp = uscp + row;
						*usp = row + 1700*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;

		case NT_I8 | NT_UNSIGNED:
			ucp0 = (unsigned char*)dataStartPtr;					// Pointer to the start of all wave data.
			for(layer=0; layer<numLayers; layer++) {
				uclp = ucp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
				for(column=0; column<numColumns; column++) {
					if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
						result = -1;								// User aborted.
						break;
					}
					uccp = uclp + column*pointsPerColumn;			// Pointer to start of data for this column.
					for(row=0; row<numRows; row++) {
						ucp = uccp + row;
						*ucp = row + 1800*column + 1000000*layer;
					}
				}
				if (result != 0)
					break;											// User abort.
			}
			break;
		
		default:	// Unknown data type - possible in a future version of Igor.
			HSetState((Handle)wavH, hState);
			return NT_FNOT_AVAIL;
			break;
	}
	
	HSetState((Handle)wavH, hState);
	WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.
	
	return result;
}

/*	WAFill3DWavePointMethod()
	
	This example shows how to access the data in a multi-dimensional wave
	using a slower but very easy access method.

	See the top of the file for instructions on how to invoke this function
	from Igor Pro 3.0 or later.
	
	By using the MDSetNumericWavePointValue routine to store into the wave, instead of
	accessing the wave directly, we relieve ourselves of the need to worry about
	the data type of the wave, at the cost of running more slowly.
*/

/*#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned
struct WAFill3DWavePointMethodParams {
	waveHndl w;
	double result;
};
typedef struct WAFill3DWavePointMethodParams WAFill3DWavePointMethodParams;
#include "XOPStructureAlignmentReset.h"

static int
WAFill3DWavePointMethod(WAFill3DWavePointMethodParams*  p)
{
	waveHndl wavH;
	int waveType;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];				// Used to pass the row, column and layer to MDSetNumericWavePointValue.
	double value[2];							// Contains, real/imaginary parts but we use the real only.
	long numRows, numColumns, numLayers;
	long row, column, layer;
	int result;
	
	p->result = 0;				// The Igor function result is always zero.
	
	wavH = p->w;
	if (wavH == NIL)
		return NOWAV;

	waveType = WaveType(wavH);
	if (waveType & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	
	if (result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes))
		return result;
	
	if (numDimensions != 3)
		return NEEDS_3D_WAVE;

	numRows = dimensionSizes[0];
	numColumns = dimensionSizes[1];
	numLayers = dimensionSizes[2];
	
	MemClear(indices, sizeof(indices));			// Unused indices must be zero.
	result = 0;
	for(layer=0; layer<numLayers; layer++) {
		indices[2] = layer;
		for(column=0; column<numColumns; column++) {
			if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
				result = -1;								// User aborted.
				break;
			}
			indices[1] = column;
			for(row=0; row<numRows; row++) {
				indices[0] = row;
				value[0] = row + 1000*column + 1000000*layer;
				if (result = MDSetNumericWavePointValue(wavH, indices, value)) {
					WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.
					return result;
				}
			}
		}
		if (result != 0)
			break;
	}

	WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.
	
	return result;
}


/*	WAFill3DWaveStorageMethod()
	
	This example shows how to access the data in a multi-dimensional wave
	using the temp storage method. It is fast and easy but requires enough
	memory for a temporary double-precision copy of the wave data.

	See the top of the file for instructions on how to invoke this function
	from Igor Pro 3.0 or later.
*/

/*#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned
struct WAFill3DWaveStorageMethodParams {
	waveHndl w;
	double result;
};
typedef struct WAFill3DWaveStorageMethodParams WAFill3DWaveStorageMethodParams;
#include "XOPStructureAlignmentReset.h"

static int
WAFill3DWaveStorageMethod(WAFill3DWaveStorageMethodParams* p)
{
	waveHndl wavH;
	int waveType;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	/*	Pointer terminology
		dp0 points to start of char data.
		dlp points to start of char data for a layer.
		dcp points to start of char data for a column.
		dp points to a particular point of char data.
	*/
/*	double *dp0, *dlp, *dcp, *dp;			// Pointers used for double data.
	long numRows, numColumns, numLayers;
	long row, column, layer;
	long pointsPerColumn, pointsPerLayer;
	long numBytes;
	double* dPtr;
	int result, result2;
	
	p->result = 0;							// The Igor function result is always zero.
	
	wavH = p->w;
	if (wavH == NIL)
		return NOWAV;

	waveType = WaveType(wavH);
	if (waveType & NT_CMPLX)
		return NO_COMPLEX_WAVE;
	if (waveType==TEXT_WAVE_TYPE)
		return NUMERIC_ACCESS_ON_TEXT_WAVE;
	
	if (result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes))
		return result;
	
	if (numDimensions != 3)
		return NEEDS_3D_WAVE;

	numRows = dimensionSizes[0];
	numColumns = dimensionSizes[1];
	numLayers = dimensionSizes[2];
	pointsPerColumn = numRows;
	pointsPerLayer = pointsPerColumn*numColumns;
	
	numBytes = WavePoints(wavH) * sizeof(double);			// Bytes needed for copy
	//	This example doesn't support complex waves.
	//	if (isComplex)
	//		numBytes *= 2;
	dPtr = (double*)NewPtr(numBytes);
	if (dPtr==NIL)
		return NOMEM;
	
	if (result = MDGetDPDataFromNumericWave(wavH, dPtr)) {	// Get a copy of the wave data.
		DisposePtr((Ptr)dPtr);
		return result;
	}
	
	result = 0;
	dp0 = dPtr;												// Pointer to the start of all wave data.
	for(layer=0; layer<numLayers; layer++) {
		dlp = dp0 + layer*pointsPerLayer;					// Pointer to start of data for this layer.
		for(column=0; column<numColumns; column++) {
			if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
				result = -1;								// User aborted.
				break;
			}
			dcp = dlp + column*pointsPerColumn;				// Pointer to start of data for this column.
			for(row=0; row<numRows; row++) {
				dp = dcp + row;
				*dp = row + 1000*column + 1000000*layer;
			}
		}
		if (result != 0)
			break;
	}
	
	if (result2 = MDStoreDPDataInNumericWave(wavH, dPtr)) {	// Store copy in the wave.
		DisposePtr((Ptr)dPtr);
		return result2;
	}
	
	DisposePtr((Ptr)dPtr);
	WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.
	
	return result;
}

static int
DoAppendAndPrepend(Handle textH, Handle prependStringH, Handle appendStringH)
{
	long textHLength;
	long appendStringHLength;
	long prependStringHLength;

	textHLength = GetHandleSize(textH);
	prependStringHLength = GetHandleSize(prependStringH);
	appendStringHLength = GetHandleSize(appendStringH);
	
	SetHandleSize(textH, textHLength + prependStringHLength + appendStringHLength);
	if (MemError())
		return NOMEM;
	memmove(*textH+prependStringHLength, *textH, textHLength);								// Make room for prependString.
	memcpy(*textH, *prependStringH, prependStringHLength);									// Prepend prependString.
	memcpy(*textH+textHLength+prependStringHLength, *appendStringH, appendStringHLength);	// Append appendString.
	return 0;
}

/*	WAModifyTextWave()
	
	This example shows how to access the data in a multi-dimensional text wave.

	See the top of the file for instructions on how to invoke this function
	from Igor Pro 3.0 or later.
*/

/*#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned
struct WAModifyTextWaveParams {
		Handle appendStringH;			// String to be appended to each wave point.
		Handle prependStringH;			// String to be prepended to each wave point.
		waveHndl w;
		double result;
};
typedef struct WAModifyTextWaveParams WAModifyTextWaveParams;
#include "XOPStructureAlignmentReset.h"

static int
WAModifyTextWave(WAModifyTextWaveParams* p)
{
	waveHndl wavH;
	int waveType;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long indices[MAX_DIMENSIONS];				// Used to pass the row, column and layer to MDSetTextWavePointValue.
	long numRows, numColumns, numLayers, numChunks;
	long row, column, layer, chunk;
	Handle textH;
	int result;
	
	result = 0;

	textH = NewHandle(0L);			// Handle used to pass text wave characters to Igor.
	if (textH == NIL) {
		result = NOMEM;
		goto done;
	}
	
	if (p->prependStringH == NIL) {
		result = USING_NULL_STRVAR;		// The user called the function with an uninitialized string variable.
		goto done;
	}
	
	if (p->appendStringH == NIL) {
		result = USING_NULL_STRVAR;		// The user called the function with an uninitialized string variable.
		goto done;
	}
	
	wavH = p->w;
	if (wavH == NIL) {
		result = NOWAV;					// The user called the function with a missing wave or uninitialized wave reference variable.
		goto done;
	}

	waveType = WaveType(wavH);
	if (waveType!=TEXT_WAVE_TYPE) {
		result = TEXT_ACCESS_ON_NUMERIC_WAVE;
		goto done;
	}
	
	if (result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes))
		goto done;

	numRows = dimensionSizes[0];
	numColumns = dimensionSizes[1];
	if (numColumns==0)
		numColumns = 1;
	numLayers = dimensionSizes[2];
	if (numLayers==0)
		numLayers = 1;
	numChunks = dimensionSizes[3];
	if (numChunks==0)
		numChunks = 1;
	
	MemClear(indices, sizeof(indices));			// Clear unused indices.
	result = 0;
	for(chunk=0; chunk<numChunks; chunk++) {
		indices[3] = chunk;
		for(layer=0; layer<numLayers; layer++) {
			indices[2] = layer;
			for(column=0; column<numColumns; column++) {
				if (gCallSpinProcess && SpinProcess()) {		// Spins cursor and allows background processing.
					result = -1;								// User aborted.
					break;
				}
				indices[1] = column;
				for(row=0; row<numRows; row++) {
					indices[0] = row;
					if (result = MDGetTextWavePointValue(wavH, indices, textH))
						goto done;
					if (result = DoAppendAndPrepend(textH, p->prependStringH, p->appendStringH))
						goto done;
					if (result = MDSetTextWavePointValue(wavH, indices, textH))
						goto done;
				}
			}
			if (result != 0)
				break;
		}
		if (result != 0)
			break;
	}
	
done:
	if (wavH != NIL)
		WaveHandleModified(wavH);				// Inform Igor that we have changed the wave.
	if (textH != NIL)
		DisposeHandle(textH);
	if (p->prependStringH)
		DisposeHandle(p->prependStringH);		// We need to get rid of input parameters.
	if (p->appendStringH)
		DisposeHandle(p->appendStringH);		// We need to get rid of input parameters.
	p->result = 0;								// The Igor function result is always zero.
	
	return(result);
}
*/
/*	RegisterFunction()
	
	Igor calls this at startup time to find the address of the
	XFUNCs added by this XOP. See XOP manual regarding "Direct XFUNCs".
*/
static long
RegisterFunction()
{
	int funcIndex;

	funcIndex = GetXOPItem(0);		// Which function is Igor asking about?
	switch (funcIndex) {
		case 0:						// MicroArrayFitFunc(wave)
			return((long)MicroArrayFitFunc);
			break;
/*		case 1:						// WAFill3DWaveDirectMethod(wave)
			return((long)WAFill3DWaveDirectMethod);
			break;
		case 2:						// WAFill3DWavePointMethod(wave)
			return((long)WAFill3DWavePointMethod);
			break;
		case 3:						// WAFill3DWaveStorageMethod(wave)
			return((long)WAFill3DWaveStorageMethod);
			break;
		case 4:						// WAModifyTextWave(wave, prependString, appendString)
			return((long)WAModifyTextWave);
			break;*/
	}
	return(NIL);
}

/*	DoFunction()
	
	Igor calls this when the user invokes one if the XOP's XFUNCs
	if we returned NIL for the XFUNC from RegisterFunction. In this
	XOP, we always use the direct XFUNC method, so Igor will never call
	this function. See XOP manual regarding "Direct XFUNCs".
*/
static int
DoFunction()
{
	int funcIndex;
	void *p;				// Pointer to structure containing function parameters and result.
	int err;

	funcIndex = GetXOPItem(0);		// Which function is being invoked ?
	p = (void *)GetXOPItem(1);		// Get pointer to params/result.
	switch (funcIndex) {
		case 0:						// MicroArrayFitFunc(wave)
			err = MicroArrayFitFunc(p);
			break;
/*		case 1:						// WAFill3DWaveDirectMethod(wave)
			err = WAFill3DWaveDirectMethod(p);
			break;
		case 2:						// WAFill3DWavePointMethod(wave)
			err = WAFill3DWavePointMethod(p);
			break;
		case 3:						// WAFill3DWaveStorageMethod(wave)
			err = WAFill3DWaveStorageMethod(p);
			break;
		case 4:						// WAModifyTextWave(wave, prependString, appendString)
			err = WAModifyTextWave(p);
			break;*/
	}
	return(err);
}

/*	XOPEntry()

	This is the entry point from the host application to the XOP for all messages after the
	INIT message.
*/
static void
XOPEntry(void)
{	
	long result = 0;

	switch (GetXOPMessage()) {
		case FUNCTION:						// Our external function being invoked ?
			result = DoFunction();
			break;

		case FUNCADDRS:
			result = RegisterFunction();
			break;
	}
	SetXOPResult(result);
}

/*	main(ioRecHandle)

	This is the initial entry point at which the host application calls XOP.
	The message sent by the host must be INIT.
	main() does any necessary initialization and then sets the XOPEntry field of the
	ioRecHandle to the address to be called for future messages.
*/
HOST_IMPORT void
main(IORecHandle ioRecHandle)
{	
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	if (igorVersion < 500)				// Requires Igor Pro 5.00 or later.
		SetXOPResult(OLD_IGOR);			// OLD_IGOR is defined in WaveAccess.h and there are corresponding error strings in WaveAccess.r and WaveAccessWinCustom.rc.
	else
		SetXOPResult(0L);
}

