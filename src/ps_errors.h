#ifndef PS_ERRORS_H
#define PS_ERRORS_H

int ps_errno;
/** 
 * Header file defining error codes
 */
enum  {
	PS_OK,
	PS_ERROR_DIMENSION_MISMATCH,
	PS_ERROR_INDEX_OUT_OF_BOUNDS,
	PS_ERROR_UNEXPECTED_EOF, 
	PS_ERROR_FILE_NOT_FOUND
};

#endif