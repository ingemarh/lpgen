/* 
 * This function is needed in the Macintosh MEX-files to very
 * quickly round extended numbers on MC68040. This routine is 
 * a lot faster than 68040 FINTR/FINTRZ instruction.
 *
 */
 
#ifdef ANSI_C 
	int ffloor(double x)
#else
	int ffloor(x)
	double x;
#endif
	{
	int *ii;
	char *ic;
    extended temp;
	
	ic=&temp;
    ii=ic+6;
	temp = x + 140739635838976; /* 2^47+2^31 */
	return( *ii - 2147483648 ); /* 2^31 */
	}
