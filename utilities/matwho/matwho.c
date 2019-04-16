/*
 * MATWHO - Lists the variables in the specified .MAT file.
 *    MATWHO is a MEX file wrapper around the MEX function 'matGetDir'
 *    provided in the 'MAT-File Library' by The MathWorks. It lists the
 *    variables stored in a .MAT file, similar to using the 'who' function
 *    with the location set to '-file', but much faster, particularly on
 *    large .MAT files. MATWHO does not support pattern matching on the
 *    variable names, it only returns the complete list of variables.
 *    It is equivalent to calling: c = who('-file','filename');
 *
 *    This function resulted from a discussion on comp.soft-sys.matlab.
 *    My thanks to Friedrich and James Tursa for their contributions and feedback.
 *    http://www.mathworks.com/matlabcentral/newsreader/view_thread/316593
 *
 * USAGE:
 *   c = matwho('filename')
 *
 * INPUT:
 *   filename - Name of the '.mat' file to read.
 *
 * OUTPUT:
 *   c = Cell array of strings that correspond to each variable name.
 *
 * See also WHO
 *
 * Author: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
 * Last Modified: $Date: 2012-04-24 12:47:59 -0400 (Tue, 24 Apr 2012) $
 * $Id: matwho.c 3887 2012-04-24 16:47:59Z bkraus $
 */

#include "mex.h"
#include "mat.h"

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    int ndir, i;
    char *fname;
    const char **dir;
    MATFile *pmat;
  
    /* Check that inputs and outputs are OK */
    if(nrhs != 1)
        mexErrMsgIdAndTxt("matwho:usage","Wrong number of input arguments.\nUsage: c = matwho(filename)");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt("matwho:usage","Too many output arguments.\nUsage: c = matwho(filename)");
    else if(!mxIsChar(prhs[0]))
        mexErrMsgIdAndTxt("matwho:usage","First argument must be a filename (string).\nUsage: c = matwho(filename)");
    
    /* Copy the first input argument to the filename. */
    fname = mxArrayToString(prhs[0]);
    
    /* Open file for reading */
    pmat = matOpen(fname, "r");
    if (pmat == NULL)
        mexErrMsgIdAndTxt("matwho:fileerror","Error opening file: %s\n", fname);

    /* Read mat file directory */
    dir = (const char **)matGetDir(pmat, &ndir);
    if (dir == NULL) {
        matClose(pmat);
        mexErrMsgIdAndTxt("matwho:fileerror","Error reading directory of file: %s\n", fname);
    } else {
        plhs[0] = mxCreateCellMatrix(ndir,1);
        for (i=0; i < ndir; i++)
            mxSetCell(plhs[0], i, mxCreateString(dir[i]));
    }
    mxFree(dir);
    mxFree(fname);
    
    /* Close file */
    if (matClose(pmat) != 0)
        mexWarnMsgIdAndTxt("matwho:fileerror","Error closing file: %s\n", fname);
}
