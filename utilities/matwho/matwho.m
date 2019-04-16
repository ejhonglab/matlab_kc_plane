function c = matwho(filename)
% MATWHO - Lists the variables in the specified .MAT file.
%    MATWHO is a MEX file wrapper around the MEX function 'matGetDir'
%    provided in the 'MAT-File Library' by The MathWorks. It lists the
%    variables stored in a .MAT file, similar to using the 'who' function
%    with the location set to '-file', but much faster, particularly on
%    large .MAT files. MATWHO does not support pattern matching on the
%    variable names, it only returns the complete list of variables.
%    It is equivalent to calling: c = who('-file','filename');
%
%    This function resulted from a discussion on comp.soft-sys.matlab.
%    My thanks to Friedrich and James Tursa for their contributions and feedback.
%    http://www.mathworks.com/matlabcentral/newsreader/view_thread/316593
%
% USAGE:
%   c = matwho('filename')
%
% INPUT:
%   filename - Name of the '.mat' file to read.
%
% OUTPUT:
%   c = Cell array of strings that correspond to each variable name.
%
% See also WHO

% Author: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Last Modified: $Date: 2012-06-22 15:44:50 -0400 (Fri, 22 Jun 2012) $
% $Id: matwho.m 4089 2012-06-22 19:44:50Z bkraus $
cfilename = [mfilename '.c'];
wstr = sprintf('\nUsing native (and much slower) ''who'' function instead.');
if(exist(cfilename,'file')==2)
    warning([mfilename, ':UncompiledMex'], '''%s'' needs to be compiled as a MEX function.%s',cfilename,wstr);
else warning([mfilename, ':NoMex'], 'Cannot find ''matwho'' MEX function.%s',wstr);
end

% Call native who with the '-file' argument.
c = who('-file',filename);
