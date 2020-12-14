function [newdata] = modifynewplane(planedata,hordim, vertdim, name, Parentname, mydir,energy, fieldsize, filecount, plane, normplane)
%%
% Copyright (C) 2020 Tomppa Pakarinen, tomppa.pakarinen@pshp.fi


% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/
%%

%% Computes needed values for new plane structure


newdata.(['NEW',num2str(filecount)]).datatype = '2D';
newdata.(['NEW',num2str(filecount)]).Plane = plane;
newdata.(['NEW',num2str(filecount)]).directory = mydir;
newdata.(['NEW',num2str(filecount)]).parentname = Parentname;
% plane name XY (convention) is reversed compared to XZ and ZY -> need to flip
if strcmp(plane,'XY')
    plane = flip(plane);
end
newdata.(['NEW',num2str(filecount)]).name = name;
newdata.(['NEW',num2str(filecount)]).Energy = energy(end);
newdata.(['NEW',num2str(filecount)]).FieldSize = fieldsize;
newdata.(['NEW',num2str(filecount)]).([lower(plane(1)), 'pos']) = vertdim;
newdata.(['NEW',num2str(filecount)]).([lower(plane(2)), 'pos']) = hordim;
newdata.(['NEW',num2str(filecount)]).dataNON = planedata;
newdata.(['NEW',num2str(filecount)]).dataMAX = planedata./max(planedata);
newdata.(['NEW',num2str(filecount)]).dataTOVAL = planedata/normplane;

[~, caxidxx] = min(abs(hordim));
[~, caxidxy] = min(abs(vertdim));

caxval = planedata(caxidxy, caxidxx);
newdata.(['NEW',num2str(filecount)]).dataCAX = planedata/caxval;
% Filling values for 2D manual normalization, now just replaced with CAX... Should be included
newdata.(['NEW',num2str(filecount)]).dataMAN = newdata.(['NEW',num2str(filecount)]).dataCAX;


truncate = 0;
% Truncate size for display, if data is too big (e.g. > 100/dim)



% Truncate data for display

maxdimsize = 120;
xtrunc = hordim;
ytrunc = vertdim;

% Truncate size for display, if data is too big
if length(hordim) > maxdimsize
    xstep = floor(length(hordim)/maxdimsize);
    xtrunc = hordim(1:xstep:end);
    truncate = 1;
else
    xstep = 1;
end

% Truncate size for display, if data is too big
if length(vertdim) > maxdimsize
    ystep = floor(length(vertdim)/maxdimsize);
    ytrunc = vertdim(1:ystep:end);
    truncate = 1;
else
    ystep = 1;
end

if truncate == 1
    
    
    disp('Data truncated to speedup display');
%     
%     if ~mod(length(ydim), ystep)
%        addy = 2;
%        ydim = ydim(2:end);
%     else
%        addy = 1;
%     end
%     if ~mod(length(xdim), ystep)
%        addx = 2;
%        xdim = xdim(2:end)
%     else
%        addx = 1;
%     end
    
    

    VdoseCAX = newdata.(['NEW',num2str(filecount)]).dataCAX(1:ystep:end, 1:xstep:end);
    VdoseMAX = newdata.(['NEW',num2str(filecount)]).dataMAX(1:ystep:end, 1:xstep:end);
    VdoseNON = newdata.(['NEW',num2str(filecount)]).dataNON(1:ystep:end, 1:xstep:end);
    VdoseMAN  = newdata.(['NEW',num2str(filecount)]).dataMAN(1:ystep:end, 1:xstep:end);
    size(newdata.(['NEW',num2str(filecount)]).dataCAX)
    newdata.(['NEW',num2str(filecount)]).DisplaydataCAX = VdoseCAX;
    newdata.(['NEW',num2str(filecount)]).DisplaydataMAX = VdoseMAX;
    newdata.(['NEW',num2str(filecount)]).DisplaydataNON = VdoseNON;
    newdata.(['NEW',num2str(filecount)]).DisplaydataMAN = VdoseMAN;
    newdata.(['NEW',num2str(filecount)]).DisplaydataTOVAL = VdoseNON/normplane;
    newdata.(['NEW',num2str(filecount)]).(['Display', ([lower(plane(1)), 'pos'])]) = ytrunc;
    newdata.(['NEW',num2str(filecount)]).(['Display', ([lower(plane(2)), 'pos'])]) = xtrunc;
    newdata.(['NEW',num2str(filecount)]).truncate = 1;
    
    
else
    newdata.(['NEW',num2str(filecount)]).DisplaydataCAX = newdata.(['NEW',num2str(filecount)]).dataCAX;
    newdata.(['NEW',num2str(filecount)]).DisplaydataMAX = newdata.(['NEW',num2str(filecount)]).dataMAX;
    newdata.(['NEW',num2str(filecount)]).DisplaydataNON = newdata.(['NEW',num2str(filecount)]).dataNON;
    newdata.(['NEW',num2str(filecount)]).DisplaydataMAN = newdata.(['NEW',num2str(filecount)]).dataMAN;
    newdata.(['NEW',num2str(filecount)]).DisplaydataTOVAL = newdata.(['NEW',num2str(filecount)]).dataNON/normplane;
    newdata.(['NEW',num2str(filecount)]).(['Display', ([lower(plane(1)), 'pos'])]) = vertdim;
    newdata.(['NEW',num2str(filecount)]).(['Display', ([lower(plane(2)), 'pos'])]) = hordim;
    newdata.(['NEW',num2str(filecount)]).truncate = 0;
end


