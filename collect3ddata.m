function [newdata] = collect3ddata(name, Parentname, energy, fieldsize, tempmat, datatype, datadirect, filecount, normalvaldist, normalvalperc, normtoval)
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

% Write header info
newdata.(['NEW',num2str(filecount)]).name = name;
newdata.(['NEW',num2str(filecount)]).parentname = Parentname;
newdata.(['NEW',num2str(filecount)]).energy = energy(end);
newdata.(['NEW',num2str(filecount)]).fieldsize = fieldsize;
newdata.(['NEW',num2str(filecount)]).pos = tempmat(:,1);
newdata.(['NEW',num2str(filecount)]).datatype = datatype;
newdata.(['NEW',num2str(filecount)]).directory = datadirect;

idx2 = [];
% Interpolate data to find normalization constants. Interpolation
datainterp = [];
midax = ceil(length(tempmat(:,1))/2);
try
    interv = abs(tempmat(midax,1) - tempmat(midax-1,1));
catch
    disp('Error. Data has too few slices to get slice interval/thickness');
end
interpdens = interv/10;

datainterp(:,2) = interp1(tempmat(:,1), tempmat(:,2), tempmat(1,1):interpdens:tempmat(end,1));
datainterp(:,1) = tempmat(1,1):interpdens:tempmat(end,1);
%newdata.(['NEW',num2str(filecount)]).interpolated  = datainterp;

% CAX normalization
try
    [~,idx2] = min(abs(datainterp(:,1)));
    tempmat2 = tempmat(:,2)/datainterp(idx2,2);
    newdata.(['NEW',num2str(filecount)]).dataCAX = tempmat2;
catch
   try
        disp('Cannot find 0-point from the data');
        newdata.(['NEW',num2str(filecount)]).dataCAX = tempmat(:,2);  
   catch
       disp('Error in CAX normalization')
       newdata.(['NEW',num2str(filecount)]).dataCAX = tempmat(:,2);  
   end
end
datainterpcax = datainterp;
% Normalize interpolated data also to CAX (If cax is defined)
if ~isempty(idx2)
    datainterpcax(:,2) = datainterp(:,2)/datainterp(idx2,2);
end
newdata.(['NEW',num2str(filecount)]).interpolatedCAX  = datainterpcax;
newdata.(['NEW',num2str(filecount)]).interpolated  = datainterp;

% MAX normalization
normvalMAX = max(tempmat(:,2));
tempmat2 = tempmat(:,2)./normvalMAX;
newdata.(['NEW',num2str(filecount)]).dataMAX = tempmat2; 

% No normalization
newdata.(['NEW',num2str(filecount)]).dataNON = tempmat(:,2);
% Normalization to a value (if existing)
try
    newdata.(['NEW',num2str(filecount)]).dataTOVAL = tempmat(:,2)/normtoval;
    disp('TOVAL normalization done');
catch
    newdata.(['NEW',num2str(filecount)]).dataTOVAL = tempmat(:,2);
    disp('Normalization to a value failed');
end
% Manual normalization (no interpolation -> finds the closest distance value)

[~, idxdist] = min(abs(datainterp(:,1) - normalvaldist));

newdata.(['NEW',num2str(filecount)]).dataMAN = normalvalperc.*tempmat(:,2)/datainterp(idxdist,2);

newdata.(['NEW',num2str(filecount)]).dataMEAN = normalvalperc.*tempmat(:,2)/datainterp(idxdist,2);

filecount = filecount + 1;






