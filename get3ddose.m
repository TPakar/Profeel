function [dose, error, xbounds, ybounds, zbounds] = get3ddose(~, filename)
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
%% Collects .3ddose data as per instructed in STATDOSE documentation: I Kawrakow, E Mainegra-Hing, DWO Rogers, F Tessier, BRB Walters. 
%  The EGSnrc Code System: Monte Carlo simulation of electron and photon transport. 
%  Technical Report PIRS-701, National Research Council Canada (2017). 

dlgTitle = 'Valid question';
dlgQuestion = 'Load error data? (2xslower)! No will return 0 error matrix';
geterrorData = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');


fid = fopen(filename,'r');
%
% Check the number of voxels in xyz-dicections

len = fscanf(fid, '%i', 3);


% Read boudaries
xbounds = fscanf(fid, '%f', len(1)+1);
ybounds = fscanf(fid, '%f', len(2)+1);
zbounds = fscanf(fid, '%f', len(3)+1);


% Define array sizes
dslice = zeros(len(1), len(2));
eslice = dslice;
dose = zeros(len(2), len(1), len(3));
error = dose;

% Read the dose and error values

% 
for k = 1:len(3) % load the dose matrix
    dslice = fscanf(fid,'%e',[len(1), len(2)]);
    dose(:,:,k) = dslice';
end

if strcmp(geterrorData, 'Yes')
    for k = 1:len(3) % load the error matrix
        eslice = fscanf(fid,'%e',[len(1), len(2)]);
        error(:,:,k) = eslice';
    end
end

toc