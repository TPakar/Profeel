function [PTWdata] = collectPTW(ptwdir, varargin)

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
%% Collects and normalized PTW data to structures (.txt export files)



if ~isempty(varargin)
   normalvaldist = varargin{1};
   normalvalperc = varargin{2}/100;
else
   normalvalperc = 1;
   normalvaldist = 0;
   disp('Ei manuaalista normalisointia');
end
% mm to cm
divv = 10;
% interpolation density (points per every datapoint)
interpsteps = 10;
filecount = 1;
datatype = 'profile';
tempmat = [];
tempheader = "d";
for k = 1:length(ptwdir)
if ptwdir(k).bytes > 1
    
    fid = fopen([ptwdir(1).folder, '\', ptwdir(k).name]);
    line = fgetl(fid);
    while (line ~= -1)                      % Read until end of the file
        if contains(line,'mm')
           divv = 10;
        elseif contains(line,'cm')
           divv = 1;
        end

        tempheader1 = line;
        tempheader = strsplit(tempheader1, ' ');
        name = ptwdir(k).name;
        if contains(lower(name),'pdd')
            datatype = 'pdd';
        else
            datatype = 'profile';
        end
        energy = tempheader([find(contains(tempheader, 'Energy')),find(contains(tempheader, 'Energy'))+2]);
        energy = strrep(energy,',', '.');
        fieldsize = tempheader([find(contains(tempheader, 'Field'))+5, find(contains(tempheader, 'Field'))+7]);
        fieldsize = strrep(fieldsize,',', '.');
        line = fgetl(fid);
        i = 1;
        tempmat = [];
        
        try
        while ~contains(line, 'Depth')
            if i > 1
                line = fgetl(fid);
                
            end
            if ~isempty(line) && ~contains(line, 'Depth')
                
                line = strrep(line, ',', '.');

                % Split the line with different delimiters
                templine = strsplit(line,' ');
                templine2 = strsplit(line, '\t');
                if length(templine) > length(templine2)
                    tempmat(i,1) = str2double(templine{1});
                    tempmat(i,2) = str2double(templine{2});
                else
                    tempmat(i,1) = str2double(templine2{1});
                    tempmat(i,2) = str2double(templine2{2});
                end
                i = i+1;
            end
            
        end
        catch
           disp('PTW valmis'); 
        end

        PTWdata.(['PTW',num2str(filecount)]).name = name;
        PTWdata.(['PTW',num2str(filecount)]).energy = energy(end);
        PTWdata.(['PTW',num2str(filecount)]).fieldsize = fieldsize;
        PTWdata.(['PTW',num2str(filecount)]).pos = tempmat(:,1)/divv;
        PTWdata.(['PTW',num2str(filecount)]).datatype = datatype;
        
        PTWdata.(['PTW',num2str(filecount)]).dataunit = 'Unknown';
        % Interpolate data to find normalization constants. Interpolation
        % density of 0.01 evalueted as sufficient
        datainterp = [];
        
        midax = floor(length(PTWdata.(['PTW',num2str(filecount)]).pos)/2);
        interv = abs(PTWdata.(['PTW',num2str(filecount)]).pos(midax) - PTWdata.(['PTW',num2str(filecount)]).pos(midax-1));
        interpdens = interv/interpsteps; 
        datainterp(:,2) = interp1(tempmat(:,1)/divv, tempmat(:,2), tempmat(1,1)/divv:interpdens:tempmat(end,1)/divv);
        datainterp(:,1) = tempmat(1,1)/divv:interpdens:tempmat(end,1)/divv;
        disp(['Interpolation set to ', num2str(interpdens), 'cm']);
        
        % CAX normalization
        idx2 = [];
        if ~contains(datatype,'pdd')
            [val,idx2] = min(abs(datainterp(:,1)));
            
            % Warn about 2 interp steps deviations from CAX
            if val > abs(interv*2)
               disp('Data does not have values near zero! Normalization to CAX might not be correct'); 
            end
            tempmat2 = tempmat(:,2)./datainterp(idx2,2);
            PTWdata.(['PTW',num2str(filecount)]).dataCAX = tempmat2;      
        else
           try
                disp('PDD-data -> CAX normalization skipped -> No normalization');
                PTWdata.(['PTW',num2str(filecount)]).dataCAX = tempmat(:,2);  
           catch
               disp('Error in CAX noramlization, no zero-point in the data -> No normalization')
               PTWdata.(['PTW',num2str(filecount)]).dataCAX = tempmat(:,2);
           end
        end
        % Normalize also interpolated data to CAX
        datainterpcax = [];
        try
            datainterpcax(:,2) = datainterp(:,2)/datainterp(idx2,2);
            datainterpcax(:,1) = datainterp(:,1);
        catch
            disp('Interpolated data not normalized. PDD?'); 
            datainterpcax(:,2) = datainterp(:,2);
            datainterpcax(:,1) = datainterp(:,1);
        end
        
        PTWdata.(['PTW',num2str(filecount)]).interpolated  = datainterp;
        PTWdata.(['PTW',num2str(filecount)]).interpolatedCAX  = datainterpcax;
        
        % MAX normalization
        normvalMAX = max(tempmat(:,2));
        tempmat2 = tempmat(:,2)./normvalMAX;
        PTWdata.(['PTW',num2str(filecount)]).dataMAX = tempmat2; 
        
        % No normalization
        PTWdata.(['PTW',num2str(filecount)]).dataNON = tempmat(:,2);
        
        % Manual normalization (no interpolation -> finds the closest distance value)
        try
            [~, idxdist] = min(abs(datainterp(:,1) - normalvaldist));
        catch 
            disp('Error in PTW normalization. Set distance value out of rahge: Set to half maximum');
            normalvaldist = max(datainterp(:,1))/2;
            [~, idxdist] = min(abs(datainterp(:,1) - normalvaldist));
        end
        PTWdata.(['PTW',num2str(filecount)]).dataMAN = normalvalperc.*tempmat(:,2)/datainterp(idxdist,2);

        filecount = filecount + 1;
    end
        fclose(fid);
end 
end
%save('testfilePTW.mat','PTWdata');


