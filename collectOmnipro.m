function [omnidata] = collectOmnipro(omnidir, varargin)
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
%% Collects and normalizes omniprodata to structures

omnidata = struct;
filecount = 1;
flipy = false;
flipx = false;
divv = 1;
if ~isempty(varargin)
   normalvaldist = varargin{1}/10;
   normalvalperc = varargin{2}/100;
   caxcorr = varargin{3};
else
   normalvalperc = 100;
   normalvaldist = 0;
end

for k = 1:length(omnidir)
   if omnidir(k).bytes > 1
       linecount = 1;
       line = 'asd';
       % Some random big number
       collectdata = 1000000000;
       ycount = 1;
       fid = fopen([omnidir(1).folder, '\', omnidir(k).name]);
       im = '';
       % Initialize energy and fieldsize, since they are required later
       omnidata.(['omni', num2str(filecount)]).Energy = 'Unknown';
       omnidata.(['omni', num2str(filecount)]).FieldSize = 'Unknown';
       lengthunit = 'cm';
       
       while ~contains(line, '</asciibody>')
           if contains(line, 'Image')
              imtemp = strsplit(line);
              nameidx = find(contains(imtemp, 'Image'));
              im = cell2mat(imtemp(nameidx+2:end));
              omnidata.(['omni', num2str(filecount)]).name = [omnidir(k).name]; %[omnidata.(['omni', num2str(filecount)]).name,' ', im];
             % omnidata.(['omni', num2str(filecount+1)]).name = [omnidata.(['omni', num2str(filecount+1)]).name,' ', im];
              %omnidata.(['omni', num2str(filecount+2)]).name = [omnidata.(['omni', num2str(filecount+2)]).name,' ', im];
           end
           line = fgetl(fid);
           
           % Separator
           if contains(line, 'Separator:')
               if contains(line, '[TAB]')
                    separator = '\t';
               elseif contains(line, ',')
                    separator = ',';
               else
                    splitted = strsplit(line, ':');
                    separator = strrep(splitted{end}, ' ', '');
               end
           end
           
           % Matrix size
           if contains(line, 'No. of Columns:')
               splittedline = strsplit(line, ':');
               line = strrep(splittedline{end},' ', '');
               matcolumns = str2double(line);
           end
           if contains(line, 'No. of Rows:')
               splittedline = strsplit(line, ':');
               line = strrep(splittedline{end},' ', '');
               matrows = str2double(line);
           end
           % Axis unit
           if contains(line, 'cm')
               divv = 1;
               lengthunit = '[cm]';
           elseif contains(line, 'mm')
               divv = 10;
               lengthunit = '[mm]';
           end

           % Collect the header
           if contains(line, 'File Name:')
               if ~contains(line, '\')
                   splittedline = strsplit(line, ':');
                   omnidata.(['omni', num2str(filecount)]).name = [strrep(splittedline{end}, ' ', ''), '_XY axis ', im];
                   %omnidata.(['omni', num2str(filecount+1)]).name = [strrep(splittedline{end}, ' ', ''), '_X axis ', im];
                   %omnidata.(['omni', num2str(filecount+2)]).name = [strrep(splittedline{end}, ' ', ''), '_Y axis ', im];
               else
                   splittedline = strsplit(line, '\');
                   omnidata.(['omni', num2str(filecount)]).name = [strrep(splittedline{end}, ' ', ''), ' XY axis ', im];
                  % omnidata.(['omni', num2str(filecount+1)]).name = [strrep(splittedline{end}, ' ', ''), ' X axis ', im];
                   %omnidata.(['omni', num2str(filecount+2)]).name = [strrep(splittedline{end}, ' ', ''), ' Y axis ', im];
               end

           elseif contains(line, ':')
                %splittedline = {};
                splittedline = strsplit(line, ':');
                temp = strrep(splittedline{1}, ' ', '');
                omnidata.(['omni', num2str(filecount)]).(temp) = strrep(splittedline{2}, ' ', '');
               % omnidata.(['omni', num2str(filecount+1)]).(temp) = strrep(splittedline{2}, ' ', '');
               % omnidata.(['omni', num2str(filecount+2)]).(temp) = strrep(splittedline{2}, ' ', '');

           end
           
           
           
           if contains(line, ['X', lengthunit])
               % Collect x positions
               collectdata = linecount + 2;
               datamatrix = zeros(matrows,matcolumns);
               %xaxis = zeros(32,1);
               ypos = zeros(matrows,1);
               % Delete spaces
               temp = strrep(line, ' ','');
               temp2 = strsplit(temp, separator);
               % Split data
               xpos = str2double(temp2(2:end-1));
               
           elseif linecount >= collectdata && ~contains(line, '</asciibody>')
               % Collect y positions
               temp = strrep(line, ' ','');
               temp2 = strsplit(temp, separator);
               ypos(ycount) = str2double(temp2(1,1));
               
               % Collect matrix data
               datamatrix(ycount,:) = str2double(temp2(2:end-1));         
               ycount = ycount + 1;    
           end
           linecount = linecount + 1;
       end
       

       if mod(matrows,2) == 0 && mod(matcolumns,2) == 0
            xave = (datamatrix(matrows/2,:) + datamatrix(matrows/2+1,:))/2;
            yave = (datamatrix(:, matcolumns/2) + datamatrix(:, matcolumns/2+1))/2;
           elseif mod(matrows,2) == 0 && mod(matcolumns,2) ~= 0
                xave = (datamatrix(matrows/2, :) + datamatrix(matrows/2+1, :))/2;
                yave = datamatrix(:, floor(matcolumns/2)+1);
           elseif mod(matcolumns,2) == 0 && mod(matrows,2) ~= 0
                xave = datamatrix(floor(matrows/2)+1, :);
                yave = (datamatrix(:, matcolumns/2) + datamatrix(:, matcolumns/2+1))/2;
           else
               xave = datamatrix(floor(matrows/2)+1, :);
               yave = datamatrix(:, floor(matcolumns/2)+1);
       end
       %disp('Mid axis profiles created automatically (Accuracy same as the measurement accuracy)');
       
       
       % Interpolate data to find normalization constants. 0.01 interpolation
       % value evaluated as sufficient, though arbitrary
        datainterpX = [];
        datainterpY = [];
        % For X
        % Check that the data direction is not flipped
         if xpos(end)<xpos(1)
             xpos = flip(xpos);
             xave = flip(xave);
             flipx = true;
         end
        
        
        xpos = xpos/divv;
        ypos = ypos/divv;
         
        midax = floor(length(xpos)/2);
        interv = abs(xpos(midax) - xpos(midax-1));
        interpdens = interv/10; 
         
        datainterpX(:,2) = interp1(xpos, xave, xpos(1):interpdens:xpos(end));
        datainterpX(:,1) = xpos(1):interpdens:xpos(end);
        disp(['Interpolation set to ', num2str(interpdens), 'cm']);
        %omnidata.(['omni',num2str(filecount+1)]).interpolated  = datainterpX; % In cm
        % For Y
        % Check that the data direction is not flipped
        if ypos(end)<ypos(1)
            ypos = flip(ypos);
            yave = flip(yave);
            flipy = true;
        end
        midax = floor(length(ypos)/2);
        interv = abs(ypos(midax) - ypos(midax-1));
        interpdens = interv/10; 
        
        datainterpY(:,2) = interp1(ypos, yave, ypos(1):interpdens:ypos(end));
        datainterpY(:,1) = ypos(1):interpdens:ypos(end);
        
       % Different normalizations
       
       % normalized to CAX (0-point from data)
       % X axis
       try
           [~, idxcaxx] = min(abs(datainterpX(:,1)));
           normvalCAXx = datainterpX(idxcaxx, 2);
       catch
            normvalCAXx = datainterpX(round(length(datainterpX(:,2))/2),2);
            disp('Error while searching mid axis (X) CAX value form Omnipro data -> CAX set to axis middle point');
       end
       % Y axis
       try
           [~, idxcaxy] = min(abs(datainterpY(:,1)));
           normvalCAXy = datainterpY(idxcaxy, 2);
       catch
            normvalCAXy = datainterpY(round(length(datainterpY(:,2))/2),2);
            disp('Error while searching mid axis (Y) CAX value form Omnipro data -> CAX set to axis middle point');
       end
       % Matrix
       try
           [devx, idxx] = min(abs(xpos));
           [devy, idxy] = min(abs(ypos));
           normvalCAXm = datamatrix(idxy, idxx);
       catch
           [sizey, sizex] = size(datamatrix);
           normvalCAXm = datamatrix(round(sizey/2), round(sizex/2));
           disp('Error in finding 2D 0-point. 2D data CAX normalized to midpoint');
       end
       

       %normval = (xave(16) + xave(17))/2;
       
    
       CAXnormalizedMdata = datamatrix/normvalCAXm;

       % Normalized to MAX
       normval = max(xave);
       normval = max(yave);

       normvalMAX = max(max(datamatrix));
       MAXnormalizedMdata = datamatrix/normvalMAX;
               
       % No normalization
       UNnormalizedXdata = xave;
       UNnormalizedYdata = yave;
       UNnormalizedMdata = datamatrix;
       
       
       % Manual
       % Normalization constant
       % For X
       [~, idxdistx] = min(abs(datainterpX(:,1) - normalvaldist));
       MANnormalizedXdata = normalvalperc*xave(:)/datainterpX(idxdistx,2);
        
       % Manual normalization (finds the closest distance value <- interpolated before
       % For Y
       [~, idxdisty] = min(abs(datainterpY(:,1) - normalvaldist));
       MANnormalizedYdata = normalvalperc*yave(:)/datainterpY(idxdisty,2);
       
       
       % Manual normalization
       try
            MANnormalizedMdata = normalvalperc*datamatrix/((datainterpX(idxdistx,2) + datainterpY(idxdisty,2))/2);
       catch
            disp('Error in initial 2D manual normalization. Renormalize manually');
            MANnormalizedMdata = normalvalperc*datamatrix;
       end

       

       omnidata.(['omni', num2str(filecount)]).dataCAX = flip(CAXnormalizedMdata);
       omnidata.(['omni', num2str(filecount)]).dataMAX = flip(MAXnormalizedMdata);
       omnidata.(['omni', num2str(filecount)]).dataNON = flip(UNnormalizedMdata);
       omnidata.(['omni', num2str(filecount)]).dataMAN = flip(MANnormalizedMdata);
       
       % Position vectors
       omnidata.(['omni', num2str(filecount)]).xpos = xpos';
       omnidata.(['omni', num2str(filecount)]).ypos = ypos;
      % omnidata.(['omni', num2str(filecount+1)]).pos = xpos';
      % omnidata.(['omni', num2str(filecount+2)]).pos = ypos;
 
       % 2D
       %[Xq, Yq] = meshgrid(xpos(1):0.01:xpos(end),ypos(1):0.01:ypos(end));
       %datamatrix1 = interp2(xpos,ypos,CAXnormalizedMdata, Xq, Yq);
        maxdimsize = 250;
        xtrunc = xpos;
        ytrunc = ypos;
        truncate = 0;
        

        % Truncate size for display, if data is too big
        if length(xpos) > maxdimsize
            xstep = floor(length(xpos)/maxdimsize);
            xtrunc = xpos(1:xstep:length(xpos));
            truncate = 1;
        else
            xstep = 1;
        end

        % Truncate size for display, if data is too big
        if length(ypos) > maxdimsize
            ystep = floor(length(ypos)/maxdimsize);
            ytrunc = ypos(1:ystep:length(ypos));
            truncate = 1;
        else
            ystep = 1;
        end
       % Truncate size for display, if data is too big (e.g. > 100/dim)
        if truncate == 1
            disp('Data truncated to speedup display');
            
            VdoseCAX = CAXnormalizedMdata(1:ystep:length(CAXnormalizedMdata(:,1)), 1:xstep:length(CAXnormalizedMdata(1,:)));
            VdoseMAX = MAXnormalizedMdata(1:ystep:length(MAXnormalizedMdata(:,1)), 1:xstep:length(MAXnormalizedMdata(1,:)));
            VdoseNON = UNnormalizedMdata(1:ystep:length(UNnormalizedMdata(:,1)), 1:xstep:length(UNnormalizedMdata(1,:)));
            VdoseMAN = MANnormalizedMdata(1:ystep:length(MANnormalizedMdata(:,1)), 1:xstep:length(MANnormalizedMdata(1,:)));

            omnidata.(['omni', num2str(filecount)]).DisplaydataCAX = flip(VdoseCAX);
            omnidata.(['omni', num2str(filecount)]).DisplaydataMAX = flip(VdoseMAX);
            omnidata.(['omni', num2str(filecount)]).DisplaydataNON = flip(VdoseNON);
            omnidata.(['omni', num2str(filecount)]).DisplaydataMAN = flip(VdoseMAN);
            omnidata.(['omni', num2str(filecount)]).Displayxpos = xtrunc';
            omnidata.(['omni', num2str(filecount)]).Displayypos = ytrunc;
            omnidata.(['omni', num2str(filecount)]).truncate = 1;
        else
            omnidata.(['omni', num2str(filecount)]).DisplaydataCAX = flip(CAXnormalizedMdata);
            omnidata.(['omni', num2str(filecount)]).DisplaydataMAX = flip(MAXnormalizedMdata);
            omnidata.(['omni', num2str(filecount)]).DisplaydataNON = flip(UNnormalizedMdata);
            omnidata.(['omni', num2str(filecount)]).DisplaydataMAN = flip(MANnormalizedMdata);
            omnidata.(['omni', num2str(filecount)]).Displayxpos = xpos';
            omnidata.(['omni', num2str(filecount)]).Displayypos = ypos;
            omnidata.(['omni', num2str(filecount)]).truncate = 0;
        end
        
        
       
       datainterpcaxX = datainterpX;
       datainterpcaxY = datainterpY;
       
       
       datainterpcaxX(:,2) = datainterpX(:,2)/normvalCAXx;
       datainterpcaxY(:,2) = datainterpY(:,2)/normvalCAXy;
       
      % omnidata.(['omni', num2str(filecount+1)]).interpolatedCAX = datainterpcaxY;
      % omnidata.(['omni', num2str(filecount+2)]).interpolatedCAX = datainterpcaxX;
       
%        try
%         [caxdevX, symX, homX, dmaxX, dminX, devX, FWX, penRX, penLX] = fieldparams(datainterpcaxX(:,2), datainterpcaxX(:,1), 'omni', caxcorr);
%         [caxdevY, symY, homY, dmaxY, dminY, devY, FWY, penRY, penLY] = fieldparams(datainterpcaxY(:,2), datainterpcaxY(:,1), 'omni', caxcorr);
%         paramsx = [caxdevX, symX, homX, dmaxX, dminX, devX, FWX, penRX, penLX];
%        paramsy = [caxdevY, symY, homY, dmaxY, dminY, devY, FWY, penRY, penLY];
%        
%        omnidata.(['omni', num2str(filecount+1)]).params = round(paramsx,4);
%        omnidata.(['omni', num2str(filecount+2)]).params = round(paramsy,4);
%        catch
%           disp('Profile fieldparameters not computed'); 
%        end
       
       
       
       % Datatypes
       omnidata.(['omni', num2str(filecount)]).datatype = '2D';
       %omnidata.(['omni', num2str(filecount+1)]).datatype = 'profile';
       %omnidata.(['omni', num2str(filecount+2)]).datatype = 'profile';
       
       fclose(fid);
       filecount = filecount + 1; %+ 3;
   end 
end


%save('Omnidata.mat','omnidata')

