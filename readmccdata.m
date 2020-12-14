function [mccdat] = readmccdata(files, normalvaldist, normalvalperc)
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
mccdat = struct;
interpsteps = 10;

filecount = 1;
for i = 1:length(files)
   fid = fopen(files{i});
   % Get the first line
   line = fgetl(fid);
   % Loop until the end of the file and collect different files on the way
   while ~contains(line, "END_SCAN_DATA")
      if contains(line, 'BEGIN_SCAN')
          temp = strsplit(files{i},'\');
          tempnam = strrep(temp{end},'mcc', '');
          mccdat.(['mcc', num2str(filecount)]).name = tempnam;
          depthok = 0;
          crossplaneok = 0;
          inplaneok = 0;
          curvetype = '';
            while ~contains(line, 'END_SCAN')
               % Get header info
               line = fgetl(fid);
               % Get positions
               if contains(line, 'SCAN_DEPTH') || contains(line, 'SCAN_OFFAXIS_INPLANE') || contains(line, 'SCAN_OFFAXIS_CROSSPLANE') || contains(line, 'SCAN_CURVETYPE')
                  splithere = find(line == '='); 
                  linestart = line(1:splithere-1);
                  linestart = regexprep(linestart, ' ', '');
                  linestart = regexprep(linestart, '\t', '');
                  %try
                  switch linestart
                      case 'SCAN_DEPTH'
                            depth = num2str(str2double(line(splithere+1:end))/10);
                            depth = regexprep(depth, ' ', '');
                            depth = regexprep(depth, '\t', '');
                            mccdat.(['mcc', num2str(filecount)]).scandepth = depth;
                            depthok = 1;
                      case 'SCAN_OFFAXIS_INPLANE'
                            inplane = num2str(str2double(line(splithere+1:end))/10);
                            inplane = regexprep(inplane, ' ', '');
                            inplane = regexprep(inplane, '\t', '');
                            mccdat.(['mcc', num2str(filecount)]).scanoffaxis_inplane =  inplane;
                            inplaneok = 1;
                      case 'SCAN_OFFAXIS_CROSSPLANE'
                            crossplane = num2str(str2double(line(splithere+1:end))/10);
                            crossplane = regexprep(crossplane, ' ', '');
                            crossplane = regexprep(crossplane, '\t', '');
                            mccdat.(['mcc', num2str(filecount)]).scanoffaxis_crossplane =  crossplane;
                            crossplaneok = 1;
                      case 'SCAN_CURVETYPE'
                          if contains(line, 'PDD')
                            curvetype = 'pdd';
                          elseif contains(line, 'PROFILE')
                              curvetype = 'profile';
                          else
                              curvetype = line(splithere+1:end);
                          end
                  end
                  %catch
                  %   disp([linestart, ' not found from the header']); 
                 % end
               end
               
               
               % Get curve type
               if contains(line, 'BEGIN_DATA')
                   if inplaneok == 0
                       inplane = 'NaN';
                   end
                   if crossplaneok == 0
                       crossplane = 'NaN';
                   end
                   if depthok == 0
                       depth = 'NaN';
                   end
                   % add scan number to the name
                   mccdat.(['mcc', num2str(filecount)]).name = [mccdat.(['mcc', num2str(filecount)]).name,'_X:', crossplane, '_Y:', inplane,'_Z:', depth,'_SCAN_',num2str(filecount)];
                   mccdat.(['mcc', num2str(filecount)]).datatype = curvetype;
%                    temp = strsplit(line,'=');
%                    if contains(lower(temp{2}), 'profile')
%                         mccdat.(['mcc', num2str(filecount)]).datatype = 'profile';
%                    else
%                         mccdat.(['mcc', num2str(filecount)]).datatype = lower(temp{2});
%                    end
               end
               
               
               % Get data units
               if contains(line, 'MEAS_UNIT')
                  splithere = find(line == '='); 
                  dataunit = line(splithere+1:end);
                  dataunit = regexprep(dataunit, ' ', '');
                  dataunit = regexprep(dataunit, '\t', '');
                  mccdat.(['mcc', num2str(filecount)]).dataunit = dataunit;
               end
               
               
               % Get energy
               if contains(line, 'ENERGY')
                   temp = strsplit(line,'=');
                   mccdat.(['mcc', num2str(filecount)]).Energy = temp{2};
               end 
               % Get SSD
                if contains(line, 'SSD')
                   temp = strsplit(line,'=');
                   % Distance unit in mm -> cm
                   mccdat.(['mcc', num2str(filecount)]).SSD = str2double(temp{2})/10;
               end 
               % Get nominal field size (from inplane axis)
                if contains(line, 'INPLANE_FIELDSIZE')
                   temp = strsplit(line,'=');
                   % Distance unit in mm -> cm
                   mccdat.(['mcc', num2str(filecount)]).fieldsize = temp{2};
                end 
               if contains(line, 'BEGIN_DATA') && ~contains(line, 'SCAN')
                  line = fgetl(fid);
                  tempmat = [];
                  reference = [];
                  while ~contains(line,'END_DATA') 
                     templine = strsplit(line,' ');
                     templine2 = strsplit(line, '\t');
                     valfound = 0;
                     if length(templine) > length(templine2)
                        % Delete empty cells
                        cleanedline = templine(~cellfun('isempty',templine));                       
                     else
                        cleanedline = templine2(~cellfun('isempty',templine2));       
                     end
                       
                     % Data has 3 columns (1 = pos, 2 = data, 3 = reference)
                     if length(cleanedline) == 3
                         tempmat(end+1,1) = str2double(cleanedline{1});
                         tempmat(end,2) = str2double(cleanedline{2}); 
                         reference(end+1,1) = str2double(cleanedline{3});
                     else
                         tempmat(end+1,1) = str2double(cleanedline{1});
                         tempmat(end,2) = str2double(cleanedline{2}); 
                     end
                     

%                         templine = strsplit(line,' ');
%                         templine2 = strsplit(line, '\t');
%                         valfound = 0;
%                         if length(templine) > length(templine2)
%                             while valfound == 0
%                                 if ~isempty(templine{kk})
%                                     tempmat(end+1,1) = str2double(templine{kk});
%                                     valfound = 1;
%                                 end
%                                 kk = kk + 1;
%                             end
%                             tempmat(end,2) = str2double(templine{end});
%                         else
%                             kk = 1;
%                             while valfound == 0
%                                 if ~isempty(templine2{kk})
%                                     tempmat(end+1,1) = str2double(templine2{kk});
%                                     valfound = 1;
%                                 end
%                                 kk = kk + 1;
%                             end
%                             tempmat(end,2) = str2double(templine2{end});
%                         end
                  line = fgetl(fid);
                  end
                  % In mm
                  mccdat.(['mcc', num2str(filecount)]).pos = tempmat(:,1)/10;
                  
                    % Compute different normalizations and field parameters
                    % Interpolate data to find normalization constants. Interpolation
                    % density of 0.01 evalueted as sufficient
                    datainterp = [];
                    
                    midax = floor(length(mccdat.(['mcc', num2str(filecount)]).pos)/2);
                    interv = abs(mccdat.(['mcc', num2str(filecount)]).pos(midax) - mccdat.(['mcc', num2str(filecount)]).pos(midax-1));
                    interpdens = interv/interpsteps; 
                    
                    datainterp(:,2) = interp1(tempmat(:,1)/10, tempmat(:,2), tempmat(1,1)/10:interpdens:tempmat(end,1)/10);
                    datainterp(:,1) = tempmat(1,1)/10:interpdens:tempmat(end,1)/10;
                    disp(['Interpolation set to ', num2str(interpdens), 'cm']);
                    
                    % CAX normalization
                    idx2 = [];
                    if ~contains(mccdat.(['mcc', num2str(filecount)]).datatype,'pdd')
                        idx2 = find((datainterp(:,1) == 0));
                        tempmat2 = tempmat(:,2)./datainterp(idx2,2);
                        mccdat.(['mcc',num2str(filecount)]).dataCAX = tempmat2;      
                    else
                       try
                            disp('PDD-data -> CAX normalization not done');
                            mccdat.(['mcc',num2str(filecount)]).dataCAX = tempmat(:,2);  
                       catch
                           disp('Error in CAX normalization. No zero-point. -> No CAX normalization')
                           mccdat.(['mcc',num2str(filecount)]).dataCAX = tempmat(:,2);  
                       end
                    end
                    % Normalize also interpolated data to CAX
                    try
                        datainterpcax = datainterp;
                        datainterpcax(:,2) = datainterp(:,2)/datainterp(idx2,2);
                        
                    catch
                        disp('Interpolated data not normalized');
                        datainterpcax = datainterp;
                        datainterpcax(:,2) = datainterp(:,2);
                        
                    end
                    mccdat.(['mcc',num2str(filecount)]).interpolated  = datainterp;
                    mccdat.(['mcc',num2str(filecount)]).interpolatedCAX  = datainterpcax;
                    
                    % MAX normalization
                    normvalMAX = max(tempmat(:,2));
                    tempmat2 = tempmat(:,2)./normvalMAX;
                    mccdat.(['mcc',num2str(filecount)]).dataMAX = tempmat2; 

                    % No normalization
                    mccdat.(['mcc',num2str(filecount)]).dataNON = tempmat(:,2);

                    % Manual normalization (no interpolation -> finds the closest distance value)

                    [~, idxdist] = min(abs(datainterp(:,1) - normalvaldist));

                    mccdat.(['mcc',num2str(filecount)]).dataMAN = normalvalperc.*tempmat(:,2)/datainterp(idxdist,2);
                    if ~isempty(reference)
                        mccdat.(['mcc',num2str(filecount)]).reference = reference/100;
                        disp('Note: Reference is not in right units');
                    end
                    filecount = filecount + 1;
               end

            end
           line = fgetl(fid);
      end
   
   end
   
   
   fclose(fid);
end
save('mcctest.mat', 'mccdat');

