function [dta, Drelative, Dabsolute, dtaA, DTApercent2, DTApercent3, Dpassperc] = DtaDComp(refstruct, measstruct, dosecrit, dtacrit, resolution, interval, dtamethod, gmask)
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
%% Computes DTA and dose difference values from given data (vectors same length)


measured = measstruct.data;
reference = refstruct.data;
DTApercent2 = [];
DTApercent3 = [];
%%
% DTAperc: Finds absolute difference vector between each value of reference
% and measured data -> finds all values that meet the dose criterion -> finds
% the minimum distance
matsize = [];
failedcountDTA = 0;
failedcountDTAA = 0;

% dta criterion to cm
dtacrit = dtacrit/10;

% 1D
% Check number of dimensions
dimnum = 0;
if size(reference,1) == 1 || size(reference,2) == 1
    dimnum = 1;
else
    dimnum = 2;
end


if dimnum == 1

% Brute force (interpolate points between points and find the closest)
len = length(measured);
dta = ones(1,len);
dtaA = ones(1,len);

if dtamethod(2) == 1
% Interpolates reference data further so that dose-equivalent points can be
% found with sufficient resolution
interp_ref = interp1(reference, 0:resolution:len);


% Find the closest value with resolution of 1/10000 from sampling accuracy.
for i = 1:len
  indices2 = find(abs(interp_ref - measured(i)) <= 0.0001);
  try
    % Find the minimum distance
    dta(i) = min(abs(indices2*resolution - i))*interval;
    if dta(i) > 2
       % Truncate high dta values e.g. to 2 cm;
       dta(i) = 2; 
    end
  catch
    % If no matching dose is found in whole reference vector -> insert NaN
    dta(i) = NaN;
  end
  
  if dta(i) > dtacrit || isnan(dta(i))
     failedcountDTA = failedcountDTA + 1;
  end
  %
end
%%
dta = dta';
DTApercent2 = (len - failedcountDTA)/len*100;

end

% 'Analytical' DTA computation. Find consequent points lower and higher than
% current value -> find intersection of the line defined by the points and
% yref = constant
if dtamethod(1) == 1
    dtaA = ones(1,len)*2;
    for i = 1:length(measured)
        % Find all the consequent points surrounding the current reference
        % point dose
        temp = [];
        [~, lowervalsidx] = find(reference <= measured(i));
        for j = 1:length(lowervalsidx)
           if lowervalsidx(j) + 1 <= length(reference)
               if reference(lowervalsidx(j)+1) >= measured(i)
                   % Find DTA for the point
                   % Create a line
              
                   y2 = reference(lowervalsidx(j)+1);
                   y1 = reference(lowervalsidx(j));
                   x1 = refstruct.pos(lowervalsidx(j));
                   x2 = refstruct.pos(lowervalsidx(j)+1);
                   x = x1;
                   
                   y0 = y1 - (y2-y1)/(x2-x1)*x;
                    
                   temp(end+1) = (measured(i) - y0)*(x2-x1)/(y2-y1);

               end
           end
            if lowervalsidx(j) - 1 > 0
               if reference(lowervalsidx(j) - 1) >= measured(i)
                   % Find DTA for the point
                   % Create a line
                   y2 = reference(lowervalsidx(j));
                   y1 = reference(lowervalsidx(j)-1);
                   x1 = refstruct.pos(lowervalsidx(j)-1);
                   x2 = refstruct.pos(lowervalsidx(j));
                   x = x1;
                   
                   y0 = y1 - (y2-y1)/(x2-x1)*x;

                   % Solve intersection with the ref point line
                   temp(end+1) = (measured(i) - y0)*(x2-x1)/(y2-y1);
               end
            end
        end
        try
            minval = min(abs(temp - measstruct.pos(i)));
            if minval > 2
                minval = 2;
            end
            dtaA(i) = minval;
        catch
           dtaA(i) = NaN;
        end
        if dtaA(i) > dtacrit || isnan(dtaA(i))
             failedcountDTAA = failedcountDTAA + 1;
        end
    end

DTApercent3 = (len - failedcountDTAA)/len*100;
end

%% Dose difference
Drelative = (measured - reference)./reference;    %./(reference*dosecrit/100);
Dabsolute = (measured - reference);
[Dpassperc, ~] = find(abs(Drelative)*100 < dosecrit);

Dpassperc = sum(Dpassperc)/length(Drelative)*100;

% 2D case
else
    % 'Brute force' (interpolate points between datapoints and find the closest values)
    [leny, lenx] = size(reference);
    dta = ones(leny,1);
    dtaA = ones(1,lenx);

    % Find dose values inside a circle, size of the seach distance
    % Brute force
    if dtamethod(2) == 1 || dtamethod(1) == 1 
        radiusx = round(1/interval(2));
        radiusy = round(1/interval(1));
        for i = 1:lenx
                for j = 1:leny
                        % Only compute if dose in referencedata exceeds the
                        % set threshold
                        if gmask(j,i) == 1
                            if abs(reference(j,i) - measured(j,i)) <= 0.0001
                                dta(j,i) = 0;
                            else
                                if i-radiusx <= 0
                                    negxdim = i-1;
                                else
                                    negxdim = radiusx;
                                end
                                if i+radiusx > lenx
                                    posxdim = lenx-i-1;
                                else
                                    posxdim = radiusx;
                                end
                                if j-radiusy <= 0
                                    negydim = j-1;
                                else
                                    negydim = radiusy; 
                                end
                                if j+radiusy > leny
                                    posydim = leny-j-1;
                                else
                                    posydim = radiusy;
                                end

                                try
                                    grid = reference(j-negydim:j+posydim, i-negxdim:i+posxdim);
                                catch
                                   disp('Computation failed');
                                   return;
                                end
                                diff = grid - measured(j,i);

                                [y,x] = find(abs(diff) < 0.0001);
                                % Equal values not found -> find the any
                                % neihbouring values which are higher and lower
                                % than the reference
                                novalues = 1; 

                                if ~isempty(y) && ~isempty(x)
                                    dta(j,i) = sqrt(min((abs(y-negydim)*interval(1)).^2 + (abs(x-negxdim)*interval(2)).^2));
                                    if dta(j,i) > 2
                                       dta(j,i) = 2; 
                                    end
                                    if dta(j,i) > dtacrit || isnan(dtaA(i))
                                        failedcountDTAA = failedcountDTAA + 1;
                                    end
                                else
                                    % Check the values around the reference
                                    [ysmaller, xsmaller] = find(diff < 0);
                                    [ylarger, xlarger] = find(diff > 0);
                                    if isempty(ysmaller) || isempty(ylarger)
                                        novalues = 1;
                                    else
                                        novalues  = 0;
                                        % Find the distance where the value
                                        % changes sign
                                        tempvals = [];
                                        for k = 1:length(ysmaller)
                                            for l = 1:length(ylarger)
                                                % If both y and x coordinates
                                                % are next to each other ->
                                                % this is the sport where the
                                                % difference surface changes
                                                % sigm
                                                if abs(ysmaller(k) - ylarger(l)) <= 1 && abs(xsmaller(k) - xlarger(l)) <= 1
                                                    % compute dta (max error = resolution/2, use 'analystical' solution or interpolation for more accurate results)
                                                    tempvals(end+1) = (abs(ysmaller(k)-negydim)*interval(1)).^2 + (abs(xsmaller(k)-negxdim)*interval(2)).^2;

                                                end
                                            end
                                        end
                                        dta(j,i) = sqrt(min(tempvals));
                                    end
                                end

                                if novalues == 1
                                    dta(j,i) = 2;
                                end

                                if dta(j,i) > dtacrit
                                   failedcountDTA = failedcountDTA + 1; 
                                end
                            end
                        else
                            dta(j,i) = NaN;
                        end
                end
        end
        
        matsize = sum(sum(gmask));
        DTApercent2 = (matsize - failedcountDTA)/matsize*100;
    end
    % Analytical    
    if dtamethod(1) == 1
        % Work in progress (?)
        disp('Analytical 2D DTA computation does not exist yet. DTA computed using brute force method');
    end
    
    %% Dose difference
    [gidxy, gidxx] = find(gmask);
    [gidxydel, gidxxdel] = find(gmask == 0);
    indexdel = sub2ind([size(gmask,1), size(gmask,2)], gidxydel, gidxxdel);
    
    index2 = sub2ind([size(gmask,1), size(gmask,2)], gidxy, gidxx);
    Drelative = zeros(size(gmask,1),size(gmask,2));
    Drelative(index2) = (measured(index2) - reference(index2))./reference(index2);    %./(reference*dosecrit/100);
    if isempty(matsize)
        matsize = sum(sum(gmask));
    end
    
    % Set indices to NaN, which do not exceed the threshold
    Drelative(indexdel) = NaN;
    
    
    Dabsolute = (measured - reference);
    Dabsolute(indexdel) = NaN;
    [Dpasspercy, ~] = find(abs(Drelative)*100 < dosecrit*2);
    
    % Number of passes
    Dpasssize = length(Dpasspercy);
    Dpassperc = Dpasssize/matsize*100;
end





