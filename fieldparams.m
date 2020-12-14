function [caxdev, sym1, hom, dmax1, dmin1, dev, fieldsize, penR, penL] = fieldparams(mydata, mypos, datatype, caxcorr)
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
%% Computes the field parameters sym, hom, dmax, dmin, pen, grad, dev
% Input data must be CAX normalized

% if ~isempty(varargin)
% 
%     caxcorr = varargin{2};
% end

% Fieldsize = 50% CAX (IEC)

% Find closest values to 50% CAX until closest values are found
% symmetrically (Finds closest values to 50% from each side)


% Find the last position with 0.5 value

lastvalfound = false;

% Find the maximum value
[maxval, idxmax] = max(mydata);


% Find the first position with 0.5 value
[~, idx1] = min(abs(mydata - 0.5*maxval));
moddata = mydata;
moddata(idx1) = [];
% Find the 'right hand side' 0.5*max value
while lastvalfound == false
    [~, idx2] = min(abs(moddata - 0.5*maxval));
    if (idx2 > idxmax && idxmax > idx1) || (idx2 < idxmax && idxmax < idx1)
        lastvalfound = true;
    end
    moddata(idx2) = [];
end


% Use the mid axis of the first and 2nd 0.5 values as evaluated CAX value
% (For profiles the penumbras are relatively symmetrical, so the symmetrical axis between 0.5*CAX and 0.5*MAX will in most cases be small)
offset = round(abs((idx2 - idx1))/2);
if idx1 < idx2
    caxevalpos = idx1 + offset;
else
    caxevalpos = idx1 - offset;
end

CAXeval = mydata(caxevalpos);


lastvalfound = false;
% Repeat the fieldsize idx search for 0.5*CAX values
% Find the first position with 0.5 value
[~, idx1] = min(abs(mydata - 0.5*CAXeval));
moddata = mydata;
moddata(idx1) = [];
% Find the 'right hand side' 0.5*max value
while lastvalfound == false
    [~, idx2] = min(abs(moddata - 0.5*CAXeval));
    if (idx2 > caxevalpos && caxevalpos > idx1) || (idx2 < caxevalpos && caxevalpos < idx1)
        lastvalfound = true;
    end
    moddata(idx2) = [];
end



% Next work only if the data is around 0-point
% if mypos(idx1) < 0
%     idx2 = idx1;
%     while mypos(idx2) <= 0
%         [~, idx2] = min(abs(moddata - 0.5));
%         moddata(idx2) = [];
%     end
% else
%     idx2 = idx1;
%     while mypos(idx2) >= 0
%         [~, idx2] = min(abs(moddata - 0.5));
%         moddata(idx2) = [];
%     end
% end


fieldsize = abs(mypos(idx2) - mypos(idx1));

%fieldsize = abs(mypos(idx1)) + abs(mypos(idx2));

% CAX deviation

% CAX deviation calculated by using the fieldsize

% if mypos(idx2) > 0 && mypos(idx1) > 0
%     caxdev = min(mypos(idx2) - mypos(idx1));
% else
%     caxdev = min(mypos(idx2) + mypos(idx1));
% end





% Find index for the wanted field region (0.8*field size, size >= 10cm, 0.6*field size when size < 10cm)

% if measurement == 'omni'
%     if fieldsize >= 10.0
%         [val, idxFSpos] = min(abs(mypos - fieldsize*0.4));
%         [val, idxFSneg] = min(abs(mypos - -fieldsize*0.4));
%     else
%         [val, idxFSpos] = min(abs(mypos - fieldsize*0.6));
%         [val, idxFSneg] = min(abs(mypos - fieldsize*0.6));
%     end
% else

% Flattened region according to the IEC definition of symmetry


if fieldsize <= 10.0
    [~, idxFSpos] = min(abs(mypos - (mypos(caxevalpos) + (fieldsize/2 - 0.5))));
    [~, idxFSneg] = min(abs(mypos - (mypos(caxevalpos) - (fieldsize/2 - 0.5))));
    if fieldsize < 5.0
       disp('Field too small! Fieldsize calculated using the IEC values for field sizes between 5.0 and 10.0 cm. Results may not be accurate'); 
    else
        disp('Fieldsize OK')
    end
elseif fieldsize <= 30
    [~, idxFSpos] = min(abs(mypos - (mypos(caxevalpos) + (fieldsize/2 - fieldsize*0.1))));
    [~, idxFSneg] = min(abs(mypos - (mypos(caxevalpos) - (fieldsize/2 - fieldsize*0.1))));
    disp('Fieldsize OK')
else
    [~, idxFSpos] = min(abs(mypos - (mypos(caxevalpos) + (fieldsize/2 - 1.5))));
    [~, idxFSneg] = min(abs(mypos - (mypos(caxevalpos) - (fieldsize/2 - 1.5))));
    disp('Fieldsize OK')
end



moddata = mydata;
[~, idx20pos] = min(abs(mydata - 0.20*CAXeval));
[~, idx80pos] = min(abs(mydata - 0.80*CAXeval));

try
    % PenumbraR (20%/80%) (IEC). Same method as for 50% CAX
    while idx20pos <= caxevalpos || idx80pos <= caxevalpos
        if idx20pos <= caxevalpos
            moddata(idx20pos) = [];
            [~, idx20pos] = min(abs(moddata - 0.2*CAXeval));
        end
        if idx80pos <= caxevalpos
            moddata(idx80pos) = [];
            [~, idx80pos] = min(abs(moddata - 0.8*CAXeval));
        end
    end
    penR = abs(mypos(idx20pos) - mypos(idx80pos));
catch
   disp('Error in penumbra calculation! CAX deviation may be too high');
   penR = 0;
end

moddata = mydata;
% PenumbraL (20%/80%) (IEC)
[~, idx20neg] = min(abs(moddata - 0.2*CAXeval));
[~, idx80neg] = min(abs(moddata - 0.8*CAXeval));

try
    while idx20neg >= caxevalpos || idx80neg >= caxevalpos
       if idx20neg >= caxevalpos
            moddata(idx20neg) = [];
            [~, idx20neg] = min(abs(moddata - 0.2*CAXeval));
        end
        if idx80neg >= caxevalpos
            moddata(idx80neg) = [];
            [~, idx80neg] = min(abs(moddata - 0.8*CAXeval));
        end
    end
    penL = abs(mypos(idx80neg) - mypos(idx20neg));
catch
   disp('Error in penumbra calculation! CAX deviation may be too high');
   penL = 0;
end


caxdev = round(min(mypos(idx2), mypos(idx1)) + fieldsize/2,4);

if caxcorr == 1
    mypos = mypos + caxdev;
    % Renormalize
    [~, caxidx] = min(abs(mypos));
    mydata = mydata/mydata(caxidx);
end


% Symmetry [sym = (d(x)/d(-x))max (IEC)

try
if mod(length(mypos),2) == 0
    % Odd data
    countuntil = idxFSneg + (idxFSpos-idxFSneg)/2;
    [val, idxZero] = min(abs(mypos));
    if val < 0
       caxval =  (mydata(idxZero) + mydata(idxZero+1))./2;
    else
       caxval =  (mydata(idxZero) + mydata(idxZero-1))./2;
    end
else 
    % Even data
    countuntil = idxFSneg + (idxFSpos-idxFSneg)/2 - 1;
    [~, idxZero] = min(abs(mypos));
    caxval = mydata(idxZero);
end
catch
   disp('Error in computing symmetry.');
   disp(lasterror);
end

dmax1 = max(mydata(idxFSneg:idxFSpos))*100;
dmin1 = min(mydata(idxFSneg:idxFSpos))*100;

% Dev (IEC)
try
    cmp(1) = abs(dmin1/100 - caxval)/caxval;
    cmp(2) = abs(dmax1/100 - caxval)/caxval;
    dev = max(cmp)*100;
catch
   disp('Error in computing deviation'); 
end

ratioarr = zeros(length(idxFSneg:countuntil),1);
for i = idxFSneg:countuntil
   val1 = mydata(i);
   flipped = flip(mydata);
   val2 = flipped(i);
   cmp(1) = abs(val1/val2);
   cmp(2) = abs(val2/val1);
   
   ratioarr(i,1) = max(cmp);
end

% Final symmetry value (IEC)
sym1 = max(ratioarr)*100;

% Homogeneity or flatness (IEC)
hom = dmax1/dmin1*100;


% Test plot on a new figure
% figure(3);
% hold on
% plot(mypos, mydata);
% scatter(mypos(idxmax), maxval);
% scatter(mypos(idx1), mydata(idx1), 'g');
% scatter(mypos(idx2), mydata(idx2), 'k');
% scatter(mypos(caxevalpos), mydata(caxevalpos), 'c');



