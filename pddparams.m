function [R100, R80, R50, D100, D200, surfdose, J1020] = pddparams(pddData, pos)
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
% Computes PDD parameters to PDD data according to IEC

% Check that the data is not in descending order
if pos(1) > pos(2)
   %pddData = flip(pddData);
   pos = flip(pos);
   disp('Data order changed to ascending for PDD parameter computation');
end

% R100 = depth of maximum dose
[maxdose, maxidx] = max(pddData);
R100 = pos(maxidx);

% R80 = depth of 80% dose (further from the surface than max point)
[~, idx80] = min(abs(pddData(maxidx:end) - 0.8*maxdose));
idx80 = idx80 + maxidx-1;
R80 = pos(idx80);

% R50 = depth of 50% dose (furher from the surface than the max point)
[~, idx50] = min(abs(pddData(maxidx:end) - 0.5*maxdose));
idx50 = idx50 + maxidx-1;
R50 = pos(idx50);

% D100 = percentage depth dose at 100mm in ratio to dose at the maximum
% Ds = surface dose (0.5mm IEC)
[~, surfdoseidx] = min(abs(pos - 0.05));
surfdose = pddData(surfdoseidx)/maxdose*100;
[~, depth100idx] = min(abs(pos - 10));
D100 = pddData(depth100idx)/maxdose*100;

% D200 = percentage depth dose at 200mm in ratio to dose at the maximum
[~, depth200idx] = min(abs(pos - 20));
D200 = pddData(depth200idx)/maxdose*100;

% J10/20 is a quality parameter defined as the ratio of D100/D200
J1020 = D100/D200;

% Testplot parameter positions
% figure(4);
% hold on
% plot(pos, pddData);
% scatter(pos(maxidx), pddData(maxidx),'r');
% scatter(pos(idx80), pddData(idx80),'g');
% scatter(pos(idx50), pddData(idx50),'b');
% scatter(pos(depth100idx), pddData(depth100idx),'k');
% scatter(pos(depth200idx), pddData(depth200idx), 'c');

hold off
% If cannot compute all the parameters, print an error message
if sum(isnan([R100, R80, R50, D100, D200, J1020])) > 1
   disp('Error in computing PDD parameters. The data does not seem to be a PDD nor a Profile. Keeping PDD parameters'); 
end