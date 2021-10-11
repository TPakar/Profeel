function [] = writeOmnifile(filename, varargin)

%% This function exports 1D and 2D data to omnipro-export form
fileID = fopen(filename,'w');

% VARARGIN
% {filename, imagename, length unit, separator, unit, plane, planeposition ,rowsize,
% colsize, energy, datatype, datafactor, number of bodies (set to 1), data, colpos, rowpos}

%%
% Print header info
disp('writing Omnifile');

prtstring = "<opimrtascii>";
fprintf(fileID, '%s\n\n', prtstring);

prtstring = "<asciiheader>";
fprintf(fileID, '%s\n', prtstring);

separator = varargin{4};
if strcmp(separator, '[TAB]')
   separator = "\t"; 
end

headernames = {"File Name:","Image Name:", "Length Unit:", "Separator:","Data Unit:", "Plane:", "Plane Position:", "No. of Rows:","No. of Columns:" ...
    ,"Energy:","Data Type:","Data Factor:", "Number of Bodies:"};
for i = 1:length(headernames)
   if i == 1
       fprintf(fileID, '%s', "File Version:");
       for j = 1:(20-length(char("File Version:")))
           % Fill in spaces to 20 as in omni export
           fprintf(fileID, '%s', " ");
       end
       fprintf(fileID, '%d\n', 3);
   end
    
   if ~isempty(varargin{i}) && ~contains(headernames{i}, "Plane Position")
       prtstring = headernames{i};
       % Add data headline
       
       fprintf(fileID, '%s', prtstring);
       for j = 1:(20-length(char(headernames{i})))
           % Fill in spaces to 20 as in omni export
           fprintf(fileID, '%s', " ");
       end
       % Add the header value
       if isnumeric(varargin{i})
          fprintf(fileID, '%g\n', varargin{i});
       else
          if strcmp(prtstring, 'Separator:')
            fprintf(fileID, '"%s"\n', varargin{i});
          elseif strcmp(prtstring, 'Data Type:') && sum(strcmp(varargin{i}, {'Absolute','Abs. Dose','Relative','Rel. Dose', 'Fluence'})) < 1
              fprintf(fileID, '%s\n', 'Relative');
          elseif strcmp(prtstring, 'Data Unit:') && sum(strcmp(varargin{i}, {'mGy','cGy', 'rad','Gy','1/10%'})) < 1
              fprintf(fileID, '%s\n', '1/10%');
          elseif strcmp(prtstring, 'Energy:')
              templist = {'MV', 'KV', 'MeV'};
              idx = [];
              for k = 1:length(templist)
                 if contains(lower(varargin{i}), lower(templist{k}))
                     idx = k;
                 end
              end
              if ~isempty(idx)
                  fprintf(fileID, '%g ', str2double(varargin{i}(1:end-length(templist{idx}))));
                  fprintf(fileID, '%s\n', templist{idx});
              else
                  fprintf(fileID, '%d ', 1);
                  fprintf(fileID, '%s\n', 'MeV');
                  disp('Energy unknown -> set to 1 MeV');
              end
          else    
            fprintf(fileID, '%s\n', varargin{i});
          end
       end
   end
end
% End of header
prtstring = "</asciiheader>";
fprintf(fileID, '%s\n\n\n', prtstring);
prtstring = "<asciibody>";
fprintf(fileID, '%s\n', prtstring);
if ~isempty(varargin{7})
    prtstring = append(headernames{7},"     ", varargin{7}(1:end-2), " ", varargin{7}(end-1:end));
    fprintf(fileID, '%s\n\n', prtstring);
end

% Column positions
if strcmp(varargin{6},'YZ')
    if strcmp(separator, '\t')
        prtstringcol = ['Y[', varargin{3},']    '];
        prtstringrow = ['Z[', varargin{3},'] '];
    else
        prtstringcol = ['Y[', varargin{3},']  ', separator];
        prtstringrow = ['Z[', varargin{3},'] '];
    end
elseif strcmp(varargin{6},'XZ')
    if strcmp(separator, '\t')
        prtstringcol = ['X[', varargin{3},']    '];
        prtstringrow = ['Z[', varargin{3},'] '];
    else
        prtstringcol = ['X[', varargin{3},']   ', separator];
        prtstringrow = ['Z[', varargin{3},'] '];
    end
else
    if strcmp(separator, '\t')
        prtstringcol = ['X[', varargin{3},']    '];
        prtstringrow = ['Y[', varargin{3},'] '];
    else
        prtstringcol = ['X[', varargin{3},']   ', separator];
        prtstringrow = ['Y[', varargin{3},'] '];
    end
end

fprintf(fileID, '%s%c', prtstringcol);
if strcmp(separator, '\t')
    fprintf(fileID, '\t');
end


for i = 1:varargin{9}
    if strcmp(separator, '\t')
        fprintf(fileID, '%.4g \t', varargin{15}(i));
    else
        fprintf(fileID, '%.4g %s', varargin{15}(i), separator);
    end
end
% Row positions
fprintf(fileID, '\n%s%c', prtstringrow);
%if strcmp(separator, '\t')
%    fprintf(fileID, '\t');
%end
fprintf(fileID, '\n');

for i = 1:length(varargin{16})
    if strcmp(separator, '\t')
        fprintf(fileID, '%.4g \t', varargin{16}(i,1));
    else
        fprintf(fileID, '%.4g %s', varargin{16}(i,1), separator);
        testarr(i,1) = varargin{16}(i,1);
    end
    
    for j = 1:size(varargin{14},2)
        if strcmp(separator, '\t')
            fprintf(fileID, '%.4g \t', varargin{14}(i,j));
        else
            testarr(i,j+1) = varargin{14}(i,j);
            fprintf(fileID, '%.4g %s', varargin{14}(i,j), separator);
        end
    end
    fprintf(fileID, '\n');
end
%save('testarr', 'testarr');

prtstring = "</asciibody>";
fprintf(fileID, '%s', prtstring);

fprintf(fileID, '\n');

prtstring = "</opimrtascii>";
fprintf(fileID, '%s', prtstring);

disp('File ready!');
