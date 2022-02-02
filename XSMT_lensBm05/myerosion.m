function filtered = myerosion(imageIN, valthreshold,varargin)
% erosion filter to remove some of the impulse noise impulse noise

% mandatory field : imageIN = the image 
%                   valthreshold = the value above which one the pixel must
%                   be changed (for speckle correlation it is better to
%                   change the minimum of them
%
% optional field : a numerical inter for the size of the median filter
% (default is 5)
%                       the option 'silent'for the silent mode


filtsz =  [5 5];
mutemode = 0;


if ~isempty(varargin)
    for k = 1:length(varargin)
        if isnumeric(varargin{k})
            filtsz = [varargin{k} varargin{k}];       
        elseif  ischar(varargin{k}) && strcmp(varargin{k},'silent')
            mutemode = 1;
        end;
    end;
end;



if isempty(imageIN),        error('Empty image');       end;
if isempty(valthreshold),   error('Empty threshold');   end;

[m,n] = size(imageIN);

pts2 = sub2ind([m n],[1 1 1 1 2 2 (m-1) (m-1) m m m m],[1 2 (n-1) n 1 n 1 n 1 2 (n-1) n]);
val2 = imageIN(pts2);

%if filtsz(1) > 5
    im1med = medfilt2(imageIN,filtsz,'symmetric');   
% else 
%     im1med = cvmedfilt2(imageIN,filtsz(1));   
% end;
    
diffIm = abs(imageIN - im1med);
diffIm(pts2) = 0;

pts1 = find(diffIm > valthreshold);

imageIN(pts1) = im1med(pts1);
imageIN(pts2) = val2;

filtered = imageIN;
if ~mutemode
    disp(['Number of points corrected with myFilt = ' num2str(numel(pts1)) ' over ' num2str(numel(imageIN)) ' (' num2str(numel(pts1)/numel(imageIN)*100) ' %)']);
    disp('');
    fprintf('\n');
end;