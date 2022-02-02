

% original imagfe path
path2 = '/data/bm05/inhouse/thomas/180625lenses/det_dist/';%'/data/bm05/inhouse/thomas/180501lenses/det_dist';
% export saving folder
saveFolder = '/data/bm05/inhouse/thomas/180625lenses/processing/disto/';%'/data/bm05/inhouse/thomas/180501lenses/processing/distoTest/';

%% %%%%%%%%%%%%%%%%%%%%%%%%%dont touch after this %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FlatfieldImagePath =   [];%'/data/bm05/inhouse/thomas/180423lenses/det_dist5/flat';%[];%'/data/bm05/inhouse/thomas/171211lenses/processing/detdisto';%xst/xst_meshdet1/';% [empty] for no darkfield
FlatfieldImageName1 = 'avim';%%'distoFlat.edf';%'img0001.edf';%'avim';%'distoFlat.edf';%'ipp10.TIF';%[empty] for no file
 


savePlaceName   = [saveFolder 'distograd.mat'];
saveplace = [saveFolder 'disto_det.mat'];


orientationscan = 'vertical';   % 'vertical'or 'horizontal' => you will have to run this script twice for complete characterization distortion
%DetectorDistortion
orientationscan = 'horizontal'; % 'vertical'or 'horizontal' => you will have to run this script twice for complete characterization distortion
DetectorDistortion




IntegDistoMaps