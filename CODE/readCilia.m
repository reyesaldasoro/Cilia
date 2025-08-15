function CiliaVolume = readCilia(currFile)

b                                       = imfinfo(currFile);
numSlices                               = numel(b);

% determine the number of channels, it may be that there are 2, 3 or 4
% channels depending on the acquisition. The ImageDescription field has
% data about the channels, use this to guide the loop
% Two channels
              % <Channel AcquisitionMode="WideField" Color="255" EmissionWavelength="450" ExcitationWavelength="405" ID="Channel:0" Name="DAPI" PinholeSize="22.988506">
              %    <DetectorSettings Binning="1x1" ID="Detector:0" />
              % </Channel>
              % <Channel AcquisitionMode="WideField" Color="65280" EmissionWavelength="525" ExcitationWavelength="488" ID="Channel:1" Name="EGFP" PinholeSize="22.988506">
              %    <DetectorSettings Binning="1x1" ID="Detector:0" />
              % </Channel>
% Three channels
             % <Channel AcquisitionMode="WideField" Color="255" EmissionWavelength="450" ExcitationWavelength="405" ID="Channel:0" Name="DAPI" PinholeSize="38.314175">
             %     <DetectorSettings Binning="1x1" ID="Detector:0" />
             %  </Channel>
             %  <Channel AcquisitionMode="WideField" Color="65280" EmissionWavelength="525" ExcitationWavelength="488" ID="Channel:1" Name="EGFP" PinholeSize="38.314175">
             %     <DetectorSettings Binning="1x1" ID="Detector:0" />
             %  </Channel>
             %  <Channel AcquisitionMode="WideField" Color="16711680" EmissionWavelength="595" ExcitationWavelength="561" ID="Channel:2" Name="Alexa Fluor 555 goat anti-mouse IgG antibody/pH 7.2" PinholeSize="38.314175">
             %     <DetectorSettings Binning="1x1" ID="Detector:0" />
             %  </Channel>
% Four channels
              % <Channel AcquisitionMode="WideField" Color="255" EmissionWavelength="450" ExcitationWavelength="405" ID="Channel:0" Name="DAPI" PinholeSize="38.314175">
              %    <DetectorSettings Binning="1x1" ID="Detector:0" />
              % </Channel>
              % <Channel AcquisitionMode="WideField" Color="65280" EmissionWavelength="525" ExcitationWavelength="488" ID="Channel:1" Name="EGFP" PinholeSize="38.314175">
              %    <DetectorSettings Binning="1x1" ID="Detector:0" />
              % </Channel>
              % <Channel AcquisitionMode="WideField" Color="16711680" EmissionWavelength="595" ExcitationWavelength="561" ID="Channel:2" Name="Alexa Fluor 555 goat anti-mouse IgG antibody/pH 7.2" PinholeSize="38.314175">
              %    <DetectorSettings Binning="1x1" ID="Detector:0" />
              % </Channel>
              % <Channel AcquisitionMode="LaserScanningConfocalMicroscopy" ContrastMethod="Brightfield" ID="Channel:3" Name="TD" PinholeSize="38.314175">
              %    <DetectorSettings Binning="1x1" ID="Detector:0" />
              % </Channel>


switch numel(strfind(b(1).ImageDescription,'ID="Channel'))
  case 2
      % Two channels, EGFP and DAPI
      CiliaVolume(b(1).Width,b(1).Height,3,numSlices/2) = 0;
      for k2=1:2:numSlices
          %tempImage           = double(imread(currFile,k2));
          CiliaVolume(:,:,1,1+((k2)-1)/2)     = 0;
          CiliaVolume(:,:,2,1+((k2)-1)/2)     = double(imread(currFile,k2+1)); % Cilia in the green
          CiliaVolume(:,:,3,1+((k2)-1)/2)     = double(imread(currFile,k2+0)); % DAPI in third to be blue
      end

  case 3
      % No Brightfield included, 3 channels
      CiliaVolume(b(1).Width,b(1).Height,3,numSlices/3) = 0;
      for k2=1:3:numSlices
          %tempImage           = double(imread(currFile,k2));
          CiliaVolume(:,:,1,1+((k2)-1)/3)     = double(imread(currFile,k2+2));  
          CiliaVolume(:,:,2,1+((k2)-1)/3)     = double(imread(currFile,k2+1));
          CiliaVolume(:,:,3,1+((k2)-1)/3)     = double(imread(currFile,k2+0));  % DAPI in third to be blue
      end
  case 4
      % Four channels 3 fluorescent plus brightfield
      CiliaVolume(b(1).Width,b(1).Height,3,numSlices/4) = 0;
      for k2=1:4:numSlices
          %tempImage           = double(imread(currFile,k2));
          CiliaVolume(:,:,1,1+((k2)-1)/4)     = double(imread(currFile,k2+2));
          CiliaVolume(:,:,2,1+((k2)-1)/4)     = double(imread(currFile,k2+1));
          CiliaVolume(:,:,3,1+((k2)-1)/4)     = double(imread(currFile,k2+0));  % DAPI in third to be blue
          CiliaVolume(:,:,4,1+((k2)-1)/4)     = double(imread(currFile,k2+3));
      end
  otherwise
      disp('Error: number of channels must be between 2 and 4')

end


% if isempty(strfind(b(1).ImageDescription,'ID="Channel:3'))
%     % No Brightfield included, 3 channels
%     CiliaVolume(b(1).Width,b(1).Height,3,numSlices/3) = 0;
%     for k2=1:3:numSlices
%         %tempImage           = double(imread(currFile,k2));
%         CiliaVolume(:,:,1,1+((k2)-1)/3)     = double(imread(currFile,k2+2));
%         CiliaVolume(:,:,2,1+((k2)-1)/3)     = double(imread(currFile,k2+1));
%         CiliaVolume(:,:,3,1+((k2)-1)/3)     = double(imread(currFile,k2+0));
%     end
% else
%     CiliaVolume(b(1).Width,b(1).Height,3,numSlices/4) = 0;
%     for k2=1:4:numSlices
%         %tempImage           = double(imread(currFile,k2));
%         CiliaVolume(:,:,1,1+((k2)-1)/4)     = double(imread(currFile,k2+2));
%         CiliaVolume(:,:,2,1+((k2)-1)/4)     = double(imread(currFile,k2+1));
%         CiliaVolume(:,:,3,1+((k2)-1)/4)     = double(imread(currFile,k2+0));
%         CiliaVolume(:,:,4,1+((k2)-1)/4)     = double(imread(currFile,k2+0));
%     end
% 
% end