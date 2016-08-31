function [circlesVec]=findPlates(Img, board)
%% [circlesVec]=findPlates(Img, board)
% -------------------------------------------------------------------------
% Purpose: To find the coordinates and the radius of each plate.
%
% Description: The location of the plates is more or less known, and the
%       constants XsRange, YsRange and RRange represent the approximate
%       location. Using findCirc, the exact location is found for each
%       plate, and this is the return value of the function.
%
% Arguments: Img - A grayscale image, taken by the scanner at a resolution
%       of 300 dpi.
%       board - (1) Delerin, (2) without holes,  (3) with holes
%
% Returns: CirclesVec - A 6x3 matrix. Each row is a plate, 
%       the columns are x,y of the centre and the radius.
% -------------------------------------------------------------------------
% Irit Levin. 03.2007
% -------------------------------------------------------------------------
% updated: Irit Levin 04.09 - boads definitions in sparate files.
%% constants
%CEN = [717, 651; 717, 1774; 717, 2896; 1902, 652; 1899, 1775; 1898, 2896];
disp([datestr(now)   '   Find Plates']);
% switch board
%     case 1
%         % Delerin Scanning cell
%         XsRange = [ 710: 725;...
%                     710: 725;...
%                     710: 725;...
%                    1890:1905;...
%                    1890:1905;...
%                    1890:1905];
% 
%         YsRange = [ 643: 658;...
%                    1766:1781;...
%                    2888:2903;...
%                     643: 658;...
%                    1767:1782;...
%                    2888:2903];
%                
%         RRange  = [508:520];
%         
%     case 2
%         % Aluminum Scanning cell without holes
%         XsRange = [ 656: 676;...
%                     658: 678;...
%                     657: 677;...
%                    1895:1915;...
%                    1898:1918;...
%                    1901:1921];
% 
% %         YsRange = [ 600: 620;...
% %                    1759:1779;...
% %                    2916:2936;...
% %                     598: 618;...
% %                    1756:1776;...
% %                    2914:2934];
%         YsRange = [ 592: 612;...
%                    1749:1769;...
%                    2906:2926;...
%                     588: 608;...
%                    1746:1766;...
%                    2905:2925];
% 
%         %RRange  = [508:520];
%         RRange  = [508:528];
%         
%     case 3
%         % Aluminum Scanning cell with holes
%         XsRange = [ 673: 693;...
%                     673: 693;...
%                     672: 692;...
%                    1932:1952;...
%                    1931:1951;...
%                    1932:1952];
% 
% %         YsRange = [ 583: 603;...
% %                    1743:1763;...
% %                    2901:2921;...
% %                     585: 605;...
% %                    1740:1760;...
% %                    2901:2921];
%           YsRange = [ 575: 595;...
%                      1733:1753;...
%                      2892:2912;...
%                       575: 595;...
%                      1732:1752;...
%                      2891:2911];
%         %RRange  = [512:528];
%         RRange  = [508:528];
%         
%     case 4
%         % New Delerin With pores
%         XsRange = [ 682: 702;...
%                     685: 705;...
%                     688: 708;...
%                    1946:1966;...
%                    1948:1968;...
%                    1949:1969];
% 
%           YsRange = [ 592: 612;...
%                      1774:1794;...
%                      2958:2978;...
%                       591: 611;...
%                      1774:1794;...
%                      2955:2975];
%                  
%         RRange  = [508:528];
%         
%     otherwise
%         error('choose: 1-delerin, 2-no holes, 3-withholes');
% end
%        
[XsRange, YsRange, RRange]=getPlatesDefinition(board);

%% finding the x,y,r of each plate
  figure; imshow(Img);%
  hold on%
progressBar = waitbar(0,'please wait...');

for i=1:6
    msg = sprintf('finding plate %d/6', i);
    waitbar((i-1)/6, progressBar , msg);
    
    [res maxRs]=findCirc(RRange,XsRange(i,:),YsRange(i,:),Img,0.2,1);
    
    circlesVec(i,:) = [res(1,1), res(1,2), res(1,3)];
      circle([circlesVec(i,1),circlesVec(i,2)],circlesVec(i,3) ,500,'r-');%
      circle([circlesVec(i,1),circlesVec(i,2)],440 ,500,'y-');%
end
close(progressBar);