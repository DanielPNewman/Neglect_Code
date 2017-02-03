
%--Two level Folder Hierarchy checker--
%This script allows you to check the file path naming hierarchy in
%experiment data which were saved following a naming standard
%for instance, Direction Dots experiment data and BlueLight Data
%Both have a pattern in their names according to different
%experiment conditions and subject number,
%for example  DD1M/DD1M1.egg filepath refers to Direction Dots (DD),
%block 1 of subject 1 from Monash university(M) eeg data.
%This pattern can be written as DD#M/DD#M$, where # and $ are substituted
%by numbers, # for subjects and $ for blocks (two level hierarchy)
%This script can check up to 2 levels of hierarchy

%Example of two level hierarchy file folders


%DirectionDots_data
%  | - DD1M
%       | - DD1M1.egg
%       | - DD1M2.eeg
%       | - DD1M3.eeg
%       ...
%       | - DD1M$.eeg (where $ is the block number)
%  | - DD2M
%       | - DD2M1.egg
%       | - DD2M2.eeg
%       | - DD2M3.eeg
%       ...
%       | - DD2MM$.eeg (where $ is the block number)

%...

%  | - DD#M (where # is the subject number)
%       | - DD#M1.egg
%       | - DD#M2.eeg
%       | - DD#M3.eeg
%       ...
%       | - DD#M$.eeg (where $ is the block number)

%NOTE: For one level hierarchy search, just set the
%variable Total Folders Level1 to 1 and take care when
%specifying datafolder variable
%If for example, it is necessary just to evaluate
%block files (second level) in direction dots experiment data
%in the DD1M folder, datafolder variable will be 
%'G:\DirectionDots_data\DD1M\', pattern will be DD1M$.eeg
%and TotalFoldersLevel1 will be 1.


%Cleaning routine
clear
close all
clc

%----SET PARAMETERS HERE----

%Insert here the main folder
datafolder = 'S:\R-MNHS-SPP\Bellgrove-data\11.Megan_ONeill\EEG_Behav_Data Files\Healthy_Older\';

%Insert here the folder path pattern as a string
%for example DD#M/DD#M$.eeg (direction dots experiment)
%this means that the script will search files with the file path pattern
%substituting all "#" (subjects in the exampple)
%and "$" (Blocks in the example) symbols for numbers
%"#" for FIRST hierarchy level and "$" for SECOND hierarchy level.
pattern='DD#M\DD#M$.eeg';


%Insert there the number of folders of the first level of hierarchy
TotalFoldersLevel1=15;

%Insert there the number of folders of the second level of hierarchy
TotalFoldersLevel2=12;


%----MAIN SCRIPT----

%Search code
for j=1:TotalFoldersLevel1 %Search first level
   
   fprintf('Folder: %i\n',j) ;
   for i=1:TotalFoldersLevel2 %Search second level
       
       %Start and reset path 
       path=[];
       for k=1:length(pattern)
            
            if(pattern(k)~='#' &&  pattern(k)~='$') %keep letters in the path

                path=[path pattern(k)];

            elseif(pattern(k)=='#')

                path=[path num2str(j)]; %Substitute # for the first level folder number

            else

                path=[path num2str(i)]; %Substitute $ for the second level folder number

            end
            
            
        end
       
       if(exist([datafolder path],'file')) %Evaluate file existence
          disp([path ' FOUND!']);
       else
          disp([path ' --NOT FOUND !!--']);
       end
       
       
   end
   fprintf('\n\n');
      
end

fprintf('The "--NOT FOUND !!--" files may have one of the following problems:\n1) The folder is missing a file\n2) The file is not following the name pattern\n -> Wrong: DD1M.egg, Right:DD1M1.egg\n3) Double extension in files, for example:\n -> Wrong: DD1M1eeg.egg, Right:DD1M1.egg\nFor any of the problems above or any other\nyou should look the folder where are the files and figure out any mistake.\n');