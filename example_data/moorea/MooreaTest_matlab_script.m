%           This program was written By Craig Nelson
%           Program was written 9 May 2009 11am - 9 pm
%           Revised July 2011 for Jim Sickman on behalf of Ted Swift
%           Revised Aug 2014 for Murage Ngatia DWR
%           Revised Jul 2015 for Stu Goldberg CMORE
%           Revised Jan 2016 for Zach Quinlan CMORE
%           Revised Sep 2016 for Zach Quinlan CMORE
%           Revised Oct 2016 for Zach Quinlan/BobWhittier CMORE
%           Revised Oct 2017 for Zach Quinlan (Added Lignin for Linda)
%           Revised and uploaded Nov 2016 by Zach Quinlan for future version control
%           Read in csv files from from Fluoromax/Aqualog / FluorEssence software
%           Collate the files into a database
%           Normalize each sample to the Raman Peak of water (average of two DIW samples selected by user)
%           Subtract averaged DIW blank samples defined by user
%           Produce a subset database excluding DIW blank samples 
%           Generate a Summary Table
%           Allow easy reference and graphing of samples in database
%           Adjustable for variations between runs in slit widths or run parameters
%           Generates a data structure suitable for use with DOMFluor (Stedmon & Bro)


%           VARIABLES USED

%          Data sets remaining after run:

%          MainData is the data structure analagous to "OriginalData" in DOMFluor
%          SummaryDataMatrix is a spreadsheet of all key values at the end of the
%          analysis, this is exported to excel format as well

%          EEMainTrimmed gives the final Raman calibrated and blank-subtracted data for all samples including blanks with Rayleigh scatter cut out
%          EEMainTrimmedRawBlanks gives same as EEMainTrimmed except blanks have not had blanks subtracted
%          EEMainSubset gives the same as EEMainTrimmed but without DIW blanks and without removal of Rayleigh scatter - designed for export to DOMFluor
%          EEMainSubsetTrimmed gives the same as EEMainTrimmed but without DIW blanks
%          EEMainSubsetTrimmedStacked and Relative turn the 3-D matrix into a 2-D matrix with each row an Ex-Em peak for ordination
%          OrdinationRowMatrix and Relative are exported to excel and are designed for use in clustering and ordination downstream applications
%          EEMTrimmedDatabase is the final database with all row and column labels which is then exported to large spreadsheet

%          Em is just the values of the Emission wavelengths
%          Ex is the values of the Excitation wavelengths

%          samplepaths and shortsamplenames give the list of the final path and file names for the subset of data delivered to the EEMainSubset and MainData sets

%           To graph any given sample, use the following command, replacing # with a
%           sample number (row numbers in your input spreadsheet) and/or replacing
%           EEMainTrimmed with Blank1Trimmed or Blank2Trimmed. If you use EEMainSubset use
%           the row numbers in the variable shortsamplenames.

%           contourf(Ex,Em,squeeze(EEMainTrimmed(#,:,:)));colorbar;

%           A fancier graph can be made with a similar command...
%           surfc(Ex,Em,squeeze(EEMainTrimmed(#,:,:)));shading interp;  colorbar;

%          blanks should be graphed on a Raman Unit scale of (0-0.01 or 0.02) using the command
%          contourf(Ex,Em,squeeze(EEMainTrimmedwithBlanks(#,:,:)),(0:0.002:0.02));colorbar;

%          Other variables used (Should be comprehensive)

%          ExStep gives the Excitation steps in nm
%          Emstep gives the Emission steps in nm
%          Slitwidth gives the slit widths used in nm
%           BlankAverage gives the average Blank subtracted from each EEM, Raman calibrated
%           clear w; %           dummy
%           clear s; %           dummy    
%           clear m; %           dummy
%           clear temp; %          dummy
%           clear Range; %          used to define the range of csv files from which dlmread grabs data   
%           clear temp2; %          dummy
%           clear temp3; %          dummy
%           clear Blank1, Blank2; %          trimmed blanks
%           clear BlankDatabase1, BlankDatabase2, Database; %          are all extended reads from the original csv files with all values
%           clear EEMlist; the imported cell variable containing the spreadsheet of samples and paths
%           clear EEMain; original data files including blanks trimmed
%           clear BlankRaman; raman calibrated blanks
%           clear EEMainRaman;raman calibrated samples
%           clear EEMainBlkSubt; Blank-subtracted raman calibrated samples
%           clear EEMfilename; %          input file of the excel sheet containing filenamesand paths to csv files
%           clear SummaryDataNumsCell; %          dummy to structure final summary table
%           clear SummaryDataNums; %          dummy to structure final summary table
%           clear SummaryDataRows; %          dummy to structure final summary table
%           clear SummaryDataColumns; %          dummy to structure final summary table
%           clear ExEm; %           dummy to structure ordination tables
%           clear OrdinationColumns; %          dummy to structure ordination tables
%           clear nSampleOrig; %          total number of samples including blanks
%           clear N; %          dummy
%           clear S; %          dummy
%           clear n; %          dummy
%           clear i; %          dummy
%           clear RamanArea; Raman peak areas Ex: 350nm from 371 to 428 nm emission
%           clear ExEnd; user input end of excitation reads
%           clear ExStart; user input start of Excitation reads
%           clear EmStart; user input start of emission reads
%           clear EmEnd; user input end of excitation reads
%           clear EmArray; calculated values of emission reads
%           clear ExArray; calculated values of excitation reads from 
%           clear EEMainSubsetTrimmedStacked and StackedRelative are designed for use in clustering and ordination analyses and list all Ex-Em pair peaks in a list
%           clear Blank1Trimmed and Blank2Trimmed give the two unique blanks averaged and subtracted from each EEM, Raman calibrated


clear; % Clears the workspace of all datasets and variables
close all; % closes all open files

%           User entered information about the dataset

%          Note that tab-delimited sample listing file should have 6 columns as follows 
%           1) A unique Sample Name
%           2) the word Blank or Sample
%           3) full path to sample file
%           4) full path to first blank
%           5) full path to second blank
%           6) full path to absorbance scan

%          In this area you will input all of your run parameters for the set
%           in the future you can hard-code these values.

EEMfilename='InputFileVariri.tsv'; % input('What is the name of the Sample listing tab-delimited text file (in single quotes, including extension)?\n?: ');
disp('Note: column 2 should list text Blank or Sample, column 1 should have unique informative names of samples,');
disp('      columns 3, 4, and 5 should list complete path names for each sample and respective blanks,');
disp('      and column 6 lists complete path names for the absorbance scan.');

Slitwidth=5;
ExStart=240;
ExEnd=500;
ExStep=5;
EmStart=248; % these values are designed to estimate the range of emissions matching the number of Aqualog export rows...
EmEnd=824.6; %...assuming there is a EmStep distance between each row or each emission value.
EmStep=4.65; %On the aqualog this is dictated by the number of pixels assigned to integrate 4.65nm = 8 pixels

%          This section calculates the array size and is specific to format of delimited text files: CSVs or DATs or TXTs 
%ExArray=ExEnd-ExStart;%           these develop integer steps for Excit/Emiss vals
%EmArray=EmEnd-EmStart;
Ex=ExEnd:-ExStep:ExStart; %          these develop axis values for Ex/Em, example 500 495 490 ... 245 240
Em=EmStart:EmStep:EmEnd;
Range=[1 1 length(Em) length(Ex)]; %          this selects the appropriate range of the scans to be read in by dlmread
                                    % Typical Aqualog scans for Nelson Lab are 53 columns (Ex) and 125 rows (Em)
%          Generate Matlab Matrices from Excel Workbook
%          ********* T Swift could not find "xlsread" in the Octave manual *********

fileID=fopen(EEMfilename); EEMlist=textscan(fileID,'%s %s %s %s %s %s'); fclose(fileID); % reads in the list of samples
Rowcount=size(EEMlist{1},1);
UniqueIDs=EEMlist{1}(2:Rowcount,1);
SampleType=EEMlist{2}(2:Rowcount,1);
SampleEEMFile=EEMlist{3}(2:Rowcount,1);
SampleBlank1File=EEMlist{4}(2:Rowcount,1);
SampleBlank2File=EEMlist{5}(2:Rowcount,1);
SampleAbsFile=EEMlist{6}(2:Rowcount,1);

Samplecount=Rowcount-1; % starts at 2 to skip the column headings
for i=1:Samplecount 
    Database(i,:,:)=dlmread(char(SampleEEMFile(i,1)),'\t',Range); %          inputs each sample EEM into a master database
    BlankDatabase1(i,:,:)=dlmread(char(SampleBlank1File(i,1)),'\t',Range); %          inputs two blanks for each sample into a blank database
    BlankDatabase2(i,:,:)=dlmread(char(SampleBlank2File(i,1)),'\t',Range);
    Blank1(i,:,:)=BlankDatabase1(i,:,1:1:size(BlankDatabase1,3));%          places relevant region of the worksheet into variable Blank
    Blank2(i,:,:)=BlankDatabase2(i,:,1:1:size(BlankDatabase2,3)); %          places relevant region of the worksheet into variable Blank
    BlankAverage(i,:,:)=(Blank1(i,:,:)+Blank2(i,:,:))/2; %          averages 2 blanks for each sample
    EEMain(i,:,:)=Database(i,:,1:1:size(Database,3)); % places relevant region of the worksheet into variable EEMain, can be modified for Varian exports to account for the fact that every other column is a list of emissions
end;
clear i;

% - Set up Absorbance Databases - should find a way to initialize at 0 for compatibility with datasets without Absorbance scans
Abs=(ExEnd:-ExStep:ExStart); % Absorbance should be integer steps already, so no need for analogous ABSArray as above for Ex/Em
ABSRange=[0 1 length(Abs)-1 1]; % This selects the correct area of the absorbance scans to input by dlmread
for i=2:Samplecount
   if strcmp(char(SampleType(i,1)),'Blank')==0 % if the sample is NOT a blank...
    ABSDatabase(i,:)=dlmread(char(SampleAbsFile(i,1)),'\t',ABSRange); % inputs each sample Absorbance scan into a separate database
    for a=1:length(Em)
        for b=1:length(Ex)
            ABSMain(i,a,b)=10^(0.5*(ABSDatabase(i,find(Abs<Em(a),1,'first')))+ABSDatabase(i,find(Abs==Ex(b)))); 
            % assigns the IFE correction factor of antilog mean Aex + Aem for each Ex-Em pair of the EEM following Kothawala et al 2013 L&O Methods
        end;
    end;
   end;
end;
   

%           Apply Inner Filter Correction (first), Raman scaling (second) and Blank subtraction (third) to each EEM subset by date

nSampleOrig=size(EEMain,1); 
for i=1:nSampleOrig
    RamanArea(i)=EmStep*sum(BlankAverage(i,find(Em<=380,1,'last'):find(Em>=426,1,'first'),find(Ex==350,1,'last'))); % sums the area in the averaged blanks under the Raman Scatter peak according to Kothawala et al 2013
    EEMainIFE(i,:,:)=times(EEMain(i,:,:),ABSMain(i,:,:));% Inner Filter Effect correction using Absorbance Scan Factors built above
    EEMainRaman(i,:,:)=EEMainIFE(i,:,:)/RamanArea(1,i); %          divides each EEM by the respective Raman area
    BlankRaman(i,:,:)=BlankAverage(i,:,:)/RamanArea(1,i); %          divides each averaged blank by the respective Raman area
    EEMainBlkSubt(i,:,:)=EEMainRaman(i,:,:)-BlankRaman(i,:,:); %          subtracts the averaged blank from each sample
end
clear i;

EEMainSubset=EEMainBlkSubt; %          makes a dummy variable to prep for removing blanks from the dataset
SampleTypetemp=SampleType;
samplepaths=SampleEEMFile;
shortsamplenames=UniqueIDs;
FinalSampleIDs=UniqueIDs;
S='Blank: Removed from final dataset';
i=1;
m=1;
while i<=size(EEMainSubset,1) %          strips the blanks out of the dataset but keeps track of their new positions
   if strcmp(char(SampleTypetemp(i,1)),'Blank')==1 % if the sample is a blank...
        EEMainSubset(i,:,:)=[]; %          removes blanks from EEMain to make new array EEMainSubset to input to DOMFluor
        samplepaths(i)=[]; %          removes blank file paths
        shortsamplenames(i)=[]; %          removes blank IDs 
        SampleTypetemp(i)=[]; %          removes blanks from EEMlist for tracking purposes
        FinalSampleIDs(m)=strcat(UniqueIDs(m,1),':',cellstr(S)); %          makes complete list of sample names with blanks tagged as removed
        NewSamplePosition(m)=0; %          makes a complete list of new sample positions in the EEMainSubset with zeros for blanks
        m=m+1; %          counter for values that are blanks or samples
    else
        NewSamplePosition(m)=i; %          makes a complete list of new sample positions in the EEMainSubset with zeros for blanks
        i=i+1; %          counter only for values that are samples
        m=m+1; %          counter for values that are blanks or samples
    end;
end;
clear i;
clear m;
clear S;

%          Remove the Rayleigh Scattering where Ex>=Em: NOTE can also be done in
%          PARAFAC using the new DOMFluor toolbox. Here we do this to allow for
%          summary statistics on the total DOM fluorescence.

EEMainTrimmed=EEMainBlkSubt; %          sets up a working version of the blank subtracted data to trim
EEMainTrimmedwithBlanks=EEMainRaman; %          sets up a working version of the Raman calibrated data to trim (blanks not subtracted)
for i=1:nSampleOrig
    Blank1Trimmed(i,:,:)=Blank1(i,:,:)/RamanArea(1,i); %          sets up a working version of the Blank1
    Blank2Trimmed(i,:,:)=Blank2(i,:,:)/RamanArea(1,i); %          sets up a working version of the Blank2
end;
EEMainTrimmed(EEMainTrimmed<0)=0;                                             %          Converts negatives to zeros
EEMainTrimmedwithBlanks(EEMainTrimmedwithBlanks<0)=0;
Blank1Trimmed(Blank1Trimmed<0)=0;  
Blank2Trimmed(Blank2Trimmed<0)=0;  
for i=1:round(length(Ex)) %          basically this section removes the rayleigh scatter by cutting away based on slit widths
    w=Em(1)-Ex(i)-(2*Slitwidth);
    s=w*-1/EmStep;
    if w < 0
        EEMainTrimmed(:,1:s,i)=0;
        EEMainTrimmedwithBlanks(:,1:s,i)=0;
        Blank1Trimmed(:,1:s,i)=0;
        Blank2Trimmed(:,1:s,i)=0;
    end;
end;

%       Export a complete Trimmed Excel database of EEMs: all samples are raman calibrated, 
%       all samples have average blanks subtracted except Blanks which are just raw, raman-calibrated data.
%   
warning off MATLAB:xlswrite:AddSheet;
for i=1:size(EEMainTrimmed,1)
    if strcmp(char(SampleType(i,1)),'Blank')==1
        EEMainTrimmedRawBlanks(i,:,:)=EEMainTrimmedwithBlanks(i,:,:); %          if the sample was a DIW blank insert data without blanks subtracted
    else
        EEMainTrimmedRawBlanks(i,:,:)=EEMainTrimmed(i,:,:); %          if the sample was a sample EEM then use the blank-subtracted data
    end;
    midarray=num2cell(cat(2,Em',squeeze(EEMainTrimmedRawBlanks(i,:,:))));
    toparrayrow=[SampleEEMFile(i,1) num2cell(Ex)];
    EEMTrimmedDatabase(i,:,:)=cat(1,toparrayrow,midarray);
    xlswrite('EEMDatabase.xlsx',squeeze(EEMTrimmedDatabase(i,:,:)),char(UniqueIDs(i,1))); %          Slows down processing and can be removed if a database is not necessary
    %csvwrite('EEMDatabase.xlsx',squeeze(EEMTrimmedDatabase(i,:,:)),char(UniqueIDs(i,1)));
    if i==size(EEMainTrimmed,1) disp('Raman-calibrated and blank subtracted samples and Raman-calibrated DIW blanks written to EEMDatabase.xls as sheets'); else end;
end;

%   Calculate indices and peak maxima
[c index] = min(abs(Em-410));NearestEm410=Em(index); %for M:C index and CobleM
[c index] = min(abs(Em-445));NearestEm445=Em(index); %for M:C index and CobleC

[c index] = min(abs(Em-380));NearestEm380=Em(index); %for BIX index
[c index] = min(abs(Em-430));NearestEm430=Em(index); %for BIX index

[c index] = min(abs(Em-434));NearestEm434=Em(index);% for HIX index
[c index] = min(abs(Em-480));NearestEm480=Em(index);% for HIX index
[c index] = min(abs(Em-300));NearestEm300=Em(index);% for HIX index
[c index] = min(abs(Em-346));NearestEm346=Em(index);% for HIX index

[c index] = min(abs(Em-470));NearestEm470=Em(index);% for FI
[c index] = min(abs(Em-520));NearestEm520=Em(index);% for FI

[c index] = min(abs(Em-450));NearestEm450=Em(index);% for CobleA
[c index] = min(abs(Em-305));NearestEm305=Em(index);% for CobleB
[c index] = min(abs(Em-340));NearestEm340=Em(index);% for CobleT
[c index] = min(abs(Em-299));NearestEm299=Em(index);% for Fpeak
[c index] = min(abs(Em-509));NearestEm509=Em(index);% for Stedmon D
[c index] = min(abs(Em-435));NearestEm435=Em(index);% for Optical Brighteners
[c index] = min(abs(Em-510));NearestEm510=Em(index);% for dieselBandII
[c index] = min(abs(Em-510));NearestEm510=Em(index);% for Petroleum
[c index] = min(abs(Em-360));NearestEm360=Em(index);% for Lignin

for i=1:size(EEMainTrimmed,1)
    CobleA(1,i)=EEMainTrimmed(i,find(Em==NearestEm450,1),find(Ex==320,1)); % picks Excite 320 Emission 450 (Coble Peak A)
    CobleB(1,i)=EEMainTrimmed(i,find(Em==NearestEm305,1),find(Ex==275,1)); % picks Excite 275 Emission 305 (Coble Peak B)
    CobleT(1,i)=EEMainTrimmed(i,find(Em==NearestEm340,1),find(Ex==275,1)); % picks Excite 275 Emission 340 (Coble Peak T)
    CobleM(1,i)=EEMainTrimmed(i,find(Em==NearestEm410,1),find(Ex==310,1)); % picks Excite 310 Emission 410 (Coble Peak M)
    CobleC(1,i)=EEMainTrimmed(i,find(Em==NearestEm445,1),find(Ex==345,1)); % picks Excite 345 Emission 445 (Coble Peak C)
    Fpeak(1,i)=EEMainTrimmed(i,find(Em==NearestEm299,1),find(Ex==240,1)); % picks Excite 240 Emission 299 (F peak phenylalanine)
    StedmonD(1,i)=EEMainTrimmed(i,find(Em==NearestEm509,1),find(Ex==390,1)); % picks Excite 390 Emission 509 (Stedmon peak D)
    OpticalBrighteners(1,i)=EEMainTrimmed(i,find(Em==NearestEm435,1),find(Ex==360,1)); % picks Excite 435 Emission 360
    DieselBandII(1,i)=EEMainTrimmed(i,find(Em==NearestEm510,1),find(Ex==410,1));%picks Excite 410, Emission 550
    Petroleum(1,i)=EEMainTrimmed(i,find(Em==NearestEm510,1),find(Ex==270,1));%Picks Excite 270, Emission 500
    Lignin(1,i)=EEMainTrimmed(i,find(Em==NearestEm360,1),find(Ex==240,1));% Picks Excite 240, Emission 360 (Lignin)
    MtoC(1,i)=CobleM(1,i)/CobleC(1,i);
    BIX(1,i)=EEMainTrimmed(i,find(Em==NearestEm380,1),find(Ex==310,1))/EEMainTrimmed(i,find(Em==NearestEm430,1),find(Ex==310,1));
    HIX(1,i)=sum(EEMainTrimmed(i,find(Em==NearestEm434,1):find(Em==NearestEm480,1),find(Ex==255,1)))/sum(EEMainTrimmed(i,find(Em==NearestEm300,1):find(Em==NearestEm346,1),find(Ex==255,1)));
    FI(1,i)=EEMainTrimmed(i,find(Em==NearestEm470,1),find(Ex==370,1))/EEMainTrimmed(i,find(Em==NearestEm520,1),find(Ex==370,1));
end;

%          Prep Data structure for use in Stedmon and Bro 2008 DOM FLuor walkthrough
MainData=struct('Ex',Ex','Em',Em','X',EEMainSubset,'nEx',length(Ex),'nEm',length(Em),'nSample',size(EEMainSubset,1),'XBackup',EEMainSubset);
disp('Final data structure for use in DOMFluor exported to structure array MainData');
% 
% %   Start DOMFluor
%DOMFLuor % Checks that the package is ready to start
%MainData % Lists the contents of the data structure to be input to DOMFluor
[CutData]=EEMCut(MainData,40,0,15,0,''); % Trims the dataset
%PlotSurfby4(1,CutData,'R.U.'); % Plots sequential plots of the EEMs
[Test1]=OutlierTest(CutData,2,2,6,'Yes','No'); % Examines the structure of PARAFAC models up to 6 models using steps 2nm emission and 5nm excite.
%PlotLL(Test1,2);
%PlotLL(Test1,3);
%PlotLL(Test1,4);
%PlotLL(Test1,5);
%PlotLL(Test1,6);
%EvalModelSurf(Test1,6);
[AnalysisData]=SplitData(Test1);
[AnalysisData]=SplitHalfAnalysis(AnalysisData,(4:6),'SplitHalf.mat');
SplitHalfValidation(AnalysisData,'1-2',4);
SplitHalfValidation(AnalysisData,'3-4',4);
SplitHalfValidation(AnalysisData,'1-2',5);
SplitHalfValidation(AnalysisData,'3-4',5);
SplitHalfValidation(AnalysisData,'1-2',6);
SplitHalfValidation(AnalysisData,'3-4',6);
[AnalysisData]=RandInitAnal(AnalysisData,4,10);
[AnalysisData]=RandInitAnal(AnalysisData,5,10);
[AnalysisData]=RandInitAnal(AnalysisData,6,10);
TCC(AnalysisData.Model4,AnalysisData.Split(1).Fac_4);
TCC(AnalysisData.Model4,AnalysisData.Split(2).Fac_4);
TCC(AnalysisData.Model4,AnalysisData.Split(3).Fac_4);
TCC(AnalysisData.Model4,AnalysisData.Split(4).Fac_4);
TCC(AnalysisData.Model5,AnalysisData.Split(1).Fac_5);
TCC(AnalysisData.Model5,AnalysisData.Split(2).Fac_5);
TCC(AnalysisData.Model5,AnalysisData.Split(3).Fac_5);
TCC(AnalysisData.Model5,AnalysisData.Split(4).Fac_5);
TCC(AnalysisData.Model6,AnalysisData.Split(1).Fac_6);
TCC(AnalysisData.Model6,AnalysisData.Split(2).Fac_6);
TCC(AnalysisData.Model6,AnalysisData.Split(3).Fac_6);
TCC(AnalysisData.Model6,AnalysisData.Split(4).Fac_6);
%PlotLL(AnalysisData,4);
%EvalModelSurf(AnalysisData,4);
[PFMax4,PFEmLoad4,PFExLoad4]=ModelOut(AnalysisData,4,'ParafacResults4Component.xlsx');
[PFMax5,PFEmLoad5,PFExLoad5]=ModelOut(AnalysisData,5,'ParafacResults5Component.xlsx');
[PFMax6,PFEmLoad6,PFExLoad6]=ModelOut(AnalysisData,6,'ParafacResults6Component.xlsx');

%          Summarize and export EEM index and PARAFAC Data
SummaryDataRows=cellstr(char('Raman peak area of blanks','Total blank1 RU','Total Blank2 RU','Total sample RU after blank subtraction','Position in final database','M:C','BIX','HIX','FI','Ultra Violet Humic-like','Marine Humic-like','Visible Humic-like','Tryptophan-like','Tyrosine-like','Phenylalanine-like','Fulvic Acid-like','Optical Brighteners','Diesel Band II','Petroleum-like','Lignin-like','PARAFAC1','PARAFAC2','PARAFAC3','PARAFAC4','PARAFAC5','PARAFAC6'));
SummaryDataColumns=cat(2,cellstr(char('Sample Name')),FinalSampleIDs');
for i=1:size(EEMainTrimmed,1)
    SummaryDataNums(1,i)=RamanArea(:,i);
    SummaryDataNums(2,i)=sum(sum(Blank1Trimmed(i,:,:)));
    SummaryDataNums(3,i)=sum(sum(Blank2Trimmed(i,:,:)));
    SummaryDataNums(4,i)=sum(sum(EEMainTrimmed(i,:,:)));
    SummaryDataNums(5,i)=NewSamplePosition(1,i);
    SummaryDataNums(6,i)=MtoC(1,i);
    SummaryDataNums(7,i)=BIX(1,i);
    SummaryDataNums(8,i)=HIX(1,i);
    SummaryDataNums(9,i)=FI(1,i);
    SummaryDataNums(10,i)=CobleA(1,i);
    SummaryDataNums(11,i)=CobleM(1,i);
    SummaryDataNums(12,i)=CobleC(1,i);
    SummaryDataNums(13,i)=CobleT(1,i);
    SummaryDataNums(14,i)=CobleB(1,i);
    SummaryDataNums(15,i)=Fpeak(1,i);
    SummaryDataNums(16,i)=StedmonD(1,i);
    SummaryDataNums(17,i)=OpticalBrighteners(1,i);
    SummaryDataNums(18,i)=DieselBandII(1,i);
    SummaryDataNums(19,i)=Petroleum(1,i);
    SummaryDataNums(20,i)=Lignin(1,i);
    if NewSamplePosition(1,i)>0
        SummaryDataNums(21,i)=PFMax6(NewSamplePosition(1,i),1);
        SummaryDataNums(22,i)=PFMax6(NewSamplePosition(1,i),2);
        SummaryDataNums(23,i)=PFMax6(NewSamplePosition(1,i),3);
        SummaryDataNums(24,i)=PFMax6(NewSamplePosition(1,i),4);
        SummaryDataNums(25,i)=PFMax6(NewSamplePosition(1,i),5);
        SummaryDataNums(26,i)=PFMax6(NewSamplePosition(1,i),6);
    else
        SummaryDataNums(21,i)=0;
        SummaryDataNums(22,i)=0;
        SummaryDataNums(23,i)=0;
        SummaryDataNums(24,i)=0;
        SummaryDataNums(25,i)=0;
        SummaryDataNums(26,i)=0;

    end;
 end;
SummaryDataNumsCell=num2cell(SummaryDataNums);
SummaryDataMatrix=cat(2,SummaryDataRows,SummaryDataNumsCell);
SummaryDataMatrix=cat(1,SummaryDataColumns,SummaryDataMatrix);
SummaryDataMatrix=SummaryDataMatrix';
xlswrite('SummaryDataMatrix.xlsx', SummaryDataMatrix); %
%csvwrite('SummaryDataMatrix.csv', SummaryDataMatrix); 
%disp ('Information on peaks, indices and sample numbers exported to SummaryDataMatrix.csv');

%          Prep data for use in Nelson Hierarchical Analyses and ordination

EEMainSubsetTrimmed=CutData.X; % Use DOMFluor to generate a trimmed matrix
%EEMainSubsetTrimmed=EEMainSubset; %          generate a reduced subset with rayleigh scatter cutaway
EEMainSubsetTrimmed(EEMainSubsetTrimmed<0)=0;  %          Converts negatives to zeros
EEMainSubsetTrimmed(isnan(EEMainSubsetTrimmed)) = 0 ;
clear i;
% for i=1:round(length(Ex)) %          basically this section removes the rayleigh scatter by cutting away based on slit widths
%     w=Em(1)-Ex(i)-(2*Slitwidth);
%     s=w*-1/EmStep;
%     if w < 0
%         EEMainSubsetTrimmed(:,1:s,i)=0;
%     end;
% end;

for i=1:size(EEMainSubsetTrimmed,1) %          stack all data so that each column is a sample with all ExEm as rows
    EEMainSubsetTrimmedStacked(:,i)=reshape(squeeze(EEMainSubsetTrimmed(i,:,:)),length(Ex)*length(Em),1);
    EEMainSubsetTrimmedStackedRelative(:,i)=EEMainSubsetTrimmedStacked(:,i)/sum(EEMainSubsetTrimmedStacked(:,i));
end;

for i=1:round(length(Em)) %          generate a row label list with excitation and emission values
    for m=1:length(Ex)
    ExEm(i,m)=cellstr(['Ex' int2str(Ex(1,m)) 'Em' int2str(Em(1,i))]);
    end;
end;
ExEm=reshape(ExEm,length(Ex)*length(Em),1); %          generates a list of Ex-Em peak names in format Ex250Em320
OrdinationColumns=cat(1,cellstr('ExEmPeak'),shortsamplenames); %          generate a column label list
OrdinationRowMatrix=cat(2,ExEm,num2cell(EEMainSubsetTrimmedStacked)); %           attach row labels to main matrices
OrdinationRowMatrixRelative=cat(2,ExEm,num2cell(EEMainSubsetTrimmedStackedRelative));
xlswrite('OrdinationMain.xlsx',cat(1,OrdinationColumns',OrdinationRowMatrix)); %          export complete matrices to excel
disp ('2-dimensional data suitable for clustering/ordination exported to OrdinationMain.xlsx');
xlswrite('OrdinationMainRelative.xlsx',cat(1,OrdinationColumns',OrdinationRowMatrixRelative));
disp('2-dimensional data relativized to total peak area for clustering/ordination exported to OrdinationMainRelative.xlsx');

%          Do some quick distance ordination and clustering algorithms with the data
RelEuclidean=pdist(EEMainSubsetTrimmedStackedRelative','euclidean');
RawEuclidean=pdist(EEMainSubsetTrimmedStacked','euclidean');
RelCityBlock=pdist(EEMainSubsetTrimmedStackedRelative','cityblock');
RawCityBlock=pdist(EEMainSubsetTrimmedStacked','cityblock');

subplot;
subplot(2,4,1); 
    dendrogram(linkage(RelEuclidean,'average'),'orientation','left','labels',shortsamplenames, 'colorthreshold', 'default');
    title('Relative Euclidean Distance');
subplot(2,4,2); 
    dendrogram(linkage(RawEuclidean,'average'),'orientation','left','labels',shortsamplenames, 'colorthreshold', 'default');
    title('Raw Euclidean Distance');
subplot(2,4,3); 
    dendrogram(linkage(RelCityBlock,'average'),'orientation','left','labels',shortsamplenames, 'colorthreshold', 'default');
    title('Relative City Block Distance');
subplot(2,4,4); 
    dendrogram(linkage(RawCityBlock,'average'),'orientation','left','labels',shortsamplenames, 'colorthreshold', 'default');
    title('Raw City Block Distance');

[RelMDSEuclid, stresseuclid]=mdscale(RelEuclidean,2);
subplot(2,4,5:6); scatter(RelMDSEuclid(:,1),RelMDSEuclid(:,2)); title('MDS of Relative Euclidean Distance');
[RelMDSCityBlock, stresscityblock]=mdscale(RelCityBlock,2);
subplot(2,4,7:8);scatter(RelMDSCityBlock(:,1),RelMDSCityBlock(:,2)); title('MDS of Relative City Block Distance');
disp('Please click on points in the MDS which you wish to label then press a key');
gname(shortsamplenames);



%          Clean up clutter by removing unused variables
clear w;
clear s;    
clear m;
clear temp;
clear Range;    
clear temp2;
clear temp3;
clear BlankAverage;
clear Blank1;
clear Blank2;
clear BlankDatabase1;
clear BlankDatabase2;
clear Database;
clear EEMlist;
clear EEMain;
clear BlankRaman;
clear EEMainRaman;
clear EEMainBlkSubt;
clear EEMfilename;
clear SummaryDataNumsCell;
clear SummaryDataNums;
clear SummaryDataRows;
clear SummaryDataColumns;
clear nSampleOrig;
clear n;
clear i;
clear RamanArea;
clear Slitwidth;
clear ExEnd;
clear ExStart;
clear EmStart;
clear EmEnd;
clear EmStep;
clear ExStep;
clear EmArray;
clear ExArray;
clear ExEm;
clear OrdinationColumns;
clear Blank1Trimmed;
clear Blank2Trimmed;
clear EEMainTrimmedwithBlanks;
clear midarray;
clear toparrayrow;
clear ABSMain;
clear ABSRange;
clear Abs;
clear a;
clear b;






