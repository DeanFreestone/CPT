% this function is the short time IF histogram

% inputs
% ~~~~
% IF = vector or matrix of IFs from CPT EMD
% IA = vector or matrix of IAs from CPT EMD
% Fs = Sampling rate (Hz)
% NBins = Number of bins for the histogram
% f_0 = Start frequency
% f_1 = Max frequency
% WindowSize = Time duration to bins the data in seconds
% Overlap = Percentage overlap in for consectutive windows

function [t_out Freq_Vector IF_image IA_IF_image] = STIFH(IF,IA,Fs,NBins,f_0,f_1,WindowSize,Overlap)

[NSamples NIMFs] = size(IF);    
if NSamples < NIMFs                 % make sure the IF and IA are the right way around
    IF = IF';
    IA = IA';
    temp = NSamples;
    NSamples = NIMFs;
    NIMFs = temp;
end

IA_IF = IA.*IF;
Freq_Vector = linspace(f_0,f_1,NBins);          % this is how the data will be binned

WindowSamples = round(WindowSize*Fs);
WindowIncrement = round((1-Overlap/100)*WindowSamples);
NWindows = floor(  NSamples / WindowIncrement  );

t_out = (WindowSamples:WindowIncrement:NSamples)/Fs;

IF_image = zeros(NBins,NWindows);               % initilise for speed
IA_IF_image = zeros(NBins,NWindows);
n_IF = zeros(NBins,NIMFs);
n_IA_IF = zeros(NBins,NIMFs);
m=1;
for n=WindowSamples+1:WindowIncrement:NSamples
    IFWindow = IF(n-WindowSamples:n,:);                 % take a window of data
    IA_IFWindow = IA_IF(n-WindowSamples:n,:);
    for nn=1:NIMFs
        n_IF(:,nn) = hist(IFWindow(:,nn),Freq_Vector);
        n_IA_IF(:,nn) = hist(IA_IFWindow(:,nn),Freq_Vector);        
    end
    IF_image(:,m) = sum(n_IF,2);
    IA_IF_image(:,m) = sum(n_IA_IF,2);
    m=m+1;
end

