%TIME SCALING
%Ths function is all about time scaling(stretching) and implements two
%different implementations...one in the time domain and one in the
%frequency domain(then converted back to the time-domain). For the latter
%implementation, a second version is given...just downsampled (decimated to
%be exact). More information will be in the paper, but this switch/case
%model works pretty easily. Please read for the variable arguments...

%I obtained alot of information from the DAFX book...edited by Udo Zolzer
%(on the works cited page)

%I.E. call function: final_dst1('phasevo','Poly Original.wav',64,'output.wav', 512); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = final_dst1(method,filename,hops_signal,OUTfile,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%overall input checks
if (~ischar(filename)) || (~ischar(OUTfile)) || (~ischar(method))
    error('filename and OUTfile must be character arrays'); 
end
if (~isnumeric(hops_signal))
    error('method and hops_signal must be numerical data'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,fs,nBits] = wavread(filename); %read .wav file...shorter samples recommended
x = sum(x,2); %make the signal mono
length(x); %length, samples, of the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%default parameter...
N = 2048; %analysis block size...window the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method %switch between types of time scaling techniques...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sola is a time-domain based algorithm for time stretching
%sola input checks...variable arguments...read carefully
case 'sola'
if nargin > 6
    error('Too many input arguments for this method'); 
end
if nargin < 4
    error('Not enough input arguments');
end
if (isempty(varargin))
    tsalpha = .5; 
    overinterval = 256; 
    disp('tsalpha = .5 and overinterval = 256'); 
else
    tsalpha = varargin{1};
    disp('certain tsalphas may cause indices to become out of range');
    disp('usual parameter: .25 <= tsalpha <= 2'); 
    overinterval = varargin{2};
end
if hops_signal > N 
    error('N must be > than hop size'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hops = 256; %analysis hop size
noframes = ceil((length(x))/(hops_signal)); %number of frames...hops of the signal  

%zeropad original signal
zeropad = ((noframes)*hops_signal) + N; 
zerosig = zeros(zeropad,1); 
zerosig(1:length(x)) = x; 

hopshift = floor(hops_signal * tsalpha); %shift of hop size...incorporating time scaling
overshift = floor(overinterval*tsalpha/2); %length of overlap between shifted hop and time-scaled signal

if (hopshift > N) || (hopshift > N - overshift)
    error('hopshift must not exceed N or (N - overshift)'); 
end

y = zerosig(1:N,1); %copy first frame to time-scaled signal
max_val = zeros(length(noframes)  - 1,1); %value of max lag
max_idx = zeros(length(noframes)  - 1,1); %index of max lag

tic;
%segment the time-based signal into N-length windows...cross-correlate with
%scaled (part of) signal to find max lag within overlap interval-in order to find 
%point of best similarity...at point of max similarity the overlapping
%blocks/locations are weighted by fade-in and fade-out functions...then
%what's left is added to signal...and added on each iteration
for i = 1:noframes - 1
     windex = zerosig(i*hops_signal+1 : N+i*hops_signal); %index of frames...
     %starting with second
     segment = xcorr(windex(1:overshift),y(i*hopshift:i*hopshift+(overshift-1),1)); 
     [max_val(i,1),max_idx(i,1)] = max(segment); 
     fadeout = 1:(-1/(length(y) - (i*hopshift-(overshift-1)+max_idx(i,1) - 1))):0;
     fadein = 0:(1/(length(y) - (i*hopshift-(overshift-1)+max_idx(i,1) - 1))):1;
     fadetail = y((i*hopshift-(overshift-1))+max_idx(i,1)- 1:length(y),1).*fadeout'; 
     fadebegin = windex(1:length(fadein)).*fadein'; crossfade = fadetail + fadebegin; 
     %output vector
     y = [y(1:i*hopshift-overshift+max_idx(i,1) - 1,1); crossfade; windex(length(fadein)+1:N)]; 
end

%normalize between 1 and - 1
y = y / max(abs(y));

%sound(y,fs) %play the soundfile

%time plots...for orig. signal and time-scaled signal
timex = (0:1/fs:(length(x) / fs)); 
timey = (0:1/fs:(length(y) / fs)); 

%plotting functions
figure(1); 
subplot(2,1,1);
plot(timex(2:end),x);
xlabel('time'); ylabel('amplitude'); title(filename); 
subplot(2,1,2);
plot(timey(2:end),y);
xlabel('time'); ylabel('amplitude'); 
title({'timescaled'; filename}, 'Color','r','fontsize',16); 
toc; 
%write out a .wav time-scaled file
wavwrite(y,fs,nBits,OUTfile); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%phase vocoder-type of block-by-block approach (fft/ifft)-time scaling in
%frequency domain...time-frequency implementation with integer ratio
%phase vocoder input checks and variable arguments
case 'phasevo' 
if nargin > 5
    error('Too man input arguments for this method'); 
end
if nargin < 4
    error('Not enough input arguments');
end
if (isempty(varargin))
    outputinterval = 256; 
    disp('outputinterval = 256'); 
else
    outputinterval = varargin{1};
end
if hops_signal > N 
    error('N must be > than hop size'); 
end
%hops = 64;

noframes = ceil((length(x))/(hops_signal)); %number of frames...hops of the signal
    
w  = window(@blackman ,N); %apply window type 
stretch = outputinterval / hops_signal;  %determin time stretch ratio/factor
    
%zeropad the signal
zeropad = ((noframes)*hops_signal) + N; 
zerosig = zeros(zeropad,1); 
zerosig(1:length(x)) = x;
    
y = zeros(N + ceil(length(zerosig)*stretch), 1); %pre-allocated output vector
length(y); 

tic;  
%segment signal into N-length windows....apply window type....extract mag
%and phase from fft(shifted to center)...multiply to calculate the
%stretched signal(in the freq domain), then take the inverse fft to obtain
%a time-based version of the scaled signal(output)...then just keep overlap and adding to
%the vector on every iteration
counter = 0; 
for k = 1:noframes
    windex = zerosig((k-1)*hops_signal+1 :(k-1)*hops_signal + N);  %index of frames...
    windex = windex .* w; 
    FFT = fft(fftshift(windex)); 
    mag = abs(FFT); 
    phase = angle(FFT); 
    calculation = (mag .* exp(sqrt(-1)*stretch*phase)); 
    newsig = fftshift(real(ifft(calculation))) .* w;
    y(counter + 1:counter + N, 1) = y(counter + 1:counter + N, 1) + newsig; 
    length(y); 
    counter = counter + outputinterval; 
    if (hops_signal * (k-1)) > (length(zerosig) - N)
        break; 
    else 
        continue; 
    end
end

%normalize signal (from first window to end of output); 
y = y(N+1:length(y)) / max(abs(y));
length(y); 

%time plots...for orig. signal and time-scaled signal
timex = (0:1/fs:(length(x) / fs)); 
timey = (0:1/fs:(length(y) / fs)); 

%sound(y,fs); %play sound

%plotting functions
figure(1); 
subplot(2,1,1);
plot(timex(2:end),x); 
xlabel('time'); ylabel('amplitude'); title(filename); 
subplot(2,1,2);
plot(timey(2:end),y);
xlabel('time'); ylabel('amplitude'); 
title({'timescaled'; filename}, 'Color','r','fontsize',16);
toc; 
%write out a .wav time-scaled file
wavwrite(y,fs,nBits,OUTfile); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%downsampled version of above code...using decimation filter...all comments
%apply from above case unless otherwise noted
case 'phasevodown'
%decimate uses lowpass filter and downsamples/resamples...R = 2....1/R*fs
x = decimate(x,2,'FIR');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 5
    error('Too man input arguments for this method'); 
end
if nargin < 4
    error('Not enough input arguments');
end
if (isempty(varargin))
    outputinterval = 256; 
    disp('outputinterval = 256'); 
else
    outputinterval = varargin{1};
end
if hops_signal > N 
    error('N must be > than hop size'); 
end
%hops = 64;

noframes = ceil((length(x))/(hops_signal)); %number of frames...hops of the signal
    
w  = window(@hanning ,N);
stretch = outputinterval / hops_signal;  
    
zeropad = ((noframes)*hops_signal) + N; 
zerosig = zeros(zeropad,1); 
zerosig(1:length(x)) = x;
    
y = zeros(N + ceil(length(zerosig)*stretch), 1); 
length(y); 

tic;  
counter = 0; 
for k = 1:noframes
    windex = zerosig((k-1)*hops_signal+1 :(k-1)*hops_signal + N);  %index of frames...
    windex = windex .* w; 
    FFT = fft(fftshift(windex)); 
    mag = abs(FFT); 
    phase = angle(FFT); 
    ft = (mag .* exp(sqrt(-1)*stretch*phase)); 
    newsig = fftshift(real(ifft(ft))) .* w;
    y(counter + 1:counter + N, 1) = y(counter + 1:counter + N, 1) + newsig; 
    length(y); 
    counter = counter + outputinterval; 
    if (hops_signal * (k-1)) > (length(zerosig) - N)
        break; 
    else 
        continue; 
    end
end

y = y(N+1:length(y)) / max(abs(y));
length(y); 

%time plots...for orig. signal and time-scaled signal
timex = (0:1/fs:(length(x) / fs)); 
timey = (0:1/fs:(length(y) / fs));

%sound(y,fs/2); %play sound (at downsampled freq)

%plotting functions
figure(1); 
subplot(2,1,1);
plot(timex(2:end),x); 
xlabel('time'); ylabel('amplitude'); title(filename); 
subplot(2,1,2);
plot(timey(2:end),y);
xlabel('time'); ylabel('amplitude'); 
title({'timescaled'; filename}, 'Color','r','fontsize',16); 
toc; 
%write out a .wav time-scaled file...at downsampled rate
wavwrite(y,fs/2,nBits,OUTfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%////zeeshan