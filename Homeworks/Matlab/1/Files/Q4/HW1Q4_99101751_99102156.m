clear;
clc;

Fs = 44100 ;
filename = 'MySong.wav';

Notes = ["G", "G", "A#", "D#", "D", ...
 "G", "G", "A#", "D", "C", ...
 "G", "G", "G", "G", "G", "G#", ...
 "G#", "G#", "G#", "G#", "G", "G"];


NoteDurations = [330, 330, 490, 490, 790, ...
 330, 330, 490, 490, 750, ...
 330, 330, 330, 490, 490, 700, ...
 330, 330, 330, 490, 490, 750];

Song = songmaker(Notes,NoteDurations) ;
player = audioplayer(Song,Fs) ; % make audio from output of songmaker

play(player)  % play audio

audiowrite(filename,Song,Fs) ; % write audio
%% functions
function y = generateNote(freq,duration,alpha)
N = floor(44100/freq) ;
M = floor((duration*44100)/1000) ;
y = zeros(M,1) ;
x = rand(N,1).*2 -1 ;
y(1:N) = x ;
y(N+1) = alpha.*(y(1))./2 ;
for i = (N+2):M
    y(i) = alpha.*(y(i - N) + y(i - N - 1))./2 ;
end
y = transpose(y) ;
end
function freq = noteFreq(note)
if note == "C"
    freq = 261.63 ;
elseif note == "C#"
    freq = 277.18 ;
elseif note == "D"
    freq = 293.66 ;
elseif note == "D#"
    freq = 311.13 ;
elseif note == "E"
    freq = 329.63 ;
elseif note == "F"
    freq = 349.23 ;
elseif note == "F#"
    freq = 369.99 ;
elseif note == "G"
    freq = 392.00 ;
elseif note == "A"
    freq = 440 ;
elseif note == "G#"
    freq = 415.30 ;
elseif note == "A#"
    freq = 466.16 ;
elseif note == "B"
    freq = 493.88 ;
else
    freq = 0 ;
end
end
function song = songmaker(notes,noteDurations)
i = size(notes,2) ;
song = [] ;
for j = 1:i
    freq = noteFreq(notes(j)) ;
    song = [song,generateNote(freq,noteDurations(j),0.995)] ;
end
end