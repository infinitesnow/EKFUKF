[file, sr] = audioread('440hz.wav');
base_freq = 440 %hz
w = get_pulse(base_freq,sr);
phase = 40;
scale = 1/25;
plot(0:1000,file(1:1001),0:1000,sin((phase:1000+phase)*w)*scale)

function w = get_pulse(freq,sr)
    freq = 1/sr*freq    % in samples/sec
    w = 2*pi*freq       % compute w
end