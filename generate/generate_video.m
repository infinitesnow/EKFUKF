function generate_video(filename,movie)
if(~exist('generated_video','dir'))
    mkdir('./generated_video')
end

v = VideoWriter(strcat('./generated_video/',filename,'.mp4'),'MPEG-4');
v.Quality = 95;
v.FrameRate = 24;
open(v);
writeVideo(v,movie);
close(v);

end