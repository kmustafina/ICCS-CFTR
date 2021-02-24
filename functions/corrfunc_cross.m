function [G] = corrfunc(imgser1, imgser2);

% Kamila Mustafina (based on David Kolin's coffunc.m)
 
G = zeros(size(imgser1));

set(gcbf,'pointer','watch');

   h = waitbar(0,'Calculating 2D cross-correlation functions...');
   
for z=1:size(imgser1,3)
    G(:,:,z) = ((fftshift(real(ifft2(fft2(double(imgser1(:,:,z))).*conj(fft2(double(imgser2(:,:,z))))))))/(mean(mean(imgser1(:,:,z)))*mean(mean(imgser2(:,:,z)))*size(imgser1,1)*size(imgser1,2))) - 1;
    if ishandle(h)
     waitbar(z/size(imgser1,3),h)
     else
         break
     end
end

if ishandle(h)
 close(h)
end

set(gcbf,'pointer','arrow');
%}
