%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gerchberg and Saxton Algorithm 
%author - Hisay Lama
%email ID - hisaylama@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,  clc,
close all,
%Input beam 
M = 800; %side length (pixel unit)
N = 600; % side length (pixel unit)
psize = 20e-6; %pxiel size (m)
x=([0.5:1:M-0.5] - M/2);
y =([0.5:1:N-0.5] - N/2);
[X,Y]=meshgrid(x,y);
w=250; % Radius of the incident beam (px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input beam projected on diffracted optical element
%u_in = circ((sqrt((X-M/4).^2 + (Y-N/4).^2))./w);
u_in = circ((sqrt(X.^2 + Y.^2))/w); %Ampliture of the input beam
I_in = abs(u_in.^2);
figure(1), %irradiance image
imagesc(x,y,I_in);
xlabel('x (m)'); ylabel('y (m)');
colormap('gray');
%axis square;
axis xy;
%I = double(I_in);
PH = 2.*pi.*(rand([N,M])); %Generation of the random phase
%u_in = u_in.*exp(-1i.*PH); %Resultant input beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the target image with only amplitude i.e. without phase
S_L = input('Enter the number of spots =  ');
P = input('Enter the inter-spot separation in pixel =  ');
w = 1; %radius of the spot (in pixel unit) 
u_target = Multiple_Spot(M,N,S_L,P);
x=([0.5:1:M-0.5] - M/2);
y =([0.5:1:N-0.5] - N/2);
[X,Y]=meshgrid(x,y);
I_target = abs(u_target.^2);
figure,
imagesc(x,y,I_target);
%axis off;
axis xy;
colormap('gray'); xlabel('x (px)'); ylabel('y (px)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting of the GS algorithm iteration for geenrating phase mask
f = input('focal length = '); %focal length in pixel units
z = 0; %distance from the focus
lambda = 1064e-9; %input('lambda = '); %light wavelength [m]
k = 2.*pi/lambda; %wave vector
[fx,fy] = meshgrid((1/psize).*[-M/2:1:M/2-1]/M,(1/psize)*[-N/2:1:N/2-1]/N);%frequency size
x = lambda*f*fx;
y = lambda*f*fy;
Phase_hot = PH;
E_inv = u_in;
for n = 1:200 %number of iteration
    Efocus = exp(1i*k*(2*f+z + x.^2 + y.^2))/(1i*lambda*f).*fftshift(fft2(abs(E_inv).*exp(-1j*Phase_hot).*exp(-1i*pi*z/(lambda*f^2)*(X.^2+Y.^2))));
    phase = mod(angle(Efocus), 2*pi);     
    E_in = abs(u_target).*exp(-1j*(phase));
    E_inv = ifftshift(ifft2(ifftshift(E_in)));
    Phase_hot = mod(angle(E_inv), 2*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the phase mask
%phase_z = 2*pi*(x.^2 + y.^2).*z/(lambda*f^2);
Phase = mod(Phase_hot, 2*pi);
figure,
imshow(mat2gray((Phase)./(2*pi)))%+ 1.2*pi
%title('Phase Plot')
axis equal tight off
box off
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the resulaant image after introducing the 
R = fftshift(ifft2(fftshift(exp(-1j*(angle(E_inv))))));
figure; imshow(mat2gray(abs(R)));
colorbar()
%imshow(mat2gray(abs(I1)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 