% Practica 5: Interpolacion

clear all
clc

% Carga de imagen
img = imread('pentagon256x256.tif');
img = im2double(img);
%img = rgb2gray(img);

% 1. Sobremuestreo espacial de la imagen con factores T = 2x2 y T = 4x4.

% Sobremuestreo 2x2
img_2x2 = sobremuestreo(img,2);

% Sobremuestreo 4x4 
img_4x4 = sobremuestreo(img,4);

% DFTs de la imagen original y las imagenes sobremuestreadas
DFT_original = fftshift(fft2(img));
DFT_2x2 = fftshift(fft2(img_2x2));
DFT_4x4 = fftshift(fft2(img_4x4));

figure('Name','Punto 1. Sobremuestreo de la imagen');
subplot(2,3,1), imagesc(img), title('Imagen original');
subplot(2,3,2), imagesc(img_2x2), title('Sobremuestreo con T = 2x2');
subplot(2,3,3), imagesc(img_4x4), title('Sobremuestreo con T = 4x4');
subplot(2,3,4), imagesc(log(1.0 + abs(DFT_original))), title('DFT imagen original');
subplot(2,3,5), imagesc(log(1.0 + abs(DFT_2x2))), title('DFT sobremuestreo con T = 2x2');
subplot(2,3,6), imagesc(log(1.0 + abs(DFT_4x4))), title('DFT sobremuestreo con T = 4x4');
colormap('gray');
%% Intercalando ceros DFT

z = sobremuestreo(DFT_2x2,2);
iz = ifft2( fftshift(z) );
Fz = [1 1];
Fz = Fz.' * Fz;
imgFz = imfilter(iz, Fz, 'conv');

figure('Name','DFT con ceros intercalados');
subplot(2,3,1), imagesc(img), title('Imagen original');
subplot(2,3,2), imagesc(img_2x2), title('Imagen sobremuestreo T = 2x2');
subplot(2,3,3), imagesc(log(1.0 + abs(DFT_2x2))), title('DFT sobremuestreo con T = 2x2');
subplot(2,3,4), imagesc(log(1.0 + abs(z))), title('DFT T = 2x2 con ceros intercalados');
subplot(2,3,5), imagesc(abs(iz)), title('Inversa de la DFT con los ceros intercalados');
subplot(2,3,6), imagesc(imgFz), title('Imagen interpolada DFT con los ceros intercalados');
colormap('gray');
%%
% 2. Interpolacion espacial de la imagen con factores T = 2x2 y T = 4x4.

% -------------------------------------------------------------------------
% Interpolacion de orden cero

% Filtro de orden cero con factor T = 2x2 
F02 = [1 1];
F02 = F02.' * F02;
imgF02 = imfilter(img_2x2, F02, 'conv');

% Filtro de orden cero con factor T = 4x4
F04 = [1 1 1 1];
F04 = F04.' * F04;
imgF04 = imfilter(img_4x4, F04, 'conv');

% -------------------------------------------------------------------------
% Interpolacion lineal (orden 1)

% Filtro de orden uno con factor T = 2x2
F12 = [0.5 1.0 0.5];
F12 = F12.' * F12;
imgF12 = imfilter(img_2x2, F12, 'conv');

% Filtro de orden uno con factor T = 4x4
F14 = [0.25 0.5 0.75 1.0 0.75 0.5 0.25];
F14 = F14.' * F14;
imgF14 = imfilter(img_4x4, F14, 'conv');

% -------------------------------------------------------------------------
% Interpolacion cubica (orden 3)

% Filtro de orden tres con factor T = 2x2
F32 = [0 -1/16 0 9/16 1  9/16 0 -1/16 0];
F32 = F32.' * F32;
imgF32 = imfilter(img_2x2, F32, 'conv');

% Filtro de orden tres con factor T = 4x4
F34 = [0 -3/128 -1/16 -9/128 0 29/128 9/16 111/128 1 111/128 9/16 29/128 0 -9/128 -1/16 -3/128 0];
F34 = F34.' * F34;
imgF34 = imfilter(img_4x4, F34, 'conv');

% -------------------------------------------------------------------------
% Despliegue de imagenes y espectros
figure('Name','Punto 2. Imagenes con interpolacion espacial');
subplot(2,3,1), imagesc(imgF02), title('Img interpolacion orden cero T = 2x2');
subplot(2,3,4), imagesc(imgF04), title('Img interpolacion orden cero T = 4x4');
subplot(2,3,2), imagesc(imgF12), title('Img interpolacion lineal T = 2x2');
subplot(2,3,5), imagesc(imgF14), title('Img interpolacion lineal T = 4x4');
subplot(2,3,3), imagesc(imgF32), title('Img interpolacion cubica T = 2x2');
subplot(2,3,6), imagesc(imgF34), title('Img interpolacion cubica T = 4x4');
colormap('gray');

figure('Name','Punto 2. Espectros de la interpolacion espacial');
subplot(2,3,1), imagesc( log(1.0 + abs( fftshift(fft2(imgF02)) )) ), title('Espectro interp orden cero T = 2x2');
subplot(2,3,4), imagesc( log(1.0 + abs( fftshift(fft2(imgF04)) )) ), title('Espectro interp orden cero T = 4x4');
subplot(2,3,2), imagesc( log(1.0 + abs( fftshift(fft2(imgF12)) )) ), title('Espectro interp lineal T = 2x2');
subplot(2,3,5), imagesc( log(1.0 + abs( fftshift(fft2(imgF14)) )) ), title('Espectro interp lineal T = 4x4');
subplot(2,3,3), imagesc( log(1.0 + abs( fftshift(fft2(imgF32)) )) ), title('Espectro interp cubica T = 2x2');
subplot(2,3,6), imagesc( log(1.0 + abs( fftshift(fft2(imgF34)) )) ), title('Espectro interp cubica T = 4x4');
colormap('gray');

%%
% 3. Interpolacion en frecuencia.

% Zero padding a la DFT de la imagen original obtenida anteriormente
zp2x2 = padarray(DFT_original, [128,128], 0, 'both');
zp4x4 = padarray(DFT_original, [384,384], 0, 'both');

% Inversa de la DFT de la imagen con zero padding
izp2x2 = ifft2( fftshift(zp2x2) );
izp4x4 = ifft2( fftshift(zp4x4) );


figure('Name','Punto 3. Imagen Original y Magnitud de la DFT');
subplot(1,2,1), imagesc(img), title('Imagen original');
subplot(1,2,2), imagesc(log(abs(DFT_original))), title('DFT imagen original');
colormap('gray');


figure('Name','Punto 3. Magnitudes DFT con ceros');
subplot(3,3,1), imagesc(log(abs(DFT_original))), xlim([0 1024]), ylim([0 1024]), title('DFT imagen original');
subplot(3,3,2), imagesc(log(abs(zp2x2))), xlim([0 1024]), ylim([0 1024]), title('DFT img con padding T = 2');
subplot(3,3,3), imagesc(log(abs(zp4x4))), xlim([0 1024]), ylim([0 1024]), title('DFT img con padding T = 4');
colormap('gray');


figure('Name','Punto 3. Espectros');
subplot(2,4,1), imagesc( log(1.0 + abs( fftshift(fft2(imgF02)) )) ), title('Espectro interp orden cero T = 2x2');
subplot(2,4,5), imagesc( log(1.0 + abs( fftshift(fft2(imgF04)) )) ), title('Espectro interp orden cero T = 4x4');
subplot(2,4,2), imagesc( log(1.0 + abs( fftshift(fft2(imgF12)) )) ), title('Espectro interp lineal T = 2x2');
subplot(2,4,6), imagesc( log(1.0 + abs( fftshift(fft2(imgF14)) )) ), title('Espectro interp lineal T = 4x4');
subplot(2,4,3), imagesc( log(1.0 + abs( fftshift(fft2(imgF32)) )) ), title('Espectro interp cubica T = 2x2');
subplot(2,4,7), imagesc( log(1.0 + abs( fftshift(fft2(imgF34)) )) ), title('Espectro interp cubica T = 4x4');
subplot(2,4,4), imagesc(log(abs(zp2x2))), title('DFT img con padding T = 2');
subplot(2,4,8), imagesc(log(abs(zp4x4))), title('DFT img con padding T = 4');
colormap('gray');

figure('Name','Punto 3. Inversa de la DFT');
subplot(2,3,1), imagesc(img), title('Imagen original');
subplot(2,3,2), imagesc( abs(izp2x2) ), title('Interpolacion T = 2');
subplot(2,3,3), imagesc( abs(izp4x4) ), title('Interpolacion T = 4');
colormap('gray');

%%
% 4. Comparacion de los resultados obtenidos

figure('Name','Punto 4. Comparacion resultados para T = 2x2');
% Punto 1
subplot(3,3,[1 4]), imagesc(img_2x2), title('Sobremuestreo espacial');

% Punto 2
subplot(3,3,2), imagesc(imgF02), title('Interpolacion de orden cero');
subplot(3,3,5), imagesc(imgF12), title('Interpolacion lineal');
subplot(3,3,8), imagesc(imgF32), title('Interpolacion cubica');

% Punto 3
subplot(3,3,6), imagesc(abs(izp2x2)), title('Interpolacion en frecuencia');
colormap('gray');

figure('Name','Punto 4. Comparacion resultados para T = 4x4');
% Punto 1
subplot(3,3,[1 4]), imagesc(img_4x4), title('Sobremuestreo espacial');

% Punto 2
subplot(3,3,2), imagesc(imgF04), title('Interpolacion de orden cero');
subplot(3,3,5), imagesc(imgF14), title('Interpolacion lineal');
subplot(3,3,8), imagesc(imgF34), title('Interpolacion cubica');

% Punto 3
subplot(3,3,6), imagesc(abs(izp4x4)), title('Interpolacion en frecuencia');
colormap('gray');

%% 
% Declaracion de funciones utilizadas

% Funcion de sobremuestreo de la imagen por un factor T
function imaux = sobremuestreo(imagen,T)
    tamy = size(imagen,1);
    tamx = size(imagen,2);
    imauxY = zeros(tamy*T,tamx);
    imaux = zeros(tamy*T,tamx*T);
    for i=0:tamy-1
        imauxY((i*T)+1,:) = imagen(i+1,:);
    end
    for j=0:tamx-1
        imaux(:,(j*T)+1) = imauxY(:,j+1);
    end
    return
end