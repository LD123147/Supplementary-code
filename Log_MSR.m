clc
clear
close all

%% Open an image selection interface
uiopen
obj = get(gca,'Children');      
pics = get(obj,'CData'); 
x_limit = xlim ;
y_limit = ylim ;
[xx,yy] = size(pics);
xd = 0.05 ; 
yd = 0.05 ; %%step
xlabel=x_limit(1):xd:x_limit(2);
ylabel=y_limit(1):yd:y_limit(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%-------fitting-------%%%%%%%%%%%%%%%%%%%%%%%%%%
 [nn,mm]=size(pics);
    for i = 1 : nn
       for j = 1 : mm
          if isnan(pics(i,j)) 
             pics(i ,j) = 0 ; 
          end
       end
    end
    pics = ( pics - min(min(pics)) ) / ( max(max(pics)) - min(min(pics)) ) ;

    x0 = zeros(nn,1) ;
    y0 = zeros(nn,1) ;

    y0_1 = zeros(nn,floor(mm/10)+1) ;

      for i = 1 : nn
            for j = 1 : 10 : mm
                y0_1(i,floor(j/10)+1) = pics(i,j) ;
            end
      end  

      for i = 1 : nn
          for j = 1 : floor(mm/10)+1
              x0(i) = i ;
              y0(i) = y0(i) + y0_1(i,j);
          end
      end

    y0 = y0/floor(mm/10)+1 ;

    f = fittype('a1.*x.^(-a2)'); % Form of fitting function
    fit1 = fit(x0,y0,f,'StartPoint',[x0(1) y0(1)]);
    a1 = fit1.a1; % The value of a
    a2 = abs(fit1.a2); % The value of b
    fdata = feval(fit1,x0); % Calculate y using a fitting function

    aa = min(fdata);
    if fdata(1) > fdata(end)
        fdata_1 =  fdata + 2 * abs(aa) + 1 ;
    else
        fdata_1 =  -fdata + 2 * abs(aa) + 1 ;
    end

    fdata_2 =  zeros(nn,1) ;
    for i = 1 : nn
       fdata_2(i) = fdata_1(1) - fdata_1(i)  ; 
    end

    fdata_3 = zeros(nn,1);
    for i = 1 : nn
       fdata_3(i) = fdata_2(i) / a2 * 10 ; 
    end

    fdata_4 =  zeros(nn,1) ;
    for i = 1 : nn
        fdata_4(i) =  fdata_3(i) + 1 ;
    end

    fdata_5 =  zeros(nn,1) ;
    for i = 1 : nn
        fdata_5(i) = fdata_4(i) - ( i / floor(nn/10) ) ;
    end
    
    for i = 1 : nn
        if fdata_5(i) <= fdata_5(50)
            fdata_5(i) = fdata_5(50) ;
        end
    end

    figure
    plot(x0,fdata_5,'r');

            pics_l = zeros(nn,mm);

            for i = 1 : nn
                pics_l(i,:) = pics(i,:) * fdata_5(i) ;
            end
            pics_l = ( pics_l - min(min(pics_l)) ) / ( max(max(pics_l)) - min(min(pics_l)) ) ;
            figure;
            imagesc(xlabel,ylabel,pics_l);  
            colorbar;
            colormap(parula);
            daspect([1 1 1]);
%             axis equal;
            title(sprintf('Log compensation'));
    %         if saves_if == 1
    %             saveas(gcf,dir0_3,'fig');
    %         end
    %%
    [rows,cols] = size(pics_l) ;

    for i = 1 : rows
       for j = 1 : cols
          if isnan(pics_l(i,j))
             pics_l(i,j) = 0 ; 
          end
       end
    end

    r = double(pics_l) ;
    g = double(pics_l) ;
    b = double(pics_l) ;
    I = cat(3,r,g,b) ; 

    I_r = double(I(:,:,1));
    I_g = double(I(:,:,2));
    I_b = double(I(:,:,3));

    I_r_log = log(I_r+1);
    I_g_log = log(I_g+1);
    I_b_log = log(I_b+1);

    Rfft1 = fft2(I_r);
    Gfft1 = fft2(I_g);
    Bfft1 = fft2(I_b);[m,n] = size(I_r);
    sigma1 = 15;
    sigma2 = 80;
    sigma3 = 200;
    f1 = fspecial('gaussian', [m, n], sigma1);
    f2 = fspecial('gaussian', [m, n], sigma2);
    f3 = fspecial('gaussian', [m, n], sigma3);
    efft1 = fft2(double(f1));
    efft2 = fft2(double(f2));
    efft3 = fft2(double(f3));

    D_r1 = ifft2(Rfft1.*efft1);
    D_g1 = ifft2(Gfft1.*efft1);
    D_b1 = ifft2(Bfft1.*efft1);
    D_r_log1 = log(D_r1 + 1);
    D_g_log1 = log(D_g1 + 1);
    D_b_log1 = log(D_b1 + 1);
    R1 = I_r_log - D_r_log1;
    G1 = I_g_log - D_g_log1;
    B1 = I_b_log - D_b_log1;

    D_r2 = ifft2(Rfft1.*efft2);
    D_g2 = ifft2(Gfft1.*efft2);
    D_b2 = ifft2(Bfft1.*efft2);
    D_r_log2 = log(D_r2 + 1);
    D_g_log2 = log(D_g2 + 1);
    D_b_log2 = log(D_b2 + 1);
    R2 = I_r_log - D_r_log2;
    G2 = I_g_log - D_g_log2;
    B2 = I_b_log - D_b_log2;

    D_r3 = ifft2(Rfft1.*efft3);
    D_g3 = ifft2(Gfft1.*efft3);
    D_b3 = ifft2(Bfft1.*efft3);
    D_r_log3 = log(D_r3 + 1);
    D_g_log3 = log(D_g3 + 1);
    D_b_log3 = log(D_b3 + 1);
    R3 = I_r_log - D_r_log3;
    G3 = I_g_log - D_g_log3;
    B3 = I_b_log - D_b_log3;

    R = 0.1*R1 + 0.4*R2 + 0.5*R3;
    G = 0.1*G1 + 0.4*G2 + 0.5*G3;
    B = 0.1*B1 + 0.4*B2 + 0.5*B3;

    R = exp(R);
    R = double(R);
    MIN = min(min(R)); 
    MAX = max(max(R));
    R = (R - MIN)/(MAX - MIN);
    R = adapthisteq(R);

    G = exp(G);
    G = double(G);
    MIN = min(min(G)); 
    MAX = max(max(G));
    G = (G - MIN)/(MAX - MIN);
    G = adapthisteq(G);

    B = exp(B);
    B = double(B);
    MIN = min(min(B)); 
    MAX = max(max(B));
    B = (B - MIN)/(MAX - MIN);
    B = adapthisteq(B);

    J_1 = cat(3, R, G, B);

    J = rgb2gray(J_1) ;

    figure ;
    imagesc(xlabel,ylabel,J) ;
%     axis xy;
    colorbar ;
    colormap(gray);
    daspect([1 1 1]);
%     axis equal
    title(sprintf('Log+MSR compensation'));
