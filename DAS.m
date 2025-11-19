% % clear;
% clc;    %%%%%2 places data read to modify
% tic;
% close all;
%parameters   ---data filtering
tic;
fhigh=1; 
flow=0.01;   %in MHz
filtrue=0;  % if you want to enable the filter, value is 1, otherwise value is 0.

%parameters   --0-data collection
name='Normal_with_reflector(52.9).lvm';            % name of your data
n0=10000;                        % data length for each detector
nd=180;                         % total number of detectors
delay=0;                        % number of points delayed when you record the                              
                                % ultrasound signal
deg=0/2*n0;
%parameters   ---reconstruction
xc=0; yc=0;  
x1=xc-25; x2=xc+25;            % image area in X direction
y1=yc-25; y2=yc+25;            % image area in Y direction
xd=0.05;  yd=0.05;              % image resolution selection. When searching   
                                % reconstruction radius, suggested value is 0.2. After 
                                % getting the reconstruction radius, 0.02  
                                % recommended 
R0=53;  R0d=0.05;  ntry=1;       % R0 is reconstruction radius (the distance 
                                % betw5een image center and ultrasound detectors). R0d is  
                                % the searching step for R0, ntry is the searching times.

v=1.495;                        % ultrasound velocity 
ntime=1/60;                     % ntime=1/f; f is sampling rate (in MHz)

%%%%load data;      
aa=-load(name);
aa=circshift(aa,deg);
result=zeros(n0,nd);
% for i=1:nd
%    result(:,i)=aa(n0*(i-1)+1:n0*i);
% end

for i=1:180
   result(:,i)=aa(n0*(i-1)+1:n0*i);
end


%%%%data filtering
if(filtrue==1)
fhighn=floor(fhigh*n0*ntime);
flown=floor(flow*n0*ntime);
if(fhighn>n0/2)
    fhighn=floor(n0/2)+1; 
end
if(flown<2)
    flown=2;
end
sprintf('high=%d, low=%d',fhighn,flown)
for i=1:nd    
    temp1=zeros(n0,1);
    temp=result(:,i);
    temp=fft(temp);
    for j=flown:fhighn
        temp1(j)=temp(j);
    end
    temp1=ifft(temp1);
    temp1=real(temp1);
    result(:,i)=temp1;
end
end
sprintf('filtering complete')

%%%%image reconstruction and show
nx=floor((x2-x1)/xd)+1;
ny=floor((y2-y1)/yd)+1;
xlabel=x1:xd:x2;
ylabel=y1:yd:y2;
for ntryy=1:ntry
    pics=zeros(ny,nx);
    R=R0+R0d*(ntryy-1);
    for xx=1:nx
        x=x1+(xx-1)*xd;
        for yy=1:ny
            y=y1+(yy-1)*yd;
            for i=1:nd
                if((x-xc)*(x-xc)+(y-yc)*(y-yc)<50*50)
                    detx=R*cos((i-2.8)*2*pi/180);
                    dety=R*sin((i-2.8)*2*pi/180);
                    dis2=(x-detx)^2+(y-dety)^2;
                    dis=floor(sqrt(dis2)/v/ntime)+1;
                    pics(yy,xx)=pics(yy,xx)+result(dis-delay,i);
                end
            end
        end
    end
    figure
    imagesc(xlabel,ylabel,pics);
    colorbar;
    colormap(gray);
    daspect([1 1 1]);
    title(sprintf('%03f mm',R));
    %figure;
end
toc