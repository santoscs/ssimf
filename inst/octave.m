% function omega = ifndq(vimf, dt)
%
%
% INPUT:   
%          vimf:        an IMF;
%          dt:          time interval of the imputted data
% OUTPUT:
%          omega:       instantanesous frequency, which is 2*PI/T, where T
%                       is the period of an ascillation
% NOTE:
%     this is a function to calculate instantaneous based on EMD method--
%     normalize the absolute values and find maximum envelope for 5 times
%     then calculate the Quadrature ,Phase angle,then take difference to them
%     finally,the instantaneous frequency values of an IMF is found. 
%
%Reference:  
%
%
%  code writer:Zhaohua Wu,mailbox:zhwu@cola.iges.org
%  footnote:S.C.Su 2009/05/14
%
% 1.set initial parameters
% 2.find absolute values
% 3.find the spline envelope for AM,take those out-loop start
%   4.Normalize the envelope out (for 5 times)
% 3.find the spline envelope for AM,take those out-loop end 
% 5.flip back those negative values after AM been removed    
% 6.Calculate the quadrature values
% 7.Calculate the differece of the phase angle and give +/- sign
% 8.create a algorithm to remove those outliner
% 9.remove those marked outliner
%10.use cubic spline to smooth the instantaneous frequency values 
%11.return the values back
%
%
% Association:  those procedure of HHT need this code
%1.EMD 
%2.EEMD
%
% Concerned function: no
%
%

function omega = ifndq(vimf, dt)
%
%1.set initial parameters
 Nnormal=5;%number of spline envelope normalization for AM
 rangetop=0.90; %the threshold of outliner remove for instantaneous frequency values
 vlength = max( size(vimf) );
 vlength_1 = vlength -1;

%2.find absolute values
 for i=1:vlength,
     abs_vimf(i)=vimf(i);
     if abs_vimf(i) < 0
         abs_vimf(i)=-vimf(i);
     end
 end

%3.find the spline envelope for AM,take those out-loop start
 for jj=1:Nnormal,
     [spmax, spmin, flag]=extrema(abs_vimf);
     dd=1:1:vlength;
     upper= spline(spmax(:,1),spmax(:,2),dd);
 
%4.Normalize the envelope out 
     for i=1:vlength,
         abs_vimf(i)=abs_vimf(i)/upper(i);
     end
 end
%3.find the spline envelope for AM,take those out-loop end 

%5.flip back those negative values after AM been removed
 for i=1:vlength,
     nvimf(i)=abs_vimf(i);
     if vimf(i) < 0;
         nvimf(i)=-abs_vimf(i);
     end
 end

%6.Calculate the quadrature values
 for i=1:vlength,
     dq(i)=sqrt(1-nvimf(i)*nvimf(i));
 end

%7.Calculate the differece of the phase angle and give +/- sign
 for i=2:vlength_1,
     devi(i)=nvimf(i+1)-nvimf(i-1);
     if devi(i)>0 & nvimf(i)<1
         dq(i)=-dq(i);
     end
 end

%8.create a algorithm to remove those outliner
rangebot=-rangetop;     
 for i=2:(vlength-1),
     if nvimf(i)>rangebot & nvimf(i) < rangetop
        %good original value,direct calculate instantaneous frequency  
         omgcos(i)=abs(nvimf(i+1)-nvimf(i-1))*0.5/sqrt(1-nvimf(i)*nvimf(i));
     else
        %bad original value,direct set -9999,mark them 
         omgcos(i)=-9999;
     end
 end
 omgcos(1)=-9999;
 omgcos(vlength)=-9999;

%9.remove those marked outliner
 jj=1;
 for i=1:vlength,
     if omgcos(i)>-1000
         ddd(jj)=i;
         temp(jj)=omgcos(i);
         jj=jj+1;
     end
 end

%10.use cubic spline to smooth the instantaneous frequency values 
 temp2=spline(ddd,temp,dd); 
 omgcos=temp2;

%11.return the values back
 for i=1:vlength,
     omega(i)=omgcos(i);
 end
 pi2=pi*2;
 omega=omega/dt;
end


%  function [spmax, spmin, flag]= extrema(in_data)
%
% This is a utility program for cubic spline envelope,
%   the code is to  find out max values and max positions
%                            min values and min positions
%    (then use matlab function spline to form the spline)
%
%   function [spmax, spmin, flag]= extrema(in_data)
%
% INPUT:
%       in_data: Inputted data, a time series to be sifted;
% OUTPUT:
%       spmax: The locations (col 1) of the maxima and its corresponding
%              values (col 2)
%       spmin: The locations (col 1) of the minima and its corresponding
%              values (col 2)
%
% NOTE:
%      EMD uses Cubic Spline to be the Maximun and Minimum Envelope for
%        the data.Besides finding spline,end points should be noticed. 
%
%References:  ? which paper?
% 
%
%
% code writer: Zhaohua Wu. 
% footnote:S.C.Su
%
% There are two seperste loops in this code .
% part1.-- find out max values and max positions 
%          process the start point and end point  
% part2.-- find out min values and max positions 
%          process the start point and end point  
% Those parts are similar.
%
% Association:eemd.m
% this function ususally used for finding spline envelope
%
% Concerned function: no
%                     (all matlab internal function)

function [spmax, spmin, flag]= extrema(in_data)

flag=1;
dsize=length(in_data);

%part1.--find local max value and do end process

%start point 
%spmax(1,1)-the first 1 means first point max value,the second 1 means first index
%spmax(1,2)-the first 1 means first point max value,the second 2 means first index
%spmax(1,1)-for position of max 
%spmax(1,2)-for value    of max

spmax(1,1) = 1;
spmax(1,2) = in_data(1);

%Loop --start find max by compare the values 
%when [ (the jj th value > than the jj-1 th value ) AND (the jj th value > than the jj+1 th value )
%the value jj is the position of the max
%the value in_data (jj) is the value of the max
%do the loop by index-jj
%after the max value is found,use index -kk to store in the matrix
%kk=1,the start point
%the last value of kk ,the end point 

jj=2;
kk=2;
while jj<dsize,
    if ( in_data(jj-1)<=in_data(jj) & in_data(jj)>=in_data(jj+1) )
        spmax(kk,1) = jj;
        spmax(kk,2) = in_data (jj);
        kk = kk+1;
    end
    jj=jj+1;
end

%end point
spmax(kk,1)=dsize;
spmax(kk,2)=in_data(dsize);

%End point process-please see reference about spline end effect
%extend the slpoe of neighbor 2 max value ---as extend value
%original value of end point -----as original value
%compare extend and original value 

if kk>=4
    slope1=(spmax(2,2)-spmax(3,2))/(spmax(2,1)-spmax(3,1));
    tmp1=slope1*(spmax(1,1)-spmax(2,1))+spmax(2,2);
    if tmp1>spmax(1,2)
        spmax(1,2)=tmp1;
    end

    slope2=(spmax(kk-1,2)-spmax(kk-2,2))/(spmax(kk-1,1)-spmax(kk-2,1));
    tmp2=slope2*(spmax(kk,1)-spmax(kk-1,1))+spmax(kk-1,2);
    if tmp2>spmax(kk,2)
        spmax(kk,2)=tmp2;
    end
else
    flag=-1;
end

%these 4 sentence seems useless.
msize=size(in_data);
dsize=max(msize);
xsize=dsize/3;
xsize2=2*xsize;


%part2.--find local min value and do end process
%the syntax are all similar with part1.
%here-explan with beginning local max-find upper starting envelope
%the end process procedure-find out the neighbor 2 local extrema value
%connect those 2 local extrema and extend the line to the end
%make judgement with 1).line extend value  2).original data value
%the bigger value is chosen for upper envelope end control point

%local max 
spmin(1,1) = 1;
spmin(1,2) = in_data(1);
jj=2;
kk=2;
while jj<dsize,
    if ( in_data(jj-1)>=in_data(jj) & in_data(jj)<=in_data(jj+1))
        spmin(kk,1) = jj;
        spmin(kk,2) = in_data (jj);
        kk = kk+1;
    end
    jj=jj+1;
end


%local min
spmin(kk,1)=dsize;
spmin(kk,2)=in_data(dsize);

if kk>=4
    slope1=(spmin(2,2)-spmin(3,2))/(spmin(2,1)-spmin(3,1));
    tmp1=slope1*(spmin(1,1)-spmin(2,1))+spmin(2,2);
    if tmp1<spmin(1,2)
        spmin(1,2)=tmp1;
    end

    slope2=(spmin(kk-1,2)-spmin(kk-2,2))/(spmin(kk-1,1)-spmin(kk-2,1));
    tmp2=slope2*(spmin(kk,1)-spmin(kk-1,1))+spmin(kk-1,2);
    if tmp2<spmin(kk,2)
        spmin(kk,2)=tmp2;
    end
else
    flag=-1;
end

flag=1;

end


%   function [sigline, logep] = significance(imfs, percenta)
%
%	that is used to obtain the "percenta" line based on Wu and
%	Huang (2004).
%
%   NOTE:   For this program to work well, the minimum data length is 36
%
%   INPUT:
%	    percenta: a parameter having a value between 0.0 ~ 1.0, e.g., 0.05 
%                 represents 95% confidence level (upper bound); and 0.95 
%                 represents 5% confidence level (lower bound) 
%       imfs:     the true IMFs from running EMD code. The first IMF must
%                 be included for it is used to obtain the relative mean
%                 energy for other IMFs. The trend is not included.
%   OUTPUT:
%       sigline:  a two column matrix, with the first column the natural
%                 logarithm of mean period, and the second column the
%                 natural logarithm of mean energy for significance line
%       logep:    a two colum matrix, with the first column the natural
%                 logarithm of mean period, and the second column the
%                 natural logarithm of mean energy for all IMFs
%
% References can be found in the "Reference" section.
%
% The code is prepared by Zhaohua Wu. For questions, please read the "Q&A" section or
% contact
%   zwu@fsu.edu
%
function [sigline, logep] = significance(imfs, percenta)

nDof = length(imfs(:,1));
pdMax = fix(log(nDof))+1;

pdIntv = linspace(1,pdMax,100);
yBar = -pdIntv;

for i=1:100,
    yUpper(i)=0;
    yLower(i)= -3-pdIntv(i)*pdIntv(i);
end

for i=1:100,
    sigline(i,1)=pdIntv(i);
    
    yPos=linspace(yUpper(i),yLower(i),5000);
    dyPos=yPos(1)-yPos(2);
    yPDF=dist_value(yPos,yBar(i),nDof);
    
    sum = 0.0;
    for jj=1:5000,
        sum = sum + yPDF(jj);
    end
    
    jj1=0;
    jj2=1;
    psum1=0.0;
    psum2=yPDF(1);
    pratio1=psum1/sum;
    pratio2=psum2/sum;
    
    while pratio2 < percenta,
        jj1=jj1+1;
        jj2=jj2+1;
        psum1=psum1+yPDF(jj1);
        psum2=psum2+yPDF(jj2);
        pratio1=psum1/sum;
        pratio2=psum2/sum;
        yref=yPos(jj1);
    end
    sigline(i,2) = yref + dyPos*(pratio2-percenta)/(pratio2-pratio1);
    sigline(i,2) = sigline(i,2) + 0.066*pdIntv(i) + 0.12;
end
sigline=1.4427*sigline;

columns=length(imfs(1,:));
for i=1:columns,
    logep(i,2)=0;
    logep(i,1)=0;
    for j=1:nDof,
        logep(i,2)=logep(i,2)+imfs(j,i)*imfs(j,i);
    end
    logep(i,2)=logep(i,2)/nDof;
end

sfactor=logep(1,2);
for i=1:columns,
    logep(i,2)=0.5636*logep(i,2)/sfactor;  % 0.6441
end

for i=1:3,
    [spmax, spmin, flag]= extrema(imfs(:,i));
    temp=length(spmax(:,1))-1;
    logep(i,1)=nDof/temp;
end
for i=4:columns,
    omega=ifndq(imfs(:,i),1);
    sumomega=0;
    for j=1:nDof,
        sumomega=sumomega+omega(j);
    end
    logep(i,1)=nDof*2*pi/sumomega;
end
logep=1.4427*log(logep);
end

%   function PDF = dist_value(yPos, yBar, nDof)
%
%   INPUT:
%          yPos: An input array at which PDF values are calculated---y value
%                yPos is an 1D 5000 pt matrix for y value in interval-[yUpper,yLower]
%	         yBar: The mean value of y ---------------------------exp(yBar)=E-bar
%	         nDof: The number of degree of freedom----------------Ndof=Npt
%   OUTPUT:
%           PDF: a normalized output array -about the PDF distribution     
%
%   NOTE:   This is a utility program being called by "confidenceLine.m".
%            this code calculate PDF value under yPos range 
%            within the mean value of y(y-bar) of a chi-square distribution
%            here main job is calculating equation(3.4)--PDF formula from the reference paper
%
% References:            
%        'A study of the characteristics of white noise using the empirical mode decomposition method' 
%        Zhaohua,Wu and Norden E. Huang
%        Proc. R. Soc. Lond. A (2004) 460,1597-1611
%
% code writer: zhwu@cola.iges.org
% footnote:S.C.Su 2009/05/31
%
%   1.start to form equation(3.4)--PDF formula
%   2.calculate PDF value for every y value
%   3.to ensure the converge of the calculation,divide by rscale
%
% Association: 
%   this function is called by confidenceLine.m 
%    calculate the PDF value  for 'sigline' 
%
% Concerned function: confidenceLine.m 
%                     the others are matlab functions.  
%

function PDF = dist_value(yPos, yBar, nDof)

ylen = length(yPos);

%1.start to form equation(3.4)--PDF formula
eBar = exp(yBar);%E-bar
evalue=exp(yPos);%E

%2.calculate PDF value for every y value 
for i=1:ylen,
    tmp1 = evalue(i)/eBar-yPos(i);%calculate---tmp1=(E/E-bar)-y
    tmp2 = -tmp1*nDof*eBar/2;%calculate--------tmp2=(-1/2)*E-bar*Ndof*tmp1
    tmp3(i) = 0.5*nDof*eBar*log(nDof) + tmp2;
end

%3.to ensure the converge of the calculation,divide by rscale
%   because the PDF is a relative value,not a absolute value
rscale = max(tmp3);

tmp4 = tmp3 - rscale;%minus means divide,after we take expontial 
PDF= exp(tmp4);
end


%   function [logep] = statistic(imfs)
%
%	that is used to obtain the "percenta" line based on Wu and
%	Huang (2004).
%
%   NOTE:   For this program to work well, the minimum data length is 36
%
%   INPUT:
%	    percenta: a parameter having a value between 0.0 ~ 1.0, e.g., 0.05 
%                 represents 95% confidence level (upper bound); and 0.95 
%                 represents 5% confidence level (lower bound) 
%       imfs:     the true IMFs from running EMD code. The first IMF must
%                 be included for it is used to obtain the relative mean
%                 energy for other IMFs. The trend is not included.
%   OUTPUT:
%       sigline:  a two column matrix, with the first column the natural
%                 logarithm of mean period, and the second column the
%                 natural logarithm of mean energy for significance line
%       logep:    a two colum matrix, with the first column the natural
%                 logarithm of mean period, and the second column the
%                 natural logarithm of mean energy for all IMFs
%
% References can be found in the "Reference" section.
%
% The code is prepared by Zhaohua Wu. For questions, please read the "Q&A" section or
% contact
%   zwu@fsu.edu
%
function [logep] = statistic(imfs)

nDof = length(imfs(:,1));

columns=length(imfs(1,:));
for i=1:columns,
    logep(i,2)=0;
    logep(i,1)=0;
    for j=1:nDof,
        logep(i,2)=logep(i,2)+imfs(j,i)*imfs(j,i);
    end
    logep(i,2)=logep(i,2)/nDof;
end

sfactor=logep(1,2);
for i=1:columns,
    logep(i,2)=0.5636*logep(i,2)/sfactor;  % 0.6441
end

for i=1:3,
    [spmax, spmin, flag]= extrema(imfs(:,i));
    temp=length(spmax(:,1))-1;
    logep(i,1)=nDof/temp;
end
for i=4:columns,
    omega=ifndq(imfs(:,i),1);
    sumomega=0;
    for j=1:nDof,
        sumomega=sumomega+omega(j);
    end
    logep(i,1)=nDof*2*pi/sumomega;
end
logep=1.4427*log(logep);
end


%   function [sigline] = criticalvalue(imfs, percenta)
%
%	that is used to obtain the "percenta" line based on Wu and
%	Huang (2004).
%
%   NOTE:   For this program to work well, the minimum data length is 36
%
%   INPUT:
%	    percenta: a parameter having a value between 0.0 ~ 1.0, e.g., 0.05 
%                 represents 95% confidence level (upper bound); and 0.95 
%                 represents 5% confidence level (lower bound) 
%       imfs:     the true IMFs from running EMD code. The first IMF must
%                 be included for it is used to obtain the relative mean
%                 energy for other IMFs. The trend is not included.
%   OUTPUT:
%       sigline:  a two column matrix, with the first column the natural
%                 logarithm of mean period, and the second column the
%                 natural logarithm of mean energy for significance line
%
% References can be found in the "Reference" section.
%
% The code is prepared by Zhaohua Wu. For questions, please read the "Q&A" section or
% contact
%   zwu@fsu.edu
%
function [sigline] = criticalvalue(n, percenta)

nDof = n;
pdMax = fix(log(nDof))+1;

pdIntv = linspace(1,pdMax,100);
yBar = -pdIntv;

for i=1:100,
    yUpper(i)=0;
    yLower(i)= -3-pdIntv(i)*pdIntv(i);
end

for i=1:100,
    sigline(i,1)=pdIntv(i);
    
    yPos=linspace(yUpper(i),yLower(i),5000);
    dyPos=yPos(1)-yPos(2);
    yPDF=dist_value(yPos,yBar(i),nDof);
    
    sum = 0.0;
    for jj=1:5000,
        sum = sum + yPDF(jj);
    end
    
    jj1=0;
    jj2=1;
    psum1=0.0;
    psum2=yPDF(1);
    pratio1=psum1/sum;
    pratio2=psum2/sum;
    
    while pratio2 < percenta,
        jj1=jj1+1;
        jj2=jj2+1;
        psum1=psum1+yPDF(jj1);
        psum2=psum2+yPDF(jj2);
        pratio1=psum1/sum;
        pratio2=psum2/sum;
        yref=yPos(jj1);
    end
    sigline(i,2) = yref + dyPos*(pratio2-percenta)/(pratio2-pratio1);
    sigline(i,2) = sigline(i,2) + 0.066*pdIntv(i) + 0.12;
end
sigline=1.4427*sigline;
end

