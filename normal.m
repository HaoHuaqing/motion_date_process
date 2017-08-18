clc
clear

i=1;
num='1.5';
subplot(2,3,6);
str=strcat('disturbance=',num);
hhq_str=strcat('hhq_',num,'.csv');
lxd_str=strcat('lxd_',num,'.csv');
xfl_str=strcat('xfl_',num,'.csv');
hzx_str=strcat('hzx_',num,'.csv');

hhq = csvread(hhq_str, 2, 0);
stickposition=hhq(:,1);
actualposition=hhq(:,2);
roadpositation=hhq(:,3);
a1=actualposition-roadpositation;

lxd = csvread(lxd_str, 2, 0);
stickposition=lxd(:,1);
actualposition=lxd(:,2);
roadpositation=lxd(:,3);
a2=actualposition-roadpositation;

xfl = csvread(xfl_str, 2, 0);
stickposition=xfl(:,1);
actualposition=xfl(:,2);
roadpositation=xfl(:,3);
a3=actualposition-roadpositation;

hzx = csvread(hzx_str, 2, 0);
stickposition=hzx(:,1);
actualposition=hzx(:,2);
roadpositation=hzx(:,3);
a4=actualposition-roadpositation;

data=[a1;a2;a3;a4];
if i>0
    [f,xc]=ecdf(data);
    ecdfhist(f,xc,100);
    h = findobj(gca,'Type','patch');
    h.FaceColor = [.5 .5 .5];
    h.EdgeColor = 'w';
    ylim([0 0.5])
    t = title(str);
    set(t, 'FontSize', 20);
    set(gca,'FontSize',16)
    x=-10:0.05:10;
    y=normpdf(x,mean(data),std(data));
    hold on 
    plot(x,y,'k','LineWidth',2);
    [h,p,stats]=chi2gof(data,'nbins',20)
end

if i==0
    [f,xc]=ecdf(data);
    ecdfhist(f,xc,100);
    h = findobj(gca,'Type','patch');
    h.FaceColor = [.5 .5 .5];
    h.EdgeColor = 'w';
    pdf_normmixture = @(data,p,mu1,mu2,sigma1,sigma2)p*normpdf(data,mu1,sigma1) + (1-p)*normpdf(data,mu2,sigma2);
    pStart = .5;
    muStart = quantile(data,[0.25 0.75])
    sigmaStart = sqrt(var(data) - .25*diff(muStart).^2)
    start = [pStart muStart sigmaStart sigmaStart];
    lb = [0 -Inf -Inf 0 0];
    ub = [1 Inf Inf Inf Inf];
    paramEsts = mle(data, 'pdf',pdf_normmixture, 'start',start, 'lower',lb, 'upper',ub)
    statset('mlecustom')
    options = statset('MaxIter',100000, 'MaxFunEvals',200000);
    paramEsts = mle(data, 'pdf',pdf_normmixture, 'start',start,'lower',lb, 'upper',ub, 'options',options)
    x=-10:0.05:10;
    y = pdf_normmixture(x,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
    z1 = paramEsts(1)*normpdf(x,paramEsts(2),paramEsts(4));
    z2 = (1-paramEsts(1))*normpdf(x,paramEsts(3),paramEsts(5));
    hold on;
    plot(x,z1,'--b','LineWidth',2)
    plot(x,z2,'--b','LineWidth',2)
    plot(x,y,'k','LineWidth',2)
    ylim([0 0.5])
    t = title(str);
    set(t, 'FontSize', 20);
    set(gca,'FontSize',16)
    [h,p,stats]=chi2gof(data,'nbins',20)
end
%disp(paramEsts(1)*(1-paramEsts(1)))