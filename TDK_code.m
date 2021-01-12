%feed (m)
f=0.1e-3;
%N rotation
Nr=100;
%nose radius (m)
r=0.6e-3;
%highest number of circles in a possible intersec.
Nip=ceil(2*r/f);
%amplitude (m)
A=0.130e-4;
% A=0.230e-5;

%wp revolution(rpm)
n=150;
%wp ang. frequency (rad/s)
Om=2*pi*n/60;
%vibration frequency (rad/s)
om=0.484*Om;
% om=0.968*Om;
om=1/3*Om;
% om=0.0692*Om;
%rotation index
is=[0:1:Nr-1];

%phase
Nph=360;
ph=linspace(0,2*pi,Nph);
%arc resolution
Na=500;
    %index start
    ins=10;
    %index end
    ine=10;
%global solutions 
xaG=zeros((Na+1)*( Nr-(ins+ine)),Nph);
zaG=zeros((Na+1)*( Nr-(ins+ine)),Nph);
for k=1:Nph
    %axial movement of the tool
    zi=is*f+ph(k)/(2*pi)*f;
    %radial movement of the tool
    xi=A*sin((2*pi*(om/Om)*(is-1))+ph(k));
    %forward matrix
    xif=hankel(xi(1:Nip+1).',[xi(Nip+1:end) nan*zeros(1,Nip)]);
    zif=hankel(zi(1:Nip+1).',[zi(Nip+1:end) nan*zeros(1,Nip)]);
    %backward matrix
    xip=toeplitz([xi(1) nan*zeros(1,Nip)].',xi);
    zip=toeplitz([zi(1) nan*zeros(1,Nip)].',zi);
    %intersections, size(x1i)==[1 Nr-1]
    
    %forward intersections
    
    d=(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2;
    
    z1if=(repmat(zif(1,1:end-1),[Nip 1])+zif(2:end,1:end-1))/2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).*sqrt(-d.*(-4*r^2+d))./(2*d);
    x1if=((repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*d-sqrt(-d.*(-4*r^2+d)).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)))./(2*d);
    
    z2if=(repmat(zif(1,1:end-1),[Nip 1])+zif(2:end,1:end-1))/2-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).*sqrt(-d.*(-4*r^2+d))./(2*d);
    x2if=((repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*d+sqrt(-d.*(-4*r^2+d)).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)))./(2*d);
    
    %
    % z1if=(sqrt(-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2.*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2)))./(2*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2))+ (repmat(zif(1,1:end-1),[Nip 1])+zif(2:end,1:end-1))/2;
    % x1if=((repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).^2-2.*(repmat(xif(1,1:end-1),[Nip 1]).*(xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1))+(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2)+(sqrt(-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2))).*(-repmat(zif(1,1:end-1),[Nip 1])+zif(2:end,1:end-1))./(2*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2));
    %
    % z2if=-(sqrt(-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2.*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2)))./(2*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2))+ (repmat(zif(1,1:end-1),[Nip 1])+zif(2:end,1:end-1))/2;
    % x2if=((repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).^2-2.*(repmat(xif(1,1:end-1),[Nip 1]).*(xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1))+(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2)+(sqrt(-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2))).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1))./(2*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2));
    
    %hypermatrix definition
    xhif=cat(3,x1if,x2if);
    zhif=cat(3,z1if,z2if);
    xhif(abs(imag(xhif))>0)=nan;                % could not use 'nan' operator for neglecting imaginary part
    zhif(abs(imag(zhif))>0)=nan;
    
    xhif=real(xhif);
    zhif=real(zhif);
    %select minimum of x
    [xIf,indxIf]=min(xhif,[],3);                                %did not work without defining i(error with is(1:end-1))
    zhif=permute(zhif,[3 1 2]);
    zIf=zhif(repmat(is(1:end-1),[Nip 1])*2*Nip+repmat((1:Nip).'-1,[1 Nr-1])*2+indxIf);
    %select real solutions
    
    
    
    %backward intersections
    d=(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2;
    
    z1ip=(repmat(zip(1,1:end-1),[Nip 1])+zip(2:end,1:end-1))/2+(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).*sqrt(-d.*(-4*r^2+d))./(2*d);
    x1ip=((repmat(xip(1,1:end-1),[Nip 1])+xip(2:end,1:end-1)).*d-sqrt(-d.*(-4*r^2+d)).*(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)))./(2*d);
    
    z2ip=(repmat(zip(1,1:end-1),[Nip 1])+zip(2:end,1:end-1))/2-(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).*sqrt(-d.*(-4*r^2+d))./(2*d);
    x2ip=((repmat(xip(1,1:end-1),[Nip 1])+xip(2:end,1:end-1)).*d+sqrt(-d.*(-4*r^2+d)).*(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)))./(2*d);
    
    %
    % z1ip=(sqrt(-(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2.*((repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2)))./(2*((repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2))+ (repmat(zip(1,1:end-1),[Nip 1])+zip(2:end,1:end-1))/2;
    % x1ip=((repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).^2-2.*(repmat(xif(1,1:end-1),[Nip 1]).*(xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1))+(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2)+(sqrt(-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2))).*(-repmat(zif(1,1:end-1),[Nip 1])+zif(2:end,1:end-1))./(2*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2));
    %
    %
    % z2ip=-(sqrt(-(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2.*((repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2)))./(2*((repmat(xip(1,1:end-1),[Nip 1])-xip(2:end,1:end-1)).^2+(repmat(zip(1,1:end-1),[Nip 1])-zip(2:end,1:end-1)).^2))+ (repmat(zip(1,1:end-1),[Nip 1])+zip(2:end,1:end-1))/2;
    % x2ip=((repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).^2-2.*(repmat(xif(1,1:end-1),[Nip 1]).*(xif(2:end,1:end-1)).*(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1))+(repmat(xif(1,1:end-1),[Nip 1])+xif(2:end,1:end-1)).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2)+(sqrt(-(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2).*(-4*r^2+(repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2))).*(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1))./(2*((repmat(xif(1,1:end-1),[Nip 1])-xif(2:end,1:end-1)).^2+(repmat(zif(1,1:end-1),[Nip 1])-zif(2:end,1:end-1)).^2));
    
    
    %hypermatrix definition
    xhip=cat(3,x1ip,x2ip);
    zhip=cat(3,z1ip,z2ip);
    xhip(abs(imag(xhip))>0)=nan;                % could not use 'nan' operator for neglecting imaginary part
    zhip(abs(imag(zhip))>0)=nan;
    
    xhip=real(xhip);
    zhip=real(zhip);
    %select minimum of x
    [xIp,indxIp]=min(xhip,[],3);                                %did not work without defining i(error with is(1:end-1))
    zhip=permute(zhip,[3 1 2]);
    zIp=zhip(repmat(is(1:end-1),[Nip 1])*2*Nip+repmat((1:Nip).'-1,[1 Nr-1])*2+indxIp);
    %select real solutions
    
    
    %calculating all arcs roughly, 1:forward, 3:past
    aifp=permute(repmat(zIf,[1 1 Nip])-permute(repmat(zIp,[1 1 Nip]),[3 2 1]),[1 3 2]);
    %1:forward, 2:past
    %effective arc
    [ae1,temp1]=min(aifp,[],1);
    [ae2,temp2]=min(aifp,[],2);
    
    
    [~,inde2]=min(ae1,[],2);
    [ae,inde1]=min(ae2,[],1);
    %m2n1
    ae=permute(ae,[1 3 2]);
    %effective intersections
    %forward
    xIfe=xIf(is(1:end-1)*Nip+permute(inde1,[1 3 2]));
    zIfe=zIf(is(1:end-1)*Nip+permute(inde1,[1 3 2]));
    xIfe(ae<0)=nan;
    zIfe(ae<0)=nan;
    %past
    xIpe=xIp(is(1:end-1)*Nip+permute(inde2,[1 3 2]));
    zIpe=zIp(is(1:end-1)*Nip+permute(inde2,[1 3 2]));
    xIpe(ae<0)=nan;
    zIpe(ae<0)=nan;
    %centre
    xie=xi(1:end-1);
    zie=zi(1:end-1);
    xie(ae<0)=nan;
    zie(ae<0)=nan;

    %initial and exit angles
    alfa=atan((zie(ins:end-ine)-zIpe(ins:end-ine))./(xie(ins:end-ine)-xIpe(ins:end-ine)));
    betta=atan((zIfe(ins:end-ine)-zie(ins:end-ine))./(xie(ins:end-ine)-xIfe(ins:end-ine)));

    ast=1:Na;
    
    da=(betta+alfa)/Na;
    
    phia=repmat(3*pi/2-alfa,[Na 1])+repmat(da,[Na 1]).*repmat(ast.',[1 Nr-(ins+ine)]);
    
    %
    za=reshape([r*cos(phia)+repmat(zie(ins:end-ine),[Na 1]);nan*ones(1, Nr-(ins+ine));],[(Na+1)*( Nr-(ins+ine)) 1]);
    xa=reshape([r*sin(phia)+repmat(xie(ins:end-ine),[Na 1]);nan*ones(1, Nr-(ins+ine));],[(Na+1)*( Nr-(ins+ine)) 1]);
    xaG(:,k)=xa;
    zaG(:,k)=za;
    k
end

figure;
surf(zaG.',repmat(ph,[size(zaG,1),1]).',xaG.'+r,'edgecolor','none');
view(0,90);
xlabel('axial direction, z(m)');
ylabel('phase, \phi (rad)');

%
%
%
%
% %check flyover
%
%
% figure;
% plot(zI,xI+r,'.')
%
% %
% % alfa=0:pi/20:2*pi
% % for j=1:1:20
% % a=r*cos(alfa)+zi(j)
% % b=r*sin(alfa)+xi(j)
% % plot(a,b)
% % hold on
% % end
%
% hold on;
% alfa=atan((zi(2:end-1)-zI(1:end-1))./(xi(2:end-1)-xI(1:end-1)));
% betta=atan((zI(2:end)-zi(2:end-1))./(xi(2:end-1)-xI(1:end-1)));
%
% Na=500;
% ast=1:Na;
% da=(betta+alfa)/Na;
%
% phia=repmat(3*pi/2-alfa,[Na 1])+repmat(da,[Na 1]).*repmat(ast.',[1 Nr-2]);
%
% za=reshape([r*cos(phia)+repmat(zi(2:end-1),[Na 1]);nan*ones(1,Nr-2);],[(Na+1)*(Nr-2) 1]);
% xa=reshape([r*sin(phia)+repmat(xi(2:end-1),[Na 1]);nan*ones(1,Nr-2);],[(Na+1)*(Nr-2) 1]);
%
% hold on;plot(za,xa+r,'k-');
%
%
%
% xa=alfa, betta
% % fi=0:pi/500:2*pi;
%
% for u=2:1:Nr-1
%
% fi=(3*pi/2-alfa(u-1)):(betta(u-1)+alfa(u-1))/Na:(3*pi/2+betta(u-1))
% a=r*cos(fi)+zi(u)
% b=r*sin(fi)+xi(u)
%     plot(a,b+r)
%
%     hold on
% end
%
%
