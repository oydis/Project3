clear all; close all;

L = 100;
pc = 0.59275;
experiments = 50;

p = linspace(pc-0.2,pc,10);
P = zeros(1,length(p));

s_max = L*L;
s_array = logspace(log10(1),log10(s_max),50);
N_s = zeros(length(p),length(s_array)-1);

for i=1:length(p)
    i
    for k=1:experiments
        r = rand(L,L);
        z = r<p(i);
        [lw,num] = bwlabel(z,4);

        img = label2rgb(lw,'jet','k','shuffle');
        s = regionprops(lw,'BoundingBox');
        bbox = cat(1,s.BoundingBox);
    
        index_list = [];
        index2 = [];
        
        for j=1:length(bbox(:,1))
            if (bbox(j,3)==L) || (bbox(j,4)==L)  % Finding index for spanning cluster
            %if bbox(j,4)==Lx
                index_list = [index_list j];
            else
                index2 = [index2 j];
            end
        end
        
        s = regionprops(lw,'Area');
        area = cat(1,s.Area);
        
        for j=1:length(index_list)
            %if index_list(j)~=0
            P(i) = P(i) + area(index_list(j))/(L*L)/experiments;
            %end
        end
        
        for j=1:length(index2)
            for l=1:length(s_array)-1
                ds = s_array(l+1) - s_array(l);
                if (s_array(l) < area(index2(j))) && (area(index2(j)) < s_array(l+1))
                    N_s(i,l) = N_s(i,l) + 1.0/experiments/ds;
                end
            end
        end
    end
end

% plot(p,P,'-')
% xlabel('p')
% ylabel('P')
% axis([0 1 0 1])

n = N_s./(L*L);

figure()
loglog(s_array(1:end-1),n(1,:),'r', ... 
    s_array(1:end-1),n(2,:),'g', ...
    s_array(1:end-1),n(3,:),'b', ...
    s_array(1:end-1),n(4,:),'m', ...
    s_array(1:end-1),n(5,:),'c', ...
    s_array(1:end-1),n(6,:),'k', ...
    s_array(1:end-1),n(7,:),'r--', ...
    s_array(1:end-1),n(8,:),'g--', ...
    s_array(1:end-1),n(9,:),'b--', ...
    s_array(1:end-1),n(10,:),'m--' ...
    )
legend('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
xlabel('cluster size s')
ylabel('cluster number density n(s,p)')
% 
% figure()
% semilogx(s_array(1:end-1),n(1,:),'r', ... 
%     s_array(1:end-1),n(2,:),'g', ...
%     s_array(1:end-1),n(3,:),'b', ...
%     s_array(1:end-1),n(4,:),'m', ...
%     s_array(1:end-1),n(5,:),'c', ...
%     s_array(1:end-1),n(6,:),'k', ...
%     s_array(1:end-1),n(7,:),'r--', ...
%     s_array(1:end-1),n(8,:),'g--', ...
%     s_array(1:end-1),n(9,:),'b--', ...
%     s_array(1:end-1),n(10,:),'m--' ...
%     )
% legend('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
% xlabel('log(cluster size s)')
% ylabel('cluster number density n(s,p)')

% Finding s_xi

n_c = zeros(length(p),length(s_array)-1);
F_half = zeros(length(p),1);

for i=1:length(p)
    n_c(i,:) = n(i,:).*(s_array(1:end-1).^(1.87));
    F_half(i) = 0.5*max(log10(n_c(i,:)));
end

n_c2 = zeros(length(p)-1,length(s_array)-1);
F_half2 = zeros(length(p)-1,1);

for i=1:length(p)-1
    n_c2(i,:) = n(i,:)./n(10,:);
    F_half2(i) = 0.5*max(log10(n_c2(i,:)));
end

figure()
loglog(s_array(1:end-1),n_c(1,:),'r', ... 
    s_array(1:end-1),n_c(2,:),'g', ...
    s_array(1:end-1),n_c(3,:),'b', ...
    s_array(1:end-1),n_c(4,:),'m', ...
    s_array(1:end-1),n_c(5,:),'c', ...
    s_array(1:end-1),n_c(6,:),'k', ...
    s_array(1:end-1),n_c(7,:),'r--', ...
    s_array(1:end-1),n_c(8,:),'g--', ...
    s_array(1:end-1),n_c(9,:),'b--', ...
    s_array(1:end-1),n_c(10,:),'m--' ...
    )
legend('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
xlabel('cluster size s')
ylabel('n(s,p)*s**tau')
    
% figure()
% semilogx(s_array(1:end-1),n_c(1,:),'r', ... 
%     s_array(1:end-1),n_c(2,:),'g', ...
%     s_array(1:end-1),n_c(3,:),'b', ...
%     s_array(1:end-1),n_c(4,:),'m', ...
%     s_array(1:end-1),n_c(5,:),'c', ...
%     s_array(1:end-1),n_c(6,:),'k', ...
%     s_array(1:end-1),n_c(7,:),'r--', ...
%     s_array(1:end-1),n_c(8,:),'g--', ...
%     s_array(1:end-1),n_c(9,:),'b--', ...
%     s_array(1:end-1),n_c(10,:),'m--' ...
%     )
% legend('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
% xlabel('log(cluster size s)')
% ylabel('cluster number density n(s,p)*s^tau')

figure()
loglog(s_array(1:end-1),n_c2(1,:),'r', ... 
    s_array(1:end-1),n_c2(2,:),'g', ...
    s_array(1:end-1),n_c2(3,:),'b', ...
    s_array(1:end-1),n_c2(4,:),'m', ...
    s_array(1:end-1),n_c2(5,:),'c', ...
    s_array(1:end-1),n_c2(6,:),'k', ...
    s_array(1:end-1),n_c2(7,:),'r--', ...
    s_array(1:end-1),n_c2(8,:),'g--', ...
    s_array(1:end-1),n_c2(9,:),'b--' ...
    )
legend('p1','p2','p3','p4','p5','p6','p7','p8','p9')
xlabel('cluster size s')
ylabel('n(s,p)/n(s,pc)')
  
% length av s = 30

% s_xi = [s_array(13) s_array(14) s_array(15) s_array(16) s_array(17) ...
%     s_array(18) s_array(20) s_array(22) s_array(26)];

% length av s = 50

s_xi = [s_array(22) s_array(23) s_array(24) s_array(25) s_array(28) ...
     s_array(30) s_array(33) s_array(36) s_array(41)];

figure()
plot(p(1:end-1),s_xi,'-o')
xlabel('p')
ylabel('characteristic size s')

X = log10(s_xi)';
Y = log10(abs(p(1:end-1)-pc))';

ft = fit(X,Y,'poly1');

figure()
plot(X,Y,'-o',X,ft(X))
xlabel('log(p-pc)')
ylabel('log(characteristic size s)')
legend('simulation','linear fit')


% Better way

x_array = zeros(length(p),length(s_array)-1);

for i=1:length(p)
    x_array(i,:) = s_array(1:end-1).*abs(p(i)-pc).^(1/0.4);
end
    
figure()
loglog(x_array(1,:),n_c(1,:),'r', ... 
    x_array(2,:),n_c(2,:),'g', ...
    x_array(3,:),n_c(3,:),'b', ...
    x_array(4,:),n_c(4,:),'m', ...
    x_array(5,:),n_c(5,:),'c', ...
    x_array(6,:),n_c(6,:),'k', ...
    x_array(7,:),n_c(7,:),'r--', ...
    x_array(8,:),n_c(8,:),'g--', ...
    x_array(9,:),n_c(9,:),'b--', ...
    x_array(10,:),n_c(10,:),'m--' ...
    )
legend('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
xlabel('|p-pc|**(-1/sigma)')
ylabel('s**tau*n(s,p)')

