clear all; close all;

%L = 100;
L = 2.^[4 5 6 7 8 9];
pc = 0.59275;
experiments = 50;

%p = linspace(pc-0.05,pc+0.05,10);
%P = zeros(1,length(p));

s_max = L(6)*L(6);
s_array = logspace(log10(1),log10(s_max),30);
%N_s = zeros(length(p),length(s_array)-1);

N_s = zeros(length(L),length(s_array)-1);

for i=1:length(L)
    i
    for k=1:experiments
        r = rand(L(i),L(i));
        z = r<pc;
        [lw,num] = bwlabel(z,4);

        img = label2rgb(lw,'jet','k','shuffle');
        s = regionprops(lw,'BoundingBox');
        bbox = cat(1,s.BoundingBox);
    
        index_list = [];
        index2 = [];
        
        for j=1:length(bbox(:,1))
            if (bbox(j,3)==L(i)) || (bbox(j,4)==L(i))  % Finding index for spanning cluster
            %if bbox(j,4)==Lx
                index_list = [index_list j];
            else
                index2 = [index2 j];
            end
        end
        
        s = regionprops(lw,'Area');
        area = cat(1,s.Area);
        
        
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

n = zeros(length(L),length(s_array)-1);

for i=1:length(L)
    n(i,:) = N_s(i,:)./(L(i)*L(i));
end

X = log10(s_array(2:end-6)');
Y = log10(n(6,2:end-5)');

ft = fit(X, Y, 'poly1');

loglog(s_array(1:end-1),n(1,:),'r', ... 
    s_array(1:end-1),n(2,:),'g', ...
    s_array(1:end-1),n(3,:),'b', ...
    s_array(1:end-1),n(4,:),'m', ...
    s_array(1:end-1),n(5,:),'c', ...
    s_array(1:end-1),n(6,:),'k' ...
    )
legend('L=16','L=32','L=64','L=128','L=256','L=512')
xlabel('cluster size s')
ylabel('cluster number density n(s,p)')

figure()
plot(log10(s_array(1:end-1)),log10(n(6,:)),X,ft(X))

figure()
semilogx(s_array(1:end-1),n(1,:),'r', ... 
    s_array(1:end-1),n(2,:),'g', ...
    s_array(1:end-1),n(3,:),'b', ...
    s_array(1:end-1),n(4,:),'m', ...
    s_array(1:end-1),n(5,:),'c', ...
    s_array(1:end-1),n(6,:),'k' ...
    )
legend('L=16','L=32','L=64','L=128','L=256','L=512')
xlabel('log(cluster size s)')
ylabel('cluster number density n(s,p)')


