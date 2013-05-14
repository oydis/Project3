clear all; close all;

Lx = 100;
Ly = 100;
pc = 0.59275;
experiments = 50;

%p = linspace(0.01,1,200);
p = linspace(pc,1,1000);
P = zeros(1,length(p));

for i=1:length(p)
    i
    for k=1:experiments
        r = rand(Lx,Ly);
        z = r<p(i);
        [lw,num] = bwlabel(z,4);

        img = label2rgb(lw,'jet','k','shuffle');
        s = regionprops(lw,'BoundingBox');
        bbox = cat(1,s.BoundingBox);
    
        index_list = [];
    
        for j=1:length(bbox(:,1))
            if (bbox(j,3)==Ly) || (bbox(j,4)==Lx)  % Finding index for spanning cluster
            %if bbox(j,4)==Lx
                index_list = [index_list j];
            end
        end
        
        s = regionprops(lw,'Area');
        area = cat(1,s.Area);
        
        for j=1:length(index_list)
            %if index_list(j)~=0
            P(i) = P(i) + area(index_list(j))/(Lx*Ly)/experiments;
            %end
        end
    end
end

plot(p,P,'-')
xlabel('p')
ylabel('P')
axis([0 1 0 1.15])

X = log10(p-pc)';
Y = log10(P)';

for i=1:length(X)
    test_x = isinf(X(i));
    test_y = isinf(Y(i));
    if (test_x==0) && (test_y==0)
        break
    end
end 

ft = fit(X(i:end), Y(i:end), 'poly1');

figure()
plot(X,Y,'-',X,ft(X))
legend('simulation','fitted linear curve')
xlabel('log10(p-pc)')
ylabel('log10(P)')

