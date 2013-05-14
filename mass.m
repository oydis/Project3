clear all; close all;

L = 2.^[4 5 6 7 8 9 10 11];
pc = 0.59275;
experiments = 50;

P = zeros(length(L),1);

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
                index_list = [index_list j];
            end
        end
        
        s = regionprops(lw,'Area');
        area = cat(1,s.Area);
        
        for j=1:length(index_list)
            P(i) = P(i) + area(index_list(j))/(L(i)*L(i))/experiments;
        end
    end
end


M = zeros(length(L)-1,1);

for i=1:length(L)-1
    M(i) = P(i)*L(i)*L(i);
end

X = log10(L(1:end-1)');
Y = log10(M);
ft = fit(X,Y,'poly1');

figure()
plot(X,Y,'-o',X,ft(X))
xlabel('log(L)')
ylabel('log(M)')
legend('simulation','fitted line')
axis([1.2 3.2 1.5 6])



