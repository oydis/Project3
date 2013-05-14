%
% exwalk.m
%
% Example of use of the walk routine
% Generate spanning cluster (l-r spanning)
lx = 64;
ly = 64;
L = 64;
%L = [100 110 120 130 140 150 160 170 180 190 200];
%p = 0.585;
pc = 0.59275;
%p = pc;
p = linspace(pc-0.05,pc+0.1,100);

P = zeros(1,length(p));
M = zeros(1,length(L));
experiments = 50;

for m=1:length(p)
%for m=1:length(L)
    m
    for m2=1:experiments 
        m2
        ncount = 0;
        perc = [];
        while (size(perc ,1)==0)
            ncount = ncount + 1;
            if (ncount >1000)
                'hei'
                return
            end
            %z=rand(L(m),L(m))<p;
            z=rand(L,L)<p(m);
            [lw,num]=bwlabel(z,4);
            %perc_x = intersect(lw(1,:),lw(L(m) ,:)); % label av spanning cluster
            perc_x = intersect(lw(1,:),lw(L ,:));
            perc = find(perc_x >0); % index til spanning cluster
            spanning_index = perc_x(perc);
        end
        s = regionprops(lw,'Area');
        clusterareas = cat(1,s.Area);
%         maxarea = max(clusterareas);
%         i = find(clusterareas==maxarea);
        %spanning_index
        maxarea = clusterareas(spanning_index);
        maxarea_1 = max(maxarea);
        i = find(clusterareas==maxarea_1);
        zz = lw == i;
        % zz now contains the spanning cluster
        %imagesc(zz);
        % Display spanning cluster
        % Run walk on this cluster
        [l,r] = walk(zz);
        zzz = l.*r;
        % Find points where both l and r are non-zero
        zadd = zz + zzz;
    %     figure()
    %     subplot(2,2,1), imagesc(zz);
    %     subplot(2,2,2), imagesc(zadd);
    %     subplot(2,2,3), imagesc(zzz >0);
    %     subplot(2,2,4), imagesc(l+r>0);


        % Mass

        for k1=1:length(zzz)
            for k2=1:length(zzz)
                if zzz(k1,k2) > 0;
                    P(m) = P(m) + 1.0/(L*L)/experiments;
                    %M(m) = M(m) + 1.0/experiments;
                end
            end
        end
    end
end


% plot(L,M)
% xlabel('L')
% ylabel('Mass')
% 
% X = log10(L);
% Y = log10(M);
% 
% ft = fit(X',Y','poly1');
% 
% plot(X,Y,'-o',X,ft(X))
% xlabel('log(L)')
% ylabel('log(M)')
% legend('simulation','linear fit')

plot(p-pc,P,'-o')
xlabel('p-pc')
ylabel('P (singly connected bonds)')



