%
% exflow.m
%
clear all; clf;
% First , find the backbone
% Generate spanning cluster (l-r spanning)
%lx = 100;
%ly = 100;
lx = [60 80 100 120 140 160 180 200 220 240 260 280 300];
ly = [60 80 100 120 140 160 180 200 220 240 260 280 300];
pc = 0.59275;
p1 = pc;
%p1 = linspace(pc-0.05,pc+0.1,100);
experiments = 100;

sc = zeros(length(lx),1);
bb = zeros(length(lx),1);
de = zeros(length(lx),1);

conductivity = zeros(length(lx),1);
%conductivity = zeros(length(p1),1);

for m1=1:length(lx)
%for m1=1:length(p)
    m1
    for m2=1:experiments
        ncount = 0;
        perc = [];
        while (size(perc ,1)==0)
            ncount = ncount + 1;
            if (ncount >1000)
                'hei'
                return
            end
            %z=rand(lx,ly)<p1;
            z=rand(lx(m1),ly(m1))<p1;
            %z=rand(lx,ly)<p1(m1);
            [lw,num]=bwlabel(z,4);
            %perc_x = intersect(lw(1,:),lw(lx ,:));
            perc_x = intersect(lw(1,:),lw(lx(m1) ,:));
            perc = find(perc_x >0);
            spanning_index = perc_x(perc);
        end
        s = regionprops(lw,'Area');
        clusterareas = cat(1,s.Area);
        maxarea = clusterareas(spanning_index);
        maxarea_1 = max(maxarea);
        i = find(clusterareas==maxarea_1);
        zz = lw == i;
        % zz now contains the spanning cluster
        % Transpose
        zzz = zz';
        % Generate bond lattice from this
        g = sitetobond(zzz);
        % Generate conductivity matrix
        [p c_eff] = FIND_COND(g,lx(m1),ly(m1));
        % Transform this onto a nx x ny lattice
        x = coltomat(full(p),lx(m1),ly(m1));
        P = x.*zzz;
        g1 = g(:,1);
        g2 = g(:,2);
        z1 = coltomat(g1,lx(m1),ly(m1));
        z2 = coltomat(g2,lx(m1),ly(m1));
        % Plotting
    %     subplot(2,2,1), imagesc(zzz);
    %     title('Spanning cluster')
    %     axis equal
    %     subplot(2,2,2), imagesc(P);
    %     title('Pressure');
    %     axis equal
        f2 = zeros(lx(m1),ly(m1));
        for iy = 1:ly(m1)-1
            f2(:,iy) = (P(:,iy) - P(:,iy+1)).*z2(:,iy);
        end
        f1 = zeros(lx(m1),ly(m1));
        for ix = 1:lx(m1)-1
            f1(ix ,:) = (P(ix ,:) - P(ix+1,:)).*z1(ix ,:);
        end
        % Find the sum of absolute fluxes into each site
        fn = zeros(lx(m1),ly(m1));
        fn = fn + abs(f1);
        fn = fn + abs(f2);
        fn(:,2:ly(m1)) = fn(:,2:ly(m1)) + abs(f2(:,1:ly(m1) -1));
        fn(:,1) = fn(:,1) + abs((P(:,1) - 1.0).*(zzz(:,1)));
        fn(2:lx(m1) ,:) = fn(2:lx(m1) ,:) + abs(f1(1:lx(m1) -1,:));
    %     subplot(2,2,3), imagesc(fn);
    %     title('Flux');
    %     axis equal
        %limit = max(max(fn))-0.0001; % ?
        limit = 0.00001;
        zfn = fn>limit;
        zbb = (zzz + 2*zfn);
        zbb = zbb/max(max(zbb));
    %     subplot(2,2,4), imagesc(zbb);
    %     title('BB and DE');
    %     axis equal

        sc(m1) = sc(m1) + ...
            length(find(fn > (max(max(fn))-0.00001)))/experiments;
        bb_temp = length(find(zfn==1));
        bb(m1) = bb(m1) + bb_temp/experiments;
        de(m1) = de(m1) + (length(find(zzz==1)) - bb_temp)/experiments;

        conductivity(m1) = conductivity(m1) + sum(fn(:,1))/experiments;
    end
end

subplot(2,1,1)
plot(lx,sc,'-x')
xlabel('L')
ylabel('M')
legend('Singly connected bonds')
subplot(2,1,2)
plot(lx,bb,'-x',lx,de,'-x')
xlabel('L')
ylabel('M')
legend('Backbone','Dangling ends')

Y1 = log10(sc);
Y2 = log10(bb);
Y3 = log10(de);
X = log10(lx');

ft1 = fit(X,Y1,'poly1');
ft2 = fit(X,Y2,'poly1');
ft3 = fit(X,Y3,'poly1');

figure()
%subplot(2,1,1)
plot(X,Y1,'-x',X,ft1(X))
xlabel('log(L)')
ylabel('log(M)')
legend('Singly connected bonds','Linear fit')
%axis([1.75 2.5 1.2 2.3])
%subplot(2,1,2)
figure()
plot(X,Y2,'-x',X,ft2(X),X,Y3,'-x',X,ft3(X))
xlabel('log(L)')
ylabel('log(M)')
legend('Backbone','bb linear fit','Dangling ends','de linear fit')
%axis([1.75 2.8 2.5 5])

figure()
subplot(2,1,1)
plot(lx,conductivity,'-x')
xlabel('L')
ylabel('Conductivity')

YC = log10(conductivity);
ftC = fit(X,YC,'poly1');

subplot(2,1,2)
plot(X,YC,'-x',X,ftC(X))
xlabel('log(L)')
ylabel('log(conductivity)')
legend('simulation','fitted line')







