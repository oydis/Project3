clear all; close all;

% Finding data

L = 50;

p = linspace(0.5,0.65,80);
PI = zeros(length(p),1);
experiments = 100;

for i=1:length(p)
    i
    for k=1:experiments
        r = rand(L,L);
        z = r<p(i);
        [lw,num] = bwlabel(z,4);

        img = label2rgb(lw,'jet','k','shuffle');
        s = regionprops(lw,'BoundingBox');
        bbox = cat(1,s.BoundingBox);
    
        test = 0;
    
        for j=1:length(bbox(:,1))
            if (bbox(j,3)==L) || (bbox(j,4)==L)  % Finding index for spanning cluster
                test = 1;
            end
        end
        
        if (test == 1)
            PI(i) = PI(i) + 1.0/experiments;
        end
    end
end

plot(p,PI)
xlabel('p')
ylabel('pi')


% Results

L_array = [25 50 100 200 400 800];
p_03 = [0.55 0.57 0.577 0.585 0.588 0.5905];
p_08 = [0.61 0.6 0.598 0.595 0.5949 0.5943];

figure()
plot(L_array,p_03,'-o',L_array,p_08,'-o')
xlabel('L')
ylabel('p(pi)')
legend('p(pi=0.3)','p(pi=0.8)')

X = log10(L_array)';
Y = log10(p_08 - p_03)';

ft = fit(X,Y,'poly1');

figure()
plot(X,Y,'-o',X,ft(X))
legend('simulation','linear fit')
xlabel('log(L)')
ylabel('log(p(pi=0.8) - p(pi=0.3))')

L_nu = L_array.^(-1/1.33);

ft2 = fit(L_nu',p_03','poly1');
ft3 = fit(L_nu',p_08','poly1');

L_nu2 = linspace(0,0.9,100);

figure()
plot(L_nu,p_03,'-o',L_nu,p_08,'-o',L_nu2,ft2(L_nu2),L_nu2,ft3(L_nu2))
xlabel('L**(-1/nu)')
ylabel('p(pi)')
legend('p(pi=0.3)','p(pi=0.8)','linear fit p(pi=0.3)','linear fit p(pi=0.8)')
axis([0 0.09 0.54 0.63])

% Data collapse to find Phi

L_phi = L_array.^(1/1.33);
pc = 0.59275;





