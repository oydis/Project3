clear all; close all;

z = rand(1e6,1).^(-3+1);

Z = linspace(min(z),max(z)+1,200);
P = zeros(length(Z),1);
f = zeros(length(Z),1);

hist(z,50)

for i=1:length(Z)
    for j=1:length(z)
        if i==1
            if z(j) < Z(i)
                f(i) = f(i) + 1.0/length(z);
            end
        else
            if (z(j) > Z(i-1)) && (z(j) < Z(i))
                f(i) = f(i) + 1.0/length(z);
            end
        end
    end
end
figure()
plot(Z,f)
xlabel('Z')
ylabel('Distribution f')

% for i=1:length(Z)
%     for j=1:length(z)
%         if z(j)<Z(i)
%             P(i) = P(i) + 1.0;
%         end
%     end
% end

for i=1:length(Z)
    if i==1
        P(i) = f(i);
    else        
        P(i) = P(i-1) + f(i); 
    end
end

figure()
plot(Z,P)
xlabel('Z')
ylabel('Cumulative distribution P')
     


% Logarithmic binning

Z_log = logspace(log10(min(z)),log10(max(z)),200);
P_log = zeros(length(Z_log)-1,1);

% for i=1:length(Z_log)
%     if i==1
%         dz = Z(i);
%     else
%         dz = Z_log(i) - Z_log(i-1);
%     end
%     for j=1:length(z)
%         if z(j)<Z_log(i)
%             P_log(i) = P_log(i) + 1.0/dz; 
%         end
%     end
% end

f_log = zeros(length(Z_log)-1,1);

for i=1:length(Z_log)-1
%     if i==1
%         dz = Z_log(i);
%     else
    dz = Z_log(i+1) - Z_log(i);
%     end
    for j=1:length(z)
%         if i==1
%             if z(j) < Z_log(i)
%                 f_log(i) = f_log(i) + 1.0/dz/length(z);
%             end
%         else
        if (Z_log(i) < z(j)) && (z(j) < Z_log(i+1))
            f_log(i) = f_log(i) + 1.0/length(z)/dz;
%             end
        end
    end
end

%Y = zeros(length(Z_log),1);

for i=1:length(Z_log)-1
    
    dz = Z_log(i+1) - Z_log(i);
    
    if i==1
        P_log(i) = f_log(i)/sum(f_log);
    else        
        P_log(i) = P_log(i-1) + f_log(i)/sum(f_log); 
    end
end

%P_log = P_log./sum(f_log);

figure()
plot(log10(Z_log(1:end-1)),P_log)
xlabel('log Z')
ylabel('Cumulative distribution P')

figure()
plot(log10(Z_log(1:end-1)),f_log)
xlabel('log Z')
ylabel('Actual distribution f')

Y = 1 - P_log; % ?

ft = fit(log(Z_log(1:end-1)'), log(Y), 'poly1');

figure()
plot(log(Z_log(1:end-1)),log(Y),log(Z_log(1:end-1)),ft(log(Z_log(1:end-1))))
legend('simulation','fitted line')
xlabel('log Z')
ylabel('log(1-P)')

