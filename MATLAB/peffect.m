N = 1024;
n = log2(N);

p = 0.5;
pval = zeros(1,N);
cval = zeros(1,N);
pval(1) = p;
for i = 2:n+1
    for j = 1:2^(i-2)
        cval(2*j - 1) = pval(j)^2;
        cval(2*j) = 2*pval(j) - pval(j)^2;
    end
    pval = cval;
end
% disp(y);
y = fliplr(1 - cval);
index = find(y>0.6);
disp(index);
figure(1)
plot(y,'.');
xlabel('Channel Index');
ylabel('Symmetric Chanel Capacity');
xlim([1 N]);

figure(2)
plot(sort(y),'.');
xlabel('Channel Index in sorted order');
ylabel('Symmetric Chanel Capacity');
xlim([1 N]);