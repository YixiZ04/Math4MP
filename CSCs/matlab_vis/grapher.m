nsim = 6;

[poptot,pops,nec,act] = dcode_simfiles(nsim);


%%
p1 = [];
p2 = [];
p3 = [];
p4 = [];
p5 = [];

for i=1:10
    a = pops{i};
    p1 = [p1 sum(a(:,:,:,1), 'all')];
    p2 = [p2 sum(a(:,:,:,2), 'all')];
    p3 = [p3 sum(a(:,:,:,3), 'all')];
    p4 = [p4 sum(a(:,:,:,4), 'all')];
    p5 = [p5 sum(a(:,:,:,5), 'all')];
end

n = length(p1);
ts = 1:1:n;
length(ts)

plot(ts, p1)
hold on
plot(ts, p2)
plot(ts, p3)
plot(ts, p4)
plot(ts, p5)
legend('CSC', 'prog-1', 'prog-2', 'hija prog-1', 'hija prog-2')