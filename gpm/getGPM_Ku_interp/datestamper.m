function time = datestamper(t)

t1 = num2str(t(1));
t2 = num2str(t(2));
t3 = num2str(t(3));
t4 = num2str(t(4));

if t(5)<10
    t5 = ['0', num2str(t(5))];
else
    t5 = num2str(t(5));
end

if t(6)<10
    t6 = ['0', num2str(floor(t(6)))];
else
    t6 = num2str(floor(t(6)));
end

time = [t3, '/', t2, '/', t1, '   ', t4, ':', t5, ':', t6];

fprintf('\t')
disp(time)