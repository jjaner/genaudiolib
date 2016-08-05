
function v = sigmf(x, ac)
 a = ac(1);
c = ac(2);
v = 1 ./ (1 + exp(-a * (x - c)));
end