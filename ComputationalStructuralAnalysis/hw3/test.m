% dy/dx + y^2 = 4x

syms y(x);

ode = diff(y,x) + y^2 == 4*x;

ySol1(x) = dsolve(ode)




cond1 = real(y(0)) == 0;
cond2 = real(y(1)) == 0;
conds = [cond1 cond2];

ySol(x) = dsolve(ode,conds)
% ySol = simplify(ySol)

% ySol2(x) = dsolve(ode,conds)
% ySol2 = simplify(ySol2)

figure(1);
plot_y(1)=ezplot(real(ySol(x)));
xlim([0 1]);

saveas(figure(1),'test.png');

quit;
