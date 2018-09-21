% Continius Linear Model

a21 = (mb*lb+mw*l)*g/(Ib+mw*l^2);
a22 = -Cb/(Ib+mw*l^2);
a23 = Cw/(Ib+mw*l^2);

a31 = -(mb*lb+mw*l)*g/(Ib+mw*l^2);
a32 = Cb/(Ib+mw*l^2);
a33 = Cw*(Ib+Iw+mw*l^2)/(Iw*(Ib+mw*l^2));

A = [1 0 1; a21 a22 a23; a31 a32 a33]

b21 = -Km/(Ib+mw*l^2);
b31 = Km*(Ib+Iw+mw*l^2)/(Iw*(Ib+mw*l^2));

B = [0; b21; b31];

% Discrete time linear model

% A = sym(A)

[eigenVectorsA, eigenValuesA] = eig(A)

[V,J] = jordan(A)

lambda = sym('lambda');

A-lambda*eye(3)

pol_l = det(A-lambda*eye(3))

solve(lambda^3+27==0,lambda)