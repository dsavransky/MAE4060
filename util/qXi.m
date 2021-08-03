function Xi = qXi(q)

Xi = [q(4)*eye(3) + skew(q(1:3)); -q(1:3).'];
