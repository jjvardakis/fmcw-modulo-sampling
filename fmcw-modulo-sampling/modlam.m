function z = modlam(y, lam)
%MODLAM  Centred modulo nonlinearity of the paper: M_lam(y) = 2*lam*( y/(2*lam) - floor((y+lam)/(2*lam)) )
%   Output lies in [-lam, lam).  Same size as y (here N x 1).
%   Equivalently z = y - eps with eps in 2*lam*Z  (Eq. 2).
z = 2*lam*( y/(2*lam) - floor((y + lam)/(2*lam)) );
end
