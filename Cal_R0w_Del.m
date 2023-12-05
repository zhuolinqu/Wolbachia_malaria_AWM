function R0w = Cal_R0w_Del(P)
% solve Delta(R0w) = 0
% Delta = b^2 - 4ac with
%   a = vu/vw;
%   b = (1 - ci)/R0w + (vu/vw - 1);
%   c = (1 - R0w)/R0w;

undersqrt = sqrt(P.ci*(-1+P.vw)*P.vw^2*(-1+P.ci*P.vw));

R0wm = P.vw + P.ci*P.vw - 2*P.ci*P.vw^2 - 2 *undersqrt;
R0wp = P.vw + P.ci*P.vw - 2*P.ci*P.vw^2 + 2 *undersqrt;

if P.vw ==1
    R0wm = 1-P.ci;
    R0wp = R0wm;
end

R0w = max(R0wm,R0wp);

b = (1 - P.ci)/R0w + (P.vu/P.vw - 1);

if b > 0 
    R0w = 1;
end

% Calculating a and c as well
%{
if P.vw==1
    a = 0;
else
    a = P.vu/P.vw;
end
c = (1 - R0w)/R0w;
%}