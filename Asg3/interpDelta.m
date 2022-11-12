% Q3b
function delta = interpDelta(delta_n,S_n, S)
% input: column vectors
%       delta_n, S_n are layers of delta and S at a time
% output: column vectors
%       delta: interpolated delta 

delta = zeros(size(S));
ub = max(S_n);
lb = min(S_n);

idx = lb <= S & S <= ub;
% normal linear interpolation for S that is in range
delta(idx) = interp1(S_n, delta_n, S(idx));

I = find(~idx);
% disp(I)

% note that the index array better to be a row vector, otherwise we have
% different behavior 

% use S is out of range, use nearest S_n's delta
for i = I'
    [~, idx2] = min(abs(S_n - S(i)));
    delta(i) = delta_n(idx2);
end

end

