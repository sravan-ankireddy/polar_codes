function out = cyclic_shift(del, mat)
    out = circshift(mat,-1*del);
end