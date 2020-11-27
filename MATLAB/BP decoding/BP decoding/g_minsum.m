function g = g_minsum(beta_bit, LLR1, LLR2)

    g = ((1 - 2*beta_bit) .* LLR1) + LLR2;
    
end