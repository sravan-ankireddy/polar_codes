function f = f_minsum(LLR1, LLR2)

	f = sign(LLR1).*sign(LLR2).*min(abs(LLR1),abs(LLR2));

end