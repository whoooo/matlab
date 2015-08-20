function n = normalize(input, norm_factor)

n = input/max(abs(input)).*norm_factor;



