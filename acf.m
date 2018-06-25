function acf = acf(samples, m )

    samples = samples -mean(samples);
    auto= xcorr(samples,'coeff');
    l = (length(auto)+1)/2;
    acf = auto(l:l+m);
end

