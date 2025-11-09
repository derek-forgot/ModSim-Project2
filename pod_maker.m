function [M] = pod_maker(n_population,pod_size)
    
    %Set these up to be used in the loop
    M = zeros(n_population, n_population);
    M_pod_n = ones(pod_size, pod_size) - eye(pod_size);
    
    for i = 1:(ceil(n_population/pod_size))+1
        slice_begin = (i - 1) * pod_size + 1;
        slice_end = slice_begin + pod_size - 1;
        if slice_end > n_population
            slice_end = n_population;
            M_pod_n = ones((slice_end-slice_begin)+1, (slice_end-slice_begin)+1) - eye((slice_end-slice_begin)+1);
        end
        M(slice_begin:slice_end, slice_begin:slice_end) = M_pod_n;
    end
end

