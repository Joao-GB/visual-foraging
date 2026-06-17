function opt = dva_n_pixel_set_options(var)
    opt = struct('dist', 45, 'pixel_dva_ratio', NaN, 'width', NaN, 'res', NaN);
    optionNames = fieldnames(opt);
    for pair = reshape(var, 2, [])
        inpName = pair{1};
        if any(strcmp(inpName, optionNames))
            opt.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name', inpName);
        end
    end
end
