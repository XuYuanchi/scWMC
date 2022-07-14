function[normed_data] = process(data)

libsize  = sum(data,2);
normed_data = bsxfun(@rdivide, data, libsize) * median(libsize);
