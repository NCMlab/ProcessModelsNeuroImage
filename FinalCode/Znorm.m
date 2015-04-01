function out = Znorm(in)
out = (in - mean(in))./std(in);