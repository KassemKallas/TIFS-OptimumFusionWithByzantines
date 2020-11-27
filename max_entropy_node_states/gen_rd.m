function [sample]= gen_rd(N, nb_samples)

% The required inputes for this function are : N:the network size,
% nb_samples: how many samples from the specified distribution we want to
% generate.


%K.K.: This funtion generates a random integer according to weighted
%discrete distribution
%where each possible outcome may have different percentage to come out.


% This check if N is odd or even just to count for the cases when N/2 -1 is
% not an integer
if mod(N,2)
    % it is odd
a = 0:floor(N/2);     
else 
    % it is even
a = 0:((N/2)-1);   
end

for cf = 0:1:length(a)-1
all_nb_config(cf+1) = nchoosek(N,cf);
end
all_config = sum(all_nb_config);

for ix = 0:length(a)-1 % to generate a weight for each possibility 

    w(ix+1) = nchoosek(N,ix) / (all_config); % this generate the distribution, then for each possible Byzantine nb (0->N/2 -1)
    %gives its correspondant weight

end
sample = a( sum( bsxfun(@ge, rand(nb_samples,1), cumsum(w./sum(w))), 2) + 1 ); % this looks at the vector w and return a vector of size nb_samples
% where these values are generated according to weights in w.


% if you want to see the wieghts assigned to each possibility change the
% [sample] to [sample,w]
% then to test this function, generate in the command window sufficient nb
% of samples by this function and then by using tabulate(samples) you will
% see the percentage of each outcome.


%============================Example ====================================
% change [sample] to [sample,w]
% In the command window
% [s w] = gen_rd(10, 5000); % 5000 is the nb. of samples you can set to
% what you want, 1 till millions
% now:
% tabulate(s)

%========================================================================

end