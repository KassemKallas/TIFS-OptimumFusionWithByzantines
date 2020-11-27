function [all_perms] = permute_matrix(N,M)

%This function will take all possible index permutations of a matrix
%binom(n,k) and then after it will add to each permutation the set of
%missing indexes in the required order. that is because nchoosek function
%is not able to completely do that but, we can exploit its retrievals and w
% we can add for each permutation the set of indexes missing, then we can
% take each row of indeces and re-arrange the matrix accordingly to have a
% matrix permutation.



all_idx = 1:1:N;
bino = nchoosek(1:1:N,M);

test = zeros(size(bino,1),N);
for x = 1:length(bino)

missing(x,:) = setdiff(all_idx,bino(x,:));

end

for i=1:size(test,1) 

    row = bino(i,:);
    
     for j=1:size(test,2)
         
         if j <= size(bino,2)
             
             test(i,j) = row(j);
            
         else
             
             test(i,j) = missing(i,j-M);
             
         end
    
     end 
end

all_perms = test;

end