function [MD, match_order, new_ca2] = calc_max_dice_adj(ca1, ca2)
%Function to match two cluster solutions by maximizing the average dice
%coefficient between pairs of clusters. This uses the function perms.m to
%loop through all possible permutations of the ca2 labels, so it is not
%really practical if you are trying to match more than 11 clusters.
%Usage: [MD matched_d relabeled_ca2] = calc_max_dice_bs(ca1, ca2)
%where
%       ca1 & ca2 - vectors of cluster labels to match & compare
%       MD - maximum average dice coefficient between matched labels
%       matched_d - dice coefficient between each pair of matched clusters
%       relabeled_ca2 - vector of cluster labels for input 2 that match
%       labels in ca1

assert(length(ca1) == length(ca2));
if(size(ca1,2) == length(ca1))
    ca1 = ca1';
end
if(size(ca2,2) == length(ca2))
    ca2 = ca2';
end
n = length(ca1);

ca1_unique = nonzeros(unique(ca1));
ca2_unique = nonzeros(unique(ca2));

% check the integrity of ca2
if length(ca1_unique) ~= length(ca2_unique)
    error('The clustering ca2 is not consistent with ca1.');
end;

c = length(ca1_unique);

% distribution of ca2
M2 = double(repmat(ca2,1,c) == repmat(ca2_unique',n,1));

all_possible = perms(1:c);
all_dice = zeros(size(ca1));
new_ca2 = zeros(size(ca1));

for ip = 1:length(all_possible)
    
    new_ca2 = M2*all_possible(ip, :)';
    all_dice(ip) = nnz(new_ca2 == ca1)/nnz(ca1);

end %ip

[MD, MDI] = max(all_dice);
match_order = all_possible(MDI, :)';
new_ca2 = M2*match_order;






