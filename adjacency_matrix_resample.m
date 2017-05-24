function adjacency = adjacency_matrix_resample(ca, write_flg, outname)
%Usage
%   adjacency_matrix(ca)
%where
%   ca - cluster labels saved in memory
%   write_flg - if =1, eta_sq will be saved to a .mat
%   outfile - is the path/name of the .mat to save
%
% MBN modified February 23, 2012 - instead of reading in cluster
% information, cluster labels are passed as an input argument


adjacency = zeros(length(ca));
for ivox = 1:length(ca)
    adjacency(ivox, logical(ca == ca(ivox))) = 1;
end

%Write out adjacency matrix as .mat
if(write_flg)
    save(outname, 'adjacency');
end


        
