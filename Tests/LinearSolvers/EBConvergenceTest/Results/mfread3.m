% Function for reading in multifab data, assumes single grid
% Contains ghost cells 
% Written by Steven Reeves, August 6th 2018
% Center for Computational Science and Engineering 
% Lawrence Berkeley National Laboratory

function [fab, dims] = mfread3(infile)
    fopen(infile,'rb'); 
    dat = fread(infile, 'double'); 
    lo = int32([dat(1) dat(3) dat(5)]);
    hi = int32([dat(2) dat(4) dat(6)]); 
    dims = [lo(1)+1 + hi(1) + 2, lo(2)+1 + hi(2) + 2, lo(3)+1 + hi(3) + 2]; %gives the size of the array with 2 ghost cells
    dat(1:6) = []; 
    fab = reshape(dat, dims(1), dims(2), dims(3));
    fclose(infile); 
end
    



