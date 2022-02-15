% Matlab/Octave script for illustrating Convergence of the EB Elliptic Solver. 
% Reads in the .mat files generated from the the MultiFabs from the Solver
% Generates Error Scaling by testing simulation against phi_true = x/sqrt(x^2 + y^2) 
% Written by Steven Reeves, August 6th 2018
% Center for Computational Science and Engineering
% Lawrence Berkeley National Laboratory 

    num = [4 5 6 7 8]; %actually log2 num
    [fab16  ,  dim16] = mfread3( 'phi-0_16.mat' ); 
    [fab32  ,  dim32] = mfread3( 'phi-0_32.mat' ); 
    [fab64  ,  dim64] = mfread3( 'phi-0_64.mat' ); 
    [fab128 , dim128] = mfread3('phi-0_128.mat' ); 
    [fab256 , dim256] = mfread3('phi-0_256.mat' );

    x16 = zeros(1,dim16(1)); 
    x16(2:end-1) = ([0:15] + 0.5)*1/16; 
    x16(end) = 1; 
    y16 = x16; 
    z16 = x16; 
    [X16, Y16, Z16] = ndgrid(x16, y16, z16); 
    f16 = (X16 + Z16 - 1.0)./sqrt((X16 - 0.5).^2 + (Y16 - 0.5).^2 + (Z16 - 0.5).^2);
    for i = 1:dim16(1)
        for j = 1:dim16(2)
            for k = 1:dim16(3)
              if(fab16(i, j, k) == 0) %Detects covered cells
                f16(i,j,k) = 0; 
	      end
            end
        end
    end
    Err16 = abs(f16 - fab16);
    clear X16 Y16 x16 y16 Z16 z16; 

   
    x32 = zeros(1,dim32(1)); 
    x32(2:end-1) = ([0:31] + 0.5)*1/32; 
    x32(end) = 1; 
    y32 = x32; 
    z32 = x32; 
    [X32, Y32, Z32] = ndgrid(x32, y32 ,z32); 
    f32 = (X32 + Z32 - 1.0)./sqrt((X32 - 0.5).^2 + (Y32 - 0.5).^2 + (Z32 - 0.5).^2);
    for i = 1:dim32(1)
        for j = 1:dim32(2)
	   for k = 1:dim32(3)
            if(fab32(i,j, k) == 0) %Detects covered cells
                f32(i,j,k) = 0; 
            end
           end
        end
    end
    Err32 = abs(f32 - fab32);
    clear X32 Y32 Z32 x32 y32 z32; 
    
    x64 = zeros(1,dim64(1)); 
    x64(2:end-1) = ([0:63] + 0.5)*1/64; 
    x64(end) = 1; 
    y64 = x64; 
    z64 = x64;
    [X64, Y64, Z64] = ndgrid(x64, y64, z64); 
    f64 = (X64 + Z64 - 1.0)./sqrt((X64 - 0.5).^2 + (Y64 - 0.5).^2 + (Z64 - 0.5).^2);
    for i = 1:dim64(1)
        for j = 1:dim64(2)
	  for k = 1:dim64(3)
            if(fab64(i,j,k) == 0) 
                f64(i,j,k) = 0; 
            end
          end
        end
    end

    Err64 = abs(f64 - fab64); 
    clear X64 Y64 Z64 x64 y64 z64;

    x128 = zeros(1,dim128(1));
    x128(2:end-1) = ([0:127] + 0.5)*1/128; 
    x128(end) = 1; 
    y128 = x128; 
    z128 = x128;
    [X128, Y128, Z128] = ndgrid(x128, y128, z128); 
    f128 = (X128 + Z128 - 1.0)./sqrt((X128 - 0.5).^2 + (Y128 - 0.5).^2 + (Z128 - 0.5).^2);
    for i = 1:dim128(1)
        for j = 1:dim128(2)
	  for k = 1:dim128(3)
            if(fab128(i,j,k) == 0) 
                f128(i,j,k) = 0; 
            end
          end
        end
    end

    Err128 = abs(f128 - fab128); 
    clear X128 Y128 Z128 x128 y128 z128;

    x256 = zeros(1,dim256(1)); 
    x256(2:end-1) = ([0:255] + 0.5)*1/256; 
    x256(end) = 1; 
    y256 = x256;
    z256 = x256; 
    [X256, Y256, Z256] = ndgrid(x256, y256, z256); 
    f256 = (X256 + Z256 - 1.0)./sqrt((X256 - 0.5).^2 + (Y256 - 0.5).^2 + (Z256 - 0.5).^2);
    for i = 1:dim256(1)
        for j = 1:dim256(2)
	   for k = 1:dim256(3)
            if(fab256(i,j,k) == 0) 
                f256(i,j,k) = 0; 
            end
	   end
        end
    end

    Err256 = abs(f256 - fab256); 
    clear X256 Y256 Z256 x256 y256 z256;

    err(1) = sum(sum(sum(Err16)))/(16*16*16);
    err(2) = sum(sum(sum(Err32)))/(32*32*32); 
    err(3) = sum(sum(sum(Err64)))/(64*64*64); 
    err(4) = sum(sum(sum(Err128)))/(128*128*128); 
    err(5) = sum(sum(sum(Err256)))/(256*256*256); 
    err = log2(err); 
    
    P = polyfit(num, err, 1); 
    scale = P(1)*num + P(2); 
    figure(1) 
    plot(num, err, '.', 'markersize', 10)
    hold on 
    plot(num, scale, 'k-', 'linewidth', 2)
    
    xlabel('Log2(N)', 'fontsize', 14); 
    ylabel('Log2(||sym - true||)','fontsize', 14); 
    legend('Error', 'O(dx^2)')
    saveas(1,"EBConvergence3D.png"); 

            
