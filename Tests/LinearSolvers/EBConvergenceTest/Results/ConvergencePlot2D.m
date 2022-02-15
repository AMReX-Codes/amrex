% Matlab/Octave script for illustrating Convergence of the EB Elliptic Solver. 
% Reads in the .mat files generated from the the MultiFabs from the Solver
% Generates Error Scaling by testing simulation against phi_true = x/sqrt(x^2 + y^2) 
% Written by Steven Reeves, August 6th 2018
% Center for Computational Science and Engineering
% Lawrence Berkeley National Laboratory 

    num = [5 6 7 8 9 10]; %actually log2 num
    [fab32  ,  dim32] = mfread2( 'phi-0_32.mat' ); 
    [fab64  ,  dim64] = mfread2( 'phi-0_64.mat' ); 
    [fab128 , dim128] = mfread2('phi-0_128.mat' ); 
    [fab256 , dim256] = mfread2('phi-0_256.mat' );
    [fab512 , dim512] = mfread2('phi-0_512.mat' ); 
    [fab1024,dim1024] = mfread2('phi-0_1024.mat'); 
   
    x32 = zeros(1,dim32(1)); 
    x32(2:end-1) = ([0:31] + 0.5)*1/32; 
    x32(end) = 1; 
    y32 = x32; 
    [X32, Y32] = meshgrid(x32,y32); 
    f32 = (X32 - 0.5)./sqrt((X32 - 0.5).^2 + (Y32 - 0.5).^2);
    for i = 1:dim32(1)
        for j = 1:dim32(2)
            if(fab32(i,j) == 0) %Detects covered cells
                f32(i,j) = 0; 
            end
        end
    end
    Err32 = abs(f32 - fab32);
    clear X32 Y32 x32 y32; 
    
    x64 = zeros(1,dim64(1)); 
    x64(2:end-1) = ([0:63] + 0.5)*1/64; 
    x64(end) = 1; 
    y64 = x64; 
    [X64, Y64] = meshgrid(x64, y64); 
    f64 = (X64 - 0.5)./sqrt((X64 - 0.5).^2 + (Y64 - 0.5).^2); 
    for i = 1:dim64(1)
        for j = 1:dim64(2)
            if(fab64(i,j) == 0) 
                f64(i,j) = 0; 
            end
        end
    end

    Err64 = abs(f64 - fab64); 
    clear X64 Y64 x64 y64;

    x128 = zeros(1,dim128(1));
    x128(2:end-1) = ([0:127] + 0.5)*1/128; 
    x128(end) = 1; 
    y128 = x128; 
    [X128, Y128] = meshgrid(x128, y128); 
    f128 = (X128 - 0.5)./sqrt((X128 - 0.5).^2 + (Y128 - 0.5).^2); 
    for i = 1:dim128(1)
        for j = 1:dim128(2)
            if(fab128(i,j) == 0) 
                f128(i,j) = 0; 
            end
        end
    end

    Err128 = abs(f128 - fab128); 
    clear X128 Y128 x128 y128;

    x256 = zeros(1,dim256(1)); 
    x256(2:end-1) = ([0:255] + 0.5)*1/256; 
    x256(end) = 1; 
    y256 = x256; 
    [X256, Y256] = meshgrid(x256, y256); 
    f256 = (X256 - 0.5)./sqrt((X256 - 0.5).^2 + (Y256 - 0.5).^2); 
    for i = 1:dim256(1)
        for j = 1:dim256(2)
            if(fab256(i,j) == 0) 
                f256(i,j) = 0; 
            end
        end
    end

    Err256 = abs(f256 - fab256); 
    clear X256 Y256 x256 y256;

    x512 = zeros(1,dim512(1)); 
    x512(2:end-1) = ([0:511] + 0.5)*1/512; 
    x512(end) = 1; 
    y512 = x512; 
    [X512, Y512] = meshgrid(x512, y512); 
    f512 = (X512 - 0.5)./sqrt((X512 - 0.5).^2 + (Y512 - 0.5).^2); 
    for i = 1:dim512(1)
        for j = 1:dim512(2)
            if(fab512(i,j) == 0) 
                f512(i,j) = 0; 
            end
        end
    end

    Err512 = abs(f512 - fab512); 
    clear X512 Y512 x512 y512;

    x1024 = zeros(1,dim1024(1)); 
    x1024(2:end-1) = ([0:1023] + 0.5)*1/1024; 
    x1024(end) = 1; 
    y1024 = x1024; 
    [X1024, Y1024] = meshgrid(x1024, y1024); 
    f1024 = (X1024 - 0.5)./sqrt((X1024 - 0.5).^2 + (Y1024 - 0.5).^2); 
    for i = 1:dim1024(1)
        for j = 1:dim1024(2)
            if(fab1024(i,j) == 0) 
                f1024(i,j) = 0; 
            end
        end
    end

    Err1024 = abs(f1024 - fab1024); 
    clear X1024 Y1024 x1024 y1024;

    err1(1) = sum(sum(Err32))/(32*32); 
    err1(2) = sum(sum(Err64))/(64*64); 
    err1(3) = sum(sum(Err128))/(128*128); 
    err1(4) = sum(sum(Err256))/(256*256); 
    err1(5) = sum(sum(Err512))/(512*512); 
    err1(6) = sum(sum(Err1024))/(1024*1024);
    
    err2(1) = sqrt(sum(sum((f32-fab32).^2))/(32*32));
    err2(2) = sqrt(sum(sum((f64-fab64).^2))/(64*64)); 
    err2(3) = sqrt(sum(sum((f128-fab128).^2))/(128*128)); 
    err2(4) = sqrt(sum(sum((f256-fab256).^2))/(256*256)); 
    err2(5) = sqrt(sum(sum((f512-fab512).^2))/(512*512)); 
    err2(6) = sqrt(sum(sum((f1024-fab1024).^2))/(1024*1024));

    errm(1) = max(max(Err32)); 
    errm(2) = max(max(Err64)); 
    errm(3) = max(max(Err128)); 
    errm(4) = max(max(Err256)); 
    errm(5) = max(max(Err512)); 
    errm(6) = max(max(Err1024)); 

    err1 = log2(err1); 
    err2 = log2(err2); 
    errm = log2(errm); 
    
    P1 = polyfit(num, err1, 1); 
    P2 = polyfit(num, err2, 1);
    Pm = polyfit(num, errm, 1); 
    scale1 = P1(1)*num + P1(2); 
    scale2 = P2(1)*num + P2(2);
    scalem = Pm(1)*num + Pm(2);
    figure(1) 
    plot(num, err1, '.', 'markersize', 10)
    hold on 
    plot(num, scale1, '-', 'linewidth', 2)
    plot(num, err2, 'k.', 'markersize', 10)
    plot(num, scale2, 'k-', 'linewidth', 2)
    plot(num, errm, 'r.', 'markersize', 10)
    plot(num, scalem, 'r-', 'linewidth', 2)
    xlabel('Log2(N)', 'fontsize', 14); 
    ylabel('Log2(Err)','fontsize', 14); 
    legend('1-norm', num2str(P1(1)), '2-norm', num2str(P2(1)), 'max-norm', num2str(Pm(1)))
    saveas(1,"EBConvergence.png"); 

            
