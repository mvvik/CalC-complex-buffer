%
%   Read and plot 4D field dump produced by "plot binary Field fileName"
%          CalC version 7.4 * Victor Matveev * April 18, 2013
%-------------------------------------------------------------------------


f  = fopen('fileName.dat', 'rb'); % Read the dump file

G  = fread(f, 1, 'int');          % Read geometry (lowest 2 bits = # of dimensions)

N_r  = fread(f, 1, 'int');        % Read number of radial nodes in 4D dump
N_th = fread(f, 1, 'int');        % Read number of angular phi nodes in 4D dump

R     = fread(f, N_r, 'double');    % Radial coordinates
Theta = fread(f, N_th, 'double');   % Angular coordinates (theta = 0)

[r, theta] = meshgrid(R, Theta);

N = N_r * N_th;               % Number of nodes (N_th = 1)
state = 0;

figure(1);

while ~state
    t = fread(f, 1, 'double'); 
    if feof(f) 
        return;  
    end;
    A = fread(f, N, 'double');  % Read the data
    B = reshape(A, N_r, N_th); % Reshape into 3D array
    Ca = squeeze(B(:,:,1))'; % Reduce to a 2D array
    [Z X Ca] = pol2cart(theta, r, Ca);
    contourf(X, Z, Ca);          
    title([' time = ', num2str(t), ' ms' ]);
    xlabel('X (\mum)'); ylabel('Z (\mum)');
    colorbar;
    drawnow;
    pause(0.2);
    state = feof(f);
end

fclose(f);                   % close the file
