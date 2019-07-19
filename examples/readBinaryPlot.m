%
%   Read and plot 4D field dump produced by "plot binary Field fileName"
%          CalC version 7.4 * Victor Matveev * April 18, 2013
%
%    Easily modifiable for any geometry or cross-section (just remove
%              extra fread statements in lower dimensionality)
%-------------------------------------------------------------------------

f  = fopen('fileName.dat', 'rb'); % Read the dump file

G  = fread(f, 1, 'int');        % Read geometry (lowest 2 bits = # of dimensions)

N1 = fread(f, 1, 'int');        % Read number of x-nodes in binary plot
N2 = fread(f, 1, 'int');        % Read number of y-nodes in binary plot
N3 = fread(f, 1, 'int');        % Read number of z-nodes in binary plot

X1 = fread(f, N1, 'double');    % x-axis coordinates
X2 = fread(f, N2, 'double');    % y-axis coordinates
X3 = fread(f, N3, 'double');    % z-axis coordinates

[X Y]  = meshgrid(X1, X2);

zNode  = 2;                     % Choose your z-slice
zCoord = X3(zNode);             % z-coordinate of the z-slize

N = N1 * N2 * N3;               % Number of nodes
state = 0;

while ~state
    t = fread(f, 1, 'double');  % Read Time
    if feof(f) 
        return;  
    end;
    A = fread(f, N, 'double');  % Read the data
    B = reshape(A, N1, N2, N3); % Reshape into 3D array
    Z = squeeze(B(:,:,zNode))'; % Take the z-slice, NOTE THE TRANSPOSITION!
    contourf(X, Y, Z);          % plot the z-slice

    title([ 'Z = ', num2str(zCoord), ' \mum; time = ', num2str(t), ' ms' ]);

    xlabel('X (\mum)'); ylabel('Y (\mum)');
    colorbar;
    drawnow;
    pause(0.2);
    state = feof(f);
end

fclose(f);                   % close the file