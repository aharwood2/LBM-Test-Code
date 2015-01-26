function f = Stream(f,N,M,c)

% Create temporary grid to prevent overwriting useful populations
fnew = zeros(size(f));

% Stream one lattice site at a time
for i = 1:N
    for j = 1:M

        for n = 1:8
            
            % Compute destination coordinates
            dest_x = 1 + mod(i-1+c(1,n)+N,N);
            dest_y = 1 + mod(j-1+c(2,n)+M,M);
            
            % Stream values
            fnew(dest_y,dest_x,n) = f(j,i,n);


        end
    end
end

fnew(:,:,9) = f(:,:,9); % Rest particle distributions
f = fnew;

end