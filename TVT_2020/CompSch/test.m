D = 1:0.1:30;
N = zeros(1,length(D));
for i = 1:length(D)
   N(i) =  MY_number_of_antennas_FDLens(D(i),pi,pi);
end