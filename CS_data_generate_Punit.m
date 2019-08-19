function [data_matrix] = CS_data_generate_Punit(mean_c,var_c,total_no_of_points,dimension)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    data_matrix=zeros(total_no_of_points,dimension);
    l=1;
    
    while l<=size(data_matrix,1)
     for k=1:dimension
        data_matrix(l,k)=mean_c+var_c*randn(1);
     end
        l=l+1;
    end
end

