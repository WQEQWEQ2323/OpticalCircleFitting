function O_A2_Plot=F_UnZerosFilling(Matrix,Target_Length)
[a,b]=size(Matrix);
O_A2_Plot=Matrix(a/2-Target_Length/2+1:a/2+Target_Length/2,b/2-Target_Length/2+1:b/2+Target_Length/2);
end