% QL rectangular matrices differentiation 
% 
% Modern square-root Kalman filter implementations use QL factorization 
% at each iteration step of the filter: given a matrix A (pre-array), 
% assume that we apply an orthogonal operator Q to the A (i.e. QA=L) 
% so as to get the lower-triangular post-array L.
% 
% In general, the A and L has the following structure: 
% one iteration of the array SR filters has the form
%    Q [A11 A12;   =  [0   L12;
%       A21 A22]       L21 L22], 
% where Q is orthogonal such that L21 (s by s) is lower triangular matrix.
% 
% This code allows for computing the post-array L and its derivative L'_{\theta} 
% NOTE. However, we still do not know how to find the block L12. Thus, if it is 
% required, then we put NaN for these values. 
%
% References: (see Lemma 1 and Algorithm 1 from the paper)
%        Kulikova M.V., Tsyganova J.V. (2015) "Constructing numerically stable 
%        Kalman filter-based algorithms for gradient-based adaptive filtering", 
%        International Journal of Adaptive Control and Signal Processing, 
%        29(11):1411-1426. DOI http://dx.doi.org/10.1002/acs.2552 
%
% Authors: Maria Kulikova:  kulikova dot maria at yahoo dot com     
% ------------------------------------------------------------------- 
% Input:
%     Pre_array      - given pre-array A;
%     s              - size of the block (L_21) that should be triangularized;
%     Diff_pre_array - given derivative of the pre-array A'_{\theta}; 
% Output:
%     Post_array      - get post-array L (block lower triangular);
%     Diff_post_array - get the derivative of the post-array L'_{\theta};
%
% ------------------------------------------------------------------- 
function [Post_array,Diff_post_array,Orthog] = Diff_QL(Pre_array,s,Diff_pre_array)

[sk,sl] = size(Pre_array); k = sk-s; l = sl-s;       % sizes of the corresponding blocks
[Orthog,Post_array] = ql(Pre_array(k+1:end,1:s));    % Q is any orthogonal rotation that
                                                     % lower-triangularizes the first
                                                     % (s+k) by (s) block of the pre-array
if l~=0
 Post_array2 = Orthog'*Pre_array(1:end,end); % compute L12, L22 blocks
 Post_array  = [Post_array, Post_array2];    % the full post-array
end;

Diff_post_array = [];
if nargin>2
  L21 = Post_array(k+1:k+s,1:s);         % blocks of the post-array
  L12 = Post_array(1:k,s+1:s+l);         % blocks of the post-array
  L22 = Post_array(k+1:s+k,s+1:s+l);     % blocks of the post-array
   
  Q_applied_DiffA = Orthog'*Diff_pre_array; % apply Q to the pre-array derivatives
  XX = Q_applied_DiffA(1:k,1:s);            % notations
  YY = Q_applied_DiffA(k+1:s+k,1:s);
  VV = Q_applied_DiffA(k+1:s+k,s+1:s+l);
  
  split_product = YY/L21; % compute the matrix product
                          % Note, MatLab treats correctly the empty blocks
   
  factor_U      = triu(split_product,1);      % U part
  factor_L      = tril(split_product,-1);     % L part
  factor_D      = diag(diag(split_product));  % D part
                                          
  Diff_L21 = (factor_L+factor_D+factor_U')*L21;
  Diff_L22 = (factor_U'-factor_U)*L22+L21'\XX'*L12+VV;
  Diff_L11 = zeros(k,k);     
  Diff_L12 = ones(k,l)*NaN;  % The proposed method does not allow for 
                             % computing L12'_{\theta}. We set this block to NaN. 
                             % Fortunately, in the most filtering methods
                             % this block is of no interest
    
 Diff_post_array = [Diff_L11, Diff_L12; Diff_L21, Diff_L22];
 
end;
end

%%%%%%     additianal functions    %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q,L] = ql(A);
% Compute QR and manipulate to get QL
[Q,R] = qr(fliplr(A));
Q = fliplr(Q);
L = rot90(R,2);
end

