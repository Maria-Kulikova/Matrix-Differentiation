% Source code for reproducing the numerical experiment (Example 1) in 
%
%    Kulikova M.V., Tsyganova J.V. (2015) "Constructing numerically stable 
%    Kalman filter-based algorithms for gradient-based adaptive filtering", 
%    International Journal of Adaptive Control and Signal Processing, 
%    29(11):1411-1426. DOI http://dx.doi.org/10.1002/acs.2552 
%
% Authors: Maria Kulikova:  kulikova dot maria at yahoo dot com     
% ------------------------------------------------------------------- 

clc; clear all; close all; 

theta = sym('theta','real'); % parameter (symbolic)
parameters = [theta];
disp('Source code for reproducing the numerical experiment (Example 1, Table 1)'); 
disp('Given: the pre-array A in the following form:');
PreArray = [theta^5/20,   theta^4/8,  theta^3/6,  theta^3/3;
            theta^4/8,   theta^3/3, theta^2/2,  theta^2/2;
            theta^3/6,   theta^2/2, theta,     1;]
disp('and the derivative of the pre-array A (with respect to theta):');
Diff_PreArray = diff(PreArray,theta)

disp('Compute: the post-array R such that R=QA and the derivative of R at the point theta = 2');
disp('         In the post-array R only the first main 3 by 3 block is upper triangular');
disp('________________________________________________________________________________________');
disp('ALGORITHM');
disp('<press any key to run the algorithm and to see the result>');
pause

theta = 2;  % system parameter and its current value
disp('The pre-array A at the point theta = 2:');
PreArray = [theta^5/20,   theta^4/8,  theta^3/6,  theta^3/3;
            theta^4/8,   theta^3/3, theta^2/2,  theta^2/2;
            theta^3/6,   theta^2/2, theta,     1;]
disp('The derivative of the pre-array A at the point theta = 2:');
Diff_PreArray  = double(subs(Diff_PreArray,parameters,2)) % substitute the parameter values
disp('________________________________________________________________________________________');
disp('Results');
disp('<press any key to run the algorithm and to see the result>');
disp('The computed post-array R (at the point theta = 2):');
[PostArray,Diff_PostArray,orth] = Diff_QR(PreArray,3,Diff_PreArray);
disp(PostArray);
disp('The orthogonal transformation (QA=R):');
disp(orth);
disp('The computed derivative of the post-array R (at the point theta = 2):');
disp(Diff_PostArray);
disp('________________________________________________________________________________________');
disp('1. Accuracy of factorization: norm( A^T*A - R^T*R ,Inf) ');
disp(norm(PreArray'*PreArray - PostArray'*PostArray,inf));

disp('2. Accuracy of derivative computation: norm( (A^T*A)^{\prime} - (R^T*R)^{\prime} ,Inf) ');
disp(norm(PreArray'*Diff_PreArray+Diff_PreArray'*PreArray-PostArray'*Diff_PostArray-Diff_PostArray'*PostArray,inf));

