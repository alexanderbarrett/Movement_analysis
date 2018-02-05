%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Matthieu Puigt   @ : Matthieu.Puigt@ast.obs-mip.fr
%
% Date : 09/06/2008
%
% Goal : 
% Demonstration program of the LI-TIFROM methods.

clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%                LINEAR INSTANTANEOUS MIXTURE PROCESS                       %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% We mix three sources
N=3;


s(1,:)=wavread('source1.wav')';
s(2,:)=wavread('source2.wav')';
s(3,:)=wavread('source3.wav')';

% fs is the sample frequency rate
fs=20000;

% Choice of the mixing matrix
A=[1, 0.9, 0.7; 0.8, 1, 1];

x=A*s;

P_obs=length(x(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%                                                        % 
%                   SEPARATING SYSTEM                    %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choice of the LI-TIFROM parameters
nb_win =10;%
nb_samp_in_win=128;%
overlap=0.75;%

[B1, coord1, compt1, Ns_found1, error1, X, f, t]=li_tifrom(x, N, nb_samp_in_win, overlap, nb_win, fs);
[B2, coord2, compt2, Ns_found2, error2, X, f, t]=li_tifrom_sym(x, N, nb_samp_in_win, overlap, nb_win, fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%                                                        % 
%              SEPARATION PERFORMANCE                    %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of the true contributions
[xc]=contributions(A,s);

% Computation of the theoretical matrix to be identified
B_theo=zeros(size(A));
for i=P_obs:-1:1
    B_theo(i,:)=A(i,:)./A(1,:);
end % for i=N:-1:1

B_theo1=[B_theo(:,1),B_theo(:,2),B_theo(:,3)];
B_theo2=[B_theo(:,1),B_theo(:,3),B_theo(:,2)];
B_theo3=[B_theo(:,2),B_theo(:,1),B_theo(:,3)];
B_theo4=[B_theo(:,2),B_theo(:,3),B_theo(:,1)];
B_theo5=[B_theo(:,3),B_theo(:,1),B_theo(:,2)];
B_theo6=[B_theo(:,3),B_theo(:,2),B_theo(:,1)];

% Performance achieved by the LI-TIFROM methods

% Computations of the Frobenius norm of the difference between theoretical
% and actual mixing matrices

norm1=min([ norm(B1-B_theo1, 'fro'), norm(B1-B_theo2, 'fro'), norm(B1-B_theo3, 'fro'), norm(B1-B_theo4, 'fro'), norm(B1-B_theo5, 'fro'), norm(B1-B_theo6, 'fro')]);
norm2=min([ norm(B2-B_theo1, 'fro'), norm(B2-B_theo2, 'fro'), norm(B2-B_theo3, 'fro'), norm(B2-B_theo4, 'fro'), norm(B2-B_theo5, 'fro'), norm(B2-B_theo6, 'fro')]);

fprintf('-----\nPerformance achieved by the LI-TIFROM methods :\n')
if Ns_found1==N
    fprintf('-----\n"Asymmetrical" LI-TIFROM method :\n')
    fprintf('Frobenius norm : %.1e\n', norm1)
else
    fprintf('The "Asymmetrical" method did not estimate the mixing matrix\n')
end % if Ns_found1==N

if Ns_found2==N
    fprintf('-----\n"Symmetrical" LI-TIFROM method :\n')
    fprintf('Frobenius norm : %.1e\n', norm2)
else
    fprintf('The "symmetrical" method did not estimate the mixing matrix\n')
end % if Ns_found2==N