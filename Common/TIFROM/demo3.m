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


% We mix four sources
N=4;


s(1,:)=wavread('source1.wav')';
s(2,:)=wavread('source2.wav')';
s(3,:)=wavread('source3.wav')';
s(4,:)=wavread('source4.wav')';

% fs is the sample frequency rate
fs=20000;

% Choice of the mixing matrix
alpha=0.9;
A=[1, alpha, alpha^2, alpha^3; ...
   alpha, 1, alpha, alpha^2; ...
   alpha^2, alpha, 1, alpha; ...
   alpha^3, alpha^2, alpha, 1];

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

% Theoretical output
theo_out=inv(B_theo)*x;

% Performance achieved by the LI-TIFROM method

% Computation of the SIR
if Ns_found1==N

    y1=inv(B1)*x;

	SIR_in_dB=[zeros(1,N)];
    SIR_out_dB_basic=[zeros(1,N)];
    SIRI_dB_basic=[zeros(1,N)];  
    
    for index=1:N
        
        [SIR_in_dB1(index),SIR_out_dB1(index),SIRI_dB1(index)]=sir(xc(:,index,index),x(index,:)',y1(1,:),theo_out(index,:));
        [SIR_in_dB2(index),SIR_out_dB2(index),SIRI_dB2(index)]=sir(xc(:,index,index),x(index,:)',y1(2,:),theo_out(index,:));
        [SIR_in_dB3(index),SIR_out_dB3(index),SIRI_dB3(index)]=sir(xc(:,index,index),x(index,:)',y1(3,:),theo_out(index,:));
        [SIR_in_dB4(index),SIR_out_dB4(index),SIRI_dB4(index)]=sir(xc(:,index,index),x(index,:)',y1(4,:),theo_out(index,:));        
        SIRI_dB_basic(index)=max([SIRI_dB1(index),SIRI_dB2(index),SIRI_dB3(index),SIRI_dB4(index)]);
        SIR_out_dB_basic(index)=max([SIR_out_dB1(index),SIR_out_dB2(index),SIR_out_dB3(index),SIR_out_dB4(index)]);
        SIR_in_dB(index)=max([SIR_in_dB1(index),SIR_in_dB2(index),SIR_in_dB3(index),SIR_in_dB4(index)]);
        
    end % for index=1:N
        
end % if Ns_found1==N

if Ns_found2==N

    y2=inv(B2)*x;
    
    SIR_out_dB_improved=[zeros(1,N)];
    SIRI_dB_improved=[zeros(1,N)];
    
    
    % Theoretical output
    theo_out=inv(B_theo)*x;

            
    for index=1:N
        
        [SIR_in_dB1(index),SIR_out_dB1(index),SIRI_dB1(index)]=sir(xc(:,index,index),x(index,:)',y2(1,:),theo_out(index,:));
        [SIR_in_dB2(index),SIR_out_dB2(index),SIRI_dB2(index)]=sir(xc(:,index,index),x(index,:)',y2(2,:),theo_out(index,:));
        [SIR_in_dB3(index),SIR_out_dB3(index),SIRI_dB3(index)]=sir(xc(:,index,index),x(index,:)',y2(3,:),theo_out(index,:));
        [SIR_in_dB4(index),SIR_out_dB4(index),SIRI_dB4(index)]=sir(xc(:,index,index),x(index,:)',y2(4,:),theo_out(index,:));        
        SIRI_dB_improved(index)=max([SIRI_dB1(index),SIRI_dB2(index),SIRI_dB3(index),SIRI_dB4(index)]);
        SIR_out_dB_improved(index)=max([SIR_out_dB1(index),SIR_out_dB2(index),SIR_out_dB3(index),SIR_out_dB4(index)]);
        
    end % for index=1:N
        
end % if Ns_found2==N

% Performance achieved by the LI-TIFROM method
fprintf('-----\nPerformance achieved by the LI-TIFROM methods :\n')
fprintf('SIR^in : %.1f dB\n', mean(SIR_in_dB))
if Ns_found1==N
    fprintf('-----\n"Asymmetrical" LI-TIFROM method :\n')
    fprintf('SIR^out : %.1f dB\n', mean(SIR_out_dB_basic))
    fprintf('SIRI : %.1f dB\n', mean(SIRI_dB_basic))
else
    fprintf('The "Asymmetrical" method did not estimate the mixing matrix')
end % if Ns_found1==N

if Ns_found2==N
    fprintf('-----\n"Symmetrical" LI-TIFROM method :\n')
    fprintf('SIR^out : %.1f dB\n', mean(SIR_out_dB_improved))
    fprintf('SIRI : %.1f dB\n', mean(SIRI_dB_improved))
else
    fprintf('The "symmetrical" method did not estimate the mixing matrix')
end % if Ns_found2==N