function [B, coord, compt, Ns_found, error, X, f, t]=li_tifrom_sym(x, N_sources, nb_samp_in_win, overlap, nb_win, fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Matthieu Puigt   @ : Matthieu.Puigt@ast.obs-mip.fr
%
% Date : 09/06/2008
%
% Goal : 
% This method is the "symmetrical" LI-TIFROM method resp. proposed and 
% described in : 
%
% [1] Y. Deville, M. Puigt, B. Albouy: Time-frequency blind signal
% separation: extended methods, performance evaluation for speech sources,
% Proceedings of the International Joint Conference on Neural Networks
% (IJCNN 2004), IEEE Catalog Number: 04CH37541C, ISBN: 0-7803-8360-5, pp.
% 255-260, Budapest, Hungary, July 25-29, 2004.
%
% [2] M. Puigt, Y. Deville: Time-frequency ratio-based blind separation 
% methods for attenuated and time-delayed sources, Mechanical Systems and
% Signal Processing, vol. 19, issue 6, pp. 1348-1379, November 2005.
%
% This is an improved version of the LI-TIFROM method resp. proposed and 
% detailed in:
%
% [3] F. Abrard, Y. Deville: A time-frequency blind signal separation
% method applicable to underdetermined mixtures of dependent sources,
% Signal Processing, vol. 85, issue 7, pp. 1389-1403, July 2005.
%
% [4] Y. Deville, M. Puigt, Temporal and time-frequency correlation-based
% blind source separation methods. Part I : Determined and underdetermined 
% linear instantaneous mixtures, Signal Processing, Volume 87, Issue 3, pp.
% 374-407, March 2007. 
% 
%
% !!!! You need the Matlab Signal Processing toolboxs to run this
% algorithm. 
%
%
% The input data are :
%    'x' which is the 2D matrix of the observations:
%       -> the first indice corresponds to the index of the source,
%       -> the second one corresponds to the temporal sample of the
%       sources.
%       
%    'N_sources' which is the number of sources. 
%    'nb_samp_in_win' which is the number of samples per STFT.
%    'overlap' which is the temporal overlap between two adjacent STFTs.
%    'nb_win' which is the number of adjacent STFTs in the analysis zones.
%    'fs' which is the sample frequency.
%
% The output data are:
%    'B' which is the estimated mixing matrix.
%    'coord' which is the 2D matrix of the coordinates of the TF zones used
%    to find the columns of B.
%    'compt' which indicates the number of TF zones used to estimate B.
%    'Ns_found' which corresponds to the number of sources which have been 
%    found.
%    'error' which indicates the reason why the program stops:
%       -> if error==0, the program found an estimate of the mixing matrix
%       (Ns_found==N_sources),
%       -> if error==1, the variance of the TF ratios is too high: we
%       assume that there is no TF single-source zone anymore,
%       -> if error==2, we have explored all the TF zones but we did not
%       find all the columns of B (one source did not occur alone).
%
%    'X' is the 3D matrix of STFTs (see the Matlab specgram's help
%    to know more about it). 
%    'f' is the vector of frequencies at which the STFTs are computed.
%    't' is the vector of time samples at which the STFTs are computed.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Stage 0 : Implementation of parameters in STFTs' computation      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6, fs = 1; end
if nargin < 5, nb_win = 10; end
if nargin < 4, overlap = 0.75; end
if nargin < 3, nb_samp_in_win = 128; end

P_obs=length(x(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 Stage 1 : Computation of the STFTs                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Computation of the number of samples which overlap in STFTs
nb_overlap= round(overlap * nb_samp_in_win);

% Computation of the STFTs
for i=1:P_obs
    [X(i,:,:), f, t]=specgram(x(i,:), nb_samp_in_win, fs, nb_samp_in_win, nb_overlap);
end % for i=1:P_obs


ratio_nb_col=length(X(1,1,:));
ratio_nb_row=length(X(1,:,1));

% Computation of the ratios of STFTs.
% In [3,4], we only compute one ratio (resp. ratio_inv in [3] and
% ratio in [4]). Here, we deal with both in order to have better
% performance (see [1,2] for explanation).
ratio=zeros(P_obs,ratio_nb_row,ratio_nb_col);
ratio_inv=zeros(P_obs,ratio_nb_row,ratio_nb_col);
for i=2:P_obs
    ratio(i,:,:)=X(i,:,:)./X(1,:,:);
    ratio_inv(i,:,:)=X(1,:,:)./X(i,:,:); % Inverse ratio
end % for i=1:P_obs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Stage 2 : Detection of TF single-source zones and identification      %
%                            of the columns of B                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization of the analysis zones  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_zones=2*floor(ratio_nb_col/nb_win)-1;


% Computation of the variance of the ratios
variance = zeros(P_obs,ratio_nb_row,nb_zones);
variance_inv = zeros(P_obs,ratio_nb_row,nb_zones);
for i=1:ratio_nb_row
    for j=1:nb_zones
        for k=2:P_obs
            
            variance(k,i,j)=var( ratio(k,i,1+round(nb_win/2)*(j-1):round(nb_win/2)*(j+1)) );
            variance_inv(k,i,j)=var( ratio_inv(k,i,1+round(nb_win/2)*(j-1):round(nb_win/2)*(j+1)) );
  
        end % for k=2:P_obs
    end % for j=1:nb_zones
end % for i=1:ratio_nb_row

% We do not study the variance of ratios X_1 / X_1 !
variance(1,:,:)=1e9;
variance_inv(1,:,:)=1e9;

% We deal with the mean of the variance (see [1,2] for explanations)
if P_obs>2  
    
    mean_variance(:,:)=mean(variance(2:end,:,:));
    mean_variance_inv(:,:)=mean(variance_inv(2:end,:,:));
    
else mean_variance(:,:)=variance(2,:,:);
     mean_variance_inv(:,:)=variance_inv(2,:,:);
    
end % if P_obs>2 

Ns_found=0; % We did not find any column of B for the moment.

% test_stop is a Boolean which indicates when stopping the algorithm.
test_stop=1;

% var_thresh is the maximum authorized variance for a TF zone to be
% "really" single-source. It allows us to speed up the computations, by
% occulting some TF analysis zones. In our first papers, the value of this
% criterion was set to 1.5e-2. We now omit it by puting a high value!
var_thresh=1e9;

% In this method, we compute the distance between a new column of B and
% each previous found one. If this distance is above a user-defined 
% threshold, we keep this new column. We call this threshold c_thresh.
c_thresh=0.15;

compt=0;
error=0;

while test_stop==1

    [m1, index1] = min(mean_variance);
    [p1, k1]=min(m1);
    [m2, index2] = min(mean_variance_inv);
    [p2, k2]=min(m2);
    % p1 (resp. p2) corresponds to the minimum of the vector m1 (resp. m2),
    % i.e. the lowest value of the matrix mean_variance (resp. 
    % mean_variance_inv). 
    % k1 (resp. k2) is the temporal index of the minimum of mean_variance
    % (resp. mean_variance_inv).
    % index1(k1) (resp. index2(k2)) is the frequency index of the minimum
    % of mean_variance (resp. mean_variance_inv).
    
    
    [p,method]=min([p1,p2]);
    test_method=method-1;
    % We look for the minimum between the minimum of the variance
    % of the ratio X(i,:,:)./X(1,:,:) and the minimum of the variance of
    % the inverse ratio X(1,:,:)./X(i,:,:). See [1,2] for explanation.
    
    B(1,Ns_found+1)=1;
    
    if test_method==0 
        
        % This is the case p1 == p.
        coord(Ns_found+1,1)=index1(k1);
        coord(Ns_found+1,2)=k1;
        
        for i=2:P_obs
            % The elements of the possible new column of B are equal to the
            % ratio X(i,:,:)./X(1,:,:)
            B(i,Ns_found+1) = real(mean( ratio(i, index1(k1),1+round(nb_win/2)*(k1-1):1+round(nb_win/2)*(k1-1)+nb_win-1) ) );       
        end % for i=2:P_obs
        
    else
        % This is the case p == p2.
        coord(Ns_found+1,1)=index2(k2);
        coord(Ns_found+1,2)=k2;
        
         for i=2:P_obs
            % The elements of the possible new column of B are equal to the
            % inverse ratio X(1,:,:)./X(i,:,:)
            B(i,Ns_found+1) = 1/real(mean( ratio_inv(i, index2(k2),1+round(nb_win/2)*(k2-1):1+round(nb_win/2)*(k2-1)+nb_win-1) ) );
       
         end %fin for i=2:P_obs       
         
    end % if test_method==0
        
           
    compt=compt+1;
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          We have to know if this column is a new column of B        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    previous_found=0;
    counter = Ns_found;
    
    while (counter>0) & (previous_found==0) 
        
        if test_method==0
            dist=abs( B( 2:end, Ns_found+1)- B(2:end , counter) );
            
        else
            dist=abs( 1./B( 2:end, Ns_found+1)- 1./B(2:end , counter) );
            
        end % if test_method==0
            
        if max(dist)< c_thresh
           previous_found =1;    % This tentative column is a previous 
                                 % found column.
           
        else counter=counter-1; % Comparison of this tentative column to 
                                % another previous found column.
                     
        end % if max(dist)< c_thresh
         
        
        
    end % while (counter>0) & (previous_found==0)
             
 
    % Comparison of the variance of this TF zone and var_thresh
  
    if p > var_thresh  % The variance is too high
        error = 1;
        fprintf('The variance is too high.\nThere is no more TF single-source analysis zones.\n')
        fprintf('We just found %g column(s) of the unknown mixing matrix',Ns_found)
        B=B(:,1:Ns_found);
        test_stop = 0;
      
    else   if (test_stop==1) & (previous_found==0)
               Ns_found=Ns_found+1;                    
           end % if (test_stop==1) & (previous_found==0)
                
           if Ns_found == N_sources % We have find N_sources columns of B !
              test_stop =0;
              
           else if ( m1(k1)== max( max( mean_variance) ) ) & ( m2(k2)== max( max( mean_variance_inv) ) )
                     % We did not find N_sources columns of B and we have studied
                     % all TF analysis zones.
                     
                     error=2;
                     fprintf('We have studied all the TF analysis zones.\n')
                     fprintf('We just found %g column(s) of the unknown mixing matrix',Ns_found)
                     B=B(:,1:Ns_found);                     
                     test_stop = 0;
                        
                else if test_method==0
                          mean_variance(index1(k1), k1) = max(max(mean_variance));
                          mean_variance_inv(index1(k1), k1) = max(max(mean_variance));
                          
                      else mean_variance_inv(index2(k2), k2) = max(max(mean_variance_inv));
                           mean_variance(index2(k2), k2) = max(max(mean_variance_inv));
                          
                      end % if test_method==0
                      
                end % if m(k)== max(max(moy_variance)))         
                      
           end % if Ns_found == N_sources                

    end % if p > var_thresh          
            
end % while test_stop==1
        
