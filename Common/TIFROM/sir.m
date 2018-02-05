function [SIR_in_dB,SIR_out_dB,SIRI_dB]=sir(contribution_j,observation_i,output,theo_out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Matthieu Puigt   @ : Matthieu.Puigt@ast.obs-mip.fr
%
% Date : 09/06/2008
%
% Goal : Computation of the SIR


P_signal=mean( contribution_j.^2 );
P_interf=mean( (observation_i-contribution_j).^2 );

SIR_in=P_signal/P_interf;
SIR_in_dB=10*log10(SIR_in);
    
P_theo_out=mean((theo_out).^2);
P_interf_out=mean((output-theo_out).^2);

SIR_out=P_theo_out/P_interf_out;
SIR_out_dB=10*log10(SIR_out);

SIRI=SIR_out/SIR_in;
SIRI_dB=10*log10(SIRI);