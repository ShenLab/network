function demoMICRF()

% add MICRF and UGM functions into current work space
addpath(genpath(pwd))

% MICRF inputs and outputs information
help MICRF

% MICRF with one node score file
nodefile='data/Inputs/hotnet_inputmeta.txt';
out=MICRF(nodefile); 

% MICRF with one selected network
netfile='HPRD';
out=MICRF(nodefile,netfile); 

% MICRF with users given network
netfile='data/Inputs/Co_PrePPI_3.txt';
out=MICRF(nodefile,netfile); 

% MICRF with output file
netfile='CoPrePPI';
outputfile='MICRFtest.txt';
out=MICRF(nodefile,netfile,outputfile); 

% MICRF with different non-risk prior
pi0=0.96;
[out,w]=MICRF(nodefile,netfile,outputfile,pi0);

end