% Cell-wide Settings
% The cell-wide parameters are grouped into a single structure enb. A
% number of the functions used in this example require a subset of the
% parameters specified below. In this example we use the configuration
% according to the RMC R.14 FDD specified in TS 36.101 Annex A.3.4 which
% uses 50 RB, 4 port, 'SpatialMux' transmission scheme, '16QAM' symbol
% modulation, 2 codewords and a code rate of 1/2.
enb.NDLRB = 50;                 % Number of resource blocks
enb.CellRefP = 4;               % Cell-specific reference signal ports
enb.NCellID = 0;                % Cell ID
enb.CyclicPrefix = 'Normal';    % Normal cyclic prefix
enb.CFI = 2;                    % Length of control region
enb.DuplexMode = 'FDD';         % FDD duplex mode
enb.TDDConfig = 1;              % Uplink/Downlink configuration (TDD only)
enb.SSC = 4;                    % Special subframe configuration (TDD only)
enb.NSubframe = 0;              % Subframe number

% Transport/Physical channel settings for ease of use the DL-SCH and PDSCH
% channel specific settings are specified in a parameter structure pdsch.
% For the R.14 FDD RMC, there are two codewords, so the modulation scheme
% is specified as a cell array containing the modulation schemes of both
% codewords. If configuring for one codeword, the modulation scheme can be
% a character vector or a cell array with character vectors.
% It is also important to configure the TrBlkSizes parameter to have the
% correct number of elements as the intended number of codewords. The
% number of soft bits for the rate matching stage is decided by the UE
% category as specified in TS 36.306 Table 4.1-1. In this example, the
% transport block size is looked up from tables in TS 36.101 Annex A.3.4.
% This can also be done by using the lteRMCDL function for R.14 RMC.

% DL-SCH Settings
TrBlkSizes = [11448; 11448];    % 2 elements for 2 codeword transmission
pdsch.RV = [0 0];               % RV for the 2 codewords
pdsch.NSoftbits = 1237248;      % No of soft channel bits for UE category 2
% PDSCH Settings
pdsch.TxScheme = 'SpatialMux';  % Transmission scheme used
pdsch.Modulation = {'16QAM','16QAM'}; % Symbol modulation for 2 codewords
pdsch.NLayers = 2;              % Two spatial transmission layers
pdsch.NTxAnts = 2;              % Number of transmit antennas
pdsch.RNTI = 1;                 % The RNTI value
pdsch.PRBSet = (0:enb.NDLRB-1)';% The PRBs for full allocation
pdsch.PMISet = 0;               % Precoding matrix index
pdsch.W = 1;                    % No UE-specific beamforming
% Only required for 'Port5', 'Port7-8', 'Port8' and 'Port7-14' schemes
if any(strcmpi(pdsch.TxScheme,{'Port5','Port7-8','Port8', 'Port7-14'}))
    pdsch.W = transpose(lteCSICodebook(pdsch.NLayers,pdsch.NTxAnts,[0 0]));
end



% Random number initialization for creating random transport block(s)
rng('default');

% Convert the modulation scheme char array or cell array to string array
% for uniform processing
 pdsch.Modulation = string(pdsch.Modulation);

% Get the number of codewords from the number of transport blocks
nCodewords = numel(TrBlkSizes);

% Generate the transport block(s)
trBlk = cell(1,nCodewords); % Initialize the codeword(s)
for n=1:nCodewords
    trBlk{n} = randi([0 1],TrBlkSizes(n),1);
end
% Get the physical channel bit capacity required for rate matching from
% ltePDSCHIndices info output
[~,pdschInfo] = ltePDSCHIndices(enb,pdsch,pdsch.PRBSet);

% Define a structure array with parameters for lteRateMatchTurbo
chs = pdsch;
chs(nCodewords) = pdsch; % For 2 codewords, the array has two elements
% Initialize the codeword(s)
cw = cell(1,nCodewords);
for n=1:nCodewords
    % CRC addition for the transport block
    crccoded = lteCRCEncode(trBlk{n},'24A');
    % Code block segmentation returns a cell array of code block segments
    % with filler bits and type-24B CRC appended as required
    blksegmented = lteCodeBlockSegment(crccoded);
    % Channel coding returns the turbo coded segments in a cell array
    chencoded = lteTurboEncode(blksegmented);

    % Bundle the parameters in structure chs for rate matching as the
    % function requires both cell-wide and channel specific parameters
    chs(n).Modulation = pdsch.Modulation{n};
    chs(n).DuplexMode = enb.DuplexMode;
    chs(n).TDDConfig = enb.TDDConfig;
    % Calculate number of layers for the codeword
    if n==1
        chs(n).NLayers = floor(pdsch.NLayers/nCodewords);
    else
        chs(n).NLayers = ceil(pdsch.NLayers/nCodewords);
    end
    % Rate matching returns a codeword after sub-block interleaving, bit
    % collection and bit selection and pruning defined for turbo encoded
    % data and merging the cell array of code block segments
    cw{n} = lteRateMatchTurbo(chencoded,pdschInfo.G(n),pdsch.RV(n),chs(n));
end



% Initialize the modulated symbols
modulated = cell(1,nCodewords);
for n=1:nCodewords
   % Generate the scrambling sequence
   scramseq = ltePDSCHPRBS(enb,pdsch.RNTI,n-1,length(cw{n}));
   % Scramble the codewords
   scrambled = xor(scramseq,cw{n});
   % Symbol modulate the scrambled codewords
   modulated{n} = lteSymbolModulate(scrambled,pdsch.Modulation{n});
end
% Layer mapping results in a (symbols per layer)-by-NLayers matrix
layermapped = lteLayerMap(pdsch,modulated);
% Precoding results in a (symbols per antenna)-by-NTxAnts matrix
precoded = lteDLPrecode(enb, pdsch, layermapped);
% Apply beamforming optionally (W should be 1 or identity if no beamforming)
pdschsymbols = precoded*pdsch.W;



% Deprecoding (pseudo-inverse based) returns (Number of symbols)-by-NLayers matrix
if (any(strcmpi(pdsch.TxScheme,{'Port5' 'Port7-8' 'Port8' 'Port7-14'})))
    rxdeprecoded=pdschsymbols*pinv(pdsch.W);
else
    rxdeprecoded = lteDLDeprecode(enb,pdsch,pdschsymbols);
end
% Layer demapping returns a cell array containing one or two codewords. The
% number of codewords is deduced from the number of modulation scheme
% character vectors
layerdemapped = lteLayerDemap(pdsch,rxdeprecoded);

% Initialize the recovered codewords
cws = cell(1,nCodewords);
for n=1:nCodewords
    % Soft demodulation of received symbols
    demodulated = lteSymbolDemodulate(layerdemapped{n},pdsch.Modulation{n},'Soft');
    % Scrambling sequence generation for descrambling
    scramseq = ltePDSCHPRBS(enb,pdsch.RNTI,n-1,length(demodulated),'signed');
    % Descrambling of received bits
    cws{n} = demodulated.*scramseq;
end



% Initialize the received transport block and CRC
rxTrBlk = cell(1,nCodewords);
crcError = zeros(1,nCodewords);
for n=1:nCodewords
    % Rate recovery stage also allows combining with soft information for
    % the HARQ process, using the input cbsbuffers. For the first
    % transmission of the transport block, the soft buffers are initialized
    % as empty. For retransmissions, the parameter cbsbuffers should be the
    % soft information from the previous transmission
    cbsbuffers = [];  % Initial transmission of the HARQ process
    % Rate recovery returns a cell array of turbo encoded code blocks
    raterecovered = lteRateRecoverTurbo(cws{n},TrBlkSizes,pdsch.RV(n),chs(n),cbsbuffers);
    NTurboDecIts = 5; % Number of turbo decoding iteration cycles
    % Turbo decoding returns a cell array of decoded code blocks
    turbodecoded = lteTurboDecode(raterecovered,NTurboDecIts);
    % Code block desegmentation concatenates the input code block segments
    % into a single output data block, after removing any filler and
    % type-24B CRC bits that may be present
    [blkdesegmented,segErr] = lteCodeBlockDesegment(turbodecoded,(TrBlkSizes+24));
    % CRC decoding returns the transport block after checking for CRC error
    [rxTrBlk{n},crcError(n)] = lteCRCDecode(blkdesegmented,'24A');
end



