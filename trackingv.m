function [trackResults, channel]= trackingv(fid, channel,trackRes,navSolutions,eph,activeChnList,svTimeTable, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record for I
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================
tic
kpt=1e-3;%kalman filter process update time
kmt=1e-3;%kalman filter measurement update time/s
pdi=1e-3;%integration time /s
tracklengthall=4000/(pdi*1000);% vector tracking total length /ms
%vector tracking section length /ms, must can be exactly divided by tracklengthall
tracklength=4000/(pdi*1000);
%vector tracking section #, must be integer
sectno=fix(tracklengthall/tracklength);

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, tracklength);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, tracklength);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, tracklength);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, tracklength);
trackResults.I_E            = zeros(1, tracklength);
trackResults.I_L            = zeros(1, tracklength);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, tracklength);
trackResults.Q_P            = zeros(1, tracklength);
trackResults.Q_L            = zeros(1, tracklength);

% Loop discriminators
trackResults.dllDiscr       = inf(1, tracklength);
trackResults.dllDiscrFilt   = inf(1, tracklength);
trackResults.pllDiscr       = inf(1, tracklength);
trackResults.pllDiscrFilt   = inf(1, tracklength);


%C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(tracklength/settings.CNo.VSMinterval));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(tracklength/settings.CNo.VSMinterval));

trackResults.CNo.PRMValue=0; %To avoid error message when
trackResults.CNo.PRMIndex=0; %tracking window is closed before completion.

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

% codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%% init
%vector tracking start time in dataset /ms
StartTime=15000;
%# of available channels
NumChan=length(activeChnList);%settings.numberOfChannels;
%one step transition matrix
F=diag([0,0,0,0,0,0,1,0]);
F(1,4)=kpt;
F(2,5)=kpt;
F(3,6)=kpt;
F(7,8)=kpt;
%process noise covariance matrix
Qw=diag([2e1,2e1,2e1,3e0,3e0,3e0,1e-7,1e-1]);

%measurement noise covariance matrix, set the init value empirically
R(1:NumChan,1:NumChan)=1500*eye(NumChan);%3000*eye(2*NumChan); % dimension=number of available channels
R(NumChan+1:2*NumChan,NumChan+1:2*NumChan)=9e2*eye(NumChan);%8e1*eye(NumChan);%

%kalman gain matrix
K=zeros(8,2*NumChan);
%estimate error matrix
P0=diag([1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1,1e0]);
%one step prediction error matrix
P=zeros(8,8);
%one step prediction of states
X_next=zeros(8,tracklength);
%clock drift record
recordddt0=zeros(1,tracklength);
%state estimates
X_est=X_next;
%states initial value
X0=[0,0,0,0,0,0,0,0]';
dt=X0(7);%clock bias
ddt=X0(8);%clock drift error
ddt0=-navSolutions.ddt(StartTime*settings.navSolRate/1000);%intial of clock drift
H=zeros(2*NumChan,8);%measurement matrix
Z=zeros(2*NumChan,tracklength);%mearurements
%initial values of position and velocity from SLL
pos0=[navSolutions.X(StartTime*settings.navSolRate/1000),navSolutions.Y(StartTime*settings.navSolRate/1000),navSolutions.Z(StartTime*settings.navSolRate/1000)];
vel0=1*[navSolutions.Vx(StartTime*settings.navSolRate/1000),navSolutions.Vy(StartTime*settings.navSolRate/1000),navSolutions.Vz(StartTime*settings.navSolRate/1000)];

% adaptive parameters for R
cnt=1;
%set the lengh of consequence data needed to get the adaptive R value
lastn=500;
mat2=zeros(2*NumChan,lastn);

dSv=zeros(1,NumChan);%delta sat position in 1ms
dPlos=zeros(1,NumChan);%delta receiver position in 1ms
Vs=zeros(1,NumChan);%relative velocity between receiver and sat
dVlos=zeros(1,NumChan);%delta receiver velocity in 1ms

carrError=zeros(1,NumChan);
carrErrorold=zeros(1,NumChan);
codeError=zeros(1,NumChan);
codeErrorold=zeros(1,NumChan);

carrFreq=zeros(1,NumChan);
codeFreq=zeros(1,NumChan);
initsample=zeros(1,NumChan);
remCodePhase=zeros(1,NumChan);
codePhaseStep=zeros(1,NumChan);


%record all R matrix
recordR=zeros(2*NumChan,tracklength);
recordtxTime=zeros(NumChan,tracklength);

%init for code phase/freq, carrier freq, transmittime
for channelNr = 1:NumChan%settings.numberOfChannels

    % Only process if PRN is non zero (acquisition was successful)
        trackResults(activeChnList(channelNr)).PRN     = trackRes(1,activeChnList(channelNr)).PRN;
        carrFreq(1,channelNr)=trackRes(1,activeChnList(channelNr)).carrFreq(StartTime);
        codeFreq(1,channelNr)=trackRes(1,activeChnList(channelNr)).codeFreq(StartTime);

        initsample(1,channelNr)=ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(StartTime));
        initsampleforcode(1,channelNr)=ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(StartTime-1));
        remCodePhase(1,channelNr)=(initsampleforcode(1,channelNr)-trackRes(1,activeChnList(channelNr)).absoluteSample(StartTime-1))...
            /settings.samplingFreq*codeFreq(1,channelNr);
        codePhaseStep(1,channelNr) = codeFreq(1,channelNr) / settings.samplingFreq;
        tTime=...
            findTransTime(initsample(channelNr),activeChnList,svTimeTable,trackRes);
        transmitTime(activeChnList(channelNr))=tTime(activeChnList(channelNr));
        remCarrPhase(1,channelNr)=0;
 
        %C/No computation
%         vsmCnt(channelNr)  = 0;
        if (settings.CNo.enableVSM==1)
            CNo='Calculating...';
        else
            CNo='Disabled';
        end
        % Get a vector with the C/A code sampled 1x/chip
        caCode0 = generateCAcode(trackRes(1,activeChnList(channelNr)).PRN);
        % Then make it possible to do early and late versions
        caCode(channelNr,:) =[caCode0(1023) caCode0 caCode0(1)];

end % for channelNr

        transmitTime0=transmitTime;
        blksize0=pdi*settings.samplingFreq*ones(1,NumChan);
%         blksize=pdi*settings.samplingFreq*ones(1,NumChan);
        blksize=ceil(pdi*settings.samplingFreq*ones(1,NumChan));
        samplepos=initsample;
        mininit=min(initsample);
        minpos=mininit;
        satClkCorr=0;
        
%% Start processing channels ==============================================
for sectcnt=1:sectno
        %=== Process the number of specified code periods =================
        for loopCnt =  1:tracklength
            %% Set up all the code phase tracking information -------------------------
            % Define index into early code vector
            transmitTime(activeChnList)=transmitTime(activeChnList)-...
                (-1*dt)/settings.c+blksize/settings.samplingFreq;
            

            [satPositionsall, satClkCorrall] = satpos([transmitTime(transmitTime>0),transmitTime0(transmitTime0>0)], ...
                    [trackRes(activeChnList).PRN,trackRes(activeChnList).PRN],eph);
            satPositions=satPositionsall(:,1:NumChan);
            satClkCorr=satClkCorrall(1:NumChan);
            satPositions0=satPositionsall(:,end/2+1:end);
            satClkCorr0=satClkCorrall(NumChan+1:2*NumChan);


            fseek(fid, ...
                dataAdaptCoeff*(0*settings.skipNumberOfSamples + ...
                minpos-1), ...
                'bof');
            [rawSignal0, samplesRead0] = fread(fid, ...
                    dataAdaptCoeff*(max(samplepos)-minpos+max(blksize)), settings.dataType);
            for m=1:NumChan
                %% Read next block of data ------------------------------------------------
                rawSignal= rawSignal0((samplepos(m)-minpos)*dataAdaptCoeff+1:(samplepos(m)-minpos+blksize(1,m))*dataAdaptCoeff)';

                if (dataAdaptCoeff==2)
%                     rawSignal1=rawSignal(1:2:end);%xin20120704
                    rawSignal1=rawSignal(1:2:end-1);%xin20120704
                    rawSignal2=rawSignal(2:2:end);
                    rawSignal = rawSignal1 + 1i .* rawSignal2;  %transpose vector
                end

                lx=satPositions(1,m)-pos0(1);
                ly=satPositions(2,m)-pos0(2);
                lz=satPositions(3,m)-pos0(3);
                norm_a=sqrt(lx*lx+ly*ly+lz*lz);
                a=[lx;ly;lz]/norm_a;
                H(m,:)=[-a(1),-a(2),-a(3),0,0,0,-1,0];
                H(m+NumChan,:)=[0,0,0,+a(1),+a(2),+a(3),0,1];
                dSv(m)=...
                    [satPositions(1,m)-satPositions0(1,m),...
                    satPositions(2,m)-satPositions0(2,m),...
                    satPositions(3,m)-satPositions0(3,m)]*a;

                Vs(m)=...
                    ([satPositions(4,m),...
                    satPositions(5,m),...
                    satPositions(6,m)]-vel0)*a;
                estPecef=pos0;
                estVecef=vel0;

                dPlos(m)=1*(1*X0(1:3)'+kmt*estVecef+estPecef-pos0)*a;
                dVlos(m)=1*X0(4:6)'*a;%X0(4:6)'*a;%
                
                codeFreq(1,m)=settings.codeFreqBasis*(1-1*(ddt0)/settings.c-Vs(m)/settings.c);
                codePhaseStep(1,m) = codeFreq(1,m) / settings.samplingFreq;
                remCodePhase(1,m) =remCodePhase(1,m) + (1*dt)/settings.c*codeFreq(1,m)-...
                    (dSv(m)-dPlos(m))/settings.c*codeFreq(1,m)+(blksize(1,m)-blksize0(1,m)).*codePhaseStep(1,m);


%% Local Code(E,P,L) Generation
                % Local Code generation - ready for action, 1 over 2
                tcode       = (remCodePhase(1,m)-earlyLateSpc) : ...
                                codePhaseStep(1,m) : ...
                                ((blksize(1,m)-1)*codePhaseStep(1,m)+remCodePhase(1,m)-earlyLateSpc);
                tcode2      = ceil(tcode) + 1;
                earlyCode   = caCode(m,tcode2);

                tcode       = (remCodePhase(1,m)+earlyLateSpc) : ...
                                codePhaseStep(1,m) : ...
                                ((blksize(1,m)-1)*codePhaseStep(1,m)+remCodePhase(1,m)+earlyLateSpc);
                tcode2      = ceil(tcode) + 1;
                lateCode    = caCode(m,tcode2);

                tcode       = remCodePhase(1,m) : ...
                                codePhaseStep(1,m) : ...
                                ((blksize(1,m)-1)*codePhaseStep(1,m)+remCodePhase(1,m));
                tcode2      = ceil(tcode) + 1;
                promptCode  = caCode(m,tcode2);
                %% Generate the carrier frequency to mix the signal to baseband -----------
                
                carrFreq(1,m)=+1*settings.IF+ddt0/settings.c*1575.42e6+Vs(m)/settings.c*1575.42e6;
                time    = (0:blksize(1,m)) ./ settings.samplingFreq;

                % Get the argument to sin/cos functions
                trigarg = ((carrFreq(1,m) * 2.0 * pi) .* time) + remCarrPhase(1,m);
                remCarrPhase(1,m) = rem(trigarg(blksize(1,m)+1), (2 * pi));

                % Finally compute the signal to mix the collected data to bandband
                carrsig = exp(1i .* trigarg(1:blksize(1,m)));

                %% Generate the six standard accumulated values ---------------------------
                % First mix to baseband
                qBasebandSignal = real(carrsig .* rawSignal);
                iBasebandSignal = imag(carrsig .* rawSignal);

                % Now get early, late, and prompt values for each
                I_E = sum(earlyCode  .* iBasebandSignal);
                Q_E = sum(earlyCode  .* qBasebandSignal);
                I_P = sum(promptCode .* iBasebandSignal);
                Q_P = sum(promptCode .* qBasebandSignal);
                I_L = sum(lateCode   .* iBasebandSignal);
                Q_L = sum(lateCode   .* qBasebandSignal);
                

                %% Find PLL error and update carrier NCO ----------------------------------

                % Implement carrier loop discriminator

                if (loopCnt==1 && sectcnt==1)
                    IP1(1,m)=I_P;
                    QP1(1,m)=Q_P;
                    carrErrorold(1,m)=carrError(1,m);
                    carrError(1,m)=0;
                else
                    dot=IP1(1,m)*I_P+QP1(1,m)*Q_P;
                    cross=IP1(1,m)*Q_P-I_P*QP1(1,m);
                    carrErrorold(1,m)=carrError(1,m);

                    carrError(1,m) = cross*sign(dot)/(2*pi*(I_P*I_P+Q_P*Q_P));

                    IP1(1,m)=I_P;
                    QP1(1,m)=Q_P;
                end
 
                
                Z(m+NumChan,loopCnt)=(carrErrorold(1,m)+(carrError(1,m)-carrErrorold(1,m))...
                    /blksize(1,m)*(blksize(1,m)-mod((samplepos(1,m)-minpos),blksize(1,m))))/pdi/1575.42e6*settings.c;

                if abs(Z(m+NumChan,loopCnt))>100
                    Z(m+NumChan,loopCnt)=0;%sign(Z(m+NumChan,loopCnt))*200;
                end


                trackResults(activeChnList(m)).carrFreq(loopCnt) = carrFreq(1,m)+(-dVlos(m))/settings.c*1575.42e6;

                %% Find DLL error and update code NCO -------------------------------------

                codeErrorold(1,m)=codeError(1,m);
                codeError(1,m) = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                    2/(sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

                Z(m,loopCnt)=(codeErrorold(1,m)+(codeError(1,m)-codeErrorold(1,m))...
                    /blksize(1,m)*(blksize(1,m)-mod((samplepos(1,m)-minpos),blksize(1,m))))/codeFreq(1,m)*settings.c;

                % Modify code freq based on NCO command
                trackResults(activeChnList(m)).codeFreq(loopCnt) = codeFreq(1,m);
                trackResults(activeChnList(m)).remCodePhase(loopCnt)=remCodePhase(1,m);
                %% Record various measures to show in postprocessing ----------------------
                % Record sample number (based on 8bit samples)
%                 trackResults(channelNr).absoluteSample(loopCnt) = (ftell(fid))/dataAdaptCoeff- remCodePhase/codePhaseStep;

                trackResults(activeChnList(m)).dllDiscr(loopCnt)       = codeError(1,m);
                trackResults(activeChnList(m)).I_E(loopCnt) = I_E;
                trackResults(activeChnList(m)).I_P(loopCnt) = I_P;
                trackResults(activeChnList(m)).I_L(loopCnt) = I_L;
                trackResults(activeChnList(m)).Q_E(loopCnt) = Q_E;
                trackResults(activeChnList(m)).Q_P(loopCnt) = Q_P;
                trackResults(activeChnList(m)).Q_L(loopCnt) = Q_L;

                if (settings.CNo.enableVSM==1)
                    if (rem(loopCnt,settings.CNo.VSMinterval)==0)
                        CNoValue=CNoVSM(trackResults(activeChnList(m)).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                            trackResults(activeChnList(m)).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
                        trackResults(activeChnList(m)).CNo.VSMValue(loopCnt/settings.CNo.VSMinterval)=CNoValue;
                        trackResults(activeChnList(m)).CNo.VSMIndex(loopCnt/settings.CNo.VSMinterval)=loopCnt;

                    end
                end

                trackResults(activeChnList(m)).status  = channel(channelNr).status;
                trackResults(activeChnList(m)).blksize(loopCnt)  = blksize(1,m);
            end

            %Kalman filter
            P=F*P0*F'+Qw;
            K=P*H'/(H*P*H'+R);
            P0=(eye(8)-K*H)*P;
            
            X_next(:,loopCnt)=F*X0;
            beta=(Z(:,loopCnt)-H*X_next(:,loopCnt));
            alpha=K*beta;
           
            X_est(:,loopCnt)=X_next(:,loopCnt)+alpha;
            X0=X_est(:,loopCnt);
            res2=beta;
            mat2(:,cnt)=res2;
            cnt=cnt+1;
            %if collected 'lastn' measurements, start over again
            if cnt==lastn
                cnt=1;
            end
            if (loopCnt>lastn ||sectcnt>1)

                recordR(:,loopCnt)=var(mat2,0,2);
                tmpR=diag(recordR(:,loopCnt));
                %update R adaptively
                R(1:NumChan,1:NumChan)=tmpR(1:NumChan,1:NumChan);
                R(NumChan+1:2*NumChan,NumChan+1:2*NumChan)=tmpR(NumChan+1:2*NumChan,NumChan+1:2*NumChan);

            end
            vel0=(-1)*X0(4:6)'+estVecef;
            pos0=1*X0(1:3)'+estPecef;
            
            Result.pos(loopCnt,:)=pos0';
            pos0=pos0+kmt*vel0;
            Result.vel(loopCnt,:)=vel0';

            dt=X0(7);
            ddt=X0(8);
            ddt0=ddt0+ddt;
            recordddt0(loopCnt)=ddt0;
            transmitTime0=transmitTime;
            recordtxTime(:,loopCnt)=transmitTime(activeChnList)';
            
            blksize = ceil((settings.codeLength-remCodePhase) ./ codePhaseStep);
            
            samplepos=samplepos+blksize;
            minpos=min(samplepos);
            
        end % for loopCnt

        disp('   Saving VLL results to file "vllresult.mat"')
%         fname=['vllresult',num2str(StartTime/1000+0*tracklength/1000*(sectcnt-1)),...
%             '~',num2str(StartTime/1000+tracklength/1000)];
        save('vllresult', ...
            'Result', 'settings', 'StartTime','trackResults', 'tracklength','activeChnList',...
            'recordR','recordddt0','P0','samplepos','blksize','recordtxTime');
%%
% if draw figures, uncomment them
        for i=1:length(activeChnList)%frequency tracking results figures
            figure
            plot(trackRes(1,activeChnList(i)).carrFreq(StartTime:StartTime+tracklength-1))
            grid on
            hold on
            plot(trackResults(1,activeChnList(i)).carrFreq,'r')
            hold off
            title(['Channel',num2str(activeChnList(i)),'PRN ',num2str(trackResults(activeChnList(i)).PRN)])
        end
        figure %X,Y,Z position results figures
        plot(navSolutions.X(StartTime/1000*settings.navSolRate:(StartTime+tracklength)/1000*settings.navSolRate-1)-navSolutions.X(StartTime/1000*settings.navSolRate))
        hold on
        plot(Result.pos(1:100:end,1)'-Result.pos(1,1),'r')
        hold off
        grid on
        figure
        plot(navSolutions.Y(StartTime/1000*settings.navSolRate:(StartTime+tracklength)/1000*settings.navSolRate-1)-navSolutions.Y(StartTime/1000*settings.navSolRate))
        hold on
        plot(Result.pos(1:100:end,2)'-Result.pos(1,2),'r')
        hold off
        grid on
        figure
        plot(navSolutions.Z(StartTime/1000*settings.navSolRate:(StartTime+tracklength)/1000*settings.navSolRate-1)-navSolutions.Z(StartTime/1000*settings.navSolRate))
        hold on
        plot(Result.pos(1:100:end,3)'-Result.pos(1,3),'r')
        hold off
        grid on

        figure %X,Y,Z velocity results figures
        plot(navSolutions.Vx(StartTime/1000*settings.navSolRate:(StartTime+tracklength)/1000*settings.navSolRate-1))
        hold on
        plot(Result.vel(1:100:end,1)','r')
        hold off
        grid on
        figure
        plot(navSolutions.Vy(StartTime/1000*settings.navSolRate:(StartTime+tracklength)/1000*settings.navSolRate-1))
        hold on
        plot(Result.vel(1:100:end,2)','r')
        hold off
        grid on
        figure
        plot(navSolutions.Vz(StartTime/1000*settings.navSolRate:(StartTime+tracklength)/1000*settings.navSolRate-1))
        hold on
        plot(Result.vel(1:100:end,3)','r')
        hold off
        grid on
        StartTime=StartTime+tracklength;
end

toc
disp('Vector tracking over')
