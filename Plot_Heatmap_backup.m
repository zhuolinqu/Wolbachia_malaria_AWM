clearvars; close all; clc

Baseline_params_malaria
P = Baseline_params_gambiae(P);
tic
vw_vec = .5:.02:1;
alpha_vec = 0:.05:0.3; %0:.1:1;
ci_vec = 0.5:.02:1;%:.1:1; %0:.1:1;

outputM = NaN(length(ci_vec),length(vw_vec),length(alpha_vec));
outputW = NaN(length(ci_vec),length(vw_vec),length(alpha_vec));
outputMw1 = NaN(length(ci_vec),length(vw_vec),length(alpha_vec));

for j1 = 1:length(vw_vec)
    P.vw = vw_vec(j1);
    P.vu = 1 - P.vw;
    for j2 = 1:length(ci_vec)
        P.ci = ci_vec(j2);
        for j3 = 1:length(alpha_vec)
            P.alpha = alpha_vec(j3);

            phiW_list = linspace(1,50,200);
            %phiW_list = P.phiW;

            R0_M = NaN(length(phiW_list),3);
            R0_W = NaN(length(phiW_list),1);
            for k = 1:length(phiW_list)
                P.phiW = phiW_list(k);

                G0w = P.vw*P.bf*P.phiW/P.mufw;
                G0u = P.bf*P.phiU/P.mufu;
                R0w = P.vw*P.mufu*P.phiW/(P.mufw*P.phiU);

                a = P.vu/P.vw;
                b = (1-P.ci)/R0w+(P.vu/P.vw-1);
                c = (1-R0w)/R0w;

                rp = (-b + sqrt(b^2-4*a*c))/(2*a);
                rn = (-b - sqrt(b^2-4*a*c))/(2*a);

                for j = 1:3
                    if j==1
                        r = rp;
                        if r >=0
                            Sustar = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
                            Swstar = r*P.mufu/P.mufw*Sustar;
                        else
                            Sustar = NaN;
                            Swstar = NaN;
                        end
                    elseif j ==2
                        r = rn;
                        if r >= 0
                            Sustar = P.Kf/(1+r*P.mufu/P.mufw)*(1-1/G0w);
                            Swstar = r*P.mufu/P.mufw*Sustar;
                        else
                            Sustar = NaN;
                            Swstar = NaN;
                        end
                    elseif j==3
                        Sustar = P.Kf*(1-1/G0u);
                        Swstar = 0;
                    end


                    R0MH = P.bm*P.betaM*...
                        (P.alpha*Swstar/(Sustar+Swstar)*P.sigma/(P.mufw+P.sigma)*1/P.mufw+...
                        Sustar/(Sustar+Swstar)*P.sigma/(P.mufu+P.sigma)*1/P.mufu);

                    R0HM = P.bh*P.h*(P.h+P.muH)*...
                        (P.rho*P.betaD*1/(P.rD+P.muD+P.muH)+...
                        (1-P.rho)*P.betaA*1/(P.rA+P.muH)+...
                        P.rho*(1-P.phi)*P.rD/(P.rD+P.muD+P.muH)*P.betaA*1/(P.rA+P.muH));

                    R0_M(k,j) = sqrt(R0MH*R0HM);
                    R0_W(k) = R0w;
                end
            end
            if max(R0_W) <1
                disp('Issue with R0_W')
                keyboard
            end
            

            
            id = find(imag(R0_M(:,1))~=0,1,'last');
            
            if isempty(id)
                outM = NaN;
                outW = NaN;
            else
                if id < length(phiW_list)
                    outM = real(R0_M(id));
                    outW = real(R0_W(id));
                else
                    outM = NaN;
                    outW = NaN;
                end
            end

            %%%% WARNING - this is not correct %%%%
            clear id2
            if outM > 1
                R0_M2 = R0_M(id:end,1);
                R0_W2 = R0_W(id:end,1);
                id2 = find(R0_M2<1,1,'first');
                if isempty(id2)
                    outMw1 = NaN;
                else
                    if id < length(phiW_list)
                        outMw1 = real(R0_W2(id2));
                    else
                        outMw1 = NaN;
                    end
                end
            else
                outMw1 = NaN;
            end


            %param(j2,j1,j3) = [P.vw P.ci P.alpha]
            outputM(j2,j1,j3) = outM;
            outputW(j2,j1,j3) = outW;
            outputMw1(j2,j1,j3) = outMw1;
            

            
            %{
figure(2)
plot(R0_W,R0_M1,'linewidth',2)
hold on
set(gca,'fontsize',18)
xlabel('R0_W')
ylabel('R_0^M')
legend('EE+','EE-','DFE')
grid on
            %}
            %disp(R0_M)
            R0_M1 = R0_M;
            R0_M1(1:id-1,1:2)= NaN;

            %{
            figure(1)
            plot(R0_W,R0_M1(:,1),'b','linewidth',2)
            hold on
            plot(R0_W,R0_M1(:,2),'b--','linewidth',2)
            plot(R0_W,R0_M1(:,3),'r','linewidth',2)
            set(gca,'fontsize',18)
            xlabel('R_0^W')
            ylabel('R_0^M')
            legend('EE+','EE-','DFE')
            grid on
            pause(.1)
            %}
            


        end

            id3 = find(outputM(j2,j1,:))
            if outputM(j2,j1,j3) > 0.9 & outputM(j2,j1,j3) < 1.1
                outsurf(j2,j1) = P.alpha;
            else 
                outsurf(j2,j1) = NaN;
            end

    end
end
outputMw1(outputMw1>1) = NaN;
toc
%%
%{
figure(5)
plot(vw_vec,outputM(:,:,1),'linewidth',2)
set(gca,'fontsize',14)
xlabel('v_w')
ylabel('R_0^M')
%xlim([0 1])
%}
%%

%aid = 1;

if length(alpha_vec)<=4
    aidr = 2;
    aidc = 2;
elseif length(alpha_vec)<=9
    aidr = 3;
    aidc = 3;
else
    disp('Error in plotting - too many alpha values')
end
%%
for aid = 1:length(alpha_vec)

figure(20)
subplot(aidr,aidc,aid)
h = imagesc(vw_vec,ci_vec,squeeze(outputM(:,:,aid)));
set(h,'AlphaData',~isnan(squeeze(outputM(:,:,aid))))
set(gca,'Ydir','normal','fontsize',14)
xlabel('v_w')
ylabel('c_i')
colorbar
sgtitle('R_0^M, EE+/EE- appearance')

hold on
[C,h] = contour(vw_vec,ci_vec,squeeze(outputM(:,:,aid)),[1 1],'k','linewidth',2);
clabel(C,h)
hold off
caxis([.5 2])
title(sprintf('alpha = %0.1f',alpha_vec(aid)))

%
figure(30)
subplot(aidr,aidc,aid)
h = imagesc(vw_vec,ci_vec,squeeze(outputW(:,:,aid)));
set(h,'AlphaData',~isnan(squeeze(outputW(:,:,aid))))
set(gca,'Ydir','normal','fontsize',14)
xlabel('v_w')
ylabel('c_i')
colorbar
sgtitle('R_0^W, EE+/EE- appearance')
hold on
[C,h] = contour(vw_vec,ci_vec,squeeze(outputW(:,:,aid)),'k','linewidth',2);
clabel(C,h)
hold off
caxis([0 1])
title(sprintf('alpha = %0.1f',alpha_vec(aid)))
%
figure(40)
subplot(aidr,aidc,aid)
h = imagesc(vw_vec,ci_vec,squeeze(outputMw1(:,:,aid)));
set(h,'AlphaData',~isnan(squeeze(outputMw1(:,:,aid))))
hold on
[C,h] = contour(vw_vec,ci_vec,squeeze(outputMw1(:,:,aid)),.5:.1:.9,'k','linewidth',2);
clabel(C,h)
hold off

set(gca,'Ydir','normal','fontsize',14)
xlabel('v_w')
ylabel('c_i')
colorbar
sgtitle('R_0^W, cross R_0^M=1')
caxis([.5 1])
title(sprintf('alpha = %0.2f',alpha_vec(aid)))
end
toc