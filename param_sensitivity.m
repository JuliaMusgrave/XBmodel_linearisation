%

% This function performs a sensitivity analysis on the best fit to
% XBmodel_linear.m
% 
% It looks at the complex modulus at the optimal solution and on 10 % 
% either side of each parameter to look at the sensitivity and effect of 
% each of the parameters
%
% Inputs: 
%       - params: fitting parameters used in the model, at the optimum value
%       - data: 1D array of complex modulus the model was fit to
%       - freqs: 1D array of the sampled frequencies 
%       - supp: boolean input indicating whether Figure S3 from the supplementary
%       should also be plotted
%
% Author: Julia Musgrave
% Date: May 2022


function [] = param_sensitivity(params,data,freqs,supp)

no_ps=length(params);

% looking at 100%, 110%, and 90% of each parameter value
dp=[1 1.1 0.9 1]; 

ps={'{\itk}_1' '{\itk}_2' '{\itk''}_3' '{\itk}_{-1}' '{\itk}_{-2}' '{\itK}' '\phi_{\itx}' '\phi_{\itv}' '\phi_{\its}' '\phi_{\itl}'};
load MLcolours.mat
colours={'k',blue,red,'k'};

%plotting prep
if supp
w=0.4;
h=0.16;
fs=logspace(-1,2,100);
figure('Units', 'normalized' ,'OuterPosition',  [0.3, 0.05, 0.4, 0.9])
end

% cycling through all the parameters
for i=1:no_ps
    
    x=params;

    if supp
    ip=i;
    if i>5
        ip=i-5;
    end
    yp=1.01-ip*0.19;
    end

    for j=1:4 
        %varying each by 10 % from optimal
        x(i)=params(i)*dp(j);
        RMSE(i,j)=XBmodel_linear(x,data,freqs);
        [~,Y]=XBmodel_linear(x);

        if supp
        subplot('Position',[0.075 yp w h])
        semilogx(fs,real(Y),'Color',colours{j},'LineWidth',1)
        hold on
                
        subplot('Position',[0.575 yp w h])
        semilogx(fs,imag(Y),'Color',colours{j},'LineWidth',1)
        hold on
        if j==3 && mod(i,5)==1 
             legend('Model fit','10% increase','10% decrease','AutoUpdate','off','Box','off','Location','northwest')
        end
        end
    end
        
        % formatting each axis
        if supp
        subplot('Position',[0.075 yp w h])
        ylabel('Elastic Mod (MPa)','FontSize',12)
        xlim([0.1 100])
        box off
        xticklabels({'0.1' '1' '10' '100'})
        annotation('textbox',[.08 yp+0.055 .1 .1],'String',ps{i},'FitBoxToText','on','FontWeight','bold');
        text(0.04,max(ylim),char(63+2*i),'FontSize',14,'FontWeight','bold')
        if mod(i,5)==0
        xlabel('Frequency (Hz)','FontSize',12)
        end
        hold off

        subplot('Position',[0.575 yp w h])
        ylabel('Viscous Mod (MPa)','FontSize',12)
        box off
        xlim([0.1 100])
        xticklabels({'0.1' '1' '10' '100'})
        text(0.04,max(ylim),char(64+2*i),'FontSize',14,'FontWeight','bold')
        if mod(i,5)==0
        xlabel('Frequency (Hz)','FontSize',12)
        end
        if i==5
        figure('Units', 'normalized' ,'OuterPosition',  [0.3, 0.05, 0.4, 0.9])
        end
        end


end

% box plot figure
figure('Name','Sensitivity analysis (Fig 8)','Units', 'normalized' ,'OuterPosition',  [0.3, 0.3, 0.22, 0.3])
rmse_perc=(RMSE-RMSE(1,1))/RMSE(1,1)*100;
x=categorical(ps);
x=reordercats(x,ps);
bar(x,mean(rmse_perc(:,2:3),2),'w','EdgeColor','black')
box off
ylabel('Change in RMSE (%)','FontSize',12)


end

