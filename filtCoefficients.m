% compute filter coefficients for rank order filter
% code include several parameters to change filter behavior

% pcsm_expType:  1 for exponental decrease of contributions 
%                2 for Gaussian shaped attenuation
% pcsm_sigma: sigma parameter (variance) for gaussian function,
%             and controls the slope for exponential attenuation
% divs_expo:  for higher values (>1), amplify the coefficients 
% divsSum: sum of coefficients
% rmside: damp cells in the opposite side if PfPc!=50, higher-more damping
%         choose between [>=0, <1] , = 0 means no effect
% PfWindow: filter length ( (window_width)*(window_height) )

divsSum = 2; % sum of coefficients
pcsm_expType = 2; % for gaussian shape = 2, for exponential = 1
PfPc = 50; % center of the filter (percent length of the window size) => 50% median
divs_expo = 1; % amplifying exponential
pcsm_sigma = 0.001; % sigma parameter for gaussian attenuation (1 is exponential)
PfWindow = 9; % number of filter coefficients = number of cells in the window = (x*z)
rmSide = 0;

if length(pcsm_sigma) == 1
    pcsm_sigma = [pcsm_sigma, 1];
end
if pcsm_expType==1
    pcsm_sigma = 1;
end

pcoeff = zeros(PfWindow) ;

% pc is the center of the filter coefficients 50% in standard
pc =  round(   ( (PfWindow-1) * PfPc / 100 )+1);
divs = zeros(1,PfWindow);
for v1 = 1:PfWindow
    if pcsm_expType == 1 
        divs(v1) =  1 / floor( 1./pcsm_sigma(1) * ...
            (abs(v1-pc) + 1) ) ^ divs_expo;
    elseif pcsm_expType == 2 
        divs(v1) = exp ( -( (v1-pc) / pcsm_sigma(2)  ) ...
            ^ 2 / (2 * pcsm_sigma(1) ^ 2) ) ^ divs_expo;
    end
end


if rmSide > 0  && PfPc~=50
  if PfPc < 50
          divs( ( floor(pc + 1) ) : end ) = ...
              divs( ( floor(pc + 1) ) : end ) * (1 - rmSide);
  elseif PfPc > 50
          divs( 1 : ( floor(pc - 1) ) ) = ...
              divs( 1 : ( floor(pc - 1) ) ) * (1 - rmSide);
  end
end

divs = (divsSum/sum(divs)) * divs;
disp( num2str(divs,'%2.6f, ') )

xtickLabels = split ( ...
              num2str( ...
              ( ( (1:PfWindow)-1)/(PfWindow-1) ) *100 ,' %4.1f' ) ...
              );
x = (-(PfWindow-1)/2): ((PfWindow-1)/2);
h = plot(x, divs,'-s','lineWidth', 0.8, 'markerSize',3, ...
    'DisplayName',['\sigma=' num2str(pcsm_sigma(1) )]); hold on
ylabel('coefficient')
xlabel('percentile (%)')
xticks(x);
set(gca,'xtickLabels',xtickLabels);
ylim([0 divsSum])
grid on
legend