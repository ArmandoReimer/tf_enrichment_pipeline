% prefix_cell = {'2019-03-26-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW',...
%                '2019-03-20-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW_01',...
%                '2019-03-20-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW_02',...
%                };
% prefix_cell = {'2019-03-21-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW_03',...
%                '2019-03-14-Dl_Venus_hbP2P_MCPmCherry_Leica_Zoom2_7uW14uW_01',...
%                '2019-04-25-Dl_Venus_hbP2P_MCPmCherry_Leica_Zoom2_7uW14uW_05'};
prefix_cell = {'2019-03-04-Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW_14uW_02'};


for i = 1:numel(prefix_cell)
    disp('fitting...')
    tic
    Prefix = prefix_cell{i}
    fit3DGaussiansToAllSpots(Prefix);
    toc
    disp('done.')
end

