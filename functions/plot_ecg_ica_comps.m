%% Plot ECG-ICA Comps
% works in conjunction with ecg_ica_corr()
% takes rejV and cV as input + a path where the pictures should be saved
% todo: fix time scale fro sec to ms

function plot_ecg_ica_comps(EEG, rejV, cV, img_save_path, subjid)

    %Plot ICA Comps
    FigH = figure('Position', [475,142,1098,739]);
    set(FigH,'Visible','off')
    
    if mod(sum(rejV),2) == 0 
        cmp_to_plot = sum(rejV);
    else
        cmp_to_plot = (sum(rejV)+1);
    end
    
       
    %plot components
    cmpIdx = find(rejV);
    cmpPosi = [1:2:sum(rejV)+15]; %Positions of the Topoplots
    timePosi = [2:2:sum(rejV)+16]; %Positions of the ERP 
    
    for nplt = 1:sum(rejV)
        
        %Plot Topoplot
        subplot(cmp_to_plot,2,cmpPosi(nplt));
        topoplot(EEG.icawinv(:,cmpIdx(nplt)), EEG.chanlocs);
        title(['Component ' num2str(cmpIdx(nplt)) ' Corr ' num2str(cV(cmpIdx(nplt)))]) ;
        
        %Plot ERP & single trials
        subplot(cmp_to_plot,2,timePosi(nplt));
        icaacttmp  = eeg_getdatact(EEG, 'component', cmpIdx(nplt));
        offset     = nan_mean(icaacttmp(:));
        era        = nan_mean(squeeze(icaacttmp)')-offset;
        era_limits = get_era_limits(era);
        erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
            '', 3, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','erp_vltg_ticks',era_limits);       
    end
    
    saveas(gcf, [img_save_path, subjid, '_ECG_Corr_comp.png'],'png');
end