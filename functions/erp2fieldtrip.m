function ft_struct = erp2fieldtrip(ERP, ft_struct, event_id)
    % ERP2FIELDTRIP Converts ERPlab ERP data to FieldTrip timelockaverage structure.
    %
    %   Description: 
    %   Takes the average Event-Related Potentials (ERPs) 
    %   from an ERPlab structure and inputs them
    %   into a FieldTrip timelockaverage structure.
    %
    %   Inputs:
    %       ERP       - ERPlab structure containing ERP data, including fields
    %                   such as bindata, binerror, and ntrials.accepted.
    %       ft_struct - FieldTrip structure to be populated with ERP data.
    %       event_id  - String identifier for the event of interest.
    %
    %   Outputs:
    %       ft_struct - Updated FieldTrip structure containing the average ERP
    %                   data, variance, and degrees of freedom for the specified
    %                   event.
    %
    %   Notes:
    %       - The function handles ERP data with 31 or 32 channels specifically,
    %         and attempts to process other channel counts using a try-catch block.
    %       - Ensure that the ERP structure is correctly formatted with the
    %         necessary fields before calling this function.

    % get relevant data
    event_idx = find(strcmp(sprintf('.{%s}', event_id),{ERP.EVENTLIST.bdf.expression}));
    
    % if 31 or 32 channels, extracts the average and variance of the ERP data for the specified event
    if size(ERP.bindata,1) == 32
        ft_struct.avg = ERP.bindata(1:32,:,event_idx );
        ft_struct.var = ERP.binerror(1:32,:,event_idx ).^2*ERP.ntrials.accepted(event_idx); % get the variance from the ERP SEM
        ft_struct.dof = ones(size(ft_struct.avg ));
    elseif size(ERP.bindata,1) == 31
        ft_struct.avg = ERP.bindata(1:31,:,event_idx );
        ft_struct.var = ERP.binerror(1:31,:,event_idx ).^2*ERP.ntrials.accepted(event_idx); % get the variance from the ERP SEM
        ft_struct.dof = ones(size(ft_struct.avg ));
    else
        % extract data for all available channels
        try
            ft_struct.avg = ERP.bindata(1: size(ERP.bindata,1),:,event_idx );
            ft_struct.var = ERP.binerror(1: size(ERP.bindata,1),:,event_idx ).^2*ERP.ntrials.accepted(event_idx); % get the variance from the ERP SEM
            ft_struct.dof = ones(size(ft_struct.avg ));
        catch
        end
    end
end