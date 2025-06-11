function [] = create_dirs(folderList, subjid)
    % CREATE_DIRS - Creates new subject or project directories if needed.
    % This helps ensure that output folders exist before saving data.

    if exist('subjid', 'var')
        for i = 1:length(folderList)
            if not(isfolder([folderList{i}, subjid]))
                mkdir([folderList{i}, subjid])
            end
        end
    else
    
        for i = 1:length(folderList)
            if not(isfolder([folderList{i}]))
                mkdir([folderList{i}])
            end
        end
    end

end