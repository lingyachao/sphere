function [DATA_DIR,type] = find_full_id(DATA_ROOT_DIR, id)

    folder_list = dir(DATA_ROOT_DIR);
    for i = 1 : length(folder_list)
        folder_name = folder_list(i).name;

        if ~isempty(strfind(folder_name, id))
            DATA_DIR = [DATA_ROOT_DIR folder_name '/'];

            if ~isempty(strfind(folder_name, 'brain'))
                type = 'brain';
            else
                type = 'sphere';
            end
        end
    end

end

