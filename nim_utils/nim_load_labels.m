function nim = nim_load_labels(nim)
% NIM_LOAD_LABELS Load parcellation labels into the nim structure
%
% Arguments:
%   nim - The nim structure
%
% Returns:
%   nim - Updated nim structure with atlas_labels field

arguments
    nim
end

labelfile = "";

if isfield(nim, 'parcellation_mask_file')
    [filepath, filename, ~] = fileparts(nim.parcellation_mask_file);

    % Priority 1: Check for XML labels file (from generate_parcellation_from_mni)
    xml_file = fullfile(filepath, 'atlas_labels.xml');
    if isfile(xml_file)
        labelfile = xml_file;
    else
        % Priority 2: Try the standard naming convention from nim_preprocessing
        % Replace _M with '' and add _atlas_labels.mat
        base_filename = strrep(filename, '_M', '');
        potential_mat_file = fullfile(filepath, [base_filename '_atlas_labels.mat']);

        if isfile(potential_mat_file)
            labelfile = potential_mat_file;
        else
            % Fallback: look for any atlas labels file in the directory
            potential_mat_file = fullfile(filepath, 'sample_atlas_labels.mat');
            if isfile(potential_mat_file)
                labelfile = potential_mat_file;
            end
        end

        % If still not found, search for any *_atlas_labels.mat file
        if ~isfile(labelfile)
            files = dir(fullfile(filepath, '*_atlas_labels.mat'));
            if ~isempty(files)
                labelfile = fullfile(filepath, files(1).name);
                fprintf('Found atlas labels file: %s\n', labelfile);
            end
        end
    end
end

% Check if file exists
if ~isfile(labelfile)
    warning('Label file not found for mask: %s. Atlas labels will not be available.', nim.parcellation_mask_file);
    return;
end

fprintf('Loading parcellation labels from %s\n', labelfile);

% Initialize the label map if it doesn't exist
if ~isfield(nim, 'atlas_labels') || ~isfield(nim.atlas_labels, 'map')
    nim.atlas_labels.map = containers.Map('KeyType', 'double', 'ValueType', 'char');
end
nim.atlas_labels.file = labelfile;

% Determine file type and load accordingly
[~, ~, ext] = fileparts(labelfile);

if strcmp(ext, '.xml')
    % Load XML file (from FSL atlases)
    try
        xml_content = fileread(labelfile);
        dom = xmlread(labelfile);
        labels = dom.getElementsByTagName('label');

        for i = 0:labels.getLength()-1
            label = labels.item(i);
            index_attr = label.getAttribute('index');
            if ~isempty(index_attr)
                label_index = str2double(char(index_attr));
                label_name = char(label.getTextContent());
                label_name = strtrim(label_name);
                nim.atlas_labels.map(label_index) = label_name;
            end
        end
        fprintf('Loaded atlas labels from XML file (%d labels)\n', nim.atlas_labels.map.Count);
    catch ex
        warning('Error loading XML file: %s', ex.message);
    end
else
    % Load MAT file
    try
        loaded_data = load(labelfile);
        if isfield(loaded_data, 'atlas_labels')
            % Copy fields from the loaded atlas_labels
            field_names = fieldnames(loaded_data.atlas_labels);
            for i = 1:length(field_names)
                nim.atlas_labels.(field_names{i}) = loaded_data.atlas_labels.(field_names{i});
            end
            fprintf('Loaded atlas labels from MAT file (%d labels)\n', nim.atlas_labels.map.Count);
        else
            warning('MAT file does not contain atlas_labels structure');
        end
    catch ex
        warning(ex.identifier, 'Error loading MAT file: %s', ex.message);
    end
end

fprintf('Loaded %d parcellation labels\n', nim.atlas_labels.map.Count);
end
