function labels = nim_load_atlas_labels(atlas_type)
% nim_load_atlas_labels: Load atlas labels from the corresponding XML file
%
% Arguments:
%   atlas_type - Type of atlas ('HarvardOxford' or 'JHU')
%
% Returns:
%   labels - Structure containing index-to-name mappings for the atlas

% Get FSL directory
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
  error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define XML file paths based on atlas type
if strcmpi(atlas_type, 'JHU')
  % For JHU atlas, check multiple possible locations
  potential_files = {
    fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-labels-1mm.xml'),
    fullfile(fsl_path, 'data/atlases/JHU-labels.xml'),
    '/Users/12salty/Documents/research-chun/fsl/data/atlases/JHU-labels.xml' % Direct path to the file you shared
    };
  
  xml_file = '';
  for i = 1:length(potential_files)
    if exist(potential_files{i}, 'file')
      xml_file = potential_files{i};
      fprintf('Found JHU atlas XML file: %s\n', xml_file);
      break;
    end
  end
elseif strcmpi(atlas_type, 'JHU-tract')
  % For JHU-tract atlas, check multiple possible locations
  potential_files = {
    fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.xml'),
    fullfile(fsl_path, 'data/atlases/JHU-labels.xml'),
    '/Users/12salty/Documents/research-chun/fsl/data/atlases/JHU-tracts.xml' % Direct path to the file you shared
    };
  
  xml_file = '';
  for i = 1:length(potential_files)
    if exist(potential_files{i}, 'file')
      xml_file = potential_files{i};
      fprintf('Found JHU-tract atlas XML file: %s\n', xml_file);
      break;
    end
  end
else
  % Default: Harvard-Oxford atlas
  potential_files = {
    fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-Cortical.xml'),
    fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.xml'),
    fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-Subcortical.xml')
    };
  
  xml_file = '';
  for i = 1:length(potential_files)
    if exist(potential_files{i}, 'file')
      xml_file = potential_files{i};
      fprintf('Found Harvard-Oxford atlas XML file: %s\n', xml_file);
      break;
    end
  end
end

% Check if XML file exists
if isempty(xml_file) || ~exist(xml_file, 'file')
  warning('Atlas label file not found for %s atlas. Using numeric labels instead.', atlas_type);
  labels = struct('map', containers.Map('KeyType', 'double', 'ValueType', 'char'));
  labels.atlas_type = atlas_type;
  return;
end

% Parse XML file
try
  % Create a map to store the label mappings
  label_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
  
  % Use xmlread to parse the XML file
  xDoc = xmlread(xml_file);
  
  % Get all label nodes
  label_nodes = xDoc.getElementsByTagName('label');
  fprintf('Found %d label nodes in XML file\n', label_nodes.getLength);
  
  % Process each label node
  for i = 0:label_nodes.getLength-1
    node = label_nodes.item(i);
    
    % Get index (parcel ID)
    index_attr = node.getAttributes.getNamedItem('index');
    if ~isempty(index_attr)
      index = str2double(char(index_attr.getValue));
      
      % Get name (the text content of the label node)
      name = char(node.getTextContent);
      name = strtrim(name); % Remove any leading/trailing whitespace
      
      % Store in map
      if ~isnan(index) && ~isempty(name)
        label_map(index) = name;
      end
    end
  end
  
  % Return the labels structure
  labels = struct('map', label_map);
  labels.atlas_type = atlas_type;
  fprintf('Successfully loaded %d atlas labels for %s atlas\n', label_map.Count, atlas_type);
  
catch ex
  warning(ex.identifier, 'Error parsing atlas XML file: %s. Using numeric labels instead.', ex.message);
  disp(['Error details: ' getReport(ex)]);
  labels = struct('map', containers.Map('KeyType', 'double', 'ValueType', 'char'));
  labels.atlas_type = atlas_type;
end
end
