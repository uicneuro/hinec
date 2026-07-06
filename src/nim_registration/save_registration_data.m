function save_registration_data(registration_data, options)
% save_registration_data: Save registration data to file

output_file = registration_data.output_file;

fprintf('Saving registration data...\n');

try
    save(output_file, 'registration_data', '-v7.3');

    file_info = dir(output_file);
    fprintf('  ✓ Registration data saved: %s (%.1f MB)\n', output_file, file_info.bytes/1024/1024);

catch ME
    error('Failed to save registration data: %s', ME.message);
end

end
