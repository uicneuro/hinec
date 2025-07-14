# HINEC: HIgh-order NEural Connectivity - Guidelines

## Running Commands
- Add paths: `addpaths`
- Process data: `main <data_location> <output_file.mat>`
- Run complete workflow: Edit and execute `runhinec.m`
- Visualization: `nim_plotparcelall(nim)` or `nim_plotparcellation(nim)`
- Data loading: `load <filename>.mat`

## Code Style Guidelines
- **Naming**: Use camelCase for functions (e.g., `nim_preprocessing`, `nim_fa`)
- **Input validation**: Use MATLAB's `arguments` block for parameter validation
- **Indentation**: 2 spaces for consistent formatting
- **Documentation**: Comment each function with purpose and parameters
- **Error handling**: Use `error()` for fatal errors, informative messages
- **File organization**: Group related functions in logical directories (e.g., `nim_utils/`, `nim_calculation/`)
- **Output feedback**: Use `disp()` or `fprintf()` to provide processing status

## Requirements
- MATLAB with Image Processing and Statistics and Machine Learning Toolboxes
- SPM12 in root directory
- FSL (must be initialized before use)