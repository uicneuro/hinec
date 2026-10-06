classdef TestExternalAtlas < matlab.unittest.TestCase
    properties
        Folder
        Reference
        Atlas
        Labels
        Data
    end
    methods(TestMethodSetup)
        function setup(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(genpath(fullfile(root, 'src')));
            addpath(root);
            tc.Folder = tempname; mkdir(tc.Folder);
            tc.addTeardown(@()rmdir(tc.Folder, 's'));
            tc.Reference = fullfile(tc.Folder, 'dwi.nii');
            niftiwrite(ones(3,4,5,2,'single'), tc.Reference);
            tc.Atlas = fullfile(tc.Folder, 'atlas.nii');
            tc.Data = zeros(3,4,5,'int32');
            tc.Data(1,:,:) = 7; tc.Data(2,:,:) = 42;
            niftiwrite(tc.Data, tc.Atlas);
            tc.Labels = fullfile(tc.Folder, 'names.tsv');
            tc.writeLabels(sprintf('index\tname\n0\tBackground\n7\tCallosal subset\n42\tOther region\n'));
        end
    end
    methods(Test)
        function stagesBothNiftiFormatsAndPreservesIDs(tc)
            for compressed = [false true]
                niftiwrite(tc.Data, tc.Atlas, 'Compressed', compressed);
                source = tc.Atlas;
                if compressed, source = [source '.gz']; end
                dest = fullfile(tc.Folder, 'parcellation.nii.gz');
                [P, labels] = nim_prepare_external_atlas(source, tc.Reference, dest, tc.Labels);
                tc.verifyEqual(P, tc.Data);
                tc.verifyEqual(niftiread(dest), tc.Data);
                tc.verifyEqual(labels.map(7), 'Callosal subset');
                tc.verifyEqual(labels.atlas_type, 'external');
                tc.verifyEqual(niftiinfo(dest).Transform.T, niftiinfo(source).Transform.T);
                fid = fopen(dest); c = onCleanup(@()fclose(fid));
                tc.verifyEqual(fread(fid,2,'uint8'), [31;139]);
                clear c;
            end
        end
        function resolvesCustomNamesAndNumericIDs(tc)
            [P,L] = nim_prepare_external_atlas(tc.Atlas, tc.Reference, '', tc.Labels);
            nim = struct('parcellation_mask',P,'atlas_labels',L,'atlas_type','external');
            by_name = nim_roi_mask(nim, {'Callosal subset'});
            by_id = nim_roi_mask(nim, {7});
            tc.verifyEqual(by_name, tc.Data == 7);
            tc.verifyEqual(by_name, by_id);
        end
        function generatesNumericNamesWithoutHumanFallback(tc)
            [~,L] = nim_prepare_external_atlas(tc.Atlas, tc.Reference);
            tc.verifyEqual(L.map(7), 'Region_7');
            tc.verifyEqual(L.map(42), 'Region_42');
            tc.verifyEqual(double(L.map.Count), 2);
        end
        function rejectsShiftedAffineDespiteMatchingDimensions(tc)
            info = niftiinfo(tc.Atlas); T = info.Transform.T; T(4,1) = T(4,1)+1;
            info.Transform = affine3d(T); info.TransformName = 'Sform';
            niftiwrite(tc.Data, tc.Atlas, info);
            tc.verifyError(@()nim_prepare_external_atlas(tc.Atlas, tc.Reference), 'nim:externalAtlasGrid');
        end
        function rejectsDifferentDimensions(tc)
            niftiwrite(ones(4,4,5,'int32'),tc.Atlas);
            tc.verifyError(@()nim_prepare_external_atlas(tc.Atlas,tc.Reference), 'nim:externalAtlasGrid');
        end
        function rejectsProbabilitiesNegativeAndNonfiniteLabels(tc)
            for value = [.5 -1 NaN Inf double(intmax('int32'))+1]
                data = double(tc.Data); data(1) = value;
                niftiwrite(data,tc.Atlas);
                tc.verifyError(@()nim_prepare_external_atlas(tc.Atlas,tc.Reference), 'nim:externalAtlasValues');
            end
        end
        function rejectsEmptyAtlas(tc)
            niftiwrite(zeros(3,4,5,'int32'),tc.Atlas);
            tc.verifyError(@()nim_prepare_external_atlas(tc.Atlas,tc.Reference), 'nim:externalAtlasEmpty');
        end
        function rejectsIncompleteAndAmbiguousLabelTables(tc)
            cases = {sprintf('index\tname\n7\tOnly one\n'), ...
                sprintf('index\tname\n7\tA\n7\tB\n42\tC\n'), ...
                sprintf('index\tname\n7\tSame\n42\tsame\n'), ...
                sprintf('id\tname\n7\tA\n42\tB\n'), ...
                sprintf('index\tname\n7\tA\n42\t\n')};
            for i = 1:numel(cases)
                tc.writeLabels(cases{i});
                tc.verifyError(@()nim_prepare_external_atlas(tc.Atlas,tc.Reference,'',tc.Labels), 'nim:externalAtlasLabels');
            end
        end
        function explicitMissingFilesDoNotFallBack(tc)
            tc.verifyError(@()nim_prepare_external_atlas('missing.nii',tc.Reference), 'nim:externalAtlasMissing');
            tc.verifyError(@()nim_prepare_external_atlas(tc.Atlas,tc.Reference,'','missing.tsv'), 'nim:externalAtlasLabels');
        end
        function preprocessingExternalBranchNeedsNoFslOrHumanAtlas(tc)
            old = getenv('FSLDIR'); tc.addTeardown(@()setenv('FSLDIR',old)); setenv('FSLDIR','');
            options = struct('atlas_file',tc.Atlas,'atlas_labels_file',tc.Labels, ...
                'use_t1_registration',true,'t1_available',true);
            [path, labelpath] = preproc_atlas_resampling(tc.Reference,tc.Folder, ...
                fullfile(tc.Folder,'sample'),'jhu',options);
            tc.verifyEqual(niftiread(path),tc.Data);
            saved = load(labelpath);
            tc.verifyEqual(saved.atlas_labels.map(42),'Other region');
            tc.verifyEqual(saved.atlas_labels.source,tc.Atlas);
        end
        function configAndOverridesPreserveExternalPaths(tc)
            configfile = fullfile(tc.Folder,'external.yml');
            fid = fopen(configfile,'w');
            fprintf(fid,'preprocessing:\n  atlas_file: Some/MacaqueAtlas.nii.gz\n  atlas_labels_file: Some/Names.tsv\n');
            fclose(fid);
            config = load_config_yaml(configfile);
            tc.verifyEqual(config.preprocessing.atlas_file,'Some/MacaqueAtlas.nii.gz');
            tc.verifyEqual(config.preprocessing.atlas_labels_file,'Some/Names.tsv');
            config = nim_config_apply_overrides(config,{'preprocessing.atlas_file=Other/Atlas.nii'});
            tc.verifyEqual(config.preprocessing.atlas_file,'Other/Atlas.nii');
        end
        function cachedOutputCannotSilentlyIgnoreNewAtlas(tc)
            configfile = fullfile(tc.Folder,'external.yml');
            fid = fopen(configfile,'w'); fprintf(fid,'preprocessing:\n  atlas_file: supplied.nii\n'); fclose(fid);
            config = load_config_yaml(configfile);
            out = fullfile(tc.Folder,'cached.mat'); dummy=1; save(out,'dummy');
            tc.verifyError(@()main('unused',out,config),'main:externalAtlasCachedOutput');
        end
        function labelFileRequiresAtlas(tc)
            configfile = fullfile(tc.Folder,'external.yml');
            fid=fopen(configfile,'w');fprintf(fid,'preprocessing:\n  atlas_labels_file: names.tsv\n');fclose(fid);
            config = load_config_yaml(configfile);
            tc.verifyError(@()main('unused','unused.mat',config),'nim:externalAtlasRequired');
        end
        function mainRawAndCachedRoutesKeepExternalAnatomy(tc)
            % Exercise main's raw/cached routing, run-directory staging, flags,
            % and saved names while stubbing expensive image estimation only.
            olddir = pwd; oldpath = path;
            tc.addTeardown(@()restore_environment(olddir,oldpath));
            cd(tc.Folder);
            for d = {'src/nim_preprocessing','src/nim_plots','src/nim_utils', ...
                    'src/nim_calculation','src/nim_parcellation','src/nim_tractography', ...
                    'src/nim_registration','tests','lib/spm12','lib/bfgs','data','stubs'}
                mkdir(d{1});
            end
            tc.writeStub('nim_read',sprintf('function n=nim_read(varargin)\nn=struct(''FA'',ones(3,4,5));\nend\n'));
            for name = {'nim_dt_spd','nim_eig','nim_fa'}
                tc.writeStub(name{1},sprintf('function n=%s(n)\nend\n',name{1}));
            end
            for name = {'nim_registration','nim_parcellation','nim_parcellation_registered','nim_load_labels'}
                tc.writeStub(name{1},sprintf('function varargout=%s(varargin)\nerror(''test:humanFallback'',''Unexpected anatomy fallback'');\nend\n',name{1}));
            end
            tc.writeStub('nim_preprocessing',sprintf([ ...
                'function nim_preprocessing(prefix,options)\n' ...
                'assert(~options.use_t1_registration);\n' ...
                'copyfile([prefix ''_raw.nii.gz''],[prefix ''.nii.gz'']);\n' ...
                'preproc_atlas_resampling([prefix ''.nii.gz''],fileparts(prefix),prefix,options.atlas_type,options);\nend\n']));
            addpath(fullfile(tc.Folder,'stubs'),'-begin');
            configfile = fullfile(tc.Folder,'external.yml');
            fid=fopen(configfile,'w');
            fprintf(fid,['preprocessing:\n  atlas_file: %s\n  atlas_labels_file: %s\n' ...
                '  t1_available: true\n  use_t1_registration: false\n  register_to_mni: false\n'],tc.Atlas,tc.Labels);
            fclose(fid); config=load_config_yaml(configfile);
            for cached = [false true]
                inputdir = fullfile(tc.Folder,sprintf('case%d',cached)); mkdir(inputdir);
                prefix = fullfile(inputdir,'sample');
                niftiwrite(ones(3,4,5,2,'single'),[prefix '_raw.nii'],'Compressed',true);
                niftiwrite(ones(3,4,5,'single'),[prefix '_T1.nii'],'Compressed',true);
                if cached, copyfile([prefix '_raw.nii.gz'],[prefix '.nii.gz']); end
                info = struct('run_dir',fullfile(inputdir,'run'),'run_id',sprintf('external_case%d',cached));
                for key = {'intermediate_dir','output_dir','logs_dir','tractography_dir'}
                    info.(key{1})=fullfile(info.run_dir,key{1});mkdir(info.(key{1}));
                end
                % This stale human sidecar must not replace external names.
                fid=fopen(fullfile(info.intermediate_dir,'atlas_labels.xml'),'w');
                fprintf(fid,'<atlas><label index="7">Wrong human region</label></atlas>');fclose(fid);
                main(prefix,fullfile(inputdir,'result.mat'),config,info);
                saved=load(fullfile(info.output_dir,'result.mat'));
                tc.verifyEqual(saved.nim.parcellation_mask,tc.Data);
                tc.verifyEqual(saved.nim.atlas_labels.map(7),'Callosal subset');
                tc.verifyEqual(saved.nim.atlas_type,'external');
                tc.verifyTrue(isfile(saved.nim.parcellation_mask_file));
                tc.verifyTrue(isfile(saved.nim.atlas_labels_file));
            end
        end
    end
    methods(Access=private)
        function writeLabels(tc, text)
            fid=fopen(tc.Labels,'w');fprintf(fid,'%s',text);fclose(fid);
        end
        function writeStub(tc, name, text)
            fid=fopen(fullfile(tc.Folder,'stubs',[name '.m']),'w');fprintf(fid,'%s',text);fclose(fid);
        end
    end
end

function restore_environment(folder, saved_path)
cd(folder);path(saved_path);
end
