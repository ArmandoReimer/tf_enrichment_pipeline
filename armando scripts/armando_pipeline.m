
dataType = '1Dg11_2xDl';
dataType = '1DgS2_2xDl';

addpath(genpath('S:\Armando\tf_enrichment_pipeline'));

main01_compile_traces(dataType,'S:\Armando\Dropbox\DorsalSyntheticsDropbox\','firstNC', 12)

main02_sample_local_protein(dataType,'S:\Armando\Dropbox\DorsalSyntheticsDropbox\',...
    'firstNC', 12, 'masking_method', 'kSnakeCirles', 'ignoreQC', true, 'write_snip_flag', true,...
 'shouldSegmentNuclei', true, 'max_rad_um', 6)

main04_make_exploratory_figs(dataType,'S:\Armando\Dropbox\DorsalSyntheticsDropbox\',...
    'firstNC', 12)