function [RawResultsPath, DataPath, FigureRoot] =   header_function(DropboxFolder, project) 
%     RawResultsPath = [DropboxFolder 'LocalEnrichmentResults\'];
    RawResultsPath = [DropboxFolder, filesep];
    DataPath = [DropboxFolder, filesep, 'ProcessedEnrichmentData\' project '\'];
    FigureRoot = [DropboxFolder, filesep, 'LocalEnrichmentFigures\PipelineOutput'];

    mkdir(FigureRoot);