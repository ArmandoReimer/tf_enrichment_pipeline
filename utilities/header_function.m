function [RawResultsPath, DataPath, FigureRoot] =   header_function(DropboxFolder, project) 
%     RawResultsPath = [DropboxFolder 'LocalEnrichmentResults\'];
    RawResultsPath = DropboxFolder;
    DataPath = [DropboxFolder 'ProcessedEnrichmentData\' project '\'];
    FigureRoot = [DropboxFolder 'LocalEnrichmentFigures\PipelineOutput\'];
    mkdir(FigureRoot);