function [mat, edf, aux] = foragingLoad(subj, ses, folder, prm)
% Importa e retorna os arquivos mat e edf associados ao sujeito e sessão
    aux = [sprintf('%02d', subj) '_' sprintf('%02d', ses)];
    st  = {[prm.matPreffix aux], [prm.edfPreffix aux]};
    ext = {prm.matExtension,  prm.edfExtension};

    fileNames = foragingFindFiles(folder, st, ext);

    mat = load(fileNames{1});
    oldFolder = pwd;
    [newFolder, edfFile, edfExt] = fileparts(fileNames{2});
    cd(newFolder);
    edf = edfmex([edfFile edfExt]);
    cd(oldFolder);
end
