
function getdatapath()::String

    dataset_alias = "AI4DBiological-Systems_NMR_data" # don't use spaces or 'strange' symbols like commas, colons, etc.
    archive_file_name = "nmr_physical_parameters_dataset.tar.gz" # the filename on the data repository that we download.
    url = "https://zenodo.org/record/8174261/files/nmr_physical_parameters_dataset.tar.gz?download=1"

    register(DataDep("$dataset_alias",
        """
        Dataset: Selected BMRB 1D 1H NMR data and physical chemistry values compiled from literature
        Author: Roy Chih Chung Wang
        License: [Creative Commons Attribution Non Commercial Share Alike 4.0 International](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
        DOI: 10.5281/zenodo.8174261

        Description:
        The experiments were downloaded from BMRB,
        The physical chemistry parameters data and sample configuration-related files were collected, organized, converted or manually entered in as JSON files by Roy.
        The physical chemistry parameters, i.e., chemical shift and J-coupling values, are from two sources:
        - GISSMO: (https://gissmo.bmrb.io/ acccessed Mar. 2023). See Dashti, et. al. DOI: 10.1021/acs.analchem.8b02660 and DOI: 10.1021/acs.analchem.7b02884 for more details.
        - From Govindaraju, et. al.'s work: DOI: 10.1002/1099-1492(200005)13:3<129::AID-NBM619>3.0.CO;2-V

        Please cite the data sources and this data repository if you find the contents helpful for your work. See the Zenodo DOI entry for more description.
        """,
        url
    ));

    #readdir(datadep"AI4DBiological-Systems NMR data") # have to manually type out the alias. Does not allow string variable substitution.
    local_dataset_archive_path = @datadep_str("$dataset_alias") # call the actual macro to allow string variable substitution.

    # extract archive, then delete. Do this only if archive file still exists.

    root_data_path = joinpath(local_dataset_archive_path, "contents")
    
    if isfile(joinpath(local_dataset_archive_path, archive_file_name))
        t = @task begin; ispath(root_data_path) || mkpath(root_data_path); end
        schedule(t); wait(t)
    
        #t = @task begin; Tar.extract(joinpath(local_dataset_archive_path, archive_file_name), root_data_path); end
        t = @task begin; extractuncompress(joinpath(local_dataset_archive_path, archive_file_name), root_data_path); end
        schedule(t); wait(t)
        rm(joinpath(local_dataset_archive_path, archive_file_name)) # delete the archive file.
    end

    return root_data_path

    # # return root_data_path. however, this unpacks in the current working directory!
    # archive_file_path = joinpath(local_dataset_archive_path, archive_file_name)
    # if isfile(archive_file_path)
    #     DataDeps.unpack(archive_file_path)
    # end
    #return local_dataset_archive_path
end

function extractuncompress(src_path, dest_path)
    tar_gz = open(src_path)
    tar = GzipDecompressorStream(tar_gz)
    dir = Tar.extract(tar, dest_path)
    close(tar)
    
    return dir
end