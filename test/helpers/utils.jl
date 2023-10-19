
############ serialization testing

function saveloadBSON(file_path, S)

    BSON.bson(file_path, S)
    return BSON.load(file_path)
end

function saveloadJSON(file_path, S)

    NMRHamiltonian.saveasJSON(file_path, S)
    return JSON3.read(read(file_path))
end

function cleanup(BSON_filename, JSON_filename, clean_up_files)
    
    if clean_up_files
        if isfile(BSON_filename)
            rm(BSON_filename)
        end
        
        if isfile(JSON_filename)
            rm(JSON_filename)
        end

    end

    return nothing
end