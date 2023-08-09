#include "metadata_map.h"

void metadata_map::build_map(std::unordered_map<std::string, std::string> *meta_map, std::string metadata_fname)
{
    std::vector<std::string> metadata_options;
    boost::split(metadata_options, metadata_fname, boost::is_any_of(","));
    std::ifstream metadata_file(metadata_options[0], std::ios_base::in);

    if(!metadata_file.is_open())
    {
        throw std::runtime_error("File could not be opened. Verify metadata file exists.\n");
    }
    if(metadata_options.size() < 3)
    {
        throw std::runtime_error("Missing required specifications for meta flag. Must include metadata file name, "
                                "sequence name column, and species identification column.\n");
    }

    std::string line;
    std::getline(metadata_file, line);
    std::for_each(++metadata_options.begin(), metadata_options.end(), [&](std::string header)
    {
        // verify metadata file header
        if(line.find(header) == std::string::npos)
        {
            throw std::runtime_error("Header '" + header + "' could not be found in metadata file. "
                                    "Verify names being entered to names in metadata file.\n");
        }
    });

    std::vector<std::string> metadata_header_row;
    std::size_t name_index = 0;
    std::size_t spec_index = 0;
    boost::split(metadata_header_row, line, boost::is_any_of("\t"));
    std::size_t count = 0;
    for(const auto& column_val : metadata_header_row)
    {
        if(column_val.compare(metadata_options[1]) == 0)
            name_index = count;
        if(column_val.compare(metadata_options[2]) == 0)
            spec_index = count;

        count++;
    }

    std::size_t missing_spec_id_count = 0;
    std::vector<std::string> metadata_row;
    std::ofstream ex_seqs_log(
        "excluded_protein_sequences.txt", std::ios::out
    );
    while(std::getline(metadata_file, line))
    {
        boost::split(metadata_row, line, boost::is_any_of("\t"));

        // check sequence has a species ID
        if (!metadata_row.at(spec_index).empty())
        {
            meta_map->emplace(metadata_row.at(name_index), metadata_row.at(spec_index));
        }
        // otherwise, assume no species ID
        else
        {
            ex_seqs_log << metadata_row.at(name_index) << "\n";
            missing_spec_id_count += 1;
        }
    }

    std::cout << "WARNING: " << missing_spec_id_count
        << " sequences in metadata file did not have a value for their"
        << " \"SpeciesID\" column - they have been excluded from this run.\n"
        << "Please review the aforementioned sequences in"
        << " \"excluded_protein_sequences.txt\"\n";
}

std::string metadata_map::get_id(std::string sequence_data, std::unordered_map<std::string, std::string> *meta_map)
{
    std::unordered_map<std::string, std::string>::const_iterator name = meta_map->find(sequence_data);

    if(name != meta_map->end())
    {
        return name->second;
    }
    else
    {
        return "";
    }
}

