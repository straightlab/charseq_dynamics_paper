using ArgParse
using BridgeFinder
#

function parse_commandline()
    s = ArgParseSettings(description= "Splice .FR reads for CHarSeq data")
    @add_arg_table s begin
        "--matchParamsRNA", "-r"
            help = "Maximum number of mismatches, min number of matches, and max lookup for RNA"
            arg_type = String
            default = "2:10:20"
		"--matchParamsDNA", "-d"
            help = "Maximum number of mismatches, min number of matches, and max lookup for RNA"
            arg_type = String
            default = "2:10:20"
		"--maxreadl", "-l"
            help = "Maximum read length"
            arg_type = Int
            default = 294
        "--verbose", "-v"
            help = "verbose, prints progress"
            action = :store_true
        "--checkpairs", "-c"
            help = "make sure reads are properly paired in input files"
            action = :store_true
        "fq_prefix"
            help = "prefix of input files"
            required = true
		"fq_spliced_prefix"
            help = "prefix of output files for spliced reads"
            required = true
		"fq_unspliced_prefix"
            help = "prefix of output files for unspliced reads"
            required = true


    end
    parsed_args = parse_args(s)
end

function main()
    dict = parse_commandline()
	splice_mates(dict["fq_prefix"],dict["fq_spliced_prefix"],dict["fq_unspliced_prefix"],dict["matchParamsRNA"],dict["matchParamsDNA"],dict["maxreadl"],dict["checkpairs"], dict["verbose"])
end

main()
