using ArgParse
using TreeTools
using TreeKnit

function parse_cmd()
    s = ArgParseSettings()
    @add_arg_table! s begin
        #"method"
        #help= "method of treeknit"
        #arg_type = String
        #required = true
        "tree_S"
        help = "path to tree S file"
        arg_type = String
        required = true
        "tree_M"
        help = "path to tree M file"
        arg_type = String
        required = true
        "tree_L"
        help = "path to tree L file"
        arg_type = String
        required = true
        "--output_S"
        help = "output path for tree S resolved"
        arg_type = String
        default = "results/tree_S_resolved.nwk"
        "--output_M"
        help = "output path for tree M resolved"
        arg_type = String
        default = "results/tree_M_resolved.nwk"
        "--output_L"
        help = "output path for tree L resolved"
        arg_type = String
        default = "results/tree_L_resolved.nwk"
    end
    return parse_args(s)
end

function main()
	args = parse_cmd()
    t_S = read_tree(args["tree_S"], label="S")
    t_M = read_tree(args["tree_M"], label="M")
    t_L = read_tree(args["tree_L"], label="L")
    MCCs = run_treeknit!([t_S, t_M, t_L], OptArgs(3;method=:better_trees))

    write_newick(args["output_S"], t_S)
    write_newick(args["output_M"], t_M)
    write_newick(args["output_L"], t_L)
end

main()
