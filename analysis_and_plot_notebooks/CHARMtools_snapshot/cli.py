import argparse

import clean_leg
import clean_splicing
import clean_isolated
import sep_clean
import clean3

def cli():
    parser = argparse.ArgumentParser(prog="CHARMtools", description="Functions for data-analysis in HiRES etc. projects")
    subcommands = parser.add_subparsers(title="These are sub-commands",metavar="command")
#--------- clean_leg sub command ---
    clean_leg_arg = subcommands.add_parser(
                            "clean_leg",
                            help="clean promiscuous legs that contacts with multiple legs")
    clean_leg_arg.set_defaults(handle=clean_leg.cli)
    clean_leg_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            nargs=1)
    clean_leg_arg.add_argument(
                            "-t", "--thread",
                            type=int,
                            dest="thread",
                            action="store",
                            default=4,
                            help="set thread number")
    clean_leg_arg.add_argument(
                            "-d","--distance",
                            dest="max_distance",
                            metavar="MAX_DISTANCE",
                            type=int,
                            action="store",
                            default=1000,
                            help="max distance to calculate adjacent legs"
    )
    clean_leg_arg.add_argument(
                            "-n","--count",
                            metavar="MAX_COUNT",
                            dest="max_count",
                            type=int,
                            action="store",
                            default=10,
                            help="number threshold of adjacent legs"
    )                   
    clean_leg_arg.add_argument("-o", "--output", 
                            dest="out_name", action="store",
                            metavar="OUTPUT_FILE",
                            required=True,
                            help="set output file name")
 # ---------#clean_splicing sub command ---
    clean_splicing_arg = subcommands.add_parser(
                            "clean_splicing", 
                            help="clean exon splicing from mRNA in contact file")
    clean_splicing_arg.set_defaults(handle=clean_splicing.cli)
    clean_splicing_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1)     
    clean_splicing_arg.add_argument(
                            "-r", "--reference", 
                            dest="gtf_filename",
                            type = str,
                            action="store", 
                            help="annotation gtf file", 
                            required=True)
    clean_splicing_arg.add_argument(
                            "-o", "--output", 
                            dest="out_name",
                            metavar="OUTPUT_FILE",
                            required=True, 
                            help="output file name", 
                            action="store")
    clean_splicing_arg.add_argument(
                            "-t", "--thread",
                            type=int,
                            dest="num_thread",
                            action="store",
                            default=4,
                            help="set thread number")

#--------- clean3 subcommand ------
    clean3_arg = subcommands.add_parser(
                            "clean3",
                            help="clean 3dg particles poorly supported by contacts"
    )
    clean3_arg.set_defaults(handle=clean3.cli)
    clean3_arg.add_argument(
                            "-i","--input",
                            dest="filename",
                            metavar="STRUCTURE_FILE",
                            help=".3dg/.xyz format structure file to clean",
                            type=str,
                            required=True
    )
    clean3_arg.add_argument(
                            "-r", "--reference",
                            dest="ref_filename",
                            metavar="PAIRS",
                            help=".pairs format contact file",
                            type=str,
                            required=True
    )
    clean3_arg.add_argument(
                            "-o", "--output",
                            dest="output",
                            metavar="CLEANED",
                            help="file name of the output cleaned structure file",
                            type=str,
                            required=True
    )
    clean3_arg.add_argument(
                            "-q", "--quantile",
                            dest="quantile",
                            metavar="QUANTILE",
                            help="quantile of particles to remove",
                            type=int,
                            default=0.06
    )
    clean3_arg.add_argument(
                            "-d", "--distance",
                            dest="distance",
                            metavar="DISTANCE",
                            help="max distance (bp) from a contact leg to a 3D genome particle",
                            type=int,
                            default=500_000
    )
    
#--------- clean_isolate subcommand ------
    clean_isolated_arg = subcommands.add_parser(
                            "clean_isolated",
                            help="remove isolated contacts according to L-0.5 distance"
    )
    clean_isolated_arg.set_defaults(handle=clean_isolated.cli)
    clean_isolated_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1
    )
    clean_isolated_arg.add_argument(
                            "-t","--thread",
                            dest="thread",
                            type=int,
                            help="set thread number",
                            default=4
    )     
    clean_isolated_arg.add_argument(
                            "-m","--dense",
                            dest="dense",
                            type=int,
                            help="number of contacts in proximity",
                            default=5)
    clean_isolated_arg.add_argument(
                            "-d","--distance",
                            dest="distance",
                            type=int,
                            help="check contacts in what L-0.5 range",
                            default=10000000) 
    clean_isolated_arg.add_argument(
                            "-o","--output",
                            dest="output_file",
                            action="store",
                            metavar="OUTPUT_FILE",
                            required=True,
                            help = "output file name",
                            type=str
    )

# --------- sep_clean subcommand ---------
    sep_clean_arg = subcommands.add_parser(
                            "sep_clean",
                            help = "seperate homologous chromosome(add (mat)/(pat) for chra chrb colomns), clean isolated contacts again. \
                                generate one more hickit compatible output file.\
                                    works with hickit imputed pairs file"
    )
    sep_clean_arg.set_defaults(handle=sep_clean.cli)
    sep_clean_arg.add_argument(
                            dest="filename",
                            metavar="INPUT_FILE",
                            help="input filename",
                            nargs=1
    )
    sep_clean_arg.add_argument(
                            "-n", "--num_thread",
                            dest="num_thread",
                            help="number of thread use",
                            default="4"
    )
    sep_clean_arg.add_argument(
                            "-o1", "--output1",
                            dest="output_file1",
                            help="output file path for .dip.pairs (the chra(mat) format) file",
                            action="store",
                            required=True
    )
    sep_clean_arg.add_argument(
                            "-o2", "--output2",
                            dest="output_file2",
                            help="output file path for .pairs (the hickit -b default format) file",
                            required=True
    )
    sep_clean_arg.add_argument(
                            "-m","--dense",
                            dest="dense",
                            type=int,
                            help="number of contacts in proximity",
                            default=5
    )
    sep_clean_arg.add_argument(
                            "-d","--distance",
                            dest="distance",
                            type=int,
                            help="check contacts in what L-0.5 range",
                            default=10000000
    ) 

    args = parser.parse_args()
    #print(args.replace_switch)
    #print(args.out_name)
    #print(args.filenames)
    if hasattr(args, "handle"):
        args.handle(args)
    else:
        parser.print_help()
if __name__ == "__main__":
    cli()