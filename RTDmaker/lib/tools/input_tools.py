import os
import sys
import glob
import shutil


def check_input(args):

    # Arguments order and grouping:
    # srQC.main(args.ass_dir, args.ref_dir,                       # Annotations to analyze
    #           args.sj_dir, args.sjreads_th,                     # SJ-QC related args
    #           args.genome_fl, args.fastq_dir, args.abund_th,    # Abundance-QC related args
    #           args.len_th, args.antlen_th,                      # Structure-QC related args
    #           args.skip_set, args.add_set,                      # Category-selection related args
    #           args.size_th, args.keep_set,                      # Memory-related arguments
    #           paths_dt, args.prefix, args.outname, logfile)     # Output related args

    # 1) Check the annotations provided by the user
    if not args.ass_dir and not args.ref_dir:
        sys.exit(f'ERROR: You must specify either assemblies or references to analyze, both are missing')

    if args.ass_dir:
        ass_files = sorted(glob.glob(os.path.join(args.ass_dir, "**/*.gtf"), recursive=True))
        if not ass_files:
            sys.exit(f'ERROR: The specified assembly directory ({args.ass_dir}) seems to be empty')

    if args.ref_dir:
        ref_files = sorted(glob.glob(os.path.join(args.ref_dir, "**/*.gtf"), recursive=True))
        if not ref_files:
            sys.exit(f'ERROR: The specified reference directory ({args.ref_dir}) seems to be empty')

    # 2) Check SJ-QC args (check if user is not skipping this analysis)
    if args.ass_dir and "SJ" not in args.skip_set:
        if not args.ass_dir:
            sys.exit(f'ERROR: No assembly directory specified')

        sj_files = sorted(glob.glob(os.path.join(args.sj_dir, "**/*SJ.out.tab"), recursive=True))
        if not sj_files:
            sys.exit(f'ERROR: The specified SJ directory ({args.sj_dir}) seems to be empty')

        if not args.sjreads_th or len(args.sjreads_th) != 2:
            sys.exit(f'ERROR: Incorrect number of values specified for --SJ-reads')

        # This also to convert the values to floats
        try:
            args.sjreads_th = tuple(map(float, args.sjreads_th))
        except ValueError:
            sys.exit(f"ERROR: One of --SJ-reads values is not a valid number: {args.sjreads_th}")

        for val in args.sjreads_th:
            if val <= 0:
                sys.exit(f'ERROR: Invalid value ({val}) specified --SJ-reads. Value must be positive.')

    # This check if the program exist, but it doesn't check if the program runs correctly, thus anal can still fail
    is_tool = lambda name: shutil.which(name) is not None

    # 3) Check Abundance-QC args (check if user is not skipping this analysis)
    if args.ass_dir and "quant" not in args.skip_set:
        # Check if Salmon and Gffread are installed in the environment
        if is_tool("salmon") is False:
            sys.exit(f'ERROR: Salmon is not installed in the environment.')

        if is_tool("gffread") is False:
            sys.exit(f'ERROR: GffRead is not installed in the environment.')

        if not args.genome_fl:
            sys.exit(f'ERROR: No Genome Fasta file specified')

        if not args.genome_fl.endswith(".fa") and not args.genome_fl.endswith(".fasta"):
            sys.exit(f'ERROR: Genome file must be in Fasta format (prefix must be .fa or .fasta)')

        if not os.path.exists(args.genome_fl):
            sys.exit(f'ERROR: Genome file does not exit: {args.genome_fl}')

        if not args.ass_dir:
            sys.exit(f'ERROR: No assembly directory specified')

        if not args.fastq_dir:
            sys.exit(f'ERROR: No Fastq directory specified')

        fastq_files = []
        for ext in ['fq.gz', 'fq', 'fastq.gz', 'fastq']:
            for fl in sorted(glob.glob(os.path.join(args.fastq_dir, f"**/*.{ext}"), recursive=True)):
                fastq_files.append(fl)

        if not fastq_files:
            sys.exit(f'ERROR: The specified Fastq directory ({args.sj_dir}) seems to be empty')

        # This also to convert the values to floats
        try:
            args.abund_th = tuple(map(float, args.abund_th))
        except ValueError:
            sys.exit(f"ERROR: One of --tpm values is not a valid number: {args.abund_th}")

        if not args.abund_th or len(args.abund_th) != 2:
            sys.exit(f'ERROR: Incorrect number of values specified for --tpm')

        for val in args.abund_th:
            if val <= 0:
                sys.exit(f'ERROR: Invalid value ({val}) specified --tpm. Value must be positive.')

    # GffRead now is also required at the end of the pipeline
    if is_tool("gffread") is False:
        sys.exit(f'ERROR: GffRead is not installed in the environment.')

    # 4) Check the values for the Structural-QC
    if not args.len_th or not 0 <= args.len_th <= 1:
        sys.exit(f'ERROR: Incorrect value for fragmentary check (--fragment-len): {args.len_th}. '
                 f'Value must be between 0 and 1.')

    if not args.antlen_th or not 0 <= args.antlen_th <= 1:
        sys.exit(f'ERROR: Incorrect value for antisense fragmentary check (--antisense-len): {args.antlen_th}. '
                 f'Value must be between 0 and 1.')

    # 5) Check optional models by user
    if "None" in args.add_set:
        if len(args.add_set) != 1:
            print(f"WARNING: 'None' was specified as one of multiple arguments for --add. Thus, all of the optional models will be removed.")
        args.add_set = set()
    else:
        args.add_set = set(args.add_set)

    if "None" in args.keep_set:
        if len(args.keep_set) != 1:
            sys.exit(f"ERROR: The 'None' option is mutually exclusive with the other options for the --keep argument.")
        args.keep_set = set()
    else:
        args.keep_set = set(args.keep_set)

    # 6) Check output related args
    if not args.outpath:
        sys.exit(f'ERROR: Path for the output folder is missing')

    if not args.outname:
        sys.exit(f'ERROR: Name for the output files is missing')

    if not args.prefix:
        args.prefix = args.outname

    # Filename should not contain an extension as it is used as a prefix for additional files
    if args.outname.endswith(".gtf"):
        args.outname = args.outname.replace(".gtf", "")

    if args.size_th < 1:
        sys.exit(f'ERROR: The maximum file size allocated for temporary files (--ram) cannot be less than 1 Gb')
    elif args.size_th > 8:
        sys.exit(f'ERROR: The maximum file size allocated for temporary files (--ram) cannot be more than 8 Gb')
    else:
        pass

    # Sanitize paths
    # Convert relative paths to absolute, and convert path to the type used by local OS system (Linux vs Windows paths)
    if args.ass_dir:
        args.ass_dir = os.path.abspath(os.path.normpath(args.ass_dir))

    if args.ref_dir:
        args.ref_dir = os.path.abspath(os.path.normpath(args.ref_dir))

    if args.sj_dir:
        args.sj_dir = os.path.abspath(os.path.normpath(args.sj_dir))

    if args.fastq_dir:
        args.fastq_dir = os.path.abspath(os.path.normpath(args.fastq_dir))

    if args.outpath:
        args.outpath = os.path.abspath(os.path.normpath(args.outpath))

    return args


def create_project_dir(outpath):

    # 3) Generate the whole project directory structure and create a dictionary to handle the multiple subfolder paths
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    # Generate sanitized sub-folders paths
    dir_out = outpath
    dir_inter = os.path.abspath(os.path.normpath(os.path.join(dir_out, "intermediary")))
    dir_removed = os.path.abspath(os.path.normpath(os.path.join(dir_out, "removed")))
    dir_report = os.path.abspath(os.path.normpath(os.path.join(dir_out, "report")))

    dir_log = os.path.abspath(os.path.normpath(os.path.join(dir_report, "logfiles")))
    dir_chm = os.path.abspath(os.path.normpath(os.path.join(dir_inter, "additional")))
    dir_high = os.path.abspath(os.path.normpath(os.path.join(dir_inter, "abundant")))
    dir_quant = os.path.abspath(os.path.normpath(os.path.join(dir_inter, "quantification")))

    dir_low = os.path.abspath(os.path.normpath(os.path.join(dir_removed, "low_abundant")))

    # Create the paths and and track their location to facilitate their handling inside the code
    paths = [("outpath", dir_out), ("inter", dir_inter), ("removed", dir_removed), ("report", dir_report),
             ("additional", dir_chm), ("low_abundant", dir_low), ("abundant", dir_high), ("quant", dir_quant),
             ("logs", dir_log)]

    paths_dt = {}
    for (dir_name, dir_path) in paths:
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        paths_dt[dir_name] = dir_path

    return paths_dt
