#!/usr/bin/env Rscript 

# fetch GEO/SRA fastq files given a GEO/SRA accession
# generating shell scripts and run later
# supported sample_ids: GSM and SRR
# author: Ning Zhang (nzhang@stowers.org)



################################################### args ##################################################
# version
version = function() {
	return('0.1.5')
}

# args
get_opts = function() {
	suppressPackageStartupMessages(library(optparse))

	# args
	opts = list(
		make_option(c('-s','--sample_sheet'), type='character', default=NULL, metavar='file',
			help='Sample sheet csv file (required)'),
		make_option(c('-r','--requester'), type='character', default=NULL, metavar='string',
			help='Requester initial (required)'),
		make_option(c('-l','--lab'), type='character', default=NULL, metavar='string',
			help='Lab name (required)'),
		make_option(c('--out'), type='character', default='/n/core/Bioinformatics/PublicData', metavar='folder',
			help='Output folder (default: /n/core/Bioinformatics/PublicData)'),
		make_option(c('--threads'), type='numeric', default=5, metavar='int',
			help=paste('Maximum number of samples to download simultaneously (default: 5).',
			'Your connection may be blocked by admin if setting this number too high')),
		make_option(c('--temp'), type='character', default='/n/ngs/tools/secundo/scripts/temp_SRR_files', metavar='folder',
			help=paste('Cache Folder for NCBI SRA toolkit (default=/n/ngs/tools/secundo/scripts/temp_SRR_files/)',
			'Change it accordingly if your vdb-config was set up differently')), 
		make_option(c('--api_key'), type='character', default='75249f4227ae94555d3a37564ee87f0b7f08', metavar='key', 
			help='NCBI API key (default: 75249f4227ae94555d3a37564ee87f0b7f08)'), 
		make_option(c('--download_mode'), action='store_true', default=TRUE, 
			help='Turn on downloading mode (generating shell scripts and downloading) (default: enabled)'), 
		make_option(c('--dryrun_mode'), action='store_false', default=FALSE, dest='download_mode', 
			help='Turn on dry-run mode (only generating shell scripts without downloading) (default: disabled)')
	)

	opt_parser = OptionParser(option_list=opts, prog='Rscript sra_download.r',
		description=c(
			sprintf('\nSRA data download utility (version: %s).', version()),
			'Currently supported sample IDs: GSM and SRR.',
			'Contact Ning Zhang (nzhang@stowers.org) for any questions or comments.'
		)
	)
	opts = parse_args(opt_parser)

	# check args
	if(is.null(opts$sample_sheet)) {
		print_help(opt_parser)
		stop('Required option --sample_sheet is not set', call.=FALSE)
	}

	if(is.null(opts$lab)) {
		print_help(opt_parser)
		stop('Required option --lab is not set', call.=FALSE)
	}

	if(is.null(opts$requester)) {
		print_help(opt_parser)
		stop('Required option --requester is not set', call.=FALSE)
	}

	return(opts)
}
###########################################################################################################


############################################### sample info ###############################################
load_sample_info = function(sample_sheet_file, requester, lab) {
	suppressPackageStartupMessages(library(readxl))

	# load sample sheet
	# sample_sheet = read.xls(sample_sheet_file)
	sample_sheet = read_xlsx(sample_sheet_file)
	colnames(sample_sheet) = c('request_id', 'sample_id', 'sample_name', 'expr_type','read_type',
	'read_length', 'ref_genome','comments')


	# check expr_type
	if(!all(as.character(unique(sample_sheet['expr_type'])[[1]]) %in% c('Stranded RNA-seq', 'Non-stranded RNA-seq'))) {
		message('Incorrect sample sheet format:')
		stop('Invalid experiment types detected! Allowed values for expr_type column: Stranded RNA-seq or Non-stranded RNA-seq',
		call.=FALSE)
	}

	# check read_type
	if(!all(as.character(unique(sample_sheet['read_type'])[[1]]) %in% c('Single', 'Paired'))) {
		message('Incorrect sample sheet format:')
		stop('Invalid read types detected! Allowed values for read_type column: Single or Paired',
		call.=FALSE)
	}

	sample_info = list()

	sample_info$requester = requester
	sample_info$lab = lab
	sample_info$sample_sheet = sample_sheet

	return(sample_info)
}
###########################################################################################################


############################################ geo/sra download #############################################
# GSM to SRR accession
# one gsm id could be mapped to multiple srr ids
gsm_to_srr = function(gsm_id, api_key) {
	suppressPackageStartupMessages(library(rentrez))

	# returning matched srr ids
	srr_ids = NA

	# entrez search
	search = entrez_search(db='sra', term=gsm_id, api_key=api_key)
	num_match = length(search$ids)

	if(num_match == 0) { # no match found
		message(sprintf('[%s] No entrez match found for sample id: %s', Sys.time(), gsm_id))
	} else if(num_match > 1) { # multiple matches found
		message(sprintf('[%s] Multiple entrez matches found for sample id: %s', Sys.time(), gsm_id))
	} else { # unique match found
		# message(sprintf('[%s] Unique entrez match found for sample id: %s', Sys.time(), gsm_id))

		if(startsWith(gsm_id, 'GSM')) { # gsm
			summary = read.table(text=entrez_fetch(db='sra', id=search$ids[1], api_key=api_key, rettype='runinfo'), 
				sep=',', header=TRUE)
			srr_ids = as.character(summary$Run)
			message(sprintf('[%s] SRR IDs found: %s => %s', Sys.time(), gsm_id, paste(srr_ids, collapse=', ')))
		} else if(startsWith(gsm_id, 'SRR')) { # srr
			message(sprintf('[%s] Input ID is already SRR ID: %s', Sys.time(), gsm_id))
			srr_ids = gsm_id
		} else { # invalid format
			message(sprintf('[%s] Incorrect input ID: %s', Sys.time(), gsm_id))
		}
	}

	# sleep for 0.4 second due to ncbi api rate limit
	# 10 queries per second if api is used
	# 3 queries per second otherwise
	Sys.sleep(0.4)

	return(srr_ids)
}

# generate commands for fastq file download
download_fastq_cmds = function(gsm_id, srr_id, out_folder, sra_cache_folder) {
	# download sra file
	cmd1 = sprintf(
		'prefetch %s -O %s',
		srr_id, sra_cache_folder # I added an output folder here, or it will download to ./ 
	)

	# move sra file to output folder
	cmd2 = sprintf(
		'mv %s %s',
		file.path(sra_cache_folder, sprintf('%s/%s.sra', srr_id, srr_id)), # This was causing problems, so i added a SRR folder to the path
		file.path(out_folder, srr_id)
	)

	# validate downloaded sra file
	cmd3 = sprintf(
		'vdb-validate %s', 
		file.path(out_folder, srr_id, sprintf('%s.sra', srr_id))
	)

	# extract fastq file
	# results will be gzipped
	# paired-end sample will generate two files: _1.fastq and _2.fastq
	# single-end sample will generate one file:  _1.fastq
	cmd4 = sprintf(
		'fastq-dump --gzip --outdir %s --split-files %s',
		file.path(out_folder, srr_id),
		file.path(out_folder, srr_id, sprintf('%s.sra', srr_id))
	)

	# return cmds
	return(sprintf('# %s: %s\n%s\n%s\n%s\n%s\n', gsm_id, srr_id, cmd1, cmd2, cmd3, cmd4))
}

# wrapper
download_cmds_wrapper = function(request_id, sample_id, out_dir_base, sra_cache_folder, api_key) {
	run_info = data.frame()

	# out dir
	# create if not exists
	out_dir = file.path(out_dir_base, request_id, sample_id)
	if(!dir.exists(out_dir)) {
		dir.create(out_dir, recursive=TRUE)
	}

	# to SRR ID if necessary
	srr_ids = gsm_to_srr(sample_id, api_key)

	if(any(!is.na(srr_ids))) { # found srr id
		for(srr_id in srr_ids) {
			# create sub folder under out_dir if not exists
			if(!dir.exists(file.path(out_dir, srr_id))) {
				dir.create(file.path(out_dir, srr_id))
			}

			# generate commands
			cmds = download_fastq_cmds(sample_id, srr_id, out_dir, sra_cache_folder)
			run_info = rbind(run_info, data.frame(request_id=request_id, sample_id=sample_id, srr_id=srr_id, cmds=cmds))
		}
	} else { # srr id not found
		cmds = sprintf('# %s\n', sample_id)
		run_info = rbind(run_info, data.frame(request_id=request_id, sample_id=sample_id, srr_id=srr_ids, cmds=cmds))
	}

	return(run_info)
}

download_wrapper = function() {

}
###########################################################################################################

######################################### analysis sample sheet ###########################################
# helper function for get_sample_info_analysis
# get fastq path
get_fq_path = function(x, out_folder) {
	fq_path = NA

	request_id = x[1]
	sample_id = x[2]
	srr_id = x[3]
	read_type = x[4]

	if(!is.na(srr_id)) {
		if(read_type == 'Paired') {
			fq_1_path = file.path(out_folder, request_id, sample_id, srr_id, sprintf('%s_1.fastq.gz', srr_id))
			fq_2_path = file.path(out_folder, request_id, sample_id, srr_id, sprintf('%s_2.fastq.gz', srr_id))

			fq_path = sprintf('%s,%s', fq_1_path, fq_2_path)
		} else if (read_type == 'Single') {
			fq_path = file.path(out_folder, request_id, sample_id, srr_id, sprintf('%s_1.fastq.gz', srr_id))
		}
	} else {
		if(read_type == 'Paired') {
			fq_1_path = 'NA'
			fq_2_path = 'NA'

			fq_path = sprintf('%s,%s', fq_1_path, fq_2_path)
		} else if (read_type == 'Single') {
			fq_path = 'NA'
		}
	}

	return(fq_path)
}

# generate sample sheet for analysis
get_sample_info_analysis = function(sample_info, run_info, out_dir) {
	suppressPackageStartupMessages(library(data.table))

	out_sample_info = sample_info$sample_sheet

	# add additional cols
	out_sample_info = merge(out_sample_info, run_info)

	out_sample_info['lab'] = rep(sample_info$lab, nrow(out_sample_info))
	out_sample_info['requester'] = rep(sample_info$requester, nrow(out_sample_info))

	out_sample_info['fq_path'] = apply(
		out_sample_info[c('request_id', 'sample_id', 'srr_id', 'read_type')], 1,
		get_fq_path, out_folder=out_dir)

	# re-order cols
	out_sample_info = out_sample_info[c('fq_path', 'request_id', 'sample_id', 'sample_name', 'srr_id',
		'expr_type', 'read_type', 'read_length', 'ref_genome', 'lab','requester')]

	# split fastq for paired end
	# code could be improved
	out_sample_sheet = data.frame()

	for(i in seq(1, nrow(out_sample_info))) {
		if(out_sample_info['read_type'][i, ] == 'Paired') { # paired end
			fq_1_path = strsplit(out_sample_info[i, 'fq_path'], ',')[[1]][1]
			fq_2_path = strsplit(out_sample_info[i, 'fq_path'], ',')[[1]][2]

			out_sample_sheet = rbind(out_sample_sheet,
				data.frame(
					fq_path=fq_1_path,
					request_id=out_sample_info['request_id'][i, ],
					sample_id=out_sample_info['sample_id'][i, ],
					srr_id=out_sample_info['srr_id'][i, ],
					sample_name=out_sample_info['sample_name'][i, ],
					expr_type=out_sample_info['expr_type'][i, ],
					read=1,
					read_type=out_sample_info['read_type'][i, ],
					read_length=out_sample_info['read_length'][i, ] + 1,
					ref_genome=out_sample_info['ref_genome'][i, ],
					lab=out_sample_info['lab'][i, ],
					requester=out_sample_info['requester'][i, ]
				)
			)

			out_sample_sheet = rbind(out_sample_sheet,
				data.frame(
					fq_path=fq_2_path,
					request_id=out_sample_info['request_id'][i, ],
					sample_id=out_sample_info['sample_id'][i, ],
					srr_id=out_sample_info['srr_id'][i, ],
					sample_name=out_sample_info['sample_name'][i, ],
					expr_type=out_sample_info['expr_type'][i, ],
					read=2,
					read_type=out_sample_info['read_type'][i, ],
					read_length=out_sample_info['read_length'][i, ] + 1,
					ref_genome=out_sample_info['ref_genome'][i, ],
					lab=out_sample_info['lab'][i, ],
					requester=out_sample_info['requester'][i, ]
				)
			)
		} else if(out_sample_info['read_type'][i, ] == 'Single') { # single end
			out_sample_sheet = rbind(out_sample_sheet,
				data.frame(
					fq_path=out_sample_info['fq_path'][i, ],
					request_id=out_sample_info['request_id'][i, ],
					sample_id=out_sample_info['sample_id'][i, ],
					srr_id=out_sample_info['srr_id'][i, ],
					sample_name=out_sample_info['sample_name'][i, ],
					expr_type=out_sample_info['expr_type'][i, ],
					read=1,
					read_type=out_sample_info['read_type'][i, ],
					read_length=out_sample_info['read_length'][i, ] + 1,
					ref_genome=out_sample_info['ref_genome'][i, ],
					lab=out_sample_info['lab'][i, ],
					requester=out_sample_info['requester'][i, ]
				)
			)
		}
	}

	return(out_sample_sheet)
}
###########################################################################################################


################################################## main ###################################################
main = function() {
	suppressPackageStartupMessages(library(foreach))
	suppressPackageStartupMessages(library(doParallel))

	# get opts
	opts = get_opts()

	sample_sheet_file = opts$sample_sheet
	requester = opts$requester
	lab = opts$lab
	out_dir = opts$out
	threads = opts$threads
	sra_cache_folder = opts$temp
	api_key = opts$api_key
	dry_run = !as.logical(opts$download_mode)

	if(!dry_run) {
		# clear sta cache folder
		message(sprintf('[%s] Program in downloading mode. First clear SRA cache folder: %s...', 
			Sys.time(), sra_cache_folder))
		do.call(file.remove, list(list.files(sra_cache_folder, full.names=TRUE)))
	} else {
		message(sprintf('[%s] Program in dry-run mode', Sys.time()))
	}

	# load sample info
	message(sprintf('[%s] Loading sample sheet...', Sys.time()))
	sample_info = load_sample_info(sample_sheet_file, requester, lab)

	# out dir
	out_dir = normalizePath(out_dir) # absolute path

	# create if not exists
	if(!dir.exists(out_dir)) {
		dir.create(out_dir, recursive=TRUE)
	}

	# create output sample sheet folder if not exists
	if(!dir.exists(file.path(out_dir, 'analysis_sample_sheet'))) {
		dir.create(file.path(out_dir, 'analysis_sample_sheet'), recursive=TRUE)
	}

	# subset sample sheet
	id = as.data.frame(cbind(sample_info$sample_sheet['request_id'], sample_info$sample_sheet['sample_id']))

	# generate commands
	message(sprintf('[%s] Generating commands for fastq downloading...', Sys.time()))
	cmds = data.frame()
	for(i in seq(1, nrow(id))) {
		request_id = id[i, 'request_id']
		sample_id = id[i, 'sample_id']

		cmds = rbind(cmds, download_cmds_wrapper(request_id, sample_id, out_dir, sra_cache_folder, api_key))
	}

	message(sprintf('[%s] ID mapping summary', Sys.time()))
	print(cmds[c('request_id', 'sample_id', 'srr_id')])

	# write commands to file
	# write(paste(as.character(cmds[, 'cmds']), collapse='\n'), file=file.path(out_dir, 'analysis_sample_sheet', 
	# 	sprintf('%s_%s_cmds.sh', request_id, sample_info$requester)))

	write(paste(as.character(cmds[, 'cmds']), collapse='\n'), file=file.path("./", 
		sprintf('%s_%s_cmds.sh', request_id, sample_info$requester)))

	# generating sample info for analysis
	message(sprintf('[%s] Generating analysis sample sheet...', Sys.time()))
	out_sample_info = get_sample_info_analysis(sample_info, cmds, out_dir)
	out_sample_info = out_sample_info[with(out_sample_info, order(request_id, sample_id)), ]

	# write.table(out_sample_info, file.path(out_dir, 'analysis_sample_sheet', sprintf('%s_%s_samplesheet.csv', 
	# 	request_id, sample_info$requester)), sep=',', row.names=FALSE, col.names=TRUE, quote=FALSE)
	
	write.table(out_sample_info, file.path('./', sprintf('%s_%s_samplesheet.csv', 
		request_id, sample_info$requester)), sep=',', row.names=FALSE, col.names=TRUE, quote=FALSE)

	# download fastq if not dry run
	if(!dry_run) {
		message(sprintf('[%s] Downloading fastq files...', Sys.time()))

		# set up parallel backend
		if(nrow(cmds) < threads) {
			message(paste(sprintf('[%s] Number of parallel taks (%d) is less than requested (%d).', 
				Sys.time(), nrow(cmds), threads), sprintf('Will use %d cpu processes', nrow(cmds))))
		}
		registerDoParallel(cores=min(threads, nrow(cmds)))

		# parallel processing
		run = foreach(row=iter(cmds, by='row'), .combine='rbind') %dopar% {
			request_id = as.character(row['request_id'][[1]])
			sample_id = as.character(row['sample_id'][[1]])
			srr_id = as.character(row['srr_id'][[1]])

			log = data.frame(
				request_id=as.character(row['request_id'][[1]]), 
				sample_id=as.character(row['sample_id'][[1]]), 
				srr_id=as.character(row['srr_id'][[1]]), 
				return_code=0
			)

			# skip sra download
			skip_sra = FALSE
			if(length(list.files(file.path(out_dir, request_id, sample_id, srr_id), '.sra')) > 0) {
				message(sprintf('[%s] %s/%s/%s: sra file found. Check validity', 
					Sys.time(), request_id, sample_id, srr_id))

				check_sra_cmd = sprintf('vdb-validate %s 2>&1', 
					file.path(out_dir, request_id, sample_id, srr_id, sprintf('%s.sra', srr_id)))
				check_sra_result = paste(system(check_sra_cmd, intern=TRUE), collapse=' ')

				if(grepl(sprintf('Table \'%s.sra\' is consistent', srr_id), check_sra_result)) {
					message(sprintf('[%s] %s/%s/%s: sra file is consistent. Skip.', 
						Sys.time(), request_id, sample_id, srr_id))
					skip_sra = TRUE
				}
			}
			
			# skip fastq dump
			skip_fastq = FALSE
			if(length(list.files(file.path(out_dir, request_id, sample_id, srr_id), '.fastq.gz')) > 0) {
				message(sprintf('[%s] %s/%s/%s: fastq files found. Checking validity...', 
					Sys.time(), request_id, sample_id, srr_id))

				check_fastq_cmd = sprintf('gzip -t %s && echo ok || echo bad 2>&1', 
					file.path(out_dir, request_id, sample_id, srr_id, '*.fastq.gz'))				
				check_fastq_result = system(check_fastq_cmd, intern=TRUE)

				if(check_fastq_result == 'ok') {
					message(sprintf('[%s] %s/%s/%s: fastq files are not corrupted. Skip.', 
						Sys.time(), request_id, sample_id, srr_id))
					skip_fastq = TRUE
					skip_sra = TRUE
				}
			}

			cmd_list = strsplit(as.character(row['cmds'][[1]]), '\n')[[1]]
			return_codes = c()

			for(cmd in cmd_list) {
				if(!startsWith(cmd, '#')) { # not comment line
					if (skip_sra & !skip_fastq) { # skip sra but not fastq
						if(startsWith(cmd, 'prefetch') | startsWith(cmd, 'mv') | startsWith(cmd, 'vdb-validate')) {
							r = 0
						} else { # only do fastq dump step
							message(sprintf('[%s] %s', Sys.time(), cmd))
							r = system(cmd)
							message(sprintf('[%s] %s: returns %d', Sys.time(), cmd, r))
						}
					} else if(skip_fastq) { # skip both sra and fastq
						r = 0
					} else { # new download
						message(sprintf('[%s] %s', Sys.time(), cmd))
						r = system(cmd)
						message(sprintf('[%s] %s: returns %d', Sys.time(), cmd, r))
					}

					return_codes = c(return_codes, r)
				}
			}

			# success return 0, otherwise return -1
			log['return_code'] = as.numeric(sum(return_codes^2) == 0) - 1

			log # 'return' value for foreach dopar
		}

		message(sprintf('[%s] Download summary ', Sys.time()))
		print(run)

	} else {
		message(sprintf('[%s] Program exits without downloading files since option --dryrun_mode is set...', 
			Sys.time()))
	}
}
###########################################################################################################

main()
